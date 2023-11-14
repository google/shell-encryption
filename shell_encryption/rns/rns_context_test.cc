/*
 * Copyright 2023 Google LLC.
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     https://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include "shell_encryption/rns/rns_context.h"

#include <memory>
#include <random>

#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include "absl/strings/str_cat.h"
#include "shell_encryption/integral_types.h"
#include "shell_encryption/montgomery.h"
#include "shell_encryption/ntt_parameters.h"
#include "shell_encryption/rns/testing/parameters.h"
#include "shell_encryption/rns/testing/testing_utils.h"
#include "shell_encryption/testing/parameters.h"
#include "shell_encryption/testing/status_matchers.h"
#include "shell_encryption/testing/status_testing.h"

namespace rlwe {
namespace {

// Useful typedefs.
using ModularInt16 = MontgomeryInt<Uint16>;
using ModularInt32 = MontgomeryInt<Uint32>;
using ModularInt64 = MontgomeryInt<Uint64>;
using ModularInt128 = MontgomeryInt<Uint128>;

using ::rlwe::testing::StatusIs;
using ::testing::HasSubstr;
using ::testing::SizeIs;

// Test fixtures.
template <typename T>
class RnsContextTest : public ::testing::Test {};
TYPED_TEST_SUITE(RnsContextTest, testing::ModularIntTypes);

TYPED_TEST(RnsContextTest, CreateFailsIfLogNIsNonPositive) {
  EXPECT_THAT(RnsContext<TypeParam>::Create(/*log_n=*/-1, {}, {}, 2),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`log_n` must be positive")));
  EXPECT_THAT(RnsContext<TypeParam>::Create(/*log_n=*/0, {}, {}, 2),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`log_n` must be positive")));
}

TYPED_TEST(RnsContextTest, CreateFailsIfPlaintextModulusIsZero) {
  EXPECT_THAT(RnsContext<TypeParam>::Create(1, {}, {}, /*plaintext_modulus=*/0),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`plaintext_modulus` must be positive")));
}

// Checks RnsContext can be created with the chosen constants.
TYPED_TEST(RnsContextTest, CreateSucceeds) {
  testing::RnsParameters<TypeParam> rns_params =
      testing::GetRnsParameters<TypeParam>();

  ASSERT_OK_AND_ASSIGN(
      RnsContext<TypeParam> context,
      RnsContext<TypeParam>::Create(rns_params.log_n, rns_params.qs,
                                    rns_params.ps, rns_params.t));
  EXPECT_EQ(context.LogN(), rns_params.log_n);
  EXPECT_EQ(context.NumMainPrimeModuli(), rns_params.qs.size());
  EXPECT_EQ(context.NumAuxPrimeModuli(), rns_params.ps.size());

  std::vector<const PrimeModulus<TypeParam>*> main_moduli =
      context.MainPrimeModuli();
  EXPECT_EQ(main_moduli.size(), rns_params.qs.size());
  for (int i = 0; i < rns_params.qs.size(); ++i) {
    EXPECT_NE(main_moduli[i]->ModParams(), nullptr);
    EXPECT_NE(main_moduli[i]->NttParams(), nullptr);
    EXPECT_EQ(main_moduli[i]->Modulus(), rns_params.qs[i]);
    EXPECT_EQ(main_moduli[i]->NttParams()->number_coeffs,
              (1 << rns_params.log_n));
  }

  std::vector<const PrimeModulus<TypeParam>*> aux_moduli =
      context.AuxPrimeModuli();
  EXPECT_EQ(aux_moduli.size(), rns_params.ps.size());
  for (int i = 0; i < rns_params.ps.size(); ++i) {
    EXPECT_NE(aux_moduli[i]->ModParams(), nullptr);
    EXPECT_NE(aux_moduli[i]->NttParams(), nullptr);
    EXPECT_EQ(aux_moduli[i]->Modulus(), rns_params.ps[i]);
    EXPECT_EQ(aux_moduli[i]->NttParams()->number_coeffs,
              (1 << rns_params.log_n));
  }

  EXPECT_EQ(context.PlaintextModulus(), rns_params.t);
}

// Checks RnsContext can be created even without any auxiliary modulus.
TYPED_TEST(RnsContextTest, CreateSucceedsWithEmptyAuxModuli) {
  testing::RnsParameters<TypeParam> rns_params =
      testing::GetRnsParameters<TypeParam>();
  ASSERT_OK_AND_ASSIGN(
      RnsContext<TypeParam> context,
      RnsContext<TypeParam>::Create(rns_params.log_n, rns_params.qs, /*ps=*/{},
                                    rns_params.t));
  EXPECT_EQ(context.LogN(), rns_params.log_n);
  EXPECT_EQ(context.NumMainPrimeModuli(), rns_params.qs.size());
  EXPECT_EQ(context.NumAuxPrimeModuli(), 0);

  std::vector<const PrimeModulus<TypeParam>*> main_moduli =
      context.MainPrimeModuli();
  EXPECT_EQ(main_moduli.size(), rns_params.qs.size());
  for (int i = 0; i < rns_params.qs.size(); ++i) {
    EXPECT_NE(main_moduli[i]->ModParams(), nullptr);
    EXPECT_NE(main_moduli[i]->NttParams(), nullptr);
    EXPECT_EQ(main_moduli[i]->Modulus(), rns_params.qs[i]);
    EXPECT_EQ(main_moduli[i]->NttParams()->number_coeffs,
              (1 << rns_params.log_n));
  }

  EXPECT_THAT(context.AuxPrimeModuli(), SizeIs(0));

  EXPECT_EQ(context.PlaintextModulus(), rns_params.t);
}

// Checks that constants precomputation in RnsContext are done correctly.
TYPED_TEST(RnsContextTest, CrtConstantsAreCorrect) {
  testing::RnsParameters<TypeParam> rns_params =
      testing::GetRnsParameters<TypeParam>();
  ASSERT_OK_AND_ASSIGN(
      RnsContext<TypeParam> context,
      RnsContext<TypeParam>::Create(rns_params.log_n, rns_params.qs,
                                    rns_params.ps, rns_params.t));

  std::vector<const PrimeModulus<TypeParam>*> main_moduli =
      context.MainPrimeModuli();
  std::vector<const PrimeModulus<TypeParam>*> aux_moduli =
      context.AuxPrimeModuli();

  // Expect q_mod_qs[i] == {qs[i] % qs[j] : j = 0..L}.
  auto q_mod_qs = context.MainPrimeModulusResidues();
  EXPECT_EQ(q_mod_qs.size(), rns_params.qs.size());
  for (int i = 0; i < rns_params.qs.size(); ++i) {
    EXPECT_EQ(q_mod_qs[i].NumModuli(), rns_params.qs.size());
    for (int j = 0; j < rns_params.qs.size(); ++j) {
      EXPECT_EQ(q_mod_qs[i].Component(j).ExportInt(main_moduli[j]->ModParams()),
                rns_params.qs[i] % rns_params.qs[j]);
    }
  }

  // Expect q_inv_mod_qs[i] == {qs[i]^(-1) % qs[j] : j = 0..L}.
  auto q_inv_mod_qs = context.MainPrimeModulusInverseResidues();
  EXPECT_EQ(q_inv_mod_qs.size(), rns_params.qs.size());
  for (int i = 0; i < rns_params.qs.size(); ++i) {
    EXPECT_EQ(q_inv_mod_qs[i].NumModuli(), rns_params.qs.size());
    for (int j = 0; j < rns_params.qs.size(); ++j) {
      if (i != j) {
        EXPECT_EQ(q_inv_mod_qs[i].Component(j).Mul(q_mod_qs[i].Component(j),
                                                   main_moduli[j]->ModParams()),
                  TypeParam::ImportOne(main_moduli[j]->ModParams()));
      }
    }
  }

  // Expect q_inv_mod_ps[i] = {qs[j]^(-1) % ps[i] : j = 0..L}.
  auto q_inv_mod_ps = context.MainPrimeModulusInverseAuxResidues();
  ASSERT_EQ(q_inv_mod_ps.size(), rns_params.ps.size());
  for (int i = 0; i < rns_params.ps.size(); ++i) {
    ASSERT_EQ(q_inv_mod_ps[i].NumModuli(), rns_params.qs.size());
    auto mod_params_pi = aux_moduli[i]->ModParams();
    for (int j = 0; j < rns_params.qs.size(); ++j) {
      ASSERT_OK_AND_ASSIGN(
          auto qj_mod_pi,
          TypeParam::ImportInt(rns_params.qs[j], mod_params_pi));
      // Expect q_inv_mod_ps[i][j] = qs[j]^(-1) % ps[i].
      EXPECT_EQ(q_inv_mod_ps[i].Component(j).Mul(qj_mod_pi, mod_params_pi),
                TypeParam::ImportOne(mod_params_pi));
    }
  }

  // Expect ql_mod_ps[i] = {\prod_{j=0}^l qs[j] % ps[i] : l = 0..L}.
  std::vector<TypeParam> qs_mod_ps;
  for (int j = 0; j < rns_params.ps.size(); ++j) {
    qs_mod_ps.push_back(TypeParam::ImportOne(aux_moduli[j]->ModParams()));
  }
  for (int i = 0; i < rns_params.qs.size(); ++i) {
    ASSERT_OK_AND_ASSIGN(auto ql_mod_ps,
                         context.MainLeveledModulusAuxResidues(i));
    ASSERT_EQ(ql_mod_ps.NumModuli(), rns_params.ps.size());
    for (int j = 0; j < rns_params.ps.size(); ++j) {
      auto mod_params_pj = aux_moduli[j]->ModParams();
      ASSERT_OK_AND_ASSIGN(
          auto qi_mod_pj,
          TypeParam::ImportInt(rns_params.qs[i], mod_params_pj));
      EXPECT_EQ(ql_mod_ps.Component(j),
                qs_mod_ps[j].MulInPlace(qi_mod_pj, mod_params_pj));
    }
  }

  // Expect t_mod_qs = t % qs[i].
  auto t_mod_qs = context.PlaintextModulusMainResidues();
  EXPECT_EQ(t_mod_qs.size(), rns_params.qs.size());
  for (int i = 0; i < rns_params.qs.size(); ++i) {
    EXPECT_EQ(t_mod_qs[i].ExportInt(main_moduli[i]->ModParams()),
              rns_params.t % rns_params.qs[i]);
  }
}

TYPED_TEST(RnsContextTest, CrtConstantsAreCorrectForArbitraryPlaintextModulus) {
  // Use a different plaintext modulus in this test. The constants should still
  // be generated correctly.
  testing::RnsParameters<TypeParam> rns_params =
      testing::GetRnsParameters<TypeParam>();
  rns_params.t += 1;

  ASSERT_OK_AND_ASSIGN(
      RnsContext<TypeParam> context,
      RnsContext<TypeParam>::Create(rns_params.log_n, rns_params.qs,
                                    rns_params.ps, rns_params.t));

  std::vector<const PrimeModulus<TypeParam>*> main_moduli =
      context.MainPrimeModuli();
  std::vector<const PrimeModulus<TypeParam>*> aux_moduli =
      context.AuxPrimeModuli();

  // Expect t_mod_qs = t % qs[i].
  auto t_mod_qs = context.PlaintextModulusMainResidues();
  EXPECT_EQ(t_mod_qs.size(), rns_params.qs.size());
  for (int i = 0; i < rns_params.qs.size(); ++i) {
    EXPECT_EQ(t_mod_qs[i].ExportInt(main_moduli[i]->ModParams()),
              rns_params.t % rns_params.qs[i]);
  }
}

TYPED_TEST(RnsContextTest,
           MainLeveledModulusAuxResiduesErrorsIfLevelIsOutOfRange) {
  testing::RnsParameters<TypeParam> rns_params =
      testing::GetRnsParameters<TypeParam>();
  ASSERT_OK_AND_ASSIGN(
      RnsContext<TypeParam> context,
      RnsContext<TypeParam>::Create(rns_params.log_n, rns_params.qs,
                                    rns_params.ps, rns_params.t));

  // `level` is negative.
  EXPECT_THAT(context.MainLeveledModulusAuxResidues(/*level=*/-1),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`level` must be non-negative")));

  // `level` is larger than the number of main moduli - 1.
  EXPECT_THAT(
      context.MainLeveledModulusAuxResidues(/*level=*/rns_params.qs.size()),
      StatusIs(
          absl::StatusCode::kInvalidArgument,
          HasSubstr(absl::StrCat("`level` must be non-negative and at most ",
                                 rns_params.qs.size() - 1))));
}

TYPED_TEST(RnsContextTest,
           MainPrimeModulusComplementsErrorsIfLevelIsOutOfRange) {
  // Create a RNS context with just main prime moduli.
  testing::RnsParameters<TypeParam> rns_params =
      testing::GetRnsParameters<TypeParam>();
  ASSERT_OK_AND_ASSIGN(
      RnsContext<TypeParam> context,
      RnsContext<TypeParam>::Create(rns_params.log_n, rns_params.qs, /*ps=*/{},
                                    rns_params.t));

  // `level` is negative.
  EXPECT_THAT(context.MainPrimeModulusComplements(/*level=*/-1),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`level` must be non-negative")));

  // `level` is larger than the number of main moduli - 1.
  EXPECT_THAT(
      context.MainPrimeModulusComplements(/*level=*/rns_params.qs.size()),
      StatusIs(
          absl::StatusCode::kInvalidArgument,
          HasSubstr(absl::StrCat("`level` must be non-negative and at most ",
                                 rns_params.qs.size() - 1))));
}

TYPED_TEST(RnsContextTest,
           MainPrimeModulusCrtFactorsErrorsIfLevelIsOutOfRange) {
  // Create a RNS context with just main prime moduli.
  testing::RnsParameters<TypeParam> rns_params =
      testing::GetRnsParameters<TypeParam>();
  ASSERT_OK_AND_ASSIGN(
      RnsContext<TypeParam> context,
      RnsContext<TypeParam>::Create(rns_params.log_n, rns_params.qs, /*ps=*/{},
                                    rns_params.t));

  // `level` is negative.
  EXPECT_THAT(context.MainPrimeModulusCrtFactors(/*level=*/-1),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`level` must be non-negative")));

  // `level` is larger than the number of main moduli - 1.
  EXPECT_THAT(
      context.MainPrimeModulusCrtFactors(/*level=*/rns_params.qs.size()),
      StatusIs(
          absl::StatusCode::kInvalidArgument,
          HasSubstr(absl::StrCat("`level` must be non-negative and at most ",
                                 rns_params.qs.size() - 1))));
}

TYPED_TEST(RnsContextTest,
           MainPrimeModulusComplementResiduesErrorsIfLevelIsOutOfRange) {
  // Create a RNS context with just main prime moduli.
  testing::RnsParameters<TypeParam> rns_params =
      testing::GetRnsParameters<TypeParam>();
  ASSERT_OK_AND_ASSIGN(
      RnsContext<TypeParam> context,
      RnsContext<TypeParam>::Create(rns_params.log_n, rns_params.qs, /*ps=*/{},
                                    rns_params.t));

  // `level` is negative.
  EXPECT_THAT(context.MainPrimeModulusComplementResidues(/*level=*/-1),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`level` must be non-negative")));

  // `level` is larger than the number of main moduli.
  EXPECT_THAT(context.MainPrimeModulusComplementResidues(
                  /*level=*/rns_params.qs.size()),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr(absl::StrCat(
                           "`level` must be non-negative and at most ",
                           rns_params.qs.size() - 1))));
}

// Checks RnsContext can be created to support BFV encoding.
TYPED_TEST(RnsContextTest, CreateForBfvSucceeds) {
  testing::RnsParameters<TypeParam> rns_params =
      testing::GetRnsParameters<TypeParam>();
  // `CreateForBfv` currently supports only an odd plaintext modulus in order to
  // instantiate Montgomery integer modulo the plaintext modulus.
  if (rns_params.t % 2 == 0) {
    rns_params.t += 1;
  }

  ASSERT_OK_AND_ASSIGN(
      RnsContext<TypeParam> context,
      RnsContext<TypeParam>::CreateForBfv(rns_params.log_n, rns_params.qs,
                                          rns_params.ps, rns_params.t));
  EXPECT_EQ(context.LogN(), rns_params.log_n);
  EXPECT_EQ(context.NumMainPrimeModuli(), rns_params.qs.size());
  EXPECT_EQ(context.NumAuxPrimeModuli(), rns_params.ps.size());

  std::vector<const PrimeModulus<TypeParam>*> main_moduli =
      context.MainPrimeModuli();
  EXPECT_EQ(main_moduli.size(), rns_params.qs.size());
  for (int i = 0; i < rns_params.qs.size(); ++i) {
    EXPECT_NE(main_moduli[i]->ModParams(), nullptr);
    EXPECT_NE(main_moduli[i]->NttParams(), nullptr);
    EXPECT_EQ(main_moduli[i]->Modulus(), rns_params.qs[i]);
    EXPECT_EQ(main_moduli[i]->NttParams()->number_coeffs,
              (1 << rns_params.log_n));
  }
  std::vector<const PrimeModulus<TypeParam>*> aux_moduli =
      context.AuxPrimeModuli();
  EXPECT_EQ(aux_moduli.size(), rns_params.ps.size());
  for (int i = 0; i < rns_params.ps.size(); ++i) {
    EXPECT_NE(aux_moduli[i]->ModParams(), nullptr);
    EXPECT_NE(aux_moduli[i]->NttParams(), nullptr);
    EXPECT_EQ(aux_moduli[i]->Modulus(), rns_params.ps[i]);
    EXPECT_EQ(aux_moduli[i]->NttParams()->number_coeffs,
              (1 << rns_params.log_n));
  }

  EXPECT_EQ(context.PlaintextModulus(), rns_params.t);
  EXPECT_NE(context.PlaintextModulusParams().ModParams(), nullptr);
  EXPECT_EQ(context.PlaintextModulusParams().Modulus(), rns_params.t);
}

TYPED_TEST(RnsContextTest, PlaintextModulusConstantsForBfvAreCorrect) {
  testing::RnsParameters<TypeParam> rns_params =
      testing::GetRnsParameters<TypeParam>();
  // `CreateForBfv` currently supports only an odd plaintext modulus in order to
  // instantiate Montgomery integer modulo the plaintext modulus.
  if (rns_params.t % 2 == 0) {
    rns_params.t += 1;
  }

  ASSERT_OK_AND_ASSIGN(
      RnsContext<TypeParam> context,
      RnsContext<TypeParam>::CreateForBfv(rns_params.log_n, rns_params.qs,
                                          rns_params.ps, rns_params.t));
  auto mod_params_t = context.PlaintextModulusParams().ModParams();
  ASSERT_NE(mod_params_t, nullptr);

  // Expect qs_mod_t[i] = q_i (mod t), q_invs_mod_t[i] = q_i^(-1) (mod t),
  // ql_invs_mod_t[l] = Q_l^(-1) (mod t) = \prod_{i=0}^l q_i^(-1) (mod t).
  auto qs_mod_t = context.MainPrimeModulusPlaintextResidues();
  auto q_invs_mod_t = context.MainPrimeModulusInversePlaintextResidues();
  auto ql_invs_mod_t = context.MainLeveledModulusInversePlaintextResidue();
  ASSERT_EQ(qs_mod_t.size(), rns_params.qs.size());
  ASSERT_EQ(q_invs_mod_t.size(), rns_params.qs.size());
  ASSERT_EQ(ql_invs_mod_t.size(), rns_params.qs.size());

  auto ql_inv_mod_t = TypeParam::ImportOne(mod_params_t);  // Q_l^(-1) (mod t).
  for (int i = 0; i < rns_params.qs.size(); ++i) {
    EXPECT_EQ(qs_mod_t[i].ExportInt(mod_params_t),
              rns_params.qs[i] % rns_params.t);

    ASSERT_OK_AND_ASSIGN(
        TypeParam qi_mod_t,
        TypeParam::ImportInt(rns_params.qs[i] % rns_params.t, mod_params_t));
    EXPECT_EQ(qi_mod_t.Mul(q_invs_mod_t[i], mod_params_t),
              TypeParam::ImportOne(mod_params_t));

    ASSERT_OK_AND_ASSIGN(TypeParam qi_inv_mod_t,
                         qi_mod_t.MultiplicativeInverseFast(mod_params_t));
    ql_inv_mod_t.MulInPlace(qi_inv_mod_t, mod_params_t);
    EXPECT_EQ(ql_invs_mod_t[i], ql_inv_mod_t);
  }
}

////////////////////////////////////////////////////////////////////////////////
// Tests for supporting finite field encoding
////////////////////////////////////////////////////////////////////////////////

// Test fixtures.
template <typename T>
class RnsContextForFiniteFieldEncodingTest : public ::testing::Test {};
TYPED_TEST_SUITE(RnsContextForFiniteFieldEncodingTest,
                 testing::ModularIntTypesForFiniteFieldEncoding);

TYPED_TEST(RnsContextForFiniteFieldEncodingTest,
           CreateForBgvFiniteFieldEncodingFailsIfPlaintextModulusIsEven) {
  EXPECT_THAT(
      RnsContext<TypeParam>::CreateForBgvFiniteFieldEncoding(
          1, {}, {}, /*plaintext_modulus=*/0),
      StatusIs(absl::StatusCode::kInvalidArgument,
               HasSubstr("`plaintext_modulus` cannot be an even number")));
  EXPECT_THAT(
      RnsContext<TypeParam>::CreateForBgvFiniteFieldEncoding(
          1, {}, {}, /*plaintext_modulus=*/2),
      StatusIs(absl::StatusCode::kInvalidArgument,
               HasSubstr("`plaintext_modulus` cannot be an even number")));
  EXPECT_THAT(
      RnsContext<TypeParam>::CreateForBgvFiniteFieldEncoding(
          1, {}, {}, /*plaintext_modulus=*/6),
      StatusIs(absl::StatusCode::kInvalidArgument,
               HasSubstr("`plaintext_modulus` cannot be an even number")));
}

TYPED_TEST(RnsContextForFiniteFieldEncodingTest,
           CreateForBgvFiniteFieldEncodingSucceeds) {
  for (auto const& rns_params :
       testing::GetRnsParametersForFiniteFieldEncoding<TypeParam>()) {
    ASSERT_OK_AND_ASSIGN(
        RnsContext<TypeParam> context,
        RnsContext<TypeParam>::CreateForBgvFiniteFieldEncoding(
            rns_params.log_n, rns_params.qs, rns_params.ps, rns_params.t));
    EXPECT_EQ(context.LogN(), rns_params.log_n);
    EXPECT_EQ(context.NumMainPrimeModuli(), rns_params.qs.size());
    EXPECT_EQ(context.NumAuxPrimeModuli(), rns_params.ps.size());

    std::vector<const PrimeModulus<TypeParam>*> main_moduli =
        context.MainPrimeModuli();
    EXPECT_EQ(main_moduli.size(), rns_params.qs.size());
    for (int i = 0; i < rns_params.qs.size(); ++i) {
      EXPECT_NE(main_moduli[i]->ModParams(), nullptr);
      EXPECT_NE(main_moduli[i]->NttParams(), nullptr);
      EXPECT_EQ(main_moduli[i]->Modulus(), rns_params.qs[i]);
      EXPECT_EQ(main_moduli[i]->NttParams()->number_coeffs,
                (1 << rns_params.log_n));
    }
    std::vector<const PrimeModulus<TypeParam>*> aux_moduli =
        context.AuxPrimeModuli();
    EXPECT_EQ(aux_moduli.size(), rns_params.ps.size());
    for (int i = 0; i < rns_params.ps.size(); ++i) {
      EXPECT_NE(aux_moduli[i]->ModParams(), nullptr);
      EXPECT_NE(aux_moduli[i]->NttParams(), nullptr);
      EXPECT_EQ(aux_moduli[i]->Modulus(), rns_params.ps[i]);
      EXPECT_EQ(aux_moduli[i]->NttParams()->number_coeffs,
                (1 << rns_params.log_n));
    }

    EXPECT_EQ(context.PlaintextModulus(), rns_params.t);
    const PrimeModulus<TypeParam>& modulus_t = context.PlaintextModulusParams();
    EXPECT_NE(modulus_t.ModParams(), nullptr);
    EXPECT_EQ(modulus_t.Modulus(), rns_params.t);
    EXPECT_NE(modulus_t.NttParams(), nullptr);
    EXPECT_EQ(modulus_t.NttParams()->number_coeffs, (1 << rns_params.log_n));
  }
}

TYPED_TEST(RnsContextForFiniteFieldEncodingTest,
           CreateForBfvFiniteFieldEncodingFailsIfPlaintextModulusIsEven) {
  EXPECT_THAT(
      RnsContext<TypeParam>::CreateForBfvFiniteFieldEncoding(
          1, {}, {}, /*plaintext_modulus=*/0),
      StatusIs(absl::StatusCode::kInvalidArgument,
               HasSubstr("`plaintext_modulus` cannot be an even number")));
  EXPECT_THAT(
      RnsContext<TypeParam>::CreateForBfvFiniteFieldEncoding(
          1, {}, {}, /*plaintext_modulus=*/2),
      StatusIs(absl::StatusCode::kInvalidArgument,
               HasSubstr("`plaintext_modulus` cannot be an even number")));
  EXPECT_THAT(
      RnsContext<TypeParam>::CreateForBfvFiniteFieldEncoding(
          1, {}, {}, /*plaintext_modulus=*/6),
      StatusIs(absl::StatusCode::kInvalidArgument,
               HasSubstr("`plaintext_modulus` cannot be an even number")));
}

TYPED_TEST(RnsContextForFiniteFieldEncodingTest,
           CreateForBfvFiniteFieldEncodingSucceeds) {
  for (auto const& rns_params :
       testing::GetRnsParametersForFiniteFieldEncoding<TypeParam>()) {
    ASSERT_OK_AND_ASSIGN(
        RnsContext<TypeParam> context,
        RnsContext<TypeParam>::CreateForBfvFiniteFieldEncoding(
            rns_params.log_n, rns_params.qs, rns_params.ps, rns_params.t));
    EXPECT_EQ(context.LogN(), rns_params.log_n);
    EXPECT_EQ(context.NumMainPrimeModuli(), rns_params.qs.size());
    EXPECT_EQ(context.NumAuxPrimeModuli(), rns_params.ps.size());

    std::vector<const PrimeModulus<TypeParam>*> main_moduli =
        context.MainPrimeModuli();
    EXPECT_EQ(main_moduli.size(), rns_params.qs.size());
    for (int i = 0; i < rns_params.qs.size(); ++i) {
      EXPECT_NE(main_moduli[i]->ModParams(), nullptr);
      EXPECT_NE(main_moduli[i]->NttParams(), nullptr);
      EXPECT_EQ(main_moduli[i]->Modulus(), rns_params.qs[i]);
      EXPECT_EQ(main_moduli[i]->NttParams()->number_coeffs,
                (1 << rns_params.log_n));
    }
    std::vector<const PrimeModulus<TypeParam>*> aux_moduli =
        context.AuxPrimeModuli();
    EXPECT_EQ(aux_moduli.size(), rns_params.ps.size());
    for (int i = 0; i < rns_params.ps.size(); ++i) {
      EXPECT_NE(aux_moduli[i]->ModParams(), nullptr);
      EXPECT_NE(aux_moduli[i]->NttParams(), nullptr);
      EXPECT_EQ(aux_moduli[i]->Modulus(), rns_params.ps[i]);
      EXPECT_EQ(aux_moduli[i]->NttParams()->number_coeffs,
                (1 << rns_params.log_n));
    }

    EXPECT_EQ(context.PlaintextModulus(), rns_params.t);
    const PrimeModulus<TypeParam>& modulus_t = context.PlaintextModulusParams();
    EXPECT_NE(modulus_t.ModParams(), nullptr);
    EXPECT_EQ(modulus_t.Modulus(), rns_params.t);
    EXPECT_NE(modulus_t.NttParams(), nullptr);
    EXPECT_EQ(modulus_t.NttParams()->number_coeffs, (1 << rns_params.log_n));
  }
}

TYPED_TEST(RnsContextForFiniteFieldEncodingTest,
           PlaintextModulusConstantsForBfvAreCorrect) {
  for (auto const& rns_params :
       testing::GetRnsParametersForFiniteFieldEncoding<TypeParam>()) {
    ASSERT_OK_AND_ASSIGN(
        RnsContext<TypeParam> context,
        RnsContext<TypeParam>::CreateForBfv(rns_params.log_n, rns_params.qs,
                                            rns_params.ps, rns_params.t));
    auto mod_params_t = context.PlaintextModulusParams().ModParams();
    ASSERT_NE(mod_params_t, nullptr);

    // Expect qs_mod_t[i] = q_i (mod t), q_invs_mod_t[i] = q_i^(-1) (mod t),
    // ql_invs_mod_t[l] = Q_l^(-1) (mod t) = \prod_{i=0}^l q_i^(-1) (mod t).
    auto qs_mod_t = context.MainPrimeModulusPlaintextResidues();
    auto q_invs_mod_t = context.MainPrimeModulusInversePlaintextResidues();
    auto ql_invs_mod_t = context.MainLeveledModulusInversePlaintextResidue();
    ASSERT_EQ(qs_mod_t.size(), rns_params.qs.size());
    ASSERT_EQ(q_invs_mod_t.size(), rns_params.qs.size());
    ASSERT_EQ(ql_invs_mod_t.size(), rns_params.qs.size());

    auto ql_inv_mod_t =
        TypeParam::ImportOne(mod_params_t);  // Q_l^(-1) (mod t).
    for (int i = 0; i < rns_params.qs.size(); ++i) {
      EXPECT_EQ(qs_mod_t[i].ExportInt(mod_params_t),
                rns_params.qs[i] % rns_params.t);

      ASSERT_OK_AND_ASSIGN(
          TypeParam qi_mod_t,
          TypeParam::ImportInt(rns_params.qs[i] % rns_params.t, mod_params_t));
      EXPECT_EQ(qi_mod_t.Mul(q_invs_mod_t[i], mod_params_t),
                TypeParam::ImportOne(mod_params_t));

      ASSERT_OK_AND_ASSIGN(TypeParam qi_inv_mod_t,
                           qi_mod_t.MultiplicativeInverseFast(mod_params_t));
      ql_inv_mod_t.MulInPlace(qi_inv_mod_t, mod_params_t);
      EXPECT_EQ(ql_invs_mod_t[i], ql_inv_mod_t);
    }
  }
}

}  // namespace
}  // namespace rlwe
