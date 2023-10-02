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

#include "shell_encryption/rns/crt_interpolation.h"

#include <vector>

#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include "absl/numeric/int128.h"
#include "absl/random/random.h"
#include "absl/status/status.h"
#include "shell_encryption/int256.h"
#include "shell_encryption/integral_types.h"
#include "shell_encryption/montgomery.h"
#include "shell_encryption/rns/rns_context.h"
#include "shell_encryption/rns/testing/parameters.h"
#include "shell_encryption/status_macros.h"
#include "shell_encryption/testing/parameters.h"
#include "shell_encryption/testing/status_matchers.h"
#include "shell_encryption/testing/status_testing.h"

namespace rlwe {
namespace {

using ::rlwe::testing::StatusIs;
using ::testing::HasSubstr;

// Test fixture to pre-compute constants for RNS moduli.
template <typename ModularInt>
class CrtInterpolationTest : public ::testing::Test {
 protected:
  using ModularIntParams = typename ModularInt::Params;
  using Integer = typename ModularInt::Int;
  using BigInteger = typename testing::CompositeModulus<Integer>::value_type;

  void SetUp() override {
    testing::RnsParameters<ModularInt> rns_params =
        testing::GetRnsParameters<ModularInt>();
    ASSERT_OK_AND_ASSIGN(auto rns_context, RnsContext<ModularInt>::Create(
                                               rns_params.log_n, rns_params.qs,
                                               rns_params.ps, rns_params.t));
    rns_context_ =
        std::make_unique<const RnsContext<ModularInt>>(std::move(rns_context));
    moduli_ = rns_context_->MainPrimeModuli();
  }

  Integer SmallestPrimeModulus() const {
    Integer q{0};
    for (auto& modulus : moduli_) {
      if (modulus->Modulus() < q) {
        q = modulus->Modulus();
      }
    }
    return q;
  }

  // Generate coefficient vectors containing 0 modulo each prime modulus.
  std::vector<std::vector<ModularInt>> GenerateZeros() const {
    int num_coeffs = 1 << rns_context_->LogN();
    std::vector<std::vector<ModularInt>> coeff_vectors;
    for (auto& modulus : moduli_) {
      std::vector<ModularInt> coeffs(
          num_coeffs, ModularInt::ImportZero(modulus->ModParams()));
      coeff_vectors.push_back(std::move(coeffs));
    }
    return coeff_vectors;
  }

  std::unique_ptr<const RnsContext<ModularInt>> rns_context_;
  std::vector<const PrimeModulus<ModularInt>*> moduli_;
};

TYPED_TEST_SUITE(CrtInterpolationTest, testing::ModularIntTypes);

// Returns a vector of random integers in [0, bound).
template <typename IntegralType>
std::vector<IntegralType> RandomCoeffs(int n, IntegralType bound) {
  absl::BitGen bitgen;
  std::vector<IntegralType> coeffs;
  for (int i = 0; i < n; ++i) {
    IntegralType coeff = absl::Uniform<IntegralType>(bitgen, 0, bound);
    coeffs.push_back(std::move(coeff));
  }
  return coeffs;
}

// Specialized for bound up to 128 bits. Note that this is not meant to be a
// uniform sampler.
template <>
std::vector<Uint128> RandomCoeffs(int n, Uint128 bound) {
  Uint64 bound_hi = absl::Uint128High64(static_cast<absl::uint128>(bound));
  Uint64 bound_lo = absl::Uint128Low64(static_cast<absl::uint128>(bound));
  absl::BitGen bitgen;
  std::vector<Uint128> coeffs;
  for (int i = 0; i < n; ++i) {
    if (bound_hi == 0) {
      Uint128 coeff =
          static_cast<Uint128>(absl::Uniform<Uint64>(bitgen, 0, bound_lo));
      coeffs.push_back(std::move(coeff));
    } else {
      Uint128 coeff =
          static_cast<Uint128>(absl::Uniform<Uint64>(bitgen, 0, bound_hi));
      coeff <<= 64;
      Uint128 extra =
          static_cast<Uint128>(absl::Uniform<Uint64>(bitgen, 0, bound_lo));
      coeff += extra;
      coeffs.push_back(std::move(coeff));
    }
  }
  return coeffs;
}

#ifdef ABSL_HAVE_INTRINSIC_INT128
template <>
std::vector<absl::uint128> RandomCoeffs(int n, absl::uint128 bound) {
  Uint64 bound_hi = absl::Uint128High64(bound);
  Uint64 bound_lo = absl::Uint128Low64(bound);
  absl::BitGen bitgen;
  std::vector<absl::uint128> coeffs;
  for (int i = 0; i < n; ++i) {
    if (bound_hi == 0) {
      absl::uint128 coeff = static_cast<absl::uint128>(
          absl::Uniform<Uint64>(bitgen, 0, bound_lo));
      coeffs.push_back(std::move(coeff));
    } else {
      absl::uint128 coeff = static_cast<absl::uint128>(
          absl::Uniform<Uint64>(bitgen, 0, bound_hi));
      coeff <<= 64;
      absl::uint128 extra = static_cast<absl::uint128>(
          absl::Uniform<Uint64>(bitgen, 0, bound_lo));
      coeff += extra;
      coeffs.push_back(std::move(coeff));
    }
  }
  return coeffs;
}
#endif

// Specialized for bound up to 256 bits. Note that this is not meant to be a
// uniform sampler.
template <>
std::vector<uint256> RandomCoeffs(int n, uint256 bound) {
  absl::BitGen bitgen;
  std::vector<uint256> coeffs;
  for (int i = 0; i < n; ++i) {
    if (Uint256High128(bound) == 0) {
      uint256 coeff = static_cast<uint256>(
          absl::Uniform<absl::uint128>(bitgen, 0, Uint256Low128(bound)));
      coeffs.push_back(std::move(coeff));
    } else {
      uint256 coeff = static_cast<uint256>(
          absl::Uniform<absl::uint128>(bitgen, 0, Uint256High128(bound)));
      coeff <<= 128;
      uint256 extra = static_cast<uint256>(
          absl::Uniform<absl::uint128>(bitgen, 0, Uint256Low128(bound)));
      coeff += extra;
      coeffs.push_back(std::move(coeff));
    }
  }
  return coeffs;
}

TYPED_TEST(CrtInterpolationTest, CrtInterpolationFailsIfCoeffsAreEmpty) {
  using Integer = typename TypeParam::Int;
  using BigInteger = typename testing::CompositeModulus<Integer>::value_type;
  EXPECT_THAT((CrtInterpolation<TypeParam, BigInteger>({}, {}, {}, {})),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`coeff_vectors` must not be empty")));
}

TYPED_TEST(CrtInterpolationTest, CrtInterpolationFailsIfWrongNumberOfModuli) {
  using Integer = typename TypeParam::Int;
  using BigInteger = typename testing::CompositeModulus<Integer>::value_type;
  // Generate vectors containing 0 modulo each prime modulus.
  std::vector<std::vector<TypeParam>> coeff_vectors = this->GenerateZeros();
  std::vector<BigInteger> modulus_hats =
      RnsModulusComplements<TypeParam, BigInteger>(this->moduli_);
  int num_moduli = this->moduli_.size();
  int level = num_moduli - 1;
  ASSERT_OK_AND_ASSIGN(std::vector<TypeParam> modulus_hat_invs,
                       this->rns_context_->MainPrimeModulusCrtFactors(level));
  EXPECT_THAT((CrtInterpolation<TypeParam, BigInteger>(
                  coeff_vectors,
                  absl::MakeSpan(this->moduli_).subspan(1, num_moduli - 1),
                  modulus_hats, modulus_hat_invs)),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr(absl::StrCat("`moduli` must contain ",
                                              num_moduli, " RNS moduli"))));

  EXPECT_THAT((CrtInterpolation<TypeParam, BigInteger>(
                  coeff_vectors, this->moduli_,
                  absl::MakeSpan(modulus_hats).subspan(0, num_moduli - 1),
                  modulus_hat_invs)),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr(absl::StrCat("`modulus_hats` must contain ",
                                              num_moduli, " elements"))));

  EXPECT_THAT(
      (CrtInterpolation<TypeParam, BigInteger>(
          coeff_vectors, this->moduli_, modulus_hats,
          absl::MakeSpan(modulus_hat_invs).subspan(0, num_moduli - 1))),
      StatusIs(absl::StatusCode::kInvalidArgument,
               HasSubstr(absl::StrCat("`modulus_hat_invs` must contain ",
                                      num_moduli, " elements"))));
}

TYPED_TEST(CrtInterpolationTest,
           CrtInterpolationIsCorrectForSmallCoefficients) {
  using Integer = typename TypeParam::Int;
  using BigInteger = typename testing::CompositeModulus<Integer>::value_type;

  int num_coeffs = 1 << this->rns_context_->LogN();
  int num_moduli = this->moduli_.size();

  // 0 <= a < q_i for all prime moduli q_i iff the RNS representation wrt
  // (q_0,..,q_L) is (a,..,a).
  Integer q = this->SmallestPrimeModulus();
  std::vector<Integer> coeffs = RandomCoeffs<Integer>(num_coeffs, q);
  std::vector<std::vector<TypeParam>> coeff_vectors(num_moduli);
  for (int i = 0; i < num_moduli; ++i) {
    for (int j = 0; j < num_coeffs; ++j) {
      ASSERT_OK_AND_ASSIGN(
          TypeParam coeff_qi,
          TypeParam::ImportInt(coeffs[j], this->moduli_[i]->ModParams()));
      coeff_vectors[i].push_back(std::move(coeff_qi));
    }
  }

  std::vector<BigInteger> modulus_hats =
      RnsModulusComplements<TypeParam, BigInteger>(this->moduli_);
  ASSERT_OK_AND_ASSIGN(
      std::vector<TypeParam> modulus_hat_invs,
      this->rns_context_->MainPrimeModulusCrtFactors(num_moduli - 1));
  ASSERT_OK_AND_ASSIGN(
      std::vector<BigInteger> interpolated_coeffs,
      (CrtInterpolation<TypeParam, BigInteger>(
          coeff_vectors, this->moduli_, modulus_hats, modulus_hat_invs)));
  EXPECT_EQ(interpolated_coeffs.size(), num_coeffs);
  for (int i = 0; i < interpolated_coeffs.size(); ++i) {
    EXPECT_EQ(interpolated_coeffs[i], static_cast<BigInteger>(coeffs[i]));
  }
}

TYPED_TEST(CrtInterpolationTest,
           CrtInterpolationIsCorrectForLargeCoefficients) {
  using Integer = typename TypeParam::Int;
  using BigInteger = typename testing::CompositeModulus<Integer>::value_type;

  // Generate a vector of random values in [0, Q) for modulus Q = q_0 * .. * q_l
  int log_n = this->rns_context_->LogN();
  int num_coeffs = 1 << log_n;
  int num_moduli = this->moduli_.size();
  BigInteger big_modulus = testing::CompositeModulus<Integer>::Value();
  std::vector<BigInteger> big_coeffs =
      RandomCoeffs<BigInteger>(num_coeffs, big_modulus);

  // Represent values in `big_coeffs` using RNS system wrt moduli q_0, .., q_l
  std::vector<std::vector<TypeParam>> coeff_vectors(num_moduli);
  for (int i = 0; i < num_moduli; ++i) {
    BigInteger qi = static_cast<BigInteger>(this->moduli_[i]->Modulus());
    for (int j = 0; j < num_coeffs; ++j) {
      Integer coeff_qi_value = static_cast<Integer>(big_coeffs[j] % qi);
      ASSERT_OK_AND_ASSIGN(
          TypeParam coeff_qi,
          TypeParam::ImportInt(coeff_qi_value, this->moduli_[i]->ModParams()));
      coeff_vectors[i].push_back(std::move(coeff_qi));
    }
  }

  // Interpolate the RNS representation, and we should get back values in
  // `big_coeffs` (mod Q).
  std::vector<BigInteger> modulus_hats =
      RnsModulusComplements<TypeParam, BigInteger>(this->moduli_);
  ASSERT_OK_AND_ASSIGN(
      std::vector<TypeParam> modulus_hat_invs,
      this->rns_context_->MainPrimeModulusCrtFactors(num_moduli - 1));
  ASSERT_OK_AND_ASSIGN(
      std::vector<BigInteger> interpolated_coeffs,
      (CrtInterpolation<TypeParam, BigInteger>(
          coeff_vectors, this->moduli_, modulus_hats, modulus_hat_invs)));

  EXPECT_EQ(interpolated_coeffs.size(), num_coeffs);
  for (int j = 0; j < interpolated_coeffs.size(); ++j) {
    EXPECT_EQ(interpolated_coeffs[j], big_coeffs[j]);
  }
}

}  // namespace
}  // namespace rlwe
