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

#include "shell_encryption/rns/rns_polynomial.h"

#include <cmath>
#include <cstddef>
#include <memory>
#include <utility>
#include <vector>

#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include "absl/log/check.h"
#include "absl/random/random.h"
#include "absl/status/status.h"
#include "absl/strings/str_cat.h"
#include "shell_encryption/integral_types.h"
#include "shell_encryption/montgomery.h"
#include "shell_encryption/polynomial.h"
#include "shell_encryption/rns/crt_interpolation.h"
#include "shell_encryption/rns/rns_context.h"
#include "shell_encryption/rns/rns_modulus.h"
#include "shell_encryption/rns/rns_serialization.pb.h"
#include "shell_encryption/rns/testing/parameters.h"
#include "shell_encryption/rns/testing/testing_utils.h"
#include "shell_encryption/status_macros.h"
#include "shell_encryption/testing/parameters.h"
#include "shell_encryption/testing/status_matchers.h"
#include "shell_encryption/testing/status_testing.h"
#include "shell_encryption/testing/testing_prng.h"

namespace rlwe {
namespace {

using Prng = testing::TestingPrng;
using ::rlwe::testing::StatusIs;
using ::testing::HasSubstr;

constexpr Uint64 kMaxCoeff = 1024;  // Upper bound for random coeff values.

// Test fixture for RnsPolynomial. All tests (unless explicitly mentioned) use
// polynomials in Z[X]/(Q, X^N+1), where N = 2^rns_context_.log_n and Q is the
// product of main prime moduli in rns_context_. See testing/parameters.h for
// concrete values used to instantiate these tests.
template <typename ModularInt>
class RnsPolynomialTest : public ::testing::Test {
  using ModularIntParams = typename ModularInt::Params;

 protected:
  void SetUp() override {
    testing::RnsParameters<ModularInt> rns_params =
        testing::GetRnsParameters<ModularInt>();
    auto rns_context = RnsContext<ModularInt>::Create(
        rns_params.log_n, rns_params.qs, rns_params.ps, rns_params.t);
    CHECK(rns_context.ok());
    rns_context_ = std::make_unique<const RnsContext<ModularInt>>(
        std::move(rns_context.value()));
    moduli_ = rns_context_->MainPrimeModuli();

    prng_ = std::make_unique<testing::TestingPrng>(0);
  }

  // Returns a coefficient vector of length 2^log_n, representing values in
  // [0, max_coeff) using mod-q values where q is specified by `mod_params`.
  // The coefficients are not pseudorandom but only for testing operations on
  // polynomials.
  absl::StatusOr<std::vector<ModularInt>> SampleCoeffs(
      int log_n, const ModularIntParams* mod_params,
      Uint64 max_coeff = kMaxCoeff) const {
    absl::BitGen bitgen;
    std::vector<ModularInt> coeffs;
    for (int i = 0; i < (1 << log_n); ++i) {
      Uint64 coeff = absl::Uniform<Uint64>(bitgen, 0, max_coeff);
      RLWE_ASSIGN_OR_RETURN(
          ModularInt coeff_mod,
          ModularInt::ImportInt(static_cast<typename ModularInt::Int>(coeff),
                                mod_params));
      coeffs.push_back(std::move(coeff_mod));
    }
    return coeffs;
  }

  // Returns a coefficient vector of length 2^log_n, representing random
  // values in {-1, 0, 1} using mod-q values.
  absl::StatusOr<std::vector<ModularInt>> SampleTernaryCoeffs(
      int log_n, const ModularIntParams* mod_params) const {
    absl::BitGen bitgen;
    std::vector<ModularInt> coeffs;
    for (int i = 0; i < (1 << log_n); ++i) {
      int8_t x = absl::Uniform<int8_t>(absl::IntervalClosed, bitgen, -1, 1);
      typename ModularInt::Int x_mod_q =
          x == -1 ? mod_params->modulus - 1 : (x == 0 ? 0 : 1);
      RLWE_ASSIGN_OR_RETURN(ModularInt coeff_mod,
                            ModularInt::ImportInt(x_mod_q, mod_params));
      coeffs.push_back(std::move(coeff_mod));
    }
    return coeffs;
  }

  // Returns a RnsPolynomial in Z[X]/(Q, X^N+1) with random coefficients.
  absl::StatusOr<RnsPolynomial<ModularInt>> SampleRnsPolynomial(
      absl::Span<const PrimeModulus<ModularInt>* const> moduli) const {
    int log_n = rns_context_->LogN();
    std::vector<std::vector<ModularInt>> a_coeff_vectors;
    for (auto modulus_qi : moduli) {
      RLWE_ASSIGN_OR_RETURN(std::vector<ModularInt> coeffs_qi,
                            SampleCoeffs(log_n, modulus_qi->ModParams()));
      a_coeff_vectors.push_back(std::move(coeffs_qi));
    }
    return RnsPolynomial<ModularInt>::Create(std::move(a_coeff_vectors),
                                             /*is_ntt=*/true);
  }

  absl::StatusOr<RnsPolynomial<ModularInt>> SampleRnsPolynomial() const {
    return SampleRnsPolynomial(this->moduli_);
  }

  std::unique_ptr<const RnsContext<ModularInt>> rns_context_;
  std::vector<const PrimeModulus<ModularInt>*> moduli_;
  std::unique_ptr<Prng> prng_;
};

TYPED_TEST_SUITE(RnsPolynomialTest, testing::ModularIntTypes);

TYPED_TEST(RnsPolynomialTest, CreateFailsIfCoeffVectorsIsEmpty) {
  EXPECT_THAT(
      RnsPolynomial<TypeParam>::Create(/*coeff_vectors=*/{}, /*is_ntt=*/false),
      StatusIs(absl::StatusCode::kInvalidArgument,
               HasSubstr("`coeff_vectors` cannot be empty")));
}

TYPED_TEST(RnsPolynomialTest, CreateFailsIfNonPowerOfTwoNumberOfCoefficients) {
  std::vector<TypeParam> coeffs(
      3, TypeParam::ImportZero(this->moduli_[0]->ModParams()));
  EXPECT_THAT(RnsPolynomial<TypeParam>::Create(
                  /*coeff_vectors=*/{std::move(coeffs)}, /*is_ntt=*/false),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`coeff_vectors` must contain vectors of "
                                 "length a power of two")));
}

TYPED_TEST(RnsPolynomialTest, CreateFailsIfCoeffVectorsHaveDifferentLengths) {
  if (this->moduli_.size() < 2) {
    // This test requires at least two prime moduli.
    return;
  }
  std::vector<TypeParam> coeffs_q0(
      8, TypeParam::ImportZero(this->moduli_[0]->ModParams()));
  std::vector<TypeParam> coeffs_q1(
      2, TypeParam::ImportZero(this->moduli_[1]->ModParams()));
  EXPECT_THAT(
      RnsPolynomial<TypeParam>::Create(
          /*coeff_vectors=*/{std::move(coeffs_q0), std::move(coeffs_q1)},
          /*is_ntt=*/false),
      StatusIs(
          absl::StatusCode::kInvalidArgument,
          HasSubstr("`coeff_vectors` must contain vectors of equal length")));
}

TYPED_TEST(RnsPolynomialTest, CreateIsCorrect) {
  // Create vectors consisting of 1 wrt each prime modulus.
  int log_n = this->rns_context_->LogN();
  int num_coeffs = 1 << log_n;
  std::vector<std::vector<TypeParam>> coeff_vectors;
  for (auto const& modulus : this->moduli_) {
    std::vector<TypeParam> coeffs(num_coeffs,
                                  TypeParam::ImportOne(modulus->ModParams()));
    coeff_vectors.push_back(std::move(coeffs));
  }
  ASSERT_EQ(coeff_vectors.size(), this->moduli_.size());  // sanity check

  // Create a RnsPolynomial using `coeff_vectors` as its coefficients.
  ASSERT_OK_AND_ASSIGN(auto a, RnsPolynomial<TypeParam>::Create(
                                   coeff_vectors, /*is_ntt=*/false));
  EXPECT_EQ(a.LogN(), log_n);
  EXPECT_EQ(a.NumCoeffs(), num_coeffs);
  EXPECT_EQ(a.NumModuli(), this->moduli_.size());
  EXPECT_EQ(a.Coeffs(), coeff_vectors);
  EXPECT_FALSE(a.IsNttForm());

  // Create another RnsPolynomial using `coeff_vectors` as its NTT coefficients.
  ASSERT_OK_AND_ASSIGN(
      auto b, RnsPolynomial<TypeParam>::Create(coeff_vectors, /*is_ntt=*/true));
  EXPECT_EQ(b.LogN(), log_n);
  EXPECT_EQ(b.NumCoeffs(), num_coeffs);
  EXPECT_EQ(b.NumModuli(), this->moduli_.size());
  EXPECT_EQ(b.Coeffs(), coeff_vectors);
  EXPECT_TRUE(b.IsNttForm());
}

TYPED_TEST(RnsPolynomialTest, CreateEmptyGivesCorrectValues) {
  auto empty = RnsPolynomial<TypeParam>::CreateEmpty();
  EXPECT_EQ(empty.LogN(), 0);
  EXPECT_EQ(empty.NumCoeffs(), 1);
  EXPECT_EQ(empty.NumModuli(), 0);
  EXPECT_TRUE(empty.Coeffs().empty());
}

TYPED_TEST(RnsPolynomialTest, CreateZeroFailsIfLogNIsNotPositive) {
  EXPECT_THAT(RnsPolynomial<TypeParam>::CreateZero(/*log_n=*/-1, {}),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`log_n` must be positive")));
  EXPECT_THAT(RnsPolynomial<TypeParam>::CreateZero(/*log_n=*/0, {}),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`log_n` must be positive")));
}

TYPED_TEST(RnsPolynomialTest, CreateZeroFailsIfModulliIsEmpty) {
  EXPECT_THAT(RnsPolynomial<TypeParam>::CreateZero(this->rns_context_->LogN(),
                                                   /*moduli=*/{}),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`moduli` must not be empty")));
}

TYPED_TEST(RnsPolynomialTest, CreateOneFailsIfLogNIsNotPositive) {
  for (int log_n = -1; log_n <= 0; ++log_n) {
    EXPECT_THAT(RnsPolynomial<TypeParam>::CreateOne(log_n, {}),
                StatusIs(absl::StatusCode::kInvalidArgument,
                         HasSubstr("`log_n` must be positive")));
  }
}

TYPED_TEST(RnsPolynomialTest, CreateOneFailsIfModulliIsEmpty) {
  EXPECT_THAT(RnsPolynomial<TypeParam>::CreateOne(this->rns_context_->LogN(),
                                                  /*moduli=*/{}),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`moduli` must not be empty")));
}

TYPED_TEST(RnsPolynomialTest, CreateZeroIsCorrect) {
  int log_n = this->rns_context_->LogN();
  int num_coeffs = 1 << log_n;
  ASSERT_OK_AND_ASSIGN(
      auto zero, RnsPolynomial<TypeParam>::CreateZero(log_n, this->moduli_));
  EXPECT_EQ(zero.LogN(), log_n);
  EXPECT_EQ(zero.NumCoeffs(), num_coeffs);
  EXPECT_EQ(zero.NumModuli(), this->moduli_.size());
  EXPECT_TRUE(zero.IsNttForm());

  // The coefficients of the polynomial 0 should be all 0.
  const std::vector<std::vector<TypeParam>>& coeff_vectors = zero.Coeffs();
  ASSERT_EQ(coeff_vectors.size(), this->moduli_.size());
  for (int i = 0; i < this->moduli_.size(); ++i) {
    ASSERT_EQ(coeff_vectors[i].size(), num_coeffs);
    for (int j = 0; j < num_coeffs; ++j) {
      EXPECT_EQ(coeff_vectors[i][j],
                TypeParam::ImportZero(this->moduli_[i]->ModParams()));
    }
  }
}

TYPED_TEST(RnsPolynomialTest, CreateOneIsCorrect) {
  int log_n = this->rns_context_->LogN();
  int num_coeffs = 1 << log_n;
  ASSERT_OK_AND_ASSIGN(
      auto one, RnsPolynomial<TypeParam>::CreateOne(log_n, this->moduli_));
  EXPECT_EQ(one.LogN(), log_n);
  EXPECT_EQ(one.NumCoeffs(), num_coeffs);
  EXPECT_EQ(one.NumModuli(), this->moduli_.size());
  EXPECT_TRUE(one.IsNttForm());

  // The coefficients of the polynomial 1 should be all 1.
  const std::vector<std::vector<TypeParam>>& coeff_vectors = one.Coeffs();
  ASSERT_EQ(coeff_vectors.size(), this->moduli_.size());
  for (int i = 0; i < this->moduli_.size(); ++i) {
    ASSERT_EQ(coeff_vectors[i].size(), num_coeffs);
    for (int j = 0; j < num_coeffs; ++j) {
      EXPECT_EQ(coeff_vectors[i][j],
                TypeParam::ImportOne(this->moduli_[i]->ModParams()));
    }
  }
}

TYPED_TEST(RnsPolynomialTest, SampleUniformFailsIfLogNIsNotPositive) {
  EXPECT_THAT(RnsPolynomial<TypeParam>::template SampleUniform<Prng>(
                  /*log_n=*/-1, this->prng_.get(), this->moduli_),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`log_n` must be positive")));
  EXPECT_THAT(RnsPolynomial<TypeParam>::template SampleUniform<Prng>(
                  /*log_n=*/0, this->prng_.get(), this->moduli_),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`log_n` must be positive")));
}

TYPED_TEST(RnsPolynomialTest, SampleUniformFailsIfPrngIsNull) {
  EXPECT_THAT(RnsPolynomial<TypeParam>::template SampleUniform<Prng>(
                  this->rns_context_->LogN(), /*prng=*/nullptr, this->moduli_),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`prng` must not be null")));
}

TYPED_TEST(RnsPolynomialTest, SampleUniformFailsIfModuliIsEmpty) {
  EXPECT_THAT(RnsPolynomial<TypeParam>::template SampleUniform<Prng>(
                  this->rns_context_->LogN(), this->prng_.get(),
                  /*moduli=*/{}),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`moduli` must not be empty")));
}

TYPED_TEST(RnsPolynomialTest, SampleUniformSucceeds) {
  int log_n = this->rns_context_->LogN();
  ASSERT_OK_AND_ASSIGN(RnsPolynomial<TypeParam> polynomial,
                       RnsPolynomial<TypeParam>::SampleUniform(
                           log_n, this->prng_.get(), this->moduli_));
  EXPECT_TRUE(polynomial.IsNttForm());
  EXPECT_EQ(polynomial.LogN(), log_n);
  EXPECT_EQ(polynomial.NumCoeffs(), (1 << log_n));
  EXPECT_EQ(polynomial.NumModuli(), this->moduli_.size());
}

TYPED_TEST(RnsPolynomialTest, ConvertFromPolynomialCoeffsFailsIfNotPowerOf2) {
  constexpr int k_num_coeffs = 3;  // not a power of 2
  const typename TypeParam::Params* mp = this->moduli_[0]->ModParams();
  std::vector<TypeParam> coeffs_q(k_num_coeffs, TypeParam::ImportZero(mp));
  EXPECT_THAT(RnsPolynomial<TypeParam>::ConvertFromPolynomialCoeffs(
                  coeffs_q, mp, this->moduli_),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`coeffs_q` must have length a power of 2")));
}

TYPED_TEST(RnsPolynomialTest,
           ConvertFromPolynomialCoeffsFailsIfModParamsIsNull) {
  int log_n = this->rns_context_->LogN();
  const typename TypeParam::Params* mp = this->moduli_[0]->ModParams();
  std::vector<TypeParam> coeffs_q(1 << log_n, TypeParam::ImportZero(mp));
  EXPECT_THAT(RnsPolynomial<TypeParam>::ConvertFromPolynomialCoeffs(
                  coeffs_q,
                  /*mod_params_q=*/nullptr, this->moduli_),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`mod_params_q` must not be null")));
}

TYPED_TEST(RnsPolynomialTest, ConvertFromPolynomialCoeffsFailsIfModuliIsEmpty) {
  int log_n = this->rns_context_->LogN();
  const typename TypeParam::Params* mp = this->moduli_[0]->ModParams();
  std::vector<TypeParam> coeffs_q(1 << log_n, TypeParam::ImportZero(mp));
  EXPECT_THAT(RnsPolynomial<TypeParam>::ConvertFromPolynomialCoeffs(
                  coeffs_q, mp, /*moduli=*/{}),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`moduli` must not be empty")));
}

TYPED_TEST(RnsPolynomialTest, ConvertFromPolynomialCoeffsSucceeds) {
  int log_n = this->rns_context_->LogN();
  auto mod_params_q0 = this->moduli_[0]->ModParams();
  ASSERT_OK_AND_ASSIGN(std::vector<TypeParam> coeffs_q0,
                       this->SampleCoeffs(log_n, mod_params_q0));
  ASSERT_OK_AND_ASSIGN(RnsPolynomial<TypeParam> polynomial,
                       RnsPolynomial<TypeParam>::ConvertFromPolynomialCoeffs(
                           coeffs_q0, mod_params_q0, this->moduli_));
  EXPECT_TRUE(polynomial.IsNttForm());
  EXPECT_EQ(polynomial.NumModuli(), this->moduli_.size());
  EXPECT_EQ(polynomial.NumCoeffs(), 1 << log_n);
  EXPECT_EQ(polynomial.LogN(), log_n);

  // Coefficients of `polynomial` should be RNS representation of values
  // in `coeffs_q0`.
  ASSERT_OK(polynomial.ConvertToCoeffForm(this->moduli_));
  for (int i = 0; i < this->moduli_.size(); ++i) {
    auto mod_params_qi = this->moduli_[i]->ModParams();
    auto coeffs_qi = polynomial.Coeffs()[i];
    ASSERT_EQ(coeffs_qi.size(), coeffs_q0.size());
    for (int j = 0; j < coeffs_q0.size(); ++j) {
      if (mod_params_q0->modulus < mod_params_qi->modulus) {
        auto x = coeffs_qi[j].ExportInt(mod_params_qi);
        ASSERT_OK_AND_ASSIGN(auto x_mod_q0,
                             TypeParam::ImportInt(x, mod_params_q0));
        auto y_mod_q0 = coeffs_q0[j];
        EXPECT_EQ(x_mod_q0, y_mod_q0);
      } else {
        auto x_mod_qi = coeffs_qi[j];
        auto y = coeffs_q0[j].ExportInt(mod_params_q0);
        ASSERT_OK_AND_ASSIGN(auto y_mod_qi,
                             TypeParam::ImportInt(y, mod_params_qi));
        EXPECT_EQ(x_mod_qi, y_mod_qi);
      }
    }
  }
}

TYPED_TEST(RnsPolynomialTest,
           ConvertBalancedFromPolynomialCoeffsFailsIfNotPowerOf2) {
  constexpr int k_num_coeffs = 3;  // not a power of 2
  const typename TypeParam::Params* mp = this->moduli_[0]->ModParams();
  std::vector<TypeParam> coeffs_q(k_num_coeffs, TypeParam::ImportZero(mp));
  EXPECT_THAT(RnsPolynomial<TypeParam>::ConvertBalancedFromPolynomialCoeffs(
                  coeffs_q, mp, this->moduli_),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`coeffs_q` must have length a power of 2")));
}

TYPED_TEST(RnsPolynomialTest,
           ConvertBalancedFromPolynomialCoeffsFailsIfModParamsIsNull) {
  int log_n = this->rns_context_->LogN();
  const typename TypeParam::Params* mp = this->moduli_[0]->ModParams();
  std::vector<TypeParam> coeffs_q(1 << log_n, TypeParam::ImportZero(mp));
  EXPECT_THAT(RnsPolynomial<TypeParam>::ConvertBalancedFromPolynomialCoeffs(
                  coeffs_q,
                  /*mod_params_q=*/nullptr, this->moduli_),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`mod_params_q` must not be null")));
}

TYPED_TEST(RnsPolynomialTest,
           ConvertBalancedFromPolynomialCoeffsFailsIfModuliIsEmpty) {
  int log_n = this->rns_context_->LogN();
  const typename TypeParam::Params* mp = this->moduli_[0]->ModParams();
  std::vector<TypeParam> coeffs_q(1 << log_n, TypeParam::ImportZero(mp));
  EXPECT_THAT(RnsPolynomial<TypeParam>::ConvertBalancedFromPolynomialCoeffs(
                  coeffs_q, mp, /*moduli=*/{}),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`moduli` must not be empty")));
}

TYPED_TEST(RnsPolynomialTest, ConvertBalancedFromPolynomialSucceeds) {
  int log_n = this->rns_context_->LogN();
  // Find the smallest prime modulus q to sample a polynomial with balanced
  // coefficients in [-q/2, q/2].
  const typename TypeParam::Params* mod_params_q0 =
      this->moduli_[0]->ModParams();
  for (auto& modulus : this->moduli_) {
    if (mod_params_q0->modulus > modulus->Modulus()) {
      mod_params_q0 = modulus->ModParams();
    }
  }

  ASSERT_OK_AND_ASSIGN(std::vector<TypeParam> coeffs_q0,
                       this->SampleCoeffs(log_n, mod_params_q0));
  ASSERT_OK_AND_ASSIGN(
      RnsPolynomial<TypeParam> polynomial,
      RnsPolynomial<TypeParam>::ConvertBalancedFromPolynomialCoeffs(
          coeffs_q0, mod_params_q0, this->moduli_));
  EXPECT_TRUE(polynomial.IsNttForm());
  EXPECT_EQ(polynomial.NumModuli(), this->moduli_.size());
  EXPECT_EQ(polynomial.NumCoeffs(), 1 << log_n);
  EXPECT_EQ(polynomial.LogN(), log_n);

  // Coefficients of `polynomial` should agree with `coeffs_q0` in representing
  // balanced integer values in [-q/2, q/2].
  typename TypeParam::Int q0 = mod_params_q0->modulus;
  typename TypeParam::Int q0_half = q0 >> 1;
  ASSERT_OK(polynomial.ConvertToCoeffForm(this->moduli_));
  for (int i = 0; i < this->moduli_.size(); ++i) {
    auto mod_params_qi = this->moduli_[i]->ModParams();
    typename TypeParam::Int qi = mod_params_qi->modulus;
    typename TypeParam::Int qi_half = qi >> 1;
    auto coeffs_qi = polynomial.Coeffs()[i];
    ASSERT_EQ(coeffs_qi.size(), coeffs_q0.size());
    for (int j = 0; j < coeffs_q0.size(); ++j) {
      auto x = coeffs_qi[j].ExportInt(mod_params_qi);
      auto y = coeffs_q0[j].ExportInt(mod_params_q0);
      bool is_x_neg = x > qi_half;
      bool is_y_neg = y > q0_half;
      typename TypeParam::Int x_abs = is_x_neg ? qi - x : x;
      typename TypeParam::Int y_abs = is_y_neg ? q0 - y : y;
      EXPECT_EQ(is_x_neg, is_y_neg);
      EXPECT_EQ(x_abs, y_abs);
    }
  }
}

TYPED_TEST(RnsPolynomialTest, ConvertToCoeffFormFailsIfAlreadyInCoeffForm) {
  // Create a RNS polynomial in coefficient form.
  ASSERT_OK_AND_ASSIGN(auto a, RnsPolynomial<TypeParam>::CreateZero(
                                   this->rns_context_->LogN(), this->moduli_,
                                   /*is_ntt=*/false));
  EXPECT_THAT(a.ConvertToCoeffForm(this->moduli_),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("Polynomial already in Coefficient form")));
}

TYPED_TEST(RnsPolynomialTest, ConvertToCoeffFormFailsIfWrongNumberOfModuli) {
  int num_moduli = this->moduli_.size();
  ASSERT_GE(num_moduli, 1);
  ASSERT_OK_AND_ASSIGN(auto a, RnsPolynomial<TypeParam>::CreateZero(
                                   this->rns_context_->LogN(), this->moduli_,
                                   /*is_ntt=*/true));
  EXPECT_THAT(a.ConvertToCoeffForm(
                  absl::MakeSpan(this->moduli_).subspan(0, num_moduli - 1)),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr(absl::StrCat("`moduli` must contain ",
                                              num_moduli, " RNS moduli"))));
}

TYPED_TEST(RnsPolynomialTest, ConvertToCoeffFormSuccesds) {
  // RNS polynomial 1 in Ntt form.
  ASSERT_OK_AND_ASSIGN(auto a, RnsPolynomial<TypeParam>::CreateOne(
                                   this->rns_context_->LogN(), this->moduli_));
  ASSERT_OK(a.ConvertToCoeffForm(this->moduli_));
  EXPECT_FALSE(a.IsNttForm());
  // The coefficient form of the polynomial 1 should have coefficient
  // vectors [1,0,..0]
  int expected_num_moduli = this->moduli_.size();
  int expected_num_coeffs = 1 << this->rns_context_->LogN();
  ASSERT_EQ(a.NumModuli(), expected_num_moduli);
  ASSERT_EQ(a.NumCoeffs(), expected_num_coeffs);
  for (int i = 0; i < expected_num_moduli; ++i) {
    auto mod_params_qi = this->moduli_[i]->ModParams();
    ASSERT_EQ(a.Coeffs()[i].size(), expected_num_coeffs);
    EXPECT_EQ(a.Coeffs()[i][0], TypeParam::ImportOne(mod_params_qi));
    for (int j = 1; j < expected_num_coeffs; ++j) {
      EXPECT_EQ(a.Coeffs()[i][j], TypeParam::ImportZero(mod_params_qi));
    }
  }
}

TYPED_TEST(RnsPolynomialTest, ConvertToNttFormFailsIfAlreadyInNttForm) {
  // Create a RNS polynomial in NTT form.
  ASSERT_OK_AND_ASSIGN(auto a, RnsPolynomial<TypeParam>::CreateZero(
                                   this->rns_context_->LogN(), this->moduli_,
                                   /*is_ntt=*/true));
  EXPECT_THAT(a.ConvertToNttForm(this->moduli_),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("Polynomial already in NTT form")));
}

TYPED_TEST(RnsPolynomialTest, ConvertToNttFormFailsIfWrongNumberOfModuli) {
  int num_moduli = this->moduli_.size();
  ASSERT_GE(num_moduli, 1);
  ASSERT_OK_AND_ASSIGN(auto a, RnsPolynomial<TypeParam>::CreateZero(
                                   this->rns_context_->LogN(), this->moduli_,
                                   /*is_ntt=*/false));
  EXPECT_THAT(a.ConvertToNttForm(
                  absl::MakeSpan(this->moduli_).subspan(0, num_moduli - 1)),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr(absl::StrCat("`moduli` must contain ",
                                              num_moduli, " RNS moduli"))));
}

TYPED_TEST(RnsPolynomialTest, ConvertToNttFormSuccesds) {
  // Random polynomial in NTT form.
  ASSERT_OK_AND_ASSIGN(RnsPolynomial<TypeParam> a, this->SampleRnsPolynomial());
  ASSERT_TRUE(a.IsNttForm());

  // Coefficient form of the same polynomial.
  RnsPolynomial<TypeParam> b = a;
  ASSERT_OK(b.ConvertToCoeffForm(this->moduli_));
  ASSERT_FALSE(b.IsNttForm());

  // Convert back to NTT form and we should have the same representation.
  ASSERT_OK(b.ConvertToNttForm(this->moduli_));
  ASSERT_TRUE(b.IsNttForm());
  EXPECT_EQ(a, b);
}

////////////////////////////////////////////////////////////
// Additions
////////////////////////////////////////////////////////////

TYPED_TEST(RnsPolynomialTest, NegateFailsIfWrongNumberOfModuli) {
  ASSERT_OK_AND_ASSIGN(auto a, RnsPolynomial<TypeParam>::CreateOne(
                                   this->rns_context_->LogN(), this->moduli_));
  int num_moduli = this->moduli_.size();
  ASSERT_GE(num_moduli, 1);
  EXPECT_THAT(
      a.NegateInPlace(absl::MakeSpan(this->moduli_).subspan(0, num_moduli - 1)),
      StatusIs(absl::StatusCode::kInvalidArgument,
               HasSubstr(absl::StrCat("`moduli` must contain ", num_moduli,
                                      " RNS moduli"))));
}

TYPED_TEST(RnsPolynomialTest, NegateIsCorrect) {
  ASSERT_OK_AND_ASSIGN(RnsPolynomial<TypeParam> a, this->SampleRnsPolynomial());
  ASSERT_OK_AND_ASSIGN(RnsPolynomial<TypeParam> b, a.Negate(this->moduli_));
  // a + (-a) should be 0.
  ASSERT_OK_AND_ASSIGN(RnsPolynomial<TypeParam> sum, a.Add(b, this->moduli_));
  ASSERT_OK_AND_ASSIGN(auto zero,
                       RnsPolynomial<TypeParam>::CreateZero(
                           this->rns_context_->LogN(), this->moduli_));
  EXPECT_EQ(sum, zero);
}

TYPED_TEST(RnsPolynomialTest, AddFailsIfThatContainsDifferentNumberOfVectors) {
  int num_moduli = this->moduli_.size();
  if (num_moduli <= 1) {
    return;  // This test requires at least two prime moduli.
  }

  // Create two polynomials with respect to different number of prime moduli.
  ASSERT_OK_AND_ASSIGN(auto a, RnsPolynomial<TypeParam>::CreateOne(
                                   this->rns_context_->LogN(), this->moduli_));
  ASSERT_OK_AND_ASSIGN(
      auto b, RnsPolynomial<TypeParam>::CreateOne(
                  this->rns_context_->LogN(),
                  absl::MakeSpan(this->moduli_).subspan(0, num_moduli - 1)));
  EXPECT_THAT(
      a.AddInPlace(b, this->moduli_),
      StatusIs(absl::StatusCode::kInvalidArgument,
               HasSubstr(absl::StrCat("`that` must contain ", num_moduli,
                                      " coefficient vectors"))));
}

TYPED_TEST(RnsPolynomialTest, AddFailsIfWrongNumberOfModuli) {
  ASSERT_OK_AND_ASSIGN(auto a, RnsPolynomial<TypeParam>::CreateOne(
                                   this->rns_context_->LogN(), this->moduli_));
  int num_moduli = this->moduli_.size();
  ASSERT_GE(num_moduli, 1);
  EXPECT_THAT(
      a.AddInPlace(a, absl::MakeSpan(this->moduli_).subspan(0, num_moduli - 1)),
      StatusIs(absl::StatusCode::kInvalidArgument,
               HasSubstr(absl::StrCat("`moduli` must contain ", num_moduli,
                                      " RNS moduli"))));
}

TYPED_TEST(RnsPolynomialTest, AddIsCorrect) {
  ASSERT_OK_AND_ASSIGN(auto zero,
                       RnsPolynomial<TypeParam>::CreateZero(
                           this->rns_context_->LogN(), this->moduli_));
  ASSERT_OK_AND_ASSIGN(RnsPolynomial<TypeParam> a, this->SampleRnsPolynomial());
  ASSERT_OK_AND_ASSIGN(RnsPolynomial<TypeParam> b, this->SampleRnsPolynomial());

  // a + 0 == a
  ASSERT_OK_AND_ASSIGN(RnsPolynomial<TypeParam> sum0,
                       a.Add(zero, this->moduli_));
  EXPECT_EQ(sum0, a);

  // a + b == b + a
  ASSERT_OK_AND_ASSIGN(RnsPolynomial<TypeParam> sum1, a.Add(b, this->moduli_));
  ASSERT_OK_AND_ASSIGN(RnsPolynomial<TypeParam> sum2, b.Add(a, this->moduli_));
  EXPECT_EQ(sum1, sum2);
}

TYPED_TEST(RnsPolynomialTest, SubFailsIfThatContainsDifferentNumberOfVectors) {
  int num_moduli = this->moduli_.size();
  if (num_moduli <= 1) {
    return;  // This test requires at least two prime moduli.
  }

  // Create two polynomials with respect to different number of prime moduli.
  ASSERT_OK_AND_ASSIGN(auto a, RnsPolynomial<TypeParam>::CreateOne(
                                   this->rns_context_->LogN(), this->moduli_));
  ASSERT_OK_AND_ASSIGN(
      auto b, RnsPolynomial<TypeParam>::CreateOne(
                  this->rns_context_->LogN(),
                  absl::MakeSpan(this->moduli_).subspan(0, num_moduli - 1)));
  EXPECT_THAT(
      a.SubInPlace(b, this->moduli_),
      StatusIs(absl::StatusCode::kInvalidArgument,
               HasSubstr(absl::StrCat("`that` must contain ", num_moduli,
                                      " coefficient vectors"))));
}

TYPED_TEST(RnsPolynomialTest, SubFailsIfWrongNumberOfModuli) {
  ASSERT_OK_AND_ASSIGN(auto a, RnsPolynomial<TypeParam>::CreateOne(
                                   this->rns_context_->LogN(), this->moduli_));
  int num_moduli = this->moduli_.size();
  ASSERT_GE(num_moduli, 1);
  EXPECT_THAT(
      a.SubInPlace(a, absl::MakeSpan(this->moduli_).subspan(0, num_moduli - 1)),
      StatusIs(absl::StatusCode::kInvalidArgument,
               HasSubstr(absl::StrCat("`moduli` must contain ", num_moduli,
                                      " RNS moduli"))));
}

TYPED_TEST(RnsPolynomialTest, SubIsCorrect) {
  ASSERT_OK_AND_ASSIGN(auto zero,
                       RnsPolynomial<TypeParam>::CreateZero(
                           this->rns_context_->LogN(), this->moduli_));
  ASSERT_OK_AND_ASSIGN(RnsPolynomial<TypeParam> a, this->SampleRnsPolynomial());
  ASSERT_OK_AND_ASSIGN(RnsPolynomial<TypeParam> a_neg, a.Negate(this->moduli_));
  // Negation is performed without Barrett reduction, so it's possible to end up
  // with two MontgomeryInt with different representations (e.g. 0 vs modulus)
  // but for the same modular integer value (in this case, 0 (mod modulus)).
  // This might incorrectly fail the equality check which compares the
  // Montgomery representation for performance consideration.
  // So we add zero to a_neg to force Barrett reduction on it.
  ASSERT_OK(a_neg.AddInPlace(zero, this->moduli_));

  // a - 0 = a
  ASSERT_OK_AND_ASSIGN(RnsPolynomial<TypeParam> diff0,
                       a.Sub(zero, this->moduli_));
  EXPECT_EQ(diff0, a);

  // 0 - a = -a
  ASSERT_OK_AND_ASSIGN(RnsPolynomial<TypeParam> diff1,
                       zero.Sub(a, this->moduli_));
  EXPECT_EQ(diff1, a_neg);

  // a - b = -(b - a);
  ASSERT_OK_AND_ASSIGN(RnsPolynomial<TypeParam> b, this->SampleRnsPolynomial());
  ASSERT_OK_AND_ASSIGN(RnsPolynomial<TypeParam> diff2, a.Sub(b, this->moduli_));
  ASSERT_OK_AND_ASSIGN(RnsPolynomial<TypeParam> diff3, b.Sub(a, this->moduli_));
  ASSERT_OK_AND_ASSIGN(RnsPolynomial<TypeParam> diff3_neg,
                       diff3.Negate(this->moduli_));
  // Same as above, force Barrett reduction on the negated values.
  ASSERT_OK(diff3_neg.AddInPlace(zero, this->moduli_));
  EXPECT_EQ(diff2, diff3_neg);
}

////////////////////////////////////////////////////////////
// Multiplications
////////////////////////////////////////////////////////////

// Negative test for MulInPlace(int_value, moduli) when moduli does not contain
// the desired number of modulus parameters.
TYPED_TEST(RnsPolynomialTest,
           MulInPlaceWithIntScalarFailsIfWrongNumberOfModuli) {
  ASSERT_OK_AND_ASSIGN(auto a, RnsPolynomial<TypeParam>::CreateOne(
                                   this->rns_context_->LogN(), this->moduli_));
  int num_moduli = this->moduli_.size();
  ASSERT_GE(num_moduli, 1);
  EXPECT_THAT(
      a.MulInPlace(0, absl::MakeSpan(this->moduli_).subspan(0, num_moduli - 1)),
      StatusIs(absl::StatusCode::kInvalidArgument,
               HasSubstr(absl::StrCat("`moduli` must contain ", num_moduli,
                                      " RNS moduli"))));
}

TYPED_TEST(RnsPolynomialTest, MulInPlaceWithIntScalarIsCorrect) {
  int num_coeffs = 1 << this->rns_context_->LogN();
  ASSERT_OK_AND_ASSIGN(RnsPolynomial<TypeParam> a, this->SampleRnsPolynomial());

  // a * 1 == a
  ASSERT_OK_AND_ASSIGN(RnsPolynomial<TypeParam> b, a.Mul(1, this->moduli_));
  ASSERT_EQ(b.NumCoeffs(), num_coeffs);
  ASSERT_EQ(b.NumModuli(), this->moduli_.size());
  for (int i = 0; i < this->moduli_.size(); ++i) {
    EXPECT_EQ(b.Coeffs()[i], a.Coeffs()[i]);
  }

  // a * 0 == 0
  ASSERT_OK_AND_ASSIGN(RnsPolynomial<TypeParam> c, a.Mul(0, this->moduli_));
  ASSERT_EQ(c.NumCoeffs(), num_coeffs);
  ASSERT_EQ(c.NumModuli(), this->moduli_.size());
  for (int i = 0; i < this->moduli_.size(); ++i) {
    TypeParam zero = TypeParam::ImportZero(this->moduli_[i]->ModParams());
    for (int j = 0; j < num_coeffs; ++j) {
      EXPECT_EQ(c.Coeffs()[i][j], zero);
    }
  }

  // a * 2 == a + a
  ASSERT_OK_AND_ASSIGN(RnsPolynomial<TypeParam> d, a.Mul(2, this->moduli_));
  ASSERT_EQ(d.NumCoeffs(), num_coeffs);
  ASSERT_EQ(d.NumModuli(), this->moduli_.size());
  for (int i = 0; i < this->moduli_.size(); ++i) {
    auto mod_params_qi = this->moduli_[i]->ModParams();
    for (int j = 0; j < num_coeffs; j++) {
      TypeParam expected = a.Coeffs()[i][j];
      expected.AddInPlace(expected, mod_params_qi);
      EXPECT_EQ(d.Coeffs()[i][j], expected);
    }
  }

  // Multiply a with prime q_i should make it vanish (mod q_i).
  for (int i = 0; i < this->moduli_.size(); ++i) {
    auto mod_params_qi = this->moduli_[i]->ModParams();
    ASSERT_OK_AND_ASSIGN(RnsPolynomial<TypeParam> e,
                         a.Mul(mod_params_qi->modulus, this->moduli_));
    TypeParam zero = TypeParam::ImportZero(mod_params_qi);
    for (int j = 0; j < num_coeffs; ++j) {
      EXPECT_EQ(e.Coeffs()[i][j], zero);
    }
  }
}

// Mul() should be equivalent to MulInPlace().
TYPED_TEST(RnsPolynomialTest, MulWithIntScalarIsCorrect) {
  ASSERT_OK_AND_ASSIGN(RnsPolynomial<TypeParam> a, this->SampleRnsPolynomial());
  for (int scalar_value = 0; scalar_value < 3; scalar_value++) {
    auto scalar = static_cast<typename TypeParam::Int>(scalar_value);
    RnsPolynomial<TypeParam> b = a;
    ASSERT_OK(b.MulInPlace(scalar, this->moduli_));
    ASSERT_OK_AND_ASSIGN(RnsPolynomial<TypeParam> c,
                         a.Mul(scalar, this->moduli_));
    EXPECT_EQ(b, c);
  }
}

// Negative test for MulInPlace(rns_int, moduli) when `rns_int` does not contain
// the correct number of modular integers.
TYPED_TEST(RnsPolynomialTest,
           MulInPlaceWithRnsScalarFailsIfWrongNumberOfModularIntegers) {
  ASSERT_OK_AND_ASSIGN(auto a, RnsPolynomial<TypeParam>::CreateOne(
                                   this->rns_context_->LogN(), this->moduli_));
  int num_moduli = this->moduli_.size();
  ASSERT_GE(num_moduli, 1);
  std::vector<TypeParam> one_mod_qs;
  std::transform(this->moduli_.begin(), this->moduli_.end(),
                 std::back_inserter(one_mod_qs), [](auto& moduli_it) {
                   return TypeParam::ImportOne(moduli_it->ModParams());
                 });
  one_mod_qs.pop_back();  // remove the last modular integer
  RnsInt<TypeParam> one_missing_last_modulus{std::move(one_mod_qs)};
  EXPECT_THAT(
      a.MulInPlace(one_missing_last_modulus, this->moduli_),
      StatusIs(absl::StatusCode::kInvalidArgument,
               HasSubstr(absl::StrCat("`scalar_mod_qs` must contain ",
                                      num_moduli, " modular integers"))));
}

// Negative test for MulInPlace(rns_int, moduli) when `moduli` does not contain
// the correct number of modulus parameters.
TYPED_TEST(RnsPolynomialTest,
           MulInPlaceWithRnsScalarFailsIfWrongNumberOfModuli) {
  ASSERT_OK_AND_ASSIGN(auto a, RnsPolynomial<TypeParam>::CreateOne(
                                   this->rns_context_->LogN(), this->moduli_));
  int num_moduli = this->moduli_.size();
  ASSERT_GE(num_moduli, 1);
  std::vector<TypeParam> one_mod_qs;
  std::transform(this->moduli_.begin(), this->moduli_.end(),
                 std::back_inserter(one_mod_qs), [](auto& moduli_it) {
                   return TypeParam::ImportOne(moduli_it->ModParams());
                 });
  RnsInt<TypeParam> one{std::move(one_mod_qs)};
  EXPECT_THAT(
      a.MulInPlace(one,
                   absl::MakeSpan(this->moduli_).subspan(0, num_moduli - 1)),
      StatusIs(absl::StatusCode::kInvalidArgument,
               HasSubstr(absl::StrCat("`moduli` must contain ", num_moduli,
                                      " RNS moduli"))));
}

TYPED_TEST(RnsPolynomialTest, MulInPlaceWithRnsScalarIsCorrect) {
  ASSERT_OK_AND_ASSIGN(RnsPolynomial<TypeParam> polynomial0,
                       this->SampleRnsPolynomial());

  // Multiply with a prime modulus q_i in RNS form, the result should be
  // the same as multiplying the polynomial with q_i as an integer.
  for (int i = 0; i < this->moduli_.size(); ++i) {
    typename TypeParam::Int qi = this->moduli_[i]->Modulus();
    std::vector<TypeParam> qi_mod_qs;
    std::transform(
        this->moduli_.begin(), this->moduli_.end(),
        std::back_inserter(qi_mod_qs), [&qi](auto& moduli_it) {
          return TypeParam::ImportInt(qi, moduli_it->ModParams()).value();
        });
    RnsInt<TypeParam> qi_rns{std::move(qi_mod_qs)};
    RnsPolynomial<TypeParam> polynomial1 = polynomial0;
    ASSERT_OK(polynomial1.MulInPlace(qi_rns, this->moduli_));

    ASSERT_OK_AND_ASSIGN(RnsPolynomial<TypeParam> expected,
                         polynomial0.Mul(qi, this->moduli_));
    EXPECT_EQ(polynomial1, expected);
  }
}

TYPED_TEST(RnsPolynomialTest,
           MulInPlaceWithRnsPolynomialFailsIfWrongNumberOfCoeffVectors) {
  int num_moduli = this->moduli_.size();
  if (num_moduli <= 1) {
    return;  // This test requires at least two prime moduli.
  }

  int log_n = this->rns_context_->LogN();
  ASSERT_OK_AND_ASSIGN(
      RnsPolynomial<TypeParam> a,
      RnsPolynomial<TypeParam>::CreateZero(log_n, this->moduli_));
  ASSERT_OK_AND_ASSIGN(
      RnsPolynomial<TypeParam> b_missing_last_modulus,
      RnsPolynomial<TypeParam>::CreateZero(
          log_n, absl::MakeSpan(this->moduli_).subspan(0, num_moduli - 1)));
  EXPECT_THAT(
      a.MulInPlace(b_missing_last_modulus, this->moduli_),
      StatusIs(absl::StatusCode::kInvalidArgument,
               HasSubstr(absl::StrCat("`that` must contain ", num_moduli,
                                      " coefficient vectors"))));
}

TYPED_TEST(RnsPolynomialTest,
           MulInPlaceWithRnsPolynomialFailsIfWrongNumberOfModuli) {
  int num_moduli = this->moduli_.size();
  if (num_moduli <= 1) {
    return;  // This test requires at least two prime moduli.
  }

  int log_n = this->rns_context_->LogN();
  ASSERT_OK_AND_ASSIGN(
      auto a, RnsPolynomial<TypeParam>::CreateZero(log_n, this->moduli_));
  ASSERT_OK_AND_ASSIGN(
      auto b, RnsPolynomial<TypeParam>::CreateZero(log_n, this->moduli_));
  EXPECT_THAT(
      a.MulInPlace(b, absl::MakeSpan(this->moduli_).subspan(0, num_moduli - 1)),
      StatusIs(absl::StatusCode::kInvalidArgument,
               HasSubstr(absl::StrCat("`moduli` must contain ", num_moduli,
                                      " RNS moduli"))));
}

TYPED_TEST(RnsPolynomialTest,
           MulInPlaceWithRnsPolynomialFailsIfThisIsNotNttForm) {
  int log_n = this->rns_context_->LogN();
  ASSERT_OK_AND_ASSIGN(RnsPolynomial<TypeParam> a_not_ntt,
                       RnsPolynomial<TypeParam>::CreateZero(
                           log_n, this->moduli_, /*is_ntt=*/false));
  ASSERT_OK_AND_ASSIGN(
      auto b, RnsPolynomial<TypeParam>::CreateZero(log_n, this->moduli_));

  EXPECT_THAT(a_not_ntt.MulInPlace(b, this->moduli_),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("RNS polynomial `this` must be in NTT form")));
}

TYPED_TEST(RnsPolynomialTest,
           MulInPlaceWithRnsPolynomialFailsIfThatIsNotNttForm) {
  int log_n = this->rns_context_->LogN();
  ASSERT_OK_AND_ASSIGN(
      auto a, RnsPolynomial<TypeParam>::CreateZero(log_n, this->moduli_));
  ASSERT_OK_AND_ASSIGN(
      auto b_not_ntt, RnsPolynomial<TypeParam>::CreateZero(log_n, this->moduli_,
                                                           /*is_ntt=*/false));

  EXPECT_THAT(a.MulInPlace(b_not_ntt, this->moduli_),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("RNS polynomial `that` must be in NTT form")));
}

TYPED_TEST(RnsPolynomialTest, MulWithRnsPolynomialIsCorrect) {
  ASSERT_OK_AND_ASSIGN(RnsPolynomial<TypeParam> a, this->SampleRnsPolynomial());

  // a * 1 == a.
  int log_n = this->rns_context_->LogN();
  ASSERT_OK_AND_ASSIGN(
      auto one, RnsPolynomial<TypeParam>::CreateOne(log_n, this->moduli_));
  ASSERT_OK_AND_ASSIGN(RnsPolynomial<TypeParam> b, a.Mul(one, this->moduli_));
  EXPECT_EQ(b, a);

  // a * 0 == 0.
  ASSERT_OK_AND_ASSIGN(
      auto zero, RnsPolynomial<TypeParam>::CreateZero(log_n, this->moduli_));
  ASSERT_OK_AND_ASSIGN(RnsPolynomial<TypeParam> c, a.Mul(zero, this->moduli_));
  EXPECT_EQ(c, zero);

  // Multiply with a monomial X^k
  constexpr int k = 5;
  int num_coeffs = 1 << log_n;
  auto mod_params_q0 = this->moduli_[0]->ModParams();
  std::vector<TypeParam> monomial_coeffs_q0(
      num_coeffs, TypeParam::ImportZero(mod_params_q0));
  monomial_coeffs_q0[k] = TypeParam::ImportOne(mod_params_q0);
  ASSERT_OK_AND_ASSIGN(RnsPolynomial<TypeParam> monomial,
                       RnsPolynomial<TypeParam>::ConvertFromPolynomialCoeffs(
                           monomial_coeffs_q0, mod_params_q0, this->moduli_));
  ASSERT_OK_AND_ASSIGN(RnsPolynomial<TypeParam> d,
                       a.Mul(monomial, this->moduli_));
  ASSERT_EQ(d.NumCoeffs(), num_coeffs);
  ASSERT_EQ(d.NumModuli(), this->moduli_.size());
  for (int i = 0; i < this->moduli_.size(); ++i) {
    auto mod_params_qi = this->moduli_[i]->ModParams();
    auto ntt_params_qi = this->moduli_[i]->NttParams();
    // Get a(X) (mod q_i) and d(X) (mod q_i) in coefficient form.
    Polynomial<TypeParam> d_qi(d.Coeffs()[i]);
    Polynomial<TypeParam> a_qi(a.Coeffs()[i]);
    std::vector<TypeParam> d_qi_coeffs =
        d_qi.InverseNtt(ntt_params_qi, mod_params_qi);
    std::vector<TypeParam> a_qi_coeffs =
        a_qi.InverseNtt(ntt_params_qi, mod_params_qi);
    // Assume a(X) = a0 + a1 * X + .. + a{n-1} * X^(n-1).
    // Then a(X) * X^k = a0 * X^k + a1 * X^(k+1) + .. + a{n-1-k} * X^(n-1)
    //                 - a{n-1-k+1} * X^0 - .. - a{n-1} * X^{k-1} (mod f(X)),
    // where f(X) = X^N + 1.
    int j = 0;
    for (; j + k < num_coeffs; ++j) {
      auto x = d_qi_coeffs[j + k];
      auto y = a_qi_coeffs[j];
      EXPECT_EQ(x, y);
    }
    for (; j < num_coeffs; ++j) {
      int index = (j + k) % num_coeffs;
      EXPECT_EQ(d_qi_coeffs[index].ExportInt(mod_params_qi),
                a_qi_coeffs[j].Negate(mod_params_qi).ExportInt(mod_params_qi));
    }
  }
}

// Mul(RnsPolynomial) should be equivalent to MulInPlace(RnsPolynomial).
TYPED_TEST(RnsPolynomialTest, MulInPlaceWithRnsPolynomialIsCorrect) {
  ASSERT_OK_AND_ASSIGN(RnsPolynomial<TypeParam> a, this->SampleRnsPolynomial());
  ASSERT_OK_AND_ASSIGN(RnsPolynomial<TypeParam> b, this->SampleRnsPolynomial());
  ASSERT_OK_AND_ASSIGN(RnsPolynomial<TypeParam> c, a.Mul(b, this->moduli_));
  RnsPolynomial<TypeParam> d = a;
  ASSERT_OK(d.MulInPlace(b, this->moduli_));
  EXPECT_EQ(d, c);
}

TYPED_TEST(RnsPolynomialTest, FusedMulAddInPlaceFailsIfWrongNumberOfModuli) {
  int num_moduli = this->moduli_.size();
  if (num_moduli <= 1) {
    return;  // This test requires at least two prime moduli.
  }

  int log_n = this->rns_context_->LogN();
  ASSERT_OK_AND_ASSIGN(
      auto f, RnsPolynomial<TypeParam>::CreateZero(log_n, this->moduli_));
  ASSERT_OK_AND_ASSIGN(
      auto h, RnsPolynomial<TypeParam>::CreateZero(log_n, this->moduli_));
  ASSERT_OK_AND_ASSIGN(
      auto g_missing_last_modulus,
      RnsPolynomial<TypeParam>::CreateZero(
          log_n, absl::MakeSpan(this->moduli_).subspan(0, num_moduli - 1)));

  EXPECT_THAT(f.FusedMulAddInPlace(g_missing_last_modulus, h, this->moduli_),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr(absl::StrCat("`a` must contain ", num_moduli,
                                              " coefficient vectors"))));
  EXPECT_THAT(f.FusedMulAddInPlace(h, g_missing_last_modulus, this->moduli_),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr(absl::StrCat("`b` must contain ", num_moduli,
                                              " coefficient vectors"))));
  EXPECT_THAT(
      f.FusedMulAddInPlace(
          h, h, absl::MakeSpan(this->moduli_).subspan(0, num_moduli - 1)),
      StatusIs(absl::StatusCode::kInvalidArgument,
               HasSubstr(absl::StrCat("`moduli` must contain ", num_moduli,
                                      " RNS moduli"))));
}

// a.FusedMulAddInPlace(b, c, moduli) should be equivalent to a += b * c.
TYPED_TEST(RnsPolynomialTest, FusedMulAddInPlaceIsCorrect) {
  ASSERT_OK_AND_ASSIGN(RnsPolynomial<TypeParam> a, this->SampleRnsPolynomial());
  ASSERT_OK_AND_ASSIGN(RnsPolynomial<TypeParam> b, this->SampleRnsPolynomial());
  ASSERT_OK_AND_ASSIGN(RnsPolynomial<TypeParam> c, this->SampleRnsPolynomial());

  RnsPolynomial<TypeParam> result = a;
  ASSERT_OK(result.FusedMulAddInPlace(b, c, this->moduli_));
  EXPECT_EQ(result.NumCoeffs(), (1 << this->rns_context_->LogN()));
  EXPECT_EQ(result.NumModuli(), this->moduli_.size());

  ASSERT_OK_AND_ASSIGN(RnsPolynomial<TypeParam> product,
                       b.Mul(c, this->moduli_));
  ASSERT_OK_AND_ASSIGN(RnsPolynomial<TypeParam> expected,
                       a.Add(product, this->moduli_));
  EXPECT_EQ(result, expected);
}

TYPED_TEST(RnsPolynomialTest, SubstituteFailsIfPowerIsNegative) {
  ASSERT_OK_AND_ASSIGN(RnsPolynomial<TypeParam> a,
                       RnsPolynomial<TypeParam>::CreateZero(
                           this->rns_context_->LogN(), this->moduli_));
  EXPECT_THAT(a.Substitute(-1, this->moduli_),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("Substitution power must be a non-negative "
                                 "odd integer less than 2*n")));
}

TYPED_TEST(RnsPolynomialTest, SubstituteFailsIfPowerIsEven) {
  ASSERT_OK_AND_ASSIGN(RnsPolynomial<TypeParam> a,
                       RnsPolynomial<TypeParam>::CreateZero(
                           this->rns_context_->LogN(), this->moduli_));
  EXPECT_THAT(a.Substitute(/*power=*/0, this->moduli_),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("Substitution power must be a non-negative "
                                 "odd integer less than 2*n")));
  EXPECT_THAT(a.Substitute(/*power=*/2, this->moduli_),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("Substitution power must be a non-negative "
                                 "odd integer less than 2*n")));
  EXPECT_THAT(a.Substitute(/*power=*/4, this->moduli_),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("Substitution power must be a non-negative "
                                 "odd integer less than 2*n")));
}

TYPED_TEST(RnsPolynomialTest, SubstituteFailsIfPowerIsOutOfBound) {
  int log_n = this->rns_context_->LogN();
  ASSERT_OK_AND_ASSIGN(
      RnsPolynomial<TypeParam> a,
      RnsPolynomial<TypeParam>::CreateZero(log_n, this->moduli_));
  int too_large_power = (1 << log_n) * 2 + 1;  // 2 * n + 1
  EXPECT_THAT(a.Substitute(too_large_power, this->moduli_),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("Substitution power must be a non-negative "
                                 "odd integer less than 2*n")));
}

TYPED_TEST(RnsPolynomialTest, SubstituteFailsIfThisIsNotNttForm) {
  ASSERT_OK_AND_ASSIGN(
      RnsPolynomial<TypeParam> a_not_ntt,
      RnsPolynomial<TypeParam>::CreateZero(this->rns_context_->LogN(),
                                           this->moduli_, /*is_ntt=*/false));
  EXPECT_THAT(a_not_ntt.Substitute(1, this->moduli_),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("RnsPolynomial must be in NTT form")));
}

TYPED_TEST(RnsPolynomialTest, SubstituteIsCorrect) {
  ASSERT_OK_AND_ASSIGN(RnsPolynomial<TypeParam> a, this->SampleRnsPolynomial());

  // First, substitute with power = 1 should be an noop.
  ASSERT_OK_AND_ASSIGN(RnsPolynomial<TypeParam> a_sub_with_identity,
                       a.Substitute(1, this->moduli_));
  ASSERT_EQ(a_sub_with_identity, a);

  // Next, substitute a with a non-trivial power and check for a (mod q_i).
  constexpr int k_power = 5;
  ASSERT_OK_AND_ASSIGN(RnsPolynomial<TypeParam> a_sub,
                       a.Substitute(k_power, this->moduli_));
  for (int i = 0; i < this->moduli_.size(); ++i) {
    // a(X^k) (mod q_i) should be the same as a_i(X^k) where a_i = a (mod q_i).
    Polynomial<TypeParam> a_i(a.Coeffs()[i]);
    ASSERT_OK_AND_ASSIGN(Polynomial<TypeParam> a_i_sub,
                         a_i.Substitute(k_power, this->moduli_[i]->NttParams(),
                                        this->moduli_[i]->ModParams()));
    Polynomial<TypeParam> a_sub_i(a_sub.Coeffs()[i]);
    EXPECT_EQ(a_sub_i, a_i_sub);
  }
}

TYPED_TEST(RnsPolynomialTest, ModReduceLsbFailsIfNumberOfModuliIsOne) {
  // Create a RnsPolynomial with only one prime moduli.
  int num_coeffs = 1 << this->rns_context_->LogN();
  std::vector<TypeParam> coeff_q0(
      num_coeffs, TypeParam::ImportZero(this->moduli_[0]->ModParams()));
  ASSERT_OK_AND_ASSIGN(
      auto a, RnsPolynomial<TypeParam>::Create({coeff_q0}, /*is_ntt=*/false));

  // We cannot apply modulus reduction on a with just one prime moduli.
  EXPECT_THAT(
      a.ModReduceLsb(this->rns_context_->PlaintextModulus(),
                     RnsInt<TypeParam>{},
                     absl::MakeSpan(this->moduli_).subspan(0, 1)),
      StatusIs(
          absl::StatusCode::kInvalidArgument,
          HasSubstr("ModReduceLsb cannot apply with only one prime modulus")));
}

TYPED_TEST(RnsPolynomialTest, ModReduceLsbFailsIfQlInvHasWrongNumberOfModuli) {
  if (this->moduli_.size() <= 1) {
    // This test requires at least two prime moduli.
    return;
  }

  ASSERT_OK_AND_ASSIGN(auto a, RnsPolynomial<TypeParam>::CreateOne(
                                   this->rns_context_->LogN(), this->moduli_));

  int level = this->moduli_.size() - 1;
  const RnsInt<TypeParam>& ql_inv_with_wrong_number_of_moduli =
      this->rns_context_->MainPrimeModulusInverseResidues()[level].Prefix(
          level - 1);
  EXPECT_THAT(a.ModReduceLsb(this->rns_context_->PlaintextModulus(),
                             ql_inv_with_wrong_number_of_moduli, this->moduli_),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr(absl::StrCat(
                           "`ql_inv` must be defined with respect to ", level,
                           " moduli"))));
}

TYPED_TEST(RnsPolynomialTest, ModReduceLsbIsCorrectForSmallCoefficients) {
  if (this->moduli_.size() <= 1) {
    // ModReduceLsb is applicable only if there are 2 or more prime moduli.
    return;
  }

  typename TypeParam::Int t = this->rns_context_->PlaintextModulus();

  // Compute the constants related to the last prime modulus q_l.
  int log_n = this->rns_context_->LogN();
  int level = this->moduli_.size() - 1;
  const RnsInt<TypeParam>& ql_inv =
      this->rns_context_->MainPrimeModulusInverseResidues()[level].Prefix(
          level);

  // For polynomial whose coefficients are smaller than a single prime modulus,
  // apply mod reduction to it should result in 0.
  absl::Span<const PrimeModulus<TypeParam>* const> moduli(this->moduli_);
  auto mod_params_q0 = moduli[0]->ModParams();
  ASSERT_OK_AND_ASSIGN(std::vector<TypeParam> coeffs_q0,
                       this->SampleCoeffs(log_n, mod_params_q0));
  ASSERT_OK_AND_ASSIGN(auto rns_poly,
                       RnsPolynomial<TypeParam>::ConvertFromPolynomialCoeffs(
                           coeffs_q0, mod_params_q0, moduli));
  ASSERT_OK_AND_ASSIGN(RnsPolynomial<TypeParam> rns_poly_reduced,
                       rns_poly.Mul(t, moduli));
  ASSERT_OK(rns_poly_reduced.ModReduceLsb(t, ql_inv, moduli));

  ASSERT_OK_AND_ASSIGN(
      auto zero_at_level_minus_one,
      RnsPolynomial<TypeParam>::CreateZero(log_n, moduli.subspan(0, level)));
  EXPECT_EQ(rns_poly_reduced, zero_at_level_minus_one);
}

TYPED_TEST(RnsPolynomialTest, ModReduceLsbIsCorrectForLargeCoefficients) {
  using Integer = typename TypeParam::Int;
  using BigInteger = typename testing::CompositeModulus<Integer>::value_type;

  if (this->moduli_.size() <= 1) {
    // ModReduceLsb is applicable only if there are 2 or more prime moduli.
    return;
  }

  // Compute the constants related to the last prime modulus q_l.
  int log_n = this->rns_context_->LogN();
  int level = this->moduli_.size() - 1;
  const RnsInt<TypeParam>& ql_inv =
      this->rns_context_->MainPrimeModulusInverseResidues()[level].Prefix(
          level);

  // If a(X) = t * e(X) + m(X), where e's coefficients are negative and with
  // small magnitude, the result of mod reduction on a(X) from Q to Q/q_l
  // should still be congruent to m(X) mod t.
  // Note that we represent a negative value e using residue Q - e (mod Q), so
  // for e with small absolute value, e (mod Q) is large and close to Q.
  Integer t = this->rns_context_->PlaintextModulus();
  absl::Span<const PrimeModulus<TypeParam>* const> moduli(this->moduli_);
  auto mod_params_q0 = moduli[0]->ModParams();
  ASSERT_OK_AND_ASSIGN(
      TypeParam neg_one,
      TypeParam::ImportInt(mod_params_q0->modulus - 1, mod_params_q0));
  std::vector<TypeParam> e_coeffs_q0((1 << log_n), neg_one);
  ASSERT_OK_AND_ASSIGN(
      auto e, RnsPolynomial<TypeParam>::ConvertBalancedFromPolynomialCoeffs(
                  e_coeffs_q0, mod_params_q0, moduli));
  ASSERT_OK_AND_ASSIGN(auto a, e.Mul(t, moduli));
  ASSERT_OK_AND_ASSIGN(
      std::vector<TypeParam> m_coeffs_q0,
      this->SampleCoeffs(log_n, mod_params_q0, static_cast<Uint64>(t)));
  ASSERT_OK_AND_ASSIGN(
      auto m, RnsPolynomial<TypeParam>::ConvertBalancedFromPolynomialCoeffs(
                  m_coeffs_q0, mod_params_q0, moduli));
  ASSERT_OK(a.AddInPlace(m, moduli));
  // Mod reduce from Q = q_0 * .. * q_L to Q/q_L.
  ASSERT_OK(a.ModReduceLsb(t, ql_inv, moduli));

  // Compute CRT interpolation of the mod-reduced RNS coefficients.
  auto moduli_reduced = moduli.subspan(0, level);
  ASSERT_OK(a.ConvertToCoeffForm(moduli_reduced));
  std::vector<BigInteger> modulus_hats =
      RnsModulusComplements<TypeParam, BigInteger>(moduli_reduced);
  ASSERT_OK_AND_ASSIGN(
      std::vector<TypeParam> modulus_hat_invs,
      this->rns_context_->MainPrimeModulusCrtFactors(level - 1));
  ASSERT_OK_AND_ASSIGN(
      std::vector<BigInteger> interpolated_coeffs,
      (CrtInterpolation<TypeParam, BigInteger>(
          a.Coeffs(), moduli_reduced, modulus_hats, modulus_hat_invs)));

  // Now let's check the mod-reduced coefficients represent m(X) (mod t).
  BigInteger q{1};
  for (auto prime_modulus : moduli_reduced) {
    q *= static_cast<BigInteger>(prime_modulus->Modulus());
  }
  BigInteger q_half = q >> 1;
  BigInteger t_big = static_cast<BigInteger>(t);
  for (int i = 0; i < (1 << log_n); ++i) {
    bool is_negative = interpolated_coeffs[i] > q_half;
    BigInteger x_abs =
        is_negative ? q - interpolated_coeffs[i] : interpolated_coeffs[i];
    BigInteger x_reminder = is_negative ? t_big - x_abs % t_big : x_abs % t_big;
    x_reminder = x_reminder % t_big;  // reduce t to 0
    BigInteger x_expected =
        static_cast<BigInteger>(m_coeffs_q0[i].ExportInt(mod_params_q0));
    EXPECT_EQ(x_reminder, x_expected);
  }
}

TYPED_TEST(RnsPolynomialTest,
           ModReduceLsbIsCorrectForArbitraryPlaintextModulus) {
  using Integer = typename TypeParam::Int;
  using BigInteger = typename testing::CompositeModulus<Integer>::value_type;

  if (this->moduli_.size() <= 1) {
    // ModReduceLsb is applicable only if there are 2 or more prime moduli.
    return;
  }

  // Use a different plaintext modulus t that doesn't satisfy q_i % t == 1 for
  // prime modulus q_i. ModReduceLsb should still be correct as long as the
  // reminder q_i % t is relatively small.
  Integer t = this->rns_context_->PlaintextModulus() + 1;

  // Compute the constants related to the last prime modulus q_l.
  int log_n = this->rns_context_->LogN();
  int level = this->moduli_.size() - 1;
  const RnsInt<TypeParam>& ql_inv =
      this->rns_context_->MainPrimeModulusInverseResidues()[level].Prefix(
          level);

  // Set a(X) = t * e(X) + m(X) (mod Q) for some error polynomial e(X) with
  // small coefficients (could be both positive and negative).
  absl::Span<const PrimeModulus<TypeParam>* const> moduli(this->moduli_);
  auto mod_params_q0 = moduli[0]->ModParams();
  ASSERT_OK_AND_ASSIGN(
      TypeParam neg_one,
      TypeParam::ImportInt(mod_params_q0->modulus - 1, mod_params_q0));
  std::vector<TypeParam> e_coeffs_q0((1 << log_n), neg_one);
  ASSERT_OK_AND_ASSIGN(
      auto e, RnsPolynomial<TypeParam>::ConvertBalancedFromPolynomialCoeffs(
                  e_coeffs_q0, mod_params_q0, moduli));
  ASSERT_OK_AND_ASSIGN(auto a, e.Mul(t, moduli));
  ASSERT_OK_AND_ASSIGN(
      std::vector<TypeParam> m_coeffs_q0,
      this->SampleCoeffs(log_n, mod_params_q0, static_cast<Uint64>(t)));
  ASSERT_OK_AND_ASSIGN(
      auto m, RnsPolynomial<TypeParam>::ConvertBalancedFromPolynomialCoeffs(
                  m_coeffs_q0, mod_params_q0, moduli));
  ASSERT_OK(a.AddInPlace(m, moduli));
  // Mod reduce from Q = q_0 * .. * q_L to Q/q_L.
  ASSERT_OK(a.ModReduceLsb(t, ql_inv, moduli));

  // Compute CRT interpolation of the mod-reduced RNS coefficients.
  auto moduli_reduced = moduli.subspan(0, level);
  ASSERT_OK(a.ConvertToCoeffForm(moduli_reduced));
  std::vector<BigInteger> modulus_hats =
      RnsModulusComplements<TypeParam, BigInteger>(moduli_reduced);
  ASSERT_OK_AND_ASSIGN(
      std::vector<TypeParam> modulus_hat_invs,
      this->rns_context_->MainPrimeModulusCrtFactors(level - 1));
  ASSERT_OK_AND_ASSIGN(
      std::vector<BigInteger> interpolated_coeffs,
      (CrtInterpolation<TypeParam, BigInteger>(
          a.Coeffs(), moduli_reduced, modulus_hats, modulus_hat_invs)));

  // Now let's check the mod-reduced coefficients represent m(X) (mod t).
  BigInteger q{1};
  for (auto prime_modulus : moduli_reduced) {
    q *= static_cast<BigInteger>(prime_modulus->Modulus());
  }
  BigInteger q_half = q >> 1;
  BigInteger t_big = static_cast<BigInteger>(t);
  for (int i = 0; i < (1 << log_n); ++i) {
    bool is_negative = interpolated_coeffs[i] > q_half;
    BigInteger x_abs =
        is_negative ? q - interpolated_coeffs[i] : interpolated_coeffs[i];
    BigInteger x_reminder = is_negative ? t_big - x_abs % t_big : x_abs % t_big;
    x_reminder = x_reminder % t_big;  // reduce t to 0
    BigInteger x_expected =
        static_cast<BigInteger>(m_coeffs_q0[i].ExportInt(mod_params_q0));
    EXPECT_EQ(x_reminder, x_expected);
  }
}

TYPED_TEST(RnsPolynomialTest, ModReduceMsbFailsIfNumberOfModuliIsOne) {
  // Create a RnsPolynomial with only one prime moduli.
  auto moduli = absl::MakeSpan(this->moduli_).subspan(0, 1);
  ASSERT_OK_AND_ASSIGN(auto a, RnsPolynomial<TypeParam>::CreateOne(
                                   this->rns_context_->LogN(), moduli));

  // We cannot apply modulus reduction on a with just one prime moduli.
  int level = this->moduli_.size() - 1;
  const RnsInt<TypeParam>& ql_inv =
      this->rns_context_->MainPrimeModulusInverseResidues()[level].Prefix(
          level);
  EXPECT_THAT(
      a.ModReduceMsb(ql_inv, moduli),
      StatusIs(
          absl::StatusCode::kInvalidArgument,
          HasSubstr("ModReduceMsb cannot apply with only one prime modulus")));
}

TYPED_TEST(RnsPolynomialTest, ModReduceMsbFailsIfWrongNumberOfModuli) {
  ASSERT_OK_AND_ASSIGN(auto a, RnsPolynomial<TypeParam>::CreateOne(
                                   this->rns_context_->LogN(), this->moduli_));
  int num_moduli = this->moduli_.size();
  int level = num_moduli - 1;
  const RnsInt<TypeParam>& ql_inv =
      this->rns_context_->MainPrimeModulusInverseResidues()[level].Prefix(
          level);
  EXPECT_THAT(
      a.ModReduceMsb(ql_inv, absl::MakeSpan(this->moduli_).subspan(0, level)),
      StatusIs(absl::StatusCode::kInvalidArgument,
               HasSubstr(absl::StrCat("`moduli` must contain ", num_moduli,
                                      " RNS moduli"))));
}

TYPED_TEST(RnsPolynomialTest, ModReduceMsbFailsIfQlInvHasWrongNumberOfModuli) {
  if (this->moduli_.size() <= 1) {
    // This test requires at least two prime moduli.
    return;
  }

  ASSERT_OK_AND_ASSIGN(auto a, RnsPolynomial<TypeParam>::CreateOne(
                                   this->rns_context_->LogN(), this->moduli_));

  int level = this->moduli_.size() - 1;
  const RnsInt<TypeParam>& ql_inv_with_wrong_number_of_moduli =
      this->rns_context_->MainPrimeModulusInverseResidues()[level].Prefix(
          level - 1);
  EXPECT_THAT(a.ModReduceMsb(ql_inv_with_wrong_number_of_moduli, this->moduli_),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr(absl::StrCat(
                           "`ql_inv` must be defined with respect to ", level,
                           " moduli"))));
}

TYPED_TEST(RnsPolynomialTest, ModReduceMsbIsCorrect) {
  using Integer = typename TypeParam::Int;
  using BigInteger = typename testing::CompositeModulus<Integer>::value_type;

  if (this->moduli_.size() <= 1) {
    // ModReduceMsb is applicable only if there are 2 or more prime moduli.
    return;
  }

  // For polynomial a (mod Q) where Q = Q' * q_l, ModReduceMsb should change
  // the polynomial to (a - (a mod q_l)) / q_l = round(a / q_l) (mod Q').
  ASSERT_OK_AND_ASSIGN(RnsPolynomial<TypeParam> a, this->SampleRnsPolynomial());
  RnsPolynomial<TypeParam> b = a;

  int level = this->moduli_.size() - 1;
  const RnsInt<TypeParam>& ql_inv =
      this->rns_context_->MainPrimeModulusInverseResidues()[level].Prefix(
          level);
  ASSERT_TRUE(b.IsNttForm());
  ASSERT_OK(b.ModReduceMsb(ql_inv, this->moduli_));
  EXPECT_EQ(b.NumModuli(), level);  // b is now modulo Q/q_l
  EXPECT_TRUE(b.IsNttForm());

  // Compute CRT interpolation of both a and mod-reduced polynomial b.
  ASSERT_OK(a.ConvertToCoeffForm(this->moduli_));
  std::vector<BigInteger> modulus_hats =
      RnsModulusComplements<TypeParam, BigInteger>(this->moduli_);
  ASSERT_OK_AND_ASSIGN(std::vector<TypeParam> modulus_hat_invs,
                       this->rns_context_->MainPrimeModulusCrtFactors(level));
  ASSERT_OK_AND_ASSIGN(
      std::vector<BigInteger> a_coeffs,
      (CrtInterpolation<TypeParam, BigInteger>(
          a.Coeffs(), this->moduli_, modulus_hats, modulus_hat_invs)));

  auto moduli_reduced = absl::MakeSpan(this->moduli_).subspan(0, level);
  ASSERT_OK(b.ConvertToCoeffForm(moduli_reduced));
  std::vector<BigInteger> modulus_hats_reduced =
      RnsModulusComplements<TypeParam, BigInteger>(moduli_reduced);
  ASSERT_OK_AND_ASSIGN(
      std::vector<TypeParam> modulus_hat_invs_reduced,
      this->rns_context_->MainPrimeModulusCrtFactors(level - 1));
  ASSERT_OK_AND_ASSIGN(std::vector<BigInteger> b_coeffs,
                       (CrtInterpolation<TypeParam, BigInteger>(
                           b.Coeffs(), moduli_reduced, modulus_hats_reduced,
                           modulus_hat_invs_reduced)));

  BigInteger ql = static_cast<BigInteger>(this->moduli_[level]->Modulus());
  BigInteger ql_half = ql >> 1;
  for (int i = 0; i < (1 << this->rns_context_->LogN()); ++i) {
    BigInteger r = a_coeffs[i] % ql;
    BigInteger rounded = r < ql_half ? a_coeffs[i] - r : a_coeffs[i] + (ql - r);
    rounded /= ql;
    EXPECT_EQ(b_coeffs[i], rounded);
  }

  // Lastly, let's check that ModReduceMsb on a polynomial in Coefficient form
  // produces the same result as if the polynomial is in NTT form.
  RnsPolynomial<TypeParam> c = a;
  ASSERT_FALSE(c.IsNttForm());
  ASSERT_OK(c.ModReduceMsb(ql_inv, this->moduli_));
  ASSERT_FALSE(c.IsNttForm());
  ASSERT_FALSE(b.IsNttForm());  // b was converted to Coeff form.
  EXPECT_EQ(c, b);
}

TYPED_TEST(RnsPolynomialTest, ExactSwitchRnsBasisFailsIfPolynomialIsInNttForm) {
  int level = this->moduli_.size() - 1;
  ASSERT_GE(level, 0);
  ASSERT_OK_AND_ASSIGN(
      auto a, RnsPolynomial<TypeParam>::CreateZero(
                  this->rns_context_->LogN(), this->moduli_, /*is_ntt=*/true));

  ASSERT_OK_AND_ASSIGN(auto q_hat_inv_mod_qs,
                       this->rns_context_->MainPrimeModulusCrtFactors(level));
  ASSERT_OK_AND_ASSIGN(
      auto q_hat_mod_ps,
      this->rns_context_->MainPrimeModulusComplementResidues(level));
  ASSERT_OK_AND_ASSIGN(
      auto q_mod_ps, this->rns_context_->MainLeveledModulusAuxResidues(level));
  EXPECT_THAT(
      a.SwitchRnsBasis(this->moduli_, this->rns_context_->AuxPrimeModuli(),
                       q_hat_inv_mod_qs, q_hat_mod_ps, q_mod_ps.RnsRep()),
      StatusIs(absl::StatusCode::kInvalidArgument,
               HasSubstr("Polynomial must be in coefficient form")));
}

TYPED_TEST(RnsPolynomialTest, ExactSwitchRnsBasisFailsIfWrongNumberOfModuli) {
  int level = this->moduli_.size() - 1;
  ASSERT_GE(level, 0);
  ASSERT_OK_AND_ASSIGN(
      auto a, RnsPolynomial<TypeParam>::CreateZero(
                  this->rns_context_->LogN(), this->moduli_, /*is_ntt=*/false));

  ASSERT_OK_AND_ASSIGN(auto q_hat_inv_mod_qs,
                       this->rns_context_->MainPrimeModulusCrtFactors(level));
  ASSERT_OK_AND_ASSIGN(
      auto q_hat_mod_ps,
      this->rns_context_->MainPrimeModulusComplementResidues(level));
  ASSERT_OK_AND_ASSIGN(
      auto q_mod_ps, this->rns_context_->MainLeveledModulusAuxResidues(level));
  EXPECT_THAT(
      a.SwitchRnsBasis(absl::MakeSpan(this->moduli_).subspan(0, level),
                       this->rns_context_->AuxPrimeModuli(), q_hat_inv_mod_qs,
                       q_hat_mod_ps, q_mod_ps.RnsRep()),
      StatusIs(absl::StatusCode::kInvalidArgument,
               HasSubstr(absl::StrCat("`this_moduli` must contain ",
                                      this->moduli_.size(), " RNS moduli"))));
}

TYPED_TEST(RnsPolynomialTest, ExactSwitchRnsBasisFailsIfEmptyOutputModuli) {
  int level = this->moduli_.size() - 1;
  ASSERT_GE(level, 0);
  ASSERT_OK_AND_ASSIGN(
      auto a, RnsPolynomial<TypeParam>::CreateZero(
                  this->rns_context_->LogN(), this->moduli_, /*is_ntt=*/false));

  ASSERT_OK_AND_ASSIGN(auto q_hat_inv_mod_qs,
                       this->rns_context_->MainPrimeModulusCrtFactors(level));
  ASSERT_OK_AND_ASSIGN(
      auto q_hat_mod_ps,
      this->rns_context_->MainPrimeModulusComplementResidues(level));
  ASSERT_OK_AND_ASSIGN(
      auto q_mod_ps, this->rns_context_->MainLeveledModulusAuxResidues(level));
  EXPECT_THAT(a.SwitchRnsBasis(this->moduli_,
                               /*output_moduli=*/{}, q_hat_inv_mod_qs,
                               q_hat_mod_ps, q_mod_ps.RnsRep()),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`output_moduli` must not be empty")));
}

TYPED_TEST(RnsPolynomialTest,
           ExactSwitchRnsBasisFailsIfWrongLevelForQHatInvModQs) {
  // This test requires the main moduli have at least two elements.
  if (this->moduli_.size() <= 1) {
    return;
  }
  int level = this->moduli_.size() - 1;
  ASSERT_GE(level, 0);
  ASSERT_OK_AND_ASSIGN(
      auto a, RnsPolynomial<TypeParam>::CreateZero(
                  this->rns_context_->LogN(), this->moduli_, /*is_ntt=*/false));

  ASSERT_OK_AND_ASSIGN(
      auto q_hat_inv_mod_qs_at_wrong_level,
      this->rns_context_->MainPrimeModulusCrtFactors(level - 1));
  ASSERT_OK_AND_ASSIGN(
      auto q_hat_mod_ps,
      this->rns_context_->MainPrimeModulusComplementResidues(level));
  ASSERT_OK_AND_ASSIGN(
      auto q_mod_ps, this->rns_context_->MainLeveledModulusAuxResidues(level));
  EXPECT_THAT(
      a.SwitchRnsBasis(this->moduli_, this->rns_context_->AuxPrimeModuli(),
                       q_hat_inv_mod_qs_at_wrong_level, q_hat_mod_ps,
                       q_mod_ps.RnsRep()),
      StatusIs(absl::StatusCode::kInvalidArgument,
               HasSubstr("`prime_q_hat_inv_mod_qs` must contain at least")));
}

TYPED_TEST(RnsPolynomialTest,
           ExactSwitchRnsBasisFailsIfWrongLevelForQHatModPs) {
  // This test requires the main moduli have at least two elements.
  if (this->moduli_.size() <= 1) {
    return;
  }
  int level = this->moduli_.size() - 1;
  ASSERT_GE(level, 0);
  ASSERT_OK_AND_ASSIGN(
      auto a, RnsPolynomial<TypeParam>::CreateZero(
                  this->rns_context_->LogN(), this->moduli_, /*is_ntt=*/false));

  ASSERT_OK_AND_ASSIGN(auto q_hat_inv_mod_qs,
                       this->rns_context_->MainPrimeModulusCrtFactors(level));
  ASSERT_OK_AND_ASSIGN(
      auto q_hat_mod_ps_at_wrong_level,
      this->rns_context_->MainPrimeModulusComplementResidues(level - 1));
  ASSERT_OK_AND_ASSIGN(
      auto q_mod_ps, this->rns_context_->MainLeveledModulusAuxResidues(level));
  EXPECT_THAT(
      a.SwitchRnsBasis(this->moduli_, this->rns_context_->AuxPrimeModuli(),
                       q_hat_inv_mod_qs, q_hat_mod_ps_at_wrong_level,
                       q_mod_ps.RnsRep()),
      StatusIs(absl::StatusCode::kInvalidArgument,
               HasSubstr("`prime_q_hat_mod_ps` must contain at least")));
}

TYPED_TEST(RnsPolynomialTest, ExactSwitchRnsBasisFailsIfWrongSizeForQModPs) {
  int level = this->moduli_.size() - 1;
  ASSERT_GE(level, 0);
  ASSERT_OK_AND_ASSIGN(
      auto a, RnsPolynomial<TypeParam>::CreateZero(
                  this->rns_context_->LogN(), this->moduli_, /*is_ntt=*/false));

  ASSERT_OK_AND_ASSIGN(auto q_hat_inv_mod_qs,
                       this->rns_context_->MainPrimeModulusCrtFactors(level));
  ASSERT_OK_AND_ASSIGN(
      auto q_hat_mod_ps,
      this->rns_context_->MainPrimeModulusComplementResidues(level));
  ASSERT_OK_AND_ASSIGN(
      auto q_mod_ps_wrong_size,
      this->rns_context_->MainLeveledModulusAuxResidues(level));
  q_mod_ps_wrong_size.zs.pop_back();
  EXPECT_THAT(a.SwitchRnsBasis(
                  this->moduli_, this->rns_context_->AuxPrimeModuli(),
                  q_hat_inv_mod_qs, q_hat_mod_ps, q_mod_ps_wrong_size.RnsRep()),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`q_mod_ps` must contain")));
}

TYPED_TEST(RnsPolynomialTest,
           ScaleAndSwitchRnsBasisFailsIfPolynomialIsInNttForm) {
  int level = this->moduli_.size() - 1;
  ASSERT_GE(level, 0);
  ASSERT_OK_AND_ASSIGN(
      auto a, RnsPolynomial<TypeParam>::CreateZero(
                  this->rns_context_->LogN(), this->moduli_, /*is_ntt=*/true));

  ASSERT_OK_AND_ASSIGN(auto q_hat_inv_mod_qs,
                       this->rns_context_->MainPrimeModulusCrtFactors(level));
  EXPECT_THAT(
      a.ScaleAndSwitchRnsBasis(
          this->moduli_, this->rns_context_->AuxPrimeModuli(), q_hat_inv_mod_qs,
          this->rns_context_->MainPrimeModulusInverseAuxResidues(),
          this->rns_context_->AuxModulusResidues()),
      StatusIs(absl::StatusCode::kInvalidArgument,
               HasSubstr("Polynomial must be in coefficient form")));
}

TYPED_TEST(RnsPolynomialTest,
           ScaleAndSwitchRnsBasisFailsIfWrongNumberOfModuli) {
  int level = this->moduli_.size() - 1;
  ASSERT_GE(level, 0);
  ASSERT_OK_AND_ASSIGN(
      auto a, RnsPolynomial<TypeParam>::CreateZero(
                  this->rns_context_->LogN(), this->moduli_, /*is_ntt=*/false));

  ASSERT_OK_AND_ASSIGN(auto q_hat_inv_mod_qs,
                       this->rns_context_->MainPrimeModulusCrtFactors(level));
  EXPECT_THAT(
      a.ScaleAndSwitchRnsBasis(
          absl::MakeSpan(this->moduli_).subspan(0, level),
          this->rns_context_->AuxPrimeModuli(), q_hat_inv_mod_qs,
          this->rns_context_->MainPrimeModulusInverseAuxResidues(),
          this->rns_context_->AuxModulusResidues()),
      StatusIs(absl::StatusCode::kInvalidArgument,
               HasSubstr(absl::StrCat("`this_moduli` must contain ",
                                      this->moduli_.size(), " RNS moduli"))));
}

TYPED_TEST(RnsPolynomialTest, ScaleAndSwitchRnsBasisFailsIfEmptyOutputModuli) {
  int level = this->moduli_.size() - 1;
  ASSERT_GE(level, 0);
  ASSERT_OK_AND_ASSIGN(
      auto a, RnsPolynomial<TypeParam>::CreateZero(
                  this->rns_context_->LogN(), this->moduli_, /*is_ntt=*/false));

  ASSERT_OK_AND_ASSIGN(auto q_hat_inv_mod_qs,
                       this->rns_context_->MainPrimeModulusCrtFactors(level));
  EXPECT_THAT(a.ScaleAndSwitchRnsBasis(
                  this->moduli_,
                  /*output_moduli=*/{}, q_hat_inv_mod_qs,
                  this->rns_context_->MainPrimeModulusInverseAuxResidues(),
                  this->rns_context_->AuxModulusResidues()),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("output_moduli` must not be empty")));
}

TYPED_TEST(RnsPolynomialTest,
           ScaleAndSwitchRnsBasisFailsIfWrongLevelForQHatInvModQs) {
  // This test requires at least two main prime moduli.
  if (this->moduli_.size() <= 1) {
    return;
  }
  int level = this->moduli_.size() - 1;
  ASSERT_GE(level, 0);
  ASSERT_OK_AND_ASSIGN(
      auto a, RnsPolynomial<TypeParam>::CreateZero(
                  this->rns_context_->LogN(), this->moduli_, /*is_ntt=*/false));

  ASSERT_OK_AND_ASSIGN(
      auto q_hat_inv_mod_qs_at_wrong_level,
      this->rns_context_->MainPrimeModulusCrtFactors(level - 1));
  EXPECT_THAT(
      a.ScaleAndSwitchRnsBasis(
          this->moduli_, this->rns_context_->AuxPrimeModuli(),
          q_hat_inv_mod_qs_at_wrong_level,
          this->rns_context_->MainPrimeModulusInverseAuxResidues(),
          this->rns_context_->AuxModulusResidues()),
      StatusIs(absl::StatusCode::kInvalidArgument,
               HasSubstr("`prime_q_hat_inv_mod_qs` must contain at least")));
}

TYPED_TEST(RnsPolynomialTest,
           ScaleAndSwitchRnsBasisFailsIfWrongSizeForQInvModPs) {
  int level = this->moduli_.size() - 1;
  ASSERT_GE(level, 0);
  ASSERT_OK_AND_ASSIGN(
      auto a, RnsPolynomial<TypeParam>::CreateZero(
                  this->rns_context_->LogN(), this->moduli_, /*is_ntt=*/false));

  ASSERT_OK_AND_ASSIGN(auto q_hat_inv_mod_qs,
                       this->rns_context_->MainPrimeModulusCrtFactors(level));
  int num_aux_moduli = this->rns_context_->AuxPrimeModuli().size();
  auto q_inv_mod_ps_wrong_size =
      this->rns_context_->MainPrimeModulusInverseAuxResidues().subspan(
          0, num_aux_moduli - 1);
  EXPECT_THAT(
      a.ScaleAndSwitchRnsBasis(
          this->moduli_, this->rns_context_->AuxPrimeModuli(), q_hat_inv_mod_qs,
          q_inv_mod_ps_wrong_size, this->rns_context_->AuxModulusResidues()),
      StatusIs(absl::StatusCode::kInvalidArgument,
               HasSubstr("`prime_q_inv_mod_ps` must contain")));
}

TYPED_TEST(RnsPolynomialTest,
           ScaleAndSwitchRnsBasisFailsIfWrongLevelForPModQs) {
  // This test requires at least two main prime moduli.
  if (this->moduli_.size() <= 1) {
    return;
  }
  int level = this->moduli_.size() - 1;
  ASSERT_GE(level, 0);
  ASSERT_OK_AND_ASSIGN(
      auto a, RnsPolynomial<TypeParam>::CreateZero(
                  this->rns_context_->LogN(), this->moduli_, /*is_ntt=*/false));

  ASSERT_OK_AND_ASSIGN(auto q_hat_inv_mod_qs,
                       this->rns_context_->MainPrimeModulusCrtFactors(level));
  auto p_mod_qs_at_wrong_level =
      this->rns_context_->AuxModulusResidues().subspan(0, level);
  EXPECT_THAT(
      a.ScaleAndSwitchRnsBasis(
          this->moduli_, this->rns_context_->AuxPrimeModuli(), q_hat_inv_mod_qs,
          this->rns_context_->MainPrimeModulusInverseAuxResidues(),
          p_mod_qs_at_wrong_level),
      StatusIs(absl::StatusCode::kInvalidArgument,
               HasSubstr("`p_mod_qs` must contain at least")));
}

TYPED_TEST(RnsPolynomialTest, ScaleAndSwitchRnsBasisToSmallerModulus) {
  using Integer = typename TypeParam::Int;
  using ModularIntParams = typename TypeParam::Params;

  constexpr Integer t = 517;

  ASSERT_OK_AND_ASSIGN(std::unique_ptr<const ModularIntParams> params_t,
                       TypeParam::Params::Create(t));
  auto modulus_t =
      absl::WrapUnique(new PrimeModulus<TypeParam>{std::move(params_t),
                                                   /*ntt_params=*/nullptr});
  std::vector<const PrimeModulus<TypeParam>*> output_moduli = {modulus_t.get()};
  std::vector<const PrimeModulus<TypeParam>*> main_moduli = this->moduli_;
  std::vector<TypeParam> qs_mod_t;      // {q_i mod t}
  std::vector<TypeParam> q_invs_mod_t;  // {q_i^(-1) mod t}
  std::vector<TypeParam> t_mod_qs;      // {t mod q_i}
  std::vector<TypeParam> t_inv_mod_qs;  // {t^(-1) mod q_i}
  for (int i = 0; i < main_moduli.size(); ++i) {
    auto mod_params_qi = main_moduli[i]->ModParams();
    ASSERT_OK_AND_ASSIGN(
        TypeParam qi_mod_t,
        TypeParam::ImportInt(mod_params_qi->modulus, modulus_t->ModParams()));
    ASSERT_OK_AND_ASSIGN(
        TypeParam qi_inv_mod_t,
        qi_mod_t.MultiplicativeInverseFast(modulus_t->ModParams()));
    q_invs_mod_t.push_back(std::move(qi_inv_mod_t));
    qs_mod_t.push_back(std::move(qi_mod_t));

    ASSERT_OK_AND_ASSIGN(TypeParam t_mod_qi,
                         TypeParam::ImportInt(t, mod_params_qi));
    t_inv_mod_qs.push_back(t_mod_qi.MultiplicativeInverse(mod_params_qi));
    t_mod_qs.push_back(std::move(t_mod_qi));
  }

  // Sample random values v mod t.
  int log_n = this->rns_context_->LogN();
  ASSERT_OK_AND_ASSIGN(std::vector<TypeParam> v_mod_t,
                       this->SampleCoeffs(log_n, modulus_t->ModParams(),
                                          static_cast<Uint64>(t)));
  // Encode v (mod t) to a = ([-Q * a] mod t) / t (mod Q).
  std::vector<TypeParam> scaled_v_mod_t = v_mod_t;
  for (int i = 0; i < main_moduli.size(); ++i) {
    ASSERT_OK(TypeParam::BatchMulInPlace(&scaled_v_mod_t, qs_mod_t[i],
                                         modulus_t->ModParams()));
  }
  for (int i = 0; i < scaled_v_mod_t.size(); ++i) {
    scaled_v_mod_t[i].NegateInPlace(modulus_t->ModParams());
  }
  ASSERT_OK_AND_ASSIGN(
      RnsPolynomial<TypeParam> a,
      RnsPolynomial<TypeParam>::ConvertFromPolynomialCoeffs(
          scaled_v_mod_t, modulus_t->ModParams(), main_moduli));
  ASSERT_OK(a.ConvertToCoeffForm(main_moduli));
  ASSERT_OK(a.MulInPlace(t_inv_mod_qs, main_moduli));

  // Compute round(t / Q * a(X)) mod t.
  int level = main_moduli.size() - 1;
  ASSERT_OK_AND_ASSIGN(std::vector<TypeParam> q_hat_inv_mod_qs,
                       this->rns_context_->MainPrimeModulusCrtFactors(level));
  RnsInt<TypeParam> q_inv_mod_t{std::move(q_invs_mod_t)};
  ASSERT_OK_AND_ASSIGN(
      RnsPolynomial<TypeParam> b,
      a.ScaleAndSwitchRnsBasis(main_moduli, output_moduli, q_hat_inv_mod_qs,
                               {q_inv_mod_t}, t_mod_qs));
  EXPECT_EQ(b.NumModuli(), output_moduli.size());

  for (int i = 0; i < (1 << log_n); ++i) {
    Integer expected = v_mod_t[i].ExportInt(modulus_t->ModParams());
    Integer actual = b.Coeffs()[0][i].ExportInt(modulus_t->ModParams());
    EXPECT_EQ(actual, expected);
  }
}

TYPED_TEST(RnsPolynomialTest, ExtendRnsBasisInPlaceHasCorrectModuliSize) {
  // Sample a random polynomial mod Q.
  ASSERT_OK_AND_ASSIGN(RnsPolynomial<TypeParam> a,
                       this->SampleRnsPolynomial(this->moduli_));
  ASSERT_OK(a.ConvertToCoeffForm(this->moduli_));

  // Extend a (mod Q) to a (mod Q*P).
  int level = this->moduli_.size() - 1;
  ASSERT_OK_AND_ASSIGN(std::vector<TypeParam> q_hat_inv_mod_qs,
                       this->rns_context_->MainPrimeModulusCrtFactors(level));
  ASSERT_OK_AND_ASSIGN(
      std::vector<RnsInt<TypeParam>> q_hat_mod_ps,
      this->rns_context_->MainPrimeModulusComplementResidues(level));
  std::vector<const PrimeModulus<TypeParam>*> aux_moduli =
      this->rns_context_->AuxPrimeModuli();
  ASSERT_OK_AND_ASSIGN(
      auto q_mod_ps, this->rns_context_->MainLeveledModulusAuxResidues(level));
  ASSERT_OK(a.ExtendRnsBasisInPlace(this->moduli_, aux_moduli, q_hat_inv_mod_qs,
                                    q_hat_mod_ps, q_mod_ps.RnsRep()));
  EXPECT_EQ(a.NumModuli(), this->moduli_.size() + aux_moduli.size());
}

TYPED_TEST(RnsPolynomialTest, ApproxSwitchRnsBasisFailsIfWrongNumberOfModuli) {
  int level = this->moduli_.size() - 1;
  ASSERT_GE(level, 0);
  ASSERT_OK_AND_ASSIGN(auto a, RnsPolynomial<TypeParam>::CreateZero(
                                   this->rns_context_->LogN(), this->moduli_));

  ASSERT_OK_AND_ASSIGN(auto q_hat_inv_mod_qs,
                       this->rns_context_->MainPrimeModulusCrtFactors(level));
  ASSERT_OK_AND_ASSIGN(
      auto q_hat_mod_ps,
      this->rns_context_->MainPrimeModulusComplementResidues(level));
  EXPECT_THAT(
      a.ApproxSwitchRnsBasis(absl::MakeSpan(this->moduli_).subspan(0, level),
                             this->rns_context_->AuxPrimeModuli(),
                             q_hat_inv_mod_qs, q_hat_mod_ps),
      StatusIs(absl::StatusCode::kInvalidArgument,
               HasSubstr(absl::StrCat("`this_moduli` must contain ",
                                      this->moduli_.size(), " RNS moduli"))));
}

TYPED_TEST(RnsPolynomialTest, ApproxSwitchRnsBasisFailsIfEmptyOutputModuli) {
  int level = this->moduli_.size() - 1;
  ASSERT_GE(level, 0);
  ASSERT_OK_AND_ASSIGN(auto a, RnsPolynomial<TypeParam>::CreateZero(
                                   this->rns_context_->LogN(), this->moduli_));

  ASSERT_OK_AND_ASSIGN(auto q_hat_inv_mod_qs,
                       this->rns_context_->MainPrimeModulusCrtFactors(level));
  ASSERT_OK_AND_ASSIGN(
      auto q_hat_mod_ps,
      this->rns_context_->MainPrimeModulusComplementResidues(level));
  EXPECT_THAT(a.ApproxSwitchRnsBasis(absl::MakeSpan(this->moduli_),
                                     /*output_moduli=*/{}, q_hat_inv_mod_qs,
                                     q_hat_mod_ps),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`output_moduli` must not be empty")));
}

TYPED_TEST(RnsPolynomialTest, ApproxSwitchRnsBasisFailsIfWrongNumberOfQHatInv) {
  int level = this->moduli_.size() - 1;
  if (level <= 0) {
    // This test is valid only when there are more than one prime moduli.
    return;
  }

  ASSERT_OK_AND_ASSIGN(auto a, RnsPolynomial<TypeParam>::CreateZero(
                                   this->rns_context_->LogN(), this->moduli_));

  ASSERT_OK_AND_ASSIGN(
      auto q_hat_inv_mod_qs_at_wrong_level,
      this->rns_context_->MainPrimeModulusCrtFactors(level - 1));
  ASSERT_OK_AND_ASSIGN(
      auto q_hat_mod_ps,
      this->rns_context_->MainPrimeModulusComplementResidues(level));
  EXPECT_THAT(
      a.ApproxSwitchRnsBasis(absl::MakeSpan(this->moduli_),
                             this->rns_context_->AuxPrimeModuli(),
                             q_hat_inv_mod_qs_at_wrong_level, q_hat_mod_ps),
      StatusIs(absl::StatusCode::kInvalidArgument,
               HasSubstr(absl::StrCat("`prime_q_hat_inv_mod_qs` must contain ",
                                      this->moduli_.size(), " elements"))));
}

TYPED_TEST(RnsPolynomialTest, ApproxSwitchRnsBasisFailsIfWrongNumberOfQHat) {
  int level = this->moduli_.size() - 1;
  if (level <= 0) {
    // This test is valid only when there are more than one prime moduli.
    return;
  }

  ASSERT_OK_AND_ASSIGN(auto a, RnsPolynomial<TypeParam>::CreateZero(
                                   this->rns_context_->LogN(), this->moduli_));

  ASSERT_OK_AND_ASSIGN(auto q_hat_inv_mod_qs,
                       this->rns_context_->MainPrimeModulusCrtFactors(level));
  ASSERT_OK_AND_ASSIGN(
      auto q_hat_mod_ps_at_wrong_level,
      this->rns_context_->MainPrimeModulusComplementResidues(level - 1));
  EXPECT_THAT(
      a.ApproxSwitchRnsBasis(absl::MakeSpan(this->moduli_),
                             this->rns_context_->AuxPrimeModuli(),
                             q_hat_inv_mod_qs, q_hat_mod_ps_at_wrong_level),
      StatusIs(absl::StatusCode::kInvalidArgument,
               HasSubstr(absl::StrCat("`prime_q_hat_mod_ps` must contain ",
                                      this->moduli_.size(), " elements"))));
}

TYPED_TEST(RnsPolynomialTest, ApproxSwitchRnsBasisIsCorrect) {
  using Integer = typename TypeParam::Int;
  using BigInteger = typename testing::CompositeModulus<Integer>::value_type;

  // ApproxSwitchRnsBasis is correct when the auxiliary modulus P is larger than
  // the main modulus Q. So for this test we define Q using a subset of q_i's
  // whose product is smaller than P.
  std::vector<const PrimeModulus<TypeParam>*> aux_moduli =
      this->rns_context_->AuxPrimeModuli();
  BigInteger p{1};
  for (auto prime_modulus : aux_moduli) {
    p *= static_cast<BigInteger>(prime_modulus->Modulus());
  }
  std::vector<const PrimeModulus<TypeParam>*> main_moduli;
  BigInteger q{1};
  for (int i = 0; i < this->moduli_.size(); ++i) {
    BigInteger qq = q * static_cast<BigInteger>(this->moduli_[i]->Modulus());
    if (qq < p) {
      q = qq;
      main_moduli.push_back(this->moduli_[i]);
    } else {
      break;
    }
  }
  ASSERT_GE(main_moduli.size(), 1);

  // Sample a random polynomial mod Q.
  int log_n = this->rns_context_->LogN();
  ASSERT_OK_AND_ASSIGN(RnsPolynomial<TypeParam> a,
                       this->SampleRnsPolynomial(main_moduli));

  // Approximately convert a(X) from mod-Q to mod-P.
  int level = main_moduli.size() - 1;
  ASSERT_OK_AND_ASSIGN(std::vector<TypeParam> q_hat_inv_mod_qs,
                       this->rns_context_->MainPrimeModulusCrtFactors(level));
  ASSERT_OK_AND_ASSIGN(
      std::vector<RnsInt<TypeParam>> q_hat_mod_ps,
      this->rns_context_->MainPrimeModulusComplementResidues(level));
  ASSERT_OK_AND_ASSIGN(RnsPolynomial<TypeParam> b,
                       a.ApproxSwitchRnsBasis(main_moduli, aux_moduli,
                                              q_hat_inv_mod_qs, q_hat_mod_ps));
  ASSERT_EQ(b.NumModuli(), aux_moduli.size());
  ASSERT_FALSE(b.IsNttForm());

  // Compute CRT interpolation of the RNS polynomial a.
  ASSERT_OK(a.ConvertToCoeffForm(main_moduli));
  auto q_hats = RnsModulusComplements<TypeParam, BigInteger>(main_moduli);
  ASSERT_OK_AND_ASSIGN(auto q_hat_invs,
                       this->rns_context_->MainPrimeModulusCrtFactors(level));
  ASSERT_OK_AND_ASSIGN(std::vector<BigInteger> a_coeffs,
                       (CrtInterpolation<TypeParam, BigInteger>(
                           a.Coeffs(), main_moduli, q_hats, q_hat_invs)));

  // Compute CRT interpolation of the RNS polynomial b.
  auto p_hats = RnsModulusComplements<TypeParam, BigInteger>(aux_moduli);
  auto p_hat_invs = this->rns_context_->AuxPrimeModulusCrtFactors();
  ASSERT_OK_AND_ASSIGN(std::vector<BigInteger> b_coeffs,
                       (CrtInterpolation<TypeParam, BigInteger>(
                           b.Coeffs(), aux_moduli, p_hats, p_hat_invs)));

  // For each coefficient a[i], ApproxSwitchRnsBasis should ensure that b[i] =
  // (a[i] + d * Q) mod P, for some small integer d. So we expect b[i] mod Q is
  // the same as a[i].
  for (int i = 0; i < (1 << log_n); ++i) {
    BigInteger expected = a_coeffs[i];  // a[i] < q < p
    BigInteger actual = b_coeffs[i] % q;
    EXPECT_EQ(actual, expected);
  }

  // Lastly, we check ApproxSwitchRnsBasis on a polynomial in Coefficient form
  // should return the same as if the input polynomial is in NTT form.
  ASSERT_FALSE(a.IsNttForm());
  ASSERT_OK_AND_ASSIGN(RnsPolynomial<TypeParam> c,
                       a.ApproxSwitchRnsBasis(main_moduli, aux_moduli,
                                              q_hat_inv_mod_qs, q_hat_mod_ps));
  ASSERT_FALSE(c.IsNttForm());
  EXPECT_EQ(c, b);
}

TYPED_TEST(RnsPolynomialTest, ApproxModReduceLsbFailsIfWrongNumberOfModuli) {
  int num_moduli = this->moduli_.size();
  ASSERT_OK_AND_ASSIGN(RnsPolynomial<TypeParam> a,
                       RnsPolynomial<TypeParam>::CreateZero(
                           this->rns_context_->LogN(), this->moduli_));
  std::vector<const PrimeModulus<TypeParam>*> aux_moduli =
      this->rns_context_->AuxPrimeModuli();
  ASSERT_OK_AND_ASSIGN(RnsPolynomial<TypeParam> b,
                       RnsPolynomial<TypeParam>::CreateZero(
                           this->rns_context_->LogN(), aux_moduli));
  EXPECT_THAT(a.ApproxModReduceLsb(
                  b, absl::MakeSpan(this->moduli_).subspan(0, num_moduli - 1),
                  aux_moduli, this->rns_context_->AuxModulusInverseResidues(),
                  this->rns_context_->AuxPrimeModulusCrtFactors(),
                  this->rns_context_->AuxPrimeModulusComplementResidues(),
                  this->rns_context_->PlaintextModulusInverseAuxResidues(),
                  this->rns_context_->AuxModulusPlaintextResidue(),
                  this->rns_context_->PlaintextModulus()),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr(absl::StrCat("`this_moduli` must contain ",
                                              num_moduli, " RNS moduli"))));
}

TYPED_TEST(RnsPolynomialTest, ApproxModReduceLsbFailsIfWrongNumberOfAuxModuli) {
  std::vector<const PrimeModulus<TypeParam>*> aux_moduli =
      this->rns_context_->AuxPrimeModuli();
  if (aux_moduli.size() <= 1) {
    // This test requires at least two auxiliary moduli.
    return;
  }

  ASSERT_OK_AND_ASSIGN(RnsPolynomial<TypeParam> a,
                       RnsPolynomial<TypeParam>::CreateZero(
                           this->rns_context_->LogN(), this->moduli_));
  ASSERT_OK_AND_ASSIGN(RnsPolynomial<TypeParam> b,
                       RnsPolynomial<TypeParam>::CreateZero(
                           this->rns_context_->LogN(), aux_moduli));
  EXPECT_THAT(
      a.ApproxModReduceLsb(
          b, this->moduli_, absl::MakeSpan(aux_moduli).subspan(0, 1),
          this->rns_context_->AuxModulusInverseResidues(),
          this->rns_context_->AuxPrimeModulusCrtFactors(),
          this->rns_context_->AuxPrimeModulusComplementResidues(),
          this->rns_context_->PlaintextModulusInverseAuxResidues(),
          this->rns_context_->AuxModulusPlaintextResidue(),
          this->rns_context_->PlaintextModulus()),
      StatusIs(absl::StatusCode::kInvalidArgument,
               HasSubstr(absl::StrCat("`aux_moduli` must contain ",
                                      aux_moduli.size(), " RNS moduli"))));
}

TYPED_TEST(RnsPolynomialTest, ApproxModReduceLsbFailsIfWrongNumberOfConstants) {
  int num_moduli = this->moduli_.size();
  auto aux_moduli = this->rns_context_->AuxPrimeModuli();
  int num_aux_moduli = aux_moduli.size();
  if (num_moduli <= 1 || num_aux_moduli <= 1) {
    // This test requires >= 2main moduli and >= 2 aux moduli.
    return;
  }

  ASSERT_OK_AND_ASSIGN(RnsPolynomial<TypeParam> a,
                       RnsPolynomial<TypeParam>::CreateZero(
                           this->rns_context_->LogN(), this->moduli_));
  ASSERT_OK_AND_ASSIGN(RnsPolynomial<TypeParam> b,
                       RnsPolynomial<TypeParam>::CreateZero(
                           this->rns_context_->LogN(), aux_moduli));
  EXPECT_THAT(a.ApproxModReduceLsb(
                  b, this->moduli_, aux_moduli,
                  this->rns_context_->AuxModulusInverseResidues().subspan(
                      0, num_moduli - 1),
                  this->rns_context_->AuxPrimeModulusCrtFactors(),
                  this->rns_context_->AuxPrimeModulusComplementResidues(),
                  this->rns_context_->PlaintextModulusInverseAuxResidues(),
                  this->rns_context_->AuxModulusPlaintextResidue(),
                  this->rns_context_->PlaintextModulus()),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr(absl::StrCat("`p_inv_mod_qs` must contain ",
                                              num_moduli, " elements"))));
  EXPECT_THAT(
      a.ApproxModReduceLsb(
          b, this->moduli_, aux_moduli,
          this->rns_context_->AuxModulusInverseResidues(),
          this->rns_context_->AuxPrimeModulusCrtFactors().subspan(
              0, num_aux_moduli - 1),
          this->rns_context_->AuxPrimeModulusComplementResidues(),
          this->rns_context_->PlaintextModulusInverseAuxResidues(),
          this->rns_context_->AuxModulusPlaintextResidue(),
          this->rns_context_->PlaintextModulus()),
      StatusIs(absl::StatusCode::kInvalidArgument,
               HasSubstr(absl::StrCat("`prime_p_hat_inv_mod_ps` must contain ",
                                      num_aux_moduli, " elements"))));

  EXPECT_THAT(
      a.ApproxModReduceLsb(
          b, this->moduli_, aux_moduli,
          this->rns_context_->AuxModulusInverseResidues(),
          this->rns_context_->AuxPrimeModulusCrtFactors(),
          this->rns_context_->AuxPrimeModulusComplementResidues().subspan(
              0, num_aux_moduli - 1),
          this->rns_context_->PlaintextModulusInverseAuxResidues(),
          this->rns_context_->AuxModulusPlaintextResidue(),
          this->rns_context_->PlaintextModulus()),
      StatusIs(absl::StatusCode::kInvalidArgument,
               HasSubstr(absl::StrCat("`prime_p_hat_mod_qs` must contain ",
                                      num_aux_moduli, " elements"))));

  EXPECT_THAT(
      a.ApproxModReduceLsb(
          b, this->moduli_, aux_moduli,
          this->rns_context_->AuxModulusInverseResidues(),
          this->rns_context_->AuxPrimeModulusCrtFactors(),
          this->rns_context_->AuxPrimeModulusComplementResidues(),
          this->rns_context_->PlaintextModulusInverseAuxResidues().subspan(
              0, num_aux_moduli - 1),
          this->rns_context_->AuxModulusPlaintextResidue(),
          this->rns_context_->PlaintextModulus()),
      StatusIs(absl::StatusCode::kInvalidArgument,
               HasSubstr(absl::StrCat("`t_inv_mod_ps` must contain ",
                                      num_aux_moduli, " elements"))));
}

// Checks that invoking `ApproxModReduceLsb()` with mod-expanded polynomials
// (i.e. generated by `ApproxModSwitchRnsBasis` to auxiliary modulus P) reduces
// modulus from Q*P to Q and keeps plaintext messages at LSB of modulus.
TYPED_TEST(RnsPolynomialTest, ApproxModReduceLsbIsCorrect) {
  using Integer = typename TypeParam::Int;
  using BigInteger = typename testing::CompositeModulus<Integer>::value_type;

  // ApproxSwitchRnsBasis is correct when the auxiliary modulus P is larger than
  // the main modulus Q. So for this test we define Q using a subset of q_i's
  // whose product is smaller than P.
  std::vector<const PrimeModulus<TypeParam>*> aux_moduli =
      this->rns_context_->AuxPrimeModuli();
  BigInteger p{1};
  for (auto prime_modulus : aux_moduli) {
    p *= static_cast<BigInteger>(prime_modulus->Modulus());
  }
  std::vector<const PrimeModulus<TypeParam>*> main_moduli;
  BigInteger q{1};
  for (int i = 0; i < this->moduli_.size(); ++i) {
    BigInteger qq = q * static_cast<BigInteger>(this->moduli_[i]->Modulus());
    if (qq < p) {
      q = qq;
      main_moduli.push_back(this->moduli_[i]);
    } else {
      break;
    }
  }
  ASSERT_GE(main_moduli.size(), 1);

  // We do the following steps to make sure ApproxModReduceLsb is correct:
  // 1. Start with m + t * e (mod Q), for some plaintext m, plaintext modulus t
  //    and small e;
  // 2. Let a = P * (m + t * e) + t * f (mod Q), i.e. scale up by P and add
  //    small error t*f;
  // 3. Let b = t * f (mod P) == P * (m + t * e) + t * f (mod P).
  // At this point the pair (a, b) represents c = P*m + P*t*e + u*Q + t*f
  // (mod Q*P) for small u. Then a.ApproxModReduceLsb(b) should update a to
  // floor(c/P) = m + t * e' (mod Q) without destroying m in the LSB of Q.
  int log_n = this->rns_context_->LogN();
  Integer t = this->rns_context_->PlaintextModulus();
  auto mod_params_q0 = main_moduli[0]->ModParams();
  ASSERT_OK_AND_ASSIGN(
      std::vector<TypeParam> m_coeffs_q0,
      this->SampleCoeffs(log_n, mod_params_q0, static_cast<Uint64>(t)));
  ASSERT_OK_AND_ASSIGN(
      auto m, RnsPolynomial<TypeParam>::ConvertBalancedFromPolynomialCoeffs(
                  m_coeffs_q0, mod_params_q0, main_moduli));
  ASSERT_OK_AND_ASSIGN(std::vector<TypeParam> e_coeffs_q0,
                       this->SampleTernaryCoeffs(log_n, mod_params_q0));
  ASSERT_OK_AND_ASSIGN(
      auto e, RnsPolynomial<TypeParam>::ConvertBalancedFromPolynomialCoeffs(
                  e_coeffs_q0, mod_params_q0, main_moduli));
  ASSERT_OK_AND_ASSIGN(RnsPolynomial<TypeParam> a, e.Mul(t, main_moduli));
  ASSERT_OK(a.AddInPlace(m, main_moduli));  // a = m + t * e (mod Q).

  // a = P * (m + t * e) (mod Q).
  int num_main_moduli = main_moduli.size();
  auto p_mod_qs =
      this->rns_context_->AuxModulusResidues().subspan(0, num_main_moduli);
  ASSERT_OK(a.MulInPlace(p_mod_qs, main_moduli));  // c = a*P (mod Q)

  // Add additional error t*f (mod Q).
  ASSERT_OK_AND_ASSIGN(std::vector<TypeParam> f_coeffs_q0,
                       this->SampleTernaryCoeffs(log_n, mod_params_q0));
  ASSERT_OK_AND_ASSIGN(
      auto f, RnsPolynomial<TypeParam>::ConvertBalancedFromPolynomialCoeffs(
                  f_coeffs_q0, mod_params_q0, main_moduli));
  ASSERT_OK(f.MulInPlace(t, main_moduli));
  ASSERT_OK(a.AddInPlace(f, main_moduli));

  // Let b = t * f (mod P).
  int level = num_main_moduli - 1;
  ASSERT_OK_AND_ASSIGN(std::vector<TypeParam> q_hat_inv_mod_qs,
                       this->rns_context_->MainPrimeModulusCrtFactors(level));
  ASSERT_OK_AND_ASSIGN(
      std::vector<RnsInt<TypeParam>> q_hat_mod_ps,
      this->rns_context_->MainPrimeModulusComplementResidues(level));
  ASSERT_OK_AND_ASSIGN(
      auto b, RnsPolynomial<TypeParam>::ConvertBalancedFromPolynomialCoeffs(
                  f_coeffs_q0, mod_params_q0, aux_moduli));
  ASSERT_OK(b.MulInPlace(t, aux_moduli));
  EXPECT_EQ(b.NumModuli(), aux_moduli.size());

  // Apply approximate modulus reduction.
  auto p_inv_mod_qs = this->rns_context_->AuxModulusInverseResidues().subspan(
      0, num_main_moduli);
  ASSERT_OK(a.ApproxModReduceLsb(
      b, main_moduli, aux_moduli, p_inv_mod_qs,
      this->rns_context_->AuxPrimeModulusCrtFactors(),
      this->rns_context_->AuxPrimeModulusComplementResidues(),
      this->rns_context_->PlaintextModulusInverseAuxResidues(),
      this->rns_context_->AuxModulusPlaintextResidue(),
      this->rns_context_->PlaintextModulus()));
  ASSERT_OK(a.ConvertToCoeffForm(main_moduli));
  ASSERT_EQ(a.NumModuli(), main_moduli.size());
  ASSERT_EQ(a.NumCoeffs(), 1 << log_n);

  // Compute CRT interpolation of a.
  std::vector<BigInteger> modulus_hats =
      RnsModulusComplements<TypeParam, BigInteger>(main_moduli);
  ASSERT_OK_AND_ASSIGN(std::vector<TypeParam> modulus_hat_invs,
                       this->rns_context_->MainPrimeModulusCrtFactors(level));
  ASSERT_OK_AND_ASSIGN(
      std::vector<BigInteger> a_coeffs,
      (CrtInterpolation<TypeParam, BigInteger>(
          a.Coeffs(), main_moduli, modulus_hats, modulus_hat_invs)));
  BigInteger q_half = q >> 1;
  BigInteger t_big = static_cast<BigInteger>(t);
  for (int i = 0; i < (1 << log_n); ++i) {
    bool is_negative = a_coeffs[i] > q_half;
    BigInteger x_abs = is_negative ? q - a_coeffs[i] : a_coeffs[i];
    BigInteger x_reminder = is_negative ? t_big - x_abs % t_big : x_abs % t_big;
    x_reminder = x_reminder % t_big;  // reduce t to 0
    BigInteger x_expected =
        static_cast<BigInteger>(m_coeffs_q0[i].ExportInt(mod_params_q0));
    EXPECT_EQ(x_reminder, x_expected);
  }
}

TYPED_TEST(RnsPolynomialTest, ApproxModReduceMsbFailsIfWrongNumberOfModuli) {
  int num_moduli = this->moduli_.size();
  ASSERT_OK_AND_ASSIGN(RnsPolynomial<TypeParam> a,
                       RnsPolynomial<TypeParam>::CreateZero(
                           this->rns_context_->LogN(), this->moduli_));
  std::vector<const PrimeModulus<TypeParam>*> aux_moduli =
      this->rns_context_->AuxPrimeModuli();
  ASSERT_OK_AND_ASSIGN(RnsPolynomial<TypeParam> b,
                       RnsPolynomial<TypeParam>::CreateZero(
                           this->rns_context_->LogN(), aux_moduli));
  EXPECT_THAT(a.ApproxModReduceMsb(
                  b, absl::MakeSpan(this->moduli_).subspan(0, num_moduli - 1),
                  aux_moduli, this->rns_context_->AuxModulusInverseResidues(),
                  this->rns_context_->AuxPrimeModulusCrtFactors(),
                  this->rns_context_->AuxPrimeModulusComplementResidues()),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr(absl::StrCat("`this_moduli` must contain ",
                                              num_moduli, " RNS moduli"))));
}

TYPED_TEST(RnsPolynomialTest, ApproxModReduceMsbFailsIfWrongNumberOfAuxModuli) {
  std::vector<const PrimeModulus<TypeParam>*> aux_moduli =
      this->rns_context_->AuxPrimeModuli();
  if (aux_moduli.size() <= 1) {
    // This test requires at least two auxiliary moduli.
    return;
  }

  ASSERT_OK_AND_ASSIGN(RnsPolynomial<TypeParam> a,
                       RnsPolynomial<TypeParam>::CreateZero(
                           this->rns_context_->LogN(), this->moduli_));
  ASSERT_OK_AND_ASSIGN(RnsPolynomial<TypeParam> b,
                       RnsPolynomial<TypeParam>::CreateZero(
                           this->rns_context_->LogN(), aux_moduli));
  EXPECT_THAT(
      a.ApproxModReduceMsb(
          b, this->moduli_, absl::MakeSpan(aux_moduli).subspan(0, 1),
          this->rns_context_->AuxModulusInverseResidues(),
          this->rns_context_->AuxPrimeModulusCrtFactors(),
          this->rns_context_->AuxPrimeModulusComplementResidues()),
      StatusIs(absl::StatusCode::kInvalidArgument,
               HasSubstr(absl::StrCat("`aux_moduli` must contain ",
                                      aux_moduli.size(), " RNS moduli"))));
}

TYPED_TEST(RnsPolynomialTest, ApproxModReduceMsbFailsIfWrongNumberOfConstants) {
  int num_moduli = this->moduli_.size();
  auto aux_moduli = this->rns_context_->AuxPrimeModuli();
  int num_aux_moduli = aux_moduli.size();
  if (num_moduli <= 1 || num_aux_moduli <= 1) {
    // This test requires >= 2main moduli and >= 2 aux moduli.
    return;
  }

  ASSERT_OK_AND_ASSIGN(RnsPolynomial<TypeParam> a,
                       RnsPolynomial<TypeParam>::CreateZero(
                           this->rns_context_->LogN(), this->moduli_));
  ASSERT_OK_AND_ASSIGN(RnsPolynomial<TypeParam> b,
                       RnsPolynomial<TypeParam>::CreateZero(
                           this->rns_context_->LogN(), aux_moduli));
  EXPECT_THAT(a.ApproxModReduceMsb(
                  b, this->moduli_, aux_moduli,
                  this->rns_context_->AuxModulusInverseResidues().subspan(
                      0, num_moduli - 1),
                  this->rns_context_->AuxPrimeModulusCrtFactors(),
                  this->rns_context_->AuxPrimeModulusComplementResidues()),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr(absl::StrCat("`p_inv_mod_qs` must contain ",
                                              num_moduli, " elements"))));
  EXPECT_THAT(
      a.ApproxModReduceMsb(
          b, this->moduli_, aux_moduli,
          this->rns_context_->AuxModulusInverseResidues(),
          this->rns_context_->AuxPrimeModulusCrtFactors().subspan(
              0, num_aux_moduli - 1),
          this->rns_context_->AuxPrimeModulusComplementResidues()),
      StatusIs(absl::StatusCode::kInvalidArgument,
               HasSubstr(absl::StrCat("`prime_p_hat_inv_mod_ps` must contain ",
                                      num_aux_moduli, " elements"))));

  EXPECT_THAT(
      a.ApproxModReduceMsb(
          b, this->moduli_, aux_moduli,
          this->rns_context_->AuxModulusInverseResidues(),
          this->rns_context_->AuxPrimeModulusCrtFactors(),
          this->rns_context_->AuxPrimeModulusComplementResidues().subspan(
              0, num_aux_moduli - 1)),
      StatusIs(absl::StatusCode::kInvalidArgument,
               HasSubstr(absl::StrCat("`prime_p_hat_mod_qs` must contain ",
                                      num_aux_moduli, " elements"))));
}

// Checks that invoking `ApproxModReduceMsb()` computes an approximate version
// of round(a / P) where a (mod Q*P). This makes sure that plaintext messages
// are kept at MSB of modulus after `ApproxModReduceMsb`.
TYPED_TEST(RnsPolynomialTest, ApproxModReduceMsbIsCorrect) {
  using Integer = typename TypeParam::Int;
  using BigInteger = typename testing::CompositeModulus<Integer>::value_type;

  // It's easier to sample error coefficients smaller than 2^32.
  constexpr size_t k_max_error_bits =
      std::min(sizeof(Integer), sizeof(Uint64)) / 2;
  constexpr Uint64 k_max_error = 1 << k_max_error_bits;

  // For correctness we need the auxiliary modulus P > the main modulus Q. So we
  // define Q using a subset of q_i's whose product is smaller than P.
  std::vector<const PrimeModulus<TypeParam>*> aux_moduli =
      this->rns_context_->AuxPrimeModuli();
  BigInteger p{1};
  for (auto prime_modulus : aux_moduli) {
    p *= static_cast<BigInteger>(prime_modulus->Modulus());
  }
  std::vector<const PrimeModulus<TypeParam>*> main_moduli;
  BigInteger q{1};
  for (int i = 0; i < this->moduli_.size(); ++i) {
    BigInteger qq = q * static_cast<BigInteger>(this->moduli_[i]->Modulus());
    if (qq < p) {
      q = qq;
      main_moduli.push_back(this->moduli_[i]);
    } else {
      break;
    }
  }
  ASSERT_GE(main_moduli.size(), 1);

  ASSERT_OK_AND_ASSIGN(RnsPolynomial<TypeParam> a,
                       this->SampleRnsPolynomial(main_moduli));

  // Convert a from mod-Q to mod-P (result in b). The pair (a, b) now is a RNS
  // polynomial modulo Q*P that approximately represents a (mod Q).
  int num_coeffs = 1 << this->rns_context_->LogN();
  int num_main_moduli = main_moduli.size();
  int level = num_main_moduli - 1;
  ASSERT_OK_AND_ASSIGN(std::vector<TypeParam> q_hat_inv_mod_qs,
                       this->rns_context_->MainPrimeModulusCrtFactors(level));
  ASSERT_OK_AND_ASSIGN(
      std::vector<RnsInt<TypeParam>> q_hat_mod_ps,
      this->rns_context_->MainPrimeModulusComplementResidues(level));

  // We check the correctness of ApproxModReduceMsb using these steps:
  // 1. Scale up a by P (mod Q), get c (mod Q);
  // 2. Add small e to c, and let b = e (mod P);
  // 3. Then compute d = c.ApproxModReduceMsb(b).
  // Then we should have d \approx a with small error. So, if a stores plaintext
  // messages in high order bits, then d should preserve such plaintext bits.
  //
  // The rational is as follows: If a' = a.ApproxSwitchRnsBasis(Q, P), then we
  // have (a, a') mod (Q, P) represents a + u*Q (mod Q*P) for small scalar u due
  // to the property of ApproxSwitchRnsBasis.
  // Now, if c = a*P + e (mod Q) and b = a'*P + e (mod P) = e (mod P), then
  // (c, b) mod (Q, P) represents a*P + u*Q*P + e (mod Q*P) = a*P + e (mod Q*P).
  // Since c.ApproxModReduceMsb(e) computes round((a*P + e) / P) (mod Q), and
  // since e is small, the result d should be close to a.
  RnsPolynomial<TypeParam> c = a;
  auto p_mod_qs =
      this->rns_context_->AuxModulusResidues().subspan(0, num_main_moduli);
  auto mod_params_q0 = main_moduli[0]->ModParams();
  ASSERT_OK(c.MulInPlace(p_mod_qs, main_moduli));  // c = a*P (mod Q)
  ASSERT_OK_AND_ASSIGN(std::vector<TypeParam> e_coeffs_q0,
                       this->SampleCoeffs(this->rns_context_->LogN(),
                                          mod_params_q0, k_max_error));
  ASSERT_OK_AND_ASSIGN(
      auto e, RnsPolynomial<TypeParam>::ConvertFromPolynomialCoeffs(
                  e_coeffs_q0, mod_params_q0, main_moduli));  // e (mod Q)
  ASSERT_OK(c.AddInPlace(e, main_moduli));  // c = a*P + e (mod Q)
  ASSERT_OK_AND_ASSIGN(
      auto b, RnsPolynomial<TypeParam>::ConvertFromPolynomialCoeffs(
                  e_coeffs_q0, mod_params_q0, aux_moduli));  // b = e (mod P)

  // Apply approximate modulus reduction.
  auto p_inv_mod_qs = this->rns_context_->AuxModulusInverseResidues().subspan(
      0, num_main_moduli);
  ASSERT_OK(c.ApproxModReduceMsb(
      b, main_moduli, aux_moduli, p_inv_mod_qs,
      this->rns_context_->AuxPrimeModulusCrtFactors(),
      this->rns_context_->AuxPrimeModulusComplementResidues()));
  ASSERT_EQ(c.NumModuli(), num_main_moduli);
  ASSERT_EQ(c.NumCoeffs(), num_coeffs);

  // Convert a, b, c to Coefficient form.
  ASSERT_OK(a.ConvertToCoeffForm(main_moduli));
  ASSERT_OK(c.ConvertToCoeffForm(main_moduli));

  // Compute CRT interpolations of a, c (mod Q) and d (mod Q*P).
  std::vector<BigInteger> modulus_hats =
      RnsModulusComplements<TypeParam, BigInteger>(main_moduli);
  ASSERT_OK_AND_ASSIGN(std::vector<TypeParam> modulus_hat_invs,
                       this->rns_context_->MainPrimeModulusCrtFactors(level));

  ASSERT_OK_AND_ASSIGN(
      std::vector<BigInteger> a_coeffs,
      (CrtInterpolation<TypeParam, BigInteger>(
          a.Coeffs(), main_moduli, modulus_hats, modulus_hat_invs)));
  ASSERT_OK_AND_ASSIGN(
      std::vector<BigInteger> c_coeffs,
      (CrtInterpolation<TypeParam, BigInteger>(
          c.Coeffs(), main_moduli, modulus_hats, modulus_hat_invs)));
  // Each coefficient of approximately mod reduced polynomial c(X) should be
  // within the range of coefficient of a(X) +/- number of main moduli.
  for (int i = 0; i < num_coeffs; ++i) {
    BigInteger x = a_coeffs[i];
    BigInteger y = c_coeffs[i];
    EXPECT_LE(x, y + num_main_moduli);
    EXPECT_LE(y, x + num_main_moduli);
  }
}

TYPED_TEST(RnsPolynomialTest, AppendCoeffVector) {
  int num_coeffs = 1 << this->rns_context_->LogN();
  int num_existing_moduli = this->moduli_.size();

  // Sample a RNS polynomial mod Q.
  ASSERT_OK_AND_ASSIGN(RnsPolynomial<TypeParam> a, this->SampleRnsPolynomial());
  ASSERT_EQ(a.NumModuli(), num_existing_moduli);
  ASSERT_EQ(a.NumCoeffs(), num_coeffs);

  // Append a coefficient vector of 1 (mod p_0) to a. Now a is modulo Q * p_0.
  std::vector<const PrimeModulus<TypeParam>*> aux_moduli =
      this->rns_context_->AuxPrimeModuli();
  ASSERT_FALSE(aux_moduli.empty());
  TypeParam zero_mod_p = TypeParam::ImportZero(aux_moduli[0]->ModParams());
  std::vector<TypeParam> coeffs_p(num_coeffs, zero_mod_p);
  RnsPolynomial<TypeParam> b = a;
  b.AppendCoeffVector(std::move(coeffs_p));
  EXPECT_EQ(b.NumModuli(), num_existing_moduli + 1);
  EXPECT_EQ(b.NumCoeffs(), num_coeffs);

  // The new coefficient vector should be the last in a.Coeffs(), and the
  // first `num_existing_moduli` coefficient vectors should remain the same.
  for (int i = 0; i < num_existing_moduli; ++i) {
    EXPECT_EQ(b.Coeffs()[i], a.Coeffs()[i]);
  }
  for (int j = 0; j < num_coeffs; ++j) {
    EXPECT_EQ(b.Coeffs()[num_existing_moduli][j], zero_mod_p);
  }
}

TYPED_TEST(RnsPolynomialTest, DetachLastCoeffVector) {
  int num_existing_moduli = this->moduli_.size();
  ASSERT_GE(num_existing_moduli, 1);

  // Sample a RNS polynomial mod Q, where the last prime moduli of Q is q_l.
  ASSERT_OK_AND_ASSIGN(RnsPolynomial<TypeParam> a, this->SampleRnsPolynomial());
  ASSERT_EQ(a.NumModuli(), num_existing_moduli);

  // Make a copy of a.
  RnsPolynomial<TypeParam> a_original = a;

  // Detach the last sub-polynomial, i.e. a (mod q_l).
  ASSERT_OK_AND_ASSIGN(std::vector<TypeParam> coeffs_ql,
                       a.DetachLastCoeffVector());
  ASSERT_EQ(a.NumModuli(), num_existing_moduli - 1);
  for (int i = 0; i < num_existing_moduli - 1; ++i) {
    EXPECT_EQ(a.Coeffs()[i], a_original.Coeffs()[i]);
  }
  EXPECT_EQ(coeffs_ql, a_original.Coeffs()[num_existing_moduli - 1]);
}

TYPED_TEST(RnsPolynomialTest, DetachLastCoeffVectorFailsIfEmptyCoeffVectors) {
  // Sample a RNS polynomial mod q_0.
  ASSERT_GE(this->moduli_.size(), 1);
  ASSERT_OK_AND_ASSIGN(RnsPolynomial<TypeParam> a,
                       RnsPolynomial<TypeParam>::CreateZero(
                           this->rns_context_->LogN(), {this->moduli_[0]}));
  // Make `a` an empty polynomial.
  ASSERT_OK(a.DetachLastCoeffVector());
  ASSERT_EQ(a.NumModuli(), 0);
  EXPECT_THAT(a.DetachLastCoeffVector(),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("RnsPolynomial is empty")));
}

TYPED_TEST(RnsPolynomialTest, ReplaceSubPolynomialAt) {
  int num_coeffs = 1 << this->rns_context_->LogN();
  int num_existing_moduli = this->moduli_.size();

  // Sample a RNS polynomial mod Q.
  ASSERT_OK_AND_ASSIGN(RnsPolynomial<TypeParam> a, this->SampleRnsPolynomial());
  ASSERT_EQ(a.NumModuli(), num_existing_moduli);
  ASSERT_EQ(a.NumCoeffs(), num_coeffs);

  // Make a copy.
  RnsPolynomial<TypeParam> a_original = a;

  // Replace the first sub-polynomial (i.e. a (mod q_0)) with new
  // coefficients.
  TypeParam zero_mod_q0 = TypeParam::ImportZero(this->moduli_[0]->ModParams());
  std::vector<TypeParam> zero_coeffs_q0(num_coeffs, zero_mod_q0);
  Polynomial<TypeParam> zero_q0(std::move(zero_coeffs_q0));
  ASSERT_OK(a.ReplaceSubPolynomialAt(0, zero_q0));
  EXPECT_EQ(a.NumModuli(), num_existing_moduli);  // should remain the same
  EXPECT_EQ(a.NumCoeffs(), num_coeffs);

  for (int j = 0; j < num_coeffs; ++j) {
    EXPECT_EQ(a.Coeffs()[0][j], zero_mod_q0);
  }
  for (int i = 1; i < num_existing_moduli; ++i) {
    EXPECT_EQ(a.Coeffs()[i], a_original.Coeffs()[i]);
  }
}

TYPED_TEST(RnsPolynomialTest, ReplaceSubPolynomialAtFailsIfIndexOutOfRange) {
  int num_coeffs = 1 << this->rns_context_->LogN();
  int num_existing_moduli = this->moduli_.size();
  ASSERT_OK_AND_ASSIGN(RnsPolynomial<TypeParam> a,
                       RnsPolynomial<TypeParam>::CreateZero(
                           this->rns_context_->LogN(), this->moduli_));
  ASSERT_EQ(a.NumModuli(), num_existing_moduli);

  // Replace the first sub-polynomial (i.e. a (mod q_0)) with new
  // coefficients.
  TypeParam zero_mod_q0 = TypeParam::ImportZero(this->moduli_[0]->ModParams());
  std::vector<TypeParam> zero_coeffs_q0(num_coeffs, zero_mod_q0);
  Polynomial<TypeParam> zero_q0(std::move(zero_coeffs_q0));

  EXPECT_THAT(
      a.ReplaceSubPolynomialAt(-1, zero_q0),
      StatusIs(absl::StatusCode::kInvalidArgument, HasSubstr("out of bound")));
  EXPECT_THAT(
      a.ReplaceSubPolynomialAt(num_existing_moduli, zero_q0),
      StatusIs(absl::StatusCode::kInvalidArgument, HasSubstr("out of bound")));
}

TYPED_TEST(RnsPolynomialTest, DeserializeFailsIfLogNIsNonPositive) {
  ASSERT_OK_AND_ASSIGN(RnsPolynomial<TypeParam> a,
                       RnsPolynomial<TypeParam>::CreateZero(
                           this->rns_context_->LogN(), this->moduli_));
  ASSERT_OK_AND_ASSIGN(SerializedRnsPolynomial serialized,
                       a.Serialize(this->moduli_));
  serialized.set_log_n(0);
  EXPECT_THAT(RnsPolynomial<TypeParam>::Deserialize(serialized, this->moduli_),
              StatusIs(absl::StatusCode::kInvalidArgument, HasSubstr("log_n")));
  serialized.set_log_n(-1);
  EXPECT_THAT(RnsPolynomial<TypeParam>::Deserialize(serialized, this->moduli_),
              StatusIs(absl::StatusCode::kInvalidArgument, HasSubstr("log_n")));
}

TYPED_TEST(RnsPolynomialTest, DeserializeFailsIfIncorrectNumberOfCoeffVectors) {
  ASSERT_OK_AND_ASSIGN(RnsPolynomial<TypeParam> a,
                       RnsPolynomial<TypeParam>::CreateZero(
                           this->rns_context_->LogN(), this->moduli_));
  ASSERT_OK_AND_ASSIGN(SerializedRnsPolynomial serialized,
                       a.Serialize(this->moduli_));
  // Append a coefficient vector to the serialized object.
  serialized.add_coeff_vectors("");
  EXPECT_THAT(RnsPolynomial<TypeParam>::Deserialize(serialized, this->moduli_),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("Number of serialized coefficient vectors")));
}

TYPED_TEST(RnsPolynomialTest, SerializeFailsIfIncorrectNumberOfModuli) {
  ASSERT_OK_AND_ASSIGN(RnsPolynomial<TypeParam> a,
                       RnsPolynomial<TypeParam>::CreateZero(
                           this->rns_context_->LogN(), this->moduli_));
  int num_moduli = this->moduli_.size();
  EXPECT_THAT(
      a.Serialize(absl::MakeSpan(this->moduli_).subspan(0, num_moduli - 1)),
      StatusIs(absl::StatusCode::kInvalidArgument, HasSubstr("moduli")));
}

TYPED_TEST(RnsPolynomialTest, SerializedDeserializes) {
  ASSERT_OK_AND_ASSIGN(RnsPolynomial<TypeParam> a, this->SampleRnsPolynomial());
  ASSERT_OK_AND_ASSIGN(SerializedRnsPolynomial serialized1,
                       a.Serialize(this->moduli_));
  EXPECT_EQ(serialized1.log_n(), a.LogN());
  EXPECT_EQ(serialized1.is_ntt(), a.IsNttForm());
  EXPECT_EQ(serialized1.coeff_vectors_size(), this->moduli_.size());

  ASSERT_OK_AND_ASSIGN(
      RnsPolynomial<TypeParam> deserialized1,
      RnsPolynomial<TypeParam>::Deserialize(serialized1, this->moduli_));
  EXPECT_EQ(deserialized1, a);

  // Change the format (Coefficient <-> NTT forms) and serialize again.
  if (a.IsNttForm()) {
    ASSERT_OK(a.ConvertToCoeffForm(this->moduli_));
  } else {
    ASSERT_OK(a.ConvertToNttForm(this->moduli_));
  }
  ASSERT_OK_AND_ASSIGN(SerializedRnsPolynomial serialized2,
                       a.Serialize(this->moduli_));
  ASSERT_OK_AND_ASSIGN(
      RnsPolynomial<TypeParam> deserialized2,
      RnsPolynomial<TypeParam>::Deserialize(serialized2, this->moduli_));
  EXPECT_EQ(deserialized2, a);
}

}  // namespace
}  // namespace rlwe
