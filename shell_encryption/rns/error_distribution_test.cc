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

#include "shell_encryption/rns/error_distribution.h"

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <memory>
#include <utility>
#include <vector>

#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include "absl/log/check.h"
#include "absl/status/status.h"
#include "absl/strings/string_view.h"
#include "shell_encryption/prng/single_thread_hkdf_prng.h"
#include "shell_encryption/rns/rns_context.h"
#include "shell_encryption/rns/rns_modulus.h"
#include "shell_encryption/rns/rns_polynomial.h"
#include "shell_encryption/rns/testing/parameters.h"
#include "shell_encryption/sampler/discrete_gaussian.h"
#include "shell_encryption/testing/parameters.h"
#include "shell_encryption/testing/status_matchers.h"
#include "shell_encryption/testing/status_testing.h"

namespace rlwe {
namespace {

using Prng = SingleThreadHkdfPrng;
using ::rlwe::testing::StatusIs;
using ::testing::HasSubstr;

constexpr int kVariance = 8;             // Common error variance.
constexpr double kBaseGaussianS = 12.8;  // For the base Gaussian distribution.
constexpr absl::string_view kPrngSeed =
    "0123456789abcdef0123456789abcdef0123456789abcdef0123456789abcdef";

// Test fixture for error distribution sampling functions.
template <typename ModularInt>
class ErrorDistributionTest : public ::testing::Test {
 protected:
  using ModularIntParams = typename ModularInt::Params;
  using Integer = typename ModularInt::Int;

  void SetUp() override {
    testing::RnsParameters<ModularInt> rns_params =
        testing::GetRnsParameters<ModularInt>();

    // Create the RNS context.
    auto rns_context = RnsContext<ModularInt>::Create(
        rns_params.log_n, rns_params.qs, rns_params.ps, rns_params.t);
    CHECK(rns_context.ok());
    rns_context_ = std::make_unique<const RnsContext<ModularInt>>(
        std::move(rns_context.value()));
    moduli_ = rns_context_->MainPrimeModuli();

    // Create a PRNG for instantiating the samplers.
    auto prng = Prng::Create(kPrngSeed.substr(0, Prng::SeedLength()));
    CHECK(prng.ok());
    prng_ = std::move(prng.value());
  }

  void SetUpGaussianSampler() {
    auto sampler = DiscreteGaussianSampler<Integer>::Create(kBaseGaussianS);
    CHECK(sampler.ok());
    dg_sampler_ = std::move(sampler.value());
  }

  // Convert a vector `coeffs` of balanced mod q integers to their signed
  // integer values, where q is given in `mod_params`. Balanced mod q
  // representation uses numbers in [0, q / 2) for positive values and numbers
  // in [q / 2, q) for negative values. We assume that `coeffs` represent
  // integers with small absolute values to fit in int32_t.
  std::vector<int32_t> ConvertToSignedIntegerValues(
      const std::vector<ModularInt>& coeffs,
      const ModularIntParams* mod_params) const {
    Integer q = mod_params->modulus;
    Integer q_half = q >> 1;
    std::vector<int32_t> values;
    values.reserve(coeffs.size());
    for (const auto& coeff : coeffs) {
      Integer coeff_bal_mod_q = coeff.ExportInt(mod_params);
      if (coeff_bal_mod_q > q_half) {  // negative value
        int32_t coeff_abs_value = static_cast<int32_t>(q - coeff_bal_mod_q);
        values.push_back(-coeff_abs_value);
      } else {  // positive value
        int32_t coeff_abs_value = static_cast<int32_t>(coeff_bal_mod_q);
        values.push_back(coeff_abs_value);
      }
    }
    return values;
  }

  std::unique_ptr<const RnsContext<ModularInt>> rns_context_;
  std::vector<const PrimeModulus<ModularInt>*> moduli_;
  std::unique_ptr<Prng> prng_;
  std::unique_ptr<const DiscreteGaussianSampler<Integer>> dg_sampler_;
};
TYPED_TEST_SUITE(ErrorDistributionTest, testing::ModularIntTypes);

// Unit tests for `SampleError()`.

TYPED_TEST(ErrorDistributionTest, SampleErrorFailsIfLogNIsNotPositive) {
  EXPECT_THAT(SampleError<TypeParam>(
                  /*log_n=*/-1, kVariance, this->moduli_, this->prng_.get()),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`log_n` must be positive")));
  EXPECT_THAT(SampleError<TypeParam>(
                  /*log_n=*/0, kVariance, this->moduli_, this->prng_.get()),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`log_n` must be positive")));
}

TYPED_TEST(ErrorDistributionTest, SampleErrorFailsIfVarianceIsNotPositive) {
  int log_n = this->rns_context_->LogN();
  EXPECT_THAT(SampleError<TypeParam>(log_n, /*variance=*/-1, this->moduli_,
                                     this->prng_.get()),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`variance` must be positive")));
  EXPECT_THAT(SampleError<TypeParam>(log_n, /*variance=*/0, this->moduli_,
                                     this->prng_.get()),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`variance` must be positive")));
}

TYPED_TEST(ErrorDistributionTest, SampleErrorFailsIfModuliIsEmpty) {
  int log_n = this->rns_context_->LogN();
  EXPECT_THAT(SampleError<TypeParam>(log_n, kVariance, /*moduli=*/{},
                                     this->prng_.get()),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`moduli` must not be empty")));
}

TYPED_TEST(ErrorDistributionTest,
           SampleErrorReturnsPolynomialWithSmallCoefficients) {
  int log_n = this->rns_context_->LogN();
  ASSERT_OK_AND_ASSIGN(RnsPolynomial<TypeParam> e,
                       SampleError<TypeParam>(log_n, kVariance, this->moduli_,
                                              this->prng_.get()));
  ASSERT_EQ(e.NumModuli(), this->moduli_.size());
  ASSERT_EQ(e.NumCoeffs(), 1 << log_n);

  // Since the error variance is sufficiently smaller than the RNS moduli qi,
  // e (mod qi) and e (mod qj) should represent the same integer values for all
  // other RNS moduli qj.
  ASSERT_OK(e.ConvertToCoeffForm(this->moduli_));
  std::vector<int32_t> e_coeffs_q0 = this->ConvertToSignedIntegerValues(
      e.Coeffs()[0], this->moduli_[0]->ModParams());
  for (int i = 1; i < this->moduli_.size(); ++i) {
    std::vector<int32_t> e_coeffs_qi = this->ConvertToSignedIntegerValues(
        e.Coeffs()[i], this->moduli_[i]->ModParams());
    EXPECT_EQ(e_coeffs_q0, e_coeffs_qi);
  }

  // The error should be small but not zero.
  int32_t e_coeff_max = 0;
  for (auto e_coeff_q0 : e_coeffs_q0) {
    int32_t e_coeff_abs = abs(e_coeff_q0);
    if (e_coeff_abs > e_coeff_max) {
      e_coeff_max = e_coeff_abs;
    }
    // Each error value should be in [-2 * variance, 2 * variance].
    EXPECT_LT(e_coeff_abs, 2 * kVariance + 1);
  }
  EXPECT_GT(e_coeff_max, 0);
}

// Tests for `SampleUniformTernary`.
TYPED_TEST(ErrorDistributionTest,
           SampleUniformTernaryFailsIfLogNIsNotPositive) {
  EXPECT_THAT(SampleUniformTernary<TypeParam>(
                  /*log_n=*/-1, this->moduli_, this->prng_.get()),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`log_n` must be positive")));
  EXPECT_THAT(SampleUniformTernary<TypeParam>(
                  /*log_n=*/0, this->moduli_, this->prng_.get()),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`log_n` must be positive")));
}

TYPED_TEST(ErrorDistributionTest, SampleUniformTernaryFailsIfModuliIsEmpty) {
  int log_n = this->rns_context_->LogN();
  EXPECT_THAT(
      SampleUniformTernary<TypeParam>(log_n, /*moduli=*/{}, this->prng_.get()),
      StatusIs(absl::StatusCode::kInvalidArgument,
               HasSubstr("`moduli` must not be empty")));
}

TYPED_TEST(ErrorDistributionTest,
           SampleUniformTernaryReturnsPolynomialWithTernaryCoefficients) {
  int log_n = this->rns_context_->LogN();
  ASSERT_OK_AND_ASSIGN(
      RnsPolynomial<TypeParam> e,
      SampleUniformTernary<TypeParam>(log_n, this->moduli_, this->prng_.get()));
  ASSERT_EQ(e.NumModuli(), this->moduli_.size());
  ASSERT_EQ(e.NumCoeffs(), 1 << log_n);
  ASSERT_TRUE(e.IsNttForm());

  // When in the coefficient form, e (mod qi) and e (mod qj) should represent
  // the same integer values in {-1, 0, 1} for all RNS moduli qi and qj.
  ASSERT_OK(e.ConvertToCoeffForm(this->moduli_));
  std::vector<int32_t> e_coeffs_q0 = this->ConvertToSignedIntegerValues(
      e.Coeffs()[0], this->moduli_[0]->ModParams());
  for (int k = 0; k < e_coeffs_q0.size(); ++k) {
    EXPECT_TRUE(e_coeffs_q0[k] <= 1 || e_coeffs_q0[k] == -1);
  }
  for (int i = 1; i < this->moduli_.size(); ++i) {
    std::vector<int32_t> e_coeffs_qi = this->ConvertToSignedIntegerValues(
        e.Coeffs()[i], this->moduli_[i]->ModParams());
    EXPECT_EQ(e_coeffs_q0, e_coeffs_qi);
  }
}

// Tests for `SampleDiscreteGaussian`.
TYPED_TEST(ErrorDistributionTest,
           SampleDiscreteGaussianFailsIfLogNIsNotPositive) {
  this->SetUpGaussianSampler();
  EXPECT_THAT(SampleDiscreteGaussian<TypeParam>(
                  /*log_n=*/-1, kBaseGaussianS, this->moduli_,
                  this->dg_sampler_.get(), this->prng_.get()),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`log_n` must be positive")));
  EXPECT_THAT(SampleDiscreteGaussian<TypeParam>(
                  /*log_n=*/0, kBaseGaussianS, this->moduli_,
                  this->dg_sampler_.get(), this->prng_.get()),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`log_n` must be positive")));
}

TYPED_TEST(ErrorDistributionTest, SampleDiscreteGaussianFailsIfModuliIsEmpty) {
  this->SetUpGaussianSampler();
  int log_n = this->rns_context_->LogN();
  EXPECT_THAT(SampleDiscreteGaussian<TypeParam>(
                  log_n, kBaseGaussianS, /*moduli=*/{}, this->dg_sampler_.get(),
                  this->prng_.get()),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`moduli` must not be empty")));
}

TYPED_TEST(ErrorDistributionTest, SampleDiscreteGaussianFailsIfSamplerIsNull) {
  int log_n = this->rns_context_->LogN();
  EXPECT_THAT(SampleDiscreteGaussian<TypeParam>(
                  log_n, kBaseGaussianS, this->moduli_, /*dg_sampler=*/nullptr,
                  this->prng_.get()),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`dg_sampler` must not be null")));
}

TYPED_TEST(ErrorDistributionTest,
           SampleDiscreteGaussianReturnsPolynomialWithBoundedCoefficients) {
  this->SetUpGaussianSampler();
  int log_n = this->rns_context_->LogN();
  ASSERT_OK_AND_ASSIGN(RnsPolynomial<TypeParam> e,
                       SampleDiscreteGaussian<TypeParam>(
                           log_n, kBaseGaussianS, this->moduli_,
                           this->dg_sampler_.get(), this->prng_.get()));
  ASSERT_EQ(e.NumModuli(), this->moduli_.size());
  ASSERT_EQ(e.NumCoeffs(), 1 << log_n);

  // We choose a small standard deviation s such that e (mod qi) and e (mod qj)
  // should represent the same integer values for all RNS moduli qi and qj.
  // Moreover, the coefficients should have bounded values.
  std::vector<int32_t> e_coeffs_q0 = this->ConvertToSignedIntegerValues(
      e.Coeffs()[0], this->moduli_[0]->ModParams());
  for (int i = 1; i < this->moduli_.size(); ++i) {
    std::vector<int32_t> e_coeffs_qi = this->ConvertToSignedIntegerValues(
        e.Coeffs()[i], this->moduli_[i]->ModParams());
    EXPECT_EQ(e_coeffs_q0, e_coeffs_qi);
  }

  auto [e_coeff_min, e_coeff_max] =
      std::minmax_element(e_coeffs_q0.begin(), e_coeffs_q0.end());
  int32_t coeff_max_bound = static_cast<int32_t>(
      DiscreteGaussianSampler<typename TypeParam::Int>::kTailBoundMultiplier *
      kBaseGaussianS);
  EXPECT_GT(*e_coeff_max, 0);
  EXPECT_LT(*e_coeff_max, coeff_max_bound);
  EXPECT_GT(*e_coeff_min, -coeff_max_bound);
  EXPECT_LT(*e_coeff_min, 0);
}

}  // namespace
}  // namespace rlwe
