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

#include "shell_encryption/dft_transformations.h"

#include <complex>
#include <memory>
#include <vector>

#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include "absl/status/status.h"
#include "absl/strings/str_cat.h"
#include "shell_encryption/context.h"
#include "shell_encryption/status_macros.h"
#include "shell_encryption/statusor.h"
#include "shell_encryption/testing/parameters.h"
#include "shell_encryption/testing/status_matchers.h"
#include "shell_encryption/testing/status_testing.h"
#include "shell_encryption/testing/testing_prng.h"

namespace rlwe {
namespace {

using ::rlwe::testing::StatusIs;
using ::testing::HasSubstr;

template <typename ModularInt>
class DftTransformationsTest : public ::testing::Test {
  using ModularIntParams = typename ModularInt::Params;

 protected:
  void SetUp() override {
    typename RlweContext<ModularInt>::Parameters params =
        testing::ContextParameters<ModularInt>::Value()[0];
    ASSERT_OK_AND_ASSIGN(context_, RlweContext<ModularInt>::Create(params));
    prng_ = std::make_unique<testing::TestingPrng>(0);
  }

  StatusOr<std::vector<ModularInt>> SampleCoeffs(
      int log_n, const ModularIntParams* mod_params) const {
    std::vector<ModularInt> coeffs;
    for (int i = 0; i < (1 << log_n); ++i) {
      RLWE_ASSIGN_OR_RETURN(ModularInt coeff,
                            ModularInt::ImportRandom(prng_.get(), mod_params));
      coeffs.push_back(std::move(coeff));
    }
    return coeffs;
  }

  std::unique_ptr<const RlweContext<ModularInt>> context_;
  std::unique_ptr<testing::TestingPrng> prng_;
};

TYPED_TEST_SUITE(DftTransformationsTest, testing::ModularIntTypes);

TYPED_TEST(DftTransformationsTest, ForwardNttFailsIfInvalidCoeffsLength) {
  // Create a vector of N+1 modular integers, where N is a power of two.
  int log_n = this->context_->GetLogN();
  std::vector<TypeParam> coeffs(
      (1 << log_n) + 1,
      TypeParam::ImportZero(this->context_->GetModulusParams()));

  EXPECT_THAT(ForwardNumberTheoreticTransform<TypeParam>(
                  coeffs, *this->context_->GetNttParams(),
                  *this->context_->GetModulusParams()),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr(absl::StrCat("`coeffs` must have contain ",
                                              (1 << log_n), " elements"))));
}

TYPED_TEST(DftTransformationsTest, ForwardNumberTheoreticTransform) {
  for (const auto& params : testing::ContextParameters<TypeParam>::Value()) {
    ASSERT_OK_AND_ASSIGN(auto context, RlweContext<TypeParam>::Create(params));
    int log_n = context->GetLogN();
    ASSERT_OK_AND_ASSIGN(TypeParam a,
                         TypeParam::ImportRandom(this->prng_.get(),
                                                 context->GetModulusParams()));

    // If coeffs = {a, 0, 0, ...}, then NTT(coeffs) = {a, a, a, ...}.
    std::vector<TypeParam> coeffs(
        (1 << log_n), TypeParam::ImportZero(context->GetModulusParams()));
    coeffs[0] = a;

    ASSERT_OK(ForwardNumberTheoreticTransform<TypeParam>(
        coeffs, *context->GetNttParams(), *context->GetModulusParams()));
    EXPECT_EQ(coeffs.size(), (1 << log_n));
    for (auto const& coeff : coeffs) {
      EXPECT_EQ(coeff, a);
    }
  }
}

TYPED_TEST(DftTransformationsTest, InverseNttFailsIfInvalidCoeffsLength) {
  // Create a vector of N+1 modular integers, where N is a power of two.
  int log_n = this->context_->GetLogN();
  std::vector<TypeParam> coeffs(
      (1 << log_n) + 1,
      TypeParam::ImportZero(this->context_->GetModulusParams()));

  EXPECT_THAT(InverseNumberTheoreticTransform<TypeParam>(
                  coeffs, *this->context_->GetNttParams(),
                  *this->context_->GetModulusParams()),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr(absl::StrCat("`coeffs` must have contain ",
                                              (1 << log_n), " elements"))));
}

TYPED_TEST(DftTransformationsTest, InverseNumberTheoreticTransform) {
  for (const auto& params : testing::ContextParameters<TypeParam>::Value()) {
    ASSERT_OK_AND_ASSIGN(auto context, RlweContext<TypeParam>::Create(params));
    int log_n = context->GetLogN();
    ASSERT_OK_AND_ASSIGN(TypeParam a,
                         TypeParam::ImportRandom(this->prng_.get(),
                                                 context->GetModulusParams()));

    // If NTT(coeffs) = {a, a, a, ...}, then coeffs = {a, 0, 0, ...}.
    std::vector<TypeParam> coeffs((1 << log_n), a);
    ASSERT_OK(InverseNumberTheoreticTransform<TypeParam>(
        coeffs, *context->GetNttParams(), *context->GetModulusParams()));
    EXPECT_EQ(coeffs.size(), (1 << log_n));
    EXPECT_EQ(coeffs[0], a);
    TypeParam zero = TypeParam::ImportZero(context->GetModulusParams());
    for (int i = 1; i < coeffs.size(); ++i) {
      EXPECT_EQ(coeffs[i], zero);
    }
  }
}

TEST(DftTransformations, IterativeHalfCooleyTukeyFailsIfNotEnoughPsis) {
  constexpr int len = 1 << 5;
  std::vector<std::complex<double>> values(len, {0, 0});
  std::vector<std::complex<double>> too_short_psis = {{1, 0}};
  EXPECT_THAT(
      IterativeHalfCooleyTukey(values, too_short_psis),
      StatusIs(absl::StatusCode::kInvalidArgument,
               HasSubstr("Not enough primitive roots in `psis_bitrev`")));
}

TEST(DftTransformations, IterativeHalfCooleyTukeyFailsIfWrongInputLength) {
  constexpr int incorrect_len = 7;  // not a power of 2.
  std::vector<std::complex<double>> values(incorrect_len, {0, 0});
  std::vector<std::complex<double>> psis(incorrect_len, {1, 0});
  EXPECT_THAT(IterativeHalfCooleyTukey(values, psis),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("size of `coeffs` must be a power of two")));
}

TEST(DftTransformations, IterativeHalfGentlemanSandeFailsIfNotEnoughPsis) {
  constexpr int len = 1 << 5;
  std::vector<std::complex<double>> values(len, {0, 0});
  std::vector<std::complex<double>> too_short_psis = {{1, 0}};
  EXPECT_THAT(
      IterativeHalfGentlemanSande(values, too_short_psis),
      StatusIs(absl::StatusCode::kInvalidArgument,
               HasSubstr("Not enough primitive roots in `psis_bitrev_inv`")));
}

TEST(DftTransformations, IterativeHalfGentlemanSandeFailsIfWrongInputLength) {
  constexpr int incorrect_len = 7;  // not a power of 2.
  std::vector<std::complex<double>> values(incorrect_len, {0, 0});
  std::vector<std::complex<double>> psis_inv(incorrect_len, {1, 0});
  EXPECT_THAT(IterativeHalfGentlemanSande(values, psis_inv),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("size of `coeffs` must be a power of two")));
}

}  // namespace
}  // namespace rlwe
