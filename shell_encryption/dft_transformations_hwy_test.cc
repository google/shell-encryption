/*
 * Copyright 2025 Google LLC.
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

#include "shell_encryption/dft_transformations_hwy.h"

#include <memory>
#include <vector>

#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include "absl/status/status.h"
#include "absl/strings/str_cat.h"
#include "hwy/targets.h"
#include "shell_encryption/context.h"
#include "shell_encryption/dft_transformations.h"
#include "shell_encryption/status_macros.h"
#include "shell_encryption/statusor.h"
#include "shell_encryption/testing/parameters.h"
#include "shell_encryption/testing/status_matchers.h"
#include "shell_encryption/testing/status_testing.h"
#include "shell_encryption/testing/testing_prng.h"

#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "shell_encryption/dft_transformations_hwy_test.cc"
#include "hwy/foreach_target.h"  // IWYU pragma: keep
#include "hwy/highway.h"

#if HWY_ONCE || HWY_IDE

namespace rlwe::internal {
namespace {

using ::rlwe::testing::StatusIs;
using ::testing::HasSubstr;

template <typename ModularInt>
class DftTransformationsHwyTest : public ::testing::Test {
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

TYPED_TEST_SUITE(DftTransformationsHwyTest, testing::ModularIntTypes);

TYPED_TEST(DftTransformationsHwyTest, ForwardNttFailsIfInvalidCoeffsLength) {
  // Create a vector of N+1 modular integers, where N is a power of two.
  int log_n = this->context_->GetLogN();
  std::vector<TypeParam> coeffs(
      (1 << log_n) + 1,
      TypeParam::ImportZero(this->context_->GetModulusParams()));

  EXPECT_THAT(ForwardNumberTheoreticTransformHwy<TypeParam>(
                  coeffs, *this->context_->GetNttParams(),
                  *this->context_->GetModulusParams()),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr(absl::StrCat("`coeffs` must have contain ",
                                              (1 << log_n), " elements"))));
}

TYPED_TEST(DftTransformationsHwyTest, ForwardNumberTheoreticTransform) {
  auto context_params = testing::ContextParameters<TypeParam>::Value();
  // Add corner test cases where log degree is very small.
  context_params.push_back({/*.modulus =*/5, /*.log_n =*/1, /*.log_t =*/1,
                            /*.variance =*/8});
  context_params.push_back({/*.modulus =*/17, /*.log_n =*/2, /*.log_t =*/1,
                            /*.variance =*/8});
  for (const auto& params : context_params) {
    ASSERT_OK_AND_ASSIGN(auto context, RlweContext<TypeParam>::Create(params));
    int log_n = context->GetLogN();

    ASSERT_OK_AND_ASSIGN(
        std::vector<TypeParam> coeffs_original,
        this->SampleCoeffs(log_n, context->GetModulusParams()));
    std::vector<TypeParam> coeffs_truth = coeffs_original;
    ASSERT_OK(ForwardNumberTheoreticTransform<TypeParam>(
        coeffs_truth, *context->GetNttParams(), *context->GetModulusParams()));
    // It's not really desirable to test different targets in a single test, but
    // I can't seem to move these to parameterized tests.
    for (auto target : hwy::SupportedAndGeneratedTargets()) {
      hwy::SetSupportedTargetsForTest(target);
      std::vector<TypeParam> coeffs = coeffs_original;
      ASSERT_OK(ForwardNumberTheoreticTransformHwy<TypeParam>(
          coeffs, *context->GetNttParams(), *context->GetModulusParams()));

      EXPECT_EQ(coeffs.size(), (1 << log_n));
      for (int i = 0; i < coeffs.size(); ++i) {
        EXPECT_EQ(coeffs[i], coeffs_truth[i]);
      }
      hwy::SetSupportedTargetsForTest(0);
    }
  }
}

}  // namespace
}  // namespace rlwe::internal

#endif
