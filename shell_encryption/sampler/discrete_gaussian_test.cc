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

#include "shell_encryption/sampler/discrete_gaussian.h"

#include <cmath>
#include <memory>
#include <vector>

#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include "absl/numeric/int128.h"
#include "absl/status/status.h"
#include "shell_encryption/integral_types.h"
#include "shell_encryption/testing/status_matchers.h"
#include "shell_encryption/testing/status_testing.h"
#include "shell_encryption/testing/testing_prng.h"

namespace rlwe {
namespace {

using ::testing::HasSubstr;
using testing::StatusIs;
#ifdef ABSL_HAVE_INTRINSIC_INT128
using TestTypes = ::testing::Types<Uint64, absl::uint128, unsigned __int128>;
#else
using TestTypes = ::testing::Types<Uint64, absl::uint128>;
#endif

// Gaussian parameter of the base sampler.
constexpr double kBaseS = 12.8;

TEST(DiscreteGaussianSampler, CreateFailsIfSBaseIsTooSmall) {
  // The Gaussian parameter of the base sampler, `s_base`, must be at least
  // sqrt(2) * kSmoothParameter, or about 7.55.
  EXPECT_THAT(DiscreteGaussianSampler<Uint64>::Create(/*s_base=*/-1),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`s_base` must be at least")));
  EXPECT_THAT(DiscreteGaussianSampler<Uint64>::Create(/*s_base=*/0),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`s_base` must be at least")));
  EXPECT_THAT(DiscreteGaussianSampler<Uint64>::Create(/*s_base=*/3),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`s_base` must be at least")));
}

TEST(DiscreteGaussianSampler, SampleWithIterationsFailsIfSIsTooSmall) {
  ASSERT_OK_AND_ASSIGN(auto sampler, DiscreteGaussianSampler<Uint64>::Create(
                                         /*s_base=*/kBaseS));
  constexpr int num_iterations = 0;
  auto prng = testing::TestingPrng(0);
  EXPECT_THAT(
      sampler->SampleWithIterations(/*s=*/kBaseS / 2.0, num_iterations, prng),
      StatusIs(absl::StatusCode::kInvalidArgument,
               HasSubstr("`s` must be at least the base s")));
}

TEST(DiscreteGaussianSampler,
     SampleWithIterationsFailsIfNumIterationsIsNegative) {
  ASSERT_OK_AND_ASSIGN(auto sampler, DiscreteGaussianSampler<Uint64>::Create(
                                         /*s_base=*/kBaseS));
  auto prng = testing::TestingPrng(0);
  EXPECT_THAT(sampler->SampleWithIterations(/*s=*/kBaseS,
                                            /*num_iterations=*/-1, prng),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`num_iterations` must be non-negative")));
}

TEST(DiscreteGaussianSampler, NumIterationsFailsIfSIsTooSmall) {
  ASSERT_OK_AND_ASSIGN(auto sampler, DiscreteGaussianSampler<Uint64>::Create(
                                         /*s_base=*/kBaseS));
  EXPECT_THAT(sampler->NumIterations(/*s=*/kBaseS / 2.0),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`s` must be at least the base s")));
}

template <typename Integer>
class DiscreteGaussianSamplerTest : public ::testing::Test {};
TYPED_TEST_SUITE(DiscreteGaussianSamplerTest, TestTypes);

// Checks that the generated samples have bounded magnitude.
TYPED_TEST(DiscreteGaussianSamplerTest, SampleHasBoundedSize) {
  using DGSampler = DiscreteGaussianSampler<TypeParam>;
  constexpr int num_samples = 1024;
  testing::TestingPrng prng(0);
  ASSERT_OK_AND_ASSIGN(auto sampler, DGSampler::Create(kBaseS));

  std::vector<double> gaussian_parameters = {kBaseS, 2 * kBaseS};
  if (sizeof(TypeParam) > 2) {
    gaussian_parameters.push_back(std::exp2(10));
  }

  for (double s : gaussian_parameters) {
    TypeParam expected_bound =
        static_cast<TypeParam>(std::ceil(DGSampler::kTailBoundMultiplier * s));
    for (int i = 0; i < num_samples; ++i) {
      ASSERT_OK_AND_ASSIGN(TypeParam x, sampler->Sample(s, prng));
      bool is_negative = x > DGSampler::kNegativeThreshold;
      if (is_negative) {
        EXPECT_LT(-x, expected_bound);
      } else {
        EXPECT_LT(x, expected_bound);
      }
    }
  }
}

// Checks that the sampler supports very large standard deviation.
TYPED_TEST(DiscreteGaussianSamplerTest, LargeStandardDeviation) {
  using DGSampler = DiscreteGaussianSampler<TypeParam>;
  constexpr int num_samples = 1024;

  if (sizeof(TypeParam) < 8) {
    GTEST_SKIP() << "Integer type must be at least 64 bits for large s";
  }

  // Use a larger base Gaussian parameter to speed up the sampling process.
  double s_base = 64 / std::sqrt(2 * M_PI);
  testing::TestingPrng prng(0);
  ASSERT_OK_AND_ASSIGN(auto sampler,
                       DiscreteGaussianSampler<TypeParam>::Create(s_base));

  std::vector<double> gaussian_parameters = {std::exp2(30), std::exp2(50)};
  for (double s : gaussian_parameters) {
    TypeParam expected_bound =
        static_cast<TypeParam>(std::ceil(DGSampler::kTailBoundMultiplier * s));
    for (int i = 0; i < num_samples; ++i) {
      ASSERT_OK_AND_ASSIGN(TypeParam x, sampler->Sample(s, prng));
      bool is_negative = x > DGSampler::kNegativeThreshold;
      if (is_negative) {
        EXPECT_LT(-x, expected_bound);
      } else {
        EXPECT_LT(x, expected_bound);
      }
    }
  }
}

}  // namespace
}  // namespace rlwe
