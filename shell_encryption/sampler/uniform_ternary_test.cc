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

#include "shell_encryption/sampler/uniform_ternary.h"

#include <memory>
#include <vector>

#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include "absl/status/status.h"
#include "shell_encryption/context.h"
#include "shell_encryption/integral_types.h"
#include "shell_encryption/montgomery.h"
#include "shell_encryption/testing/parameters.h"
#include "shell_encryption/testing/status_matchers.h"
#include "shell_encryption/testing/status_testing.h"
#include "shell_encryption/testing/testing_prng.h"

namespace rlwe {
namespace {

using ::testing::HasSubstr;
using testing::StatusIs;

constexpr int kTestingRounds = 10;

TEST(SampleFromUniformTernary, FailsIfNumCoeffsIsNotPositive) {
  using ModularInt = MontgomeryInt<Uint32>;
  constexpr Uint32 k_modulus = 17;
  ASSERT_OK_AND_ASSIGN(auto mod_params, ModularInt::Params::Create(k_modulus));
  auto prng = testing::TestingPrng(0);
  EXPECT_THAT(SampleFromUniformTernary<ModularInt>(/*num_coeffs=*/-1,
                                                   mod_params.get(), &prng),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`num_coeffs` must be positive")));
  EXPECT_THAT(SampleFromUniformTernary<ModularInt>(/*num_coeffs=*/0,
                                                   mod_params.get(), &prng),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`num_coeffs` must be positive")));
}

TEST(SampleFromUniformTernary, FailsIfNullModParams) {
  using ModularInt = MontgomeryInt<Uint32>;
  constexpr int k_num_coeffs = 32;
  auto prng = testing::TestingPrng(0);
  EXPECT_THAT(SampleFromUniformTernary<ModularInt>(
                  k_num_coeffs, /*mod_params=*/nullptr, &prng),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`mod_params` must not be null")));
}

TEST(SampleFromUniformTernary, FailsIfNullPrng) {
  using ModularInt = MontgomeryInt<Uint32>;
  constexpr Uint32 k_modulus = 17;
  constexpr int k_num_coeffs = 32;
  ASSERT_OK_AND_ASSIGN(auto mod_params, ModularInt::Params::Create(k_modulus));
  EXPECT_THAT(SampleFromUniformTernary<ModularInt>(
                  k_num_coeffs, mod_params.get(), /*prng=*/nullptr),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`prng` must not be null")));
}

template <typename ModularInt>
class SampleFromUniformTernaryTest : public ::testing::Test {};
TYPED_TEST_SUITE(SampleFromUniformTernaryTest, testing::ModularIntTypes);

TYPED_TEST(SampleFromUniformTernaryTest, CheckUpperBoundOnTernaryDistribution) {
  using Integer = typename TypeParam::Int;

  auto prng = std::make_unique<rlwe::testing::TestingPrng>(0);

  for (const auto& params :
       rlwe::testing::ContextParameters<TypeParam>::Value()) {
    ASSERT_OK_AND_ASSIGN(auto context,
                         rlwe::RlweContext<TypeParam>::Create(params));
    Integer modulus = context->GetModulus();
    for (int i = 0; i < kTestingRounds; i++) {
      ASSERT_OK_AND_ASSIGN(
          std::vector<TypeParam> coeffs,
          rlwe::SampleFromUniformTernary<TypeParam>(
              context->GetN(), context->GetModulusParams(), prng.get()));

      // Check that each coefficient is in [-1, 1]
      for (int j = 0; j < context->GetN(); j++) {
        Integer reduced = coeffs[j].ExportInt(context->GetModulusParams());
        if (reduced > (modulus >> 1)) {
          EXPECT_EQ(modulus - reduced, 1);
        } else {
          EXPECT_LE(reduced, 1);
        }
      }
    }
  }
}

}  // namespace
}  // namespace rlwe
