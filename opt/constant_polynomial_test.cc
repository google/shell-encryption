/*
 * Copyright 2021 Google LLC.
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

#include "opt/constant_polynomial.h"

#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include "context.h"
#include "polynomial.h"  // Needs to be included since Polynomial is forward declarated.
#include "testing/parameters.h"
#include "testing/status_matchers.h"
#include "testing/status_testing.h"

namespace {

using rlwe::testing::StatusIs;

template <typename ModularInt>
class ConstantPolynomialTest : public ::testing::Test {};
TYPED_TEST_SUITE(ConstantPolynomialTest, rlwe::testing::ModularIntTypes);

TYPED_TEST(ConstantPolynomialTest, Works) {
  for (const auto& params :
       rlwe::testing::ContextParameters<TypeParam>::Value()) {
    ASSERT_OK_AND_ASSIGN(auto context,
                         rlwe::RlweContext<TypeParam>::Create(params));
    std::vector<typename TypeParam::Int> constant, constant_barrett;

    for (auto length_constant : {1, 10, 1024}) {
      constant.resize(length_constant);
      for (auto length_constant_barrett : {1, 10, 1024}) {
        constant_barrett.resize(length_constant_barrett);

        if (length_constant == length_constant_barrett) {
          ASSERT_OK_AND_ASSIGN(auto p,
                               rlwe::ConstantPolynomial<TypeParam>::Create(
                                   constant, constant_barrett));
        } else {
          EXPECT_THAT(
              rlwe::ConstantPolynomial<TypeParam>::Create(constant,
                                                          constant_barrett),
              StatusIs(::absl::StatusCode::kInvalidArgument,
                       "The vectors of Int do not have the same size."));
        }
      }
    }
  }
}

}  // namespace
