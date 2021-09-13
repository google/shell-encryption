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

#include "shell_encryption/opt/lazy_polynomial.h"

#include <vector>

#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include "absl/status/status.h"
#include "shell_encryption/context.h"
#include "shell_encryption/montgomery.h"
#include "shell_encryption/polynomial.h"  // Needs to be included since Polynomial is forward declarated.
#include "shell_encryption/prng/single_thread_hkdf_prng.h"
#include "shell_encryption/statusor.h"
#include "shell_encryption/testing/parameters.h"
#include "shell_encryption/testing/status_matchers.h"
#include "shell_encryption/testing/status_testing.h"

namespace {

template <typename ModularInt>
class LazyPolynomialTest : public ::testing::Test {
 protected:
  rlwe::StatusOr<rlwe::Polynomial<ModularInt>> Sample(
      const rlwe::RlweContext<ModularInt>* context) {
    RLWE_ASSIGN_OR_RETURN(auto seed,
                          rlwe::SingleThreadHkdfPrng::GenerateSeed());
    RLWE_ASSIGN_OR_RETURN(auto prng, rlwe::SingleThreadHkdfPrng::Create(seed));
    std::vector<ModularInt> coefficients;
    coefficients.reserve(context->GetN());
    for (size_t i = 0; i < context->GetN(); i++) {
      RLWE_ASSIGN_OR_RETURN(
          auto r,
          ModularInt::ImportRandom(prng.get(), context->GetModulusParams()));
      coefficients.push_back(r);
    }
    return rlwe::Polynomial<ModularInt>::ConvertToNtt(
        std::move(coefficients), context->GetNttParams(),
        context->GetModulusParams());
  }
};
TYPED_TEST_SUITE(LazyPolynomialTest, rlwe::testing::ModularIntTypes);

TYPED_TEST(LazyPolynomialTest, CreateExportEmpty) {
  for (const auto& params :
       rlwe::testing::ContextParameters<TypeParam>::Value()) {
    ASSERT_OK_AND_ASSIGN(auto context,
                         rlwe::RlweContext<TypeParam>::Create(params));

    ASSERT_OK_AND_ASSIGN(
        auto empty,
        (rlwe::LazyPolynomial<TypeParam,
                              typename rlwe::internal::BigInt<
                                  typename TypeParam::Int>::value_type>::
             CreateEmpty(context->GetN(), context->GetModulusParams())));
    auto exported = empty.Export(context->GetModulusParams());

    EXPECT_EQ(exported, rlwe::Polynomial<TypeParam>(
                            context->GetN(), context->GetModulusParams()));
  }
}

TYPED_TEST(LazyPolynomialTest, CreateExport) {
  for (const auto& params :
       rlwe::testing::ContextParameters<TypeParam>::Value()) {
    ASSERT_OK_AND_ASSIGN(auto context,
                         rlwe::RlweContext<TypeParam>::Create(params));

    ASSERT_OK_AND_ASSIGN(auto poly1, this->Sample(context.get()));
    ASSERT_OK_AND_ASSIGN(auto poly2, this->Sample(context.get()));

    ASSERT_OK_AND_ASSIGN(
        auto empty,
        (rlwe::LazyPolynomial<
            TypeParam,
            typename rlwe::internal::BigInt<typename TypeParam::Int>::
                value_type>::Create(poly1.Coeffs(), poly2.Coeffs(),
                                    context->GetModulusParams())));
    auto exported = empty.Export(context->GetModulusParams());

    ASSERT_OK_AND_ASSIGN(auto product,
                         poly1.Mul(poly2, context->GetModulusParams()));

    EXPECT_EQ(exported, product);
  }
}

TYPED_TEST(LazyPolynomialTest, ComputeInnerProduct) {
  const size_t kNumberProducts = 10;

  for (const auto& params :
       rlwe::testing::ContextParameters<TypeParam>::Value()) {
    ASSERT_OK_AND_ASSIGN(auto context,
                         rlwe::RlweContext<TypeParam>::Create(params));

    std::vector<rlwe::Polynomial<TypeParam>> poly1s, poly2s;
    for (size_t i = 0; i < kNumberProducts; i++) {
      ASSERT_OK_AND_ASSIGN(auto poly1, this->Sample(context.get()));
      ASSERT_OK_AND_ASSIGN(auto poly2, this->Sample(context.get()));
      poly1s.push_back(std::move(poly1));
      poly2s.push_back(std::move(poly2));
    }

    rlwe::Polynomial<TypeParam> inner_product(context->GetN(),
                                              context->GetModulusParams());
    ASSERT_OK_AND_ASSIGN(
        auto lazy_inner_product,
        (rlwe::LazyPolynomial<TypeParam,
                              typename rlwe::internal::BigInt<
                                  typename TypeParam::Int>::value_type>::
             CreateEmpty(context->GetN(), context->GetModulusParams())));

    for (size_t i = 0; i < kNumberProducts; i++) {
      ASSERT_OK(inner_product.FusedMulAddInPlace(poly1s[i], poly2s[i],
                                                 context->GetModulusParams()));
      ASSERT_OK(lazy_inner_product.FusedMulAddInPlace(
          poly1s[i].Coeffs(), poly2s[i].Coeffs(), context->GetModulusParams()));
    }

    auto exported = lazy_inner_product.Export(context->GetModulusParams());
    EXPECT_EQ(exported, inner_product);
  }
}

TYPED_TEST(LazyPolynomialTest, SameSizeCheckCreation) {
  for (const auto& params :
       rlwe::testing::ContextParameters<TypeParam>::Value()) {
    ASSERT_OK_AND_ASSIGN(auto context,
                         rlwe::RlweContext<TypeParam>::Create(params));

    ASSERT_OK_AND_ASSIGN(auto poly1, this->Sample(context.get()));
    rlwe::Polynomial<TypeParam> poly2(context->GetN() + 1,
                                      context->GetModulusParams());

    EXPECT_THAT((rlwe::LazyPolynomial<
                    TypeParam,
                    typename rlwe::internal::BigInt<typename TypeParam::Int>::
                        value_type>::Create(poly1.Coeffs(), poly2.Coeffs(),
                                            context->GetModulusParams())),
                rlwe::testing::StatusIs(
                    absl::StatusCode::kInvalidArgument,
                    testing::HasSubstr("not all of the same size")));
  }
}

TYPED_TEST(LazyPolynomialTest, SameSizeCheckFusedOperation) {
  for (const auto& params :
       rlwe::testing::ContextParameters<TypeParam>::Value()) {
    ASSERT_OK_AND_ASSIGN(auto context,
                         rlwe::RlweContext<TypeParam>::Create(params));

    ASSERT_OK_AND_ASSIGN(auto poly1, this->Sample(context.get()));
    rlwe::Polynomial<TypeParam> poly2(context->GetN() + 1,
                                      context->GetModulusParams());

    ASSERT_OK_AND_ASSIGN(
        auto empty,
        (rlwe::LazyPolynomial<TypeParam,
                              typename rlwe::internal::BigInt<
                                  typename TypeParam::Int>::value_type>::
             CreateEmpty(context->GetN(), context->GetModulusParams())));

    EXPECT_THAT(empty.FusedMulAddInPlace(poly1.Coeffs(), poly2.Coeffs(),
                                         context->GetModulusParams()),
                rlwe::testing::StatusIs(
                    absl::StatusCode::kInvalidArgument,
                    testing::HasSubstr("not all of the same size")));
  }
}

}  // namespace
