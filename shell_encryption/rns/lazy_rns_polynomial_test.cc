// Copyright 2024 Google LLC
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include "shell_encryption/rns/lazy_rns_polynomial.h"

#include <memory>
#include <vector>

#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include "absl/log/check.h"
#include "absl/status/status.h"
#include "shell_encryption/rns/rns_context.h"
#include "shell_encryption/rns/rns_modulus.h"
#include "shell_encryption/rns/rns_polynomial.h"
#include "shell_encryption/rns/testing/parameters.h"
#include "shell_encryption/testing/parameters.h"
#include "shell_encryption/testing/status_matchers.h"
#include "shell_encryption/testing/status_testing.h"
#include "shell_encryption/testing/testing_prng.h"

namespace rlwe {
namespace {

using Prng = testing::TestingPrng;
using ::rlwe::testing::StatusIs;
using ::testing::HasSubstr;

template <typename ModularInt>
class LazyRnsPolynomialTest : public ::testing::Test {
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

  std::unique_ptr<const RnsContext<ModularInt>> rns_context_;
  std::vector<const PrimeModulus<ModularInt>*> moduli_;
  std::unique_ptr<Prng> prng_;
};

TYPED_TEST_SUITE(LazyRnsPolynomialTest, testing::ModularIntTypes);

template <typename ModularInt>
class LazyRnsPolynomialNegativeTest : public LazyRnsPolynomialTest<ModularInt> {
};

TYPED_TEST_SUITE(LazyRnsPolynomialNegativeTest,
                 testing::ModularIntTypesForNegativeTests);

TYPED_TEST(LazyRnsPolynomialNegativeTest, CreateFailsIfPolynomialNotInNttForm) {
  int log_n = this->rns_context_->LogN();
  ASSERT_OK_AND_ASSIGN(auto poly_not_ntt, RnsPolynomial<TypeParam>::CreateZero(
                                              log_n, this->moduli_,
                                              /*is_ntt=*/false));
  ASSERT_OK_AND_ASSIGN(
      auto poly_ntt, RnsPolynomial<TypeParam>::CreateZero(log_n, this->moduli_,
                                                          /*is_ntt=*/true));
  EXPECT_THAT(
      LazyRnsPolynomial<TypeParam>::Create(poly_not_ntt, poly_ntt,
                                           this->moduli_),
      StatusIs(absl::StatusCode::kInvalidArgument, HasSubstr("NTT form")));
  EXPECT_THAT(
      LazyRnsPolynomial<TypeParam>::Create(poly_ntt, poly_not_ntt,
                                           this->moduli_),
      StatusIs(absl::StatusCode::kInvalidArgument, HasSubstr("NTT form")));
}

TYPED_TEST(LazyRnsPolynomialNegativeTest,
           CreateFailsIfPolynomialDegreeMismatch) {
  int log_n = this->rns_context_->LogN();
  ASSERT_OK_AND_ASSIGN(
      auto poly, RnsPolynomial<TypeParam>::CreateZero(log_n, this->moduli_));
  ASSERT_OK_AND_ASSIGN(auto poly_short, RnsPolynomial<TypeParam>::CreateZero(
                                            log_n - 1, this->moduli_));
  EXPECT_THAT(
      LazyRnsPolynomial<TypeParam>::Create(poly, poly_short, this->moduli_),
      StatusIs(absl::StatusCode::kInvalidArgument,
               HasSubstr("same number of coefficients")));
  EXPECT_THAT(
      LazyRnsPolynomial<TypeParam>::Create(poly_short, poly, this->moduli_),
      StatusIs(absl::StatusCode::kInvalidArgument,
               HasSubstr("same number of coefficients")));
}

TYPED_TEST(LazyRnsPolynomialNegativeTest,
           CreateFailsIfPolynomialModuliMismatch) {
  int log_n = this->rns_context_->LogN();
  int num_moduli = this->moduli_.size();
  ASSERT_GE(num_moduli, 2) << "Must have at least two moduli.";
  ASSERT_OK_AND_ASSIGN(
      auto poly_reduced_moduli,
      RnsPolynomial<TypeParam>::CreateZero(
          log_n, absl::MakeSpan(this->moduli_).subspan(0, num_moduli - 1)));
  ASSERT_OK_AND_ASSIGN(
      auto poly, RnsPolynomial<TypeParam>::CreateZero(log_n, this->moduli_));
  EXPECT_THAT(LazyRnsPolynomial<TypeParam>::Create(poly_reduced_moduli, poly,
                                                   this->moduli_),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("must be defined wrt `moduli`")));
  EXPECT_THAT(LazyRnsPolynomial<TypeParam>::Create(poly, poly_reduced_moduli,
                                                   this->moduli_),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("must be defined wrt `moduli`")));
}

TYPED_TEST(LazyRnsPolynomialNegativeTest, CreateZeroFailsIfLogNIsNotPositive) {
  EXPECT_THAT(
      LazyRnsPolynomial<TypeParam>::CreateZero(/*log_n=*/-1, this->moduli_),
      StatusIs(absl::StatusCode::kInvalidArgument,
               HasSubstr("`log_n` must be positive")));
  EXPECT_THAT(
      LazyRnsPolynomial<TypeParam>::CreateZero(/*log_n=*/0, this->moduli_),
      StatusIs(absl::StatusCode::kInvalidArgument,
               HasSubstr("`log_n` must be positive")));
}

TYPED_TEST(LazyRnsPolynomialNegativeTest, CreateZeroFailsIfModuliIsEmpty) {
  int log_n = this->rns_context_->LogN();
  EXPECT_THAT(LazyRnsPolynomial<TypeParam>::CreateZero(log_n,
                                                       /*moduli=*/{}),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`moduli` must not be empty")));
}

TYPED_TEST(LazyRnsPolynomialNegativeTest, ExportFailsIfModuliMismatch) {
  int log_n = this->rns_context_->LogN();
  int num_moduli = this->moduli_.size();
  ASSERT_GE(num_moduli, 2) << "Must have at least two moduli.";
  ASSERT_OK_AND_ASSIGN(auto lazy, LazyRnsPolynomial<TypeParam>::CreateZero(
                                      log_n, this->moduli_));
  EXPECT_THAT(
      lazy.Export(absl::MakeSpan(this->moduli_).subspan(0, num_moduli - 1)),
      StatusIs(absl::StatusCode::kInvalidArgument,
               HasSubstr("`moduli` does not contain enough")));
}

TYPED_TEST(LazyRnsPolynomialNegativeTest,
           FusedMulAddFailsIfPolynomialNotInNttForm) {
  int log_n = this->rns_context_->LogN();
  ASSERT_OK_AND_ASSIGN(
      auto poly_ntt, RnsPolynomial<TypeParam>::CreateZero(log_n, this->moduli_,
                                                          /*is_ntt=*/true));
  ASSERT_OK_AND_ASSIGN(auto poly_not_ntt,
                       RnsPolynomial<TypeParam>::CreateZero(
                           log_n, this->moduli_, /*is_ntt=*/false));
  ASSERT_OK_AND_ASSIGN(auto lazy, LazyRnsPolynomial<TypeParam>::CreateZero(
                                      log_n, this->moduli_));
  EXPECT_THAT(
      lazy.FusedMulAddInPlace(poly_ntt, poly_not_ntt, this->moduli_),
      StatusIs(absl::StatusCode::kInvalidArgument, HasSubstr("NTT form")));
  EXPECT_THAT(
      lazy.FusedMulAddInPlace(poly_not_ntt, poly_ntt, this->moduli_),
      StatusIs(absl::StatusCode::kInvalidArgument, HasSubstr("NTT form")));
}

TYPED_TEST(LazyRnsPolynomialNegativeTest,
           FusedMulAddFailsIfPolynomialModuliMismatch) {
  int log_n = this->rns_context_->LogN();
  int num_moduli = this->moduli_.size();
  ASSERT_GE(num_moduli, 2) << "Must have at least two moduli.";
  auto moduli_reduced =
      absl::MakeSpan(this->moduli_).subspan(0, num_moduli - 1);
  ASSERT_OK_AND_ASSIGN(auto poly_reduced, RnsPolynomial<TypeParam>::CreateZero(
                                              log_n, moduli_reduced));
  ASSERT_OK_AND_ASSIGN(
      auto poly, RnsPolynomial<TypeParam>::CreateZero(log_n, this->moduli_));
  ASSERT_OK_AND_ASSIGN(auto lazy, LazyRnsPolynomial<TypeParam>::CreateZero(
                                      log_n, this->moduli_));
  EXPECT_THAT(lazy.FusedMulAddInPlace(poly, poly_reduced, this->moduli_),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("must all be defined wrt `moduli`")));
  EXPECT_THAT(lazy.FusedMulAddInPlace(poly_reduced, poly, this->moduli_),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("must all be defined wrt `moduli`")));
  EXPECT_THAT(lazy.FusedMulAddInPlace(poly, poly, moduli_reduced),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("must all be defined wrt `moduli`")));
}

TYPED_TEST(LazyRnsPolynomialNegativeTest,
           FusedMulAddFailsIfPolynomialDegreeMismatch) {
  int log_n = this->rns_context_->LogN();
  ASSERT_OK_AND_ASSIGN(
      auto poly, RnsPolynomial<TypeParam>::CreateZero(log_n, this->moduli_));
  ASSERT_OK_AND_ASSIGN(auto poly_small, RnsPolynomial<TypeParam>::CreateZero(
                                            log_n - 1, this->moduli_));
  ASSERT_OK_AND_ASSIGN(auto lazy, LazyRnsPolynomial<TypeParam>::CreateZero(
                                      log_n, this->moduli_));
  EXPECT_THAT(lazy.FusedMulAddInPlace(poly, poly_small, this->moduli_),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("must have the same number of coefficients")));
  EXPECT_THAT(lazy.FusedMulAddInPlace(poly_small, poly, this->moduli_),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("must have the same number of coefficients")));
  EXPECT_THAT(lazy.FusedMulAddInPlace(poly_small, poly_small, this->moduli_),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("must have the same number of coefficients as "
                                 "this lazy polynomial")));
}

TYPED_TEST(LazyRnsPolynomialTest, CreateExport) {
  int log_n = this->rns_context_->LogN();
  ASSERT_OK_AND_ASSIGN(RnsPolynomial<TypeParam> a,
                       RnsPolynomial<TypeParam>::SampleUniform(
                           log_n, this->prng_.get(), this->moduli_));
  ASSERT_OK_AND_ASSIGN(RnsPolynomial<TypeParam> b,
                       RnsPolynomial<TypeParam>::SampleUniform(
                           log_n, this->prng_.get(), this->moduli_));
  ASSERT_OK_AND_ASSIGN(
      auto lazy, LazyRnsPolynomial<TypeParam>::Create(a, b, this->moduli_));
  ASSERT_OK_AND_ASSIGN(RnsPolynomial<TypeParam> lazy_export,
                       lazy.Export(this->moduli_));
  ASSERT_OK_AND_ASSIGN(RnsPolynomial<TypeParam> expected,
                       a.Mul(b, this->moduli_));
  EXPECT_EQ(lazy_export, expected);
}

TYPED_TEST(LazyRnsPolynomialTest, FusedMulAddInPlace) {
  int log_n = this->rns_context_->LogN();
  int dimensions[] = {1, 2, 16, 128, 1024};
  for (int num_polys : dimensions) {
    std::vector<RnsPolynomial<TypeParam>> as, bs;
    as.reserve(num_polys);
    bs.reserve(num_polys);
    for (int i = 0; i < num_polys; ++i) {
      ASSERT_OK_AND_ASSIGN(RnsPolynomial<TypeParam> a,
                           RnsPolynomial<TypeParam>::SampleUniform(
                               log_n, this->prng_.get(), this->moduli_));
      ASSERT_OK_AND_ASSIGN(RnsPolynomial<TypeParam> b,
                           RnsPolynomial<TypeParam>::SampleUniform(
                               log_n, this->prng_.get(), this->moduli_));
      as.push_back(std::move(a));
      bs.push_back(std::move(b));
    }
    // Compute the dot product between vectors of random polynomials
    // using the lazy polynomial's `FusedMulAddInPlace` function.
    ASSERT_OK_AND_ASSIGN(auto lazy, LazyRnsPolynomial<TypeParam>::CreateZero(
                                        log_n, this->moduli_));
    for (int i = 0; i < num_polys; ++i) {
      ASSERT_OK(lazy.FusedMulAddInPlace(as[i], bs[i], this->moduli_));
    }

    // Check the results.
    ASSERT_OK_AND_ASSIGN(auto expected,
                         RnsPolynomial<TypeParam>::CreateZero(
                             log_n, this->moduli_, /*is_ntt=*/true));
    for (int i = 0; i < num_polys; ++i) {
      ASSERT_OK(expected.FusedMulAddInPlace(as[i], bs[i], this->moduli_));
    }
    ASSERT_OK_AND_ASSIGN(RnsPolynomial<TypeParam> lazy_export,
                         lazy.Export(this->moduli_));
    EXPECT_EQ(lazy_export, expected);
  }
}

}  // namespace
}  // namespace rlwe
