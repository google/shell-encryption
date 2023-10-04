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

#include "shell_encryption/rns/rns_gadget.h"

#include <cstddef>
#include <memory>
#include <vector>

#include <gmock/gmock.h>
#include <gtest/gtest.h>
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

using ::rlwe::testing::StatusIs;
using ::testing::HasSubstr;
using Prng = testing::TestingPrng;

constexpr int kLogGadgetBase = 10;

template <typename ModularInt>
class RnsGadgetTest : public ::testing::Test {
 protected:
  using ModularIntParams = typename ModularInt::Params;
  using Integer = typename ModularInt::Int;
  using RnsParams = testing::RnsParameters<ModularInt>;

  void SetUp() override {
    RnsParams params = testing::GetRnsParameters<ModularInt>();
    auto rns_context = RnsContext<ModularInt>::Create(params.log_n, params.qs,
                                                      params.ps, params.t);
    int level = params.qs.size() - 1;
    rns_context_ = std::make_unique<const RnsContext<ModularInt>>(
        std::move(rns_context.value()));
    moduli_ = rns_context_->MainPrimeModuli();
    q_hats_ = rns_context_->MainPrimeModulusComplements(level).value();
    q_hat_invs_ = rns_context_->MainPrimeModulusCrtFactors(level).value();
    log_bs_.resize(moduli_.size(), kLogGadgetBase);

    prng_ = std::make_unique<testing::TestingPrng>(0);
  }

  std::unique_ptr<const RnsContext<ModularInt>> rns_context_;

  std::vector<const PrimeModulus<ModularInt>*> moduli_;
  std::vector<ModularInt> q_hats_;      // {qi_hat mod qi}_i
  std::vector<ModularInt> q_hat_invs_;  // {qi_hat_inv}_i
  std::vector<size_t> log_bs_;

  std::unique_ptr<Prng> prng_;
};

TYPED_TEST_SUITE(RnsGadgetTest, testing::ModularIntTypes);

TYPED_TEST(RnsGadgetTest, CreateFailsIfLogNIsNotPositive) {
  EXPECT_THAT(
      RnsGadget<TypeParam>::Create(/*log_n=*/-1, this->log_bs_, this->q_hats_,
                                   this->q_hat_invs_, this->moduli_),
      StatusIs(absl::StatusCode::kInvalidArgument,
               HasSubstr("`log_n` must be positive")));
  EXPECT_THAT(
      RnsGadget<TypeParam>::Create(/*log_n=*/0, this->log_bs_, this->q_hats_,
                                   this->q_hat_invs_, this->moduli_),
      StatusIs(absl::StatusCode::kInvalidArgument,
               HasSubstr("`log_n` must be positive")));
}

TYPED_TEST(RnsGadgetTest, CreateFailsIfModuliIsEmpty) {
  EXPECT_THAT(
      RnsGadget<TypeParam>::Create(/*log_n=*/1, this->log_bs_, this->q_hats_,
                                   this->q_hat_invs_, /*moduli=*/{}),
      StatusIs(absl::StatusCode::kInvalidArgument,
               HasSubstr("`moduli` must not be empty")));
}

TYPED_TEST(RnsGadgetTest, CreateFailsIfLogBsHasWrongSize) {
  std::vector<size_t> incorrect_log_bs = this->log_bs_;
  ASSERT_FALSE(incorrect_log_bs.empty());
  incorrect_log_bs.pop_back();  // Smaller than this->moduli_
  EXPECT_THAT(
      RnsGadget<TypeParam>::Create(/*log_n=*/1, incorrect_log_bs, this->q_hats_,
                                   this->q_hat_invs_, this->moduli_),
      StatusIs(absl::StatusCode::kInvalidArgument,
               HasSubstr("`log_bs` must contain")));
}

TYPED_TEST(RnsGadgetTest, CreateFailsIfQHatsHasWrongSize) {
  int num_moduli = this->moduli_.size();
  ASSERT_GT(num_moduli, 0);
  EXPECT_THAT(RnsGadget<TypeParam>::Create(
                  /*log_n=*/1, this->log_bs_,
                  absl::MakeSpan(this->q_hats_).subspan(0, num_moduli - 1),
                  this->q_hat_invs_, this->moduli_),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`q_hats` must contain")));
}

TYPED_TEST(RnsGadgetTest, CreateFailsIfQHatInvsHasWrongSize) {
  int num_moduli = this->moduli_.size();
  ASSERT_GT(num_moduli, 0);
  EXPECT_THAT(RnsGadget<TypeParam>::Create(
                  /*log_n=*/1, this->log_bs_, this->q_hats_,
                  absl::MakeSpan(this->q_hat_invs_).subspan(0, num_moduli - 1),
                  this->moduli_),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`q_hat_invs` must contain")));
}

TYPED_TEST(RnsGadgetTest, DecomposeFailsIfInputPolynomialIsNotCoeffForm) {
  int log_n = this->rns_context_->LogN();
  ASSERT_OK_AND_ASSIGN(
      RnsGadget<TypeParam> gadget,
      RnsGadget<TypeParam>::Create(log_n, this->log_bs_, this->q_hats_,
                                   this->q_hat_invs_, this->moduli_));

  ASSERT_OK_AND_ASSIGN(RnsPolynomial<TypeParam> a,
                       RnsPolynomial<TypeParam>::CreateZero(
                           log_n, this->moduli_, /*is_ntt=*/true));
  ASSERT_TRUE(a.IsNttForm());

  EXPECT_THAT(gadget.Decompose(a, this->moduli_),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`a` must be in coefficient form")));
}

TYPED_TEST(RnsGadgetTest,
           DecomposeFailsIfInputPolynomialHasIncorrectNumberOfModuli) {
  // This test requires at least two RNS moduli.
  int num_moduli = this->moduli_.size();
  if (num_moduli <= 1) {
    return;
  }

  int log_n = this->rns_context_->LogN();
  ASSERT_OK_AND_ASSIGN(
      RnsGadget<TypeParam> gadget,
      RnsGadget<TypeParam>::Create(log_n, this->log_bs_, this->q_hats_,
                                   this->q_hat_invs_, this->moduli_));

  ASSERT_OK_AND_ASSIGN(
      RnsPolynomial<TypeParam> a,
      RnsPolynomial<TypeParam>::CreateZero(
          log_n, absl::MakeSpan(this->moduli_).subspan(0, num_moduli - 1),
          /*is_ntt=*/false));
  ASSERT_FALSE(a.IsNttForm());
  ASSERT_LT(a.NumModuli(), num_moduli);
  EXPECT_THAT(gadget.Decompose(a, this->moduli_),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`a` must be defined with respect to")));
}

TYPED_TEST(RnsGadgetTest, DecomposedDigitsCanRecoverInputPolynomial) {
  int log_n = this->rns_context_->LogN();
  ASSERT_OK_AND_ASSIGN(
      RnsGadget<TypeParam> gadget,
      RnsGadget<TypeParam>::Create(log_n, this->log_bs_, this->q_hats_,
                                   this->q_hat_invs_, this->moduli_));

  ASSERT_OK_AND_ASSIGN(auto a, RnsPolynomial<TypeParam>::SampleUniform(
                                   log_n, this->prng_.get(), this->moduli_));
  ASSERT_OK(a.ConvertToCoeffForm(this->moduli_));

  ASSERT_OK_AND_ASSIGN(auto digits, gadget.Decompose(a, this->moduli_));
  ASSERT_EQ(digits.size(), gadget.Dimension());

  // Compute the inner product <g, g^-1(a)>, which should recover polynomial a.
  ASSERT_OK_AND_ASSIGN(
      auto sum, RnsPolynomial<TypeParam>::CreateZero(log_n, this->moduli_));
  for (int i = 0; i < digits.size(); ++i) {
    ASSERT_OK(digits[i].ConvertToNttForm(this->moduli_));
    RnsPolynomial<TypeParam> x = digits[i];
    ASSERT_OK(x.MulInPlace(gadget.Component(i), this->moduli_));
    ASSERT_OK(sum.AddInPlace(x, this->moduli_));
  }
  ASSERT_OK(sum.ConvertToCoeffForm(this->moduli_));
  EXPECT_EQ(sum, a);
}

TYPED_TEST(RnsGadgetTest, DecomposeMultiplePolynomials) {
  constexpr int k_num_polys = 3;

  int log_n = this->rns_context_->LogN();
  ASSERT_OK_AND_ASSIGN(
      RnsGadget<TypeParam> gadget,
      RnsGadget<TypeParam>::Create(log_n, this->log_bs_, this->q_hats_,
                                   this->q_hat_invs_, this->moduli_));

  std::vector<RnsPolynomial<TypeParam>> polys;
  polys.reserve(k_num_polys);
  for (int i = 0; i < k_num_polys; ++i) {
    ASSERT_OK_AND_ASSIGN(auto poly,
                         RnsPolynomial<TypeParam>::SampleUniform(
                             log_n, this->prng_.get(), this->moduli_));
    ASSERT_OK(poly.ConvertToCoeffForm(this->moduli_));
    polys.push_back(std::move(poly));
  }

  ASSERT_OK_AND_ASSIGN(auto poly_digits,
                       gadget.Decompose(polys, this->moduli_));
  ASSERT_EQ(poly_digits.size(), k_num_polys);
  for (auto const& digits : poly_digits) {
    ASSERT_EQ(digits.size(), gadget.Dimension());
  }

  // For each polynomial in polys, we can recover it from its decomposition
  // `digits` = g^-1(polynomial) in poly_digits by taking the inner product
  // <g, digits>.
  for (int i = 0; i < k_num_polys; ++i) {
    ASSERT_OK_AND_ASSIGN(
        auto sum, RnsPolynomial<TypeParam>::CreateZero(log_n, this->moduli_));
    for (int j = 0; j < poly_digits[i].size(); ++j) {
      ASSERT_OK(poly_digits[i][j].ConvertToNttForm(this->moduli_));
      RnsPolynomial<TypeParam> x = poly_digits[i][j];
      ASSERT_OK(x.MulInPlace(gadget.Component(j), this->moduli_));
      ASSERT_OK(sum.AddInPlace(x, this->moduli_));
    }
    ASSERT_OK(sum.ConvertToCoeffForm(this->moduli_));
    EXPECT_EQ(sum, polys[i]);
  }
}

}  // namespace
}  // namespace rlwe
