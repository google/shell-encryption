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

#include "shell_encryption/rns/coefficient_encoder.h"

#include <vector>

#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include "absl/log/check.h"
#include "absl/status/status.h"
#include "shell_encryption/montgomery.h"
#include "shell_encryption/rns/rns_context.h"
#include "shell_encryption/rns/rns_modulus.h"
#include "shell_encryption/rns/testing/parameters.h"
#include "shell_encryption/rns/testing/testing_utils.h"
#include "shell_encryption/testing/parameters.h"
#include "shell_encryption/testing/status_matchers.h"
#include "shell_encryption/testing/status_testing.h"

namespace rlwe {

namespace {

using ::rlwe::testing::StatusIs;
using ::testing::HasSubstr;

// Test fixture for CoefficientEncoder, which encodes integers modulo t as the
// coefficients of polynomials in Z[X]/(Q, X^N+1). See testing/parameters.h for
// concrete parameters used to instantiate these tests.
template <typename ModularInt>
class CoefficientEncoderTest : public ::testing::Test {
 protected:
  using ModularIntParams = typename ModularInt::Params;
  using Integer = typename ModularInt::Int;

  void SetUp() override {
    testing::RnsParameters<ModularInt> rns_params =
        testing::GetRnsParameters<ModularInt>();
    auto rns_context = RnsContext<ModularInt>::Create(
        rns_params.log_n, rns_params.qs, rns_params.ps, rns_params.t);
    CHECK(rns_context.ok());
    rns_context_ = std::make_unique<const RnsContext<ModularInt>>(
        std::move(rns_context.value()));
    moduli_ = rns_context_->MainPrimeModuli();
  }

  // Returns the integer modulo q representing the integer `x` (assuming 0 <=
  // `x` < `max_value` < q) using balanced representation, that is, values x in
  // [0, max_value/2] are represented by x (mod q), and values x in
  // (max_value/2, max_value) are represented by
  // -(max_value - x) (mod q).
  Integer BalancedModRepFor(Integer x, Integer max_value, Integer q) {
    Integer max_value_half = max_value >> 1;
    if (x <= max_value_half) {
      return x % q;
    }
    // Now we have the second case x \in [max_value / 2, max_value).
    Integer x_abs = max_value - x;
    return (q - x_abs) % q;
  }

  std::unique_ptr<const RnsContext<ModularInt>> rns_context_;
  std::vector<const PrimeModulus<ModularInt>*> moduli_;
};

TYPED_TEST_SUITE(CoefficientEncoderTest, testing::ModularIntTypes);

TYPED_TEST(CoefficientEncoderTest, CreateFailsIfContextIsNull) {
  EXPECT_THAT(CoefficientEncoder<TypeParam>::Create(/*context=*/nullptr),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`context` must not be null")));
}

TYPED_TEST(CoefficientEncoderTest, CreateSucceeds) {
  ASSERT_OK_AND_ASSIGN(auto encoder, CoefficientEncoder<TypeParam>::Create(
                                         this->rns_context_.get()));
  EXPECT_EQ(encoder.LogN(), this->rns_context_->LogN());
  EXPECT_EQ(encoder.PlaintextModulus(), this->rns_context_->PlaintextModulus());
}

TYPED_TEST(CoefficientEncoderTest, EncodeBgvFailsIfMessageVectorIsTooLong) {
  using Integer = typename TypeParam::Int;
  ASSERT_OK_AND_ASSIGN(auto encoder, CoefficientEncoder<TypeParam>::Create(
                                         this->rns_context_.get()));
  int num_coeffs_max = 1 << this->rns_context_->LogN();
  std::vector<Integer> messages(num_coeffs_max + 1, 0);  // longer than 2^log_n
  EXPECT_THAT(encoder.EncodeBgv(messages, this->moduli_),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr(absl::StrCat("`messages` can contain at most ",
                                              num_coeffs_max, " elements"))));
}

TYPED_TEST(CoefficientEncoderTest, EncodeBgvFailsIfEmptyModuli) {
  using Integer = typename TypeParam::Int;
  ASSERT_OK_AND_ASSIGN(auto encoder, CoefficientEncoder<TypeParam>::Create(
                                         this->rns_context_.get()));
  int num_coeffs = 1 << this->rns_context_->LogN();
  std::vector<Integer> messages(num_coeffs, 0);
  EXPECT_THAT(encoder.EncodeBgv(messages, /*moduli=*/{}),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`moduli` cannot be empty")));
}

TYPED_TEST(CoefficientEncoderTest, DecodeBgvFailsIfWrongNumberOfModuli) {
  using Integer = typename TypeParam::Int;
  ASSERT_OK_AND_ASSIGN(auto encoder, CoefficientEncoder<TypeParam>::Create(
                                         this->rns_context_.get()));
  ASSERT_OK_AND_ASSIGN(auto zero,
                       RnsPolynomial<TypeParam>::CreateZero(
                           this->rns_context_->LogN(), this->moduli_));
  // Run DecodeBgv() against a smaller set of RNS moduli should result in an
  // error.
  int num_moduli = this->moduli_.size();
  EXPECT_THAT(
      encoder.DecodeBgv(
          zero, absl::MakeSpan(this->moduli_).subspan(0, num_moduli - 1)),
      StatusIs(absl::StatusCode::kInvalidArgument,
               HasSubstr(absl::StrCat("`moduli` must contain ", num_moduli,
                                      " RNS moduli"))));
}

// Checks `EncodeBgv()` works as expected when the message vector contains fewer
// than 2^log_n elements, in which case the remaining polynomial coefficients
// should be padded with 0s.
TYPED_TEST(CoefficientEncoderTest, EncodeBgvIsCorrectForShortVectors) {
  using Integer = typename TypeParam::Int;
  ASSERT_OK_AND_ASSIGN(auto encoder, CoefficientEncoder<TypeParam>::Create(
                                         this->rns_context_.get()));
  int num_coeffs = 1 << this->rns_context_->LogN();
  int num_moduli = this->moduli_.size();

  // Empty message vector
  ASSERT_OK_AND_ASSIGN(auto poly0,
                       encoder.EncodeBgv(/*messages=*/{}, this->moduli_));
  EXPECT_EQ(poly0.NumCoeffs(), num_coeffs);
  EXPECT_EQ(poly0.NumModuli(), num_moduli);
  for (int i = 0; i < num_moduli; ++i) {
    auto mod_params_qi = this->moduli_[i]->ModParams();
    auto zero = TypeParam::ImportZero(mod_params_qi);
    for (auto const& coeff : poly0.Coeffs()[i]) {
      EXPECT_EQ(coeff, zero);
    }
  }

  // A vector with two messages.
  Integer plaintext_modulus = this->rns_context_->PlaintextModulus();
  std::vector<Integer> messages =
      testing::SampleMessages<Integer>(2, plaintext_modulus);
  ASSERT_OK_AND_ASSIGN(auto poly1, encoder.EncodeBgv(messages, this->moduli_));
  EXPECT_EQ(poly1.NumCoeffs(), num_coeffs);
  EXPECT_EQ(poly1.NumModuli(), num_moduli);
  for (int i = 0; i < num_moduli; ++i) {
    auto mod_params_qi = this->moduli_[i]->ModParams();
    auto qi = mod_params_qi->modulus;
    auto zero = TypeParam::ImportZero(mod_params_qi);
    int j = 0;
    for (; j < messages.size(); ++j) {
      EXPECT_EQ(poly1.Coeffs()[i][j].ExportInt(mod_params_qi),
                this->BalancedModRepFor(messages[j], plaintext_modulus, qi));
    }
    for (; j < num_coeffs; ++j) {
      EXPECT_EQ(poly1.Coeffs()[i][j], zero);
    }
  }
}

TYPED_TEST(CoefficientEncoderTest, EncodeBgvIsCorrectForWholeVectors) {
  using Integer = typename TypeParam::Int;
  ASSERT_OK_AND_ASSIGN(auto encoder, CoefficientEncoder<TypeParam>::Create(
                                         this->rns_context_.get()));
  int num_coeffs = 1 << this->rns_context_->LogN();
  int num_moduli = this->moduli_.size();
  Integer plaintext_modulus = this->rns_context_->PlaintextModulus();

  // A vector with `num_coeffs` many messages.
  std::vector<Integer> messages =
      testing::SampleMessages<Integer>(num_coeffs, plaintext_modulus);
  ASSERT_OK_AND_ASSIGN(auto poly, encoder.EncodeBgv(messages, this->moduli_));
  EXPECT_EQ(poly.NumCoeffs(), num_coeffs);
  EXPECT_EQ(poly.NumModuli(), num_moduli);
  for (int i = 0; i < num_moduli; ++i) {
    auto mod_params_qi = this->moduli_[i]->ModParams();
    auto qi = mod_params_qi->modulus;
    for (int j = 0; j < num_coeffs; ++j) {
      EXPECT_EQ(poly.Coeffs()[i][j].ExportInt(mod_params_qi),
                this->BalancedModRepFor(messages[j], plaintext_modulus, qi));
    }
  }
}

TYPED_TEST(CoefficientEncoderTest, DecodeBgvIsCorrectForErrorlessEncoding) {
  using Integer = typename TypeParam::Int;
  Integer t = this->rns_context_->PlaintextModulus();
  ASSERT_OK_AND_ASSIGN(auto encoder, CoefficientEncoder<TypeParam>::Create(
                                         this->rns_context_.get()));

  // Sample messages mod t and encode them to a RNS polynomial.
  int num_coeffs = 1 << this->rns_context_->LogN();
  std::vector<Integer> messages =
      testing::SampleMessages<Integer>(num_coeffs, t);
  ASSERT_OK_AND_ASSIGN(auto poly, encoder.EncodeBgv(messages, this->moduli_));

  // Now decode the RNS polynomial.
  ASSERT_OK_AND_ASSIGN(auto decoded, encoder.DecodeBgv(poly, this->moduli_));
  EXPECT_EQ(decoded, messages);
}

TYPED_TEST(CoefficientEncoderTest, DecodeBgvCanRemoveErrors) {
  using Integer = typename TypeParam::Int;
  Integer t = this->rns_context_->PlaintextModulus();
  ASSERT_OK_AND_ASSIGN(auto encoder, CoefficientEncoder<TypeParam>::Create(
                                         this->rns_context_.get()));

  // Sample messages mod t and encode them to a RNS polynomial.
  int log_n = this->rns_context_->LogN();
  int num_coeffs = 1 << log_n;
  std::vector<Integer> messages =
      testing::SampleMessages<Integer>(num_coeffs, t);
  ASSERT_OK_AND_ASSIGN(auto noisy, encoder.EncodeBgv(messages, this->moduli_));

  // Sample an error polynomial and then add to the plaintext polynomial.
  RnsPolynomial<TypeParam> error =
      testing::SampleTernaryNoises<TypeParam>(log_n, this->moduli_);
  ASSERT_OK(error.MulInPlace(t, this->moduli_));
  ASSERT_OK(noisy.AddInPlace(error, this->moduli_));

  // Now decode the noisy polynomial.
  ASSERT_OK_AND_ASSIGN(auto decoded, encoder.DecodeBgv(noisy, this->moduli_));
  EXPECT_EQ(decoded, messages);
}

TYPED_TEST(CoefficientEncoderTest, EncodeBgvIsAdditiveHomomorphic) {
  using Integer = typename TypeParam::Int;
  Integer t = this->rns_context_->PlaintextModulus();
  ASSERT_OK_AND_ASSIGN(auto encoder, CoefficientEncoder<TypeParam>::Create(
                                         this->rns_context_.get()));

  // Sample two vectors of messages mod t and encode them to RNS polynomials.
  int num_coeffs = 1 << this->rns_context_->LogN();
  std::vector<Integer> xs = testing::SampleMessages<Integer>(num_coeffs, t);
  std::vector<Integer> ys = testing::SampleMessages<Integer>(num_coeffs, t);
  ASSERT_OK_AND_ASSIGN(auto poly0, encoder.EncodeBgv(xs, this->moduli_));
  ASSERT_OK_AND_ASSIGN(auto poly1, encoder.EncodeBgv(ys, this->moduli_));

  // poly0 + poly1 should encode (x + y) mod t for x, y in xs, ys.
  ASSERT_OK_AND_ASSIGN(auto sum, poly0.Add(poly1, this->moduli_));

  // Now decode sum and check results.
  ASSERT_OK_AND_ASSIGN(auto decoded_sum, encoder.DecodeBgv(sum, this->moduli_));
  ASSERT_EQ(decoded_sum.size(), num_coeffs);
  for (int i = 0; i < num_coeffs; ++i) {
    EXPECT_EQ(decoded_sum[i], (xs[i] + ys[i]) % t);
  }

  // poly0 - poly1 should encodes (x - y) mod t for x, y in xs, ys.
  ASSERT_OK_AND_ASSIGN(auto diff, poly0.Sub(poly1, this->moduli_));
  ASSERT_OK_AND_ASSIGN(auto decoded_diff,
                       encoder.DecodeBgv(diff, this->moduli_));
  ASSERT_EQ(decoded_diff.size(), num_coeffs);
  for (int i = 0; i < num_coeffs; ++i) {
    bool is_x_larger = xs[i] > ys[i];
    Integer diff_abs = is_x_larger ? xs[i] - ys[i] : ys[i] - xs[i];
    Integer diff_exp = is_x_larger ? diff_abs : t - diff_abs;
    EXPECT_EQ(decoded_diff[i], diff_exp % t);  // make sure we reduce mod t.
  }
}

TYPED_TEST(CoefficientEncoderTest,
           EncodeBgvIsAdditiveHomomorphicWithScalarMultiplication) {
  using Integer = typename TypeParam::Int;

  // Three random scalar multiplications should be sufficient.
  constexpr int k_num_scalars = 3;

  Integer t = this->rns_context_->PlaintextModulus();
  ASSERT_OK_AND_ASSIGN(auto encoder, CoefficientEncoder<TypeParam>::Create(
                                         this->rns_context_.get()));

  // Sample messages mod t and encode them to a RNS polynomial.
  int num_coeffs = 1 << this->rns_context_->LogN();
  std::vector<Integer> xs = testing::SampleMessages<Integer>(num_coeffs, t);
  ASSERT_OK_AND_ASSIGN(auto poly, encoder.EncodeBgv(xs, this->moduli_));

  // Sample some scalars mod t.
  std::vector<Integer> scalars =
      testing::SampleMessages<Integer>(k_num_scalars, t);

  for (auto scalar : scalars) {
    // scalar * poly1 should encode (x + y) mod t for x, y in xs, ys.
    ASSERT_OK_AND_ASSIGN(auto scaled, poly.Mul(scalar, this->moduli_));

    // Decode and check results.
    ASSERT_OK_AND_ASSIGN(auto decoded,
                         encoder.DecodeBgv(scaled, this->moduli_));
    ASSERT_EQ(decoded.size(), num_coeffs);
    for (int i = 0; i < num_coeffs; ++i) {
      EXPECT_EQ(decoded[i], (xs[i] * scalar) % t);
    }
  }
}

}  // namespace

}  // namespace rlwe
