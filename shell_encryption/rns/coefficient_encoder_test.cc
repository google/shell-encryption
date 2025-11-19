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

#include <memory>
#include <vector>

#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include "absl/log/check.h"
#include "absl/status/status.h"
#include "absl/strings/str_cat.h"
#include "absl/types/span.h"
#include "shell_encryption/rns/crt_interpolation.h"
#include "shell_encryption/rns/rns_context.h"
#include "shell_encryption/rns/rns_modulus.h"
#include "shell_encryption/rns/rns_polynomial.h"
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
  // (max_value/2, max_value) are represented by -(max_value - x) (mod q).
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

template <typename ModularInt>
class CoefficientEncoderNegativeTest
    : public CoefficientEncoderTest<ModularInt> {};

TYPED_TEST_SUITE(CoefficientEncoderNegativeTest,
                 testing::ModularIntTypesForNegativeTests);

TYPED_TEST(CoefficientEncoderNegativeTest, CreateFailsIfContextIsNull) {
  EXPECT_THAT(CoefficientEncoder<TypeParam>::Create(/*context=*/nullptr),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`context` must not be null")));
}

TYPED_TEST(CoefficientEncoderNegativeTest,
           EncodeBgvFailsIfMessageVectorIsTooLong) {
  using Integer = typename TypeParam::Int;
  ASSERT_OK_AND_ASSIGN(auto encoder, CoefficientEncoder<TypeParam>::Create(
                                         this->rns_context_.get()));
  int num_coeffs_max = 1 << this->rns_context_->LogN();
  std::vector<Integer> messages(num_coeffs_max + 1, 0);  // longer than 2^log_n
  EXPECT_THAT(encoder.EncodeBgv(messages, this->moduli_),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr(absl::StrCat("`messages` can contain at most ",
                                              num_coeffs_max, " elements"))));
  Integer t = this->rns_context_->PlaintextModulus();
  EXPECT_THAT(encoder.template EncodeBgv<Integer>(messages, t, this->moduli_),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr(absl::StrCat("`messages` can contain at most ",
                                              num_coeffs_max, " elements"))));
}

TYPED_TEST(CoefficientEncoderNegativeTest, EncodeBgvFailsIfEmptyModuli) {
  using Integer = typename TypeParam::Int;
  ASSERT_OK_AND_ASSIGN(auto encoder, CoefficientEncoder<TypeParam>::Create(
                                         this->rns_context_.get()));
  int num_coeffs = 1 << this->rns_context_->LogN();
  std::vector<Integer> messages(num_coeffs, 0);
  EXPECT_THAT(encoder.EncodeBgv(messages, /*moduli=*/{}),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`moduli` cannot be empty")));

  Integer t = this->rns_context_->PlaintextModulus();
  EXPECT_THAT(encoder.template EncodeBgv<Integer>(messages, t, /*moduli=*/{}),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`moduli` cannot be empty")));
}

TYPED_TEST(CoefficientEncoderNegativeTest,
           DecodeBgvFailsIfWrongNumberOfModuli) {
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

TYPED_TEST(CoefficientEncoderNegativeTest,
           DecodeBgvWithBigPlaintextFailsIfWrongNumberOfModuli) {
  using Integer = typename TypeParam::Int;
  ASSERT_OK_AND_ASSIGN(auto encoder, CoefficientEncoder<TypeParam>::Create(
                                         this->rns_context_.get()));
  ASSERT_OK_AND_ASSIGN(auto zero,
                       RnsPolynomial<TypeParam>::CreateZero(
                           this->rns_context_->LogN(), this->moduli_));
  // Run DecodeBgv() which accepts large plaintext modulus, with constants
  // whose sizes do not match the plaintext polynomial.
  int num_moduli = this->moduli_.size();
  int level = num_moduli - 1;
  Integer t = this->rns_context_->PlaintextModulus();
  std::vector<Integer> modulus_hats =
      RnsModulusComplements<TypeParam, Integer>(this->moduli_);
  ASSERT_OK_AND_ASSIGN(std::vector<TypeParam> modulus_hat_invs,
                       this->rns_context_->MainPrimeModulusCrtFactors(level));
  EXPECT_THAT(
      encoder.template DecodeBgv<Integer>(
          zero, t, absl::MakeSpan(this->moduli_).subspan(0, num_moduli - 1),
          modulus_hats, modulus_hat_invs),
      StatusIs(absl::StatusCode::kInvalidArgument,
               HasSubstr(absl::StrCat("`moduli`"))));
  EXPECT_THAT(encoder.template DecodeBgv<Integer>(
                  zero, t, this->moduli_,
                  absl::MakeSpan(modulus_hats).subspan(0, num_moduli - 1),
                  modulus_hat_invs),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr(absl::StrCat("`modulus_hats`"))));
  EXPECT_THAT(encoder.template DecodeBgv<Integer>(
                  zero, t, this->moduli_, modulus_hats,
                  absl::MakeSpan(modulus_hat_invs).subspan(0, num_moduli - 1)),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr(absl::StrCat("`modulus_hat_invs`"))));
}

TYPED_TEST(CoefficientEncoderTest, CreateSucceeds) {
  ASSERT_OK_AND_ASSIGN(auto encoder, CoefficientEncoder<TypeParam>::Create(
                                         this->rns_context_.get()));
  EXPECT_EQ(encoder.LogN(), this->rns_context_->LogN());
  EXPECT_EQ(encoder.PlaintextModulus(), this->rns_context_->PlaintextModulus());
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
    // scalar * poly should encode (scalar * x) mod t for x in xs.
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

TYPED_TEST(CoefficientEncoderTest, EncodeBgvWithLargePlaintextModulus) {
  using Integer = typename TypeParam::Int;
  using BigInteger = typename testing::CompositeModulus<Integer>::value_type;

  int modulus_bits = 0;
  for (auto modulus : this->moduli_) {
    modulus_bits += modulus->ModParams()->log_modulus;
  }
  // Let the plaintext modulus be 2^(log(Q) - 16), which should be satisfied
  // by most parameters.
  if (modulus_bits < 16) {
    GTEST_SKIP() << "Insufficient modulus bits";
  }
  BigInteger plaintext_modulus = 1;
  plaintext_modulus <<= (modulus_bits - 16);

  ASSERT_OK_AND_ASSIGN(auto encoder, CoefficientEncoder<TypeParam>::Create(
                                         this->rns_context_.get()));

  // Encode some special values with small absolute value.
  std::vector<BigInteger> xs;
  xs.push_back(1);
  xs.push_back(static_cast<BigInteger>(plaintext_modulus - 1));  // -1
  xs.push_back(static_cast<BigInteger>(plaintext_modulus - 2));  // -2
  ASSERT_OK_AND_ASSIGN(auto poly, encoder.template EncodeBgv<BigInteger>(
                                      xs, plaintext_modulus, this->moduli_));
  ASSERT_FALSE(poly.IsNttForm());
  for (int i = 0; i < this->moduli_.size(); ++i) {
    // Since values in xs are all small (in absolute values), their CRT coeffs
    // are easy to compute.
    auto mod_params_qi = this->moduli_[i]->ModParams();
    Integer qi = mod_params_qi->modulus;
    EXPECT_EQ(poly.Coeffs()[i][0].ExportInt(mod_params_qi), 1);
    EXPECT_EQ(poly.Coeffs()[i][1].ExportInt(mod_params_qi), qi - 1);
    EXPECT_EQ(poly.Coeffs()[i][2].ExportInt(mod_params_qi), qi - 2);
  }
  // Decode the polynomial.
  auto level = this->moduli_.size() - 1;
  std::vector<BigInteger> modulus_hats =
      RnsModulusComplements<TypeParam, BigInteger>(this->moduli_);
  ASSERT_OK_AND_ASSIGN(std::vector<TypeParam> modulus_hat_invs,
                       this->rns_context_->MainPrimeModulusCrtFactors(level));
  ASSERT_OK_AND_ASSIGN(auto decoded, encoder.template DecodeBgv<BigInteger>(
                                         poly, plaintext_modulus, this->moduli_,
                                         modulus_hats, modulus_hat_invs));
  ASSERT_EQ(decoded.size(), 1 << this->rns_context_->LogN());
  EXPECT_EQ(decoded[0], xs[0]);
  EXPECT_EQ(decoded[1], xs[1]);
  EXPECT_EQ(decoded[2], xs[2]);
  for (int i = 3; i < decoded.size(); ++i) {
    EXPECT_EQ(decoded[i], 0);
  }
}

////////////////////////////////////////////////////////////////////////////////
// Tests for BFV coefficient encoding/decoding
////////////////////////////////////////////////////////////////////////////////

template <typename ModularInt>
class CoefficientEncoderBfvTest : public ::testing::Test {
 protected:
  void SetUp() override {
    testing::RnsParameters<ModularInt> rns_params =
        testing::GetRnsParameters<ModularInt>();
    // BFV encoding currently only supports odd plaintext modulus.
    if (rns_params.t % 2 == 0) {
      rns_params.t++;
    }
    // Create a RnsContext suitable to use BFV encoding, in particular, this
    // creates Montgomery integer parameters for the plaintext modulus.
    auto rns_context = RnsContext<ModularInt>::CreateForBfv(
        rns_params.log_n, rns_params.qs, rns_params.ps, rns_params.t);
    CHECK(rns_context.ok());
    rns_context_ = std::make_unique<const RnsContext<ModularInt>>(
        std::move(rns_context.value()));
    main_moduli_ = rns_context_->MainPrimeModuli();
    aux_moduli_ = rns_context_->AuxPrimeModuli();
  }

  std::unique_ptr<const RnsContext<ModularInt>> rns_context_;
  std::vector<const PrimeModulus<ModularInt>*> main_moduli_;
  std::vector<const PrimeModulus<ModularInt>*> aux_moduli_;
};

TYPED_TEST_SUITE(CoefficientEncoderBfvTest, testing::ModularIntTypes);

TYPED_TEST(CoefficientEncoderBfvTest, EncodeBfvFailsIfMessageVectorIsTooLong) {
  using Integer = typename TypeParam::Int;
  ASSERT_OK_AND_ASSIGN(auto encoder, CoefficientEncoder<TypeParam>::Create(
                                         this->rns_context_.get()));
  int num_coeffs_max = 1 << this->rns_context_->LogN();
  std::vector<Integer> messages(num_coeffs_max + 1, 0);  // longer than 2^log_n
  EXPECT_THAT(encoder.EncodeBfv(messages, this->main_moduli_),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr(absl::StrCat("`messages` can contain at most ",
                                              num_coeffs_max, " elements"))));
}

TYPED_TEST(CoefficientEncoderBfvTest, EncodeBfvFailsIfEmptyModuli) {
  using Integer = typename TypeParam::Int;
  ASSERT_OK_AND_ASSIGN(auto encoder, CoefficientEncoder<TypeParam>::Create(
                                         this->rns_context_.get()));
  int num_coeffs = 1 << this->rns_context_->LogN();
  std::vector<Integer> messages(num_coeffs, 0);
  EXPECT_THAT(encoder.EncodeBfv(messages, /*moduli=*/{}),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`moduli` cannot be empty")));
}

TYPED_TEST(CoefficientEncoderBfvTest,
           EncodeBfvFailsIfNullPlaintextModulusParams) {
  using Integer = typename TypeParam::Int;

  // Create a generic RnsContext which doesn't define Montgomery parameters for
  // the plaintext modulus.
  testing::RnsParameters<TypeParam> rns_params =
      testing::GetRnsParameters<TypeParam>();
  ASSERT_OK_AND_ASSIGN(auto rns_context, RnsContext<TypeParam>::Create(
                                             rns_params.log_n, rns_params.qs,
                                             rns_params.ps, rns_params.t));
  ASSERT_OK_AND_ASSIGN(auto encoder,
                       CoefficientEncoder<TypeParam>::Create(&rns_context));

  // EncodeBfv() should fail as rns_context has null Montgomery parameters.
  int num_coeffs = 1 << this->rns_context_->LogN();
  std::vector<Integer> messages(num_coeffs, 0);
  EXPECT_THAT(encoder.EncodeBfv(messages, rns_context.MainPrimeModuli()),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("RnsContext does not contain a valid "
                                 "plaintext modulus parameters")));
}

TYPED_TEST(CoefficientEncoderBfvTest, DecodeBfvFailsIfWrongNumberOfModuli) {
  ASSERT_OK_AND_ASSIGN(auto encoder, CoefficientEncoder<TypeParam>::Create(
                                         this->rns_context_.get()));
  ASSERT_OK_AND_ASSIGN(auto zero,
                       RnsPolynomial<TypeParam>::CreateZero(
                           this->rns_context_->LogN(), this->main_moduli_));
  // Run DecodeBfv() with a smaller set of RNS moduli should return an error.
  int num_moduli = this->main_moduli_.size();
  EXPECT_THAT(
      encoder.DecodeBfv(
          zero, absl::MakeSpan(this->main_moduli_).subspan(0, num_moduli - 1)),
      StatusIs(absl::StatusCode::kInvalidArgument,
               HasSubstr(absl::StrCat("`moduli` must contain ", num_moduli,
                                      " RNS moduli"))));
}

TYPED_TEST(CoefficientEncoderBfvTest, DecodeBfvFailsIfNullPlaintextModParams) {
  testing::RnsParameters<TypeParam> rns_params =
      testing::GetRnsParameters<TypeParam>();
  // Create a RnsContext not suitable for BFV encoding, in particular, it does
  // not define Montgomery integer parameters for the plaintext modulus.
  ASSERT_OK_AND_ASSIGN(auto context, RnsContext<TypeParam>::Create(
                                         rns_params.log_n, rns_params.qs,
                                         rns_params.ps, rns_params.t));
  ASSERT_EQ(context.PlaintextModulusParams().ModParams(), nullptr);

  // Create an encoder based on the RnsContext.
  ASSERT_OK_AND_ASSIGN(auto encoder,
                       CoefficientEncoder<TypeParam>::Create(&context));

  std::vector<const PrimeModulus<TypeParam>*> moduli =
      context.MainPrimeModuli();
  ASSERT_OK_AND_ASSIGN(
      auto zero, RnsPolynomial<TypeParam>::CreateZero(context.LogN(), moduli));
  // DecodeBfv() should fail as the context doesn't define Montgomery parameters
  // for the plaintext modulus.
  EXPECT_THAT(encoder.DecodeBfv(zero, moduli),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("RnsContext does not contain valid "
                                 "plaintext modulus parameters")));
}

TYPED_TEST(CoefficientEncoderBfvTest, EncodeBfvDecodes) {
  using Integer = typename TypeParam::Int;
  ASSERT_OK_AND_ASSIGN(auto encoder, CoefficientEncoder<TypeParam>::Create(
                                         this->rns_context_.get()));

  // Sample messages mod t and encode the vector to a BFV plaintext polynomial.
  int num_coeffs = 1 << this->rns_context_->LogN();
  Integer t = this->rns_context_->PlaintextModulus();
  std::vector<Integer> messages =
      testing::SampleMessages<Integer>(num_coeffs, t);
  ASSERT_OK_AND_ASSIGN(auto poly,
                       encoder.EncodeBfv(messages, this->main_moduli_));

  // Now decode the RNS polynomial.
  ASSERT_OK_AND_ASSIGN(auto decoded,
                       encoder.DecodeBfv(poly, this->main_moduli_));
  EXPECT_EQ(decoded, messages);
}

TYPED_TEST(CoefficientEncoderBfvTest, EncodeBfvPadsShortMessageVectorWithZero) {
  using Integer = typename TypeParam::Int;
  ASSERT_OK_AND_ASSIGN(auto encoder, CoefficientEncoder<TypeParam>::Create(
                                         this->rns_context_.get()));

  int num_coeffs = 1 << this->rns_context_->LogN();
  int num_moduli = this->main_moduli_.size();

  // Empty message vector
  ASSERT_OK_AND_ASSIGN(auto poly0,
                       encoder.EncodeBfv(/*messages=*/{}, this->main_moduli_));
  ASSERT_EQ(poly0.NumCoeffs(), num_coeffs);
  ASSERT_EQ(poly0.NumModuli(), num_moduli);
  ASSERT_OK(poly0.ConvertToCoeffForm(this->main_moduli_));
  for (int i = 0; i < num_moduli; ++i) {
    auto mod_params_qi = this->main_moduli_[i]->ModParams();
    auto zero = TypeParam::ImportZero(mod_params_qi);
    for (auto const& coeff : poly0.Coeffs()[i]) {
      EXPECT_EQ(coeff, zero);
    }
  }

  // A vector with two messages.
  Integer t = this->rns_context_->PlaintextModulus();
  std::vector<Integer> messages = testing::SampleMessages<Integer>(2, t);
  ASSERT_OK_AND_ASSIGN(auto poly1,
                       encoder.EncodeBfv(messages, this->main_moduli_));
  ASSERT_EQ(poly1.NumCoeffs(), num_coeffs);
  ASSERT_EQ(poly1.NumModuli(), num_moduli);
  ASSERT_OK(poly1.ConvertToCoeffForm(this->main_moduli_));
  for (int i = 0; i < num_moduli; ++i) {
    auto mod_params_qi = this->main_moduli_[i]->ModParams();
    auto zero = TypeParam::ImportZero(mod_params_qi);
    for (int j = messages.size(); j < num_coeffs; ++j) {
      EXPECT_EQ(poly1.Coeffs()[i][j], zero);
    }
  }
  ASSERT_OK_AND_ASSIGN(auto decoded,
                       encoder.DecodeBfv(poly1, this->main_moduli_));
  ASSERT_EQ(decoded.size(), num_coeffs);
  for (int j = 0; j < messages.size(); ++j) {
    EXPECT_EQ(decoded[j], messages[j]);
  }
  for (int j = messages.size(); j < num_coeffs; ++j) {
    EXPECT_EQ(decoded[j], 0);
  }
}

TYPED_TEST(CoefficientEncoderBfvTest, DecodeBfvCanRemoveErrors) {
  using Integer = typename TypeParam::Int;
  Integer t = this->rns_context_->PlaintextModulus();
  ASSERT_OK_AND_ASSIGN(auto encoder, CoefficientEncoder<TypeParam>::Create(
                                         this->rns_context_.get()));

  // Sample messages mod t and encode them to a RNS polynomial.
  int log_n = this->rns_context_->LogN();
  int num_coeffs = 1 << log_n;
  std::vector<Integer> messages =
      testing::SampleMessages<Integer>(num_coeffs, t);
  ASSERT_OK_AND_ASSIGN(auto noisy,
                       encoder.EncodeBfv(messages, this->main_moduli_));
  ASSERT_TRUE(noisy.IsNttForm());
  ASSERT_OK(noisy.ConvertToCoeffForm(this->main_moduli_));

  // Sample an error polynomial and then add to the plaintext polynomial.
  RnsPolynomial<TypeParam> error =
      testing::SampleTernaryNoises<TypeParam>(log_n, this->main_moduli_);
  ASSERT_OK(noisy.AddInPlace(error, this->main_moduli_));

  // Now decode the noisy polynomial.
  ASSERT_OK_AND_ASSIGN(auto decoded,
                       encoder.DecodeBfv(noisy, this->main_moduli_));
  EXPECT_EQ(decoded, messages);
}

TYPED_TEST(CoefficientEncoderBfvTest, EncodeBfvIsAdditiveHomomorphic) {
  using Integer = typename TypeParam::Int;
  Integer t = this->rns_context_->PlaintextModulus();
  ASSERT_OK_AND_ASSIGN(auto encoder, CoefficientEncoder<TypeParam>::Create(
                                         this->rns_context_.get()));

  // Sample two vectors of messages mod t and encode them to RNS polynomials.
  int num_coeffs = 1 << this->rns_context_->LogN();
  std::vector<Integer> xs = testing::SampleMessages<Integer>(num_coeffs, t);
  std::vector<Integer> ys = testing::SampleMessages<Integer>(num_coeffs, t);
  ASSERT_OK_AND_ASSIGN(auto poly0, encoder.EncodeBfv(xs, this->main_moduli_));
  ASSERT_OK_AND_ASSIGN(auto poly1, encoder.EncodeBfv(ys, this->main_moduli_));

  // poly0 + poly1 should encode (x + y) mod t for all x in xs and y in ys.
  ASSERT_OK_AND_ASSIGN(auto sum, poly0.Add(poly1, this->main_moduli_));

  // Now decode sum and check results.
  ASSERT_OK_AND_ASSIGN(auto decoded_sum,
                       encoder.DecodeBfv(sum, this->main_moduli_));
  ASSERT_EQ(decoded_sum.size(), num_coeffs);
  for (int i = 0; i < num_coeffs; ++i) {
    EXPECT_EQ(decoded_sum[i], (xs[i] + ys[i]) % t);
  }

  // poly0 - poly1 should encodes (x - y) mod t for x, y in xs, ys.
  ASSERT_OK_AND_ASSIGN(auto diff, poly0.Sub(poly1, this->main_moduli_));
  ASSERT_OK_AND_ASSIGN(auto decoded_diff,
                       encoder.DecodeBfv(diff, this->main_moduli_));
  ASSERT_EQ(decoded_diff.size(), num_coeffs);
  for (int i = 0; i < num_coeffs; ++i) {
    bool is_x_larger = xs[i] > ys[i];
    Integer diff_abs = is_x_larger ? xs[i] - ys[i] : ys[i] - xs[i];
    Integer diff_exp = is_x_larger ? diff_abs : t - diff_abs;
    EXPECT_EQ(decoded_diff[i], diff_exp % t);  // make sure we reduce mod t.
  }
}

TYPED_TEST(CoefficientEncoderBfvTest,
           EncodeBfvIsAdditiveHomomorphicWithScalarMultiplication) {
  using Integer = typename TypeParam::Int;

  // Three random scalar multiplications should be sufficient.
  constexpr int k_num_scalars = 3;

  ASSERT_OK_AND_ASSIGN(auto encoder, CoefficientEncoder<TypeParam>::Create(
                                         this->rns_context_.get()));
  // Sample messages mod t and encode them to a RNS polynomial.
  Integer t = this->rns_context_->PlaintextModulus();
  int num_coeffs = 1 << this->rns_context_->LogN();
  std::vector<Integer> xs = testing::SampleMessages<Integer>(num_coeffs, t);
  ASSERT_OK_AND_ASSIGN(auto poly, encoder.EncodeBfv(xs, this->main_moduli_));

  // Sample some scalars mod t.
  std::vector<Integer> scalars =
      testing::SampleMessages<Integer>(k_num_scalars, t);
  for (auto scalar : scalars) {
    // scalar * poly should encode (scalar * x) mod t for x in xs.
    ASSERT_OK_AND_ASSIGN(auto scaled, poly.Mul(scalar, this->main_moduli_));

    // Decode and check results.
    ASSERT_OK_AND_ASSIGN(auto decoded,
                         encoder.DecodeBfv(scaled, this->main_moduli_));
    ASSERT_EQ(decoded.size(), num_coeffs);
    for (int i = 0; i < num_coeffs; ++i) {
      EXPECT_EQ(decoded[i], (xs[i] * scalar) % t);
    }
  }
}

}  // namespace

}  // namespace rlwe
