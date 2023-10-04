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

#include "shell_encryption/rns/finite_field_encoder.h"

#include <vector>

#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include "absl/log/check.h"
#include "absl/random/random.h"
#include "absl/status/status.h"
#include "shell_encryption/montgomery.h"
#include "shell_encryption/rns/crt_interpolation.h"
#include "shell_encryption/rns/rns_context.h"
#include "shell_encryption/rns/rns_modulus.h"
#include "shell_encryption/rns/testing/parameters.h"
#include "shell_encryption/rns/testing/testing_utils.h"
#include "shell_encryption/testing/status_matchers.h"
#include "shell_encryption/testing/status_testing.h"

namespace rlwe {

namespace {

using ::rlwe::testing::StatusIs;
using ::testing::HasSubstr;

// Test fixture for FiniteFieldEncoder, which encodes integers modulo t as the
// coefficients of polynomials in Z[X]/(Q, X^N+1). See testing/parameters.h for
// concrete parameters used to instantiate these tests.
template <typename ModularInt>
class FiniteFieldEncoderTest : public ::testing::Test {
 protected:
  using ModularIntParams = typename ModularInt::Params;
  using Integer = typename ModularInt::Int;

  void SetUpRnsParameters(const testing::RnsParameters<ModularInt>& params) {
    // Use parameters suitable for finite field encoding.
    auto rns_context = RnsContext<ModularInt>::Create(params.log_n, params.qs,
                                                      params.ps, params.t);
    CHECK(rns_context.ok());
    rns_context_ = std::make_unique<const RnsContext<ModularInt>>(
        std::move(rns_context.value()));
    moduli_ = rns_context_->MainPrimeModuli();
  }

  void SetUpDefaultRnsParameters() {
    auto all_params =
        testing::GetRnsParametersForFiniteFieldEncoding<ModularInt>();
    CHECK(!all_params.empty());
    SetUpRnsParameters(all_params[0]);
  }

  // Checking if CRT interpolation has enough precision
  // to successfully execute. If it doesn't, we skip the test until
  // Required precision should be (log_Q)+log(\sum_i q_i).
  bool HasInsufficientPrecisionForCrtInterp() {
    using Integer = typename ModularInt::Int;
    using BigInteger = typename testing::CompositeModulus<Integer>::value_type;
    int log_Q = 0;
    BigInteger sum = 0;
    for (auto modulus : this->moduli_) {
      // can't just compute std::ceil(std::log2(Q)) as absl::uint128 does not
      // support std::log2.
      auto modulus_q = modulus->ModParams()->modulus;
      for (BigInteger num_bits = modulus_q; num_bits > 0; num_bits >>= 1) {
        log_Q += 1;
      }
      sum += modulus_q;
    }
    int log_sum = 0;
    for (BigInteger num_bits = sum; num_bits > 0; num_bits >>= 1) {
      log_sum += 1;
    }
    return (log_Q + log_sum >= 8 * sizeof(BigInteger));
  }

  std::unique_ptr<const RnsContext<ModularInt>> rns_context_;
  std::vector<const PrimeModulus<ModularInt>*> moduli_;
};

TYPED_TEST_SUITE(FiniteFieldEncoderTest,
                 testing::ModularIntTypesForFiniteFieldEncoding);

TYPED_TEST(FiniteFieldEncoderTest, CreateFailsIfContextIsNull) {
  EXPECT_THAT(FiniteFieldEncoder<TypeParam>::Create(/*context=*/nullptr),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`context` must not be null")));
}

TYPED_TEST(FiniteFieldEncoderTest, CreateFailsIfPlaintextModulusIsEven) {
  auto all_params =
      testing::GetRnsParametersForFiniteFieldEncoding<TypeParam>();
  ASSERT_FALSE(all_params.empty());
  ASSERT_OK_AND_ASSIGN(auto context, RnsContext<TypeParam>::Create(
                                         all_params[0].log_n, all_params[0].qs,
                                         /*ps=*/{}, /*plaintext_modulus=*/2));
  EXPECT_THAT(
      FiniteFieldEncoder<TypeParam>::Create(&context),
      StatusIs(absl::StatusCode::kInvalidArgument,
               HasSubstr("Plaintext modulus cannot be an even number")));
}

TYPED_TEST(FiniteFieldEncoderTest, CreateSucceeds) {
  for (auto const& params :
       testing::GetRnsParametersForFiniteFieldEncoding<TypeParam>()) {
    this->SetUpRnsParameters(params);
    ASSERT_OK_AND_ASSIGN(auto encoder, FiniteFieldEncoder<TypeParam>::Create(
                                           this->rns_context_.get()));
    EXPECT_EQ(encoder.LogN(), this->rns_context_->LogN());
    EXPECT_EQ(encoder.PlaintextModulus(),
              this->rns_context_->PlaintextModulus());
  }
}

TYPED_TEST(FiniteFieldEncoderTest, EncodeBgvFailsIfMessageVectorIsTooLong) {
  this->SetUpDefaultRnsParameters();
  ASSERT_OK_AND_ASSIGN(auto encoder, FiniteFieldEncoder<TypeParam>::Create(
                                         this->rns_context_.get()));
  int num_coeffs_max = 1 << this->rns_context_->LogN();
  std::vector<typename TypeParam::Int> messages(num_coeffs_max + 1,
                                                0);  // longer than 2^log_n
  EXPECT_THAT(encoder.EncodeBgv(messages, this->moduli_),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr(absl::StrCat("`messages` can contain at most ",
                                              num_coeffs_max, " elements"))));
}

TYPED_TEST(FiniteFieldEncoderTest, EncodeBgvFailsIfEmptyModuli) {
  this->SetUpDefaultRnsParameters();
  ASSERT_OK_AND_ASSIGN(auto encoder, FiniteFieldEncoder<TypeParam>::Create(
                                         this->rns_context_.get()));
  int num_coeffs = 1 << this->rns_context_->LogN();
  std::vector<typename TypeParam::Int> messages(num_coeffs, 0);
  EXPECT_THAT(encoder.EncodeBgv(messages, /*moduli=*/{}),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`moduli` cannot be empty")));
}

TYPED_TEST(FiniteFieldEncoderTest, BgvEncodedResultDecodes) {
  using Integer = typename TypeParam::Int;
  for (auto const& params :
       testing::GetRnsParametersForFiniteFieldEncoding<TypeParam>()) {
    this->SetUpRnsParameters(params);
    Integer t = this->rns_context_->PlaintextModulus();
    ASSERT_OK_AND_ASSIGN(auto encoder, FiniteFieldEncoder<TypeParam>::Create(
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
}

TYPED_TEST(FiniteFieldEncoderTest, BgvEncodedResultCrtDecodes) {
  using Integer = typename TypeParam::Int;
  using BigInteger = typename testing::CompositeModulus<Integer>::value_type;
  for (auto const& params :
       testing::GetRnsParametersForFiniteFieldEncoding<TypeParam>()) {
    this->SetUpRnsParameters(params);
    if (this->HasInsufficientPrecisionForCrtInterp()) {
      continue;
    }
    Integer t = this->rns_context_->PlaintextModulus();
    int num_moduli = this->moduli_.size();
    auto level = num_moduli - 1;
    std::vector<BigInteger> modulus_hats =
        RnsModulusComplements<TypeParam, BigInteger>(this->moduli_);
    ASSERT_OK_AND_ASSIGN(std::vector<TypeParam> modulus_hat_invs,
                         this->rns_context_->MainPrimeModulusCrtFactors(level));

    ASSERT_OK_AND_ASSIGN(auto encoder, FiniteFieldEncoder<TypeParam>::Create(
                                           this->rns_context_.get()));

    // Sample messages mod t and encode them to a RNS polynomial.
    int num_coeffs = 1 << this->rns_context_->LogN();
    // Messages between [0, t)
    std::vector<Integer> messages =
        testing::SampleMessages<Integer>(num_coeffs, t);
    ASSERT_OK_AND_ASSIGN(auto poly, encoder.EncodeBgv(messages, this->moduli_));

    // Now decode the RNS polynomial.
    ASSERT_OK_AND_ASSIGN(auto decoded,
                         encoder.template DecodeBgvWithCrt<BigInteger>(
                             poly, this->moduli_, absl::MakeSpan(modulus_hats),
                             absl::MakeSpan(modulus_hat_invs)));
    EXPECT_EQ(decoded, messages);
  }
}

// Checks if passed in inconsistent-length lists of moduli then DecodeBgvWithCrt
// catches this and fails.
TYPED_TEST(FiniteFieldEncoderTest,
           DecodeBgvWithCrtFailsWithImproperModuliCaseOne) {
  using Integer = typename TypeParam::Int;
  using BigInteger = typename testing::CompositeModulus<Integer>::value_type;
  this->SetUpDefaultRnsParameters();
  int num_moduli = this->moduli_.size();
  auto level = num_moduli - 1;
  std::vector<BigInteger> modulus_hats =
      RnsModulusComplements<TypeParam, BigInteger>(this->moduli_);
  ASSERT_OK_AND_ASSIGN(std::vector<TypeParam> modulus_hat_invs,
                       this->rns_context_->MainPrimeModulusCrtFactors(level));

  ASSERT_OK_AND_ASSIGN(auto encoder, FiniteFieldEncoder<TypeParam>::Create(
                                         this->rns_context_.get()));

  int num_coeffs = 1 << this->rns_context_->LogN();
  std::vector<Integer> messages(num_coeffs, 0);
  ASSERT_OK_AND_ASSIGN(auto poly, encoder.EncodeBgv(messages, this->moduli_));

  auto status = encoder.template DecodeBgvWithCrt<BigInteger>(
      poly, this->moduli_, {}, absl::MakeSpan(modulus_hat_invs));
  EXPECT_THAT(
      status,
      StatusIs(
          absl::StatusCode::kInvalidArgument,
          HasSubstr(
              "`modulus_hats` and `modulus_hat_invs` must have the same ")));
}

TYPED_TEST(FiniteFieldEncoderTest,
           DecodeBgvWithCrtFailsWithImproperModuliCaseTwo) {
  using Integer = typename TypeParam::Int;
  using BigInteger = typename testing::CompositeModulus<Integer>::value_type;
  this->SetUpDefaultRnsParameters();
  int num_moduli = this->moduli_.size();
  auto level = num_moduli - 1;
  std::vector<BigInteger> modulus_hats =
      RnsModulusComplements<TypeParam, BigInteger>(this->moduli_);
  ASSERT_OK_AND_ASSIGN(std::vector<TypeParam> modulus_hat_invs,
                       this->rns_context_->MainPrimeModulusCrtFactors(level));

  ASSERT_OK_AND_ASSIGN(auto encoder, FiniteFieldEncoder<TypeParam>::Create(
                                         this->rns_context_.get()));

  int num_coeffs = 1 << this->rns_context_->LogN();
  std::vector<Integer> messages(num_coeffs, 0);
  ASSERT_OK_AND_ASSIGN(auto poly, encoder.EncodeBgv(messages, this->moduli_));

  auto status = encoder.template DecodeBgvWithCrt<BigInteger>(
      poly, {}, absl::MakeSpan(modulus_hats), absl::MakeSpan(modulus_hat_invs));
  EXPECT_THAT(
      status,
      StatusIs(
          absl::StatusCode::kInvalidArgument,
          HasSubstr(
              "`modulus_hats` and `modulus_hat_invs` must have the same ")));
}

TYPED_TEST(FiniteFieldEncoderTest,
           DecodeBgvWithCrtFailsWithImproperModuliCaseThree) {
  using Integer = typename TypeParam::Int;
  using BigInteger = typename testing::CompositeModulus<Integer>::value_type;
  this->SetUpDefaultRnsParameters();
  int num_moduli = this->moduli_.size();
  auto level = num_moduli - 1;
  std::vector<BigInteger> modulus_hats =
      RnsModulusComplements<TypeParam, BigInteger>(this->moduli_);
  ASSERT_OK_AND_ASSIGN(std::vector<TypeParam> modulus_hat_invs,
                       this->rns_context_->MainPrimeModulusCrtFactors(level));

  ASSERT_OK_AND_ASSIGN(auto encoder, FiniteFieldEncoder<TypeParam>::Create(
                                         this->rns_context_.get()));

  int num_coeffs = 1 << this->rns_context_->LogN();
  std::vector<Integer> messages(num_coeffs, 0);
  ASSERT_OK_AND_ASSIGN(auto poly, encoder.EncodeBgv(messages, this->moduli_));

  auto status = encoder.template DecodeBgvWithCrt<BigInteger>(
      poly, this->moduli_, absl::MakeSpan(modulus_hats), {});
  EXPECT_THAT(
      status,
      StatusIs(
          absl::StatusCode::kInvalidArgument,
          HasSubstr(
              "`modulus_hats` and `modulus_hat_invs` must have the same ")));
}

// Checks that `Encode(messages)` pads the `messages` vector with 0 entries if
// its size is smaller than number of available slots.
TYPED_TEST(FiniteFieldEncoderTest, EncodeBgvCanPadShortMessageVector) {
  using Integer = typename TypeParam::Int;
  for (auto const& params :
       testing::GetRnsParametersForFiniteFieldEncoding<TypeParam>()) {
    this->SetUpRnsParameters(params);
    ASSERT_OK_AND_ASSIGN(auto encoder, FiniteFieldEncoder<TypeParam>::Create(
                                           this->rns_context_.get()));

    // First, encode a vector with just one element.
    int num_coeffs = 1 << this->rns_context_->LogN();
    Integer t = this->rns_context_->PlaintextModulus();
    std::vector<Integer> messages1 = testing::SampleMessages<Integer>(1, t);
    ASSERT_OK_AND_ASSIGN(auto poly1,
                         encoder.EncodeBgv(messages1, this->moduli_));

    ASSERT_OK_AND_ASSIGN(auto decoded1,
                         encoder.DecodeBgv(poly1, this->moduli_));
    // Decoded vector should have size equal to the number of available slots.
    EXPECT_EQ(decoded1.size(), num_coeffs);
    EXPECT_EQ(decoded1[0], messages1[0]);
    for (int i = 1; i < num_coeffs; ++i) {
      EXPECT_EQ(decoded1[i], 0);  // Remaining slots should be filled with 0.
    }

    // Next, encode a vector with more than one elements but shorter than
    // a full length vector.
    std::vector<Integer> messages2 =
        testing::SampleMessages<Integer>(num_coeffs - 1, t);
    ASSERT_OK_AND_ASSIGN(auto poly2,
                         encoder.EncodeBgv(messages2, this->moduli_));
    ASSERT_OK_AND_ASSIGN(auto decoded2,
                         encoder.DecodeBgv(poly2, this->moduli_));
    EXPECT_EQ(decoded2.size(), num_coeffs);
    for (int i = 0; i < num_coeffs - 1; ++i) {
      EXPECT_EQ(decoded2[i], messages2[i]);
    }
    EXPECT_EQ(decoded2[num_coeffs - 1], 0);  // The last slot is filled with 0.
  }
}

TYPED_TEST(FiniteFieldEncoderTest, NoisyBgvPlaintextDecodes) {
  using Integer = typename TypeParam::Int;
  for (auto const& params :
       testing::GetRnsParametersForFiniteFieldEncoding<TypeParam>()) {
    this->SetUpRnsParameters(params);
    ASSERT_OK_AND_ASSIGN(auto encoder, FiniteFieldEncoder<TypeParam>::Create(
                                           this->rns_context_.get()));

    // Sample messages mod t and encode them to a RNS polynomial.
    int log_n = this->rns_context_->LogN();
    Integer t = this->rns_context_->PlaintextModulus();
    std::vector<Integer> messages =
        testing::SampleMessages<Integer>(1 << log_n, t);
    ASSERT_OK_AND_ASSIGN(RnsPolynomial<TypeParam> noisy,
                         encoder.EncodeBgv(messages, this->moduli_));
    ASSERT_OK(noisy.ConvertToCoeffForm(this->moduli_));

    // Sample an error polynomial and then add to the plaintext polynomial.
    RnsPolynomial<TypeParam> error =
        testing::SampleTernaryNoises<TypeParam>(log_n, this->moduli_);
    ASSERT_OK(error.MulInPlace(t, this->moduli_));
    ASSERT_OK(noisy.AddInPlace(error, this->moduli_));

    // Now decode the RNS polynomial.
    ASSERT_OK_AND_ASSIGN(auto decoded, encoder.DecodeBgv(noisy, this->moduli_));
    EXPECT_EQ(decoded, messages);
  }
}

TYPED_TEST(FiniteFieldEncoderTest, NoisyBgvPlaintextCrtDecodes) {
  using Integer = typename TypeParam::Int;
  using BigInteger = typename testing::CompositeModulus<Integer>::value_type;
  for (auto const& params :
       testing::GetRnsParametersForFiniteFieldEncoding<TypeParam>()) {
    this->SetUpRnsParameters(params);
    if (this->HasInsufficientPrecisionForCrtInterp()) {
      continue;
    }
    int num_moduli = this->moduli_.size();
    auto level = num_moduli - 1;
    std::vector<BigInteger> modulus_hats =
        RnsModulusComplements<TypeParam, BigInteger>(this->moduli_);
    ASSERT_OK_AND_ASSIGN(std::vector<TypeParam> modulus_hat_invs,
                         this->rns_context_->MainPrimeModulusCrtFactors(level));

    ASSERT_OK_AND_ASSIGN(auto encoder, FiniteFieldEncoder<TypeParam>::Create(
                                           this->rns_context_.get()));

    // Sample messages mod t and encode them to a RNS polynomial.
    int log_n = this->rns_context_->LogN();
    Integer t = this->rns_context_->PlaintextModulus();
    std::vector<Integer> messages =
        testing::SampleMessages<Integer>(1 << log_n, t);
    ASSERT_OK_AND_ASSIGN(RnsPolynomial<TypeParam> noisy,
                         encoder.EncodeBgv(messages, this->moduli_));
    ASSERT_OK(noisy.ConvertToCoeffForm(this->moduli_));

    // Sample an error polynomial and then add to the plaintext polynomial.
    RnsPolynomial<TypeParam> error =
        testing::SampleTernaryNoises<TypeParam>(log_n, this->moduli_);
    ASSERT_OK(error.MulInPlace(t, this->moduli_));
    ASSERT_OK(noisy.AddInPlace(error, this->moduli_));

    // Now decode the RNS polynomial.
    ASSERT_OK_AND_ASSIGN(auto decoded,
                         encoder.template DecodeBgvWithCrt<BigInteger>(
                             noisy, this->moduli_, absl::MakeSpan(modulus_hats),
                             absl::MakeSpan(modulus_hat_invs)));
    EXPECT_EQ(decoded, messages);
  }
}

TYPED_TEST(FiniteFieldEncoderTest, EncodeIsAdditiveHomomorphic) {
  using Integer = typename TypeParam::Int;
  for (auto const& params :
       testing::GetRnsParametersForFiniteFieldEncoding<TypeParam>()) {
    this->SetUpRnsParameters(params);
    ASSERT_OK_AND_ASSIGN(auto encoder, FiniteFieldEncoder<TypeParam>::Create(
                                           this->rns_context_.get()));

    // Sample two vectors of messages mod t and encode them to RNS polynomials.
    int num_coeffs = 1 << this->rns_context_->LogN();
    Integer t = this->rns_context_->PlaintextModulus();
    std::vector<Integer> xs = testing::SampleMessages<Integer>(num_coeffs, t);
    std::vector<Integer> ys = testing::SampleMessages<Integer>(num_coeffs, t);
    ASSERT_OK_AND_ASSIGN(auto poly0, encoder.EncodeBgv(xs, this->moduli_));
    ASSERT_OK_AND_ASSIGN(auto poly1, encoder.EncodeBgv(ys, this->moduli_));

    // poly0 + poly1 should encode (x + y) mod t for x, y in xs, ys.
    ASSERT_OK_AND_ASSIGN(auto sum, poly0.Add(poly1, this->moduli_));

    // Now decode sum and check results.
    ASSERT_OK_AND_ASSIGN(auto decoded_sum,
                         encoder.DecodeBgv(sum, this->moduli_));
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
}

TYPED_TEST(FiniteFieldEncoderTest, EncodeBgvIsMultiplicativeHomomorphic) {
  using Integer = typename TypeParam::Int;
  using BigInteger = typename TypeParam::BigInt;
  for (auto const& params :
       testing::GetRnsParametersForFiniteFieldEncoding<TypeParam>()) {
    this->SetUpRnsParameters(params);
    ASSERT_OK_AND_ASSIGN(auto encoder, FiniteFieldEncoder<TypeParam>::Create(
                                           this->rns_context_.get()));

    // Sample two vectors of messages mod t and encode them to RNS polynomials.
    int num_coeffs = 1 << this->rns_context_->LogN();
    Integer t = this->rns_context_->PlaintextModulus();
    std::vector<Integer> xs = testing::SampleMessages<Integer>(num_coeffs, t);
    std::vector<Integer> ys = testing::SampleMessages<Integer>(num_coeffs, t);
    ASSERT_OK_AND_ASSIGN(auto poly0, encoder.EncodeBgv(xs, this->moduli_));
    ASSERT_OK_AND_ASSIGN(auto poly1, encoder.EncodeBgv(ys, this->moduli_));

    // poly0 * poly1 should encode pointwise (x * y) mod t for x, y in xs, ys.
    ASSERT_OK_AND_ASSIGN(auto product, poly0.Mul(poly1, this->moduli_));

    // Now decode sum and check results.
    ASSERT_OK_AND_ASSIGN(auto decoded,
                         encoder.DecodeBgv(product, this->moduli_));
    ASSERT_EQ(decoded.size(), num_coeffs);
    for (int i = 0; i < num_coeffs; ++i) {
      BigInteger expected = static_cast<BigInteger>(xs[i]);
      expected *= static_cast<BigInteger>(ys[i]);
      expected %= static_cast<BigInteger>(t);
      EXPECT_EQ(decoded[i], expected);
    }
  }
}

}  // namespace

}  // namespace rlwe
