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

#include <cmath>
#include <memory>
#include <vector>

#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include "absl/log/check.h"
#include "absl/random/random.h"
#include "absl/status/status.h"
#include "absl/types/span.h"
#include "shell_encryption/montgomery.h"
#include "shell_encryption/rns/crt_interpolation.h"
#include "shell_encryption/rns/rns_context.h"
#include "shell_encryption/rns/rns_integer.h"
#include "shell_encryption/rns/rns_modulus.h"
#include "shell_encryption/rns/rns_polynomial.h"
#include "shell_encryption/rns/testing/parameters.h"
#include "shell_encryption/rns/testing/testing_utils.h"
#include "shell_encryption/testing/status_matchers.h"
#include "shell_encryption/testing/status_testing.h"

namespace rlwe {

namespace {

using ::rlwe::testing::StatusIs;
using ::testing::HasSubstr;

// Test fixture for FiniteFieldEncoder, which encodes integers modulo t in the
// "slots" of polynomials in Z[X]/(Q, X^N+1). See testing/parameters.h for
// concrete parameters used to instantiate these tests.
template <typename ModularInt>
class FiniteFieldEncoderTest : public ::testing::Test {
 protected:
  using ModularIntParams = typename ModularInt::Params;
  using Integer = typename ModularInt::Int;

  void SetUpRnsParameters(const testing::RnsParameters<ModularInt>& params) {
    // Use parameters suitable for finite field encoding.
    auto rns_context = RnsContext<ModularInt>::CreateForBgvFiniteFieldEncoding(
        params.log_n, params.qs, params.ps, params.t);
    CHECK(rns_context.ok());
    rns_context_ = std::make_unique<const RnsContext<ModularInt>>(
        std::move(rns_context.value()));
    main_moduli_ = rns_context_->MainPrimeModuli();
    aux_moduli_ = rns_context_->AuxPrimeModuli();
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
    for (auto modulus : this->main_moduli_) {
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
  std::vector<const PrimeModulus<ModularInt>*> main_moduli_;
  std::vector<const PrimeModulus<ModularInt>*> aux_moduli_;
};

TYPED_TEST_SUITE(FiniteFieldEncoderTest,
                 testing::ModularIntTypesForFiniteFieldEncoding);

TYPED_TEST(FiniteFieldEncoderTest, CreateFailsIfContextIsNull) {
  EXPECT_THAT(FiniteFieldEncoder<TypeParam>::Create(/*context=*/nullptr),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`context` must not be null")));
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
  EXPECT_THAT(encoder.EncodeBgv(messages, this->main_moduli_),
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
    ASSERT_OK_AND_ASSIGN(auto poly,
                         encoder.EncodeBgv(messages, this->main_moduli_));

    // Now decode the RNS polynomial.
    ASSERT_OK_AND_ASSIGN(auto decoded,
                         encoder.DecodeBgv(poly, this->main_moduli_));
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
    int num_moduli = this->main_moduli_.size();
    auto level = num_moduli - 1;
    std::vector<BigInteger> modulus_hats =
        RnsModulusComplements<TypeParam, BigInteger>(this->main_moduli_);
    ASSERT_OK_AND_ASSIGN(std::vector<TypeParam> modulus_hat_invs,
                         this->rns_context_->MainPrimeModulusCrtFactors(level));

    ASSERT_OK_AND_ASSIGN(auto encoder, FiniteFieldEncoder<TypeParam>::Create(
                                           this->rns_context_.get()));

    // Sample messages mod t and encode them to a RNS polynomial.
    int num_coeffs = 1 << this->rns_context_->LogN();
    // Messages between [0, t)
    std::vector<Integer> messages =
        testing::SampleMessages<Integer>(num_coeffs, t);
    ASSERT_OK_AND_ASSIGN(auto poly,
                         encoder.EncodeBgv(messages, this->main_moduli_));

    // Now decode the RNS polynomial.
    ASSERT_OK_AND_ASSIGN(
        auto decoded,
        encoder.template DecodeBgvWithCrt<BigInteger>(
            poly, this->main_moduli_, absl::MakeSpan(modulus_hats),
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
  int num_moduli = this->main_moduli_.size();
  auto level = num_moduli - 1;
  std::vector<BigInteger> modulus_hats =
      RnsModulusComplements<TypeParam, BigInteger>(this->main_moduli_);
  ASSERT_OK_AND_ASSIGN(std::vector<TypeParam> modulus_hat_invs,
                       this->rns_context_->MainPrimeModulusCrtFactors(level));

  ASSERT_OK_AND_ASSIGN(auto encoder, FiniteFieldEncoder<TypeParam>::Create(
                                         this->rns_context_.get()));

  int num_coeffs = 1 << this->rns_context_->LogN();
  std::vector<Integer> messages(num_coeffs, 0);
  ASSERT_OK_AND_ASSIGN(auto poly,
                       encoder.EncodeBgv(messages, this->main_moduli_));

  auto status = encoder.template DecodeBgvWithCrt<BigInteger>(
      poly, this->main_moduli_, {}, absl::MakeSpan(modulus_hat_invs));
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
  int num_moduli = this->main_moduli_.size();
  auto level = num_moduli - 1;
  std::vector<BigInteger> modulus_hats =
      RnsModulusComplements<TypeParam, BigInteger>(this->main_moduli_);
  ASSERT_OK_AND_ASSIGN(std::vector<TypeParam> modulus_hat_invs,
                       this->rns_context_->MainPrimeModulusCrtFactors(level));

  ASSERT_OK_AND_ASSIGN(auto encoder, FiniteFieldEncoder<TypeParam>::Create(
                                         this->rns_context_.get()));

  int num_coeffs = 1 << this->rns_context_->LogN();
  std::vector<Integer> messages(num_coeffs, 0);
  ASSERT_OK_AND_ASSIGN(auto poly,
                       encoder.EncodeBgv(messages, this->main_moduli_));

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
  int num_moduli = this->main_moduli_.size();
  auto level = num_moduli - 1;
  std::vector<BigInteger> modulus_hats =
      RnsModulusComplements<TypeParam, BigInteger>(this->main_moduli_);
  ASSERT_OK_AND_ASSIGN(std::vector<TypeParam> modulus_hat_invs,
                       this->rns_context_->MainPrimeModulusCrtFactors(level));

  ASSERT_OK_AND_ASSIGN(auto encoder, FiniteFieldEncoder<TypeParam>::Create(
                                         this->rns_context_.get()));

  int num_coeffs = 1 << this->rns_context_->LogN();
  std::vector<Integer> messages(num_coeffs, 0);
  ASSERT_OK_AND_ASSIGN(auto poly,
                       encoder.EncodeBgv(messages, this->main_moduli_));

  auto status = encoder.template DecodeBgvWithCrt<BigInteger>(
      poly, this->main_moduli_, absl::MakeSpan(modulus_hats), {});
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
                         encoder.EncodeBgv(messages1, this->main_moduli_));

    ASSERT_OK_AND_ASSIGN(auto decoded1,
                         encoder.DecodeBgv(poly1, this->main_moduli_));

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
                         encoder.EncodeBgv(messages2, this->main_moduli_));
    ASSERT_OK_AND_ASSIGN(auto decoded2,
                         encoder.DecodeBgv(poly2, this->main_moduli_));

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
                         encoder.EncodeBgv(messages, this->main_moduli_));
    ASSERT_OK(noisy.ConvertToCoeffForm(this->main_moduli_));

    // Sample an error polynomial and then add to the plaintext polynomial.
    RnsPolynomial<TypeParam> error =
        testing::SampleTernaryNoises<TypeParam>(log_n, this->main_moduli_);
    ASSERT_OK(error.MulInPlace(t, this->main_moduli_));
    ASSERT_OK(noisy.AddInPlace(error, this->main_moduli_));

    // Now decode the RNS polynomial.
    ASSERT_OK_AND_ASSIGN(auto decoded,
                         encoder.DecodeBgv(noisy, this->main_moduli_));
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
    int num_moduli = this->main_moduli_.size();
    auto level = num_moduli - 1;
    std::vector<BigInteger> modulus_hats =
        RnsModulusComplements<TypeParam, BigInteger>(this->main_moduli_);
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
                         encoder.EncodeBgv(messages, this->main_moduli_));
    ASSERT_OK(noisy.ConvertToCoeffForm(this->main_moduli_));

    // Sample an error polynomial and then add to the plaintext polynomial.
    RnsPolynomial<TypeParam> error =
        testing::SampleTernaryNoises<TypeParam>(log_n, this->main_moduli_);
    ASSERT_OK(error.MulInPlace(t, this->main_moduli_));
    ASSERT_OK(noisy.AddInPlace(error, this->main_moduli_));

    // Now decode the RNS polynomial.
    ASSERT_OK_AND_ASSIGN(
        auto decoded,
        encoder.template DecodeBgvWithCrt<BigInteger>(
            noisy, this->main_moduli_, absl::MakeSpan(modulus_hats),
            absl::MakeSpan(modulus_hat_invs)));
    EXPECT_EQ(decoded, messages);
  }
}

TYPED_TEST(FiniteFieldEncoderTest, EncodeBgvIsAdditiveHomomorphic) {
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
    ASSERT_OK_AND_ASSIGN(auto poly0, encoder.EncodeBgv(xs, this->main_moduli_));
    ASSERT_OK_AND_ASSIGN(auto poly1, encoder.EncodeBgv(ys, this->main_moduli_));

    // poly0 + poly1 should encode (x + y) mod t for x, y in xs, ys.
    ASSERT_OK_AND_ASSIGN(auto sum, poly0.Add(poly1, this->main_moduli_));

    // Now decode sum and check results.
    ASSERT_OK_AND_ASSIGN(auto decoded_sum,
                         encoder.DecodeBgv(sum, this->main_moduli_));
    ASSERT_EQ(decoded_sum.size(), num_coeffs);
    for (int i = 0; i < num_coeffs; ++i) {
      EXPECT_EQ(decoded_sum[i], (xs[i] + ys[i]) % t);
    }

    // poly0 - poly1 should encodes (x - y) mod t for x, y in xs, ys.
    ASSERT_OK_AND_ASSIGN(auto diff, poly0.Sub(poly1, this->main_moduli_));
    ASSERT_OK_AND_ASSIGN(auto decoded_diff,
                         encoder.DecodeBgv(diff, this->main_moduli_));
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
    ASSERT_OK_AND_ASSIGN(auto poly0, encoder.EncodeBgv(xs, this->main_moduli_));
    ASSERT_OK_AND_ASSIGN(auto poly1, encoder.EncodeBgv(ys, this->main_moduli_));

    // poly0 * poly1 should encode pointwise (x * y) mod t for x, y in xs, ys.
    ASSERT_OK_AND_ASSIGN(auto product, poly0.Mul(poly1, this->main_moduli_));

    // Now decode product and check results.
    ASSERT_OK_AND_ASSIGN(auto decoded,
                         encoder.DecodeBgv(product, this->main_moduli_));
    ASSERT_EQ(decoded.size(), num_coeffs);
    for (int i = 0; i < num_coeffs; ++i) {
      BigInteger expected = static_cast<BigInteger>(xs[i]);
      expected *= static_cast<BigInteger>(ys[i]);
      expected %= static_cast<BigInteger>(t);
      EXPECT_EQ(decoded[i], expected);
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
// BFV encoding tests
////////////////////////////////////////////////////////////////////////////////

template <typename ModularInt>
class FiniteFieldEncoderBfvTest : public ::testing::Test {
 protected:
  void SetUpRnsParameters(const testing::RnsParameters<ModularInt>& params) {
    // Use parameters suitable for finite field encoding in BFV.
    auto rns_context = RnsContext<ModularInt>::CreateForBfvFiniteFieldEncoding(
        params.log_n, params.qs, params.ps, params.t);
    CHECK(rns_context.ok());
    rns_context_ = std::make_unique<const RnsContext<ModularInt>>(
        std::move(rns_context.value()));
    main_moduli_ = rns_context_->MainPrimeModuli();
    aux_moduli_ = rns_context_->AuxPrimeModuli();
  }

  void SetUpDefaultRnsParameters() {
    auto all_params =
        testing::GetBfvParametersForFiniteFieldEncoding<ModularInt>();
    CHECK(!all_params.empty());
    SetUpRnsParameters(all_params[0]);
  }

  std::unique_ptr<const RnsContext<ModularInt>> rns_context_;
  std::vector<const PrimeModulus<ModularInt>*> main_moduli_;
  std::vector<const PrimeModulus<ModularInt>*> aux_moduli_;
};

TYPED_TEST_SUITE(FiniteFieldEncoderBfvTest,
                 testing::ModularIntTypesForFiniteFieldEncoding);

TYPED_TEST(FiniteFieldEncoderBfvTest, EncodeBfvFailsIfMessageVectorIsTooLong) {
  this->SetUpDefaultRnsParameters();
  ASSERT_OK_AND_ASSIGN(auto encoder, FiniteFieldEncoder<TypeParam>::Create(
                                         this->rns_context_.get()));
  int num_coeffs_max = 1 << this->rns_context_->LogN();
  std::vector<typename TypeParam::Int> messages(num_coeffs_max + 1,
                                                0);  // longer than 2^log_n
  EXPECT_THAT(encoder.EncodeBfv(messages, this->main_moduli_),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr(absl::StrCat("`messages` can contain at most ",
                                              num_coeffs_max, " elements"))));
}

TYPED_TEST(FiniteFieldEncoderBfvTest, EncodeBfvFailsIfEmptyModuli) {
  this->SetUpDefaultRnsParameters();
  ASSERT_OK_AND_ASSIGN(auto encoder, FiniteFieldEncoder<TypeParam>::Create(
                                         this->rns_context_.get()));
  int num_coeffs = 1 << this->rns_context_->LogN();
  std::vector<typename TypeParam::Int> messages(num_coeffs, 0);
  EXPECT_THAT(encoder.EncodeBfv(messages, /*moduli=*/{}),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`moduli` cannot be empty")));
}

TYPED_TEST(FiniteFieldEncoderBfvTest, EncodeBfvDecodes) {
  using Integer = typename TypeParam::Int;
  for (auto const& params :
       testing::GetBfvParametersForFiniteFieldEncoding<TypeParam>()) {
    this->SetUpRnsParameters(params);

    ASSERT_OK_AND_ASSIGN(auto encoder, FiniteFieldEncoder<TypeParam>::Create(
                                           this->rns_context_.get()));

    // Sample messages mod t and encode them to a RNS polynomial.
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
}

TYPED_TEST(FiniteFieldEncoderBfvTest, NoisyBfvPlaintextDecodes) {
  using Integer = typename TypeParam::Int;
  for (auto const& params :
       testing::GetBfvParametersForFiniteFieldEncoding<TypeParam>()) {
    this->SetUpRnsParameters(params);
    ASSERT_OK_AND_ASSIGN(auto encoder, FiniteFieldEncoder<TypeParam>::Create(
                                           this->rns_context_.get()));

    // Sample messages mod t and encode them to a RNS polynomial.
    int log_n = this->rns_context_->LogN();
    int num_coeffs = 1 << log_n;
    Integer t = this->rns_context_->PlaintextModulus();

    std::vector<Integer> messages =
        testing::SampleMessages<Integer>(num_coeffs, t);
    ASSERT_OK_AND_ASSIGN(RnsPolynomial<TypeParam> noisy,
                         encoder.EncodeBfv(messages, this->main_moduli_));
    ASSERT_OK(noisy.ConvertToCoeffForm(this->main_moduli_));

    // Sample an error polynomial and then add to the plaintext polynomial.
    RnsPolynomial<TypeParam> error =
        testing::SampleTernaryNoises<TypeParam>(log_n, this->main_moduli_);
    ASSERT_OK(noisy.AddInPlace(error, this->main_moduli_));

    // Now decode the RNS polynomial.
    ASSERT_OK_AND_ASSIGN(auto decoded,
                         encoder.DecodeBfv(noisy, this->main_moduli_));
    EXPECT_EQ(decoded, messages);
  }
}

TYPED_TEST(FiniteFieldEncoderBfvTest, EncodeBfvIsAdditiveHomomorphic) {
  using Integer = typename TypeParam::Int;
  for (auto const& params :
       testing::GetBfvParametersForFiniteFieldEncoding<TypeParam>()) {
    this->SetUpRnsParameters(params);
    ASSERT_OK_AND_ASSIGN(auto encoder, FiniteFieldEncoder<TypeParam>::Create(
                                           this->rns_context_.get()));

    // Sample two vectors of messages mod t and encode them to RNS polynomials.
    int num_coeffs = 1 << this->rns_context_->LogN();
    Integer t = this->rns_context_->PlaintextModulus();
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

    // poly0 - poly1 should encodes (x - y) mod t for all x in xs and y in ys.
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
}

TYPED_TEST(FiniteFieldEncoderBfvTest,
           EncodeBfvIsAdditiveHomomorphicWithScalarMultiplication) {
  using Integer = typename TypeParam::Int;

  // Three random scalar multiplications should be sufficient.
  constexpr int k_num_scalars = 3;

  for (auto const& params :
       testing::GetBfvParametersForFiniteFieldEncoding<TypeParam>()) {
    int plaintext_bits = ceil(log2(static_cast<double>(params.t)));
    int main_modulus_bits = 0;
    for (auto qi : params.qs) {
      main_modulus_bits += floor(log2(static_cast<double>(qi)));
    }
    if (main_modulus_bits < 2 * plaintext_bits) {
      continue;
    }

    this->SetUpRnsParameters(params);
    ASSERT_OK_AND_ASSIGN(auto encoder, FiniteFieldEncoder<TypeParam>::Create(
                                           this->rns_context_.get()));

    // Sample messages mod t and encode them to the slots of a RNS polynomial.
    int num_coeffs = 1 << this->rns_context_->LogN();
    Integer t = this->rns_context_->PlaintextModulus();
    std::vector<Integer> xs = testing::SampleMessages<Integer>(num_coeffs, t);
    ASSERT_OK_AND_ASSIGN(auto poly, encoder.EncodeBfv(xs, this->main_moduli_));
    ASSERT_TRUE(poly.IsNttForm());
    ASSERT_OK(poly.ConvertToCoeffForm(this->main_moduli_));

    // Sample some scalars mod t.
    std::vector<Integer> scalars =
        testing::SampleMessages<Integer>(k_num_scalars, t);
    for (auto scalar : scalars) {
      // Q/t * messages * scalar (mod Q).
      ASSERT_OK_AND_ASSIGN(auto scaled, poly.Mul(scalar, this->main_moduli_));

      // Now decode the scaled polynomial and check results.
      ASSERT_OK_AND_ASSIGN(auto decoded,
                           encoder.DecodeBfv(scaled, this->main_moduli_));
      ASSERT_EQ(decoded.size(), num_coeffs);
      for (int i = 0; i < num_coeffs; ++i) {
        EXPECT_EQ(decoded[i], (xs[i] * scalar) % t);
      }
    }
  }
}

TYPED_TEST(FiniteFieldEncoderBfvTest, EncodeBfvIsMultiplicativeHomomorphic) {
  using Integer = typename TypeParam::Int;
  using BigInteger = typename TypeParam::BigInt;
  for (auto const& params :
       testing::GetBfvParametersForFiniteFieldEncoding<TypeParam>()) {
    // For correct decoding, we will need Q >= norm(poly0 * poly1), which can be
    // bounded by N * t^2 .
    Integer t = params.t;
    int plaintext_bits = ceil(log2(static_cast<double>(t)));
    int main_modulus_bits = 0;
    for (auto qi : params.qs) {
      main_modulus_bits += floor(log2(static_cast<double>(qi)));
    }
    if (main_modulus_bits < 2 * plaintext_bits + params.log_n) {
      continue;
    }

    this->SetUpRnsParameters(params);
    ASSERT_OK_AND_ASSIGN(auto encoder, FiniteFieldEncoder<TypeParam>::Create(
                                           this->rns_context_.get()));

    // Sample two vectors of messages mod t and encode them to RNS polynomials.
    int num_coeffs = 1 << this->rns_context_->LogN();

    std::vector<Integer> xs = testing::SampleMessages<Integer>(num_coeffs, t);
    std::vector<Integer> ys = testing::SampleMessages<Integer>(num_coeffs, t);
    ASSERT_OK_AND_ASSIGN(auto poly0, encoder.EncodeBfv(xs, this->main_moduli_));
    ASSERT_OK_AND_ASSIGN(auto poly1, encoder.EncodeBfv(ys, this->main_moduli_));

    ASSERT_TRUE(poly0.IsNttForm());
    ASSERT_TRUE(poly1.IsNttForm());
    ASSERT_OK(poly0.ConvertToCoeffForm(this->main_moduli_));
    ASSERT_OK(poly1.ConvertToCoeffForm(this->main_moduli_));

    // Multiply two polynomials encoding xs and ys in slots, using the algorithm
    // proposed in https://eprint.iacr.org/2021/204. We have two polynomials
    // poly0 ~= round(Q/t * iNTT(xs)) and poly1 ~= round(Q/t * iNTT(ys)), both
    // modulo Q. Their product is computed as follows:
    // 1) Scale poly0 (mod Q) to poly0_aux = P/t * poly0 (mod P);
    // 2) Convert RNS basis to get poly0_main = P/t * poly0 (mod Q);
    // 3) Convert RNS basis to get poly1_aux = Q/t * poly1 (mod P);
    // 4) Multiply polynomials from the previous steps wrt Q and P, and we get
    //    QP/t^2 * poly0 * poly1 (mod Q) and (mod P);
    // 5) Modulus reduce the polynomial in step 4 by the auxiliary modulus P,
    //    and get Q/t * poly0 * poly1 (mod Q).

    int level = this->main_moduli_.size() - 1;
    ASSERT_OK_AND_ASSIGN(std::vector<TypeParam> q_hat_inv_mod_qs,
                         this->rns_context_->MainPrimeModulusCrtFactors(level));
    // P/t * iNTT(xs) (mod P).
    ASSERT_OK_AND_ASSIGN(
        RnsPolynomial<TypeParam> poly0_scaled_aux,
        poly0.ScaleAndSwitchRnsBasis(
            this->main_moduli_, this->aux_moduli_, q_hat_inv_mod_qs,
            this->rns_context_->MainPrimeModulusInverseAuxResidues(),
            this->rns_context_->AuxModulusResidues()));

    // P/t * iNTT(xs) (mod Q).
    ASSERT_OK_AND_ASSIGN(
        RnsPolynomial<TypeParam> poly0_scaled_main,
        poly0_scaled_aux.SwitchRnsBasis(
            this->aux_moduli_, this->main_moduli_,
            this->rns_context_->AuxPrimeModulusCrtFactors(),
            this->rns_context_->AuxPrimeModulusComplementResidues(),
            this->rns_context_->AuxModulusResidues()));

    // Q/t * iNTT(ys) (mod P).
    ASSERT_OK_AND_ASSIGN(
        std::vector<RnsInt<TypeParam>> q_hat_mod_ps,
        this->rns_context_->MainPrimeModulusComplementResidues(level));
    ASSERT_OK_AND_ASSIGN(
        auto q_mod_ps,
        this->rns_context_->MainLeveledModulusAuxResidues(level));
    ASSERT_OK_AND_ASSIGN(
        RnsPolynomial<TypeParam> poly1_aux,
        poly1.SwitchRnsBasis(this->main_moduli_, this->aux_moduli_,
                             q_hat_inv_mod_qs, q_hat_mod_ps,
                             q_mod_ps.RnsRep()));

    // PQ/t^2 * iNTT(xs) * iNTT(ys) (mod Q).
    ASSERT_OK(poly0_scaled_main.ConvertToNttForm(this->main_moduli_));
    ASSERT_OK(poly1.ConvertToNttForm(this->main_moduli_));
    ASSERT_OK_AND_ASSIGN(RnsPolynomial<TypeParam> product_scaled_main,
                         poly0_scaled_main.Mul(poly1, this->main_moduli_));
    // PQ/t^2 * iNTT(xs) * iNTT(ys) (mod P).
    ASSERT_OK(poly0_scaled_aux.ConvertToNttForm(this->aux_moduli_));
    ASSERT_OK(poly1_aux.ConvertToNttForm(this->aux_moduli_));
    ASSERT_OK_AND_ASSIGN(RnsPolynomial<TypeParam> product_scaled_aux,
                         poly0_scaled_aux.Mul(poly1_aux, this->aux_moduli_));

    // // Q/t * iNTT(xs) * iNTT(ys) (mod Q).
    RnsPolynomial<TypeParam> product = product_scaled_main;
    ASSERT_OK(product.ScaleAndReduceRnsBasisInPlace(
        product_scaled_aux, t, this->main_moduli_, this->aux_moduli_,
        this->rns_context_->AuxPrimeModulusCrtFactors(),
        this->rns_context_->AuxPrimeModulusComplementResidues(),
        this->rns_context_->AuxModulusResidues(),
        this->rns_context_->AuxModulusInverseResidues()));

    // Now decode product and check results.
    ASSERT_OK_AND_ASSIGN(auto decoded,
                         encoder.DecodeBfv(product, this->main_moduli_));
    ASSERT_EQ(decoded.size(), num_coeffs);
    for (int i = 0; i < num_coeffs; ++i) {
      BigInteger expected = static_cast<BigInteger>(xs[i]);
      expected *= static_cast<BigInteger>(ys[i]);
      expected %= static_cast<BigInteger>(t);
      EXPECT_EQ(decoded[i], expected);
    }
  }
}

TYPED_TEST(FiniteFieldEncoderBfvTest,
           EncodeBfvIsMultiplicativeHomomorphicWithUnscaledPlaintext) {
  using Integer = typename TypeParam::Int;
  using BigInteger = typename TypeParam::BigInt;
  for (auto const& params :
       testing::GetBfvParametersForFiniteFieldEncoding<TypeParam>()) {
    // For correct decoding, we will need Q >= norm(poly0 * poly1), which can be
    // bounded by N * t^2 .
    Integer t = params.t;
    int plaintext_bits = ceil(log2(static_cast<double>(t)));
    int main_modulus_bits = 0;
    for (auto qi : params.qs) {
      main_modulus_bits += floor(log2(static_cast<double>(qi)));
    }
    if (main_modulus_bits < 2 * plaintext_bits + params.log_n) {
      continue;
    }

    this->SetUpRnsParameters(params);
    ASSERT_OK_AND_ASSIGN(auto encoder, FiniteFieldEncoder<TypeParam>::Create(
                                           this->rns_context_.get()));

    // Sample two vectors of messages mod t and encode them to RNS polynomials.
    int num_coeffs = 1 << this->rns_context_->LogN();
    std::vector<Integer> xs = testing::SampleMessages<Integer>(num_coeffs, t);
    std::vector<Integer> ys = testing::SampleMessages<Integer>(num_coeffs, t);
    ASSERT_OK_AND_ASSIGN(auto poly0, encoder.EncodeBfv(xs, this->main_moduli_));
    ASSERT_OK_AND_ASSIGN(auto poly1, encoder.EncodeBfv(ys, this->main_moduli_,
                                                       /*is_scaled=*/false));
    ASSERT_TRUE(poly0.IsNttForm());
    ASSERT_TRUE(poly1.IsNttForm());

    ASSERT_OK_AND_ASSIGN(RnsPolynomial<TypeParam> product,
                         poly0.Mul(poly1, this->main_moduli_));

    // Now decode product and check results.
    ASSERT_OK_AND_ASSIGN(auto decoded,
                         encoder.DecodeBfv(product, this->main_moduli_));
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
