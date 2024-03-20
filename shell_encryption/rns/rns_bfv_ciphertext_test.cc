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

#include "shell_encryption/rns/rns_bfv_ciphertext.h"

#include <cmath>
#include <memory>
#include <utility>
#include <vector>

#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include "absl/log/check.h"
#include "absl/status/status.h"
#include "absl/strings/string_view.h"
#include "shell_encryption/prng/single_thread_hkdf_prng.h"
#include "shell_encryption/rns/coefficient_encoder.h"
#include "shell_encryption/rns/finite_field_encoder.h"
#include "shell_encryption/rns/rns_context.h"
#include "shell_encryption/rns/rns_error_params.h"
#include "shell_encryption/rns/rns_modulus.h"
#include "shell_encryption/rns/rns_polynomial.h"
#include "shell_encryption/rns/rns_secret_key.h"
#include "shell_encryption/rns/testing/parameters.h"
#include "shell_encryption/rns/testing/testing_utils.h"
#include "shell_encryption/testing/parameters.h"
#include "shell_encryption/testing/status_matchers.h"
#include "shell_encryption/testing/status_testing.h"

namespace rlwe {
namespace {

using Prng = SingleThreadHkdfPrng;
using ::rlwe::testing::StatusIs;
using ::testing::HasSubstr;

constexpr int kVariance = 8;
constexpr absl::string_view kPrngSeed =
    "0123456789abcdef0123456789abcdef0123456789abcdef0123456789abcdef";

// Tests for homomorphic operations on BFV ciphertexts.
template <typename ModularInt>
class RnsBfvCiphertextTest : public ::testing::Test {
 protected:
  void SetUp() override {
    testing::RnsParameters<ModularInt> rns_params =
        testing::GetRnsParameters<ModularInt>();
    if (rns_params.t % 2 == 0) {
      rns_params.t++;
    }

    // Create the RNS context.
    auto rns_context = RnsContext<ModularInt>::CreateForBfv(
        rns_params.log_n, rns_params.qs, rns_params.ps, rns_params.t);
    CHECK(rns_context.ok());
    rns_context_ = std::make_unique<const RnsContext<ModularInt>>(
        std::move(rns_context.value()));
    main_moduli_ = rns_context_->MainPrimeModuli();

    // Create the coefficient encoder.
    auto coeff_encoder =
        CoefficientEncoder<ModularInt>::Create(rns_context_.get());
    CHECK(coeff_encoder.ok());
    coeff_encoder_ = std::make_unique<const CoefficientEncoder<ModularInt>>(
        std::move(coeff_encoder.value()));

    // Create the error parameters.
    int log_t = floor(std::log2(static_cast<double>(rns_params.t)));
    auto error_params = RnsErrorParams<ModularInt>::Create(
        rns_params.log_n, main_moduli_, /*aux_moduli=*/{}, log_t,
        sqrt(kVariance));
    CHECK(error_params.ok());
    error_params_ = std::make_unique<const RnsErrorParams<ModularInt>>(
        std::move(error_params.value()));

    // Create a PRNG for sampling secret key.
    auto prng = Prng::Create(kPrngSeed.substr(0, Prng::SeedLength()));
    CHECK(prng.ok());
    prng_ = std::move(prng.value());
  }

  // Sample a secret key from prng.
  absl::StatusOr<RnsRlweSecretKey<ModularInt>> SampleKey() const {
    return RnsRlweSecretKey<ModularInt>::Sample(rns_context_->LogN(), kVariance,
                                                rns_context_->MainPrimeModuli(),
                                                prng_.get());
  }

  std::unique_ptr<const RnsContext<ModularInt>> rns_context_;
  std::vector<const PrimeModulus<ModularInt>*> main_moduli_;
  std::unique_ptr<const CoefficientEncoder<ModularInt>> coeff_encoder_;
  std::unique_ptr<const RnsErrorParams<ModularInt>> error_params_;
  std::unique_ptr<Prng> prng_;
};

TYPED_TEST_SUITE(RnsBfvCiphertextTest, testing::ModularIntTypes);

TYPED_TEST(RnsBfvCiphertextTest, NegatedCiphertext) {
  using Integer = typename TypeParam::Int;
  ASSERT_OK_AND_ASSIGN(RnsRlweSecretKey<TypeParam> key, this->SampleKey());

  // Sample a plaintext polynomial and encrypt it.
  Integer t = this->rns_context_->PlaintextModulus();
  std::vector<Integer> messages =
      testing::SampleMessages(1 << this->rns_context_->LogN(), t);
  ASSERT_OK_AND_ASSIGN(
      RnsBfvCiphertext<TypeParam> ciphertext,
      key.EncryptBfv(messages, this->coeff_encoder_.get(),
                     this->error_params_.get(), this->prng_.get()));

  // Homomorphically negate the ciphertext.
  ASSERT_OK_AND_ASSIGN(RnsBfvCiphertext<TypeParam> ciphertext_neg,
                       ciphertext.Negate());

  // The negated ciphertext should decrypt to negated plaintext.
  ASSERT_OK_AND_ASSIGN(
      std::vector<Integer> dec_neg_messages,
      key.DecryptBfv(ciphertext_neg, this->coeff_encoder_.get()));
  ASSERT_EQ(dec_neg_messages.size(), messages.size());
  for (int i = 0; i < messages.size(); ++i) {
    Integer expected = (t - messages[i]) % t;  // -messages[i] (mod t).
    EXPECT_EQ(dec_neg_messages[i], expected);
  }

  // The sum of two ciphertexts should decrypt to zero.
  ASSERT_OK_AND_ASSIGN(RnsBfvCiphertext<TypeParam> ciphertext_sum,
                       ciphertext + ciphertext_neg);
  ASSERT_OK_AND_ASSIGN(
      std::vector<Integer> dec_sum_messages,
      key.DecryptBfv(ciphertext_sum, this->coeff_encoder_.get()));
  ASSERT_EQ(dec_sum_messages.size(), messages.size());
  for (int i = 0; i < dec_sum_messages.size(); ++i) {
    EXPECT_EQ(dec_sum_messages[i], 0);
  }
}

TYPED_TEST(RnsBfvCiphertextTest, HomomorphicAddition) {
  using Integer = typename TypeParam::Int;
  ASSERT_OK_AND_ASSIGN(RnsRlweSecretKey<TypeParam> key, this->SampleKey());
  const CoefficientEncoder<TypeParam>* encoder = this->coeff_encoder_.get();

  // Sample two plaintext polynomials and encrypt them.
  int num_coeffs = 1 << this->rns_context_->LogN();
  Integer t = this->rns_context_->PlaintextModulus();
  std::vector<Integer> messages0 = testing::SampleMessages(num_coeffs, t);
  std::vector<Integer> messages1 = testing::SampleMessages(num_coeffs, t);
  ASSERT_OK_AND_ASSIGN(
      RnsBfvCiphertext<TypeParam> ciphertext0,
      key.EncryptBfv(messages0, this->coeff_encoder_.get(),
                     this->error_params_.get(), this->prng_.get()));
  ASSERT_OK_AND_ASSIGN(
      RnsBfvCiphertext<TypeParam> ciphertext1,
      key.EncryptBfv(messages1, this->coeff_encoder_.get(),
                     this->error_params_.get(), this->prng_.get()));

  // Homomorphically add the two ciphertexts.
  ASSERT_OK_AND_ASSIGN(RnsBfvCiphertext<TypeParam> ciphertext_sum,
                       ciphertext0 + ciphertext1);

  // The ciphertext should decrypt to the sum of messages (mod t).
  ASSERT_OK_AND_ASSIGN(std::vector<Integer> dec_messages,
                       key.DecryptBfv(ciphertext_sum, encoder));
  ASSERT_EQ(dec_messages.size(), num_coeffs);
  for (int i = 0; i < num_coeffs; ++i) {
    Integer expected = (messages0[i] + messages1[i]) % t;
    EXPECT_EQ(dec_messages[i], expected);
  }
}

TYPED_TEST(RnsBfvCiphertextTest, HomomorphicAdditionPlaintext) {
  using Integer = typename TypeParam::Int;
  ASSERT_OK_AND_ASSIGN(RnsRlweSecretKey<TypeParam> key, this->SampleKey());
  const CoefficientEncoder<TypeParam>* encoder = this->coeff_encoder_.get();

  // Encrypt random plaintext messages.
  int num_coeffs = 1 << this->rns_context_->LogN();
  Integer t = this->rns_context_->PlaintextModulus();
  std::vector<Integer> messages0 = testing::SampleMessages(num_coeffs, t);
  std::vector<Integer> messages1 = testing::SampleMessages(num_coeffs, t);
  ASSERT_OK_AND_ASSIGN(
      RnsBfvCiphertext<TypeParam> ciphertext,
      key.EncryptBfv(messages0, this->coeff_encoder_.get(),
                     this->error_params_.get(), this->prng_.get()));

  // Encode `messages1` to a plaintext polynomial and homomorphically add it to
  // `ciphertext0`.
  ASSERT_OK_AND_ASSIGN(RnsPolynomial<TypeParam> plaintext,
                       encoder->EncodeBfv(messages1, this->main_moduli_));
  ASSERT_TRUE(plaintext.IsNttForm());
  ASSERT_OK_AND_ASSIGN(RnsBfvCiphertext<TypeParam> ciphertext_sum,
                       ciphertext + plaintext);

  // This ciphertext should also decrypt to the sum of messages (mod t).
  ASSERT_OK_AND_ASSIGN(std::vector<Integer> dec_messages,
                       key.DecryptBfv(ciphertext_sum, encoder));
  ASSERT_EQ(dec_messages.size(), num_coeffs);
  for (int i = 0; i < num_coeffs; ++i) {
    Integer expected = (messages0[i] + messages1[i]) % t;
    EXPECT_EQ(dec_messages[i], expected);
  }
}

TYPED_TEST(RnsBfvCiphertextTest, HomomorphicSubtraction) {
  using Integer = typename TypeParam::Int;
  ASSERT_OK_AND_ASSIGN(RnsRlweSecretKey<TypeParam> key, this->SampleKey());
  const CoefficientEncoder<TypeParam>* encoder = this->coeff_encoder_.get();

  // Sample two plaintext polynomials and encrypt them.
  int num_coeffs = 1 << this->rns_context_->LogN();
  Integer t = this->rns_context_->PlaintextModulus();
  std::vector<Integer> messages0 = testing::SampleMessages(num_coeffs, t);
  std::vector<Integer> messages1 = testing::SampleMessages(num_coeffs, t);
  ASSERT_OK_AND_ASSIGN(
      RnsBfvCiphertext<TypeParam> ciphertext0,
      key.EncryptBfv(messages0, this->coeff_encoder_.get(),
                     this->error_params_.get(), this->prng_.get()));
  ASSERT_OK_AND_ASSIGN(
      RnsBfvCiphertext<TypeParam> ciphertext1,
      key.EncryptBfv(messages1, this->coeff_encoder_.get(),
                     this->error_params_.get(), this->prng_.get()));

  // Homomorphically subtract the two ciphertexts.
  ASSERT_OK_AND_ASSIGN(RnsBfvCiphertext<TypeParam> ciphertext_diff,
                       ciphertext0 - ciphertext1);

  // The ciphertext should decrypt to the difference of messages (mod t).
  ASSERT_OK_AND_ASSIGN(std::vector<Integer> dec_messages,
                       key.DecryptBfv(ciphertext_diff, encoder));
  ASSERT_EQ(dec_messages.size(), num_coeffs);
  for (int i = 0; i < num_coeffs; ++i) {
    Integer expected = (t + messages0[i] - messages1[i]) % t;
    EXPECT_EQ(dec_messages[i], expected);
  }
}

TYPED_TEST(RnsBfvCiphertextTest, HomomorphicSubtractionPlaintext) {
  using Integer = typename TypeParam::Int;
  ASSERT_OK_AND_ASSIGN(RnsRlweSecretKey<TypeParam> key, this->SampleKey());
  const CoefficientEncoder<TypeParam>* encoder = this->coeff_encoder_.get();

  // Encrypt random plaintext messages.
  int num_coeffs = 1 << this->rns_context_->LogN();
  Integer t = this->rns_context_->PlaintextModulus();
  std::vector<Integer> messages0 = testing::SampleMessages(num_coeffs, t);
  std::vector<Integer> messages1 = testing::SampleMessages(num_coeffs, t);
  ASSERT_OK_AND_ASSIGN(
      RnsBfvCiphertext<TypeParam> ciphertext,
      key.EncryptBfv(messages0, this->coeff_encoder_.get(),
                     this->error_params_.get(), this->prng_.get()));

  // Encode `messages1` to a plaintext polynomial and homomorphically add it to
  // `ciphertext0`.
  ASSERT_OK_AND_ASSIGN(RnsPolynomial<TypeParam> plaintext,
                       encoder->EncodeBfv(messages1, this->main_moduli_));
  ASSERT_TRUE(plaintext.IsNttForm());
  ASSERT_OK_AND_ASSIGN(RnsBfvCiphertext<TypeParam> ciphertext_diff,
                       ciphertext - plaintext);

  // This ciphertext should also decrypt to the diff of messages (mod t).
  ASSERT_OK_AND_ASSIGN(std::vector<Integer> dec_messages,
                       key.DecryptBfv(ciphertext_diff, encoder));
  ASSERT_EQ(dec_messages.size(), num_coeffs);
  for (int i = 0; i < num_coeffs; ++i) {
    Integer expected = (t + messages0[i] - messages1[i]) % t;
    EXPECT_EQ(dec_messages[i], expected);
  }
}

TYPED_TEST(RnsBfvCiphertextTest, AbsorbFailsIfPlaintextIsNotNttForm) {
  ASSERT_OK_AND_ASSIGN(auto zero, RnsPolynomial<TypeParam>::CreateZero(
                                      this->rns_context_->LogN(),
                                      this->main_moduli_, /*is_ntt=*/false));
  RnsBfvCiphertext<TypeParam> ciphertext(
      /*components=*/{zero, zero}, this->main_moduli_,
      /*power_of_s=*/1, /*error=*/0, this->error_params_.get(),
      this->rns_context_.get());

  EXPECT_THAT(ciphertext.AbsorbInPlace(zero),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`plaintext` must be in NTT form")));
}

TYPED_TEST(RnsBfvCiphertextTest, AbsorbSimpleFailsIfPlaintextLevelMismatches) {
  // This test requires at least two prime moduli.
  if (this->main_moduli_.size() < 2) {
    GTEST_SKIP() << "Test requires at least two prime moduli.";
  }

  ASSERT_OK_AND_ASSIGN(auto zero,
                       RnsPolynomial<TypeParam>::CreateZero(
                           this->rns_context_->LogN(), this->main_moduli_));
  RnsBfvCiphertext<TypeParam> ciphertext(
      /*components=*/{zero, zero}, this->main_moduli_,
      /*power_of_s=*/1, /*error=*/0, this->error_params_.get(),
      this->rns_context_.get());

  ASSERT_OK_AND_ASSIGN(auto zero_at_level_zero,
                       RnsPolynomial<TypeParam>::CreateZero(
                           this->rns_context_->LogN(),
                           absl::MakeSpan(this->main_moduli_).subspan(0, 1)));
  EXPECT_THAT(ciphertext.AbsorbInPlaceSimple(zero_at_level_zero),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`plaintext` has a mismatched level")));
}

TYPED_TEST(RnsBfvCiphertextTest, AbsorbSimpleFailsIfPlaintextIsNotNttForm) {
  ASSERT_OK_AND_ASSIGN(auto zero, RnsPolynomial<TypeParam>::CreateZero(
                                      this->rns_context_->LogN(),
                                      this->main_moduli_, /*is_ntt=*/false));
  RnsBfvCiphertext<TypeParam> ciphertext(
      /*components=*/{zero, zero}, this->main_moduli_,
      /*power_of_s=*/1, /*error=*/0, this->error_params_.get(),
      this->rns_context_.get());

  EXPECT_THAT(ciphertext.AbsorbInPlaceSimple(zero),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`plaintext` must be in NTT form")));
}

TYPED_TEST(RnsBfvCiphertextTest, MulFailsIfLevelMismatches) {
  // This test requires at least two prime moduli.
  if (this->main_moduli_.size() < 2) {
    return;
  }

  ASSERT_OK_AND_ASSIGN(auto zero,
                       RnsPolynomial<TypeParam>::CreateZero(
                           this->rns_context_->LogN(), this->main_moduli_));
  RnsBfvCiphertext<TypeParam> ciphertext0(
      /*components=*/{zero, zero}, this->main_moduli_,
      /*power_of_s=*/1, /*error=*/0, this->error_params_.get(),
      this->rns_context_.get());

  // Create a ciphertext at the lowest level.
  std::vector<const PrimeModulus<TypeParam>*> moduli_at_level_zero = {
      this->main_moduli_[0]};
  ASSERT_OK_AND_ASSIGN(auto zero_at_level_zero,
                       RnsPolynomial<TypeParam>::CreateZero(
                           this->rns_context_->LogN(), moduli_at_level_zero));
  RnsBfvCiphertext<TypeParam> ciphertext1(
      /*components=*/{zero_at_level_zero, zero_at_level_zero},
      moduli_at_level_zero,
      /*power_of_s=*/1, /*error=*/0, this->error_params_.get(),
      this->rns_context_.get());

  EXPECT_THAT(ciphertext0.Mul(ciphertext1),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`that` has a mismatched level")));
}

TYPED_TEST(RnsBfvCiphertextTest, MulFailsIfPowerOfSMismatches) {
  ASSERT_OK_AND_ASSIGN(auto zero,
                       RnsPolynomial<TypeParam>::CreateZero(
                           this->rns_context_->LogN(), this->main_moduli_));
  RnsBfvCiphertext<TypeParam> ciphertext0(
      /*components=*/{zero, zero}, this->main_moduli_,
      /*power_of_s=*/1, /*error=*/0, this->error_params_.get(),
      this->rns_context_.get());
  RnsBfvCiphertext<TypeParam> ciphertext1(
      /*components=*/{zero, zero}, this->main_moduli_,
      /*power_of_s=*/3, /*error=*/0, this->error_params_.get(),
      this->rns_context_.get());

  // Mul() should fail as two ciphertexts have different PowerOfS.
  EXPECT_THAT(
      ciphertext0.Mul(ciphertext1),
      StatusIs(
          absl::StatusCode::kInvalidArgument,
          HasSubstr("Ciphertexts must be encrypted with the same key power")));
}

TYPED_TEST(RnsBfvCiphertextTest, ModReduceFailsIfAtLevelZero) {
  // Create a ciphertext with a single prime modulus.
  std::vector<const PrimeModulus<TypeParam>*> moduli_at_level_zero = {
      this->main_moduli_[0]};
  ASSERT_OK_AND_ASSIGN(auto zero,
                       RnsPolynomial<TypeParam>::CreateZero(
                           this->rns_context_->LogN(), moduli_at_level_zero));
  RnsBfvCiphertext<TypeParam> ciphertext(
      /*components=*/{zero, zero}, moduli_at_level_zero,
      /*power_of_s=*/5, /*error=*/0, this->error_params_.get(),
      this->rns_context_.get());

  EXPECT_THAT(ciphertext.ModReduce(),
              StatusIs(absl::StatusCode::kFailedPrecondition,
                       HasSubstr("Cannot perform ModReduce with insufficient "
                                 "number of prime moduli")));
}

TYPED_TEST(RnsBfvCiphertextTest, ModReducedCiphertextDecrypts) {
  using Integer = typename TypeParam::Int;

  if (this->main_moduli_.size() < 2) {
    // There is only one prime modulus in the moduli chain, so we cannot perform
    // modulus reduction further.
    return;
  }

  ASSERT_OK_AND_ASSIGN(RnsRlweSecretKey<TypeParam> key, this->SampleKey());

  // Sample a plaintext polynomial and encrypt it wrt the full moduli chain.
  Integer t = this->rns_context_->PlaintextModulus();
  std::vector<Integer> messages =
      testing::SampleMessages(1 << this->rns_context_->LogN(), t);
  ASSERT_OK_AND_ASSIGN(
      RnsBfvCiphertext<TypeParam> ciphertext,
      key.EncryptBfv(messages, this->coeff_encoder_.get(),
                     this->error_params_.get(), this->prng_.get()));

  // ModReduce the ciphertext and get round(Q/q_L * ciphertext) mod (Q/q_L).
  int level = ciphertext.Level();
  int degree = ciphertext.Degree();
  ASSERT_OK(ciphertext.ModReduce());
  ASSERT_EQ(ciphertext.Level(), level - 1);  // level should be reduced by 1
  ASSERT_EQ(ciphertext.Degree(), degree);    // degree should not change

  // We also need to modulus reduce the secret key.
  ASSERT_OK(key.ModReduce());
  ASSERT_EQ(key.Level(), level - 1);  // key level should match the ciphertext

  ASSERT_OK_AND_ASSIGN(std::vector<Integer> dec_messages,
                       key.DecryptBfv(ciphertext, this->coeff_encoder_.get()));
  EXPECT_EQ(dec_messages, messages);
}

////////////////////////////////////////////////////////////////////////////////
// Homomorphic operations with finite field packed encoding
////////////////////////////////////////////////////////////////////////////////

template <typename ModularInt>
class RnsBfvCiphertextPackedTest : public ::testing::Test {
 protected:
  void SetUp() override {
    // Create a PRNG for sampling secret key.
    auto prng = Prng::Create(kPrngSeed.substr(0, Prng::SeedLength()));
    CHECK(prng.ok());
    prng_ = std::move(prng.value());
  }

  void SetUpBfvRnsParameters(const testing::RnsParameters<ModularInt>& params) {
    // Use parameters suitable for finite field encoding.
    auto rns_context = RnsContext<ModularInt>::CreateForBfvFiniteFieldEncoding(
        params.log_n, params.qs, params.ps, params.t);
    CHECK(rns_context.ok());
    rns_context_ = std::make_unique<const RnsContext<ModularInt>>(
        std::move(rns_context.value()));
    main_moduli_ = rns_context_->MainPrimeModuli();
    aux_moduli_ = rns_context_->AuxPrimeModuli();

    // Create the error parameters.
    int log_t = floor(std::log2(static_cast<double>(params.t)));
    auto error_params = RnsErrorParams<ModularInt>::Create(
        params.log_n, main_moduli_, /*aux_moduli=*/{}, log_t, sqrt(kVariance));
    CHECK(error_params.ok());
    error_params_ = std::make_unique<const RnsErrorParams<ModularInt>>(
        std::move(error_params.value()));

    // Create a finite field encoder.
    auto encoder = FiniteFieldEncoder<ModularInt>::Create(rns_context_.get());
    CHECK(encoder.ok());
    encoder_ = std::make_unique<const FiniteFieldEncoder<ModularInt>>(
        std::move(encoder.value()));
  }

  void SetUpDefaultBfvRnsParameters() {
    auto all_params =
        testing::GetBfvParametersForFiniteFieldEncoding<ModularInt>();
    CHECK(!all_params.empty());
    SetUpBfvRnsParameters(all_params[0]);
  }

  // Sample a secret key from prng.
  absl::StatusOr<RnsRlweSecretKey<ModularInt>> SampleKey() const {
    return RnsRlweSecretKey<ModularInt>::Sample(rns_context_->LogN(), kVariance,
                                                rns_context_->MainPrimeModuli(),
                                                prng_.get());
  }

  std::unique_ptr<const RnsContext<ModularInt>> rns_context_;
  std::vector<const PrimeModulus<ModularInt>*> main_moduli_;
  std::vector<const PrimeModulus<ModularInt>*> aux_moduli_;
  std::unique_ptr<const FiniteFieldEncoder<ModularInt>> encoder_;
  std::unique_ptr<const RnsErrorParams<ModularInt>> error_params_;
  std::unique_ptr<Prng> prng_;
};

TYPED_TEST_SUITE(RnsBfvCiphertextPackedTest,
                 testing::ModularIntTypesForFiniteFieldEncoding);

TYPED_TEST(RnsBfvCiphertextPackedTest, NegatedCiphertext) {
  using Integer = typename TypeParam::Int;
  for (auto const& params :
       testing::GetRnsParametersForFiniteFieldEncoding<TypeParam>()) {
    this->SetUpBfvRnsParameters(params);

    ASSERT_OK_AND_ASSIGN(RnsRlweSecretKey<TypeParam> key, this->SampleKey());

    // Sample a plaintext polynomial and encrypt it.
    Integer t = this->rns_context_->PlaintextModulus();
    std::vector<Integer> messages =
        testing::SampleMessages(1 << this->rns_context_->LogN(), t);
    ASSERT_OK_AND_ASSIGN(
        RnsBfvCiphertext<TypeParam> ciphertext,
        key.EncryptBfv(messages, this->encoder_.get(),
                       this->error_params_.get(), this->prng_.get()));

    // Homomorphically negate the ciphertext.
    ASSERT_OK_AND_ASSIGN(RnsBfvCiphertext<TypeParam> ciphertext_neg,
                         ciphertext.Negate());

    // The negated ciphertext should decrypt to negated plaintext.
    ASSERT_OK_AND_ASSIGN(std::vector<Integer> dec_neg_messages,
                         key.DecryptBfv(ciphertext_neg, this->encoder_.get()));
    ASSERT_EQ(dec_neg_messages.size(), messages.size());
    for (int i = 0; i < messages.size(); ++i) {
      Integer expected = (t - messages[i]) % t;  // -messages[i] (mod t).
      EXPECT_EQ(dec_neg_messages[i], expected);
    }
  }
}

TYPED_TEST(RnsBfvCiphertextPackedTest, HomomorphicAddition) {
  using Integer = typename TypeParam::Int;
  for (auto const& params :
       testing::GetRnsParametersForFiniteFieldEncoding<TypeParam>()) {
    this->SetUpBfvRnsParameters(params);

    ASSERT_OK_AND_ASSIGN(RnsRlweSecretKey<TypeParam> key, this->SampleKey());

    // Sample two plaintext polynomials and encrypt them.
    int num_coeffs = 1 << this->rns_context_->LogN();
    Integer t = this->rns_context_->PlaintextModulus();
    std::vector<Integer> messages0 = testing::SampleMessages(num_coeffs, t);
    std::vector<Integer> messages1 = testing::SampleMessages(num_coeffs, t);
    ASSERT_OK_AND_ASSIGN(
        RnsBfvCiphertext<TypeParam> ciphertext0,
        key.EncryptBfv(messages0, this->encoder_.get(),
                       this->error_params_.get(), this->prng_.get()));
    ASSERT_OK_AND_ASSIGN(
        RnsBfvCiphertext<TypeParam> ciphertext1,
        key.EncryptBfv(messages1, this->encoder_.get(),
                       this->error_params_.get(), this->prng_.get()));

    // Homomorphically add the two ciphertexts.
    ASSERT_OK_AND_ASSIGN(RnsBfvCiphertext<TypeParam> ciphertext_sum,
                         ciphertext0 + ciphertext1);

    // The ciphertext should decrypt to the sum of messages (mod t).
    ASSERT_OK_AND_ASSIGN(std::vector<Integer> dec_messages,
                         key.DecryptBfv(ciphertext_sum, this->encoder_.get()));
    ASSERT_EQ(dec_messages.size(), num_coeffs);
    for (int i = 0; i < num_coeffs; ++i) {
      Integer expected = (messages0[i] + messages1[i]) % t;
      EXPECT_EQ(dec_messages[i], expected);
    }
  }
}

TYPED_TEST(RnsBfvCiphertextPackedTest, HomomorphicAdditionPlaintext) {
  using Integer = typename TypeParam::Int;
  for (auto const& params :
       testing::GetRnsParametersForFiniteFieldEncoding<TypeParam>()) {
    this->SetUpBfvRnsParameters(params);
    ASSERT_OK_AND_ASSIGN(RnsRlweSecretKey<TypeParam> key, this->SampleKey());
    const FiniteFieldEncoder<TypeParam>* encoder = this->encoder_.get();

    // Encrypt random plaintext messages.
    int num_coeffs = 1 << this->rns_context_->LogN();
    Integer t = this->rns_context_->PlaintextModulus();
    std::vector<Integer> messages0 = testing::SampleMessages(num_coeffs, t);
    std::vector<Integer> messages1 = testing::SampleMessages(num_coeffs, t);
    ASSERT_OK_AND_ASSIGN(
        RnsBfvCiphertext<TypeParam> ciphertext,
        key.EncryptBfv(messages0, encoder, this->error_params_.get(),
                       this->prng_.get()));

    // Encode `messages1` to a plaintext polynomial and homomorphically add it
    // to `ciphertext0`.
    ASSERT_OK_AND_ASSIGN(RnsPolynomial<TypeParam> plaintext,
                         encoder->EncodeBfv(messages1, this->main_moduli_));
    ASSERT_TRUE(plaintext.IsNttForm());
    ASSERT_OK_AND_ASSIGN(RnsBfvCiphertext<TypeParam> ciphertext_sum,
                         ciphertext + plaintext);

    // This ciphertext should also decrypt to the sum of messages (mod t).
    ASSERT_OK_AND_ASSIGN(std::vector<Integer> dec_messages,
                         key.DecryptBfv(ciphertext_sum, encoder));
    ASSERT_EQ(dec_messages.size(), num_coeffs);
    for (int i = 0; i < num_coeffs; ++i) {
      Integer expected = (messages0[i] + messages1[i]) % t;
      EXPECT_EQ(dec_messages[i], expected);
    }
  }
}

TYPED_TEST(RnsBfvCiphertextPackedTest, HomomorphicSubtraction) {
  using Integer = typename TypeParam::Int;
  for (auto const& params :
       testing::GetRnsParametersForFiniteFieldEncoding<TypeParam>()) {
    this->SetUpBfvRnsParameters(params);
    ASSERT_OK_AND_ASSIGN(RnsRlweSecretKey<TypeParam> key, this->SampleKey());

    // Sample two plaintext polynomials and encrypt them.
    int num_coeffs = 1 << this->rns_context_->LogN();
    Integer t = this->rns_context_->PlaintextModulus();
    std::vector<Integer> messages0 = testing::SampleMessages(num_coeffs, t);
    std::vector<Integer> messages1 = testing::SampleMessages(num_coeffs, t);
    ASSERT_OK_AND_ASSIGN(
        RnsBfvCiphertext<TypeParam> ciphertext0,
        key.EncryptBfv(messages0, this->encoder_.get(),
                       this->error_params_.get(), this->prng_.get()));
    ASSERT_OK_AND_ASSIGN(
        RnsBfvCiphertext<TypeParam> ciphertext1,
        key.EncryptBfv(messages1, this->encoder_.get(),
                       this->error_params_.get(), this->prng_.get()));

    // Homomorphically subtract the two ciphertexts.
    ASSERT_OK_AND_ASSIGN(RnsBfvCiphertext<TypeParam> ciphertext_diff,
                         ciphertext0 - ciphertext1);

    // The ciphertext should decrypt to the difference of messages (mod t).
    ASSERT_OK_AND_ASSIGN(std::vector<Integer> dec_messages,
                         key.DecryptBfv(ciphertext_diff, this->encoder_.get()));
    ASSERT_EQ(dec_messages.size(), num_coeffs);
    for (int i = 0; i < num_coeffs; ++i) {
      Integer expected = (t + messages0[i] - messages1[i]) % t;
      EXPECT_EQ(dec_messages[i], expected);
    }
  }
}

TYPED_TEST(RnsBfvCiphertextPackedTest, HomomorphicSubtractionPlaintext) {
  using Integer = typename TypeParam::Int;
  for (auto const& params :
       testing::GetRnsParametersForFiniteFieldEncoding<TypeParam>()) {
    this->SetUpBfvRnsParameters(params);
    ASSERT_OK_AND_ASSIGN(RnsRlweSecretKey<TypeParam> key, this->SampleKey());
    const FiniteFieldEncoder<TypeParam>* encoder = this->encoder_.get();

    // Encrypt random plaintext messages.
    int num_coeffs = 1 << this->rns_context_->LogN();
    Integer t = this->rns_context_->PlaintextModulus();
    std::vector<Integer> messages0 = testing::SampleMessages(num_coeffs, t);
    std::vector<Integer> messages1 = testing::SampleMessages(num_coeffs, t);
    ASSERT_OK_AND_ASSIGN(
        RnsBfvCiphertext<TypeParam> ciphertext,
        key.EncryptBfv(messages0, this->encoder_.get(),
                       this->error_params_.get(), this->prng_.get()));

    // Encode `messages1` to a plaintext polynomial and homomorphically add it
    // to `ciphertext0`.
    ASSERT_OK_AND_ASSIGN(RnsPolynomial<TypeParam> plaintext,
                         encoder->EncodeBfv(messages1, this->main_moduli_));
    ASSERT_TRUE(plaintext.IsNttForm());
    ASSERT_OK_AND_ASSIGN(RnsBfvCiphertext<TypeParam> ciphertext_diff,
                         ciphertext - plaintext);

    // This ciphertext should also decrypt to the diff of messages (mod t).
    ASSERT_OK_AND_ASSIGN(std::vector<Integer> dec_messages,
                         key.DecryptBfv(ciphertext_diff, encoder));
    ASSERT_EQ(dec_messages.size(), num_coeffs);
    for (int i = 0; i < num_coeffs; ++i) {
      Integer expected = (t + messages0[i] - messages1[i]) % t;
      EXPECT_EQ(dec_messages[i], expected);
    }
  }
}

TYPED_TEST(RnsBfvCiphertextPackedTest, HomomorphicMulWithScalar) {
  using Integer = typename TypeParam::Int;

  // Three random scalar multiplications should be sufficient.
  constexpr int k_num_scalars = 3;

  for (auto const& params :
       testing::GetRnsParametersForFiniteFieldEncoding<TypeParam>()) {
    int plaintext_bits = ceil(log2(static_cast<double>(params.t)));
    int main_modulus_bits = 0;
    for (auto qi : params.qs) {
      main_modulus_bits += floor(log2(static_cast<double>(qi)));
    }
    if (main_modulus_bits < 2 * plaintext_bits) {
      GTEST_SKIP() << "Insufficient modulus for multiplication.";
    }

    this->SetUpBfvRnsParameters(params);
    ASSERT_OK_AND_ASSIGN(RnsRlweSecretKey<TypeParam> key, this->SampleKey());

    int log_n = this->rns_context_->LogN();
    Integer t = this->rns_context_->PlaintextModulus();
    std::vector<Integer> messages = testing::SampleMessages(1 << log_n, t);
    ASSERT_OK_AND_ASSIGN(RnsBfvCiphertext<TypeParam> ciphertext,
                         key.template EncryptBfv<FiniteFieldEncoder<TypeParam>>(
                             messages, this->encoder_.get(),
                             this->error_params_.get(), this->prng_.get()));

    // Multiply the ciphertext with a plaintext scalar.
    std::vector<Integer> scalars = testing::SampleMessages(k_num_scalars, t);
    for (int j = 0; j < k_num_scalars; ++j) {
      ASSERT_OK_AND_ASSIGN(auto ciphertext_product, ciphertext* scalars[j]);
      ASSERT_OK_AND_ASSIGN(
          std::vector<Integer> dec_messages,
          key.template DecryptBfv<FiniteFieldEncoder<TypeParam>>(
              ciphertext_product, this->encoder_.get()));
      for (int i = 0; i < (1 << log_n); ++i) {
        EXPECT_EQ(dec_messages[i], (messages[i] * scalars[j]) % t);
      }
    }
  }
}

TYPED_TEST(RnsBfvCiphertextPackedTest, HomomorphicMulWithPlaintext) {
  using Integer = typename TypeParam::Int;

  for (auto const& params :
       testing::GetBfvParametersForFiniteFieldEncoding<TypeParam>()) {
    int plaintext_bits = ceil(log2(static_cast<double>(params.t)));
    int main_modulus_bits = 0;
    for (auto qi : params.qs) {
      main_modulus_bits += floor(log2(static_cast<double>(qi)));
    }
    if (main_modulus_bits < 2 * plaintext_bits + params.log_n) {
      GTEST_SKIP() << "Insufficient modulus for multiplication.";
    }
    this->SetUpBfvRnsParameters(params);
    ASSERT_OK_AND_ASSIGN(RnsRlweSecretKey<TypeParam> key, this->SampleKey());

    int log_n = this->rns_context_->LogN();
    Integer t = this->rns_context_->PlaintextModulus();
    const FiniteFieldEncoder<TypeParam>* encoder = this->encoder_.get();

    std::vector<Integer> messages0 = testing::SampleMessages(1 << log_n, t);
    std::vector<Integer> messages1 = testing::SampleMessages(1 << log_n, t);
    ASSERT_OK_AND_ASSIGN(
        RnsBfvCiphertext<TypeParam> ciphertext,
        key.template EncryptBfv<FiniteFieldEncoder<TypeParam>>(
            messages0, encoder, this->error_params_.get(), this->prng_.get()));
    ASSERT_OK_AND_ASSIGN(auto plaintext,
                         encoder->EncodeBfv(messages1, this->main_moduli_));

    ASSERT_OK_AND_ASSIGN(auto ciphertext_product, ciphertext* plaintext);
    ASSERT_OK_AND_ASSIGN(std::vector<Integer> dec_messages,
                         key.template DecryptBfv<FiniteFieldEncoder<TypeParam>>(
                             ciphertext_product, encoder));
    for (int i = 0; i < (1 << log_n); ++i) {
      EXPECT_EQ(dec_messages[i], (messages0[i] * messages1[i]) % t);
    }
  }
}

// Tests `AbsorbSimple` which multiplies with a plaintext that encodes the
// messages without scaling them up.
TYPED_TEST(RnsBfvCiphertextPackedTest, HomomorphicMulWithUnscaledPlaintext) {
  using Integer = typename TypeParam::Int;

  for (auto const& params :
       testing::GetBfvParametersForFiniteFieldEncoding<TypeParam>()) {
    int plaintext_bits = ceil(log2(static_cast<double>(params.t)));
    int main_modulus_bits = 0;
    for (auto qi : params.qs) {
      main_modulus_bits += floor(log2(static_cast<double>(qi)));
    }
    if (main_modulus_bits < 2 * plaintext_bits + params.log_n) {
      GTEST_SKIP() << "Insufficient modulus for multiplication.";
    }
    this->SetUpBfvRnsParameters(params);
    ASSERT_OK_AND_ASSIGN(RnsRlweSecretKey<TypeParam> key, this->SampleKey());

    int log_n = this->rns_context_->LogN();
    Integer t = this->rns_context_->PlaintextModulus();
    const FiniteFieldEncoder<TypeParam>* encoder = this->encoder_.get();

    std::vector<Integer> messages0 = testing::SampleMessages(1 << log_n, t);
    std::vector<Integer> messages1 = testing::SampleMessages(1 << log_n, t);
    ASSERT_OK_AND_ASSIGN(
        RnsBfvCiphertext<TypeParam> ciphertext,
        key.template EncryptBfv<FiniteFieldEncoder<TypeParam>>(
            messages0, encoder, this->error_params_.get(), this->prng_.get()));
    ASSERT_OK_AND_ASSIGN(
        auto plaintext,
        encoder->EncodeBfv(messages1, this->main_moduli_, /*is_scaled=*/false));

    ASSERT_OK_AND_ASSIGN(auto ciphertext_product,
                         ciphertext.AbsorbSimple(plaintext));
    ASSERT_OK_AND_ASSIGN(std::vector<Integer> dec_messages,
                         key.template DecryptBfv<FiniteFieldEncoder<TypeParam>>(
                             ciphertext_product, encoder));
    for (int i = 0; i < (1 << log_n); ++i) {
      EXPECT_EQ(dec_messages[i], (messages0[i] * messages1[i]) % t);
    }
  }
}

TYPED_TEST(RnsBfvCiphertextPackedTest, HomomorphicFusedAbsorbAdd) {
  using Integer = typename TypeParam::Int;

  for (auto const& params :
       testing::GetBfvParametersForFiniteFieldEncoding<TypeParam>()) {
    int plaintext_bits = ceil(log2(static_cast<double>(params.t)));
    int main_modulus_bits = 0;
    for (auto qi : params.qs) {
      main_modulus_bits += floor(log2(static_cast<double>(qi)));
    }
    if (main_modulus_bits < 2 * plaintext_bits + params.log_n) {
      GTEST_SKIP() << "Insufficient modulus for multiplication.";
    }
    this->SetUpBfvRnsParameters(params);
    ASSERT_OK_AND_ASSIGN(RnsRlweSecretKey<TypeParam> key, this->SampleKey());

    int log_n = this->rns_context_->LogN();
    Integer t = this->rns_context_->PlaintextModulus();
    const FiniteFieldEncoder<TypeParam>* encoder = this->encoder_.get();

    std::vector<Integer> messages0 = testing::SampleMessages(1 << log_n, t);
    std::vector<Integer> messages1 = testing::SampleMessages(1 << log_n, t);
    std::vector<Integer> messages2 =
        testing::SampleMessages<Integer>(1 << log_n, 2);
    ASSERT_OK_AND_ASSIGN(
        RnsBfvCiphertext<TypeParam> ciphertext0,
        key.template EncryptBfv<FiniteFieldEncoder<TypeParam>>(
            messages0, encoder, this->error_params_.get(), this->prng_.get()));
    ASSERT_OK_AND_ASSIGN(
        RnsBfvCiphertext<TypeParam> ciphertext1,
        key.template EncryptBfv<FiniteFieldEncoder<TypeParam>>(
            messages1, encoder, this->error_params_.get(), this->prng_.get()));
    ASSERT_OK_AND_ASSIGN(
        auto plaintext,
        encoder->EncodeBfv(messages2, this->main_moduli_, /*is_scaled=*/false));

    ASSERT_OK(ciphertext0.FusedAbsorbAddInPlace(ciphertext1, plaintext));
    ASSERT_OK_AND_ASSIGN(std::vector<Integer> dec_messages,
                         key.template DecryptBfv<FiniteFieldEncoder<TypeParam>>(
                             ciphertext0, encoder));
    for (int i = 0; i < (1 << log_n); ++i) {
      Integer expected = (messages0[i] + messages1[i] * messages2[i]) % t;
      EXPECT_EQ(dec_messages[i], expected);
    }
  }
}

TYPED_TEST(RnsBfvCiphertextPackedTest, HomomorphicFusedAbsorbAddWithoutPad) {
  using Integer = typename TypeParam::Int;

  for (auto const& params :
       testing::GetBfvParametersForFiniteFieldEncoding<TypeParam>()) {
    int plaintext_bits = ceil(log2(static_cast<double>(params.t)));
    int main_modulus_bits = 0;
    for (auto qi : params.qs) {
      main_modulus_bits += floor(log2(static_cast<double>(qi)));
    }
    if (main_modulus_bits < 2 * plaintext_bits + params.log_n) {
      GTEST_SKIP() << "Insufficient modulus for multiplication.";
    }
    this->SetUpBfvRnsParameters(params);
    ASSERT_OK_AND_ASSIGN(RnsRlweSecretKey<TypeParam> key, this->SampleKey());

    int log_n = this->rns_context_->LogN();
    Integer t = this->rns_context_->PlaintextModulus();
    const FiniteFieldEncoder<TypeParam>* encoder = this->encoder_.get();

    std::vector<Integer> messages0 = testing::SampleMessages(1 << log_n, t);
    std::vector<Integer> messages1 = testing::SampleMessages(1 << log_n, t);
    std::vector<Integer> messages2 =
        testing::SampleMessages<Integer>(1 << log_n, 2);
    ASSERT_OK_AND_ASSIGN(
        RnsBfvCiphertext<TypeParam> ciphertext0,
        key.template EncryptBfv<FiniteFieldEncoder<TypeParam>>(
            messages0, encoder, this->error_params_.get(), this->prng_.get()));
    ASSERT_OK_AND_ASSIGN(
        RnsBfvCiphertext<TypeParam> ciphertext1,
        key.template EncryptBfv<FiniteFieldEncoder<TypeParam>>(
            messages1, encoder, this->error_params_.get(), this->prng_.get()));
    ASSERT_OK_AND_ASSIGN(
        auto plaintext,
        encoder->EncodeBfv(messages2, this->main_moduli_, /*is_scaled=*/false));

    // Precompute the "a" part of the resulting ciphertext.
    ASSERT_OK_AND_ASSIGN(RnsPolynomial<TypeParam> pad0,
                         ciphertext0.Component(1));
    ASSERT_OK_AND_ASSIGN(RnsPolynomial<TypeParam> pad1,
                         ciphertext1.Component(1));
    ASSERT_OK(pad0.FusedMulAddInPlace(pad1, plaintext, this->main_moduli_));

    // Fused absorb add without updating the "a" part.
    ASSERT_OK(
        ciphertext0.FusedAbsorbAddInPlaceWithoutPad(ciphertext1, plaintext));

    // Now update the "a" part.
    ASSERT_OK(ciphertext0.SetPadComponent(std::move(pad0)));

    ASSERT_OK_AND_ASSIGN(std::vector<Integer> dec_messages,
                         key.template DecryptBfv<FiniteFieldEncoder<TypeParam>>(
                             ciphertext0, encoder));
    for (int i = 0; i < (1 << log_n); ++i) {
      Integer expected = (messages0[i] + messages1[i] * messages2[i]) % t;
      EXPECT_EQ(dec_messages[i], expected);
    }
  }
}

TYPED_TEST(RnsBfvCiphertextPackedTest, HomomorphicMulWithCiphertext) {
  using Integer = typename TypeParam::Int;

  for (auto const& params :
       testing::GetBfvParametersForFiniteFieldEncoding<TypeParam>()) {
    int plaintext_bits = ceil(log2(static_cast<double>(params.t)));
    int main_modulus_bits = 0;
    for (auto qi : params.qs) {
      main_modulus_bits += floor(log2(static_cast<double>(qi)));
    }
    if (main_modulus_bits < 2 * plaintext_bits + params.log_n) {
      GTEST_SKIP() << "Insufficient modulus for multiplication.";
    }
    this->SetUpBfvRnsParameters(params);
    ASSERT_OK_AND_ASSIGN(RnsRlweSecretKey<TypeParam> key, this->SampleKey());

    int log_n = this->rns_context_->LogN();
    Integer t = this->rns_context_->PlaintextModulus();
    const FiniteFieldEncoder<TypeParam>* encoder = this->encoder_.get();

    std::vector<Integer> messages0 = testing::SampleMessages(1 << log_n, t);
    std::vector<Integer> messages1 = testing::SampleMessages(1 << log_n, t);
    ASSERT_OK_AND_ASSIGN(
        RnsBfvCiphertext<TypeParam> ciphertext0,
        key.template EncryptBfv<FiniteFieldEncoder<TypeParam>>(
            messages0, encoder, this->error_params_.get(), this->prng_.get()));
    ASSERT_OK_AND_ASSIGN(
        RnsBfvCiphertext<TypeParam> ciphertext1,
        key.template EncryptBfv<FiniteFieldEncoder<TypeParam>>(
            messages1, encoder, this->error_params_.get(), this->prng_.get()));

    ASSERT_OK_AND_ASSIGN(auto ciphertext_product, ciphertext0* ciphertext1);
    ASSERT_OK_AND_ASSIGN(std::vector<Integer> dec_messages,
                         key.template DecryptBfv<FiniteFieldEncoder<TypeParam>>(
                             ciphertext_product, encoder));
    for (int i = 0; i < (1 << log_n); ++i) {
      EXPECT_EQ(dec_messages[i], (messages0[i] * messages1[i]) % t);
    }
  }
}

}  // namespace
}  // namespace rlwe
