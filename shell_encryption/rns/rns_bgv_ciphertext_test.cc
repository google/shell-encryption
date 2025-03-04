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

#include "shell_encryption/rns/rns_bgv_ciphertext.h"

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
#include "shell_encryption/rns/rns_integer.h"
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

// Tests for homomorphic operations on BGV ciphertexts.
template <typename ModularInt>
class RnsBgvCiphertextTest : public ::testing::Test {
 protected:
  void SetUp() override {
    testing::RnsParameters<ModularInt> rns_params =
        testing::GetRnsParameters<ModularInt>();

    // Create the RNS context.
    auto rns_context = RnsContext<ModularInt>::Create(
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

TYPED_TEST_SUITE(RnsBgvCiphertextTest, testing::ModularIntTypes);

TYPED_TEST(RnsBgvCiphertextTest, NegatedCiphertext) {
  using Integer = typename TypeParam::Int;
  ASSERT_OK_AND_ASSIGN(RnsRlweSecretKey<TypeParam> key, this->SampleKey());

  // Sample a plaintext polynomial and encrypt it.
  Integer t = this->rns_context_->PlaintextModulus();
  std::vector<Integer> messages =
      testing::SampleMessages(1 << this->rns_context_->LogN(), t);
  ASSERT_OK_AND_ASSIGN(
      RnsBgvCiphertext<TypeParam> ciphertext,
      key.EncryptBgv(messages, this->coeff_encoder_.get(),
                     this->error_params_.get(), this->prng_.get()));

  // Homomorphically negate the ciphertext.
  ASSERT_OK_AND_ASSIGN(RnsBgvCiphertext<TypeParam> ciphertext_neg,
                       ciphertext.Negate());

  // The negated ciphertext should decrypt to negated plaintext.
  ASSERT_OK_AND_ASSIGN(
      std::vector<Integer> dec_neg_messages,
      key.DecryptBgv(ciphertext_neg, this->coeff_encoder_.get()));
  ASSERT_EQ(dec_neg_messages.size(), messages.size());
  for (int i = 0; i < messages.size(); ++i) {
    Integer expected = (t - messages[i]) % t;  // -messages[i] (mod t).
    EXPECT_EQ(dec_neg_messages[i], expected);
  }

  // The sum of two ciphertexts should decrypt to zero.
  ASSERT_OK_AND_ASSIGN(RnsBgvCiphertext<TypeParam> ciphertext_sum,
                       ciphertext + ciphertext_neg);
  ASSERT_OK_AND_ASSIGN(
      std::vector<Integer> dec_sum_messages,
      key.DecryptBgv(ciphertext_sum, this->coeff_encoder_.get()));
  ASSERT_EQ(dec_sum_messages.size(), messages.size());
  for (int i = 0; i < dec_sum_messages.size(); ++i) {
    EXPECT_EQ(dec_sum_messages[i], 0);
  }
}

TYPED_TEST(RnsBgvCiphertextTest, HomomorphicAddition) {
  using Integer = typename TypeParam::Int;
  ASSERT_OK_AND_ASSIGN(RnsRlweSecretKey<TypeParam> key, this->SampleKey());
  const CoefficientEncoder<TypeParam>* encoder = this->coeff_encoder_.get();

  // Sample two plaintext polynomials and encrypt them.
  int num_coeffs = 1 << this->rns_context_->LogN();
  Integer t = this->rns_context_->PlaintextModulus();
  std::vector<Integer> messages0 = testing::SampleMessages(num_coeffs, t);
  std::vector<Integer> messages1 = testing::SampleMessages(num_coeffs, t);
  ASSERT_OK_AND_ASSIGN(
      RnsBgvCiphertext<TypeParam> ciphertext0,
      key.EncryptBgv(messages0, this->coeff_encoder_.get(),
                     this->error_params_.get(), this->prng_.get()));
  ASSERT_OK_AND_ASSIGN(
      RnsBgvCiphertext<TypeParam> ciphertext1,
      key.EncryptBgv(messages1, this->coeff_encoder_.get(),
                     this->error_params_.get(), this->prng_.get()));

  // Homomorphically add the two ciphertexts.
  ASSERT_OK_AND_ASSIGN(RnsBgvCiphertext<TypeParam> ciphertext_sum,
                       ciphertext0 + ciphertext1);

  // The ciphertext should decrypt to the sum of messages (mod t).
  ASSERT_OK_AND_ASSIGN(std::vector<Integer> dec_sum_messages,
                       key.DecryptBgv(ciphertext_sum, encoder));
  ASSERT_EQ(dec_sum_messages.size(), num_coeffs);
  for (int i = 0; i < num_coeffs; ++i) {
    Integer expected = (messages0[i] + messages1[i]) % t;
    EXPECT_EQ(dec_sum_messages[i], expected);
  }

  // Encode `messages1` to a plaintext polynomial and homomorphically add it to
  // `ciphertext0`.
  ASSERT_OK_AND_ASSIGN(RnsPolynomial<TypeParam> plaintext1,
                       encoder->EncodeBgv(messages1, this->main_moduli_));
  ASSERT_OK(plaintext1.ConvertToNttForm(this->main_moduli_));
  ASSERT_OK_AND_ASSIGN(RnsBgvCiphertext<TypeParam> ciphertext_sum_plaintext,
                       ciphertext0 + plaintext1);

  // This ciphertext should also decrypt to the sum of messages (mod t).
  ASSERT_OK_AND_ASSIGN(std::vector<Integer> dec_sum_plaintext_messages,
                       key.DecryptBgv(ciphertext_sum_plaintext, encoder));
  ASSERT_EQ(dec_sum_plaintext_messages.size(), num_coeffs);
  for (int i = 0; i < num_coeffs; ++i) {
    Integer expected = (messages0[i] + messages1[i]) % t;
    EXPECT_EQ(dec_sum_plaintext_messages[i], expected);
  }

  // Create a trivial ciphertext (plaintext, 0), and use `AddInPlaceWithoutPad`.
  RnsBgvCiphertext<TypeParam> ciphertext2({plaintext1}, this->main_moduli_,
                                          /*power_of_s=*/1, /*error=*/0,
                                          this->error_params_.get());
  ASSERT_OK(ciphertext0.AddInPlaceWithoutPad(ciphertext2));
  ASSERT_OK_AND_ASSIGN(std::vector<Integer> dec_sum_without_pad_messages,
                       key.DecryptBgv(ciphertext0, encoder));
  ASSERT_EQ(dec_sum_without_pad_messages.size(), num_coeffs);
  for (int i = 0; i < num_coeffs; ++i) {
    Integer expected = (messages0[i] + messages1[i]) % t;
    EXPECT_EQ(dec_sum_without_pad_messages[i], expected);
  }
}

TYPED_TEST(RnsBgvCiphertextTest, HomomorphicSubtraction) {
  using Integer = typename TypeParam::Int;
  ASSERT_OK_AND_ASSIGN(RnsRlweSecretKey<TypeParam> key, this->SampleKey());
  const CoefficientEncoder<TypeParam>* encoder = this->coeff_encoder_.get();

  // Sample two plaintext polynomials and encrypt them.
  int num_coeffs = 1 << this->rns_context_->LogN();
  Integer t = this->rns_context_->PlaintextModulus();
  std::vector<Integer> messages0 = testing::SampleMessages(num_coeffs, t);
  std::vector<Integer> messages1 = testing::SampleMessages(num_coeffs, t);
  ASSERT_OK_AND_ASSIGN(
      RnsBgvCiphertext<TypeParam> ciphertext0,
      key.EncryptBgv(messages0, this->coeff_encoder_.get(),
                     this->error_params_.get(), this->prng_.get()));
  ASSERT_OK_AND_ASSIGN(
      RnsBgvCiphertext<TypeParam> ciphertext1,
      key.EncryptBgv(messages1, this->coeff_encoder_.get(),
                     this->error_params_.get(), this->prng_.get()));

  // Homomorphically subtract the two ciphertexts.
  ASSERT_OK_AND_ASSIGN(RnsBgvCiphertext<TypeParam> ciphertext_diff,
                       ciphertext0 - ciphertext1);

  // The ciphertext should decrypt to the difference of messages (mod t).
  ASSERT_OK_AND_ASSIGN(std::vector<Integer> dec_diff_messages,
                       key.DecryptBgv(ciphertext_diff, encoder));
  ASSERT_EQ(dec_diff_messages.size(), num_coeffs);
  for (int i = 0; i < num_coeffs; ++i) {
    Integer expected = (t + messages0[i] - messages1[i]) % t;
    EXPECT_EQ(dec_diff_messages[i], expected);
  }

  // Encode `messages1` to a plaintext polynomial and homomorphically add it to
  // `ciphertext0`.
  ASSERT_OK_AND_ASSIGN(RnsPolynomial<TypeParam> plaintext1,
                       encoder->EncodeBgv(messages1, this->main_moduli_));
  ASSERT_OK(plaintext1.ConvertToNttForm(this->main_moduli_));
  ASSERT_OK_AND_ASSIGN(RnsBgvCiphertext<TypeParam> ciphertext_diff_plaintext,
                       ciphertext0 - plaintext1);

  // This ciphertext should also decrypt to the diff of messages (mod t).
  ASSERT_OK_AND_ASSIGN(std::vector<Integer> dec_diff_plaintext_messages,
                       key.DecryptBgv(ciphertext_diff_plaintext, encoder));
  ASSERT_EQ(dec_diff_plaintext_messages.size(), num_coeffs);
  for (int i = 0; i < num_coeffs; ++i) {
    Integer expected = (t + messages0[i] - messages1[i]) % t;
    EXPECT_EQ(dec_diff_plaintext_messages[i], expected);
  }

  // Create a trivial ciphertext (plaintext, 0), and use `SubInPlaceWithoutPad`.
  RnsBgvCiphertext<TypeParam> ciphertext2({plaintext1}, this->main_moduli_,
                                          /*power_of_s=*/1, /*error=*/0,
                                          this->error_params_.get());
  ASSERT_OK(ciphertext0.SubInPlaceWithoutPad(ciphertext2));
  ASSERT_OK_AND_ASSIGN(std::vector<Integer> dec_diff_without_pad_messages,
                       key.DecryptBgv(ciphertext0, encoder));
  ASSERT_EQ(dec_diff_without_pad_messages.size(), num_coeffs);
  for (int i = 0; i < num_coeffs; ++i) {
    Integer expected = (t + messages0[i] - messages1[i]) % t;
    EXPECT_EQ(dec_diff_without_pad_messages[i], expected);
  }
}

TYPED_TEST(RnsBgvCiphertextTest, ModReduceFailsIfAtLevelZero) {
  // Create a ciphertext with a single prime modulus.
  std::vector<const PrimeModulus<TypeParam>*> moduli_at_level_zero = {
      this->main_moduli_[0]};
  ASSERT_OK_AND_ASSIGN(auto zero,
                       RnsPolynomial<TypeParam>::CreateZero(
                           this->rns_context_->LogN(), moduli_at_level_zero));
  RnsBgvCiphertext<TypeParam> ciphertext(
      /*components=*/{zero, zero}, moduli_at_level_zero,
      /*power_of_s=*/1, /*error=*/0, this->error_params_.get());

  auto q_inv_mod_qs = this->rns_context_->MainPrimeModulusInverseResidues();
  RnsInt<TypeParam> ql_inv = q_inv_mod_qs[0].Prefix(0);
  EXPECT_THAT(
      ciphertext.ModReduce(this->rns_context_->PlaintextModulus(), ql_inv),
      StatusIs(absl::StatusCode::kFailedPrecondition,
               HasSubstr("Cannot perform ModReduce with insufficient "
                         "number of prime moduli")));
}

TYPED_TEST(RnsBgvCiphertextTest, ModReducedCiphertextDecrypts) {
  using Integer = typename TypeParam::Int;
  if (this->main_moduli_.size() < 2) {
    // There is only one prime modulus in the moduli chain, so we cannot perform
    // modulus reduction further.
    GTEST_SKIP() << "Insufficient number of prime moduli for ModReduce.";
  }

  ASSERT_OK_AND_ASSIGN(RnsRlweSecretKey<TypeParam> key, this->SampleKey());

  // Sample a plaintext polynomial and encrypt it wrt the full moduli chain.
  Integer t = this->rns_context_->PlaintextModulus();
  std::vector<Integer> messages =
      testing::SampleMessages(1 << this->rns_context_->LogN(), t);
  ASSERT_OK_AND_ASSIGN(
      RnsBgvCiphertext<TypeParam> ciphertext,
      key.EncryptBgv(messages, this->coeff_encoder_.get(),
                     this->error_params_.get(), this->prng_.get()));

  // ModReduce the ciphertext and get round(Q/q_L * ciphertext) mod (Q/q_L).
  auto q_inv_mod_qs = this->rns_context_->MainPrimeModulusInverseResidues();

  int level = ciphertext.Level();
  int degree = ciphertext.Degree();
  RnsInt<TypeParam> ql_inv = q_inv_mod_qs[level].Prefix(level);
  ASSERT_OK(ciphertext.ModReduce(t, ql_inv));
  ASSERT_EQ(ciphertext.Level(), level - 1);  // level should be reduced by 1
  ASSERT_EQ(ciphertext.Degree(), degree);    // degree should not change

  // We also need to modulus reduce the secret key.
  ASSERT_OK(key.ModReduce());
  ASSERT_EQ(key.Level(), level - 1);  // key level should match the ciphertext

  ASSERT_OK_AND_ASSIGN(std::vector<Integer> dec_messages,
                       key.DecryptBgv(ciphertext, this->coeff_encoder_.get()));
  EXPECT_EQ(dec_messages, messages);
}

////////////////////////////////////////////////////////////////////////////////
// Homomorphic operations with finite field packed encoding
////////////////////////////////////////////////////////////////////////////////

template <typename ModularInt>
class RnsBgvCiphertextPackedTest : public ::testing::Test {
 protected:
  void SetUp() override {
    // Create a PRNG for sampling secret key.
    auto prng = Prng::Create(kPrngSeed.substr(0, Prng::SeedLength()));
    CHECK(prng.ok());
    prng_ = std::move(prng.value());
  }

  void SetUpBgvRnsParameters(const testing::RnsParameters<ModularInt>& params) {
    // Use parameters suitable for finite field encoding.
    auto rns_context = RnsContext<ModularInt>::CreateForBgvFiniteFieldEncoding(
        params.log_n, params.qs, params.ps, params.t);
    CHECK(rns_context.ok());
    rns_context_ = std::make_unique<const RnsContext<ModularInt>>(
        std::move(rns_context.value()));
    main_moduli_ = rns_context_->MainPrimeModuli();

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

  void SetUpDefaultBgvRnsParameters() {
    auto all_params =
        testing::GetRnsParametersForFiniteFieldEncoding<ModularInt>();
    CHECK(!all_params.empty());
    SetUpBgvRnsParameters(all_params[0]);
  }

  // Sample a secret key from prng.
  absl::StatusOr<RnsRlweSecretKey<ModularInt>> SampleKey() const {
    return RnsRlweSecretKey<ModularInt>::Sample(rns_context_->LogN(), kVariance,
                                                rns_context_->MainPrimeModuli(),
                                                prng_.get());
  }

  std::unique_ptr<const RnsContext<ModularInt>> rns_context_;
  std::vector<const PrimeModulus<ModularInt>*> main_moduli_;
  std::unique_ptr<const FiniteFieldEncoder<ModularInt>> encoder_;
  std::unique_ptr<const RnsErrorParams<ModularInt>> error_params_;
  std::unique_ptr<Prng> prng_;
};

TYPED_TEST_SUITE(RnsBgvCiphertextPackedTest,
                 testing::ModularIntTypesForFiniteFieldEncoding);

TYPED_TEST(RnsBgvCiphertextPackedTest, NegatedCiphertext) {
  using Integer = typename TypeParam::Int;
  for (auto const& params :
       testing::GetRnsParametersForFiniteFieldEncoding<TypeParam>()) {
    this->SetUpBgvRnsParameters(params);
    ASSERT_OK_AND_ASSIGN(RnsRlweSecretKey<TypeParam> key, this->SampleKey());

    // Sample a plaintext polynomial and encrypt it.
    Integer t = this->rns_context_->PlaintextModulus();
    std::vector<Integer> messages =
        testing::SampleMessages(1 << this->rns_context_->LogN(), t);
    ASSERT_OK_AND_ASSIGN(
        RnsBgvCiphertext<TypeParam> ciphertext,
        key.EncryptBgv(messages, this->encoder_.get(),
                       this->error_params_.get(), this->prng_.get()));

    // Homomorphically negate the ciphertext.
    ASSERT_OK_AND_ASSIGN(RnsBgvCiphertext<TypeParam> ciphertext_neg,
                         ciphertext.Negate());

    // The negated ciphertext should decrypt to negated plaintext.
    ASSERT_OK_AND_ASSIGN(std::vector<Integer> dec_neg_messages,
                         key.DecryptBgv(ciphertext_neg, this->encoder_.get()));
    ASSERT_EQ(dec_neg_messages.size(), messages.size());
    for (int i = 0; i < messages.size(); ++i) {
      Integer expected = (t - messages[i]) % t;  // -messages[i] (mod t).
      EXPECT_EQ(dec_neg_messages[i], expected);
    }
  }
}

TYPED_TEST(RnsBgvCiphertextPackedTest, HomomorphicAddition) {
  using Integer = typename TypeParam::Int;
  for (auto const& params :
       testing::GetRnsParametersForFiniteFieldEncoding<TypeParam>()) {
    this->SetUpBgvRnsParameters(params);

    ASSERT_OK_AND_ASSIGN(RnsRlweSecretKey<TypeParam> key, this->SampleKey());
    const FiniteFieldEncoder<TypeParam>* encoder = this->encoder_.get();

    // Sample two plaintext polynomials and encrypt them.
    int num_coeffs = 1 << this->rns_context_->LogN();
    Integer t = this->rns_context_->PlaintextModulus();
    std::vector<Integer> messages0 = testing::SampleMessages(num_coeffs, t);
    std::vector<Integer> messages1 = testing::SampleMessages(num_coeffs, t);
    ASSERT_OK_AND_ASSIGN(
        RnsBgvCiphertext<TypeParam> ciphertext0,
        key.EncryptBgv(messages0, encoder, this->error_params_.get(),
                       this->prng_.get()));
    ASSERT_OK_AND_ASSIGN(
        RnsBgvCiphertext<TypeParam> ciphertext1,
        key.EncryptBgv(messages1, encoder, this->error_params_.get(),
                       this->prng_.get()));

    // Homomorphically add the two ciphertexts.
    ASSERT_OK_AND_ASSIGN(RnsBgvCiphertext<TypeParam> ciphertext_sum,
                         ciphertext0 + ciphertext1);

    // The ciphertext should decrypt to the sum of messages (mod t).
    ASSERT_OK_AND_ASSIGN(std::vector<Integer> dec_sum_messages,
                         key.DecryptBgv(ciphertext_sum, encoder));
    ASSERT_EQ(dec_sum_messages.size(), num_coeffs);
    for (int i = 0; i < num_coeffs; ++i) {
      Integer expected = (messages0[i] + messages1[i]) % t;
      EXPECT_EQ(dec_sum_messages[i], expected);
    }

    // Encode `messages1` to a plaintext polynomial and homomorphically add it
    // to `ciphertext0`.
    ASSERT_OK_AND_ASSIGN(RnsPolynomial<TypeParam> plaintext1,
                         encoder->EncodeBgv(messages1, this->main_moduli_));
    ASSERT_OK_AND_ASSIGN(RnsBgvCiphertext<TypeParam> ciphertext_sum_plaintext,
                         ciphertext0 + plaintext1);

    // This ciphertext should also decrypt to the sum of messages (mod t).
    ASSERT_OK_AND_ASSIGN(std::vector<Integer> dec_sum_plaintext_messages,
                         key.DecryptBgv(ciphertext_sum_plaintext, encoder));
    ASSERT_EQ(dec_sum_plaintext_messages.size(), num_coeffs);
    for (int i = 0; i < num_coeffs; ++i) {
      Integer expected = (messages0[i] + messages1[i]) % t;
      EXPECT_EQ(dec_sum_plaintext_messages[i], expected);
    }
  }
}

TYPED_TEST(RnsBgvCiphertextPackedTest, HomomorphicSubtraction) {
  using Integer = typename TypeParam::Int;
  for (auto const& params :
       testing::GetRnsParametersForFiniteFieldEncoding<TypeParam>()) {
    this->SetUpBgvRnsParameters(params);
    ASSERT_OK_AND_ASSIGN(RnsRlweSecretKey<TypeParam> key, this->SampleKey());
    const FiniteFieldEncoder<TypeParam>* encoder = this->encoder_.get();

    // Sample two plaintext polynomials and encrypt them.
    int num_coeffs = 1 << this->rns_context_->LogN();
    Integer t = this->rns_context_->PlaintextModulus();
    std::vector<Integer> messages0 = testing::SampleMessages(num_coeffs, t);
    std::vector<Integer> messages1 = testing::SampleMessages(num_coeffs, t);
    ASSERT_OK_AND_ASSIGN(
        RnsBgvCiphertext<TypeParam> ciphertext0,
        key.EncryptBgv(messages0, encoder, this->error_params_.get(),
                       this->prng_.get()));
    ASSERT_OK_AND_ASSIGN(
        RnsBgvCiphertext<TypeParam> ciphertext1,
        key.EncryptBgv(messages1, encoder, this->error_params_.get(),
                       this->prng_.get()));

    // Homomorphically subtract the two ciphertexts.
    ASSERT_OK_AND_ASSIGN(RnsBgvCiphertext<TypeParam> ciphertext_diff,
                         ciphertext0 - ciphertext1);

    // The ciphertext should decrypt to the difference of messages (mod t).
    ASSERT_OK_AND_ASSIGN(std::vector<Integer> dec_diff_messages,
                         key.DecryptBgv(ciphertext_diff, encoder));
    ASSERT_EQ(dec_diff_messages.size(), num_coeffs);
    for (int i = 0; i < num_coeffs; ++i) {
      Integer expected = (t + messages0[i] - messages1[i]) % t;
      EXPECT_EQ(dec_diff_messages[i], expected);
    }

    // Encode `messages1` to a plaintext polynomial and homomorphically add it
    // to `ciphertext0`.
    ASSERT_OK_AND_ASSIGN(RnsPolynomial<TypeParam> plaintext1,
                         encoder->EncodeBgv(messages1, this->main_moduli_));
    ASSERT_OK_AND_ASSIGN(RnsBgvCiphertext<TypeParam> ciphertext_diff_plaintext,
                         ciphertext0 - plaintext1);

    // This ciphertext should also decrypt to the diff of messages (mod t).
    ASSERT_OK_AND_ASSIGN(std::vector<Integer> dec_diff_plaintext_messages,
                         key.DecryptBgv(ciphertext_diff_plaintext, encoder));
    ASSERT_EQ(dec_diff_plaintext_messages.size(), num_coeffs);
    for (int i = 0; i < num_coeffs; ++i) {
      Integer expected = (t + messages0[i] - messages1[i]) % t;
      EXPECT_EQ(dec_diff_plaintext_messages[i], expected);
    }
  }
}

TYPED_TEST(RnsBgvCiphertextPackedTest, HomomorphicMulWithPlaintext) {
  using Integer = typename TypeParam::Int;
  for (auto const& params :
       testing::GetRnsParametersForFiniteFieldEncoding<TypeParam>()) {
    this->SetUpBgvRnsParameters(params);
    int log_n = this->rns_context_->LogN();
    ASSERT_OK_AND_ASSIGN(RnsRlweSecretKey<TypeParam> key, this->SampleKey());

    Integer t = this->rns_context_->PlaintextModulus();
    std::vector<Integer> messages = testing::SampleMessages(1 << log_n, t);
    ASSERT_OK_AND_ASSIGN(RnsBgvCiphertext<TypeParam> ciphertext,
                         key.template EncryptBgv<FiniteFieldEncoder<TypeParam>>(
                             messages, this->encoder_.get(),
                             this->error_params_.get(), this->prng_.get()));

    // Multiply the ciphertext with a plaintext polynomial that encodes a
    // projection vector.
    std::vector<Integer> projections(1 << log_n, 0);
    for (int i = 0; i < (1 << log_n); i += 2) {
      projections[i] = 1;
    }
    ASSERT_OK_AND_ASSIGN(
        RnsPolynomial<TypeParam> plaintext,
        this->encoder_->EncodeBgv(projections, this->main_moduli_));

    ASSERT_OK(ciphertext.AbsorbInPlace(plaintext));

    ASSERT_OK_AND_ASSIGN(std::vector<Integer> dec_messages,
                         key.template DecryptBgv<FiniteFieldEncoder<TypeParam>>(
                             ciphertext, this->encoder_.get()));
    for (int i = 0; i < (1 << log_n); ++i) {
      EXPECT_EQ(dec_messages[i], messages[i] * projections[i]);
    }
  }
}

TYPED_TEST(RnsBgvCiphertextPackedTest, HomomorphicFusedAbsorbAdd) {
  using Integer = typename TypeParam::Int;
  for (auto const& params :
       testing::GetRnsParametersForFiniteFieldEncoding<TypeParam>()) {
    this->SetUpBgvRnsParameters(params);
    int log_n = this->rns_context_->LogN();
    ASSERT_OK_AND_ASSIGN(RnsRlweSecretKey<TypeParam> key, this->SampleKey());

    Integer t = this->rns_context_->PlaintextModulus();
    std::vector<Integer> messages0 = testing::SampleMessages(1 << log_n, t);
    std::vector<Integer> messages1 = testing::SampleMessages(1 << log_n, t);
    ASSERT_OK_AND_ASSIGN(RnsBgvCiphertext<TypeParam> ciphertext0,
                         key.template EncryptBgv<FiniteFieldEncoder<TypeParam>>(
                             messages0, this->encoder_.get(),
                             this->error_params_.get(), this->prng_.get()));
    ASSERT_OK_AND_ASSIGN(RnsBgvCiphertext<TypeParam> ciphertext1,
                         key.template EncryptBgv<FiniteFieldEncoder<TypeParam>>(
                             messages1, this->encoder_.get(),
                             this->error_params_.get(), this->prng_.get()));

    std::vector<Integer> projections(1 << log_n, 0);
    for (int i = 0; i < (1 << log_n); i += 2) {
      projections[i] = 1;
    }
    ASSERT_OK_AND_ASSIGN(
        RnsPolynomial<TypeParam> plaintext,
        this->encoder_->EncodeBgv(projections, this->main_moduli_));

    ASSERT_OK(ciphertext0.FusedAbsorbAddInPlace(ciphertext1, plaintext));

    ASSERT_OK_AND_ASSIGN(std::vector<Integer> dec_messages,
                         key.template DecryptBgv<FiniteFieldEncoder<TypeParam>>(
                             ciphertext0, this->encoder_.get()));
    for (int i = 0; i < (1 << log_n); ++i) {
      Integer expected = (messages0[i] + messages1[i] * projections[i]) % t;
      EXPECT_EQ(dec_messages[i], expected);
    }
  }
}

TYPED_TEST(RnsBgvCiphertextPackedTest,
           HomomorphicFusedAbsorbAddWithComputedPad) {
  using Integer = typename TypeParam::Int;
  for (auto const& params :
       testing::GetRnsParametersForFiniteFieldEncoding<TypeParam>()) {
    this->SetUpBgvRnsParameters(params);
    int log_n = this->rns_context_->LogN();
    ASSERT_OK_AND_ASSIGN(RnsRlweSecretKey<TypeParam> key, this->SampleKey());

    Integer t = this->rns_context_->PlaintextModulus();
    std::vector<Integer> messages0 = testing::SampleMessages(1 << log_n, t);
    std::vector<Integer> messages1 = testing::SampleMessages(1 << log_n, t);
    ASSERT_OK_AND_ASSIGN(RnsBgvCiphertext<TypeParam> ciphertext0,
                         key.template EncryptBgv<FiniteFieldEncoder<TypeParam>>(
                             messages0, this->encoder_.get(),
                             this->error_params_.get(), this->prng_.get()));
    ASSERT_OK_AND_ASSIGN(RnsBgvCiphertext<TypeParam> ciphertext1,
                         key.template EncryptBgv<FiniteFieldEncoder<TypeParam>>(
                             messages1, this->encoder_.get(),
                             this->error_params_.get(), this->prng_.get()));

    std::vector<Integer> projections(1 << log_n, 0);
    for (int i = 0; i < (1 << log_n); i += 2) {
      projections[i] = 1;
    }
    ASSERT_OK_AND_ASSIGN(
        RnsPolynomial<TypeParam> plaintext,
        this->encoder_->EncodeBgv(projections, this->main_moduli_));

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
                         key.template DecryptBgv<FiniteFieldEncoder<TypeParam>>(
                             ciphertext0, this->encoder_.get()));
    for (int i = 0; i < (1 << log_n); ++i) {
      Integer expected = (messages0[i] + messages1[i] * projections[i]) % t;
      EXPECT_EQ(dec_messages[i], expected);
    }
  }
}

TYPED_TEST(RnsBgvCiphertextPackedTest, HomomorphicMulWithCiphertext) {
  using Integer = typename TypeParam::Int;
  for (auto const& params :
       testing::GetRnsParametersForFiniteFieldEncoding<TypeParam>()) {
    this->SetUpBgvRnsParameters(params);

    int log_n = this->rns_context_->LogN();
    double expected_error_bound =
        this->error_params_->B_secretkey_encryption() *
        this->error_params_->B_secretkey_encryption() * sqrt(1 << log_n);
    double composite_modulus = 1;
    for (auto prime_modulus : this->main_moduli_) {
      composite_modulus *= static_cast<double>(prime_modulus->Modulus());
    }
    if (expected_error_bound >= composite_modulus) {
      // The test parameters are not suitable for homomorphic multiplication as
      // the error in the product ciphertext is expected to be larger than the
      // ciphertext modulus.
      continue;
    }

    ASSERT_OK_AND_ASSIGN(RnsRlweSecretKey<TypeParam> key, this->SampleKey());

    Integer t = this->rns_context_->PlaintextModulus();
    std::vector<Integer> messages0 = testing::SampleMessages(1 << log_n, t);
    std::vector<Integer> messages1 = testing::SampleMessages(1 << log_n, t);
    ASSERT_OK_AND_ASSIGN(RnsBgvCiphertext<TypeParam> ciphertext0,
                         key.template EncryptBgv<FiniteFieldEncoder<TypeParam>>(
                             messages0, this->encoder_.get(),
                             this->error_params_.get(), this->prng_.get()));
    ASSERT_OK_AND_ASSIGN(RnsBgvCiphertext<TypeParam> ciphertext1,
                         key.template EncryptBgv<FiniteFieldEncoder<TypeParam>>(
                             messages1, this->encoder_.get(),
                             this->error_params_.get(), this->prng_.get()));

    // Multiply the two ciphertexts, and the result is a degree 2 ciphertext.
    ASSERT_OK_AND_ASSIGN(RnsBgvCiphertext<TypeParam> ciphertext_product,
                         ciphertext0 * ciphertext1);
    EXPECT_EQ(ciphertext_product.Degree(), 2);
    EXPECT_EQ(ciphertext_product.Level(), ciphertext0.Level());

    // Decrypt the product ciphertext, and expect the result is the pointwise
    // product of messages0 and messages1.
    ASSERT_OK_AND_ASSIGN(std::vector<Integer> dec_messages,
                         key.template DecryptBgv<FiniteFieldEncoder<TypeParam>>(
                             ciphertext_product, this->encoder_.get()));
    for (int i = 0; i < (1 << log_n); ++i) {
      // Both messages0 and messages1 contains numbers mod t, and since the
      // plaintext modulus t is chosen to be sufficiently small, we can multiply
      // two messages without overflow the int type.
      EXPECT_EQ(dec_messages[i], (messages0[i] * messages1[i]) % t);
    }
  }
}

}  // namespace
}  // namespace rlwe
