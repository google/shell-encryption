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

#include "shell_encryption/rns/rns_bgv_public_key.h"

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
#include "shell_encryption/rns/rns_bgv_ciphertext.h"
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
constexpr PrngType kPrngType = PRNG_TYPE_HKDF;

template <typename ModularInt>
class RnsBgvPublicKeyTest : public ::testing::Test {
 protected:
  void SetUp() override {
    testing::RnsParameters<ModularInt> rns_params =
        testing::GetRnsParameters<ModularInt>();

    // Do not initialize the test environment if the parameters are too small
    // for public key encryption. The test cases should check if rns_context_
    // has been initialized before running any test.
    if (!IsParameterTooSmallForPublicKeyEncryption(rns_params)) {
      GTEST_SKIP() << "Test parameters are too small for public key tests";
    }

    // Create the RNS context.
    auto rns_context = RnsContext<ModularInt>::Create(
        rns_params.log_n, rns_params.qs, rns_params.ps, rns_params.t);
    CHECK(rns_context.ok());
    rns_context_ = std::make_unique<const RnsContext<ModularInt>>(
        std::move(rns_context.value()));
    moduli_ = rns_context_->MainPrimeModuli();

    // Create a coefficient encoder.
    auto encoder = CoefficientEncoder<ModularInt>::Create(rns_context_.get());
    CHECK(encoder.ok());
    encoder_ = std::make_unique<const CoefficientEncoder<ModularInt>>(
        std::move(encoder.value()));

    // Create the error parameters.
    int log_t = floor(std::log2(static_cast<double>(rns_params.t)));
    auto error_params = RnsErrorParams<ModularInt>::Create(
        rns_params.log_n, moduli_, /*aux_moduli=*/{}, log_t, sqrt(kVariance));
    CHECK(error_params.ok());
    error_params_ = std::make_unique<const RnsErrorParams<ModularInt>>(
        std::move(error_params.value()));

    // Create a PRNG for sampling secret key.
    auto prng = Prng::Create(kPrngSeed.substr(0, Prng::SeedLength()));
    CHECK(prng.ok());
    prng_ = std::move(prng.value());

    auto secret_key = RnsRlweSecretKey<ModularInt>::Sample(
        rns_context_->LogN(), kVariance, rns_context_->MainPrimeModuli(),
        prng_.get());
    CHECK(secret_key.ok());
    secret_key_ = std::make_unique<const RnsRlweSecretKey<ModularInt>>(
        std::move(secret_key.value()));
  }

  bool IsParameterTooSmallForPublicKeyEncryption(
      const testing::RnsParameters<ModularInt>& params) const {
    int main_modulus_bits = 0;
    for (auto q : params.qs) {
      main_modulus_bits += std::floor(std::log2(static_cast<double>(q)));
    }
    double expected_error =
        16.0 * kVariance * (1 << params.log_n) * static_cast<double>(params.t);
    double expected_error_bits = std::ceil(std::log2(expected_error));
    return expected_error_bits + 1 < main_modulus_bits;
  }

  std::unique_ptr<const RnsContext<ModularInt>> rns_context_;
  std::vector<const PrimeModulus<ModularInt>*> moduli_;
  std::unique_ptr<const CoefficientEncoder<ModularInt>> encoder_;
  std::unique_ptr<const RnsErrorParams<ModularInt>> error_params_;
  std::unique_ptr<Prng> prng_;
  std::unique_ptr<const RnsRlweSecretKey<ModularInt>> secret_key_;
};

TYPED_TEST_SUITE(RnsBgvPublicKeyTest, testing::ModularIntTypes);

template <typename ModularInt>
class RnsBgvPublicKeyNegativeTest : public RnsBgvPublicKeyTest<ModularInt> {};

TYPED_TEST_SUITE(RnsBgvPublicKeyNegativeTest,
                 testing::ModularIntTypesForNegativeTests);

TYPED_TEST(RnsBgvPublicKeyNegativeTest, CreateFailsIfVarianceIsNotPositive) {
  EXPECT_THAT(RnsBgvPublicKey<TypeParam>::Create(
                  *this->secret_key_, /*variance=*/0, kPrngType,
                  this->rns_context_->PlaintextModulus()),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`variance` must be positive")));
  EXPECT_THAT(RnsBgvPublicKey<TypeParam>::Create(
                  *this->secret_key_, /*variance=*/-1, kPrngType,
                  this->rns_context_->PlaintextModulus()),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`variance` must be positive")));
}

TYPED_TEST(RnsBgvPublicKeyNegativeTest, CreateFailsIfPrngTypeIsInvalid) {
  EXPECT_THAT(RnsBgvPublicKey<TypeParam>::Create(
                  *this->secret_key_, kVariance, PRNG_TYPE_INVALID,
                  this->rns_context_->PlaintextModulus()),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("PrngType not specified correctly")));
}

TYPED_TEST(RnsBgvPublicKeyNegativeTest, EncryptFailsWithNullEncoder) {
  using Encoder = CoefficientEncoder<TypeParam>;

  ASSERT_OK_AND_ASSIGN(RnsBgvPublicKey<TypeParam> public_key,
                       RnsBgvPublicKey<TypeParam>::Create(
                           *this->secret_key_, kVariance, kPrngType,
                           this->rns_context_->PlaintextModulus()));
  EXPECT_THAT(public_key.template Encrypt<Encoder>(
                  /*messages=*/{}, /*encoder=*/nullptr,
                  this->error_params_.get(), this->prng_.get()),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`encoder` must not be null")));
}

TYPED_TEST(RnsBgvPublicKeyNegativeTest, EncryptFailsWithNullErrorParams) {
  ASSERT_OK_AND_ASSIGN(RnsBgvPublicKey<TypeParam> public_key,
                       RnsBgvPublicKey<TypeParam>::Create(
                           *this->secret_key_, kVariance, kPrngType,
                           this->rns_context_->PlaintextModulus()));
  EXPECT_THAT(public_key.Encrypt(
                  /*messages=*/{}, this->encoder_.get(),
                  /*error_params=*/nullptr, this->prng_.get()),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`error_params` must not be null")));
}

TYPED_TEST(RnsBgvPublicKeyNegativeTest, EncryptFailsWithNullPrng) {
  ASSERT_OK_AND_ASSIGN(RnsBgvPublicKey<TypeParam> public_key,
                       RnsBgvPublicKey<TypeParam>::Create(
                           *this->secret_key_, kVariance, kPrngType,
                           this->rns_context_->PlaintextModulus()));
  EXPECT_THAT(public_key.Encrypt(
                  /*messages=*/{}, this->encoder_.get(),
                  this->error_params_.get(), /*prng=*/nullptr),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`prng` must not be null")));
}

////////////////////////////////////////////////////////////////////////////////
// Functional tests
////////////////////////////////////////////////////////////////////////////////

TYPED_TEST(RnsBgvPublicKeyTest, CreateSucceeds) {
  ASSERT_OK_AND_ASSIGN(RnsBgvPublicKey<TypeParam> public_key,
                       RnsBgvPublicKey<TypeParam>::Create(
                           *this->secret_key_, kVariance, kPrngType,
                           this->rns_context_->PlaintextModulus()));
  const int log_n = this->rns_context_->LogN();
  const int num_coeffs = 1 << log_n;
  const int num_moduli = this->moduli_.size();
  EXPECT_EQ(public_key.LogN(), log_n);
  EXPECT_EQ(public_key.Level(), num_moduli - 1);
  EXPECT_EQ(public_key.NumCoeffs(), num_coeffs);
  EXPECT_EQ(public_key.NumModuli(), num_moduli);
  const RnsPolynomial<TypeParam>& key_a = public_key.KeyA();
  const RnsPolynomial<TypeParam>& key_b = public_key.KeyB();
  EXPECT_EQ(key_a.LogN(), log_n);
  EXPECT_EQ(key_a.NumCoeffs(), num_coeffs);
  EXPECT_EQ(key_a.NumModuli(), num_moduli);
  EXPECT_EQ(key_b.LogN(), log_n);
  EXPECT_EQ(key_b.NumCoeffs(), num_coeffs);
  EXPECT_EQ(key_b.NumModuli(), num_moduli);
}

TYPED_TEST(RnsBgvPublicKeyTest, EncryptDecrypts) {
  using Integer = typename TypeParam::Int;

  // Generate a public key.
  Integer t = this->rns_context_->PlaintextModulus();
  ASSERT_OK_AND_ASSIGN(RnsBgvPublicKey<TypeParam> public_key,
                       RnsBgvPublicKey<TypeParam>::Create(
                           *this->secret_key_, kVariance, kPrngType, t));

  int num_messages = 1 << this->rns_context_->LogN();
  std::vector<Integer> messages = testing::SampleMessages(num_messages, t);
  ASSERT_OK_AND_ASSIGN(
      RnsBgvCiphertext<TypeParam> ciphertext,
      public_key.Encrypt(messages, this->encoder_.get(),
                         this->error_params_.get(), this->prng_.get()));

  ASSERT_OK_AND_ASSIGN(
      std::vector<Integer> dec_messages,
      this->secret_key_->DecryptBgv(ciphertext, this->encoder_.get()));
  EXPECT_EQ(dec_messages, messages);
}

}  // namespace
}  // namespace rlwe
