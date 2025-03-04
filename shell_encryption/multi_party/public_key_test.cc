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

#include "shell_encryption/multi_party/public_key.h"

#include <cmath>
#include <memory>
#include <utility>
#include <vector>

#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include "absl/log/check.h"
#include "absl/status/status.h"
#include "shell_encryption/multi_party/public_key_share.h"
#include "shell_encryption/multi_party/public_parameter.h"
#include "shell_encryption/multi_party/secret_key_share.h"
#include "shell_encryption/multi_party/testing/parameters.h"
#include "shell_encryption/rns/coefficient_encoder.h"
#include "shell_encryption/rns/rns_bfv_ciphertext.h"
#include "shell_encryption/rns/rns_context.h"
#include "shell_encryption/rns/rns_error_params.h"
#include "shell_encryption/rns/rns_modulus.h"
#include "shell_encryption/rns/rns_polynomial.h"
#include "shell_encryption/rns/testing/testing_utils.h"
#include "shell_encryption/testing/status_matchers.h"
#include "shell_encryption/testing/status_testing.h"
#include "shell_encryption/testing/testing_prng.h"

namespace rlwe {
namespace multi_party {
namespace {

using Prng = ::rlwe::testing::TestingPrng;
using ::rlwe::testing::SampleMessages;
using ::rlwe::testing::StatusIs;
using ::testing::HasSubstr;
using testing::ModularInt32;

constexpr int kVariance = 8;
constexpr PrngType kPrngType = PRNG_TYPE_HKDF;

template <typename ModularInt>
class PublicKeyTest : public ::testing::Test {
 protected:
  void SetUp() override {
    testing::MpaheParameters<ModularInt> params =
        testing::GetMultiPartyParameters<ModularInt>();

    // Create the RNS context.
    auto rns_context = RnsContext<ModularInt>::CreateForBfv(
        params.log_n, params.qs, params.ps, params.t);
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
    int log_t = std::floor(std::log2(static_cast<double>(params.t)));
    auto error_params = RnsErrorParams<ModularInt>::Create(
        params.log_n, moduli_, /*aux_moduli=*/{}, log_t, sqrt(kVariance));
    CHECK(error_params.ok());
    error_params_ = std::make_unique<const RnsErrorParams<ModularInt>>(
        std::move(error_params.value()));

    // Create the public parameter, which is always needed in tests.
    auto public_parameter = PublicParameter<ModularInt>::Create(
        rns_context_.get(), kVariance, kPrngType);
    CHECK(public_parameter.ok());
    public_parameter_ = std::move(public_parameter.value());

    // Create a PRNG for testing only.
    prng_ = std::make_unique<Prng>(0);
  }

  void SetupSecretKeyShare() {
    auto secret_key_share =
        SecretKeyShare<ModularInt>::Sample(rns_context_.get(), prng_.get());
    CHECK(secret_key_share.ok());
    secret_key_share_ = std::make_unique<const SecretKeyShare<ModularInt>>(
        std::move(secret_key_share.value()));
  }

  void SetupPublicKeyShare() {
    auto public_key_share = PublicKeyShare<ModularInt>::Create(
        secret_key_share_.get(), public_parameter_.get(), kPrngType);
    CHECK(public_key_share.ok());
    public_key_share_ = std::make_unique<const PublicKeyShare<ModularInt>>(
        std::move(public_key_share.value()));
  }

  RnsPolynomial<ModularInt> RawDecrypt(
      const RnsPolynomial<ModularInt>& secret_key,
      const RnsBfvCiphertext<ModularInt>& ciphertext) {
    CHECK(ciphertext.Degree() == 1) << "only degree 1 ciphertext is supported";
    RnsPolynomial<ModularInt> c0 = ciphertext.Component(0).value();
    RnsPolynomial<ModularInt> c1 = ciphertext.Component(1).value();
    auto status = c0.FusedMulAddInPlace(c1, secret_key, moduli_);
    CHECK(status.ok());
    return c0;
  }

  std::unique_ptr<const RnsContext<ModularInt>> rns_context_;
  std::vector<const PrimeModulus<ModularInt>*> moduli_;
  std::unique_ptr<const CoefficientEncoder<ModularInt>> encoder_;
  std::unique_ptr<const RnsErrorParams<ModularInt>> error_params_;
  std::unique_ptr<const PublicParameter<ModularInt>> public_parameter_;
  std::unique_ptr<const SecretKeyShare<ModularInt>> secret_key_share_;
  std::unique_ptr<const PublicKeyShare<ModularInt>> public_key_share_;
  std::unique_ptr<Prng> prng_;
};

TYPED_TEST_SUITE(PublicKeyTest, testing::ModularIntTypesForMultiParty);

class PublicKeyNegativeTest : public PublicKeyTest<ModularInt32> {};

TEST_F(PublicKeyNegativeTest, CreateFailsIfPublicParameterIsNull) {
  SetupSecretKeyShare();
  SetupPublicKeyShare();
  EXPECT_THAT(PublicKey<ModularInt32>::Create(
                  /*public_parameter=*/nullptr, {*public_key_share_.get()}),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`public_parameter` must not be null")));
}

TEST_F(PublicKeyNegativeTest, CreateFailsIfEmptyPublicKeyShares) {
  EXPECT_THAT(PublicKey<ModularInt32>::Create(public_parameter_.get(),
                                              /*public_key_shares=*/{}),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`public_key_shares` must not be empty")));
}

TEST_F(PublicKeyNegativeTest, DeserializeFailsIfPublicParameterIsNull) {
  SerializedPublicKey serialized;
  EXPECT_THAT(PublicKey<ModularInt32>::Deserialize(
                  serialized, /*public_parameter=*/nullptr),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`public_parameter` must not be null")));
}

TYPED_TEST(PublicKeyTest, Encrypt) {
  using Integer = typename TypeParam::Int;

  this->SetupSecretKeyShare();
  this->SetupPublicKeyShare();
  ASSERT_OK_AND_ASSIGN(
      PublicKey<TypeParam> public_key,
      PublicKey<TypeParam>::Create(this->public_parameter_.get(),
                                   {*this->public_key_share_.get()}));

  int num_coeffs = 1 << this->rns_context_->LogN();
  std::vector<Integer> messages = SampleMessages<Integer>(num_coeffs, 8);

  ASSERT_OK_AND_ASSIGN(
      auto ciphertext,
      public_key.Encrypt(messages, this->encoder_.get(),
                         this->error_params_.get(), this->prng_.get()));
  EXPECT_EQ(ciphertext.Degree(), 1);
  EXPECT_EQ(ciphertext.Level(), this->moduli_.size() - 1);

  // Decrypt using the only secret key share.
  RnsPolynomial<TypeParam> noisy_plaintext =
      this->RawDecrypt(this->secret_key_share_->Key(), ciphertext);
  ASSERT_OK_AND_ASSIGN(
      std::vector<Integer> decrypted,
      this->encoder_->DecodeBfv(noisy_plaintext, this->moduli_));
  EXPECT_EQ(decrypted, messages);
}

TYPED_TEST(PublicKeyTest, EncryptExplicit) {
  using Integer = typename TypeParam::Int;

  this->SetupSecretKeyShare();
  this->SetupPublicKeyShare();
  ASSERT_OK_AND_ASSIGN(
      PublicKey<TypeParam> public_key,
      PublicKey<TypeParam>::Create(this->public_parameter_.get(),
                                   {*this->public_key_share_.get()}));

  int num_coeffs = 1 << this->rns_context_->LogN();
  std::vector<Integer> messages = SampleMessages<Integer>(num_coeffs, 8);

  auto b = RnsPolynomial<TypeParam>::CreateEmpty();
  auto a = RnsPolynomial<TypeParam>::CreateEmpty();
  auto r = RnsPolynomial<TypeParam>::CreateEmpty();
  auto e = RnsPolynomial<TypeParam>::CreateEmpty();

  ASSERT_OK(public_key.EncryptExplicit(messages, this->encoder_.get(),
                                       this->error_params_.get(),
                                       this->prng_.get(), &b, &a, &r, &e));

  // Check that everything gets populated.
  EXPECT_EQ(b.NumCoeffs(), 1 << this->rns_context_->LogN());
  EXPECT_EQ(b.NumModuli(), this->rns_context_->NumMainPrimeModuli());
  EXPECT_TRUE(b.IsNttForm());
  EXPECT_EQ(a.NumCoeffs(), 1 << this->rns_context_->LogN());
  EXPECT_EQ(a.NumModuli(), this->rns_context_->NumMainPrimeModuli());
  EXPECT_TRUE(a.IsNttForm());
  EXPECT_EQ(r.NumCoeffs(), 1 << this->rns_context_->LogN());
  EXPECT_EQ(r.NumModuli(), this->rns_context_->NumMainPrimeModuli());
  EXPECT_TRUE(r.IsNttForm());
  EXPECT_EQ(e.NumCoeffs(), 1 << this->rns_context_->LogN());
  EXPECT_EQ(e.NumModuli(), this->rns_context_->NumMainPrimeModuli());
  EXPECT_TRUE(e.IsNttForm());

  // Decrypt using the only secret key share.
  ASSERT_OK(
      b.FusedMulAddInPlace(a, this->secret_key_share_->Key(), this->moduli_));
  ASSERT_OK_AND_ASSIGN(std::vector<Integer> decrypted,
                       this->encoder_->DecodeBfv(b, this->moduli_));
  EXPECT_EQ(decrypted, messages);
}

TYPED_TEST(PublicKeyTest, SerializeDeserialize) {
  this->SetupSecretKeyShare();
  this->SetupPublicKeyShare();
  ASSERT_OK_AND_ASSIGN(
      PublicKey<TypeParam> public_key,
      PublicKey<TypeParam>::Create(this->public_parameter_.get(),
                                   {*this->public_key_share_.get()}));

  ASSERT_OK_AND_ASSIGN(auto serialized, public_key.Serialize());
  ASSERT_OK_AND_ASSIGN(auto deserialized,
                       PublicKey<TypeParam>::Deserialize(
                           serialized, this->public_parameter_.get()));
  EXPECT_EQ(deserialized.ComponentB(), public_key.ComponentB());
}

}  // namespace
}  // namespace multi_party
}  // namespace rlwe
