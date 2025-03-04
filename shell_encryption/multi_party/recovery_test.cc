/*
 * Copyright 2025 Google LLC.
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

#include "shell_encryption/multi_party/recovery.h"

#include <cmath>
#include <memory>
#include <utility>
#include <vector>

#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include "absl/log/check.h"
#include "absl/status/status.h"
#include "shell_encryption/multi_party/public_key.h"
#include "shell_encryption/multi_party/public_key_share.h"
#include "shell_encryption/multi_party/public_parameter.h"
#include "shell_encryption/multi_party/secret_key_share.h"
#include "shell_encryption/multi_party/testing/parameters.h"
#include "shell_encryption/rns/coefficient_encoder.h"
#include "shell_encryption/rns/rns_context.h"
#include "shell_encryption/rns/rns_error_params.h"
#include "shell_encryption/rns/rns_modulus.h"
#include "shell_encryption/rns/rns_polynomial.h"
#include "shell_encryption/rns/testing/testing_utils.h"
#include "shell_encryption/sampler/discrete_gaussian.h"
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

constexpr int kNumParties = 3;
constexpr int kPublicKeyVariance = 8;
constexpr int kMaxValue = 72;
constexpr double kSBase = 12.8;
constexpr PrngType kPrngType = PRNG_TYPE_HKDF;

template <typename ModularInt>
class RecoverTest : public ::testing::Test {
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
        params.log_n, moduli_, /*aux_moduli=*/{}, log_t,
        sqrt(kPublicKeyVariance));
    CHECK(error_params.ok());
    error_params_ = std::make_unique<const RnsErrorParams<ModularInt>>(
        std::move(error_params.value()));

    // Create public parameter.
    auto public_parameter = PublicParameter<ModularInt>::Create(
        rns_context_.get(), kPublicKeyVariance, kPrngType);
    CHECK(public_parameter.ok());
    public_parameter_ = std::move(public_parameter.value());

    params_ = std::make_unique<const testing::MpaheParameters<ModularInt>>(
        std::move(params));

    // Create a PRNG for testing only.
    prng_ = std::make_unique<Prng>(0);
  }

  std::unique_ptr<const testing::MpaheParameters<ModularInt>> params_;
  std::unique_ptr<const RnsContext<ModularInt>> rns_context_;
  std::vector<const PrimeModulus<ModularInt>*> moduli_;
  std::unique_ptr<const CoefficientEncoder<ModularInt>> encoder_;
  std::unique_ptr<const RnsErrorParams<ModularInt>> error_params_;
  std::unique_ptr<const PublicParameter<ModularInt>> public_parameter_;
  std::unique_ptr<Prng> prng_;
};

class RecoverNegativeTest : public RecoverTest<ModularInt32> {};

TEST_F(RecoverNegativeTest, RecoverMessagesFailsIfEmptyPartialDecryptions) {
  ASSERT_OK_AND_ASSIGN(RnsPolynomial<ModularInt32> fake_ct_b,
                       RnsPolynomial<ModularInt32>::CreateZero(
                           rns_context_->LogN(), moduli_, /*is_ntt=*/true));
  EXPECT_THAT(
      RecoverMessages<ModularInt32>(/*partial_decryptions=*/{}, fake_ct_b,
                                    *public_parameter_, encoder_.get()),
      StatusIs(absl::StatusCode::kInvalidArgument,
               HasSubstr("`partial_decryptions` must not be empty")));
}

TYPED_TEST_SUITE(RecoverTest, testing::ModularIntTypesForMultiParty);

// Simulate a multi-party protocol with `kNumParties` many parties.
TYPED_TEST(RecoverTest, RecoverSucceeds) {
  using Integer = typename TypeParam::Int;

  std::vector<SecretKeyShare<TypeParam>> secret_key_shares;
  std::vector<PublicKeyShare<TypeParam>> public_key_shares;
  secret_key_shares.reserve(kNumParties);
  public_key_shares.reserve(kNumParties);
  for (int i = 0; i < kNumParties; ++i) {
    ASSERT_OK_AND_ASSIGN(auto sk_share,
                         SecretKeyShare<TypeParam>::Sample(
                             this->rns_context_.get(), this->prng_.get()));
    ASSERT_OK_AND_ASSIGN(
        auto pk_share,
        PublicKeyShare<TypeParam>::Create(
            &sk_share, this->public_parameter_.get(), kPrngType));
    secret_key_shares.push_back(std::move(sk_share));
    public_key_shares.push_back(std::move(pk_share));
  }

  ASSERT_OK_AND_ASSIGN(auto public_key,
                       PublicKey<TypeParam>::Create(
                           this->public_parameter_.get(), public_key_shares));

  // Encrypt two vectors, and then homomorphically add the two ciphertexts.
  int num_coeffs = 1 << this->rns_context_->LogN();
  std::vector<Integer> coeffs0 = SampleMessages<Integer>(num_coeffs, kMaxValue);
  std::vector<Integer> coeffs1 = SampleMessages<Integer>(num_coeffs, kMaxValue);
  ASSERT_OK_AND_ASSIGN(
      auto ciphertext0,
      public_key.Encrypt(coeffs0, this->encoder_.get(),
                         this->error_params_.get(), this->prng_.get()));
  ASSERT_OK_AND_ASSIGN(
      auto ciphertext1,
      public_key.Encrypt(coeffs1, this->encoder_.get(),
                         this->error_params_.get(), this->prng_.get()));
  ASSERT_OK_AND_ASSIGN(auto ciphertext, ciphertext0 + ciphertext1);

  // Let all secret key share holders do partial decryption.
  ASSERT_OK_AND_ASSIGN(auto component_b, ciphertext.Component(0));
  ASSERT_OK_AND_ASSIGN(auto component_a, ciphertext.Component(1));
  ASSERT_OK_AND_ASSIGN(auto dg_sampler,
                       DiscreteGaussianSampler<Integer>::Create(kSBase));
  double s_flood = this->params_->s_flood;
  std::vector<RnsPolynomial<TypeParam>> partial_decryptions;
  partial_decryptions.reserve(kNumParties);
  for (int i = 0; i < kNumParties; ++i) {
    ASSERT_OK_AND_ASSIGN(
        auto partial_decryption,
        secret_key_shares[i].PartialDecrypt(
            component_a, s_flood, dg_sampler.get(), this->prng_.get()));
    partial_decryptions.push_back(std::move(partial_decryption));
  }

  // Recover the messages from the partial decryptions.
  ASSERT_OK_AND_ASSIGN(auto decrypted,
                       RecoverMessages<TypeParam>(
                           partial_decryptions, component_b,
                           *this->public_parameter_, this->encoder_.get()));
  ASSERT_EQ(decrypted.size(), num_coeffs);
  Integer t = this->encoder_.get()->PlaintextModulus();
  for (int i = 0; i < num_coeffs; ++i) {
    Integer decrypted_coeff = decrypted[i];
    Integer expected = (coeffs0[i] + coeffs1[i]) % t;
    EXPECT_EQ(decrypted_coeff, expected);
  }
}

}  // namespace
}  // namespace multi_party
}  // namespace rlwe
