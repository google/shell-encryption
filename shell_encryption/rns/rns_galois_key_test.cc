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

#include "shell_encryption/rns/rns_galois_key.h"

#include <cmath>
#include <cstddef>
#include <memory>
#include <vector>

#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include "absl/log/check.h"
#include "absl/status/status.h"
#include "absl/strings/string_view.h"
#include "shell_encryption/rns/finite_field_encoder.h"
#include "shell_encryption/rns/rns_bfv_ciphertext.h"
#include "shell_encryption/rns/rns_bgv_ciphertext.h"
#include "shell_encryption/rns/rns_context.h"
#include "shell_encryption/rns/rns_error_params.h"
#include "shell_encryption/rns/rns_gadget.h"
#include "shell_encryption/rns/rns_modulus.h"
#include "shell_encryption/rns/rns_polynomial.h"
#include "shell_encryption/rns/rns_secret_key.h"
#include "shell_encryption/rns/testing/parameters.h"
#include "shell_encryption/rns/testing/testing_utils.h"
#include "shell_encryption/testing/status_matchers.h"
#include "shell_encryption/testing/status_testing.h"
#include "shell_encryption/testing/testing_prng.h"

namespace rlwe {
namespace {

using ::rlwe::testing::StatusIs;
using ::testing::HasSubstr;
using Prng = testing::TestingPrng;

constexpr int kVariance = 8;
constexpr PrngType kPrngType = PRNG_TYPE_HKDF;
constexpr absl::string_view kPrngSeed =
    "0123456789abcdef0123456789abcdef0123456789abcdef0123456789abcdef";

// Test fixture for RnsGaloiosKey. We use parameters that support SIMD-style
// finite field encoding. See testing/parameters.h for concrete parameters used
// to instantiate these tests.
template <typename ModularInt>
class RnsGaloisKeyTest : public ::testing::Test {
 protected:
  using ModularIntParams = typename ModularInt::Params;
  using Integer = typename ModularInt::Int;
  using RnsParams = testing::RnsParameters<ModularInt>;

  void SetUp() override { prng_ = std::make_unique<testing::TestingPrng>(0); }

  void SetUpBgvParameters(const RnsParams& params) {
    // Use parameters suitable for finite field encoding.
    auto rns_context = RnsContext<ModularInt>::CreateForBgvFiniteFieldEncoding(
        params.log_n, params.qs, params.ps, params.t);
    CHECK(rns_context.ok());
    rns_context_ = std::make_unique<const RnsContext<ModularInt>>(
        std::move(rns_context.value()));
    main_moduli_ = rns_context_->MainPrimeModuli();

    int level = params.qs.size() - 1;
    auto q_hats = rns_context_->MainPrimeModulusComplements(level).value();
    auto q_hat_invs = rns_context_->MainPrimeModulusCrtFactors(level).value();

    std::vector<size_t> log_bs(params.qs.size(), params.log_gadget_base);
    auto gadget = RnsGadget<ModularInt>::Create(params.log_n, log_bs, q_hats,
                                                q_hat_invs, main_moduli_);
    CHECK(gadget.ok());
    gadget_ = std::make_unique<const RnsGadget<ModularInt>>(
        std::move(gadget.value()));

    auto error_params = RnsErrorParams<ModularInt>::Create(
        params.log_n, main_moduli_, /*aux_moduli=*/{},
        log2(static_cast<double>(params.t)), sqrt(kVariance));
    CHECK(error_params.ok());
    error_params_ = std::make_unique<const RnsErrorParams<ModularInt>>(
        std::move(error_params.value()));
  }

  void SetUpBfvParameters(const RnsParams& params) {
    // Use parameters suitable for finite field encoding in BFV.
    auto rns_context = RnsContext<ModularInt>::CreateForBfvFiniteFieldEncoding(
        params.log_n, params.qs, params.ps, params.t);
    CHECK(rns_context.ok());
    rns_context_ = std::make_unique<const RnsContext<ModularInt>>(
        std::move(rns_context.value()));
    main_moduli_ = rns_context_->MainPrimeModuli();
    aux_moduli_ = rns_context_->AuxPrimeModuli();

    int level = params.qs.size() - 1;
    auto q_hats = rns_context_->MainPrimeModulusComplements(level).value();
    auto q_hat_invs = rns_context_->MainPrimeModulusCrtFactors(level).value();

    std::vector<size_t> log_bs(params.qs.size(), params.log_gadget_base);
    auto gadget = RnsGadget<ModularInt>::Create(params.log_n, log_bs, q_hats,
                                                q_hat_invs, main_moduli_);
    CHECK(gadget.ok());
    gadget_ = std::make_unique<const RnsGadget<ModularInt>>(
        std::move(gadget.value()));

    auto error_params = RnsErrorParams<ModularInt>::Create(
        params.log_n, main_moduli_, /*aux_moduli=*/{},
        log2(static_cast<double>(params.t)), sqrt(kVariance));
    CHECK(error_params.ok());
    error_params_ = std::make_unique<const RnsErrorParams<ModularInt>>(
        std::move(error_params.value()));
  }

  RnsParams GetDefaultBgvParameters() const {
    auto all_params =
        testing::GetRnsParametersForFiniteFieldEncoding<ModularInt>();
    CHECK(!all_params.empty());
    return all_params[0];
  }

  void SetUpDefaultBgvParameters() {
    SetUpBgvParameters(GetDefaultBgvParameters());
  }

  RnsParams GetDefaultBfvParameters() const {
    auto all_params =
        testing::GetBfvParametersForFiniteFieldEncoding<ModularInt>();
    CHECK(!all_params.empty());
    return all_params[0];
  }

  void SetUpDefaultBfvParameters() {
    SetUpBgvParameters(GetDefaultBfvParameters());
  }

  // Sample a secret key from prng.
  absl::StatusOr<RnsRlweSecretKey<ModularInt>> SampleSecretKey() const {
    return RnsRlweSecretKey<ModularInt>::Sample(rns_context_->LogN(), kVariance,
                                                main_moduli_, prng_.get());
  }

  std::unique_ptr<const RnsContext<ModularInt>> rns_context_;
  std::vector<const PrimeModulus<ModularInt>*> main_moduli_;
  std::vector<const PrimeModulus<ModularInt>*> aux_moduli_;
  std::unique_ptr<const RnsGadget<ModularInt>> gadget_;
  std::unique_ptr<const RnsErrorParams<ModularInt>> error_params_;
  std::unique_ptr<Prng> prng_;
};

TYPED_TEST_SUITE(RnsGaloisKeyTest,
                 testing::ModularIntTypesForFiniteFieldEncoding);

TYPED_TEST(RnsGaloisKeyTest, CreateFailsIfPowerIsOutOfRange) {
  this->SetUpDefaultBgvParameters();
  ASSERT_OK_AND_ASSIGN(RnsRlweSecretKey<TypeParam> secret_key,
                       this->SampleSecretKey());

  // Negative power.
  EXPECT_THAT(
      RnsGaloisKey<TypeParam>::CreateForBgv(
          secret_key, /*power=*/-1, kVariance, this->gadget_.get(),
          this->rns_context_->PlaintextModulus(), kPrngType),
      StatusIs(absl::StatusCode::kInvalidArgument,
               HasSubstr("`power` must be a non-negative odd integer")));
  // Zero power.
  EXPECT_THAT(
      RnsGaloisKey<TypeParam>::CreateForBgv(
          secret_key, /*power=*/0, kVariance, this->gadget_.get(),
          this->rns_context_->PlaintextModulus(), kPrngType),
      StatusIs(absl::StatusCode::kInvalidArgument,
               HasSubstr("`power` must be a non-negative odd integer")));
  // Positive small even power.
  EXPECT_THAT(
      RnsGaloisKey<TypeParam>::CreateForBgv(
          secret_key, /*power=*/2, kVariance, this->gadget_.get(),
          this->rns_context_->PlaintextModulus(), kPrngType),
      StatusIs(absl::StatusCode::kInvalidArgument,
               HasSubstr("`power` must be a non-negative odd integer")));
  // Positive large power that is out of range.
  int num_coeffs = 1 << this->rns_context_->LogN();
  EXPECT_THAT(
      RnsGaloisKey<TypeParam>::CreateForBgv(
          secret_key, /*power=*/2 * num_coeffs + 1, kVariance,
          this->gadget_.get(), this->rns_context_->PlaintextModulus(),
          kPrngType),
      StatusIs(absl::StatusCode::kInvalidArgument,
               HasSubstr("`power` must be a non-negative odd integer")));
}

TYPED_TEST(RnsGaloisKeyTest, CreateFailsIfVarianceIsNotPositive) {
  this->SetUpDefaultBgvParameters();
  ASSERT_OK_AND_ASSIGN(RnsRlweSecretKey<TypeParam> secret_key,
                       this->SampleSecretKey());

  EXPECT_THAT(RnsGaloisKey<TypeParam>::CreateForBgv(
                  secret_key, /*power=*/5, /*variance=*/-1, this->gadget_.get(),
                  this->rns_context_->PlaintextModulus(), kPrngType),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`variance` must be positive")));
  EXPECT_THAT(RnsGaloisKey<TypeParam>::CreateForBgv(
                  secret_key, /*power=*/5, /*variance=*/0, this->gadget_.get(),
                  this->rns_context_->PlaintextModulus(), kPrngType),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`variance` must be positive")));
}

TYPED_TEST(RnsGaloisKeyTest, CreateFailsIfGadgetIsNull) {
  this->SetUpDefaultBgvParameters();
  ASSERT_OK_AND_ASSIGN(RnsRlweSecretKey<TypeParam> secret_key,
                       this->SampleSecretKey());

  EXPECT_THAT(RnsGaloisKey<TypeParam>::CreateForBgv(
                  secret_key, /*power=*/5, kVariance, /*gadget=*/nullptr,
                  this->rns_context_->PlaintextModulus(), kPrngType),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`gadget` must not be null")));
}

TYPED_TEST(RnsGaloisKeyTest, CreateFailsIfPrngTypeIsInvalid) {
  this->SetUpDefaultBgvParameters();
  ASSERT_OK_AND_ASSIGN(RnsRlweSecretKey<TypeParam> secret_key,
                       this->SampleSecretKey());

  EXPECT_THAT(RnsGaloisKey<TypeParam>::CreateForBgv(
                  secret_key, /*power=*/5, kVariance, this->gadget_.get(),
                  this->rns_context_->PlaintextModulus(), PRNG_TYPE_INVALID),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("PrngType not specified correctly")));
}

TYPED_TEST(RnsGaloisKeyTest, CreateFromKeyComponentsFailsIfIncorrectDimension) {
  this->SetUpDefaultBgvParameters();

  // Mock up the two component vectors using zero polynomial.
  int log_n = this->rns_context_->LogN();
  int dimension = this->gadget_->Dimension();
  ASSERT_OK_AND_ASSIGN(
      RnsPolynomial<TypeParam> zero,
      RnsPolynomial<TypeParam>::CreateZero(log_n, this->main_moduli_));
  std::vector<RnsPolynomial<TypeParam>> keys_wrong_size;
  std::vector<RnsPolynomial<TypeParam>> keys_correct_size(dimension, zero);

  EXPECT_THAT(
      RnsGaloisKey<TypeParam>::CreateFromKeyComponents(
          keys_wrong_size, keys_correct_size, /*power=*/5, this->gadget_.get(),
          this->main_moduli_, kPrngSeed, kPrngType),
      StatusIs(absl::StatusCode::kInvalidArgument,
               HasSubstr("must have the same length as the gadget dimension")));
  EXPECT_THAT(
      RnsGaloisKey<TypeParam>::CreateFromKeyComponents(
          keys_correct_size, keys_wrong_size, /*power=*/5, this->gadget_.get(),
          this->main_moduli_, kPrngSeed, kPrngType),
      StatusIs(absl::StatusCode::kInvalidArgument,
               HasSubstr("must have the same length as the gadget dimension")));
}

TYPED_TEST(RnsGaloisKeyTest, DeserializeFailsIfGadgetIsNull) {
  constexpr int k_power = 5;
  this->SetUpDefaultBfvParameters();
  ASSERT_OK_AND_ASSIGN(RnsRlweSecretKey<TypeParam> secret_key,
                       this->SampleSecretKey());
  ASSERT_OK_AND_ASSIGN(
      RnsGaloisKey<TypeParam> gk,
      RnsGaloisKey<TypeParam>::CreateForBfv(secret_key, k_power, kVariance,
                                            this->gadget_.get(), kPrngType));
  ASSERT_OK_AND_ASSIGN(SerializedRnsGaloisKey serialized, gk.Serialize());

  EXPECT_THAT(RnsGaloisKey<TypeParam>::Deserialize(
                  serialized, /*gadget=*/nullptr, this->main_moduli_),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`gadget` must not be null")));
}

TYPED_TEST(RnsGaloisKeyTest, DeserializeFailsIfKeybsIsEmpty) {
  constexpr int k_power = 5;
  this->SetUpDefaultBfvParameters();
  ASSERT_OK_AND_ASSIGN(RnsRlweSecretKey<TypeParam> secret_key,
                       this->SampleSecretKey());
  ASSERT_OK_AND_ASSIGN(
      RnsGaloisKey<TypeParam> gk,
      RnsGaloisKey<TypeParam>::CreateForBfv(secret_key, k_power, kVariance,
                                            this->gadget_.get(), kPrngType));
  ASSERT_OK_AND_ASSIGN(SerializedRnsGaloisKey serialized, gk.Serialize());
  serialized.clear_key_bs();

  EXPECT_THAT(RnsGaloisKey<TypeParam>::Deserialize(
                  serialized, this->gadget_.get(), this->main_moduli_),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`key_bs` must not be empty")));
}

TYPED_TEST(RnsGaloisKeyTest, ApplyToFailsIfCiphertextHasMismatchPowerOfS) {
  constexpr int k_gk_power = 5;
  constexpr int k_ctxt_power = 1;

  this->SetUpDefaultBgvParameters();
  ASSERT_OK_AND_ASSIGN(RnsRlweSecretKey<TypeParam> secret_key,
                       this->SampleSecretKey());

  ASSERT_OK_AND_ASSIGN(
      RnsGaloisKey<TypeParam> gk,
      RnsGaloisKey<TypeParam>::CreateForBgv(
          secret_key, k_gk_power, kVariance, this->gadget_.get(),
          this->rns_context_->PlaintextModulus(), kPrngType));

  // Create a zero ciphertext (0, 0) for negative testing only.
  int log_n = this->rns_context_->LogN();
  ASSERT_OK_AND_ASSIGN(
      RnsPolynomial<TypeParam> zero,
      RnsPolynomial<TypeParam>::CreateZero(log_n, this->main_moduli_,
                                           /*is_ntt=*/false));
  RnsBgvCiphertext<TypeParam> ciphertext({zero, zero}, this->main_moduli_,
                                         k_ctxt_power, /*error=*/0,
                                         this->error_params_.get());

  // ApplyTo will fail since ciphertext has a mismatch PowerOfS.
  EXPECT_THAT(
      gk.ApplyTo(ciphertext),
      StatusIs(
          absl::StatusCode::kInvalidArgument,
          HasSubstr(
              "`ciphertext` does not have a matching substitution power")));
}

TYPED_TEST(RnsGaloisKeyTest, ApplyToFailsIfCiphertextHasMismatchRnsModuli) {
  // This test requires at least two RNS moduli.
  auto params = this->GetDefaultBgvParameters();
  if (params.qs.size() < 2) {
    return;
  }

  constexpr int k_power = 5;
  this->SetUpBgvParameters(params);

  ASSERT_OK_AND_ASSIGN(RnsRlweSecretKey<TypeParam> secret_key,
                       this->SampleSecretKey());
  ASSERT_OK_AND_ASSIGN(RnsGaloisKey<TypeParam> gk,
                       RnsGaloisKey<TypeParam>::CreateForBgv(
                           secret_key, k_power, kVariance, this->gadget_.get(),
                           this->rns_context_->PlaintextModulus(), kPrngType));

  // Create a zero ciphertext (0, 0) with a smaller set of RNS moduli.
  int log_n = this->rns_context_->LogN();
  std::vector<const PrimeModulus<TypeParam>*> small_moduli = this->main_moduli_;
  small_moduli.pop_back();
  ASSERT_OK_AND_ASSIGN(RnsPolynomial<TypeParam> zero,
                       RnsPolynomial<TypeParam>::CreateZero(log_n, small_moduli,
                                                            /*is_ntt=*/false));
  RnsBgvCiphertext<TypeParam> ciphertext({zero, zero}, small_moduli, k_power,
                                         /*error=*/0,
                                         this->error_params_.get());

  // ApplyTo will fail since ciphertrext has a mismatch RNS moduli set.
  EXPECT_THAT(
      gk.ApplyTo(ciphertext),
      StatusIs(
          absl::StatusCode::kInvalidArgument,
          HasSubstr("`ciphertext` does not have a matching RNS moduli set")));
}

TYPED_TEST(RnsGaloisKeyTest, ApplyToFailsIfCiphertextDegreeLargerThanOne) {
  constexpr int k_power = 5;

  this->SetUpDefaultBgvParameters();
  ASSERT_OK_AND_ASSIGN(RnsRlweSecretKey<TypeParam> secret_key,
                       this->SampleSecretKey());

  ASSERT_OK_AND_ASSIGN(RnsGaloisKey<TypeParam> gk,
                       RnsGaloisKey<TypeParam>::CreateForBgv(
                           secret_key, k_power, kVariance, this->gadget_.get(),
                           this->rns_context_->PlaintextModulus(), kPrngType));

  // Create a degree-2 ciphertext (0, 0, 0).
  int log_n = this->rns_context_->LogN();
  ASSERT_OK_AND_ASSIGN(
      RnsPolynomial<TypeParam> zero,
      RnsPolynomial<TypeParam>::CreateZero(log_n, this->main_moduli_,
                                           /*is_ntt=*/false));
  RnsBgvCiphertext<TypeParam> ciphertext({zero, zero, zero}, this->main_moduli_,
                                         k_power, /*error=*/0,
                                         this->error_params_.get());

  // ApplyTo will fail since ciphertrext has a mismatch degree.
  EXPECT_THAT(
      gk.ApplyTo(ciphertext),
      StatusIs(
          absl::StatusCode::kInvalidArgument,
          HasSubstr("Galois key can only apply to a ciphertext of degree 1")));
}

TYPED_TEST(RnsGaloisKeyTest, GaloisKeyForBgvCanRotateEncryptedSlots) {
  using Encoder = FiniteFieldEncoder<TypeParam>;
  using Integer = typename TypeParam::Int;

  constexpr int k_power = 5;  // rotation by one slot.

  for (auto const& params :
       testing::GetRnsParametersForFiniteFieldEncoding<TypeParam>()) {
    this->SetUpBgvParameters(params);

    ASSERT_OK_AND_ASSIGN(auto encoder,
                         Encoder::Create(this->rns_context_.get()));

    ASSERT_OK_AND_ASSIGN(RnsRlweSecretKey<TypeParam> secret_key,
                         this->SampleSecretKey());

    auto gadget = this->gadget_.get();
    auto prng = this->prng_.get();
    Integer t = this->rns_context_->PlaintextModulus();
    ASSERT_OK_AND_ASSIGN(
        RnsGaloisKey<TypeParam> gk,
        RnsGaloisKey<TypeParam>::CreateForBgv(secret_key, k_power, kVariance,
                                              gadget, t, kPrngType));
    ASSERT_EQ(gk.Gadget(), gadget);
    ASSERT_EQ(gk.Dimension(), gadget->Dimension());
    ASSERT_EQ(gk.SubstitutionPower(), k_power);

    int num_coeffs = 1 << this->rns_context_->LogN();
    std::vector<Integer> values =
        testing::SampleMessages<Integer>(num_coeffs, t);
    ASSERT_OK_AND_ASSIGN(
        RnsBgvCiphertext<TypeParam> ciphertext,
        secret_key.template EncryptBgv<Encoder>(
            values, &encoder, this->error_params_.get(), prng));

    // Apply the automorphism X -> X^k_power to the ciphertext, which is now
    // encrypted under a secret key (1, s(X^k_power)). We use the Galois key
    // to transform this ciphertext back to be encrypted under the canonical
    // secret key (1, s(X)).
    ASSERT_OK_AND_ASSIGN(auto ciphertext_sub, ciphertext.Substitute(k_power));
    ASSERT_OK_AND_ASSIGN(auto ciphertext_rot, gk.ApplyTo(ciphertext_sub));

    ASSERT_OK_AND_ASSIGN(
        auto decryptions,
        secret_key.template DecryptBgv<Encoder>(ciphertext_rot, &encoder));

    // The slots are divided into two groups, and by applying the Galois key
    // each group of the encrypted slots are rotated by 1 position.
    int group_size = num_coeffs / 2;
    std::vector<Integer> expected;
    for (int i = 0; i < group_size; ++i) {
      int index = (i + 1) % group_size;
      expected.push_back(values[index]);
    }
    for (int i = 0; i < group_size; ++i) {
      int index = group_size + (i + 1) % group_size;
      expected.push_back(values[index]);
    }

    EXPECT_EQ(decryptions, expected);
  }
}

// Checks that Galois key can be constructed using all supported PRNG types.
TYPED_TEST(RnsGaloisKeyTest, GaloisKeyCanBeInstantiatedUsingPrngs) {
  using Encoder = FiniteFieldEncoder<TypeParam>;
  using Integer = typename TypeParam::Int;

  constexpr int k_power = 5;  // rotation by one slot.

  // Using the default parameter for each integer type is sufficient.
  this->SetUpDefaultBgvParameters();
  // Keep the same secret key and a ciphertext that encrypts values in slots.
  ASSERT_OK_AND_ASSIGN(RnsRlweSecretKey<TypeParam> secret_key,
                       this->SampleSecretKey());
  ASSERT_OK_AND_ASSIGN(auto encoder, Encoder::Create(this->rns_context_.get()));
  Integer t = this->rns_context_->PlaintextModulus();
  int num_coeffs = 1 << this->rns_context_->LogN();
  std::vector<Integer> values = testing::SampleMessages<Integer>(num_coeffs, t);
  ASSERT_OK_AND_ASSIGN(
      RnsBgvCiphertext<TypeParam> ciphertext,
      secret_key.template EncryptBgv<Encoder>(
          values, &encoder, this->error_params_.get(), this->prng_.get()));

  // Apply the automorphism X -> X^k_power to the ciphertext, which is now
  // encrypted under a secret key (1, s(X^k_power)). We use the Galois key
  // to transform this ciphertext back to be encrypted under the canonical
  // secret key (1, s(X)).
  ASSERT_OK_AND_ASSIGN(auto ciphertext_sub, ciphertext.Substitute(k_power));

  for (auto prng_type : {PRNG_TYPE_CHACHA, PRNG_TYPE_HKDF}) {
    // Create a Galois key from each PRNG type, and apply it to ciphertext_sub.
    ASSERT_OK_AND_ASSIGN(
        RnsGaloisKey<TypeParam> gk,
        RnsGaloisKey<TypeParam>::CreateForBgv(
            secret_key, k_power, kVariance, this->gadget_.get(), t, prng_type));
    ASSERT_OK_AND_ASSIGN(auto ciphertext_rot, gk.ApplyTo(ciphertext_sub));

    // The resulting ciphertext now encrypts a rotated slot vector.
    ASSERT_OK_AND_ASSIGN(
        auto decryptions,
        secret_key.template DecryptBgv<Encoder>(ciphertext_rot, &encoder));

    int group_size = num_coeffs / 2;
    std::vector<Integer> expected;
    for (int i = 0; i < group_size; ++i) {
      int index = (i + 1) % group_size;
      expected.push_back(values[index]);
    }
    for (int i = 0; i < group_size; ++i) {
      int index = group_size + (i + 1) % group_size;
      expected.push_back(values[index]);
    }

    EXPECT_EQ(decryptions, expected);
  }
}

// Checks that applying two Galois keys with substitution powers a, b, resp,
// on a BGV ciphertext is equivalent to substituting X with X^(ab) on plaintext.
TYPED_TEST(RnsGaloisKeyTest, GaloisKeyForBgvComposedApplication) {
  using Encoder = FiniteFieldEncoder<TypeParam>;
  using Integer = typename TypeParam::Int;

  constexpr int k_power0 = 5;   // rotation by one slot.
  constexpr int k_power1 = 25;  // rotation by two slots.

  for (auto const& params :
       testing::GetRnsParametersForFiniteFieldEncoding<TypeParam>()) {
    this->SetUpBgvParameters(params);

    ASSERT_OK_AND_ASSIGN(auto encoder,
                         Encoder::Create(this->rns_context_.get()));

    ASSERT_OK_AND_ASSIGN(RnsRlweSecretKey<TypeParam> secret_key,
                         this->SampleSecretKey());

    auto prng = this->prng_.get();
    Integer t = this->rns_context_->PlaintextModulus();

    // Create two Galois keys with different substitution powers.
    ASSERT_OK_AND_ASSIGN(RnsGaloisKey<TypeParam> gk0,
                         RnsGaloisKey<TypeParam>::CreateForBgv(
                             secret_key, k_power0, kVariance,
                             this->gadget_.get(), t, kPrngType));
    ASSERT_OK_AND_ASSIGN(RnsGaloisKey<TypeParam> gk1,
                         RnsGaloisKey<TypeParam>::CreateForBgv(
                             secret_key, k_power1, kVariance,
                             this->gadget_.get(), t, kPrngType));

    int num_coeffs = 1 << this->rns_context_->LogN();
    std::vector<Integer> values =
        testing::SampleMessages<Integer>(num_coeffs, t);
    ASSERT_OK_AND_ASSIGN(
        RnsBgvCiphertext<TypeParam> ciphertext,
        secret_key.template EncryptBgv<Encoder>(
            values, &encoder, this->error_params_.get(), prng));

    // Apply the first automorphism to the ciphertext and then apply the first
    // Galois key.
    ASSERT_OK_AND_ASSIGN(auto ciphertext_sub0, ciphertext.Substitute(k_power0));
    ASSERT_OK_AND_ASSIGN(auto ciphertext_rot0, gk0.ApplyTo(ciphertext_sub0));
    ASSERT_EQ(ciphertext_rot0.PowerOfS(), 1);

    // Apply the second automorphism to the new ciphertext and then apply the
    // second Galois key.
    ASSERT_OK_AND_ASSIGN(auto ciphertext_sub1,
                         ciphertext_rot0.Substitute(k_power1));
    ASSERT_OK_AND_ASSIGN(auto ciphertext_rot1, gk1.ApplyTo(ciphertext_sub1));
    ASSERT_EQ(ciphertext_rot1.PowerOfS(), 1);

    // Now the last ciphertext should encrypt a slot vector that is rotated by
    // the combined offset.
    ASSERT_OK_AND_ASSIGN(
        auto decryptions,
        secret_key.template DecryptBgv<Encoder>(ciphertext_rot1, &encoder));

    int group_size = num_coeffs / 2;
    std::vector<Integer> expected;
    for (int i = 0; i < group_size; ++i) {
      int index = (i + 3) % group_size;
      expected.push_back(values[index]);
    }
    for (int i = 0; i < group_size; ++i) {
      int index = group_size + (i + 3) % group_size;
      expected.push_back(values[index]);
    }

    EXPECT_EQ(decryptions, expected);
  }
}

TYPED_TEST(RnsGaloisKeyTest, GaloisKeyForBfvCanRotateEncryptedSlots) {
  using Encoder = FiniteFieldEncoder<TypeParam>;
  using Integer = typename TypeParam::Int;

  constexpr int k_power = 5;  // rotation by one slot.

  for (auto const& params :
       testing::GetBfvParametersForFiniteFieldEncoding<TypeParam>()) {
    this->SetUpBfvParameters(params);

    ASSERT_OK_AND_ASSIGN(auto encoder,
                         Encoder::Create(this->rns_context_.get()));

    ASSERT_OK_AND_ASSIGN(RnsRlweSecretKey<TypeParam> secret_key,
                         this->SampleSecretKey());

    auto gadget = this->gadget_.get();
    auto prng = this->prng_.get();
    Integer t = this->rns_context_->PlaintextModulus();
    ASSERT_OK_AND_ASSIGN(
        RnsGaloisKey<TypeParam> gk,
        RnsGaloisKey<TypeParam>::CreateForBfv(secret_key, k_power, kVariance,
                                              gadget, kPrngType));
    ASSERT_EQ(gk.Gadget(), gadget);
    ASSERT_EQ(gk.Dimension(), gadget->Dimension());
    ASSERT_EQ(gk.SubstitutionPower(), k_power);

    int num_coeffs = 1 << this->rns_context_->LogN();
    std::vector<Integer> values =
        testing::SampleMessages<Integer>(num_coeffs, t);
    ASSERT_OK_AND_ASSIGN(
        RnsBfvCiphertext<TypeParam> ciphertext,
        secret_key.template EncryptBfv<Encoder>(
            values, &encoder, this->error_params_.get(), prng));

    // Apply the automorphism X -> X^k_power to the ciphertext, which is now
    // encrypted under a secret key (1, s(X^k_power)). We use the Galois key
    // to transform this ciphertext back to be encrypted under the canonical
    // secret key (1, s(X)).
    ASSERT_OK_AND_ASSIGN(auto ciphertext_sub, ciphertext.Substitute(k_power));
    ASSERT_OK_AND_ASSIGN(auto ciphertext_rot, gk.ApplyTo(ciphertext_sub));

    ASSERT_OK_AND_ASSIGN(
        auto decryptions,
        secret_key.template DecryptBfv<Encoder>(ciphertext_rot, &encoder));

    // The slots are divided into two groups, and by applying the Galois key
    // each group of the encrypted slots are rotated by 1 position.
    int group_size = num_coeffs / 2;
    std::vector<Integer> expected;
    for (int i = 0; i < group_size; ++i) {
      int index = (i + 1) % group_size;
      expected.push_back(values[index]);
    }
    for (int i = 0; i < group_size; ++i) {
      int index = group_size + (i + 1) % group_size;
      expected.push_back(values[index]);
    }

    EXPECT_EQ(decryptions, expected);
  }
}

// Checks that applying two Galois keys with substitution powers a, b, resp,
// on a BFV ciphertext is equivalent to substituting X with X^(ab) on plaintext.
TYPED_TEST(RnsGaloisKeyTest, GaloisKeyForBfvComposedApplication) {
  using Encoder = FiniteFieldEncoder<TypeParam>;
  using Integer = typename TypeParam::Int;

  constexpr int k_power0 = 5;   // rotation by one slot.
  constexpr int k_power1 = 25;  // rotation by two slots.

  for (auto const& params :
       testing::GetBfvParametersForFiniteFieldEncoding<TypeParam>()) {
    this->SetUpBfvParameters(params);

    ASSERT_OK_AND_ASSIGN(auto encoder,
                         Encoder::Create(this->rns_context_.get()));

    ASSERT_OK_AND_ASSIGN(RnsRlweSecretKey<TypeParam> secret_key,
                         this->SampleSecretKey());

    // Create two Galois keys with different substitution powers.
    ASSERT_OK_AND_ASSIGN(
        RnsGaloisKey<TypeParam> gk0,
        RnsGaloisKey<TypeParam>::CreateForBfv(secret_key, k_power0, kVariance,
                                              this->gadget_.get(), kPrngType));
    ASSERT_OK_AND_ASSIGN(
        RnsGaloisKey<TypeParam> gk1,
        RnsGaloisKey<TypeParam>::CreateForBfv(secret_key, k_power1, kVariance,
                                              this->gadget_.get(), kPrngType));

    int num_coeffs = 1 << this->rns_context_->LogN();
    Integer t = this->rns_context_->PlaintextModulus();
    std::vector<Integer> values =
        testing::SampleMessages<Integer>(num_coeffs, t);
    auto prng = this->prng_.get();
    ASSERT_OK_AND_ASSIGN(
        RnsBfvCiphertext<TypeParam> ciphertext,
        secret_key.template EncryptBfv<Encoder>(
            values, &encoder, this->error_params_.get(), prng));

    // Apply the first automorphism to the ciphertext and then apply the first
    // Galois key.
    ASSERT_OK_AND_ASSIGN(auto ciphertext_sub0, ciphertext.Substitute(k_power0));
    ASSERT_OK_AND_ASSIGN(auto ciphertext_rot0, gk0.ApplyTo(ciphertext_sub0));
    ASSERT_EQ(ciphertext_rot0.PowerOfS(), 1);

    // Apply the second automorphism to the new ciphertext and then apply the
    // second Galois key.
    ASSERT_OK_AND_ASSIGN(auto ciphertext_sub1,
                         ciphertext_rot0.Substitute(k_power1));
    ASSERT_OK_AND_ASSIGN(auto ciphertext_rot1, gk1.ApplyTo(ciphertext_sub1));
    ASSERT_EQ(ciphertext_rot1.PowerOfS(), 1);

    // Now the last ciphertext should encrypt a slot vector that is rotated by
    // the combined offset.
    ASSERT_OK_AND_ASSIGN(
        auto decryptions,
        secret_key.template DecryptBfv<Encoder>(ciphertext_rot1, &encoder));

    int group_size = num_coeffs / 2;
    std::vector<Integer> expected;
    for (int i = 0; i < group_size; ++i) {
      int index = (i + 3) % group_size;
      expected.push_back(values[index]);
    }
    for (int i = 0; i < group_size; ++i) {
      int index = group_size + (i + 3) % group_size;
      expected.push_back(values[index]);
    }

    EXPECT_EQ(decryptions, expected);
  }
}

TYPED_TEST(RnsGaloisKeyTest, PrecomputeRandomPadsFailsIfDigitsHaveWrongSize) {
  constexpr int k_power = 5;

  this->SetUpDefaultBgvParameters();
  ASSERT_OK_AND_ASSIGN(RnsRlweSecretKey<TypeParam> secret_key,
                       this->SampleSecretKey());
  ASSERT_OK_AND_ASSIGN(RnsGaloisKey<TypeParam> gk,
                       RnsGaloisKey<TypeParam>::CreateForBgv(
                           secret_key, k_power, kVariance, this->gadget_.get(),
                           this->rns_context_->PlaintextModulus(), kPrngType));

  // Create a vector of dimension - 1 polynomials.
  int dimension = this->gadget_->Dimension();
  int log_n = this->rns_context_->LogN();
  ASSERT_OK_AND_ASSIGN(
      RnsPolynomial<TypeParam> zero,
      RnsPolynomial<TypeParam>::CreateZero(log_n, this->main_moduli_,
                                           /*is_ntt=*/false));
  std::vector<RnsPolynomial<TypeParam>> digits_wrong_size(dimension - 1, zero);
  EXPECT_THAT(
      gk.PrecomputeRandomPad(digits_wrong_size),
      StatusIs(absl::StatusCode::kInvalidArgument,
               HasSubstr("`ciphertext_pad_digits` has incorrect size")));
}

TYPED_TEST(RnsGaloisKeyTest, ApplyWithRandomPadFailsIfMismatchPowerOfS) {
  constexpr int k_gk_power = 5;
  constexpr int k_ctxt_power = 7;

  this->SetUpDefaultBgvParameters();
  ASSERT_OK_AND_ASSIGN(RnsRlweSecretKey<TypeParam> secret_key,
                       this->SampleSecretKey());
  ASSERT_OK_AND_ASSIGN(
      RnsGaloisKey<TypeParam> gk,
      RnsGaloisKey<TypeParam>::CreateForBgv(
          secret_key, k_gk_power, kVariance, this->gadget_.get(),
          this->rns_context_->PlaintextModulus(), kPrngType));

  // Create a zero ciphertext with a different power of s.
  int log_n = this->rns_context_->LogN();
  ASSERT_OK_AND_ASSIGN(
      RnsPolynomial<TypeParam> zero,
      RnsPolynomial<TypeParam>::CreateZero(log_n, this->main_moduli_,
                                           /*is_ntt=*/false));
  RnsBgvCiphertext<TypeParam> ciphertext({zero, zero}, this->main_moduli_,
                                         k_ctxt_power, /*error=*/0,
                                         this->error_params_.get());

  // Vector of decomposed digits of the "a" part of ciphertext.
  int dimension = this->gadget_->Dimension();
  std::vector<RnsPolynomial<TypeParam>> digits(dimension, zero);

  EXPECT_THAT(
      gk.ApplyToWithRandomPad(ciphertext, digits, zero),
      StatusIs(
          absl::StatusCode::kInvalidArgument,
          HasSubstr(
              "`ciphertext` does not have a matching substitution power")));
}

TYPED_TEST(RnsGaloisKeyTest, ApplyWithRandomPadFailsIfMismatchRnsModuli) {
  // This test requires at least RNS moduli.
  auto params = this->GetDefaultBgvParameters();
  if (params.qs.size() < 2) {
    return;
  }

  constexpr int k_power = 5;
  this->SetUpBgvParameters(params);
  ASSERT_OK_AND_ASSIGN(RnsRlweSecretKey<TypeParam> secret_key,
                       this->SampleSecretKey());
  ASSERT_OK_AND_ASSIGN(RnsGaloisKey<TypeParam> gk,
                       RnsGaloisKey<TypeParam>::CreateForBgv(
                           secret_key, k_power, kVariance, this->gadget_.get(),
                           this->rns_context_->PlaintextModulus(), kPrngType));

  // Create a zero ciphertext (0, 0) with a smaller set of RNS moduli.
  int log_n = this->rns_context_->LogN();
  std::vector<const PrimeModulus<TypeParam>*> small_moduli = this->main_moduli_;
  small_moduli.pop_back();
  ASSERT_OK_AND_ASSIGN(RnsPolynomial<TypeParam> zero,
                       RnsPolynomial<TypeParam>::CreateZero(log_n, small_moduli,
                                                            /*is_ntt=*/false));
  RnsBgvCiphertext<TypeParam> ciphertext({zero, zero}, small_moduli, k_power,
                                         /*error=*/0,
                                         this->error_params_.get());

  // Vector of decomposed digits of the "a" part of ciphertext.
  int dimension = this->gadget_->Dimension();
  std::vector<RnsPolynomial<TypeParam>> digits(dimension, zero);

  EXPECT_THAT(
      gk.ApplyToWithRandomPad(ciphertext, digits, zero),
      StatusIs(
          absl::StatusCode::kInvalidArgument,
          HasSubstr("`ciphertext` does not have a matching RNS moduli set")));
}

TYPED_TEST(RnsGaloisKeyTest, GaloisKeyForBgvWithPrecomputedPads) {
  using Encoder = FiniteFieldEncoder<TypeParam>;
  using Integer = typename TypeParam::Int;

  constexpr int k_power = 5;  // rotation by one slot.

  for (auto const& params :
       testing::GetRnsParametersForFiniteFieldEncoding<TypeParam>()) {
    this->SetUpBgvParameters(params);

    ASSERT_OK_AND_ASSIGN(auto encoder,
                         Encoder::Create(this->rns_context_.get()));

    ASSERT_OK_AND_ASSIGN(RnsRlweSecretKey<TypeParam> secret_key,
                         this->SampleSecretKey());

    auto gadget = this->gadget_.get();
    int k = gadget->Dimension();
    int log_n = this->rns_context_->LogN();
    ASSERT_OK_AND_ASSIGN(
        std::vector<RnsPolynomial<TypeParam>> key_as,
        RnsGaloisKey<TypeParam>::SampleRandomPad(k, log_n, this->main_moduli_,
                                                 kPrngSeed, kPrngType));

    Integer t = this->rns_context_->PlaintextModulus();
    ASSERT_OK_AND_ASSIGN(RnsGaloisKey<TypeParam> gk,
                         RnsGaloisKey<TypeParam>::CreateWithRandomPadForBgv(
                             std::move(key_as), secret_key, k_power, kVariance,
                             gadget, t, kPrngSeed, kPrngType));
    ASSERT_EQ(gk.Gadget(), gadget);
    ASSERT_EQ(gk.Dimension(), gadget->Dimension());
    ASSERT_EQ(gk.SubstitutionPower(), k_power);

    int num_coeffs = 1 << log_n;
    std::vector<Integer> values =
        testing::SampleMessages<Integer>(num_coeffs, t);
    ASSERT_OK_AND_ASSIGN(
        RnsBgvCiphertext<TypeParam> ciphertext,
        secret_key.template EncryptBgv<Encoder>(
            values, &encoder, this->error_params_.get(), this->prng_.get()));

    // Precompute the gadget decomposition of the "a" component of `ciphertext`,
    // and also precompute the "a" component of the key-switched ciphertext.
    ASSERT_OK_AND_ASSIGN(RnsPolynomial<TypeParam> c1, ciphertext.Component(1));
    ASSERT_OK_AND_ASSIGN(RnsPolynomial<TypeParam> c1_sub,
                         c1.Substitute(k_power, this->main_moduli_));
    ASSERT_OK(c1_sub.ConvertToCoeffForm(this->main_moduli_));
    ASSERT_OK_AND_ASSIGN(std::vector<RnsPolynomial<TypeParam>> c1_sub_digits,
                         gadget->Decompose(c1_sub, this->main_moduli_));
    for (auto& c1_sub_digit : c1_sub_digits) {
      ASSERT_OK(c1_sub_digit.ConvertToNttForm(this->main_moduli_));
    }
    ASSERT_OK_AND_ASSIGN(RnsPolynomial<TypeParam> c1_rot,
                         gk.PrecomputeRandomPad(c1_sub_digits));

    // Apply the automorphism X -> X^k_power to the ciphertext with precomputed
    // "a" components of the resulting ciphertext.
    ASSERT_OK_AND_ASSIGN(auto ciphertext_sub, ciphertext.Substitute(k_power));
    ASSERT_OK_AND_ASSIGN(auto ciphertext_rot,
                         gk.ApplyToWithRandomPad(ciphertext_sub, c1_sub_digits,
                                                 std::move(c1_rot)));

    // Decrypt and check results of homomorphic rotation.
    ASSERT_OK_AND_ASSIGN(
        auto decryptions,
        secret_key.template DecryptBgv<Encoder>(ciphertext_rot, &encoder));
    int group_size = num_coeffs / 2;
    std::vector<Integer> expected;
    for (int i = 0; i < group_size; ++i) {
      int index = (i + 1) % group_size;
      expected.push_back(values[index]);
    }
    for (int i = 0; i < group_size; ++i) {
      int index = group_size + (i + 1) % group_size;
      expected.push_back(values[index]);
    }

    EXPECT_EQ(decryptions, expected);
  }
}

TYPED_TEST(RnsGaloisKeyTest, GaloisKeyForBfvWithPrecomputedPads) {
  using Encoder = FiniteFieldEncoder<TypeParam>;
  using Integer = typename TypeParam::Int;

  constexpr int k_power = 5;  // rotation by one slot.

  for (auto const& params :
       testing::GetBfvParametersForFiniteFieldEncoding<TypeParam>()) {
    this->SetUpBfvParameters(params);

    ASSERT_OK_AND_ASSIGN(auto encoder,
                         Encoder::Create(this->rns_context_.get()));

    ASSERT_OK_AND_ASSIGN(RnsRlweSecretKey<TypeParam> secret_key,
                         this->SampleSecretKey());

    auto gadget = this->gadget_.get();
    int k = gadget->Dimension();
    int log_n = this->rns_context_->LogN();
    ASSERT_OK_AND_ASSIGN(
        std::vector<RnsPolynomial<TypeParam>> key_as,
        RnsGaloisKey<TypeParam>::SampleRandomPad(k, log_n, this->main_moduli_,
                                                 kPrngSeed, kPrngType));

    Integer t = this->rns_context_->PlaintextModulus();
    ASSERT_OK_AND_ASSIGN(RnsGaloisKey<TypeParam> gk,
                         RnsGaloisKey<TypeParam>::CreateWithRandomPadForBfv(
                             std::move(key_as), secret_key, k_power, kVariance,
                             gadget, kPrngSeed, kPrngType));
    ASSERT_EQ(gk.Gadget(), gadget);
    ASSERT_EQ(gk.Dimension(), gadget->Dimension());
    ASSERT_EQ(gk.SubstitutionPower(), k_power);

    int num_coeffs = 1 << log_n;
    std::vector<Integer> values =
        testing::SampleMessages<Integer>(num_coeffs, t);
    ASSERT_OK_AND_ASSIGN(
        RnsBfvCiphertext<TypeParam> ciphertext,
        secret_key.template EncryptBfv<Encoder>(
            values, &encoder, this->error_params_.get(), this->prng_.get()));

    // Precompute the gadget decomposition of the "a" component of `ciphertext`,
    // and also precompute the "a" component of the key-switched ciphertext.
    ASSERT_OK_AND_ASSIGN(RnsPolynomial<TypeParam> c1, ciphertext.Component(1));
    ASSERT_OK_AND_ASSIGN(RnsPolynomial<TypeParam> c1_sub,
                         c1.Substitute(k_power, this->main_moduli_));
    ASSERT_OK(c1_sub.ConvertToCoeffForm(this->main_moduli_));
    ASSERT_OK_AND_ASSIGN(std::vector<RnsPolynomial<TypeParam>> c1_sub_digits,
                         gadget->Decompose(c1_sub, this->main_moduli_));
    for (auto& c1_sub_digit : c1_sub_digits) {
      ASSERT_OK(c1_sub_digit.ConvertToNttForm(this->main_moduli_));
    }
    ASSERT_OK_AND_ASSIGN(RnsPolynomial<TypeParam> c1_rot,
                         gk.PrecomputeRandomPad(c1_sub_digits));

    // Apply the automorphism X -> X^k_power to the ciphertext with precomputed
    // "a" components of the resulting ciphertext.
    ASSERT_OK_AND_ASSIGN(auto ciphertext_sub, ciphertext.Substitute(k_power));
    ASSERT_OK_AND_ASSIGN(auto ciphertext_rot,
                         gk.ApplyToWithRandomPad(ciphertext_sub, c1_sub_digits,
                                                 std::move(c1_rot)));

    // Decrypt and check results of homomorphic rotation.
    ASSERT_OK_AND_ASSIGN(
        auto decryptions,
        secret_key.template DecryptBfv<Encoder>(ciphertext_rot, &encoder));
    int group_size = num_coeffs / 2;
    std::vector<Integer> expected;
    for (int i = 0; i < group_size; ++i) {
      int index = (i + 1) % group_size;
      expected.push_back(values[index]);
    }
    for (int i = 0; i < group_size; ++i) {
      int index = group_size + (i + 1) % group_size;
      expected.push_back(values[index]);
    }

    EXPECT_EQ(decryptions, expected);
  }
}

// Check that a Galois key composed from two Galois keys with substitution power
// a and b, resp, can key switch a ciphertext with substitution power a*b % 2N.
TYPED_TEST(RnsGaloisKeyTest, ComposedGaloisKeyForBgv) {
  using Encoder = FiniteFieldEncoder<TypeParam>;
  using Integer = typename TypeParam::Int;
  constexpr int k_log_gadget_base = 5;
  constexpr int k_power0 = 5;   // rotation by one slot.
  constexpr int k_power1 = 25;  // rotation by two slots.
  constexpr int k_combined_offset = 3;
  for (auto params :
       testing::GetRnsParametersForFiniteFieldEncoding<TypeParam>()) {
    // Use a smaller gadget base to keep the composed Galois key's error small.
    params.log_gadget_base = k_log_gadget_base;
    this->SetUpBgvParameters(params);

    ASSERT_OK_AND_ASSIGN(auto encoder,
                         Encoder::Create(this->rns_context_.get()));
    ASSERT_OK_AND_ASSIGN(RnsRlweSecretKey<TypeParam> secret_key,
                         this->SampleSecretKey());

    // Create two Galois keys with different substitution powers.
    Integer t = this->rns_context_->PlaintextModulus();
    ASSERT_OK_AND_ASSIGN(RnsGaloisKey<TypeParam> gk0,
                         RnsGaloisKey<TypeParam>::CreateForBgv(
                             secret_key, k_power0, kVariance,
                             this->gadget_.get(), t, kPrngType));
    ASSERT_OK_AND_ASSIGN(RnsGaloisKey<TypeParam> gk1,
                         RnsGaloisKey<TypeParam>::CreateForBgv(
                             secret_key, k_power1, kVariance,
                             this->gadget_.get(), t, kPrngType));

    // Compose the two Galois keys.
    ASSERT_OK_AND_ASSIGN(RnsGaloisKey<TypeParam> gk_composed, gk0.Compose(gk1));
    ASSERT_EQ(gk_composed.SubstitutionPower(), k_power0 * k_power1);

    // Apply the composed Galois key on a ciphertext w/ substitution power a*b.
    int num_coeffs = 1 << this->rns_context_->LogN();
    std::vector<Integer> values =
        testing::SampleMessages<Integer>(num_coeffs, t);
    ASSERT_OK_AND_ASSIGN(
        RnsBgvCiphertext<TypeParam> ciphertext,
        secret_key.template EncryptBgv<Encoder>(
            values, &encoder, this->error_params_.get(), this->prng_.get()));
    ASSERT_OK_AND_ASSIGN(auto ciphertext_sub,
                         ciphertext.Substitute(k_power0 * k_power1));
    ASSERT_OK_AND_ASSIGN(auto ciphertext_rot,
                         gk_composed.ApplyTo(ciphertext_sub));
    ASSERT_EQ(ciphertext_rot.PowerOfS(), 1);

    // The key-switched ciphertext should encrypt a slot vector that is rotated
    // by the combined offset.
    ASSERT_OK_AND_ASSIGN(
        auto decryptions,
        secret_key.template DecryptBgv<Encoder>(ciphertext_rot, &encoder));

    int group_size = num_coeffs / 2;
    std::vector<Integer> expected;
    for (int i = 0; i < group_size; ++i) {
      int index = (i + k_combined_offset) % group_size;
      expected.push_back(values[index]);
    }
    for (int i = 0; i < group_size; ++i) {
      int index = group_size + (i + k_combined_offset) % group_size;
      expected.push_back(values[index]);
    }

    EXPECT_EQ(decryptions, expected);
  }
}

// Similar to the previous test, but generate a composed Galois key using
// precomputed "a" (aka pad) components and use the Galois key for BFV.
TYPED_TEST(RnsGaloisKeyTest, ComposedGaloisKeyForBfvWithPrecomputedPads) {
  using Encoder = FiniteFieldEncoder<TypeParam>;
  using Integer = typename TypeParam::Int;
  constexpr int k_log_gadget_base = 5;
  constexpr int k_power0 = 5;   // rotation by one slot.
  constexpr int k_power1 = 25;  // rotation by two slots.
  constexpr int k_combined_offset = 3;
  for (auto params :
       testing::GetBfvParametersForFiniteFieldEncoding<TypeParam>()) {
    // This is a very rough estimation on whether the parameter is suitable for
    // Galois key composition.
    if (params.log_gadget_base < k_log_gadget_base * 2) {
      continue;
    }

    // Use a smaller gadget base to keep the composed Galois key's error small.
    params.log_gadget_base = k_log_gadget_base;
    this->SetUpBfvParameters(params);

    ASSERT_OK_AND_ASSIGN(auto encoder,
                         Encoder::Create(this->rns_context_.get()));
    ASSERT_OK_AND_ASSIGN(RnsRlweSecretKey<TypeParam> secret_key,
                         this->SampleSecretKey());

    // Create two Galois keys with different substitution powers.
    Integer t = this->rns_context_->PlaintextModulus();
    ASSERT_OK_AND_ASSIGN(
        RnsGaloisKey<TypeParam> gk0,
        RnsGaloisKey<TypeParam>::CreateForBfv(secret_key, k_power0, kVariance,
                                              this->gadget_.get(), kPrngType));
    ASSERT_OK_AND_ASSIGN(
        RnsGaloisKey<TypeParam> gk1,
        RnsGaloisKey<TypeParam>::CreateForBfv(secret_key, k_power1, kVariance,
                                              this->gadget_.get(), kPrngType));

    // Precompute the "a" components of the composed Galois key.
    int gk_dimension = gk1.Dimension();
    const std::vector<RnsPolynomial<TypeParam>>& gk1_as = gk1.GetKeyA();
    std::vector<RnsPolynomial<TypeParam>> gk1_sub_as;
    gk1_sub_as.reserve(gk_dimension);
    for (int i = 0; i < gk_dimension; ++i) {
      ASSERT_OK_AND_ASSIGN(auto a_sub,
                           gk1_as[i].Substitute(k_power0, this->main_moduli_));
      ASSERT_OK(a_sub.ConvertToCoeffForm(this->main_moduli_));
      gk1_sub_as.push_back(std::move(a_sub));
    }
    ASSERT_OK_AND_ASSIGN(
        auto gk1_a_sub_digits,
        this->gadget_->Decompose(gk1_sub_as, this->main_moduli_));
    for (auto& digits : gk1_a_sub_digits) {
      for (auto& digit : digits) {
        ASSERT_OK(digit.ConvertToNttForm(this->main_moduli_));
      }
    }
    ASSERT_OK_AND_ASSIGN(
        auto gk_composed_as,
        RnsGaloisKey<TypeParam>::CreatePadOfComposedKey(
            gk0.GetKeyA(), gk1_a_sub_digits, this->main_moduli_));

    // Generate the composed key using the precompued "a" components.
    ASSERT_OK_AND_ASSIGN(
        RnsGaloisKey<TypeParam> gk_composed,
        gk0.ComposeWithPads(gk1, gk1_a_sub_digits, gk_composed_as));
    ASSERT_EQ(gk_composed.SubstitutionPower(), k_power0 * k_power1);

    // Apply the composed Galois key on a ciphertext w/ substitution power a*b.
    int num_coeffs = 1 << this->rns_context_->LogN();
    std::vector<Integer> values =
        testing::SampleMessages<Integer>(num_coeffs, t);
    ASSERT_OK_AND_ASSIGN(
        auto ciphertext,
        secret_key.template EncryptBfv<Encoder>(
            values, &encoder, this->error_params_.get(), this->prng_.get()));
    ASSERT_OK_AND_ASSIGN(auto ciphertext_sub,
                         ciphertext.Substitute(k_power0 * k_power1));
    ASSERT_OK_AND_ASSIGN(auto ciphertext_rot,
                         gk_composed.ApplyTo(ciphertext_sub));
    ASSERT_EQ(ciphertext_rot.PowerOfS(), 1);

    // The key-switched ciphertext should encrypt a slot vector that is rotated
    // by the combined offset.
    ASSERT_OK_AND_ASSIGN(
        auto decryptions,
        secret_key.template DecryptBfv<Encoder>(ciphertext_rot, &encoder));

    int group_size = num_coeffs / 2;
    std::vector<Integer> expected;
    for (int i = 0; i < group_size; ++i) {
      int index = (i + k_combined_offset) % group_size;
      expected.push_back(values[index]);
    }
    for (int i = 0; i < group_size; ++i) {
      int index = group_size + (i + k_combined_offset) % group_size;
      expected.push_back(values[index]);
    }
    EXPECT_EQ(decryptions, expected);
  }
}

TYPED_TEST(RnsGaloisKeyTest, SerializedDeserializes) {
  constexpr int k_power = 5;

  auto params = testing::GetRnsParametersForFiniteFieldEncoding<TypeParam>()[0];
  this->SetUpBgvParameters(params);

  ASSERT_OK_AND_ASSIGN(RnsRlweSecretKey<TypeParam> secret_key,
                       this->SampleSecretKey());
  auto gadget = this->gadget_.get();
  auto t = this->rns_context_->PlaintextModulus();
  ASSERT_OK_AND_ASSIGN(
      RnsGaloisKey<TypeParam> gk,
      RnsGaloisKey<TypeParam>::CreateForBgv(secret_key, k_power, kVariance,
                                            gadget, t, kPrngType));

  ASSERT_OK_AND_ASSIGN(SerializedRnsGaloisKey serialized, gk.Serialize());
  ASSERT_OK_AND_ASSIGN(RnsGaloisKey<TypeParam> deserialized,
                       RnsGaloisKey<TypeParam>::Deserialize(
                           serialized, gadget, this->main_moduli_));
  ASSERT_EQ(deserialized.Dimension(), gk.Dimension());
  ASSERT_EQ(deserialized.SubstitutionPower(), gk.SubstitutionPower());
  ASSERT_EQ(deserialized.GetKeyA(), gk.GetKeyA());
  ASSERT_EQ(deserialized.GetKeyB(), gk.GetKeyB());
}

}  // namespace
}  // namespace rlwe
