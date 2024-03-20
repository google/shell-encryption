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

#include "shell_encryption/rns/rns_secret_key.h"

#include <cmath>
#include <memory>
#include <utility>
#include <vector>

#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include "absl/log/check.h"
#include "absl/status/status.h"
#include "absl/strings/string_view.h"
#include "absl/types/span.h"
#include "shell_encryption/prng/single_thread_hkdf_prng.h"
#include "shell_encryption/rns/coefficient_encoder.h"
#include "shell_encryption/rns/finite_field_encoder.h"
#include "shell_encryption/rns/rns_bfv_ciphertext.h"
#include "shell_encryption/rns/rns_bgv_ciphertext.h"
#include "shell_encryption/rns/rns_context.h"
#include "shell_encryption/rns/rns_error_params.h"
#include "shell_encryption/rns/rns_modulus.h"
#include "shell_encryption/rns/rns_polynomial.h"
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

// Test fixture to take care of messy setup.
template <typename ModularInt>
class RnsRlweSecretKeyTest : public ::testing::Test {
 protected:
  void SetUp() override {
    // Create a PRNG for sampling secret key.
    auto prng = Prng::Create(kPrngSeed.substr(0, Prng::SeedLength()));
    CHECK(prng.ok());
    prng_ = std::move(prng.value());
  }

  void SetUpBgvContext() {
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
  }

  void SetUpBfvContext() {
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

TYPED_TEST_SUITE(RnsRlweSecretKeyTest, testing::ModularIntTypes);

TYPED_TEST(RnsRlweSecretKeyTest, SampleFailsIfLogNIsNotPositive) {
  this->SetUpBgvContext();
  EXPECT_THAT(
      RnsRlweSecretKey<TypeParam>::Sample(
          /*log_n=*/-1, kVariance, this->main_moduli_, this->prng_.get()),
      StatusIs(absl::StatusCode::kInvalidArgument,
               HasSubstr("`log_n` must be positive")));
  EXPECT_THAT(
      RnsRlweSecretKey<TypeParam>::Sample(
          /*log_n=*/0, kVariance, this->main_moduli_, this->prng_.get()),
      StatusIs(absl::StatusCode::kInvalidArgument,
               HasSubstr("`log_n` must be positive")));
}

TYPED_TEST(RnsRlweSecretKeyTest, SampleFailsIfVarianceIsNotPositive) {
  this->SetUpBgvContext();
  EXPECT_THAT(RnsRlweSecretKey<TypeParam>::Sample(
                  this->rns_context_->LogN(), /*variance=*/-1,
                  this->main_moduli_, this->prng_.get()),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`variance` must be positive")));
  EXPECT_THAT(RnsRlweSecretKey<TypeParam>::Sample(
                  this->rns_context_->LogN(),
                  /*variance=*/0, this->main_moduli_, this->prng_.get()),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`variance` must be positive")));
}

TYPED_TEST(RnsRlweSecretKeyTest, SampleFailsIfModuliIsEmpty) {
  this->SetUpBgvContext();
  EXPECT_THAT(
      RnsRlweSecretKey<TypeParam>::Sample(this->rns_context_->LogN(), kVariance,
                                          /*moduli=*/{}, nullptr),
      StatusIs(absl::StatusCode::kInvalidArgument,
               HasSubstr("`moduli` must contain at least one element")));
}

TYPED_TEST(RnsRlweSecretKeyTest, SampleFailsIfPrngIsNull) {
  this->SetUpBgvContext();
  EXPECT_THAT(RnsRlweSecretKey<TypeParam>::Sample(
                  this->rns_context_->LogN(), kVariance,
                  this->rns_context_->MainPrimeModuli(), /*prng=*/nullptr),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`prng` must not be null")));
}

TYPED_TEST(RnsRlweSecretKeyTest, SampleSucceeds) {
  this->SetUpBgvContext();
  ASSERT_OK_AND_ASSIGN(
      RnsRlweSecretKey<TypeParam> key,
      RnsRlweSecretKey<TypeParam>::Sample(this->rns_context_->LogN(), kVariance,
                                          this->rns_context_->MainPrimeModuli(),
                                          this->prng_.get()));
  EXPECT_EQ(key.LogN(), this->rns_context_->LogN());
  EXPECT_EQ(key.Level() + 1, this->rns_context_->NumMainPrimeModuli());
  EXPECT_EQ(key.NumCoeffs(), 1 << this->rns_context_->LogN());
  EXPECT_EQ(key.NumModuli(), this->rns_context_->NumMainPrimeModuli());

  // The secret key is sampled wrt the main moduli in the RNS context.
  absl::Span<const PrimeModulus<TypeParam>* const> key_moduli = key.Moduli();
  EXPECT_EQ(key_moduli.size(), this->main_moduli_.size());
  for (int i = 0; i < key_moduli.size(); ++i) {
    EXPECT_EQ(key_moduli[i]->Modulus(), this->main_moduli_[i]->Modulus());
  }
  EXPECT_EQ(key.Variance(), kVariance);
}

TYPED_TEST(RnsRlweSecretKeyTest, KeyCoefficientsAreBoundedAndConsistent) {
  this->SetUpBgvContext();
  ASSERT_OK_AND_ASSIGN(RnsRlweSecretKey<TypeParam> key, this->SampleKey());
  // Get a copy of key polynomial s, where the canonical secret key is (1, s).
  RnsPolynomial<TypeParam> s = key.Key();
  // Convert the key polynomial s to coefficient form.
  ASSERT_OK(s.ConvertToCoeffForm(this->main_moduli_));

  // The coefficients of the polynomial s (when lifted to Z[X]/(X^N+1)) should
  // have bounded absolute values (<= 2*variance). Furthermore, since variance
  // is smaller than all prime moduli q_i, we should have s_coeffs[h] (mod
  // q_i) and s_coeffs[h] (mod q_j) represent the same integer value for all
  // q_i, q_j
  const std::vector<std::vector<TypeParam>>& s_coeffs = s.Coeffs();
  EXPECT_EQ(s_coeffs.size(), this->main_moduli_.size());
  std::vector<typename TypeParam::Int> s_abs_coeffs(s.NumCoeffs());
  std::vector<bool> s_sign_coeffs(s.NumCoeffs());  // is s_coeffs negative?
  for (int i = 0; i < s_coeffs.size(); ++i) {
    auto s_coeffs_qi = s_coeffs[i];
    auto mod_params_qi = this->main_moduli_[i]->ModParams();
    typename TypeParam::Int qi = mod_params_qi->modulus;
    typename TypeParam::Int qi_half = qi / 2;
    for (int h = 0; h < s_coeffs_qi.size(); ++h) {
      typename TypeParam::Int s_coeff_mod_qi =
          s_coeffs_qi[h].ExportInt(mod_params_qi);
      bool is_coeff_negative = s_coeff_mod_qi > qi_half;
      typename TypeParam::Int s_abs_coeff =
          is_coeff_negative ? qi - s_coeff_mod_qi : s_coeff_mod_qi;
      EXPECT_LE(s_abs_coeff, 2 * kVariance);
      if (i == 0) {
        s_abs_coeffs[h] = s_abs_coeff;
        s_sign_coeffs[h] = is_coeff_negative;
      } else {
        EXPECT_EQ(s_abs_coeffs[h], s_abs_coeff);
        EXPECT_EQ(s_sign_coeffs[h], is_coeff_negative);
      }
    }
  }
}

TYPED_TEST(RnsRlweSecretKeyTest, EncryptPolynomialBgvFailsIfEncoderIsNull) {
  this->SetUpBgvContext();
  ASSERT_OK_AND_ASSIGN(RnsRlweSecretKey<TypeParam> key, this->SampleKey());

  ASSERT_OK_AND_ASSIGN(auto zero, RnsPolynomial<TypeParam>::CreateZero(
                                      this->rns_context_->LogN(),
                                      this->main_moduli_, /*is_ntt=*/true));
  EXPECT_THAT(key.template EncryptPolynomialBgv<CoefficientEncoder<TypeParam>>(
                  /*plaintext=*/zero, /*encoder=*/nullptr,
                  this->error_params_.get(), this->prng_.get()),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`encoder` must not be null")));
}

TYPED_TEST(RnsRlweSecretKeyTest,
           EncryptPolynomialBgvFailsIfPlaintextCoeffForm) {
  this->SetUpBgvContext();
  ASSERT_OK_AND_ASSIGN(RnsRlweSecretKey<TypeParam> key, this->SampleKey());

  ASSERT_OK_AND_ASSIGN(auto zero, RnsPolynomial<TypeParam>::CreateZero(
                                      this->rns_context_->LogN(),
                                      this->main_moduli_, /*is_ntt=*/false));
  ASSERT_EQ(zero.IsNttForm(), false);

  EXPECT_THAT(key.template EncryptPolynomialBgv<CoefficientEncoder<TypeParam>>(
                  /*plaintext=*/zero, /*encoder=*/this->coeff_encoder_.get(),
                  this->error_params_.get(), this->prng_.get()),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`plaintext` must be in NTT form")));
}

TYPED_TEST(RnsRlweSecretKeyTest, EncryptBgvFailsIfEncoderIsNull) {
  this->SetUpBgvContext();
  ASSERT_OK_AND_ASSIGN(RnsRlweSecretKey<TypeParam> key, this->SampleKey());
  EXPECT_THAT(key.template EncryptBgv<CoefficientEncoder<TypeParam>>(
                  /*messages=*/{}, /*encoder=*/nullptr,
                  this->error_params_.get(), this->prng_.get()),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`encoder` must not be null")));
}

TYPED_TEST(RnsRlweSecretKeyTest, EncryptBgvFailsIfErrorParamsIsNull) {
  this->SetUpBgvContext();
  ASSERT_OK_AND_ASSIGN(RnsRlweSecretKey<TypeParam> key, this->SampleKey());
  EXPECT_THAT(key.EncryptBgv(/*messages=*/{}, this->coeff_encoder_.get(),
                             /*error_params=*/nullptr, this->prng_.get()),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`error_params` must not be null")));
}

TYPED_TEST(RnsRlweSecretKeyTest, EncryptBgvFailsIfPrngIsNull) {
  this->SetUpBgvContext();
  ASSERT_OK_AND_ASSIGN(RnsRlweSecretKey<TypeParam> key, this->SampleKey());
  EXPECT_THAT(key.EncryptBgv(/*messages=*/{}, this->coeff_encoder_.get(),
                             this->error_params_.get(), /*prng=*/nullptr),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`prng` must not be null")));
}

TYPED_TEST(RnsRlweSecretKeyTest, DecryptBgvFailsIfCiphertextPowerOfSIsNotOne) {
  this->SetUpBgvContext();
  ASSERT_OK_AND_ASSIGN(RnsRlweSecretKey<TypeParam> key, this->SampleKey());

  // Create a ciphertext with power_of_s == 1.
  ASSERT_OK_AND_ASSIGN(auto zero,
                       RnsPolynomial<TypeParam>::CreateZero(
                           this->rns_context_->LogN(), this->main_moduli_));
  RnsBgvCiphertext<TypeParam> ciphertext(
      /*components=*/{zero, zero}, this->main_moduli_,
      /*power_of_s=*/5, /*error=*/0, this->error_params_.get());
  ASSERT_EQ(ciphertext.PowerOfS(), 5);

  EXPECT_THAT(
      key.DecryptBgv(ciphertext, this->coeff_encoder_.get()),
      StatusIs(
          absl::StatusCode::kInvalidArgument,
          HasSubstr(
              "Cannot decrypt `ciphertext` with power of s not equal to 1.")));
}

TYPED_TEST(RnsRlweSecretKeyTest,
           DecryptBgvFailsIfCiphertextHasMismatchPolyDegree) {
  this->SetUpBgvContext();
  ASSERT_OK_AND_ASSIGN(RnsRlweSecretKey<TypeParam> key, this->SampleKey());

  // Create a ciphertext with log_n that's different from the secret key.
  int log_n = this->rns_context_->LogN();
  ASSERT_OK_AND_ASSIGN(auto zero, RnsPolynomial<TypeParam>::CreateZero(
                                      log_n - 1, this->main_moduli_));
  RnsBgvCiphertext<TypeParam> ciphertext(
      /*components=*/{zero, zero}, this->main_moduli_,
      /*power_of_s=*/1, /*error=*/0, this->error_params_.get());
  ASSERT_EQ(ciphertext.LogN(), log_n - 1);

  EXPECT_THAT(key.DecryptBgv(ciphertext, this->coeff_encoder_.get()),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("Cannot decrypt `ciphertext` with a "
                                 "mismatching polynomial degree.")));
}

TYPED_TEST(RnsRlweSecretKeyTest,
           DecryptBgvFailsIfCiphertextHasMismatchNumberOfModuli) {
  this->SetUpBgvContext();
  // This test requires at least two RNS moduli.
  if (this->main_moduli_.size() <= 1) {
    return;
  }

  ASSERT_OK_AND_ASSIGN(RnsRlweSecretKey<TypeParam> key, this->SampleKey());

  // Create a ciphertext with fewer number of RNS moduli.
  int log_n = this->rns_context_->LogN();
  int level = this->main_moduli_.size() - 1;
  auto moduli_reduced = this->main_moduli_;
  moduli_reduced.pop_back();
  ASSERT_OK_AND_ASSIGN(
      auto zero, RnsPolynomial<TypeParam>::CreateZero(log_n, moduli_reduced));
  RnsBgvCiphertext<TypeParam> ciphertext(
      /*components=*/{zero, zero}, moduli_reduced,
      /*power_of_s=*/1, /*error=*/0, this->error_params_.get());
  ASSERT_EQ(ciphertext.NumModuli(), level);

  EXPECT_THAT(key.DecryptBgv(ciphertext, this->coeff_encoder_.get()),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("Cannot decrypt `ciphertext` with a "
                                 "mismatching number of prime moduli.")));
}

TYPED_TEST(RnsRlweSecretKeyTest, EncryptBgvDecrypts) {
  using Integer = typename TypeParam::Int;
  this->SetUpBgvContext();
  ASSERT_OK_AND_ASSIGN(RnsRlweSecretKey<TypeParam> key, this->SampleKey());
  int log_n = this->rns_context_->LogN();
  Integer t = this->rns_context_->PlaintextModulus();
  std::vector<Integer> messages = testing::SampleMessages(1 << log_n, t);
  ASSERT_OK_AND_ASSIGN(
      RnsBgvCiphertext<TypeParam> ciphertext,
      key.EncryptBgv(messages, this->coeff_encoder_.get(),
                     this->error_params_.get(), this->prng_.get()));

  ASSERT_OK_AND_ASSIGN(std::vector<Integer> dec_messages,
                       key.DecryptBgv(ciphertext, this->coeff_encoder_.get()));
  EXPECT_EQ(dec_messages, messages);
}

TYPED_TEST(RnsRlweSecretKeyTest, DecryptBgvOnDegreeTwoCiphertext) {
  using Integer = typename TypeParam::Int;
  this->SetUpBgvContext();
  double expected_error_bound = this->error_params_->B_secretkey_encryption() *
                                this->error_params_->B_secretkey_encryption();
  if (expected_error_bound >=
      static_cast<double>(testing::CompositeModulus<Integer>::Value())) {
    // The test parameters are not suitable for homomorphic multiplication as
    // the error in the product ciphertext is expected to be larger than the
    // ciphertext modulus. This should only apply to the test with Uint16 as
    // the underlying integer type.
    return;
  }

  ASSERT_OK_AND_ASSIGN(RnsRlweSecretKey<TypeParam> key, this->SampleKey());

  // Sample two message vectors and encrypt them.
  const CoefficientEncoder<TypeParam>* encoder = this->coeff_encoder_.get();
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

  // Multiply the two ciphertexts together to get a degree 2 ciphertext.
  ASSERT_OK_AND_ASSIGN(RnsBgvCiphertext<TypeParam> ciphertext_product,
                       ciphertext0 * ciphertext1);
  EXPECT_EQ(ciphertext_product.Degree(), 2);
  EXPECT_EQ(ciphertext_product.Level(), ciphertext0.Level());

  // Decrypt the product ciphertext.
  ASSERT_OK_AND_ASSIGN(std::vector<Integer> dec_messages,
                       key.DecryptBgv(ciphertext_product, encoder));

  // Encode the message vectors and multiply them over the plaintext space.
  ASSERT_OK_AND_ASSIGN(RnsPolynomial<TypeParam> plaintext0,
                       encoder->EncodeBgv(messages0, this->main_moduli_));
  ASSERT_OK_AND_ASSIGN(RnsPolynomial<TypeParam> plaintext1,
                       encoder->EncodeBgv(messages1, this->main_moduli_));
  ASSERT_OK(plaintext0.ConvertToNttForm(this->main_moduli_));
  ASSERT_OK(plaintext1.ConvertToNttForm(this->main_moduli_));
  ASSERT_OK_AND_ASSIGN(RnsPolynomial<TypeParam> expected_polynomial,
                       plaintext0.Mul(plaintext1, this->main_moduli_));
  ASSERT_OK_AND_ASSIGN(
      std::vector<Integer> expected_messages,
      encoder->DecodeBgv(expected_polynomial, this->main_moduli_));
  EXPECT_EQ(dec_messages, expected_messages);
}

TYPED_TEST(RnsRlweSecretKeyTest, EncryptBfvFailsIfEncoderIsNull) {
  this->SetUpBfvContext();
  ASSERT_OK_AND_ASSIGN(RnsRlweSecretKey<TypeParam> key, this->SampleKey());
  EXPECT_THAT(key.template EncryptBfv<CoefficientEncoder<TypeParam>>(
                  /*messages=*/{}, /*encoder=*/nullptr,
                  this->error_params_.get(), this->prng_.get()),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`encoder` must not be null")));
}

TYPED_TEST(RnsRlweSecretKeyTest, EncryptBfvFailsIfErrorParamsIsNull) {
  this->SetUpBfvContext();
  ASSERT_OK_AND_ASSIGN(RnsRlweSecretKey<TypeParam> key, this->SampleKey());
  EXPECT_THAT(key.EncryptBfv(/*messages=*/{}, this->coeff_encoder_.get(),
                             /*error_params=*/nullptr, this->prng_.get()),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`error_params` must not be null")));
}

TYPED_TEST(RnsRlweSecretKeyTest, EncryptBfvFailsIfPrngIsNull) {
  this->SetUpBfvContext();
  ASSERT_OK_AND_ASSIGN(RnsRlweSecretKey<TypeParam> key, this->SampleKey());
  EXPECT_THAT(key.EncryptBfv(/*messages=*/{}, this->coeff_encoder_.get(),
                             this->error_params_.get(), /*prng=*/nullptr),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`prng` must not be null")));
}

TYPED_TEST(RnsRlweSecretKeyTest, EncryptBfvFailsIfPrngPadIsNull) {
  this->SetUpBfvContext();
  ASSERT_OK_AND_ASSIGN(RnsRlweSecretKey<TypeParam> key, this->SampleKey());
  EXPECT_THAT(key.EncryptBfv(/*messages=*/{}, this->coeff_encoder_.get(),
                             this->error_params_.get(),
                             /*prng=*/this->prng_.get(), /*prng_pad=*/nullptr),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`prng_pad` must not be null")));
}

TYPED_TEST(RnsRlweSecretKeyTest, DecryptBfvFailsIfCiphertextPowerOfSIsNotOne) {
  this->SetUpBfvContext();
  ASSERT_OK_AND_ASSIGN(RnsRlweSecretKey<TypeParam> key, this->SampleKey());

  // Create a ciphertext with power_of_s == 1.
  ASSERT_OK_AND_ASSIGN(auto zero,
                       RnsPolynomial<TypeParam>::CreateZero(
                           this->rns_context_->LogN(), this->main_moduli_));
  RnsBfvCiphertext<TypeParam> ciphertext(
      /*components=*/{zero, zero}, this->main_moduli_,
      /*power_of_s=*/5, /*error=*/0, this->error_params_.get(),
      this->rns_context_.get());
  ASSERT_EQ(ciphertext.PowerOfS(), 5);

  EXPECT_THAT(
      key.DecryptBfv(ciphertext, this->coeff_encoder_.get()),
      StatusIs(
          absl::StatusCode::kInvalidArgument,
          HasSubstr(
              "Cannot decrypt `ciphertext` with power of s not equal to 1.")));
}

TYPED_TEST(RnsRlweSecretKeyTest,
           DecryptBfvFailsIfCiphertextHasMismatchPolyDegree) {
  this->SetUpBfvContext();
  ASSERT_OK_AND_ASSIGN(RnsRlweSecretKey<TypeParam> key, this->SampleKey());

  // Create a ciphertext with log_n that's different from the secret key.
  int log_n = this->rns_context_->LogN();
  ASSERT_OK_AND_ASSIGN(auto zero, RnsPolynomial<TypeParam>::CreateZero(
                                      log_n - 1, this->main_moduli_));
  RnsBfvCiphertext<TypeParam> ciphertext(
      /*components=*/{zero, zero}, this->main_moduli_,
      /*power_of_s=*/1, /*error=*/0, this->error_params_.get(),
      this->rns_context_.get());
  ASSERT_EQ(ciphertext.LogN(), log_n - 1);

  EXPECT_THAT(key.DecryptBfv(ciphertext, this->coeff_encoder_.get()),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("Cannot decrypt `ciphertext` with a "
                                 "mismatching polynomial degree.")));
}

TYPED_TEST(RnsRlweSecretKeyTest,
           DecryptBfvFailsIfCiphertextHasMismatchNumberOfModuli) {
  this->SetUpBfvContext();
  // This test requires at least two RNS moduli.
  if (this->main_moduli_.size() <= 1) {
    return;
  }

  ASSERT_OK_AND_ASSIGN(RnsRlweSecretKey<TypeParam> key, this->SampleKey());

  // Create a ciphertext with fewer number of RNS moduli.
  int log_n = this->rns_context_->LogN();
  int level = this->main_moduli_.size() - 1;
  auto moduli_reduced = this->main_moduli_;
  moduli_reduced.pop_back();
  ASSERT_OK_AND_ASSIGN(
      auto zero, RnsPolynomial<TypeParam>::CreateZero(log_n, moduli_reduced));
  RnsBfvCiphertext<TypeParam> ciphertext(
      /*components=*/{zero, zero}, moduli_reduced,
      /*power_of_s=*/1, /*error=*/0, this->error_params_.get(),
      this->rns_context_.get());
  ASSERT_EQ(ciphertext.NumModuli(), level);

  EXPECT_THAT(key.DecryptBfv(ciphertext, this->coeff_encoder_.get()),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("Cannot decrypt `ciphertext` with a "
                                 "mismatching number of prime moduli.")));
}

TYPED_TEST(RnsRlweSecretKeyTest, DecryptBfvFailsIfEncoderIsNull) {
  using Encoder = CoefficientEncoder<TypeParam>;
  this->SetUpBfvContext();
  ASSERT_OK_AND_ASSIGN(RnsRlweSecretKey<TypeParam> key, this->SampleKey());

  ASSERT_OK_AND_ASSIGN(auto zero,
                       RnsPolynomial<TypeParam>::CreateZero(
                           this->rns_context_->LogN(), this->main_moduli_));
  RnsBfvCiphertext<TypeParam> ciphertext(
      /*components=*/{zero, zero}, this->main_moduli_,
      /*power_of_s=*/1, /*error=*/0, this->error_params_.get(),
      this->rns_context_.get());

  EXPECT_THAT(key.template DecryptBfv<Encoder>(ciphertext, /*encoder=*/nullptr),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`encoder` must not be null")));
}

TYPED_TEST(RnsRlweSecretKeyTest, EncryptBfvDecrypts) {
  using Integer = typename TypeParam::Int;
  this->SetUpBfvContext();
  ASSERT_OK_AND_ASSIGN(RnsRlweSecretKey<TypeParam> key, this->SampleKey());
  int log_n = this->rns_context_->LogN();
  Integer t = this->rns_context_->PlaintextModulus();
  std::vector<Integer> messages = testing::SampleMessages(1 << log_n, t);

  ASSERT_OK_AND_ASSIGN(
      RnsBfvCiphertext<TypeParam> ciphertext,
      key.EncryptBfv(messages, this->coeff_encoder_.get(),
                     this->error_params_.get(), this->prng_.get()));

  ASSERT_OK_AND_ASSIGN(std::vector<Integer> dec_messages,
                       key.DecryptBfv(ciphertext, this->coeff_encoder_.get()));
  EXPECT_EQ(dec_messages, messages);
}

TYPED_TEST(RnsRlweSecretKeyTest, EncryptBfvWithGivenRandomPadDecrypts) {
  using Integer = typename TypeParam::Int;
  this->SetUpBfvContext();
  ASSERT_OK_AND_ASSIGN(RnsRlweSecretKey<TypeParam> key, this->SampleKey());
  int log_n = this->rns_context_->LogN();
  Integer t = this->rns_context_->PlaintextModulus();
  std::vector<Integer> messages = testing::SampleMessages(1 << log_n, t);

  // Create a PRNG for sampling the random "a" polynomial in the ciphertext.
  ASSERT_OK_AND_ASSIGN(auto prng_seed_pad, Prng::GenerateSeed());
  ASSERT_OK_AND_ASSIGN(auto prng_pad, Prng::Create(prng_seed_pad));

  ASSERT_OK_AND_ASSIGN(RnsBfvCiphertext<TypeParam> ciphertext,
                       key.EncryptBfv(messages, this->coeff_encoder_.get(),
                                      this->error_params_.get(),
                                      this->prng_.get(), prng_pad.get()));

  // Check that the "a" component in `ciphertext` is generated as expected.
  ASSERT_OK_AND_ASSIGN(auto prng_pad_check, Prng::Create(prng_seed_pad));
  ASSERT_OK_AND_ASSIGN(auto expected_pad,
                       RnsPolynomial<TypeParam>::SampleUniform(
                           log_n, prng_pad_check.get(), this->main_moduli_));
  ASSERT_OK(expected_pad.NegateInPlace(this->main_moduli_));
  ASSERT_OK_AND_ASSIGN(auto actual_pad, ciphertext.Component(1));
  EXPECT_EQ(actual_pad, expected_pad);

  ASSERT_OK_AND_ASSIGN(std::vector<Integer> dec_messages,
                       key.DecryptBfv(ciphertext, this->coeff_encoder_.get()));
  EXPECT_EQ(dec_messages, messages);
}

////////////////////////////////////////////////////////////////////////////////
// Encrypt with finite field encoder
////////////////////////////////////////////////////////////////////////////////

template <typename ModularInt>
class RnsRlweSecretKeyPackedTest : public ::testing::Test {
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

  void SetUpDefaultBgvRnsParameters() {
    auto all_params =
        testing::GetRnsParametersForFiniteFieldEncoding<ModularInt>();
    CHECK(!all_params.empty());
    SetUpBgvRnsParameters(all_params[0]);
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

  // Checking if CRT interpolation has enough precision
  // to successfully execute. If it doesn't, we skip the test until
  // Required precision should be (log_Q)+log(\sum_i q_i).
  bool HasInsufficientPrecisionForCrtInterp() const {
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
  std::unique_ptr<const FiniteFieldEncoder<ModularInt>> encoder_;
  std::unique_ptr<const RnsErrorParams<ModularInt>> error_params_;
  std::unique_ptr<Prng> prng_;
};

TYPED_TEST_SUITE(RnsRlweSecretKeyPackedTest,
                 testing::ModularIntTypesForFiniteFieldEncoding);

TYPED_TEST(RnsRlweSecretKeyPackedTest, EncryptBgvDecrypts) {
  using Integer = typename TypeParam::Int;
  for (auto const& params :
       testing::GetRnsParametersForFiniteFieldEncoding<TypeParam>()) {
    this->SetUpBgvRnsParameters(params);
    int log_n = this->rns_context_->LogN();
    ASSERT_OK_AND_ASSIGN(RnsRlweSecretKey<TypeParam> key, this->SampleKey());

    // Create a finite field encoder
    ASSERT_OK_AND_ASSIGN(auto encoder, FiniteFieldEncoder<TypeParam>::Create(
                                           this->rns_context_.get()));

    Integer t = this->rns_context_->PlaintextModulus();
    std::vector<Integer> messages = testing::SampleMessages(1 << log_n, t);
    ASSERT_OK_AND_ASSIGN(
        RnsBgvCiphertext<TypeParam> ciphertext,
        key.template EncryptBgv<FiniteFieldEncoder<TypeParam>>(
            messages, &encoder, this->error_params_.get(), this->prng_.get()));

    ASSERT_OK_AND_ASSIGN(std::vector<Integer> dec_messages,
                         key.template DecryptBgv<FiniteFieldEncoder<TypeParam>>(
                             ciphertext, &encoder));
    EXPECT_EQ(dec_messages, messages);
  }
}

TYPED_TEST(RnsRlweSecretKeyPackedTest, EncryptBgvCrtDecrypts) {
  using Integer = typename TypeParam::Int;
  using BigInteger = typename testing::CompositeModulus<Integer>::value_type;
  for (auto const& params :
       testing::GetRnsParametersForFiniteFieldEncoding<TypeParam>()) {
    this->SetUpBgvRnsParameters(params);
    if (this->HasInsufficientPrecisionForCrtInterp()) {
      continue;
    }
    int num_moduli = this->main_moduli_.size();
    auto level = num_moduli - 1;
    std::vector<BigInteger> modulus_hats =
        RnsModulusComplements<TypeParam, BigInteger>(this->main_moduli_);
    ASSERT_OK_AND_ASSIGN(std::vector<TypeParam> modulus_hat_invs,
                         this->rns_context_->MainPrimeModulusCrtFactors(level));
    int log_n = this->rns_context_->LogN();
    ASSERT_OK_AND_ASSIGN(RnsRlweSecretKey<TypeParam> key, this->SampleKey());
    // Create a finite field encoder
    ASSERT_OK_AND_ASSIGN(auto encoder, FiniteFieldEncoder<TypeParam>::Create(
                                           this->rns_context_.get()));
    Integer t = this->rns_context_->PlaintextModulus();
    std::vector<Integer> messages = testing::SampleMessages(1 << log_n, t);
    ASSERT_OK_AND_ASSIGN(
        RnsRlweCiphertext<TypeParam> ciphertext,
        key.template EncryptBgv<FiniteFieldEncoder<TypeParam>>(
            messages, &encoder, this->error_params_.get(), this->prng_.get()));

    ASSERT_OK_AND_ASSIGN(
        std::vector<Integer> dec_messages,
        (key.template DecryptBgvWithCrt<FiniteFieldEncoder<TypeParam>,
                                        BigInteger>(
            ciphertext, &encoder, absl::MakeSpan(modulus_hats),
            absl::MakeSpan(modulus_hat_invs))));
    EXPECT_EQ(dec_messages, messages);
  }
}

TYPED_TEST(RnsRlweSecretKeyPackedTest, PassWrongModuliToBgvCrtDecryption) {
  using Integer = typename TypeParam::Int;
  using BigInteger = typename testing::CompositeModulus<Integer>::value_type;
  this->SetUpDefaultBgvRnsParameters();

  int num_moduli = this->main_moduli_.size();
  auto level = num_moduli - 1;
  std::vector<BigInteger> modulus_hats =
      RnsModulusComplements<TypeParam, BigInteger>(this->main_moduli_);
  ASSERT_OK_AND_ASSIGN(std::vector<TypeParam> modulus_hat_invs,
                       this->rns_context_->MainPrimeModulusCrtFactors(level));
  int log_n = this->rns_context_->LogN();
  ASSERT_OK_AND_ASSIGN(RnsRlweSecretKey<TypeParam> key, this->SampleKey());
  // Create a finite field encoder
  ASSERT_OK_AND_ASSIGN(auto encoder, FiniteFieldEncoder<TypeParam>::Create(
                                         this->rns_context_.get()));
  Integer t = this->rns_context_->PlaintextModulus();
  std::vector<Integer> messages = testing::SampleMessages(1 << log_n, t);
  ASSERT_OK_AND_ASSIGN(
      RnsRlweCiphertext<TypeParam> ciphertext,
      key.template EncryptBgv<FiniteFieldEncoder<TypeParam>>(
          messages, &encoder, this->error_params_.get(), this->prng_.get()));

  auto status = (key.template DecryptBgvWithCrt<FiniteFieldEncoder<TypeParam>,
                                                BigInteger>(
      ciphertext, &encoder, {}, absl::MakeSpan(modulus_hat_invs)));
  EXPECT_THAT(
      status,
      StatusIs(absl::StatusCode::kInvalidArgument,
               HasSubstr("Cannot decrypt `ciphertext` with a mismatching number"
                         " of prime moduli.")));
}

TYPED_TEST(RnsRlweSecretKeyPackedTest,
           CrtDecryptBgvFailsIfCiphertextPowerOfSIsNotOne) {
  using Integer = typename TypeParam::Int;
  using BigInteger = typename testing::CompositeModulus<Integer>::value_type;
  this->SetUpDefaultBgvRnsParameters();
  int num_moduli = this->main_moduli_.size();
  auto level = num_moduli - 1;
  std::vector<BigInteger> modulus_hats =
      RnsModulusComplements<TypeParam, BigInteger>(this->main_moduli_);
  ASSERT_OK_AND_ASSIGN(std::vector<TypeParam> modulus_hat_invs,
                       this->rns_context_->MainPrimeModulusCrtFactors(level));

  ASSERT_OK_AND_ASSIGN(RnsRlweSecretKey<TypeParam> key, this->SampleKey());
  // Create a ciphertext with power_of_s != 1.
  ASSERT_OK_AND_ASSIGN(auto zero,
                       RnsPolynomial<TypeParam>::CreateZero(
                           this->rns_context_->LogN(), this->main_moduli_));
  RnsRlweCiphertext<TypeParam> ciphertext(
      /*components=*/{zero, zero}, this->main_moduli_,
      /*power_of_s=*/5, /*error=*/0, this->error_params_.get());
  ASSERT_EQ(ciphertext.PowerOfS(), 5);

  ASSERT_OK_AND_ASSIGN(auto encoder, FiniteFieldEncoder<TypeParam>::Create(
                                         this->rns_context_.get()));
  EXPECT_THAT((key.template DecryptBgvWithCrt<FiniteFieldEncoder<TypeParam>,
                                              BigInteger>(
                  ciphertext, &encoder, modulus_hats, modulus_hat_invs)),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("Cannot decrypt `ciphertext` with power of "
                                 "s not equal to 1.")));
}

TYPED_TEST(RnsRlweSecretKeyPackedTest,
           DecryptBgvWithCrtFailsIfCiphertextHasMismatchPolyDegree) {
  using Integer = typename TypeParam::Int;
  using BigInteger = typename testing::CompositeModulus<Integer>::value_type;
  this->SetUpDefaultBgvRnsParameters();
  int num_moduli = this->main_moduli_.size();
  auto level = num_moduli - 1;
  std::vector<BigInteger> modulus_hats =
      RnsModulusComplements<TypeParam, BigInteger>(this->main_moduli_);
  ASSERT_OK_AND_ASSIGN(std::vector<TypeParam> modulus_hat_invs,
                       this->rns_context_->MainPrimeModulusCrtFactors(level));

  ASSERT_OK_AND_ASSIGN(RnsRlweSecretKey<TypeParam> key, this->SampleKey());

  // Create a ciphertext with log_n that's different from the secret key.
  int log_n = this->rns_context_->LogN();
  ASSERT_OK_AND_ASSIGN(auto zero, RnsPolynomial<TypeParam>::CreateZero(
                                      log_n - 1, this->main_moduli_));
  RnsRlweCiphertext<TypeParam> ciphertext(
      /*components=*/{zero, zero}, this->main_moduli_,
      /*power_of_s=*/1, /*error=*/0, this->error_params_.get());
  ASSERT_EQ(ciphertext.LogN(), log_n - 1);
  ASSERT_OK_AND_ASSIGN(auto encoder, FiniteFieldEncoder<TypeParam>::Create(
                                         this->rns_context_.get()));
  EXPECT_THAT((key.template DecryptBgvWithCrt<FiniteFieldEncoder<TypeParam>,
                                              BigInteger>(
                  ciphertext, &encoder, modulus_hats, modulus_hat_invs)),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("Cannot decrypt `ciphertext` with a "
                                 "mismatching polynomial degree.")));
}

TYPED_TEST(RnsRlweSecretKeyPackedTest, DecryptBgvWithCrtFailsIfEncoderIsNull) {
  using Integer = typename TypeParam::Int;
  using BigInteger = typename testing::CompositeModulus<Integer>::value_type;
  this->SetUpDefaultBgvRnsParameters();
  int num_moduli = this->main_moduli_.size();
  auto level = num_moduli - 1;
  std::vector<BigInteger> modulus_hats =
      RnsModulusComplements<TypeParam, BigInteger>(this->main_moduli_);
  ASSERT_OK_AND_ASSIGN(std::vector<TypeParam> modulus_hat_invs,
                       this->rns_context_->MainPrimeModulusCrtFactors(level));
  int log_n = this->rns_context_->LogN();
  ASSERT_OK_AND_ASSIGN(RnsRlweSecretKey<TypeParam> key, this->SampleKey());
  ASSERT_OK_AND_ASSIGN(auto encoder, FiniteFieldEncoder<TypeParam>::Create(
                                         this->rns_context_.get()));
  Integer t = this->rns_context_->PlaintextModulus();
  std::vector<Integer> messages = testing::SampleMessages(1 << log_n, t);
  ASSERT_OK_AND_ASSIGN(
      RnsRlweCiphertext<TypeParam> ciphertext,
      key.template EncryptBgv<FiniteFieldEncoder<TypeParam>>(
          messages, &encoder, this->error_params_.get(), this->prng_.get()));

  EXPECT_THAT((key.template DecryptBgvWithCrt<FiniteFieldEncoder<TypeParam>,
                                              BigInteger>(
                  ciphertext, /*encoder=*/nullptr, absl::MakeSpan(modulus_hats),
                  absl::MakeSpan(modulus_hat_invs))),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`encoder` must not be null.")));
}

TYPED_TEST(RnsRlweSecretKeyPackedTest, EncryptBfvDecrypts) {
  using Integer = typename TypeParam::Int;

  for (auto const& params :
       testing::GetBfvParametersForFiniteFieldEncoding<TypeParam>()) {
    this->SetUpBfvRnsParameters(params);
    int log_n = this->rns_context_->LogN();
    Integer t = this->rns_context_->PlaintextModulus();
    ASSERT_OK_AND_ASSIGN(RnsRlweSecretKey<TypeParam> key, this->SampleKey());

    // Create a finite field encoder
    ASSERT_OK_AND_ASSIGN(auto encoder, FiniteFieldEncoder<TypeParam>::Create(
                                           this->rns_context_.get()));

    std::vector<Integer> messages = testing::SampleMessages(1 << log_n, t);
    ASSERT_OK_AND_ASSIGN(
        RnsBfvCiphertext<TypeParam> ciphertext,
        key.template EncryptBfv<FiniteFieldEncoder<TypeParam>>(
            messages, &encoder, this->error_params_.get(), this->prng_.get()));

    ASSERT_OK_AND_ASSIGN(std::vector<Integer> dec_messages,
                         key.template DecryptBfv<FiniteFieldEncoder<TypeParam>>(
                             ciphertext, &encoder));
    EXPECT_EQ(dec_messages, messages);
  }
}

}  // namespace
}  // namespace rlwe
