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

#include "shell_encryption/rns/rns_relinearization_key.h"

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
class RnsRelinKeyTest : public ::testing::Test {
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
    moduli_ = rns_context_->MainPrimeModuli();

    int level = params.qs.size() - 1;
    auto q_hats = rns_context_->MainPrimeModulusComplements(level).value();
    auto q_hat_invs = rns_context_->MainPrimeModulusCrtFactors(level).value();

    std::vector<size_t> log_bs(params.qs.size(), params.log_gadget_base);
    auto gadget = RnsGadget<ModularInt>::Create(params.log_n, log_bs, q_hats,
                                                q_hat_invs, moduli_);
    CHECK(gadget.ok());
    gadget_ = std::make_unique<const RnsGadget<ModularInt>>(
        std::move(gadget.value()));

    auto error_params = RnsErrorParams<ModularInt>::Create(
        params.log_n, moduli_, /*aux_moduli=*/{},
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
    moduli_ = rns_context_->MainPrimeModuli();

    int level = params.qs.size() - 1;
    auto q_hats = rns_context_->MainPrimeModulusComplements(level).value();
    auto q_hat_invs = rns_context_->MainPrimeModulusCrtFactors(level).value();

    std::vector<size_t> log_bs(params.qs.size(), params.log_gadget_base);
    auto gadget = RnsGadget<ModularInt>::Create(params.log_n, log_bs, q_hats,
                                                q_hat_invs, moduli_);
    CHECK(gadget.ok());
    gadget_ = std::make_unique<const RnsGadget<ModularInt>>(
        std::move(gadget.value()));

    auto error_params = RnsErrorParams<ModularInt>::Create(
        params.log_n, moduli_, /*aux_moduli=*/{},
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
                                                moduli_, prng_.get());
  }

  std::unique_ptr<const RnsContext<ModularInt>> rns_context_;
  std::vector<const PrimeModulus<ModularInt>*> moduli_;
  std::unique_ptr<const RnsGadget<ModularInt>> gadget_;
  std::unique_ptr<const RnsErrorParams<ModularInt>> error_params_;
  std::unique_ptr<Prng> prng_;
};

TYPED_TEST_SUITE(RnsRelinKeyTest,
                 testing::ModularIntTypesForFiniteFieldEncoding);

template <typename ModularInt>
class RnsRelinKeyNegativeTest : public RnsRelinKeyTest<ModularInt> {};

TYPED_TEST_SUITE(RnsRelinKeyNegativeTest,
                 testing::ModularIntTypesForNegativeTests);

TYPED_TEST(RnsRelinKeyNegativeTest, CreateFailsIfDegreeIsLessThan2) {
  constexpr absl::string_view k_error = "`degree` must be at least 2";

  this->SetUpDefaultBgvParameters();
  ASSERT_OK_AND_ASSIGN(RnsRlweSecretKey<TypeParam> secret_key,
                       this->SampleSecretKey());

  // Negative degree.
  EXPECT_THAT(RnsRelinKey<TypeParam>::CreateForBgv(
                  secret_key, /*degree=*/-1, kVariance, this->gadget_.get(),
                  this->rns_context_->PlaintextModulus(), kPrngType),
              StatusIs(absl::StatusCode::kInvalidArgument, HasSubstr(k_error)));
  // Zero degree.
  EXPECT_THAT(RnsRelinKey<TypeParam>::CreateForBgv(
                  secret_key, /*degree=*/0, kVariance, this->gadget_.get(),
                  this->rns_context_->PlaintextModulus(), kPrngType),
              StatusIs(absl::StatusCode::kInvalidArgument, HasSubstr(k_error)));
  // degree = 1.
  EXPECT_THAT(RnsRelinKey<TypeParam>::CreateForBgv(
                  secret_key, /*degree=*/1, kVariance, this->gadget_.get(),
                  this->rns_context_->PlaintextModulus(), kPrngType),
              StatusIs(absl::StatusCode::kInvalidArgument, HasSubstr(k_error)));
}

TYPED_TEST(RnsRelinKeyNegativeTest, CreateFailsIfVarianceIsNotPositive) {
  this->SetUpDefaultBgvParameters();
  ASSERT_OK_AND_ASSIGN(RnsRlweSecretKey<TypeParam> secret_key,
                       this->SampleSecretKey());
  EXPECT_THAT(
      RnsRelinKey<TypeParam>::CreateForBgv(
          secret_key, /*degree=*/2, /*variance=*/-1, this->gadget_.get(),
          this->rns_context_->PlaintextModulus(), kPrngType),
      StatusIs(absl::StatusCode::kInvalidArgument,
               HasSubstr("`variance` must be positive")));
}

TYPED_TEST(RnsRelinKeyNegativeTest, CreateFailsIfGadgetIsNull) {
  this->SetUpDefaultBgvParameters();
  ASSERT_OK_AND_ASSIGN(RnsRlweSecretKey<TypeParam> secret_key,
                       this->SampleSecretKey());
  EXPECT_THAT(RnsRelinKey<TypeParam>::CreateForBgv(
                  secret_key, /*degree=*/2, kVariance, /*gadget=*/nullptr,
                  this->rns_context_->PlaintextModulus(), kPrngType),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`gadget` must not be null")));
}

TYPED_TEST(RnsRelinKeyNegativeTest, CreateFailsIfInvalidPrngType) {
  this->SetUpDefaultBgvParameters();
  ASSERT_OK_AND_ASSIGN(RnsRlweSecretKey<TypeParam> secret_key,
                       this->SampleSecretKey());
  EXPECT_THAT(RnsRelinKey<TypeParam>::CreateForBgv(
                  secret_key, /*degree=*/2, kVariance, this->gadget_.get(),
                  this->rns_context_->PlaintextModulus(), PRNG_TYPE_INVALID),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("PrngType not specified correctly")));
}

TYPED_TEST(RnsRelinKeyNegativeTest,
           SampleRandomPadFailsIfDimensionIsNotPositive) {
  this->SetUpDefaultBgvParameters();
  ASSERT_OK_AND_ASSIGN(RnsRlweSecretKey<TypeParam> secret_key,
                       this->SampleSecretKey());
  EXPECT_THAT(RnsRelinKey<TypeParam>::SampleRandomPad(
                  /*dimension=*/0, /*degree=*/2, /*log_n=*/1, this->moduli_,
                  kPrngSeed, kPrngType),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`dimension` must be positive")));
}

TYPED_TEST(RnsRelinKeyNegativeTest, SampleRandomPadFailsIfDegreeIsLessThan2) {
  this->SetUpDefaultBgvParameters();
  ASSERT_OK_AND_ASSIGN(RnsRlweSecretKey<TypeParam> secret_key,
                       this->SampleSecretKey());
  EXPECT_THAT(RnsRelinKey<TypeParam>::SampleRandomPad(
                  /*dimension=*/1, /*degree=*/0, /*log_n=*/1, this->moduli_,
                  kPrngSeed, kPrngType),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`degree` must be at least 2")));
}

TYPED_TEST(RnsRelinKeyNegativeTest, SampleRandomPadFailsIfLogNIsNotPositive) {
  this->SetUpDefaultBgvParameters();
  ASSERT_OK_AND_ASSIGN(RnsRlweSecretKey<TypeParam> secret_key,
                       this->SampleSecretKey());
  EXPECT_THAT(RnsRelinKey<TypeParam>::SampleRandomPad(
                  /*dimension=*/1, /*degree=*/2, /*log_n=*/-1, this->moduli_,
                  kPrngSeed, kPrngType),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`log_n` must be positive")));
}

TYPED_TEST(RnsRelinKeyNegativeTest, SampleRandomPadFailsIfInvalidPrngType) {
  this->SetUpDefaultBgvParameters();
  ASSERT_OK_AND_ASSIGN(RnsRlweSecretKey<TypeParam> secret_key,
                       this->SampleSecretKey());
  EXPECT_THAT(RnsRelinKey<TypeParam>::SampleRandomPad(
                  /*dimension=*/1, /*degree=*/2, /*log_n=*/1, this->moduli_,
                  kPrngSeed, PRNG_TYPE_INVALID),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("PrngType not specified correctly")));
}

TYPED_TEST(RnsRelinKeyNegativeTest, ApplyToFailsIfCiphertextDegreeIsTooLarge) {
  constexpr int k_rk_degree = 2;

  this->SetUpDefaultBgvParameters();
  ASSERT_OK_AND_ASSIGN(RnsRlweSecretKey<TypeParam> secret_key,
                       this->SampleSecretKey());

  ASSERT_OK_AND_ASSIGN(
      RnsRelinKey<TypeParam> rk,
      RnsRelinKey<TypeParam>::CreateForBgv(
          secret_key, k_rk_degree, kVariance, this->gadget_.get(),
          this->rns_context_->PlaintextModulus(), kPrngType));

  // Create a zero degree-3 ciphertext (0, 0, 0, 0) for negative testing only.
  int log_n = this->rns_context_->LogN();
  ASSERT_OK_AND_ASSIGN(
      RnsPolynomial<TypeParam> zero,
      RnsPolynomial<TypeParam>::CreateZero(log_n, this->moduli_,
                                           /*is_ntt=*/false));
  RnsBgvCiphertext<TypeParam> ciphertext(
      {zero, zero, zero, zero}, this->moduli_,
      /*power=*/1, /*error=*/0, this->error_params_.get());

  // ApplyTo will fail since ciphertext has a mismatch PowerOfS.
  EXPECT_THAT(rk.ApplyTo(ciphertext),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`ciphertext` degree is larger")));
}

TYPED_TEST(RnsRelinKeyNegativeTest,
           ApplyToFailsIfCiphertextHasMismatchRnsModuli) {
  // This test requires at least two RNS moduli.
  testing::RnsParameters<TypeParam> params;
  for (auto const& params_candidate :
       testing::GetRnsParametersForFiniteFieldEncoding<TypeParam>()) {
    if (params_candidate.qs.size() >= 2) {
      params = params_candidate;
      break;
    }
  }
  ASSERT_NE(params.log_n, 0)
      << "Cannot find parameters with at least two prime moduli";

  constexpr int k_degree = 2;
  this->SetUpBgvParameters(params);

  ASSERT_OK_AND_ASSIGN(RnsRlweSecretKey<TypeParam> secret_key,
                       this->SampleSecretKey());
  ASSERT_OK_AND_ASSIGN(RnsRelinKey<TypeParam> rk,
                       RnsRelinKey<TypeParam>::CreateForBgv(
                           secret_key, k_degree, kVariance, this->gadget_.get(),
                           this->rns_context_->PlaintextModulus(), kPrngType));

  // Create a zero ciphertext (0, 0) with a smaller set of RNS moduli.
  int log_n = this->rns_context_->LogN();
  std::vector<const PrimeModulus<TypeParam>*> small_moduli = this->moduli_;
  small_moduli.pop_back();
  ASSERT_OK_AND_ASSIGN(RnsPolynomial<TypeParam> zero,
                       RnsPolynomial<TypeParam>::CreateZero(log_n, small_moduli,
                                                            /*is_ntt=*/false));
  RnsBgvCiphertext<TypeParam> ciphertext(
      {zero, zero}, small_moduli, /*power=*/1,
      /*error=*/0, this->error_params_.get());

  // ApplyTo will fail since ciphertext has a mismatch RNS moduli set.
  EXPECT_THAT(
      rk.ApplyTo(ciphertext),
      StatusIs(
          absl::StatusCode::kInvalidArgument,
          HasSubstr("`ciphertext` does not have a matching RNS moduli set")));
}

TYPED_TEST(RnsRelinKeyNegativeTest, ApplyToFailsIfCiphertextPowerOfSIsNot1) {
  constexpr int k_rk_degree = 2;
  constexpr int k_incorrect_power = 5;

  this->SetUpDefaultBgvParameters();
  ASSERT_OK_AND_ASSIGN(RnsRlweSecretKey<TypeParam> secret_key,
                       this->SampleSecretKey());

  ASSERT_OK_AND_ASSIGN(
      RnsRelinKey<TypeParam> rk,
      RnsRelinKey<TypeParam>::CreateForBgv(
          secret_key, k_rk_degree, kVariance, this->gadget_.get(),
          this->rns_context_->PlaintextModulus(), kPrngType));

  // Mock up a zero ciphertext (0, 0) with PowerOfS == 5.
  int log_n = this->rns_context_->LogN();
  ASSERT_OK_AND_ASSIGN(
      RnsPolynomial<TypeParam> zero,
      RnsPolynomial<TypeParam>::CreateZero(log_n, this->moduli_,
                                           /*is_ntt=*/false));
  RnsBgvCiphertext<TypeParam> ciphertext({zero, zero}, this->moduli_,
                                         k_incorrect_power, /*error=*/0,
                                         this->error_params_.get());

  // ApplyTo will fail since ciphertext has a mismatch PowerOfS.
  EXPECT_THAT(rk.ApplyTo(ciphertext),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("Relinearization key can only apply to a "
                                 "ciphertext of power 1")));
}

////////////////////////////////////////////////////////////////////////////////
// Functional tests
////////////////////////////////////////////////////////////////////////////////

TYPED_TEST(RnsRelinKeyTest, CreateForBgvSucceeds) {
  using Integer = typename TypeParam::Int;
  constexpr int k_degree = 2;

  // Using the default parameter for each integer type is sufficient.
  this->SetUpDefaultBgvParameters();
  // Keep the same secret key and a ciphertext that encrypts values in slots.
  ASSERT_OK_AND_ASSIGN(RnsRlweSecretKey<TypeParam> secret_key,
                       this->SampleSecretKey());

  auto gadget = this->gadget_.get();
  Integer t = this->rns_context_->PlaintextModulus();
  for (auto prng_type : {PRNG_TYPE_CHACHA, PRNG_TYPE_HKDF}) {
    ASSERT_OK_AND_ASSIGN(
        RnsRelinKey<TypeParam> rk,
        RnsRelinKey<TypeParam>::CreateForBgv(secret_key, k_degree, kVariance,
                                             gadget, t, prng_type));
    ASSERT_EQ(rk.Gadget(), gadget);
    ASSERT_EQ(rk.Dimension(), gadget->Dimension());
    ASSERT_EQ(rk.GetKeyA().size(), gadget->Dimension());
    ASSERT_EQ(rk.GetKeyB().size(), gadget->Dimension());
    ASSERT_EQ(rk.Degree(), k_degree);
  }
}

TYPED_TEST(RnsRelinKeyTest, CreateForBfvSucceeds) {
  constexpr int k_degree = 2;

  // Using the default parameter for each integer type is sufficient.
  this->SetUpDefaultBgvParameters();
  // Keep the same secret key and a ciphertext that encrypts values in slots.
  ASSERT_OK_AND_ASSIGN(RnsRlweSecretKey<TypeParam> secret_key,
                       this->SampleSecretKey());

  auto gadget = this->gadget_.get();
  for (auto prng_type : {PRNG_TYPE_CHACHA, PRNG_TYPE_HKDF}) {
    ASSERT_OK_AND_ASSIGN(
        RnsRelinKey<TypeParam> rk,
        RnsRelinKey<TypeParam>::CreateForBfv(secret_key, k_degree, kVariance,
                                             gadget, prng_type));
    ASSERT_EQ(rk.Gadget(), gadget);
    ASSERT_EQ(rk.Dimension(), gadget->Dimension());
    ASSERT_EQ(rk.GetKeyA().size(), gadget->Dimension());
    ASSERT_EQ(rk.GetKeyB().size(), gadget->Dimension());
    ASSERT_EQ(rk.Degree(), k_degree);
  }
}

TYPED_TEST(RnsRelinKeyTest, ApplyToBgvCiphertextOfDegree2) {
  using Encoder = FiniteFieldEncoder<TypeParam>;
  using Integer = typename TypeParam::Int;

  constexpr int k_degree = 2;

  for (auto const& params :
       testing::GetRnsParametersForFiniteFieldEncoding<TypeParam>()) {
    this->SetUpBgvParameters(params);

    int log_n = this->rns_context_->LogN();
    int n = 1 << log_n;
    double expected_error_bound =
        this->error_params_->B_secretkey_encryption() *
        this->error_params_->B_secretkey_encryption() * sqrt(n);
    double composite_modulus = 1;
    for (auto prime_modulus : this->moduli_) {
      composite_modulus *= static_cast<double>(prime_modulus->Modulus());
    }
    if (expected_error_bound >= composite_modulus) {
      // The test parameters are not suitable for homomorphic multiplication
      // as the error in the product ciphertext is expected to be larger than
      // the ciphertext modulus.
      GTEST_SKIP() << "Insufficient modulus for multiplication.";
    }

    ASSERT_OK_AND_ASSIGN(auto encoder,
                         Encoder::Create(this->rns_context_.get()));

    ASSERT_OK_AND_ASSIGN(RnsRlweSecretKey<TypeParam> secret_key,
                         this->SampleSecretKey());

    auto gadget = this->gadget_.get();
    auto prng = this->prng_.get();
    Integer t = this->rns_context_->PlaintextModulus();
    ASSERT_OK_AND_ASSIGN(
        RnsRelinKey<TypeParam> rk,
        RnsRelinKey<TypeParam>::CreateForBgv(secret_key, k_degree, kVariance,
                                             gadget, t, kPrngType));
    ASSERT_EQ(rk.Gadget(), gadget);
    ASSERT_EQ(rk.Dimension(), gadget->Dimension());
    ASSERT_EQ(rk.Degree(), k_degree);

    std::vector<Integer> values0 = testing::SampleMessages<Integer>(n, t);
    std::vector<Integer> values1 = testing::SampleMessages<Integer>(n, t);
    ASSERT_OK_AND_ASSIGN(
        RnsBgvCiphertext<TypeParam> ciphertext0,
        secret_key.template EncryptBgv<Encoder>(
            values0, &encoder, this->error_params_.get(), prng));
    ASSERT_OK_AND_ASSIGN(
        RnsBgvCiphertext<TypeParam> ciphertext1,
        secret_key.template EncryptBgv<Encoder>(
            values1, &encoder, this->error_params_.get(), prng));

    // Multiply the two ciphertexts, and the result is a degree 2 ciphertext.
    ASSERT_OK_AND_ASSIGN(RnsBgvCiphertext<TypeParam> ciphertext_product,
                         ciphertext0 * ciphertext1);
    EXPECT_EQ(ciphertext_product.Degree(), 2);
    EXPECT_EQ(ciphertext_product.Level(), ciphertext0.Level());

    // Relinearize the degree-2 ciphertext.
    ASSERT_OK_AND_ASSIGN(RnsBgvCiphertext<TypeParam> ciphertext_linear,
                         rk.ApplyTo(ciphertext_product));
    EXPECT_EQ(ciphertext_linear.Degree(), 1);
    EXPECT_EQ(ciphertext_linear.Level(), ciphertext0.Level());

    // Check for decryption results.
    ASSERT_OK_AND_ASSIGN(
        auto decryptions,
        secret_key.template DecryptBgv<Encoder>(ciphertext_linear, &encoder));
    ASSERT_EQ(decryptions.size(), n);
    for (int i = 0; i < decryptions.size(); ++i) {
      Integer expected = (values0[i] * values1[i]) % t;
      EXPECT_EQ(decryptions[i], expected);
    }
  }
}

TYPED_TEST(RnsRelinKeyTest, ApplyToBfvCiphertextOfDegree2) {
  using Encoder = FiniteFieldEncoder<TypeParam>;
  using Integer = typename TypeParam::Int;

  constexpr int k_degree = 2;

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
    this->SetUpBfvParameters(params);

    ASSERT_OK_AND_ASSIGN(auto encoder,
                         Encoder::Create(this->rns_context_.get()));

    ASSERT_OK_AND_ASSIGN(RnsRlweSecretKey<TypeParam> secret_key,
                         this->SampleSecretKey());

    auto gadget = this->gadget_.get();
    auto prng = this->prng_.get();
    Integer t = this->rns_context_->PlaintextModulus();
    ASSERT_OK_AND_ASSIGN(
        RnsRelinKey<TypeParam> rk,
        RnsRelinKey<TypeParam>::CreateForBfv(secret_key, k_degree, kVariance,
                                             gadget, kPrngType));
    ASSERT_EQ(rk.Gadget(), gadget);
    ASSERT_EQ(rk.Dimension(), gadget->Dimension());
    ASSERT_EQ(rk.Degree(), k_degree);

    int num_coeffs = 1 << this->rns_context_->LogN();
    std::vector<Integer> values0 =
        testing::SampleMessages<Integer>(num_coeffs, t);
    std::vector<Integer> values1 =
        testing::SampleMessages<Integer>(num_coeffs, t);
    ASSERT_OK_AND_ASSIGN(
        RnsBfvCiphertext<TypeParam> ciphertext0,
        secret_key.template EncryptBfv<Encoder>(
            values0, &encoder, this->error_params_.get(), prng));
    ASSERT_OK_AND_ASSIGN(
        RnsBfvCiphertext<TypeParam> ciphertext1,
        secret_key.template EncryptBfv<Encoder>(
            values1, &encoder, this->error_params_.get(), prng));

    // Multiply the two ciphertexts, and the result is a degree 2 ciphertext.
    ASSERT_OK_AND_ASSIGN(RnsBfvCiphertext<TypeParam> ciphertext_product,
                         ciphertext0 * ciphertext1);
    EXPECT_EQ(ciphertext_product.Degree(), 2);
    EXPECT_EQ(ciphertext_product.Level(), ciphertext0.Level());

    // Relinearize the degree-2 ciphertext.
    ASSERT_OK_AND_ASSIGN(RnsBfvCiphertext<TypeParam> ciphertext_linear,
                         rk.ApplyTo(ciphertext_product));
    EXPECT_EQ(ciphertext_linear.Degree(), 1);
    EXPECT_EQ(ciphertext_linear.Level(), ciphertext0.Level());

    // Check for decryption results.
    ASSERT_OK_AND_ASSIGN(
        auto decryptions,
        secret_key.template DecryptBfv<Encoder>(ciphertext_linear, &encoder));
    ASSERT_EQ(decryptions.size(), num_coeffs);
    for (int i = 0; i < decryptions.size(); ++i) {
      Integer expected = (values0[i] * values1[i]) % t;
      EXPECT_EQ(decryptions[i], expected);
    }
  }
}

}  // namespace
}  // namespace rlwe
