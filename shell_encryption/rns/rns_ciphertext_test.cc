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

#include "shell_encryption/rns/rns_ciphertext.h"

#include <cmath>
#include <memory>
#include <vector>

#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include "absl/log/check.h"
#include "absl/random/random.h"
#include "absl/status/status.h"
#include "absl/strings/string_view.h"
#include "shell_encryption/integral_types.h"
#include "shell_encryption/prng/single_thread_hkdf_prng.h"
#include "shell_encryption/rns/coefficient_encoder.h"
#include "shell_encryption/rns/rns_context.h"
#include "shell_encryption/rns/rns_error_params.h"
#include "shell_encryption/rns/rns_modulus.h"
#include "shell_encryption/rns/rns_polynomial.h"
#include "shell_encryption/rns/rns_secret_key.h"
#include "shell_encryption/rns/rns_serialization.pb.h"
#include "shell_encryption/rns/testing/parameters.h"
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
class RnsCiphertextTest : public ::testing::Test {
 protected:
  void SetUp() override {
    // Create the RNS context.
    testing::RnsParameters<ModularInt> rns_params =
        testing::GetRnsParameters<ModularInt>();
    auto rns_context = RnsContext<ModularInt>::Create(
        rns_params.log_n, rns_params.qs, rns_params.ps, rns_params.t);
    CHECK(rns_context.ok());
    rns_context_ = std::make_unique<const RnsContext<ModularInt>>(
        std::move(rns_context.value()));
    moduli_ = rns_context_->MainPrimeModuli();

    // Create the coefficient encoder.
    auto coeff_encoder =
        CoefficientEncoder<ModularInt>::Create(rns_context_.get());
    CHECK(coeff_encoder.ok());
    coeff_encoder_ = std::make_unique<const CoefficientEncoder<ModularInt>>(
        std::move(coeff_encoder.value()));

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
  }

  // Sample a secret key from prng.
  absl::StatusOr<RnsRlweSecretKey<ModularInt>> SampleKey() const {
    return RnsRlweSecretKey<ModularInt>::Sample(rns_context_->LogN(), kVariance,
                                                rns_context_->MainPrimeModuli(),
                                                prng_.get());
  }

  std::unique_ptr<const RnsContext<ModularInt>> rns_context_;
  std::vector<const PrimeModulus<ModularInt>*> moduli_;
  std::unique_ptr<const CoefficientEncoder<ModularInt>> coeff_encoder_;
  std::unique_ptr<const RnsErrorParams<ModularInt>> error_params_;
  std::unique_ptr<Prng> prng_;
};

TYPED_TEST_SUITE(RnsCiphertextTest, testing::ModularIntTypes);

TYPED_TEST(RnsCiphertextTest, AddFailsIfDegreeMismatch) {
  ASSERT_OK_AND_ASSIGN(auto zero,
                       RnsPolynomial<TypeParam>::CreateZero(
                           this->rns_context_->LogN(), this->moduli_));

  // Create a ciphertext with two components (degree = 1).
  RnsRlweCiphertext<TypeParam> ct0(/*components=*/{zero, zero}, this->moduli_,
                                   /*power_of_s=*/1, /*error=*/0,
                                   this->error_params_.get());
  ASSERT_EQ(ct0.Degree(), 1);
  // Create a ciphertext with three components (degree = 2).
  RnsRlweCiphertext<TypeParam> ct1(
      /*components=*/{zero, zero, zero}, this->moduli_,
      /*power_of_s=*/1, /*error=*/0, this->error_params_.get());
  ASSERT_EQ(ct1.Degree(), 2);

  EXPECT_THAT(ct0.AddInPlace(ct1),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`that` has a mismatched degree")));
}

TYPED_TEST(RnsCiphertextTest, AddFailsIfLevelMismatch) {
  // This test requires at least two RNS moduli.
  if (this->moduli_.size() <= 1) {
    return;
  }

  // Create a ciphertext at highest level.
  int curr_level = this->moduli_.size() - 1;
  ASSERT_OK_AND_ASSIGN(auto zero,
                       RnsPolynomial<TypeParam>::CreateZero(
                           this->rns_context_->LogN(), this->moduli_));
  RnsRlweCiphertext<TypeParam> ct0(/*components=*/{zero, zero}, this->moduli_,
                                   /*power_of_s=*/1, /*error=*/0,
                                   this->error_params_.get());
  ASSERT_EQ(ct0.Level(), curr_level);

  // Create a ciphertext with fewer RNS moduli (so at a smaller level).
  std::vector<const PrimeModulus<TypeParam>*> moduli_reduced = this->moduli_;
  moduli_reduced.pop_back();
  ASSERT_OK_AND_ASSIGN(auto zero_reduced,
                       RnsPolynomial<TypeParam>::CreateZero(
                           this->rns_context_->LogN(), moduli_reduced));
  RnsRlweCiphertext<TypeParam> ct1(
      /*components=*/{zero_reduced, zero_reduced}, moduli_reduced,
      /*power_of_s=*/1, /*error=*/0, this->error_params_.get());
  ASSERT_EQ(ct1.Level(), curr_level - 1);

  EXPECT_THAT(ct0.AddInPlace(ct1),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`that` has a mismatched level")));
  EXPECT_THAT(ct0.AddInPlaceWithoutPad(ct1),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`that` has a mismatched level")));
}

TYPED_TEST(RnsCiphertextTest, AddFailsIfPowerOfSMismatch) {
  ASSERT_OK_AND_ASSIGN(auto zero,
                       RnsPolynomial<TypeParam>::CreateZero(
                           this->rns_context_->LogN(), this->moduli_));

  // Create a ciphertext with power_of_s == 1.
  RnsRlweCiphertext<TypeParam> ct0(/*components=*/{zero, zero}, this->moduli_,
                                   /*power_of_s=*/1, /*error=*/0,
                                   this->error_params_.get());
  ASSERT_EQ(ct0.PowerOfS(), 1);

  // Create a ciphertext with power_of_s == 5.
  RnsRlweCiphertext<TypeParam> ct1(
      /*components=*/{zero, zero}, this->moduli_,
      /*power_of_s=*/5, /*error=*/0, this->error_params_.get());
  ASSERT_EQ(ct1.PowerOfS(), 5);

  EXPECT_THAT(
      ct0.AddInPlace(ct1),
      StatusIs(absl::StatusCode::kInvalidArgument,
               HasSubstr("`that` is encrypted with a different key power")));
  EXPECT_THAT(
      ct0.AddInPlaceWithoutPad(ct1),
      StatusIs(absl::StatusCode::kInvalidArgument,
               HasSubstr("`that` is encrypted with a different key power")));
}

TYPED_TEST(RnsCiphertextTest, AddPlaintextFailsIfLevelMismatch) {
  // This test requires at least two RNS moduli.
  if (this->moduli_.size() <= 1) {
    return;
  }

  // Create a ciphertext at highest level.
  int curr_level = this->moduli_.size() - 1;
  ASSERT_OK_AND_ASSIGN(auto zero,
                       RnsPolynomial<TypeParam>::CreateZero(
                           this->rns_context_->LogN(), this->moduli_));
  RnsRlweCiphertext<TypeParam> ct0(/*components=*/{zero, zero}, this->moduli_,
                                   /*power_of_s=*/1, /*error=*/0,
                                   this->error_params_.get());
  ASSERT_EQ(ct0.Level(), curr_level);

  // Create a plaintext polynomial with fewer RNS moduli (so at a lower level).
  std::vector<const PrimeModulus<TypeParam>*> moduli_reduced = this->moduli_;
  moduli_reduced.pop_back();
  ASSERT_OK_AND_ASSIGN(auto zero_reduced,
                       RnsPolynomial<TypeParam>::CreateZero(
                           this->rns_context_->LogN(), moduli_reduced));
  ASSERT_EQ(zero_reduced.NumModuli(), this->moduli_.size() - 1);
  EXPECT_THAT(ct0.AddInPlace(zero_reduced),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`plaintext` has a mismatched level")));
}

TYPED_TEST(RnsCiphertextTest, SubFailsIfDegreeMismatch) {
  ASSERT_OK_AND_ASSIGN(auto zero,
                       RnsPolynomial<TypeParam>::CreateZero(
                           this->rns_context_->LogN(), this->moduli_));

  // Create a ciphertext with two components (degree = 1).
  RnsRlweCiphertext<TypeParam> ct0(/*components=*/{zero, zero}, this->moduli_,
                                   /*power_of_s=*/1, /*error=*/0,
                                   this->error_params_.get());
  ASSERT_EQ(ct0.Degree(), 1);
  // Create a ciphertext with three components (degree = 2).
  RnsRlweCiphertext<TypeParam> ct1(
      /*components=*/{zero, zero, zero}, this->moduli_,
      /*power_of_s=*/1, /*error=*/0, this->error_params_.get());
  ASSERT_EQ(ct1.Degree(), 2);

  EXPECT_THAT(ct0.SubInPlace(ct1),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`that` has a mismatched degree")));
}

TYPED_TEST(RnsCiphertextTest, SubFailsIfLevelMismatch) {
  // This test requires at least two RNS moduli.
  if (this->moduli_.size() <= 1) {
    return;
  }

  // Create a ciphertext at highest level.
  int curr_level = this->moduli_.size() - 1;
  ASSERT_OK_AND_ASSIGN(auto zero,
                       RnsPolynomial<TypeParam>::CreateZero(
                           this->rns_context_->LogN(), this->moduli_));
  RnsRlweCiphertext<TypeParam> ct0(/*components=*/{zero, zero}, this->moduli_,
                                   /*power_of_s=*/1, /*error=*/0,
                                   this->error_params_.get());
  ASSERT_EQ(ct0.Level(), curr_level);

  // Create a ciphertext with fewer RNS moduli (so at a smaller level).
  std::vector<const PrimeModulus<TypeParam>*> moduli_reduced = this->moduli_;
  moduli_reduced.pop_back();
  ASSERT_OK_AND_ASSIGN(auto zero_reduced,
                       RnsPolynomial<TypeParam>::CreateZero(
                           this->rns_context_->LogN(), moduli_reduced));
  RnsRlweCiphertext<TypeParam> ct1(
      /*components=*/{zero_reduced, zero_reduced}, moduli_reduced,
      /*power_of_s=*/1, /*error=*/0, this->error_params_.get());
  ASSERT_EQ(ct1.Level(), curr_level - 1);

  EXPECT_THAT(ct0.SubInPlace(ct1),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`that` has a mismatched level")));
  EXPECT_THAT(ct0.SubInPlaceWithoutPad(ct1),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`that` has a mismatched level")));
}

TYPED_TEST(RnsCiphertextTest, SubPlaintextFailsIfLevelMismatch) {
  // This test requires at least two RNS moduli.
  if (this->moduli_.size() <= 1) {
    return;
  }

  // Create a ciphertext at highest level.
  int curr_level = this->moduli_.size() - 1;
  ASSERT_OK_AND_ASSIGN(auto zero,
                       RnsPolynomial<TypeParam>::CreateZero(
                           this->rns_context_->LogN(), this->moduli_));
  RnsRlweCiphertext<TypeParam> ct0(/*components=*/{zero, zero}, this->moduli_,
                                   /*power_of_s=*/1, /*error=*/0,
                                   this->error_params_.get());
  ASSERT_EQ(ct0.Level(), curr_level);

  // Create a plaintext polynomial with fewer RNS moduli (so at a lower level).
  std::vector<const PrimeModulus<TypeParam>*> moduli_reduced = this->moduli_;
  moduli_reduced.pop_back();
  ASSERT_OK_AND_ASSIGN(auto zero_reduced,
                       RnsPolynomial<TypeParam>::CreateZero(
                           this->rns_context_->LogN(), moduli_reduced));
  ASSERT_EQ(zero_reduced.NumModuli(), this->moduli_.size() - 1);

  EXPECT_THAT(ct0.SubInPlace(zero_reduced),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`plaintext` has a mismatched level")));
}

TYPED_TEST(RnsCiphertextTest, SubFailsIfPowerOfSMismatch) {
  ASSERT_OK_AND_ASSIGN(auto zero,
                       RnsPolynomial<TypeParam>::CreateZero(
                           this->rns_context_->LogN(), this->moduli_));

  // Create a ciphertext with power_of_s == 1.
  RnsRlweCiphertext<TypeParam> ct0(/*components=*/{zero, zero}, this->moduli_,
                                   /*power_of_s=*/1, /*error=*/0,
                                   this->error_params_.get());
  ASSERT_EQ(ct0.PowerOfS(), 1);

  // Create a ciphertext with power_of_s == 5.
  RnsRlweCiphertext<TypeParam> ct1(
      /*components=*/{zero, zero}, this->moduli_,
      /*power_of_s=*/5, /*error=*/0, this->error_params_.get());
  ASSERT_EQ(ct1.PowerOfS(), 5);

  EXPECT_THAT(
      ct0.SubInPlace(ct1),
      StatusIs(absl::StatusCode::kInvalidArgument,
               HasSubstr("`that` is encrypted with a different key power")));
  EXPECT_THAT(
      ct0.SubInPlaceWithoutPad(ct1),
      StatusIs(absl::StatusCode::kInvalidArgument,
               HasSubstr("`that` is encrypted with a different key power")));
}

TYPED_TEST(RnsCiphertextTest, AbsorbScalarFailsIfModuliMismatch) {
  // This test requires at least two RNS moduli.
  if (this->moduli_.size() <= 1) {
    return;
  }

  // Create a ciphertext at the highest level.
  int curr_level = this->moduli_.size() - 1;
  ASSERT_OK_AND_ASSIGN(auto zero,
                       RnsPolynomial<TypeParam>::CreateZero(
                           this->rns_context_->LogN(), this->moduli_));
  RnsRlweCiphertext<TypeParam> ct(/*components=*/{zero, zero}, this->moduli_,
                                  /*power_of_s=*/1, /*error=*/0,
                                  this->error_params_.get());
  ASSERT_EQ(ct.Level(), curr_level);

  // Multiply with a scalar wrt few number of RNS moduli.
  auto scalar = this->rns_context_->PlaintextModulusInverseMainResidues();
  EXPECT_THAT(ct.AbsorbInPlace(scalar.subspan(0, curr_level)),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr(absl::StrCat("`scalar_mod_qs` must contain ",
                                              this->moduli_.size(),
                                              " modular integers"))));
}

TYPED_TEST(RnsCiphertextTest, FusedAbsorbAddInPlaceFailsIfDegreeMismatch) {
  ASSERT_OK_AND_ASSIGN(auto zero,
                       RnsPolynomial<TypeParam>::CreateZero(
                           this->rns_context_->LogN(), this->moduli_));
  RnsRlweCiphertext<TypeParam> ct_at_degree1(
      /*components=*/{zero, zero}, this->moduli_,
      /*power_of_s=*/1, /*error=*/0, this->error_params_.get());
  RnsRlweCiphertext<TypeParam> ct_at_degree2(
      /*components=*/{zero, zero, zero}, this->moduli_,
      /*power_of_s=*/1, /*error=*/0, this->error_params_.get());
  ASSERT_EQ(ct_at_degree1.Degree(), 1);
  ASSERT_EQ(ct_at_degree2.Degree(), 2);
  EXPECT_THAT(ct_at_degree1.FusedAbsorbAddInPlace(ct_at_degree2, zero),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`ctxt` has a mismatched degree")));
}

TYPED_TEST(RnsCiphertextTest, FusedAbsorbAddInPlaceFailsIfLevelMismatch) {
  // This test requires at least two RNS moduli.
  if (this->moduli_.size() <= 1) {
    return;
  }

  std::vector<const PrimeModulus<TypeParam>*> small_moduli = this->moduli_;
  small_moduli.pop_back();

  ASSERT_OK_AND_ASSIGN(auto zero,
                       RnsPolynomial<TypeParam>::CreateZero(
                           this->rns_context_->LogN(), this->moduli_));
  ASSERT_OK_AND_ASSIGN(auto zero_lower_level,
                       RnsPolynomial<TypeParam>::CreateZero(
                           this->rns_context_->LogN(), small_moduli));
  RnsRlweCiphertext<TypeParam> ct(
      /*components=*/{zero, zero}, this->moduli_,
      /*power_of_s=*/1, /*error=*/0, this->error_params_.get());
  RnsRlweCiphertext<TypeParam> ct_lower_level(
      /*components=*/{zero, zero}, small_moduli,
      /*power_of_s=*/1, /*error=*/0, this->error_params_.get());
  ASSERT_NE(ct.Level(), ct_lower_level.Level());
  EXPECT_THAT(ct.FusedAbsorbAddInPlace(ct_lower_level, zero),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`ctxt` has a mismatched level")));
}

TYPED_TEST(RnsCiphertextTest, FusedAbsorbAddInPlaceFailsIfPowerOfSMismatch) {
  ASSERT_OK_AND_ASSIGN(auto zero,
                       RnsPolynomial<TypeParam>::CreateZero(
                           this->rns_context_->LogN(), this->moduli_));
  RnsRlweCiphertext<TypeParam> ct0(
      /*components=*/{zero, zero}, this->moduli_,
      /*power_of_s=*/1, /*error=*/0, this->error_params_.get());
  RnsRlweCiphertext<TypeParam> ct1(
      /*components=*/{zero, zero}, this->moduli_,
      /*power_of_s=*/5, /*error=*/0, this->error_params_.get());
  ASSERT_EQ(ct0.PowerOfS(), 1);
  ASSERT_EQ(ct1.PowerOfS(), 5);
  EXPECT_THAT(ct0.FusedAbsorbAddInPlace(ct1, zero),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`ctxt` has a mismatched key power")));
}

TYPED_TEST(RnsCiphertextTest,
           FusedAbsorbAddInPlaceWithoutPadFailsIfNotDegreeOne) {
  ASSERT_OK_AND_ASSIGN(auto zero,
                       RnsPolynomial<TypeParam>::CreateZero(
                           this->rns_context_->LogN(), this->moduli_));
  RnsRlweCiphertext<TypeParam> ct_at_degree1(
      /*components=*/{zero, zero}, this->moduli_,
      /*power_of_s=*/1, /*error=*/0, this->error_params_.get());
  RnsRlweCiphertext<TypeParam> ct_at_degree2(
      /*components=*/{zero, zero, zero}, this->moduli_,
      /*power_of_s=*/1, /*error=*/0, this->error_params_.get());
  ASSERT_EQ(ct_at_degree1.Degree(), 1);
  ASSERT_EQ(ct_at_degree2.Degree(), 2);

  EXPECT_THAT(
      ct_at_degree1.FusedAbsorbAddInPlaceWithoutPad(ct_at_degree2, zero),
      StatusIs(absl::StatusCode::kInvalidArgument,
               HasSubstr("only applies to degree 1 ciphertext")));
  EXPECT_THAT(
      ct_at_degree2.FusedAbsorbAddInPlaceWithoutPad(ct_at_degree1, zero),
      StatusIs(absl::StatusCode::kInvalidArgument,
               HasSubstr("only applies to degree 1 ciphertext")));
}

TYPED_TEST(RnsCiphertextTest,
           FusedAbsorbAddInPlaceWithoutPadFailsIfLevelMismatch) {
  // This test requires at least two RNS moduli.
  if (this->moduli_.size() <= 1) {
    return;
  }

  std::vector<const PrimeModulus<TypeParam>*> small_moduli = this->moduli_;
  small_moduli.pop_back();

  ASSERT_OK_AND_ASSIGN(auto zero,
                       RnsPolynomial<TypeParam>::CreateZero(
                           this->rns_context_->LogN(), this->moduli_));
  ASSERT_OK_AND_ASSIGN(auto zero_lower_level,
                       RnsPolynomial<TypeParam>::CreateZero(
                           this->rns_context_->LogN(), small_moduli));
  RnsRlweCiphertext<TypeParam> ct(
      /*components=*/{zero, zero}, this->moduli_,
      /*power_of_s=*/1, /*error=*/0, this->error_params_.get());
  RnsRlweCiphertext<TypeParam> ct_lower_level(
      /*components=*/{zero, zero}, small_moduli,
      /*power_of_s=*/1, /*error=*/0, this->error_params_.get());
  ASSERT_NE(ct.Level(), ct_lower_level.Level());
  EXPECT_THAT(ct.FusedAbsorbAddInPlaceWithoutPad(ct_lower_level, zero),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`ctxt` has a mismatched level")));
}

TYPED_TEST(RnsCiphertextTest,
           FusedAbsorbAddInPlaceWithoutPadFailsIfPowerOfSMismatch) {
  ASSERT_OK_AND_ASSIGN(auto zero,
                       RnsPolynomial<TypeParam>::CreateZero(
                           this->rns_context_->LogN(), this->moduli_));
  RnsRlweCiphertext<TypeParam> ct0(
      /*components=*/{zero, zero}, this->moduli_,
      /*power_of_s=*/1, /*error=*/0, this->error_params_.get());
  RnsRlweCiphertext<TypeParam> ct1(
      /*components=*/{zero, zero}, this->moduli_,
      /*power_of_s=*/5, /*error=*/0, this->error_params_.get());
  ASSERT_EQ(ct0.PowerOfS(), 1);
  ASSERT_EQ(ct1.PowerOfS(), 5);
  EXPECT_THAT(ct0.FusedAbsorbAddInPlaceWithoutPad(ct1, zero),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`ctxt` has a mismatched key power")));
}

TYPED_TEST(RnsCiphertextTest, SetPadComponentFailsIfNotDegreeOne) {
  ASSERT_OK_AND_ASSIGN(auto zero,
                       RnsPolynomial<TypeParam>::CreateZero(
                           this->rns_context_->LogN(), this->moduli_));
  RnsRlweCiphertext<TypeParam> ct_at_degree2(
      /*components=*/{zero, zero, zero}, this->moduli_,
      /*power_of_s=*/1, /*error=*/0, this->error_params_.get());
  ASSERT_EQ(ct_at_degree2.Degree(), 2);

  EXPECT_THAT(ct_at_degree2.SetPadComponent(zero),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("only applies to degree 1 ciphertext")));
}

TYPED_TEST(RnsCiphertextTest, DeserializeFailsIfErrorParamsIsNull) {
  if (sizeof(TypeParam) > sizeof(Uint32)) {
    GTEST_SKIP() << "Skip negative test for large integer types.";
  }

  ASSERT_OK_AND_ASSIGN(auto zero,
                       RnsPolynomial<TypeParam>::CreateZero(
                           this->rns_context_->LogN(), this->moduli_));
  RnsRlweCiphertext<TypeParam> ciphertext(
      /*components=*/{zero, zero}, this->moduli_,
      /*power_of_s=*/1, /*error=*/0, this->error_params_.get());
  ASSERT_OK_AND_ASSIGN(SerializedRnsRlweCiphertext serialized,
                       ciphertext.Serialize());
}

TYPED_TEST(RnsCiphertextTest, SerializedDeserializes) {
  constexpr int k_power_of_s = 5;
  // Create a hypothetical ciphertext.
  ASSERT_OK_AND_ASSIGN(auto c0, RnsPolynomial<TypeParam>::SampleUniform(
                                    this->rns_context_->LogN(),
                                    this->prng_.get(), this->moduli_));
  ASSERT_OK_AND_ASSIGN(auto c1, RnsPolynomial<TypeParam>::SampleUniform(
                                    this->rns_context_->LogN(),
                                    this->prng_.get(), this->moduli_));
  RnsRlweCiphertext<TypeParam> ciphertext(
      /*components=*/{std::move(c0), std::move(c1)}, this->moduli_,
      k_power_of_s, this->error_params_->B_secretkey_encryption(),
      this->error_params_.get());

  ASSERT_OK_AND_ASSIGN(SerializedRnsRlweCiphertext serialized,
                       ciphertext.Serialize());
  ASSERT_OK_AND_ASSIGN(
      RnsRlweCiphertext<TypeParam> deserialized,
      RnsRlweCiphertext<TypeParam>::Deserialize(serialized, this->moduli_,
                                                this->error_params_.get()));
  ASSERT_EQ(deserialized.Degree(), ciphertext.Degree());
  ASSERT_EQ(deserialized.Level(), ciphertext.Level());
  ASSERT_EQ(deserialized.NumCoeffs(), ciphertext.NumCoeffs());
  ASSERT_EQ(deserialized.NumModuli(), ciphertext.NumModuli());
  ASSERT_EQ(deserialized.ErrorParams(), ciphertext.ErrorParams());
  ASSERT_EQ(deserialized.PowerOfS(), ciphertext.PowerOfS());
  ASSERT_EQ(deserialized.Error(), ciphertext.Error());
  for (int i = 0; i < ciphertext.Len(); ++i) {
    ASSERT_OK_AND_ASSIGN(auto deserialized_ci, deserialized.Component(i));
    ASSERT_OK_AND_ASSIGN(auto ciphertext_ci, ciphertext.Component(i));
    EXPECT_EQ(deserialized_ci, ciphertext_ci);
  }
}

}  // namespace
}  // namespace rlwe
