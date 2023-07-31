/*
 * Copyright 2022 Google LLC.
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

#include "shell_encryption/public_key_encryption.h"

#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include "absl/status/status.h"
#include "shell_encryption/constants.h"
#include "shell_encryption/context.h"
#include "shell_encryption/montgomery.h"
#include "shell_encryption/ntt_parameters.h"
#include "shell_encryption/polynomial.h"
#include "shell_encryption/serialization.pb.h"
#include "shell_encryption/status_macros.h"
#include "shell_encryption/symmetric_encryption.h"
#include "shell_encryption/testing/parameters.h"
#include "shell_encryption/testing/status_matchers.h"
#include "shell_encryption/testing/status_testing.h"
#include "shell_encryption/testing/testing_prng.h"
#include "shell_encryption/testing/testing_utils.h"

namespace {

using ::rlwe::testing::StatusIs;
using ::testing::HasSubstr;

// Set constants.
const int kTestingRounds = 10;

// Test fixture.
template <typename ModularInt>
class PublicKeyRlweEncryptionTest : public ::testing::Test {
 public:
  // Sample a secret key.
  rlwe::StatusOr<rlwe::SymmetricRlweKey<ModularInt>> SampleSecretKey(
      const rlwe::RlweContext<ModularInt>& context) {
    RLWE_ASSIGN_OR_RETURN(std::string prng_seed,
                          rlwe::SingleThreadHkdfPrng::GenerateSeed());
    RLWE_ASSIGN_OR_RETURN(auto prng,
                          rlwe::SingleThreadHkdfPrng::Create(prng_seed));
    return rlwe::SymmetricRlweKey<ModularInt>::Sample(
        context.GetLogN(), context.GetVariance(), context.GetLogT(),
        context.GetModulusParams(), context.GetNttParams(), prng.get());
  }

  // Generate a public key.
  rlwe::StatusOr<rlwe::PublicRlweKey<ModularInt>> GeneratePublicKey(
      const rlwe::SymmetricRlweKey<ModularInt>& secret_key,
      const rlwe::RlweContext<ModularInt>& context) {
    return rlwe::PublicRlweKey<ModularInt>::Create(
        secret_key, context.GetVariance(), rlwe::PRNG_TYPE_HKDF);
  }

  // Helper function to encrypt a plaintext (represented as its coefficient
  // vector).
  rlwe::StatusOr<rlwe::SymmetricRlweCiphertext<ModularInt>> Encrypt(
      const rlwe::PublicRlweKey<ModularInt>& public_key,
      const std::vector<typename ModularInt::Int>& coeffs,
      const rlwe::RlweContext<ModularInt>& context) {
    RLWE_ASSIGN_OR_RETURN(std::vector<ModularInt> coeffs_mod,
                          rlwe::testing::ConvertToMontgomery<ModularInt>(
                              coeffs, context.GetModulusParams()));
    rlwe::Polynomial<ModularInt> plaintext_ntt =
        rlwe::Polynomial<ModularInt>::ConvertToNtt(
            coeffs_mod, context.GetNttParams(), context.GetModulusParams());

    RLWE_ASSIGN_OR_RETURN(std::string prng_seed,
                          rlwe::SingleThreadHkdfPrng::GenerateSeed());
    RLWE_ASSIGN_OR_RETURN(auto prng,
                          rlwe::SingleThreadHkdfPrng::Create(prng_seed));
    return public_key.Encrypt(plaintext_ntt, context.GetVariance(),
                              *context.GetErrorParams(), prng.get());
  }
};
TYPED_TEST_SUITE(PublicKeyRlweEncryptionTest, rlwe::testing::ModularIntTypes);

// Ensure that a public key (b, a) is an symmetric encryption of 0
TYPED_TEST(PublicKeyRlweEncryptionTest, PublicKeyIsEncryptionOfZero) {
  for (const auto& params :
       rlwe::testing::ContextParameters<TypeParam>::Value()) {
    ASSERT_OK_AND_ASSIGN(auto context,
                         rlwe::RlweContext<TypeParam>::Create(params));

    ASSERT_OK_AND_ASSIGN(rlwe::SymmetricRlweKey<TypeParam> secret_key,
                         this->SampleSecretKey(*context));
    ASSERT_OK_AND_ASSIGN(rlwe::PublicRlweKey<TypeParam> public_key,
                         this->GeneratePublicKey(secret_key, *context));
    ASSERT_EQ(public_key.Len(), context->GetN());
    ASSERT_EQ(
        public_key.PlaintextModulus().ExportInt(&public_key.ModulusParams()),
        context->GetT());

    // Create a cipherext out of the polynomials a and b in public_key
    rlwe::SymmetricRlweCiphertext<TypeParam> ciphertext(
        {public_key.GetB(), public_key.GetA()}, /*power_of_s=*/1,
        context->GetErrorParams()->B_encryption(), context->GetModulusParams(),
        context->GetErrorParams());

    // Check if the public key is an encryption of 0
    ASSERT_OK_AND_ASSIGN(std::vector<typename TypeParam::Int> decrypted,
                         rlwe::Decrypt<TypeParam>(secret_key, ciphertext));
    ASSERT_EQ(decrypted.size(), context->GetN());
    for (size_t i = 0; i < decrypted.size(); ++i) {
      EXPECT_EQ(decrypted[i], 0);
    }
  }
}

// Ensure that public key encryption produces ciphertexts that can be correctly
// decrypted.
TYPED_TEST(PublicKeyRlweEncryptionTest, EncryptionIsCorrect) {
  for (const auto& params :
       rlwe::testing::ContextParameters<TypeParam>::Value()) {
    ASSERT_OK_AND_ASSIGN(auto context,
                         rlwe::RlweContext<TypeParam>::Create(params));
    for (unsigned int i = 0; i < kTestingRounds; i++) {
      ASSERT_OK_AND_ASSIGN(rlwe::SymmetricRlweKey<TypeParam> secret_key,
                           this->SampleSecretKey(*context));
      ASSERT_OK_AND_ASSIGN(rlwe::PublicRlweKey<TypeParam> public_key,
                           this->GeneratePublicKey(secret_key, *context));

      std::vector<typename TypeParam::Int> plaintext =
          rlwe::testing::SamplePlaintext<TypeParam>(context->GetN(),
                                                    context->GetT());
      // Encrypt using the public key and decrypt using the secret key.
      ASSERT_OK_AND_ASSIGN(rlwe::SymmetricRlweCiphertext<TypeParam> ciphertext,
                           this->Encrypt(public_key, plaintext, *context));
      ASSERT_OK_AND_ASSIGN(std::vector<typename TypeParam::Int> decrypted,
                           rlwe::Decrypt<TypeParam>(secret_key, ciphertext));
      // Check correctness.
      EXPECT_EQ(plaintext, decrypted);
    }
  }
}

// Ensure that multiple public keys can be derived from the same secret key,
// and that they can all be used to encrypt plaintexts.
TYPED_TEST(PublicKeyRlweEncryptionTest, MultiplePublicKeyInstances) {
  for (const auto& params :
       rlwe::testing::ContextParameters<TypeParam>::Value()) {
    ASSERT_OK_AND_ASSIGN(auto context,
                         rlwe::RlweContext<TypeParam>::Create(params));
    for (unsigned int i = 0; i < kTestingRounds; i++) {
      ASSERT_OK_AND_ASSIGN(rlwe::SymmetricRlweKey<TypeParam> secret_key,
                           this->SampleSecretKey(*context));
      // Generate two public keys from the same secret key
      ASSERT_OK_AND_ASSIGN(rlwe::PublicRlweKey<TypeParam> public_key0,
                           this->GeneratePublicKey(secret_key, *context));
      ASSERT_OK_AND_ASSIGN(rlwe::PublicRlweKey<TypeParam> public_key1,
                           this->GeneratePublicKey(secret_key, *context));
      // The two public keys should have different values
      rlwe::Polynomial<TypeParam> pk0_a = public_key0.GetA();
      rlwe::Polynomial<TypeParam> pk0_b = public_key0.GetB();
      rlwe::Polynomial<TypeParam> pk1_a = public_key1.GetA();
      rlwe::Polynomial<TypeParam> pk1_b = public_key1.GetB();
      EXPECT_NE(pk0_a.Coeffs(), pk1_a.Coeffs());
      EXPECT_NE(pk0_b.Coeffs(), pk1_b.Coeffs());

      // Encrypt the same plaintext using the public keys
      std::vector<typename TypeParam::Int> plaintext =
          rlwe::testing::SamplePlaintext<TypeParam>(context->GetN(),
                                                    context->GetT());
      ASSERT_OK_AND_ASSIGN(rlwe::SymmetricRlweCiphertext<TypeParam> ciphertext0,
                           this->Encrypt(public_key0, plaintext, *context));
      ASSERT_OK_AND_ASSIGN(rlwe::SymmetricRlweCiphertext<TypeParam> ciphertext1,
                           this->Encrypt(public_key1, plaintext, *context));
      // The two ciphertexts should be different
      ASSERT_OK_AND_ASSIGN(rlwe::Polynomial<TypeParam> ct0_b,
                           ciphertext0.Component(0));
      ASSERT_OK_AND_ASSIGN(rlwe::Polynomial<TypeParam> ct0_a,
                           ciphertext0.Component(1));
      ASSERT_OK_AND_ASSIGN(rlwe::Polynomial<TypeParam> ct1_b,
                           ciphertext1.Component(0));
      ASSERT_OK_AND_ASSIGN(rlwe::Polynomial<TypeParam> ct1_a,
                           ciphertext1.Component(1));
      EXPECT_NE(ct0_b.Coeffs(), ct1_b.Coeffs());
      EXPECT_NE(ct0_a.Coeffs(), ct1_a.Coeffs());

      // Decryption of the two ciphertexts should both be correct
      ASSERT_OK_AND_ASSIGN(std::vector<typename TypeParam::Int> decrypted0,
                           rlwe::Decrypt<TypeParam>(secret_key, ciphertext0));
      ASSERT_OK_AND_ASSIGN(std::vector<typename TypeParam::Int> decrypted1,
                           rlwe::Decrypt<TypeParam>(secret_key, ciphertext1));
      EXPECT_EQ(decrypted0, plaintext);
      EXPECT_EQ(decrypted1, plaintext);
    }
  }
}

// Check that public-key serialization works.
TYPED_TEST(PublicKeyRlweEncryptionTest, SerializePublicKey) {
  for (const auto& params :
       rlwe::testing::ContextParameters<TypeParam>::Value()) {
    ASSERT_OK_AND_ASSIGN(auto context,
                         rlwe::RlweContext<TypeParam>::Create(params));
    for (int i = 0; i < kTestingRounds; i++) {
      ASSERT_OK_AND_ASSIGN(rlwe::SymmetricRlweKey<TypeParam> secret_key,
                           this->SampleSecretKey(*context));
      ASSERT_OK_AND_ASSIGN(rlwe::PublicRlweKey<TypeParam> public_key,
                           this->GeneratePublicKey(secret_key, *context));

      // Serialize.
      ASSERT_OK_AND_ASSIGN(rlwe::SerializedPublicRlweKey serialized,
                           public_key.Serialize());

      // Get plaintext modulus.
      ASSERT_OK_AND_ASSIGN(
          TypeParam t_mod,
          TypeParam::ImportInt(context->GetT(), context->GetModulusParams()));
      // Deserialize.
      ASSERT_OK_AND_ASSIGN(rlwe::PublicRlweKey<TypeParam> deserialized,
                           rlwe::PublicRlweKey<TypeParam>::Deserialize(
                               serialized, t_mod, context->GetModulusParams(),
                               context->GetNttParams()));

      // Check equality.
      ASSERT_EQ(public_key.Len(), deserialized.Len());
      ASSERT_EQ(public_key.PlaintextModulus(), deserialized.PlaintextModulus());
      rlwe::Polynomial<TypeParam> deserialized_a = deserialized.GetA();
      rlwe::Polynomial<TypeParam> deserialized_b = deserialized.GetB();
      rlwe::Polynomial<TypeParam> publickey_a = public_key.GetA();
      rlwe::Polynomial<TypeParam> publickey_b = public_key.GetB();
      for (int j = 0; j < deserialized.Len(); ++j) {
        // The two components in the public key must be equal
        EXPECT_EQ(publickey_a.Coeffs()[j], deserialized_a.Coeffs()[j]);
        EXPECT_EQ(publickey_b.Coeffs()[j], deserialized_b.Coeffs()[j]);
      }

      // Test that a ciphertext encrypted with the serialized key decrypts under
      // the original secret key.
      std::vector<typename TypeParam::Int> plaintext =
          rlwe::testing::SamplePlaintext<TypeParam>(context->GetN(),
                                                    context->GetT());
      ASSERT_OK_AND_ASSIGN(rlwe::SymmetricRlweCiphertext<TypeParam> ciphertext,
                           this->Encrypt(deserialized, plaintext, *context));
      ASSERT_OK_AND_ASSIGN(std::vector<typename TypeParam::Int> decrypted,
                           rlwe::Decrypt<TypeParam>(secret_key, ciphertext));
      EXPECT_EQ(decrypted, plaintext);
    }
  }
}

// Check for error handling in generating a public key.
TYPED_TEST(PublicKeyRlweEncryptionTest, InvalidCreation) {
  for (const auto& params :
       rlwe::testing::ContextParameters<TypeParam>::Value()) {
    ASSERT_OK_AND_ASSIGN(auto context,
                         rlwe::RlweContext<TypeParam>::Create(params));
    ASSERT_OK_AND_ASSIGN(rlwe::SymmetricRlweKey<TypeParam> secret_key,
                         this->SampleSecretKey(*context));
    // Generate a public key using an invalid PRNG type
    EXPECT_THAT(
        rlwe::PublicRlweKey<TypeParam>::Create(
            secret_key, context->GetVariance(), rlwe::PRNG_TYPE_INVALID),
        StatusIs(::absl::StatusCode::kInvalidArgument,
                 HasSubstr("PrngType not specified correctly.")));
  }
}

// Check for error handling in invalid deserialization.
TYPED_TEST(PublicKeyRlweEncryptionTest, InvalidDeserialization) {
  for (const auto& params :
       rlwe::testing::ContextParameters<TypeParam>::Value()) {
    ASSERT_OK_AND_ASSIGN(auto context,
                         rlwe::RlweContext<TypeParam>::Create(params));
    ASSERT_OK_AND_ASSIGN(rlwe::SymmetricRlweKey<TypeParam> secret_key,
                         this->SampleSecretKey(*context));
    ASSERT_OK_AND_ASSIGN(rlwe::PublicRlweKey<TypeParam> public_key,
                         this->GeneratePublicKey(secret_key, *context));

    // Get plaintext modulus.
    ASSERT_OK_AND_ASSIGN(
        TypeParam t_mod,
        TypeParam::ImportInt(context->GetT(), context->GetModulusParams()));
    // Serialize.
    ASSERT_OK_AND_ASSIGN(rlwe::SerializedPublicRlweKey serialized,
                         public_key.Serialize());
    // Tamper the serialization with an invalid value.
    serialized.set_prng_type(rlwe::PRNG_TYPE_INVALID);
    // Deserialize.
    EXPECT_THAT(rlwe::PublicRlweKey<TypeParam>::Deserialize(
                    serialized, t_mod, context->GetModulusParams(),
                    context->GetNttParams()),
                StatusIs(::absl::StatusCode::kInvalidArgument,
                         HasSubstr("Invalid PRNG type is specified.")));
  }
}

// Check that public-key encryption generates different ciphertexts when
// different plaintexts are encrypted.
TYPED_TEST(PublicKeyRlweEncryptionTest,
           DistinctCiphertextsForDistinctPlaintexts) {
  for (const auto& params :
       rlwe::testing::ContextParameters<TypeParam>::Value()) {
    ASSERT_OK_AND_ASSIGN(auto context,
                         rlwe::RlweContext<TypeParam>::Create(params));
    ASSERT_OK_AND_ASSIGN(rlwe::SymmetricRlweKey<TypeParam> secret_key,
                         this->SampleSecretKey(*context));
    ASSERT_OK_AND_ASSIGN(rlwe::PublicRlweKey<TypeParam> public_key,
                         this->GeneratePublicKey(secret_key, *context));

    // Encrypt two different plaintext messages. Binary values are used here
    // to be compatible with all possible plaintext moduli.
    std::vector<typename TypeParam::Int> plaintext0{0, 1, 0, 1, 1};
    std::vector<typename TypeParam::Int> plaintext1{0, 0, 1, 1, 0};
    plaintext0.resize(context->GetN(), 0);
    plaintext1.resize(context->GetN(), 0);
    ASSERT_OK_AND_ASSIGN(rlwe::SymmetricRlweCiphertext<TypeParam> ciphertext0,
                         this->Encrypt(public_key, plaintext0, *context));
    ASSERT_OK_AND_ASSIGN(rlwe::SymmetricRlweCiphertext<TypeParam> ciphertext1,
                         this->Encrypt(public_key, plaintext1, *context));

    // Sanity check: the two ciphertexts should decrypt correctly
    ASSERT_OK_AND_ASSIGN(std::vector<typename TypeParam::Int> decrypted0,
                         rlwe::Decrypt<TypeParam>(secret_key, ciphertext0));
    EXPECT_EQ(decrypted0, plaintext0);
    ASSERT_OK_AND_ASSIGN(std::vector<typename TypeParam::Int> decrypted1,
                         rlwe::Decrypt<TypeParam>(secret_key, ciphertext1));
    EXPECT_EQ(decrypted1, plaintext1);

    // Check the polynomial components of the ciphertexts are different
    ASSERT_OK_AND_ASSIGN(rlwe::Polynomial<TypeParam> ct0_b,
                         ciphertext0.Component(0));
    ASSERT_OK_AND_ASSIGN(rlwe::Polynomial<TypeParam> ct0_a,
                         ciphertext0.Component(1));
    ASSERT_OK_AND_ASSIGN(rlwe::Polynomial<TypeParam> ct1_b,
                         ciphertext1.Component(0));
    ASSERT_OK_AND_ASSIGN(rlwe::Polynomial<TypeParam> ct1_a,
                         ciphertext1.Component(1));
    EXPECT_NE(ct0_b.Coeffs(), ct1_b.Coeffs());
    EXPECT_NE(ct0_a.Coeffs(), ct1_a.Coeffs());
  }
}

// Check that public-key encryption generates different ciphertexts when
// different randomness are encrypted.
TYPED_TEST(PublicKeyRlweEncryptionTest,
           DistinctCiphertextsForDistinctRandomness) {
  for (const auto& params :
       rlwe::testing::ContextParameters<TypeParam>::Value()) {
    ASSERT_OK_AND_ASSIGN(auto context,
                         rlwe::RlweContext<TypeParam>::Create(params));
    ASSERT_OK_AND_ASSIGN(rlwe::SymmetricRlweKey<TypeParam> secret_key,
                         this->SampleSecretKey(*context));
    ASSERT_OK_AND_ASSIGN(rlwe::PublicRlweKey<TypeParam> public_key,
                         this->GeneratePublicKey(secret_key, *context));

    // Encrypt the same plaintext twice using different PRNG seeds
    std::vector<typename TypeParam::Int> plaintext{0, 1, 0, 1, 1};
    plaintext.resize(context->GetN(), 0);
    ASSERT_OK_AND_ASSIGN(std::vector<TypeParam> plaintext_mod,
                         rlwe::testing::ConvertToMontgomery<TypeParam>(
                             plaintext, context->GetModulusParams()));
    rlwe::Polynomial<TypeParam> plaintext_ntt =
        rlwe::Polynomial<TypeParam>::ConvertToNtt(plaintext_mod,
                                                  context->GetNttParams(),
                                                  context->GetModulusParams());

    ASSERT_OK_AND_ASSIGN(std::string prng_seed0,
                         rlwe::SingleThreadHkdfPrng::GenerateSeed());
    ASSERT_OK_AND_ASSIGN(std::string prng_seed1,
                         rlwe::SingleThreadHkdfPrng::GenerateSeed());
    // Sanity check: We should have different PRNG seeds.
    ASSERT_NE(prng_seed0, prng_seed1);

    ASSERT_OK_AND_ASSIGN(auto prng0,
                         rlwe::SingleThreadHkdfPrng::Create(prng_seed0));
    ASSERT_OK_AND_ASSIGN(auto prng1,
                         rlwe::SingleThreadHkdfPrng::Create(prng_seed1));
    ASSERT_OK_AND_ASSIGN(
        rlwe::SymmetricRlweCiphertext<TypeParam> ciphertext0,
        public_key.Encrypt(plaintext_ntt, context->GetVariance(),
                           *context->GetErrorParams(), prng0.get()));
    ASSERT_OK_AND_ASSIGN(
        rlwe::SymmetricRlweCiphertext<TypeParam> ciphertext1,
        public_key.Encrypt(plaintext_ntt, context->GetVariance(),
                           *context->GetErrorParams(), prng1.get()));

    // Check the polynomial components of the ciphertexts are different
    ASSERT_OK_AND_ASSIGN(rlwe::Polynomial<TypeParam> ct0_b,
                         ciphertext0.Component(0));
    ASSERT_OK_AND_ASSIGN(rlwe::Polynomial<TypeParam> ct0_a,
                         ciphertext0.Component(1));
    ASSERT_OK_AND_ASSIGN(rlwe::Polynomial<TypeParam> ct1_b,
                         ciphertext1.Component(0));
    ASSERT_OK_AND_ASSIGN(rlwe::Polynomial<TypeParam> ct1_a,
                         ciphertext1.Component(1));
    EXPECT_NE(ct0_b.Coeffs(), ct1_b.Coeffs());
    EXPECT_NE(ct0_a.Coeffs(), ct1_a.Coeffs());
  }
}

// Check that using a wrong secret key for decryption does not output a
// correct plaintext
TYPED_TEST(PublicKeyRlweEncryptionTest, IncorrectSecretKeyInDecryption) {
  for (const auto& params :
       rlwe::testing::ContextParameters<TypeParam>::Value()) {
    ASSERT_OK_AND_ASSIGN(auto context,
                         rlwe::RlweContext<TypeParam>::Create(params));
    // Sample two secret keys
    ASSERT_OK_AND_ASSIGN(rlwe::SymmetricRlweKey<TypeParam> secret_key0,
                         this->SampleSecretKey(*context));
    ASSERT_OK_AND_ASSIGN(rlwe::SymmetricRlweKey<TypeParam> secret_key1,
                         this->SampleSecretKey(*context));
    // Generate a public key from the first secret key
    ASSERT_OK_AND_ASSIGN(rlwe::PublicRlweKey<TypeParam> public_key,
                         this->GeneratePublicKey(secret_key0, *context));

    // Encrypt the plaintext using the public key
    std::vector<typename TypeParam::Int> plaintext =
        rlwe::testing::SamplePlaintext<TypeParam>(context->GetN(),
                                                  context->GetT());
    ASSERT_OK_AND_ASSIGN(rlwe::SymmetricRlweCiphertext<TypeParam> ciphertext,
                         this->Encrypt(public_key, plaintext, *context));

    // Sanity check: Decrypt using the correct secret key
    ASSERT_OK_AND_ASSIGN(std::vector<typename TypeParam::Int> correct_decrypted,
                         rlwe::Decrypt<TypeParam>(secret_key0, ciphertext));
    EXPECT_EQ(correct_decrypted, plaintext);

    // Decrypt using the other secret key, and expect incorrect output
    ASSERT_OK_AND_ASSIGN(std::vector<typename TypeParam::Int> wrong_decrypted,
                         rlwe::Decrypt<TypeParam>(secret_key1, ciphertext));
    EXPECT_NE(wrong_decrypted, plaintext);
  }
}

}  // namespace
