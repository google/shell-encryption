/*
 * Copyright 2017 Google Inc.
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

#include <cstdint>
#include <iostream>
#include <random>
#include <vector>

#include <gtest/gtest.h>
#include "constants.h"
#include "montgomery.h"
#include "ntt_parameters.h"
#include "ntt_polynomial.h"
#include "serialization.pb.h"
#include "symmetric_encryption.h"
#include "testing/utils.h"

namespace {

// Seed for the random number generator that is used to create test plaintexts.
unsigned int seed = 1;

// Set constants.
const uint64_t kDefaultLogT = 2;
const uint64_t kDefaultT = 1 << kDefaultLogT;
const uint64_t kDefaultVariance = 8;
const uint64_t kCoeffs = rlwe::kNewhopeDegreeBound;
const uint64_t kLogCoeffs = rlwe::kNewhopeLogDegreeBound;
const uint64_t kTestingRounds = 200;

// Useful typedefs.
using uint_m = rlwe::MontgomeryInt;
rlwe::MontgomeryIntParams params14(rlwe::kNewhopeLogR, rlwe::kNewhopeModulus);
using NttPolynomial = rlwe::NttPolynomial<uint_m>;
using Ciphertext = rlwe::SymmetricRlweCiphertext<uint_m>;
using Key = rlwe::SymmetricRlweKey<uint_m>;

// Test fixture.
class SymmetricRlweEncryptionTest : public ::testing::Test {
 protected:
  // Convert a vector of integers to a vector of montgomery integers.
  static std::vector<uint_m> ConvertToMontgomery(
      const std::vector<uint_m::Int>& coeffs) {
    std::vector<uint_m> output(coeffs.size(), uint_m::ImportInt(&params14, 0));
    for (unsigned int i = 0; i < output.size(); i++) {
      output[i] = uint_m::ImportInt(&params14, coeffs[i]);
    }
    return output;
  }

  // Sample a random key.
  Key SampleKey(uint_m::Int variance = kDefaultVariance,
                uint_m::Int log_t = kDefaultLogT) {
    return Key::Sample(kLogCoeffs, variance, log_t, prng_.get(), &params14,
                       ntt_params_);
  }

  // Sample a random plaintext.
  static std::vector<uint_m::Int> SamplePlaintext(uint_m::Int t = kDefaultT) {
    std::vector<uint_m::Int> plaintext(kCoeffs);
    for (unsigned int i = 0; i < kCoeffs; i++) {
      plaintext[i] = rand_r(&seed) % t;
    }

    return plaintext;
  }

  // Encrypt a plaintext.
  Ciphertext Encrypt(const Key& key,
                     const std::vector<uint_m::Int>& plaintext) const {
    auto plaintext_ntt = NttPolynomial::ConvertToNtt(
        ConvertToMontgomery(plaintext), ntt_params_);
    return rlwe::Encrypt<uint_m>(key, plaintext_ntt);
  }

  // Set up global variables.
  void SetUp() override {
    prng_.reset(new rlwe::testing::TestingPrng(1001));  // NOLINT

    // Acquire all of the NTT parameters.
    ntt_params_ = rlwe::InitializeNttParameters<uint_m>(&params14, kLogCoeffs);
  }

  rlwe::NttParameters<uint_m> ntt_params_;
  std::unique_ptr<rlwe::testing::TestingPrng> prng_;
};

// Ensure that RemoveError works correctly on negative numbers for several
// different values of t.
TEST_F(SymmetricRlweEncryptionTest, RemoveErrorNegative) {
  for (int t = 2; t < 16; t++) {
    for (uint64_t i = rlwe::kNewhopeModulus / 2 + 1; i < rlwe::kNewhopeModulus;
         i++) {
      // Create a vector that exclusively contains the value i.
      std::vector<uint_m> error_and_message(kCoeffs,
                                            uint_m::ImportInt(&params14, i));
      auto result =
          rlwe::RemoveError<uint_m>(error_and_message, params14.modulus, t);

      // Compute the expected result using signed arithmetic. Derive its
      // negative equivalent by subtracting out kNewhopeModulus and taking
      // that negative value (mod t).
      int64_t expected = (static_cast<int64_t>(i) -
                          static_cast<int64_t>(rlwe::kNewhopeModulus)) %
                         t;

      // Finally, turn any negative values into their positive equivalents
      // (mod t).
      if (expected < 0) {
        expected += t;
      }

      for (unsigned int j = 0; j < kCoeffs; j++) {
        EXPECT_EQ(expected, result[j]) << t << i;
      }
    }
  }
}

// Ensure that RemoveError works correctly on negative numbers for several
// different values of t.
TEST_F(SymmetricRlweEncryptionTest, RemoveErrorPositive) {
  for (int t = 2; t < 16; t++) {
    for (unsigned int i = 0; i <= rlwe::kNewhopeModulus / 2; i++) {
      // Create a vector that exclusively contains the value i.
      std::vector<uint_m> error_and_message(kCoeffs,
                                            uint_m::ImportInt(&params14, i));
      auto result =
          rlwe::RemoveError<uint_m>(error_and_message, params14.modulus, t);

      for (unsigned int j = 0; j < kCoeffs; j++) {
        EXPECT_EQ(i % t, result[j]);
      }
    }
  }
}

// Ensure that the encryption scheme can decrypt its own ciphertexts.
TEST_F(SymmetricRlweEncryptionTest, CanDecrypt) {
  for (unsigned int i = 0; i < kTestingRounds; i++) {
    auto key = SampleKey();
    auto plaintext = SamplePlaintext();
    auto ciphertext = Encrypt(key, plaintext);
    auto decrypted = rlwe::Decrypt<uint_m>(key, ciphertext);

    EXPECT_EQ(plaintext, decrypted);
  }
}

TEST_F(SymmetricRlweEncryptionTest, AdditivelyHomomorphic) {
  for (unsigned int i = 0; i < kTestingRounds; i++) {
    auto key = SampleKey();

    std::vector<uint_m::Int> plaintext1 = SamplePlaintext();
    std::vector<uint_m::Int> plaintext2 = SamplePlaintext();

    auto ciphertext1 = Encrypt(key, plaintext1);
    auto ciphertext2 = Encrypt(key, plaintext2);
    auto ciphertext3 = ciphertext1 + ciphertext2;

    std::vector<uint_m::Int> decrypted =
        rlwe::Decrypt<uint_m>(key, ciphertext3);

    for (unsigned int j = 0; j < plaintext1.size(); j++) {
      EXPECT_EQ((plaintext1[j] + plaintext2[j]) % kDefaultT, decrypted[j]);
    }
  }
}

TEST_F(SymmetricRlweEncryptionTest, AddToZero) {
  for (unsigned int i = 0; i < kTestingRounds; i++) {
    auto key = SampleKey();

    std::vector<uint_m::Int> plaintext = SamplePlaintext();

    Ciphertext ciphertext1(&params14);
    auto ciphertext2 = Encrypt(key, plaintext);
    auto ciphertext3 = ciphertext1 + ciphertext2;

    std::vector<uint_m::Int> decrypted =
        rlwe::Decrypt<uint_m>(key, ciphertext3);

    EXPECT_EQ(plaintext, decrypted);
  }
}

TEST_F(SymmetricRlweEncryptionTest, Absorb) {
  for (unsigned int i = 0; i < kTestingRounds; i++) {
    auto key = SampleKey();

    // Create the initial plaintexts.
    std::vector<uint_m::Int> plaintext = SamplePlaintext();
    NttPolynomial plaintext_ntt = NttPolynomial::ConvertToNtt(
        ConvertToMontgomery(plaintext), ntt_params_);
    std::vector<uint_m::Int> to_absorb = SamplePlaintext();
    NttPolynomial to_absorb_ntt = NttPolynomial::ConvertToNtt(
        ConvertToMontgomery(to_absorb), ntt_params_);

    // Create our expected value.
    NttPolynomial expected_ntt = plaintext_ntt * to_absorb_ntt;
    std::vector<uint_m::Int> expected = rlwe::RemoveError<uint_m>(
        expected_ntt.InverseNtt(ntt_params_), params14.modulus, kDefaultT);

    // Encrypt, absorb, and decrypt.
    auto ciphertext = Encrypt(key, plaintext) * to_absorb_ntt;
    std::vector<uint_m::Int> decrypted = rlwe::Decrypt<uint_m>(key, ciphertext);

    EXPECT_EQ(expected, decrypted);
  }
}

TEST_F(SymmetricRlweEncryptionTest, MultiplicativelyHomomorphic) {
  for (int i = 0; i < 200; i++) {
    // Since error builds up *very* quickly for homomorphic encryption, use
    // alternate values of t and variance.
    uint_m::Int log_t = 1;
    uint_m::Int variance = 4;

    auto key = SampleKey(variance, log_t);

    // Create the initial plaintexts.
    std::vector<uint_m::Int> plaintext1 = SamplePlaintext(1 << log_t);
    NttPolynomial plaintext1_ntt = NttPolynomial::ConvertToNtt(
        ConvertToMontgomery(plaintext1), ntt_params_);
    std::vector<uint_m::Int> plaintext2 = SamplePlaintext(1 << log_t);
    NttPolynomial plaintext2_ntt = NttPolynomial::ConvertToNtt(
        ConvertToMontgomery(plaintext2), ntt_params_);

    // Encrypt, multiply, and decrypt.
    auto ciphertext1 = Encrypt(key, plaintext1);
    auto ciphertext2 = Encrypt(key, plaintext2);
    std::vector<uint_m::Int> decrypted =
        rlwe::Decrypt<uint_m>(key, ciphertext1 * ciphertext2);

    // Create the polynomial we expect.
    NttPolynomial expected_ntt = plaintext1_ntt * plaintext2_ntt;
    std::vector<uint_m::Int> expected = rlwe::RemoveError<uint_m>(
        expected_ntt.InverseNtt(ntt_params_), params14.modulus, 1 << log_t);

    EXPECT_EQ(expected, decrypted);
  }
}

TEST_F(SymmetricRlweEncryptionTest, ManyHomomorphicAdds) {
  for (int i = 0; i < 20; i++) {
    // Sample a starting plaintext and ciphertext and create aggregators;
    std::vector<uint_m::Int> plaintext = SamplePlaintext();
    std::vector<uint_m::Int> plaintext_sum = plaintext;
    Key key = SampleKey();
    auto ciphertext_sum = Encrypt(key, plaintext);

    for (int j = 0; j < 50; j++) {
      // Sample a fresh plaintext.
      plaintext = SamplePlaintext();
      auto ciphertext = Encrypt(key, plaintext);

      // Add the new plaintext to the old plaintext.
      for (unsigned int k = 0; k < kCoeffs; k++) {
        plaintext_sum[k] += plaintext[k];
        plaintext_sum[k] %= kDefaultT;
      }

      // Add the new ciphertext to the old ciphertext.
      ciphertext_sum = ciphertext_sum + ciphertext;
      std::vector<uint_m::Int> decrypted =
          rlwe::Decrypt<uint_m>(key, ciphertext_sum);

      // Ensure the values are the same.
      EXPECT_EQ(plaintext_sum, decrypted);
    }
  }
}

TEST_F(SymmetricRlweEncryptionTest, Serialization) {
  for (int i = 0; i < 200; i++) {
    Key key = SampleKey();
    std::vector<uint_m::Int> plaintext = SamplePlaintext();
    auto ciphertext = Encrypt(key, plaintext);

    // Serialize and deserialize.
    rlwe::SerializedSymmetricRlweCiphertext serialized = ciphertext.Serialize();
    auto deserialized = Ciphertext::Deserialize(&params14, serialized);

    // Decrypt and check equality.
    std::vector<uint_m::Int> deserialized_plaintext =
        rlwe::Decrypt<uint_m>(key, deserialized);

    EXPECT_EQ(plaintext, deserialized_plaintext);
  }
}

// Test modulus switching.
TEST_F(SymmetricRlweEncryptionTest, ModulusReduction) {
  // Constants for this test.
  int log_num_coeffs = 10;
  int num_coeffs = 1 << log_num_coeffs;
  unsigned int seed = 0;

  // 30-bit modulus.
  rlwe::MontgomeryIntParams params30(rlwe::kLogR30, rlwe::kModulus30);
  auto ntt_params30 =
      rlwe::InitializeNttParameters<uint_m>(&params30, log_num_coeffs);

  // 14-bit modulus.
  auto ntt_params14 =
      rlwe::InitializeNttParameters<uint_m>(&params14, log_num_coeffs);

  for (int trials = 0; trials < 1000; trials++) {
    // Create a key.
    auto key30 = rlwe::SymmetricRlweKey<uint_m>::Sample(
        log_num_coeffs, kDefaultVariance, kDefaultLogT, prng_.get(), &params30,
        ntt_params30);
    auto key14 = key30.SwitchModulus<uint_m>(&params14, ntt_params14);

    // Create a plaintext.
    std::vector<uint_m::Int> plaintext(num_coeffs);
    std::vector<uint_m> plaintext_montgomery;
    for (uint_m::Int& coeff : plaintext) {
      coeff = rand_r(&seed) % kDefaultT;
      plaintext_montgomery.push_back(uint_m::ImportInt(&params30, coeff));
    }

    // Encrypt.
    auto plaintext_ntt = rlwe::NttPolynomial<uint_m>::ConvertToNtt(
        plaintext_montgomery, ntt_params30);
    auto ciphertext30 = rlwe::Encrypt<uint_m>(key30, plaintext_ntt);

    // Switch moduli.
    auto ciphertext14 = ciphertext30.SwitchModulus<uint_m>(
        ntt_params30, &params14, ntt_params14, kDefaultT);

    // Decrypt in the smaller modulus.
    auto decrypted = rlwe::Decrypt<uint_m>(key14, ciphertext14);

    EXPECT_EQ(plaintext, decrypted);
  }
}

}  // namespace
