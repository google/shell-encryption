// Copyright 2021 Google LLC
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

#include "oblivious_expand.h"

#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <random>
#include <vector>

#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include "absl/memory/memory.h"
#include "absl/status/status.h"
#include "constants.h"
#include "galois_key.h"
#include "montgomery.h"
#include "ntt_parameters.h"
#include "symmetric_encryption.h"
#include "testing/status_matchers.h"
#include "testing/status_testing.h"
#include "testing/testing_utils.h"

namespace rlwe {
namespace {
using ::rlwe::testing::StatusIs;
using ::testing::HasSubstr;

// Set constants. Set small enough that the test should succeed with very high
// probability.
const Uint64 kLogT = 3;
const Uint64 kVariance = 4;

// 128 coefficients is insecure, but it's small enough that ObliviousExpand
// works relatively fast.
const Uint64 kLogCoeffs = 7;
const Uint64 kCoeffs = 1 << kLogCoeffs;

const Uint64 kTotalSize =
    kCoeffs + 3;  // Test vectors larger than the number of coefficients.
const Uint64 kNumberIndices = kTotalSize / 2;  // The number of indices set to 1

unsigned int seed = 637;
std::mt19937 mt_rand(seed);

// Useful typedefs.
using uint_m = MontgomeryInt<Uint64>;
using Key = SymmetricRlweKey<uint_m>;
using Ciphertext = SymmetricRlweCiphertext<uint_m>;

class ObliviousExpandTest : public ::testing::TestWithParam<PrngType> {
 protected:
  void SetUp() override {
    prng_type_ = GetParam();

    ASSERT_OK_AND_ASSIGN(params59_, uint_m::Params::Create(kModulus59));
    ASSERT_OK_AND_ASSIGN(auto ntt_params, InitializeNttParameters<uint_m>(
                                              kLogCoeffs, params59_.get()));
    ntt_params_ =
        absl::make_unique<const NttParameters<uint_m>>(std::move(ntt_params));
    ASSERT_OK_AND_ASSIGN(
        auto error_params,
        ErrorParams<uint_m>::Create(kLogT, kVariance, params59_.get(),
                                    ntt_params_.get()));
    error_params_ = absl::make_unique<const ErrorParams<uint_m>>(error_params);
    ASSERT_OK_AND_ASSIGN(std::string prng_seed, GenerateSeed());
    ASSERT_OK_AND_ASSIGN(auto prng, CreatePrng(prng_seed));
    ASSERT_OK_AND_ASSIGN(
        auto key, Key::Sample(kLogCoeffs, kVariance, kLogT, params59_.get(),
                              ntt_params_.get(), prng.get()));
    key_ = absl::make_unique<Key>(std::move(key));
    for (int i = 0; i < kLogCoeffs; i++) {
      ASSERT_OK_AND_ASSIGN(
          auto galois_key,
          GaloisKey<uint_m>::Create(*key_, prng_type_, (kCoeffs >> i) + 1, 10));
      galois_keys_.push_back(galois_key);
    }
  }

  rlwe::StatusOr<std::string> GenerateSeed() {
    return rlwe::testing::GenerateSeed(prng_type_);
  }

  rlwe::StatusOr<std::unique_ptr<rlwe::SecurePrng>> CreatePrng(
      absl::string_view seed) {
    return rlwe::testing::CreatePrng(seed, prng_type_);
  }

  // Produces a vector of ciphertext of length (total_size / compression_factor)
  // such that it compresses, normalizes, and encrypts and the binary vector
  // with 1 at the positions specified by indices, and 0 elsewhere.
  StatusOr<std::vector<Ciphertext>> CompressAndEncryptVector(
      int total_size, const std::vector<int>& indices,
      int log_compression_factor) {
    std::vector<Ciphertext> result;

    // Create the compressed plaintexts.
    RLWE_ASSIGN_OR_RETURN(std::vector<Polynomial<uint_m>> plaintexts,
                          MakeCompressedVector<uint_m>(
                              total_size, indices, log_compression_factor,
                              params59_.get(), ntt_params_.get()));

    // When expanding the vector of ciphertexts, the coefficients will be
    // multiplied by 2^levels_of_recursion. Since the plaintext space cannot be
    // a power of two in the ObliviousExpand protocol, we pre-multiply the
    // plaintexts by 2^-log_compression_factor modulo the plaintext space.
    auto normalizer = ObliviousExpander<uint_m>::ComputeNormalizer(
        log_compression_factor, kLogT);
    RLWE_ASSIGN_OR_RETURN(auto m_normalizer,
                          uint_m::ImportInt(normalizer, key_->ModulusParams()));
    for (auto& plaintext : plaintexts) {
      RLWE_RETURN_IF_ERROR(
          plaintext.MulInPlace(m_normalizer, key_->ModulusParams()));
    }

    RLWE_ASSIGN_OR_RETURN(std::string prng_seed, GenerateSeed());
    RLWE_ASSIGN_OR_RETURN(auto prng, CreatePrng(prng_seed));

    // Encrypt entries of plaintexts.
    for (const auto& plaintext : plaintexts) {
      RLWE_ASSIGN_OR_RETURN(
          auto encrypt,
          Encrypt(*key_, plaintext, error_params_.get(), prng.get()));
      result.push_back(encrypt);
    }
    return result;
  }

  // Sample a set of number_ones distinct indices in [0, kTotalSize-1].
  std::vector<int> SampleRandomIndices(int number_ones = kNumberIndices) {
    std::vector<int> indices;
    indices.reserve(number_ones);
    while (indices.size() != number_ones) {
      int index = mt_rand() % kTotalSize;
      if (std::count(indices.begin(), indices.end(), index) == 0) {
        indices.push_back(index);
      }
    }
    return indices;
  }

  std::unique_ptr<const uint_m::Params> params59_;
  std::unique_ptr<const NttParameters<uint_m>> ntt_params_;
  std::unique_ptr<const ErrorParams<uint_m>> error_params_;
  std::unique_ptr<Key> key_;
  std::vector<GaloisKey<uint_m>> galois_keys_;

  PrngType prng_type_;
};

TEST_P(ObliviousExpandTest, CompressedVectorHasCorrectSizeAndEntries) {
  auto indices = SampleRandomIndices();
  ASSERT_OK_AND_ASSIGN(
      auto plaintexts,
      MakeCompressedVector<uint_m>(kTotalSize, indices, kLogCoeffs,
                                   params59_.get(), ntt_params_.get()));
  EXPECT_EQ(plaintexts.size(),
            ceil(kTotalSize / static_cast<double>(1 << kLogCoeffs)));

  std::vector<uint_m> expected_vector(kTotalSize,
                                      uint_m::ImportZero(params59_.get()));
  for (const auto& index : indices) {
    expected_vector[index] = uint_m::ImportOne(params59_.get());
  }

  std::vector<uint_m> coefficients_compressed_vector;
  for (const auto& polynomial : plaintexts) {
    auto coefficients =
        polynomial.InverseNtt(ntt_params_.get(), params59_.get());
    coefficients_compressed_vector.insert(coefficients_compressed_vector.end(),
                                          coefficients.begin(),
                                          coefficients.end());
  }

  EXPECT_GT(coefficients_compressed_vector.size(), kTotalSize);
  for (int i = 0; i < kTotalSize; i++) {
    EXPECT_EQ(expected_vector[i], coefficients_compressed_vector[i]);
  }
  for (int i = kTotalSize; i < coefficients_compressed_vector.size(); i++) {
    EXPECT_EQ(uint_m::ImportZero(params59_.get()),
              coefficients_compressed_vector[i]);
  }
}

TEST_P(ObliviousExpandTest, CompressionFactorOutOfRange) {
  EXPECT_THAT(MakeCompressedVector<uint_m>(kTotalSize, {kTotalSize}, kLogCoeffs,
                                           params59_.get(), ntt_params_.get()),
              StatusIs(::absl::StatusCode::kInvalidArgument,
                       HasSubstr("Index out of range")));
}

TEST_P(ObliviousExpandTest, InputIndexOutOfRange) {
  EXPECT_THAT(
      MakeCompressedVector<uint_m>(kTotalSize, {kTotalSize - 1}, kLogCoeffs + 1,
                                   params59_.get(), ntt_params_.get()),
      StatusIs(::absl::StatusCode::kInvalidArgument,
               HasSubstr("Compression factor")));
}

TEST_P(ObliviousExpandTest, ObliviousExpandProducesQueryVector) {
  // Verifies that the ObliviousExpand algorithm produces the expected vector.
  std::vector<uint_m::Int> zero(kCoeffs, 0);
  std::vector<uint_m::Int> one(zero);
  one[0] = 1;

  // Create expander.
  ASSERT_OK_AND_ASSIGN(
      auto expander,
      GaloisKeysObliviousExpander<uint_m>::Create(
          galois_keys_, kLogT, params59_.get(), ntt_params_.get()));

  // We test when kTestIndex is in the range [0, 2^3], and [0, 2^7].
  for (unsigned int log_compression_factor : {4, 6}) {
    std::vector<int> indices = SampleRandomIndices();

    // Create the compressed vector and the galois_keys.
    ASSERT_OK_AND_ASSIGN(
        auto encryptions,
        CompressAndEncryptVector(kTotalSize, indices, log_compression_factor));

    // ObliviousExpand.
    ASSERT_OK_AND_ASSIGN(
        auto res, expander->ObliviousExpand(encryptions, log_compression_factor,
                                            kTotalSize));

    // Expect that the ciphertexts at index in indices is Enc(1) and the rest
    // are Enc(0).
    EXPECT_EQ(res.size(), kTotalSize);
    for (int i = 0; i < res.size(); i++) {
      ASSERT_OK_AND_ASSIGN(std::vector<uint_m::Int> one_or_zero,
                           Decrypt(*key_, res[i]));
      // Ensure it has the correct value.
      if (std::count(indices.begin(), indices.end(), i)) {
        EXPECT_EQ(one_or_zero, one);
      } else {
        EXPECT_EQ(one_or_zero, zero);
      }
    }
  }
}

TEST_P(ObliviousExpandTest, GaloisKeysMalformed) {
  // Drop the first Galois Key.
  std::vector<GaloisKey<uint_m>> bad_galois_keys = {galois_keys_.begin() + 1,
                                                    galois_keys_.end()};
  // Create expander.
  EXPECT_THAT(GaloisKeysObliviousExpander<uint_m>::Create(
                  bad_galois_keys, kLogT, params59_.get(), ntt_params_.get()),
              StatusIs(::absl::StatusCode::kInvalidArgument,
                       HasSubstr("Substitution power")));
}

TEST_P(ObliviousExpandTest, OutputLengthOutOfRange) {
  ASSERT_OK_AND_ASSIGN(auto req, CompressAndEncryptVector(
                                     kTotalSize, {kTotalSize - 1}, kLogCoeffs));
  ASSERT_OK_AND_ASSIGN(
      auto expander,
      GaloisKeysObliviousExpander<uint_m>::Create(
          {galois_keys_[0]}, kLogT, params59_.get(), ntt_params_.get()));

  // Can only produce (kCoeffs * req.size) ciphertexts.
  EXPECT_THAT(expander->ObliviousExpand(req, /*levels_of_expand=*/1,
                                        1 + (kCoeffs * req.size())),
              StatusIs(::absl::StatusCode::kInvalidArgument,
                       HasSubstr("Output length")));
}

TEST_P(ObliviousExpandTest, ObliviousExpandWithKeyGenerator) {
  // Verifies that the server produces the expected truncated expanded query
  // vector using a single GaloisKey. The compression factor needs to be
  // <= kLogCoeffs - 1.
  std::vector<uint_m::Int> zero(kCoeffs, 0);
  std::vector<uint_m::Int> one(zero);
  one[0] = 1;

  ASSERT_OK_AND_ASSIGN(
      auto generator,
      GaloisKey<uint_m>::Create(
          *key_, prng_type_,
          GaloisGeneratorObliviousExpander<uint_m>::GetSubstitutionGenerator(),
          10));
  ASSERT_OK_AND_ASSIGN(
      auto expander, GaloisGeneratorObliviousExpander<uint_m>::Create(
                         generator, kLogT, params59_.get(), ntt_params_.get()));

  // We test when kTestIndex is in the range [0, 2^3], and [0, 2^7].
  for (unsigned int log_compression_factor : {Uint64{4}, kLogCoeffs - 1}) {
    std::vector<int> indices = SampleRandomIndices();

    // Create the compressed encryptions and the galois_keys.
    ASSERT_OK_AND_ASSIGN(
        auto encryptions,
        CompressAndEncryptVector(kTotalSize, indices, log_compression_factor));

    // ObliviousExpand.
    ASSERT_OK_AND_ASSIGN(
        auto res, expander->ObliviousExpand(encryptions, log_compression_factor,
                                            kTotalSize));

    // Expect that the ciphertexts at index in indices is Enc(1) and the rest
    // are Enc(0).
    EXPECT_EQ(res.size(), kTotalSize);
    for (int i = 0; i < res.size(); i++) {
      ASSERT_OK_AND_ASSIGN(std::vector<uint_m::Int> one_or_zero,
                           Decrypt(*key_, res[i]));
      // Ensure it has the correct value.
      if (std::count(indices.begin(), indices.end(), i)) {
        EXPECT_EQ(one_or_zero, one);
      } else {
        EXPECT_EQ(one_or_zero, zero);
      }
    }
  }
}

TEST_P(ObliviousExpandTest,
       ObliviousExpandWithKeyGeneratorFailsFullCompression) {
  // One cannot use a full compression factor when a generator is used.
  unsigned int log_compression_factor = kLogCoeffs;

  std::vector<uint_m::Int> zero(kCoeffs, 0);
  std::vector<uint_m::Int> one(zero);
  one[0] = 1;

  ASSERT_OK_AND_ASSIGN(
      auto generator,
      GaloisKey<uint_m>::Create(
          *key_, prng_type_,
          GaloisGeneratorObliviousExpander<uint_m>::GetSubstitutionGenerator(),
          10));
  ASSERT_OK_AND_ASSIGN(
      auto expander, GaloisGeneratorObliviousExpander<uint_m>::Create(
                         generator, kLogT, params59_.get(), ntt_params_.get()));

  std::vector<int> indices = SampleRandomIndices();

  // Create the compressed encryptions and the galois_keys.
  ASSERT_OK_AND_ASSIGN(
      auto encryptions,
      CompressAndEncryptVector(kTotalSize, indices, log_compression_factor));

  // ObliviousExpand fails

  EXPECT_THAT(expander->ObliviousExpand(encryptions, log_compression_factor,
                                        kTotalSize),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr(absl::StrCat(
                           "be due to using levels_of_expand = ", kLogCoeffs,
                           " when generating the ObliviousExpander."))));
}

TEST_P(ObliviousExpandTest, DefaultObliviousExpand) {
  // Verifies that the default expand is a no-op.
  ASSERT_OK_AND_ASSIGN(auto expander,
                       DefaultObliviousExpander<uint_m>::Create(
                           kLogT, params59_.get(), ntt_params_.get()));
  std::vector<int> indices = SampleRandomIndices();

  // Create the compressed vector and the galois_keys.
  ASSERT_OK_AND_ASSIGN(auto encryptions,
                       CompressAndEncryptVector(kTotalSize, indices,
                                                /*log_compression_factor=*/0));

  // ObliviousExpand.
  ASSERT_OK_AND_ASSIGN(
      auto res, expander->ObliviousExpand(encryptions, /*levels_of_expand=*/0,
                                          kTotalSize));

  // Expect that the result is encrypts the same message as the request.
  EXPECT_EQ(res.size(), kTotalSize);
  for (int i = 0; i < res.size(); i++) {
    ASSERT_OK_AND_ASSIGN(auto dec1, Decrypt(*key_, res[i]));
    ASSERT_OK_AND_ASSIGN(auto dec2, Decrypt(*key_, encryptions[i]));
    EXPECT_EQ(dec1, dec2);
  }
}

INSTANTIATE_TEST_SUITE_P(ParameterizedTest, ObliviousExpandTest,
                         ::testing::Values(rlwe::PRNG_TYPE_CHACHA,
                                           rlwe::PRNG_TYPE_HKDF));

}  // namespace
}  // namespace rlwe
