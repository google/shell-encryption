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

#include "shell_encryption/aux_relinearization_key.h"

#include <memory>
#include <random>

#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include "absl/status/status.h"
#include "shell_encryption/constants.h"
#include "shell_encryption/integral_types.h"
#include "shell_encryption/montgomery.h"
#include "shell_encryption/ntt_parameters.h"
#include "shell_encryption/polynomial.h"
#include "shell_encryption/prng/chacha_prng.h"
#include "shell_encryption/prng/hkdf_prng.h"
#include "shell_encryption/prng/single_thread_chacha_prng.h"
#include "shell_encryption/prng/single_thread_hkdf_prng.h"
#include "shell_encryption/serialization.pb.h"
#include "shell_encryption/status_macros.h"
#include "shell_encryption/symmetric_encryption.h"
#include "shell_encryption/testing/status_matchers.h"
#include "shell_encryption/testing/status_testing.h"
#include "shell_encryption/testing/testing_prng.h"
#include "shell_encryption/testing/testing_utils.h"

namespace rlwe {
namespace {

// Useful typedefs.
using ModularInt = MontgomeryInt<uint64_t>;
using Polynomial = Polynomial<ModularInt>;
using Ciphertext = SymmetricRlweCiphertext<ModularInt>;
using ModularIntParams = ModularInt::Params;
using NttParameters = NttParameters<ModularInt>;
using Key = SymmetricRlweKey<ModularInt>;
using AuxModRelinearizationKey = AuxModRelinearizationKey<ModularInt>;
using ErrorParams = ErrorParams<ModularInt>;
using ::rlwe::testing::StatusIs;
using ::testing::HasSubstr;

// Set constants.
constexpr int kLogPlaintextModulus = 1;
constexpr uint64_t kPlaintextModulus = (1ULL << kLogPlaintextModulus) + 1;  // t
constexpr int kDefaultVariance = 8;
constexpr int kLogCoeffs = 10;
constexpr int kCoeffs = 1 << kLogCoeffs;  // 1024

// We use two groups of moduli (q, p) in this test, both support NTT in the ring
// Z[X]/(X^1024 + 1, q) and Z[X]/(X^1024 + 1, p), where q is the ciphertext
// modulus and p is the auxiliary modulus.
// Note that these moduli do not offer sufficient security
// as p*q is too small for lattice dimension 1024; they are chosen only for
// testing correctness.
constexpr uint64_t kMainModulus0 = kModulus25;         // q = 33538049
constexpr uint64_t kAuxiliaryModulus0 = 134215681ULL;  // 26bit
constexpr uint64_t kMainModulus1 = kModulus44;         // q = 17592169240577
constexpr uint64_t kAuxiliaryModulus1 = 1152921504606830593ULL;  // 60bit
// Seed for the RNG
constexpr unsigned int kRandSeed = 1;

class AuxModRelinearizationKeyTest : public ::testing::TestWithParam<PrngType> {
 public:
  AuxModRelinearizationKeyTest()
      : prng_type_(GetParam()), mt_rand_(kRandSeed) {}

 protected:
  void SetUp() override {
    ASSERT_OK_AND_ASSIGN(mod_params_small_main_,
                         ModularInt::Params::Create(kMainModulus0));
    ASSERT_OK_AND_ASSIGN(mod_params_small_aux_,
                         ModularInt::Params::Create(kAuxiliaryModulus0));
    ASSERT_OK_AND_ASSIGN(mod_params_large_main_,
                         ModularInt::Params::Create(kMainModulus1));
    ASSERT_OK_AND_ASSIGN(mod_params_large_aux_,
                         ModularInt::Params::Create(kAuxiliaryModulus1));
    ASSERT_OK_AND_ASSIGN(auto ntt_params_small_main,
                         InitializeNttParameters<ModularInt>(
                             kLogCoeffs, mod_params_small_main_.get()));
    ASSERT_OK_AND_ASSIGN(auto ntt_params_small_aux,
                         InitializeNttParameters<ModularInt>(
                             kLogCoeffs, mod_params_small_aux_.get()));
    ASSERT_OK_AND_ASSIGN(auto ntt_params_large_main,
                         InitializeNttParameters<ModularInt>(
                             kLogCoeffs, mod_params_large_main_.get()));
    ASSERT_OK_AND_ASSIGN(auto ntt_params_large_aux,
                         InitializeNttParameters<ModularInt>(
                             kLogCoeffs, mod_params_large_aux_.get()));

    ntt_params_small_main_ =
        std::make_unique<const NttParameters>(std::move(ntt_params_small_main));
    ntt_params_small_aux_ =
        std::make_unique<const NttParameters>(std::move(ntt_params_small_aux));
    ntt_params_large_main_ =
        std::make_unique<const NttParameters>(std::move(ntt_params_large_main));
    ntt_params_large_aux_ =
        std::make_unique<const NttParameters>(std::move(ntt_params_large_aux));

    ASSERT_OK_AND_ASSIGN(
        auto error_params_small,
        ErrorParams::Create(kLogPlaintextModulus, kDefaultVariance,
                            mod_params_small_main_.get(),
                            ntt_params_small_main_.get()));
    ASSERT_OK_AND_ASSIGN(
        auto error_params_large,
        ErrorParams::Create(kLogPlaintextModulus, kDefaultVariance,
                            mod_params_large_main_.get(),
                            ntt_params_large_main_.get()));
    error_params_small_ =
        std::make_unique<const ErrorParams>(error_params_small);
    error_params_large_ =
        std::make_unique<const ErrorParams>(error_params_large);
  }

  StatusOr<std::string> GenerateSeed() {
    return testing::GenerateSeed(prng_type_);
  }

  StatusOr<std::unique_ptr<SecurePrng>> CreatePrng(absl::string_view seed) {
    return testing::CreatePrng(seed, prng_type_);
  }

  // Convert a vector of integers to a vector of montgomery integers.
  StatusOr<std::vector<ModularInt>> ConvertToMontgomery(
      const std::vector<ModularInt::Int>& coeffs,
      const ModularInt::Params* params) {
    std::vector<ModularInt> output(coeffs.size(),
                                   ModularInt::ImportZero(params));
    for (int i = 0; i < output.size(); ++i) {
      RLWE_ASSIGN_OR_RETURN(output[i],
                            ModularInt::ImportInt(coeffs[i], params));
    }
    return output;
  }

  // Sample a secret key.
  StatusOr<Key> SampleKey(const ModularIntParams* mod_params,
                          const NttParameters* ntt_params,
                          int log_t = kLogPlaintextModulus) {
    RLWE_ASSIGN_OR_RETURN(std::string prng_seed, GenerateSeed());
    RLWE_ASSIGN_OR_RETURN(auto prng, CreatePrng(prng_seed));
    return Key::Sample(kLogCoeffs, kDefaultVariance, log_t, mod_params,
                       ntt_params, prng.get());
  }

  // Sample a random plaintext and convert it to NTT form.
  StatusOr<Polynomial> SamplePlaintextNtt(
      const ModularIntParams* mod_params, const NttParameters* ntt_params,
      ModularInt::Int t = kPlaintextModulus) {
    std::vector<ModularInt::Int> plaintext(kCoeffs);
    for (int i = 0; i < kCoeffs; ++i) {
      plaintext[i] = mt_rand_() % t;
    }
    RLWE_ASSIGN_OR_RETURN(auto plaintext_m,
                          ConvertToMontgomery(plaintext, mod_params));
    return Polynomial::ConvertToNtt(plaintext_m, ntt_params, mod_params);
  }

  // Encrypt a plaintext.
  StatusOr<Ciphertext> Encrypt(const Key& key, const Polynomial& plaintext_ntt,
                               const ErrorParams* error_params) {
    RLWE_ASSIGN_OR_RETURN(std::string prng_seed, GenerateSeed());
    RLWE_ASSIGN_OR_RETURN(auto prng, CreatePrng(prng_seed));
    return ::rlwe::Encrypt<ModularInt>(key, plaintext_ntt, error_params,
                                       prng.get());
  }

  PrngType prng_type_;
  std::mt19937 mt_rand_;

  std::unique_ptr<const ModularIntParams> mod_params_small_main_;
  std::unique_ptr<const ModularIntParams> mod_params_small_aux_;
  std::unique_ptr<const ModularIntParams> mod_params_large_main_;
  std::unique_ptr<const ModularIntParams> mod_params_large_aux_;
  std::unique_ptr<const NttParameters> ntt_params_small_main_;
  std::unique_ptr<const NttParameters> ntt_params_small_aux_;
  std::unique_ptr<const NttParameters> ntt_params_large_main_;
  std::unique_ptr<const NttParameters> ntt_params_large_aux_;

  std::unique_ptr<const ErrorParams> error_params_small_;
  std::unique_ptr<const ErrorParams> error_params_large_;
};

// Checks for error handling when creating a relinearization key to have
// `degree` == 0.
TEST_P(AuxModRelinearizationKeyTest, ErrorIfDegreeIsZero) {
  ASSERT_OK_AND_ASSIGN(auto key, SampleKey(mod_params_small_main_.get(),
                                           ntt_params_small_main_.get()));
  EXPECT_THAT(
      AuxModRelinearizationKey::Create(key, prng_type_,
                                       /*degree=*/0,
                                       mod_params_small_aux_.get(),
                                       ntt_params_small_aux_.get()),
      StatusIs(absl::StatusCode::kInvalidArgument, HasSubstr("degree")));
}

// Checks for error handling when creating a relinearization key to have
// negative `degree`.
TEST_P(AuxModRelinearizationKeyTest, ErrorIfDegreeIsNegative) {
  ASSERT_OK_AND_ASSIGN(auto key, SampleKey(mod_params_small_main_.get(),
                                           ntt_params_small_main_.get()));
  EXPECT_THAT(
      AuxModRelinearizationKey::Create(key, prng_type_,
                                       /*degree=*/-1,
                                       mod_params_small_aux_.get(),
                                       ntt_params_small_aux_.get()),
      StatusIs(absl::StatusCode::kInvalidArgument, HasSubstr("degree")));
}

// Checks for error handling when creating a relinearization key with a null
// mod_params_main parameter.
TEST_P(AuxModRelinearizationKeyTest, ErrorIfModParamsAuxIsNull) {
  ASSERT_OK_AND_ASSIGN(auto key, SampleKey(mod_params_small_main_.get(),
                                           ntt_params_small_main_.get()));
  EXPECT_THAT(AuxModRelinearizationKey::Create(key, prng_type_,
                                               /*degree=*/2,
                                               /*mod_params_aux=*/nullptr,
                                               ntt_params_small_aux_.get()),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("mod_params_aux must not be null")));
}

// Checks for error handling when creating a relinearization key with a null
// ntt_params_main parameter.
TEST_P(AuxModRelinearizationKeyTest, ErrorIfNttParamsAuxIsNull) {
  ASSERT_OK_AND_ASSIGN(auto key, SampleKey(mod_params_small_main_.get(),
                                           ntt_params_small_main_.get()));
  EXPECT_THAT(AuxModRelinearizationKey::Create(key, prng_type_,
                                               /*degree=*/2,
                                               mod_params_small_aux_.get(),
                                               /*ntt_params_aux=*/nullptr),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("ntt_params_aux must not be null")));
}

// Checks for error handling when creating a relinearization key to have
// `substitution_power` == 0.
TEST_P(AuxModRelinearizationKeyTest, ErrorIfSubstitutionPowerIsZero) {
  ASSERT_OK_AND_ASSIGN(auto key, SampleKey(mod_params_small_main_.get(),
                                           ntt_params_small_main_.get()));
  EXPECT_THAT(AuxModRelinearizationKey::Create(key, prng_type_,
                                               /*degree=*/2,
                                               mod_params_small_aux_.get(),
                                               ntt_params_small_aux_.get(),
                                               /*substitution_power=*/0),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("substitution_power")));
}

// Checks for error handling when creating a relinearization key to have
// negative `substitution_power`.
TEST_P(AuxModRelinearizationKeyTest, ErrorIfSubstitutionPowerIsNegative) {
  ASSERT_OK_AND_ASSIGN(auto key, SampleKey(mod_params_small_main_.get(),
                                           ntt_params_small_main_.get()));
  EXPECT_THAT(AuxModRelinearizationKey::Create(key, prng_type_,
                                               /*degree=*/2,
                                               mod_params_small_aux_.get(),
                                               ntt_params_small_aux_.get(),
                                               /*substitution_power=*/-1),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("substitution_power")));
}

// Checks for error handling when creating a relinearization key with an
// invalid PrngType.
TEST_P(AuxModRelinearizationKeyTest, ErrorIfInvalidPrngType) {
  ASSERT_OK_AND_ASSIGN(auto key, SampleKey(mod_params_small_main_.get(),
                                           ntt_params_small_main_.get()));
  EXPECT_THAT(AuxModRelinearizationKey::Create(key, PRNG_TYPE_INVALID,
                                               /*degree=*/2,
                                               mod_params_small_aux_.get(),
                                               ntt_params_small_aux_.get()),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("PrngType not specified correctly")));
}

TEST_P(AuxModRelinearizationKeyTest, ErrorIfModulusMismatchInApplyTo) {
  // Create a relinearization key against the small modulus
  ASSERT_OK_AND_ASSIGN(auto key, SampleKey(mod_params_small_main_.get(),
                                           ntt_params_small_main_.get()));
  ASSERT_OK_AND_ASSIGN(auto relinearization_key,
                       AuxModRelinearizationKey::Create(
                           key, prng_type_,
                           /*degree=*/1, mod_params_small_aux_.get(),
                           ntt_params_small_aux_.get(),
                           /*substitution_power=*/5));

  AuxModRelinearizationKey another_key = relinearization_key;
  another_key = std::move(relinearization_key);

  // Encrypt a plaintext under a key wrt to the large modulus
  ASSERT_OK_AND_ASSIGN(auto key_large, SampleKey(mod_params_large_main_.get(),
                                                 ntt_params_large_main_.get()));
  ASSERT_OK_AND_ASSIGN(auto plaintext_ntt,
                       SamplePlaintextNtt(mod_params_large_main_.get(),
                                          ntt_params_large_main_.get()));
  ASSERT_OK_AND_ASSIGN(auto ciphertext, Encrypt(key_large, plaintext_ntt,
                                                error_params_large_.get()));

  // Apply the relinearization key to the ciphertext should result in an error
  EXPECT_THAT(
      relinearization_key.ApplyTo(ciphertext),
      StatusIs(
          absl::StatusCode::kInvalidArgument,
          HasSubstr(
              "Ciphertext modulus does not match key-switching key modulus")));
}

TEST_P(AuxModRelinearizationKeyTest,
       ErrorIfCiphertextDegreeIsSmallerInApplyTo) {
  // Create a relinearization key of degree 2
  ASSERT_OK_AND_ASSIGN(auto key, SampleKey(mod_params_small_main_.get(),
                                           ntt_params_small_main_.get()));
  ASSERT_OK_AND_ASSIGN(auto relinearization_key,
                       AuxModRelinearizationKey::Create(
                           key, prng_type_,
                           /*degree=*/2, mod_params_small_aux_.get(),
                           ntt_params_small_aux_.get()));
  ASSERT_EQ(relinearization_key.Degree(), 2);

  // Encrypt a plaintext and get a ciphertext of degree 1
  ASSERT_OK_AND_ASSIGN(auto plaintext_ntt,
                       SamplePlaintextNtt(mod_params_small_main_.get(),
                                          ntt_params_small_main_.get()));
  ASSERT_OK_AND_ASSIGN(auto ciphertext,
                       Encrypt(key, plaintext_ntt, error_params_small_.get()));

  // Apply the relinearization key to the ciphertext should result in an error
  EXPECT_THAT(
      relinearization_key.ApplyTo(ciphertext),
      StatusIs(absl::StatusCode::kInvalidArgument,
               HasSubstr("Ciphertext has incompatible number of components")));
}

TEST_P(AuxModRelinearizationKeyTest, ErrorIfCiphertextDegreeIsLargerInApplyTo) {
  // Create a relinearization key of degree 1
  ASSERT_OK_AND_ASSIGN(auto key, SampleKey(mod_params_small_main_.get(),
                                           ntt_params_small_main_.get()));
  ASSERT_OK_AND_ASSIGN(auto relinearization_key,
                       AuxModRelinearizationKey::Create(
                           key, prng_type_,
                           /*degree=*/1, mod_params_small_aux_.get(),
                           ntt_params_small_aux_.get(),
                           /*substitution_power=*/5));
  ASSERT_EQ(relinearization_key.Degree(), 1);

  // Encrypt a plaintext and get a ciphertext of degree 1
  ASSERT_OK_AND_ASSIGN(auto plaintext_ntt,
                       SamplePlaintextNtt(mod_params_small_main_.get(),
                                          ntt_params_small_main_.get()));
  ASSERT_OK_AND_ASSIGN(auto ciphertext,
                       Encrypt(key, plaintext_ntt, error_params_small_.get()));
  // Compute the square and get a ciphertext of degree 2
  ASSERT_OK_AND_ASSIGN(auto square, ciphertext* ciphertext);
  ASSERT_EQ(square.Len(), 3);  // length == degree + 1

  // Apply the relinearization key to the square should result in an error
  EXPECT_THAT(
      relinearization_key.ApplyTo(square),
      StatusIs(absl::StatusCode::kInvalidArgument,
               HasSubstr("Ciphertext has incompatible number of components")));
}

TEST_P(AuxModRelinearizationKeyTest,
       ErrorIfDeserializeWithZeroNumberOfComponents) {
  ASSERT_OK_AND_ASSIGN(auto key, SampleKey(mod_params_small_main_.get(),
                                           ntt_params_small_main_.get()));
  ASSERT_OK_AND_ASSIGN(auto relinearization_key,
                       AuxModRelinearizationKey::Create(
                           key, prng_type_,
                           /*degree=*/2, mod_params_small_aux_.get(),
                           ntt_params_small_aux_.get()));
  ASSERT_EQ(relinearization_key.Degree(), 2);

  // Serialize.
  ASSERT_OK_AND_ASSIGN(SerializedAuxModRelinearizationKey serialized,
                       relinearization_key.Serialize());
  ASSERT_EQ(serialized.num_components(), 1);  // num_components + 1 == degree
  // Tamper the serialization with an invalid `num_components`.
  serialized.set_num_components(0);
  // Deserialize should result in an error.
  EXPECT_THAT(
      AuxModRelinearizationKey::Deserialize(
          serialized, mod_params_small_main_.get(),
          ntt_params_small_main_.get(), mod_params_small_aux_.get(),
          ntt_params_small_aux_.get(),
          key.PlaintextModulus().ExportInt(mod_params_small_main_.get())),
      StatusIs(absl::StatusCode::kInvalidArgument,
               HasSubstr("The number of components")));
}

TEST_P(AuxModRelinearizationKeyTest,
       ErrorIfDeserializeWithNegativeNumberOfComponents) {
  ASSERT_OK_AND_ASSIGN(auto key, SampleKey(mod_params_small_main_.get(),
                                           ntt_params_small_main_.get()));
  ASSERT_OK_AND_ASSIGN(auto relinearization_key,
                       AuxModRelinearizationKey::Create(
                           key, prng_type_,
                           /*degree=*/2, mod_params_small_aux_.get(),
                           ntt_params_small_aux_.get()));
  ASSERT_EQ(relinearization_key.Degree(), 2);

  // Serialize.
  ASSERT_OK_AND_ASSIGN(SerializedAuxModRelinearizationKey serialized,
                       relinearization_key.Serialize());
  ASSERT_EQ(serialized.num_components(), 1);  // num_components + 1 == degree
  // Tamper the serialization with an invalid `num_components`.
  serialized.set_num_components(-1);
  // Deserialize should result in an error.
  EXPECT_THAT(
      AuxModRelinearizationKey::Deserialize(
          serialized, mod_params_small_main_.get(),
          ntt_params_small_main_.get(), mod_params_small_aux_.get(),
          ntt_params_small_aux_.get(),
          key.PlaintextModulus().ExportInt(mod_params_small_main_.get())),
      StatusIs(absl::StatusCode::kInvalidArgument,
               HasSubstr("The number of components")));
}

TEST_P(AuxModRelinearizationKeyTest,
       ErrorIfDeserializeWithMismatchNumberOfComponents) {
  ASSERT_OK_AND_ASSIGN(auto key, SampleKey(mod_params_small_main_.get(),
                                           ntt_params_small_main_.get()));
  ASSERT_OK_AND_ASSIGN(auto relinearization_key,
                       AuxModRelinearizationKey::Create(
                           key, prng_type_,
                           /*degree=*/2, mod_params_small_aux_.get(),
                           ntt_params_small_aux_.get()));
  ASSERT_EQ(relinearization_key.Degree(), 2);

  // Serialize.
  ASSERT_OK_AND_ASSIGN(SerializedAuxModRelinearizationKey serialized,
                       relinearization_key.Serialize());
  ASSERT_EQ(serialized.num_components(), 1);  // num_components + 1 == degree
  ASSERT_EQ(serialized.b_size(), 2);          // b_size == num_components * 2
  // Tamper the serialization with a `num_components` value that violates the
  // invariant with the size of the serialized "b" array.
  serialized.set_num_components(2);
  // Deserialize should result in an error.
  EXPECT_THAT(
      AuxModRelinearizationKey::Deserialize(
          serialized, mod_params_small_main_.get(),
          ntt_params_small_main_.get(), mod_params_small_aux_.get(),
          ntt_params_small_aux_.get(),
          key.PlaintextModulus().ExportInt(mod_params_small_main_.get())),
      StatusIs(absl::StatusCode::kInvalidArgument,
               HasSubstr("The length of serialized")));
}

TEST_P(AuxModRelinearizationKeyTest, ErrorIfDeserializeWithNullModParamsMain) {
  ASSERT_OK_AND_ASSIGN(auto key, SampleKey(mod_params_small_main_.get(),
                                           ntt_params_small_main_.get()));
  ASSERT_OK_AND_ASSIGN(auto relinearization_key,
                       AuxModRelinearizationKey::Create(
                           key, prng_type_,
                           /*degree=*/2, mod_params_small_aux_.get(),
                           ntt_params_small_aux_.get()));

  // Serialize.
  ASSERT_OK_AND_ASSIGN(SerializedAuxModRelinearizationKey serialized,
                       relinearization_key.Serialize());
  // Deserialize with a null mod_params_main argument.
  EXPECT_THAT(
      AuxModRelinearizationKey::Deserialize(
          serialized, /*mod_params_main=*/nullptr, ntt_params_small_main_.get(),
          mod_params_small_aux_.get(), ntt_params_small_aux_.get(),
          key.PlaintextModulus().ExportInt(mod_params_small_main_.get())),
      StatusIs(absl::StatusCode::kInvalidArgument,
               HasSubstr("mod_params_main must not be null")));
}

TEST_P(AuxModRelinearizationKeyTest, ErrorIfDeserializeWithNullNttParamsMain) {
  ASSERT_OK_AND_ASSIGN(auto key, SampleKey(mod_params_small_main_.get(),
                                           ntt_params_small_main_.get()));
  ASSERT_OK_AND_ASSIGN(auto relinearization_key,
                       AuxModRelinearizationKey::Create(
                           key, prng_type_,
                           /*degree=*/2, mod_params_small_aux_.get(),
                           ntt_params_small_aux_.get()));

  // Serialize.
  ASSERT_OK_AND_ASSIGN(SerializedAuxModRelinearizationKey serialized,
                       relinearization_key.Serialize());
  // Deserialize with a null ntt_params_main argument.
  EXPECT_THAT(
      AuxModRelinearizationKey::Deserialize(
          serialized, mod_params_small_main_.get(),
          /*ntt_params_main=*/nullptr, mod_params_small_aux_.get(),
          ntt_params_small_aux_.get(),
          key.PlaintextModulus().ExportInt(mod_params_small_main_.get())),
      StatusIs(absl::StatusCode::kInvalidArgument,
               HasSubstr("ntt_params_main must not be null")));
}

TEST_P(AuxModRelinearizationKeyTest, ErrorIfDeserializeWithNullModParamsAux) {
  ASSERT_OK_AND_ASSIGN(auto key, SampleKey(mod_params_small_main_.get(),
                                           ntt_params_small_main_.get()));
  ASSERT_OK_AND_ASSIGN(auto relinearization_key,
                       AuxModRelinearizationKey::Create(
                           key, prng_type_,
                           /*degree=*/2, mod_params_small_aux_.get(),
                           ntt_params_small_aux_.get()));

  // Serialize.
  ASSERT_OK_AND_ASSIGN(SerializedAuxModRelinearizationKey serialized,
                       relinearization_key.Serialize());
  // Deserialize with a null mod_params_main argument.
  EXPECT_THAT(
      AuxModRelinearizationKey::Deserialize(
          serialized, mod_params_small_main_.get(),
          ntt_params_small_main_.get(), /*mod_params_aux=*/nullptr,
          ntt_params_small_aux_.get(),
          key.PlaintextModulus().ExportInt(mod_params_small_main_.get())),
      StatusIs(absl::StatusCode::kInvalidArgument,
               HasSubstr("mod_params_aux must not be null")));
}

TEST_P(AuxModRelinearizationKeyTest, ErrorIfDeserializeWithNullNttParamsAux) {
  ASSERT_OK_AND_ASSIGN(auto key, SampleKey(mod_params_small_main_.get(),
                                           ntt_params_small_main_.get()));
  ASSERT_OK_AND_ASSIGN(auto relinearization_key,
                       AuxModRelinearizationKey::Create(
                           key, prng_type_,
                           /*degree=*/2, mod_params_small_aux_.get(),
                           ntt_params_small_aux_.get()));

  // Serialize.
  ASSERT_OK_AND_ASSIGN(SerializedAuxModRelinearizationKey serialized,
                       relinearization_key.Serialize());
  // Deserialize with a null ntt_params_main argument.
  EXPECT_THAT(
      AuxModRelinearizationKey::Deserialize(
          serialized, mod_params_small_main_.get(),
          ntt_params_small_main_.get(), mod_params_small_aux_.get(),
          /*ntt_params_aux=*/nullptr,
          key.PlaintextModulus().ExportInt(mod_params_small_main_.get())),
      StatusIs(absl::StatusCode::kInvalidArgument,
               HasSubstr("ntt_params_aux must not be null")));
}

TEST_P(AuxModRelinearizationKeyTest, ErrorIfDeserializeWithInvalidPrngType) {
  ASSERT_OK_AND_ASSIGN(auto key, SampleKey(mod_params_small_main_.get(),
                                           ntt_params_small_main_.get()));
  ASSERT_OK_AND_ASSIGN(auto relinearization_key,
                       AuxModRelinearizationKey::Create(
                           key, prng_type_,
                           /*degree=*/2, mod_params_small_aux_.get(),
                           ntt_params_small_aux_.get()));
  ASSERT_EQ(relinearization_key.Degree(), 2);

  // Serialize.
  ASSERT_OK_AND_ASSIGN(SerializedAuxModRelinearizationKey serialized,
                       relinearization_key.Serialize());
  // Tamper the serialization with an invalid `prng_type` value.
  serialized.set_prng_type(PRNG_TYPE_INVALID);
  // Deserialize should result in an error.
  EXPECT_THAT(
      AuxModRelinearizationKey::Deserialize(
          serialized, mod_params_small_main_.get(),
          ntt_params_small_main_.get(), mod_params_small_aux_.get(),
          ntt_params_small_aux_.get(),
          key.PlaintextModulus().ExportInt(mod_params_small_main_.get())),
      StatusIs(absl::StatusCode::kInvalidArgument,
               HasSubstr("Invalid PRNG type is specified")));
}

// Checks that relinearization key can correctly convert a degree-2 ciphertext
// (c0, c1, c2) to a standard degree-1 ciphertext (c0, c1).
TEST_P(AuxModRelinearizationKeyTest,
       RelinearizationIsCorrectOnCiphertextOfDegree2) {
  auto mod_params_main = mod_params_small_main_.get();
  auto ntt_params_main = ntt_params_small_main_.get();
  auto mod_params_aux = mod_params_small_aux_.get();
  auto ntt_params_aux = ntt_params_small_aux_.get();
  auto error_params = error_params_small_.get();

  constexpr int kDegree = 2;
  ASSERT_OK_AND_ASSIGN(auto key, SampleKey(mod_params_main, ntt_params_main));
  ASSERT_OK_AND_ASSIGN(
      auto relinearization_key,
      AuxModRelinearizationKey::Create(key, prng_type_, kDegree, mod_params_aux,
                                       ntt_params_aux));
  ASSERT_EQ(relinearization_key.Degree(), kDegree);

  // Sample two plaintexts.
  ASSERT_OK_AND_ASSIGN(Polynomial plaintext1_ntt,
                       SamplePlaintextNtt(mod_params_main, ntt_params_main));
  ASSERT_OK_AND_ASSIGN(Polynomial plaintext2_ntt,
                       SamplePlaintextNtt(mod_params_main, ntt_params_main));

  // Encrypt two plaintexts and then multiply the ciphertexts to get a product
  // ciphertext (c0, c1, c2). Then apply the relinearization key.
  ASSERT_OK_AND_ASSIGN(auto ciphertext1,
                       Encrypt(key, plaintext1_ntt, error_params));
  ASSERT_OK_AND_ASSIGN(auto ciphertext2,
                       Encrypt(key, plaintext2_ntt, error_params));
  ASSERT_OK_AND_ASSIGN(auto product, ciphertext1* ciphertext2);
  ASSERT_EQ(product.Len(), 3);  // product has three components c0, c1, c2
  ASSERT_OK_AND_ASSIGN(auto relinearized_product,
                       relinearization_key.ApplyTo(product));
  // The converted ciphertext should be in standard form (c0, c1) wrt the
  // main modulus, and whose power of the secret key should be 1.
  ASSERT_EQ(relinearized_product.Len(), 2);
  ASSERT_EQ(relinearized_product.ModulusParams()->modulus,
            mod_params_main->modulus);
  ASSERT_EQ(relinearized_product.PowerOfS(), 1);

  // Decrypt the relinearized product.
  ASSERT_OK_AND_ASSIGN(std::vector<ModularInt::Int> decrypted,
                       Decrypt<ModularInt>(key, relinearized_product));

  // Create the plaintext result we expect (plaintext1 * plaintext2).
  ASSERT_OK(plaintext1_ntt.MulInPlace(plaintext2_ntt, mod_params_main));
  std::vector<ModularInt::Int> expected = RemoveError<ModularInt>(
      plaintext1_ntt.InverseNtt(ntt_params_main, mod_params_main),
      mod_params_main->modulus, kPlaintextModulus, mod_params_main);
  EXPECT_EQ(decrypted, expected);
}

// Checks that a relinearization key (of degree 3) can correctly convert a
// degree-3 ciphertext (c0, c1, c2, c3) to degree-1 ciphertext (c0, c1).
TEST_P(AuxModRelinearizationKeyTest,
       RelinearizationIsCorrectOnCiphertextOfDegree3) {
  auto mod_params_main = mod_params_large_main_.get();
  auto ntt_params_main = ntt_params_large_main_.get();
  auto mod_params_aux = mod_params_large_aux_.get();
  auto ntt_params_aux = ntt_params_large_aux_.get();
  auto error_params = error_params_large_.get();

  constexpr int kDegree = 3;
  ASSERT_OK_AND_ASSIGN(auto key, SampleKey(mod_params_main, ntt_params_main));
  ASSERT_OK_AND_ASSIGN(
      auto relinearization_key,
      AuxModRelinearizationKey::Create(key, prng_type_, kDegree, mod_params_aux,
                                       ntt_params_aux));
  ASSERT_EQ(relinearization_key.Degree(), kDegree);

  // Sample three plaintexts.
  ASSERT_OK_AND_ASSIGN(Polynomial plaintext1_ntt,
                       SamplePlaintextNtt(mod_params_main, ntt_params_main));
  ASSERT_OK_AND_ASSIGN(Polynomial plaintext2_ntt,
                       SamplePlaintextNtt(mod_params_main, ntt_params_main));
  ASSERT_OK_AND_ASSIGN(Polynomial plaintext3_ntt,
                       SamplePlaintextNtt(mod_params_main, ntt_params_main));

  // Encrypt them, multiply the three ciphertexts together, and then apply the
  // relinearization key on the product.
  ASSERT_OK_AND_ASSIGN(auto ciphertext1,
                       Encrypt(key, plaintext1_ntt, error_params));
  ASSERT_OK_AND_ASSIGN(auto ciphertext2,
                       Encrypt(key, plaintext2_ntt, error_params));
  ASSERT_OK_AND_ASSIGN(auto ciphertext3,
                       Encrypt(key, plaintext3_ntt, error_params));
  ASSERT_OK_AND_ASSIGN(auto product1, ciphertext1* ciphertext2);
  ASSERT_OK_AND_ASSIGN(auto product2, product1* ciphertext3);
  ASSERT_EQ(product2.Len(), 4);  // four components c0, c1, c2, c3
  ASSERT_OK_AND_ASSIGN(auto relinearized_product,
                       relinearization_key.ApplyTo(product2));
  // The converted ciphertext should be in standard form (c0, c1) with power
  // of the secret key being 1, and wrt the main modulus.
  ASSERT_EQ(relinearized_product.Len(), 2);
  ASSERT_EQ(relinearized_product.ModulusParams()->modulus,
            mod_params_main->modulus);
  ASSERT_EQ(relinearized_product.PowerOfS(), 1);

  // Decrypt the relinearized product.
  ASSERT_OK_AND_ASSIGN(std::vector<ModularInt::Int> decrypted,
                       Decrypt<ModularInt>(key, relinearized_product));

  // Create the polynomial we expect.
  ASSERT_OK(plaintext1_ntt.MulInPlace(plaintext2_ntt, mod_params_main));
  ASSERT_OK(plaintext1_ntt.MulInPlace(plaintext3_ntt, mod_params_main));
  std::vector<ModularInt::Int> expected = RemoveError<ModularInt>(
      plaintext1_ntt.InverseNtt(ntt_params_main, mod_params_main),
      mod_params_main->modulus, kPlaintextModulus, mod_params_main);

  EXPECT_EQ(decrypted, expected);
}

// Checks that for certain parameters, even if the auxiliary modulus p is
// smaller than the main modulus q, relinearization can still generate
// correct result.
// Specifically, we set q to be a 44-bit prime, and p 26-bit prime; we then
// set the plaintext modulus t to be a 16-bit prime. When multiplying two fresh
// ciphertexts with a 5-bit error each (in l_inf coeff norm), the product has
// a 10-bit error. We then relinearize the product using the 26-bit auxiliary
// modulus p, introducing roughly a 28-bit error, and hence the resulting
// ciphertext should still be decryptable.
// This test is useful to make sure that we don't introduce unnecessary errors.
TEST_P(AuxModRelinearizationKeyTest,
       RelinearizationWithSmallAuxModIsCorrectOnCiphertextOfDegree2) {
  auto mod_params_main = mod_params_large_main_.get();
  auto ntt_params_main = ntt_params_large_main_.get();
  auto mod_params_aux = mod_params_small_aux_.get();
  auto ntt_params_aux = ntt_params_small_aux_.get();
  auto error_params = error_params_small_.get();

  constexpr int kLogT = 16;
  constexpr int kDegree = 2;
  ASSERT_OK_AND_ASSIGN(auto key,
                       SampleKey(mod_params_main, ntt_params_main, kLogT));
  ASSERT_OK_AND_ASSIGN(
      auto relinearization_key,
      AuxModRelinearizationKey::Create(key, prng_type_, kDegree, mod_params_aux,
                                       ntt_params_aux));
  ASSERT_EQ(relinearization_key.Degree(), kDegree);

  // Sample two plaintexts.
  ModularInt::Int t = key.PlaintextModulus().ExportInt(mod_params_main);
  ASSERT_OK_AND_ASSIGN(Polynomial plaintext1_ntt,
                       SamplePlaintextNtt(mod_params_main, ntt_params_main, t));
  ASSERT_OK_AND_ASSIGN(Polynomial plaintext2_ntt,
                       SamplePlaintextNtt(mod_params_main, ntt_params_main, t));

  // Encrypt two plaintexts and then multiply the ciphertexts to get a product
  // ciphertext (c0, c1, c2). Then apply the relinearization key.
  ASSERT_OK_AND_ASSIGN(auto ciphertext1,
                       Encrypt(key, plaintext1_ntt, error_params));
  ASSERT_OK_AND_ASSIGN(auto ciphertext2,
                       Encrypt(key, plaintext2_ntt, error_params));
  ASSERT_OK_AND_ASSIGN(auto product, ciphertext1* ciphertext2);
  ASSERT_EQ(product.Len(), 3);  // product has three components c0, c1, c2
  ASSERT_OK_AND_ASSIGN(auto relinearized_product,
                       relinearization_key.ApplyTo(product));
  // The converted ciphertext should be in standard form (c0, c1) wrt the
  // main modulus, and whose power of the secret key should be 1.
  ASSERT_EQ(relinearized_product.Len(), 2);
  ASSERT_EQ(relinearized_product.ModulusParams()->modulus,
            mod_params_main->modulus);
  ASSERT_EQ(relinearized_product.PowerOfS(), 1);

  // Decrypt the relinearized product.
  ASSERT_OK_AND_ASSIGN(std::vector<ModularInt::Int> decrypted,
                       Decrypt<ModularInt>(key, relinearized_product));

  // Create the plaintext result we expect (plaintext1 * plaintext2).
  ASSERT_OK(plaintext1_ntt.MulInPlace(plaintext2_ntt, mod_params_main));
  std::vector<ModularInt::Int> expected = RemoveError<ModularInt>(
      plaintext1_ntt.InverseNtt(ntt_params_main, mod_params_main),
      mod_params_main->modulus, t, mod_params_main);

  EXPECT_EQ(decrypted, expected);
}

// Checks that a key-switching key (with `substitution_power` != 1) can
// correctly convert a ciphertext encrypting m(X^j) under secret key s(X^j) to
// ciphertext encrypting m(X^j) under the secret key s = s(X).
TEST_P(AuxModRelinearizationKeyTest,
       KeySwitchingIsCorrectWithNonTrivialSubPower) {
  auto mod_params_main = mod_params_small_main_.get();
  auto ntt_params_main = ntt_params_small_main_.get();
  auto mod_params_aux = mod_params_small_aux_.get();
  auto ntt_params_aux = ntt_params_small_aux_.get();
  auto error_params = error_params_small_.get();

  constexpr int kDegree = 1;
  constexpr int kPower = 5;

  ASSERT_OK_AND_ASSIGN(auto key, SampleKey(mod_params_main, ntt_params_main));
  ASSERT_OK_AND_ASSIGN(auto switching_key, AuxModRelinearizationKey::Create(
                                               key, prng_type_, kDegree,
                                               mod_params_aux, ntt_params_aux,
                                               /*substitution_power=*/kPower));
  ASSERT_EQ(switching_key.Degree(), kDegree);
  ASSERT_EQ(switching_key.SubstitutionPower(), kPower);

  // Sample a plaintext polynomial m(X).
  ASSERT_OK_AND_ASSIGN(Polynomial plaintext_ntt,
                       SamplePlaintextNtt(mod_params_main, ntt_params_main));

  // Encrypt m(X), substitute X with X^j to get m(X^j), and apply the
  // key-switching key.
  ASSERT_OK_AND_ASSIGN(auto ciphertext,
                       Encrypt(key, plaintext_ntt, error_params));
  ASSERT_OK_AND_ASSIGN(auto ciphertext_sub,
                       ciphertext.Substitute(kPower, ntt_params_main));
  ASSERT_EQ(ciphertext_sub.PowerOfS(), kPower);  // sanity check on power
  ASSERT_OK_AND_ASSIGN(auto switched_ciphertext,
                       switching_key.ApplyTo(ciphertext_sub));
  // The converted ciphertext should be in standard form (c0, c1) with power
  // of the secret key being 1, and wrt the main modulus.
  ASSERT_EQ(switched_ciphertext.PowerOfS(), 1);
  ASSERT_EQ(switched_ciphertext.Len(), 2);
  ASSERT_EQ(switched_ciphertext.ModulusParams()->modulus,
            mod_params_main->modulus);

  // Decrypt the converted ciphertext.
  ASSERT_OK_AND_ASSIGN(std::vector<ModularInt::Int> decrypted,
                       Decrypt<ModularInt>(key, switched_ciphertext));

  // Create the polynomial m(X^j) we expect.
  ASSERT_OK_AND_ASSIGN(
      auto plaintext_sub,
      plaintext_ntt.Substitute(kPower, ntt_params_main, mod_params_main));
  std::vector<ModularInt::Int> expected = RemoveError<ModularInt>(
      plaintext_sub.InverseNtt(ntt_params_main, mod_params_main),
      mod_params_main->modulus, kPlaintextModulus, mod_params_main);

  EXPECT_EQ(decrypted, expected);
}

TEST_P(AuxModRelinearizationKeyTest, SerializeAuxModRelinearizationKey) {
  auto mod_params_main = mod_params_small_main_.get();
  auto ntt_params_main = ntt_params_small_main_.get();
  auto mod_params_aux = mod_params_small_aux_.get();
  auto ntt_params_aux = ntt_params_small_aux_.get();
  auto error_params = error_params_small_.get();

  constexpr int kDegree = 1;
  constexpr int kPower = 5;
  ASSERT_OK_AND_ASSIGN(auto key, SampleKey(mod_params_main, ntt_params_main));
  ASSERT_OK_AND_ASSIGN(
      auto relinearization_key,
      AuxModRelinearizationKey::Create(key, prng_type_, kDegree, mod_params_aux,
                                       ntt_params_aux,
                                       /*substitution_power=*/kPower));

  // Serialize and deserialize.
  ASSERT_OK_AND_ASSIGN(SerializedAuxModRelinearizationKey serialized,
                       relinearization_key.Serialize());
  ASSERT_OK_AND_ASSIGN(
      auto deserialized,
      AuxModRelinearizationKey::Deserialize(
          serialized, mod_params_main, ntt_params_main, mod_params_aux,
          ntt_params_aux, key.PlaintextModulus().ExportInt(mod_params_main)));
  ASSERT_EQ(deserialized.Degree(), kDegree);
  ASSERT_EQ(deserialized.SubstitutionPower(), kPower);

  // Create the initial plaintexts.
  ASSERT_OK_AND_ASSIGN(Polynomial plaintext_ntt,
                       SamplePlaintextNtt(mod_params_main, ntt_params_main));

  // Encrypt, substitute, apply the deserialized relinearization key.
  ASSERT_OK_AND_ASSIGN(auto ciphertext,
                       Encrypt(key, plaintext_ntt, error_params));
  ASSERT_OK_AND_ASSIGN(auto ciphertext_sub,
                       ciphertext.Substitute(kPower, ntt_params_main));
  ASSERT_OK_AND_ASSIGN(auto switched_ciphertext,
                       deserialized.ApplyTo(ciphertext_sub));
  // The converted ciphertext should be in standard form (c0, c1) with power
  // of the secret key being 1, and wrt the main modulus.
  ASSERT_EQ(switched_ciphertext.PowerOfS(), 1);
  ASSERT_EQ(switched_ciphertext.Len(), 2);
  ASSERT_EQ(switched_ciphertext.ModulusParams()->modulus,
            mod_params_main->modulus);

  // Decrypt the converted ciphertext.
  ASSERT_OK_AND_ASSIGN(std::vector<ModularInt::Int> decrypted,
                       Decrypt<ModularInt>(key, switched_ciphertext));

  // Create the polynomial we expect.
  ASSERT_OK_AND_ASSIGN(
      auto plaintext_sub,
      plaintext_ntt.Substitute(kPower, ntt_params_main, mod_params_main));
  std::vector<ModularInt::Int> expected = RemoveError<ModularInt>(
      plaintext_sub.InverseNtt(ntt_params_main, mod_params_main),
      mod_params_main->modulus, kPlaintextModulus, mod_params_main);

  EXPECT_EQ(decrypted, expected);
}

INSTANTIATE_TEST_SUITE_P(ParameterizedTest, AuxModRelinearizationKeyTest,
                         ::testing::Values(PRNG_TYPE_CHACHA, PRNG_TYPE_HKDF));

}  // namespace
}  // namespace rlwe
