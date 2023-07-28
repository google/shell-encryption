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

#include "shell_encryption/aux_galois_key.h"

#include <memory>
#include <random>

#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include "absl/status/status.h"
#include "shell_encryption/constants.h"
#include "shell_encryption/montgomery.h"
#include "shell_encryption/ntt_parameters.h"
#include "shell_encryption/polynomial.h"
#include "shell_encryption/serialization.pb.h"
#include "shell_encryption/status_macros.h"
#include "shell_encryption/symmetric_encryption.h"
#include "shell_encryption/testing/status_matchers.h"
#include "shell_encryption/testing/status_testing.h"
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
using AuxModGaloisKey = AuxModGaloisKey<ModularInt>;
using ErrorParams = ErrorParams<ModularInt>;
using ::rlwe::testing::StatusIs;
using ::testing::HasSubstr;

// Set constants.
constexpr int kLogPlaintextModulus = 1;
constexpr uint64_t kPlaintextModulus = (1ULL << kLogPlaintextModulus) + 1;  // t
constexpr int kDefaultVariance = 8;
constexpr int kLogCoeffs = 10;
constexpr int kCoeffs = 1 << kLogCoeffs;  // 1024

// We use moduli q and p that both support NTT in the ring
// Z[X]/(X^1024 + 1, q) and Z[X]/(X^1024 + 1, p), where q is the main modulus
// and p is the auxiliary modulus.
// Note that these moduli do not offer sufficient security as the lattice
// dimension 1024 is too small for the modulus q*p; they are chosen only to
// test correctness.
constexpr uint64_t kMainModulus = kModulus25;         // q = 33538049
constexpr uint64_t kAuxiliaryModulus = 134215681ULL;  // 26-bit
// Seed for the RNG
constexpr unsigned int kRandSeed = 1;

class AuxModGaloisKeyTest : public ::testing::TestWithParam<PrngType> {
 public:
  AuxModGaloisKeyTest() : prng_type_(GetParam()), mt_rand_(kRandSeed) {}

 protected:
  void SetUp() override {
    ASSERT_OK_AND_ASSIGN(mod_params_main_,
                         ModularInt::Params::Create(kMainModulus));
    ASSERT_OK_AND_ASSIGN(mod_params_aux_,
                         ModularInt::Params::Create(kAuxiliaryModulus));
    ASSERT_OK_AND_ASSIGN(auto ntt_params_main,
                         InitializeNttParameters<ModularInt>(
                             kLogCoeffs, mod_params_main_.get()));
    ASSERT_OK_AND_ASSIGN(
        auto ntt_params_aux,
        InitializeNttParameters<ModularInt>(kLogCoeffs, mod_params_aux_.get()));

    ntt_params_main_ =
        std::make_unique<const NttParameters>(std::move(ntt_params_main));
    ntt_params_aux_ =
        std::make_unique<const NttParameters>(std::move(ntt_params_aux));

    ASSERT_OK_AND_ASSIGN(
        auto error_params_small,
        ErrorParams::Create(kLogPlaintextModulus, kDefaultVariance,
                            mod_params_main_.get(), ntt_params_main_.get()));
    error_params_ = std::make_unique<const ErrorParams>(error_params_small);
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
      const ModularInt::Params* mod_params) {
    std::vector<ModularInt> output(coeffs.size(),
                                   ModularInt::ImportZero(mod_params));
    for (int i = 0; i < output.size(); ++i) {
      RLWE_ASSIGN_OR_RETURN(output[i],
                            ModularInt::ImportInt(coeffs[i], mod_params));
    }
    return output;
  }

  // Sample a secret key.
  StatusOr<Key> SampleKey() {
    RLWE_ASSIGN_OR_RETURN(std::string prng_seed, GenerateSeed());
    RLWE_ASSIGN_OR_RETURN(auto prng, CreatePrng(prng_seed));
    return Key::Sample(kLogCoeffs, kDefaultVariance, kLogPlaintextModulus,
                       mod_params_main_.get(), ntt_params_main_.get(),
                       prng.get());
  }

  // Sample a random plaintext and convert it to NTT form wrt the main modulus
  StatusOr<Polynomial> SamplePlaintextNtt() {
    std::vector<ModularInt::Int> plaintext(kCoeffs);
    for (int i = 0; i < kCoeffs; ++i) {
      plaintext[i] = mt_rand_() % kPlaintextModulus;
    }
    RLWE_ASSIGN_OR_RETURN(
        auto plaintext_m,
        ConvertToMontgomery(plaintext, mod_params_main_.get()));
    return Polynomial::ConvertToNtt(plaintext_m, ntt_params_main_.get(),
                                    mod_params_main_.get());
  }

  // Encrypt a plaintext.
  StatusOr<Ciphertext> Encrypt(const Key& key,
                               const Polynomial& plaintext_ntt) {
    RLWE_ASSIGN_OR_RETURN(std::string prng_seed, GenerateSeed());
    RLWE_ASSIGN_OR_RETURN(auto prng, CreatePrng(prng_seed));
    return ::rlwe::Encrypt<ModularInt>(key, plaintext_ntt, error_params_.get(),
                                       prng.get());
  }

  // Accessors
  const ModularIntParams* ModParamsMain() const {
    return mod_params_main_.get();
  }
  const NttParameters* NttParamsMain() const { return ntt_params_main_.get(); }
  const ModularIntParams* ModParamsAux() const { return mod_params_aux_.get(); }
  const NttParameters* NttParamsAux() const { return ntt_params_aux_.get(); }

  PrngType prng_type_;
  std::mt19937 mt_rand_;
  std::unique_ptr<const ModularIntParams> mod_params_main_;
  std::unique_ptr<const ModularIntParams> mod_params_aux_;
  std::unique_ptr<const NttParameters> ntt_params_main_;
  std::unique_ptr<const NttParameters> ntt_params_aux_;
  std::unique_ptr<const ErrorParams> error_params_;
};

// Checks for error handling when applying the Galois key on a ciphertext with
// a mismatch substitution power.
TEST_P(AuxModGaloisKeyTest, ErrorIfSubstitutionPowerMismatch) {
  ASSERT_OK_AND_ASSIGN(auto key, SampleKey());
  // Create a Galois key with substitution power 5.
  constexpr int k_power = 5;
  ASSERT_OK_AND_ASSIGN(
      auto galois_key,
      AuxModGaloisKey::Create(key, prng_type_, /*substitution_power=*/k_power,
                              ModParamsAux(), NttParamsAux()));
  ASSERT_EQ(galois_key.SubstitutionPower(), k_power);

  // Encrypt a plaintext which has default substitution power 1.
  ASSERT_OK_AND_ASSIGN(auto plaintext_ntt, SamplePlaintextNtt());
  ASSERT_OK_AND_ASSIGN(auto ciphertext, Encrypt(key, plaintext_ntt));
  ASSERT_EQ(ciphertext.PowerOfS(), 1);

  // Apply the Galois key to the ciphertext should result in an error.
  EXPECT_THAT(galois_key.ApplyTo(ciphertext),
              StatusIs(::absl::StatusCode::kInvalidArgument,
                       HasSubstr("Ciphertext PowerOfS")));
}

// Checks for error handling when applying the Galois key on a ciphertext with
// more than two components.
TEST_P(AuxModGaloisKeyTest, ErrorIfCiphertextNumberOfComponentsMismatch) {
  constexpr int k_power = 25;
  ASSERT_OK_AND_ASSIGN(auto key, SampleKey());
  ASSERT_OK_AND_ASSIGN(
      auto galois_key,
      AuxModGaloisKey::Create(key, prng_type_, /*substitution_power=*/k_power,
                              ModParamsAux(), NttParamsAux()));

  // Encrypt two plaintexts and compute their product.
  ASSERT_OK_AND_ASSIGN(auto plaintext1_ntt, SamplePlaintextNtt());
  ASSERT_OK_AND_ASSIGN(auto plaintext2_ntt, SamplePlaintextNtt());
  ASSERT_OK_AND_ASSIGN(auto ciphertext1, Encrypt(key, plaintext1_ntt));
  ASSERT_OK_AND_ASSIGN(auto ciphertext2, Encrypt(key, plaintext2_ntt));
  ASSERT_OK_AND_ASSIGN(auto product, ciphertext1* ciphertext2);
  // Substitute X with X^power in the product.
  ASSERT_OK_AND_ASSIGN(auto product_sub,
                       product.Substitute(k_power, NttParamsMain()));
  // product_sub has the desired PowerOfS.
  ASSERT_EQ(product_sub.PowerOfS(), k_power);
  // product_sub has three components (c0, c1, c2) and it is not compatible with
  // our Galois key.
  ASSERT_EQ(product_sub.Len(), 3);

  // Apply the Galois key to the substituted product ciphertext and it should
  // result in an error.
  EXPECT_THAT(
      galois_key.ApplyTo(product_sub),
      StatusIs(::absl::StatusCode::kInvalidArgument,
               HasSubstr("Ciphertext has incompatible number of components")));
}

// Checks for error handling when creating a relinearization key with a null
// mod_params_main parameter.
TEST_P(AuxModGaloisKeyTest, ErrorIfModParamsAuxIsNull) {
  ASSERT_OK_AND_ASSIGN(auto key, SampleKey());
  EXPECT_THAT(
      AuxModGaloisKey::Create(key, prng_type_,
                              /*substitution_power=*/5,
                              /*mod_params_aux=*/nullptr, NttParamsAux()),
      StatusIs(absl::StatusCode::kInvalidArgument,
               HasSubstr("mod_params_aux must not be null")));
}

// Checks for error handling when creating a relinearization key with a null
// ntt_params_main parameter.
TEST_P(AuxModGaloisKeyTest, ErrorIfNttParamsAuxIsNull) {
  ASSERT_OK_AND_ASSIGN(auto key, SampleKey());
  EXPECT_THAT(AuxModGaloisKey::Create(key, prng_type_,
                                      /*degree=*/2, ModParamsAux(),
                                      /*ntt_params_aux=*/nullptr),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("ntt_params_aux must not be null")));
}

TEST_P(AuxModGaloisKeyTest, KeySwitchedCiphertextDecrypts) {
  constexpr int k_power = 5;
  ASSERT_OK_AND_ASSIGN(auto key, SampleKey());
  ASSERT_OK_AND_ASSIGN(
      auto galois_key,
      AuxModGaloisKey::Create(key, prng_type_, /*substitution_power=*/k_power,
                              ModParamsAux(), NttParamsAux()));

  // Encrypt a plaintext m(X).
  ASSERT_OK_AND_ASSIGN(auto plaintext_ntt, SamplePlaintextNtt());
  ASSERT_OK_AND_ASSIGN(auto ciphertext, Encrypt(key, plaintext_ntt));
  // Substitute X with X^power in the ciphertext.
  ASSERT_OK_AND_ASSIGN(auto ciphertext_sub,
                       ciphertext.Substitute(k_power, NttParamsMain()));
  ASSERT_EQ(ciphertext_sub.PowerOfS(), k_power);
  // Apply the Galois key to ciphertext_sub and decrypt.
  ASSERT_OK_AND_ASSIGN(auto ciphertext_switched,
                       galois_key.ApplyTo(ciphertext_sub));
  // The switched ciphertext should be in standard form (c0, c1) with power
  // of the secret key being 1, and wrt the main modulus.
  ASSERT_EQ(ciphertext_switched.PowerOfS(), 1);
  ASSERT_EQ(ciphertext_switched.Len(), 2);
  ASSERT_EQ(ciphertext_switched.ModulusParams()->modulus,
            mod_params_main_->modulus);

  // Decrypt the switched ciphertext.
  ASSERT_OK_AND_ASSIGN(std::vector<ModularInt::Int> decrypted,
                       Decrypt<ModularInt>(key, ciphertext_switched));

  // Create the polynomial m(X^j) we expect.
  ASSERT_OK_AND_ASSIGN(
      auto plaintext_sub,
      plaintext_ntt.Substitute(k_power, NttParamsMain(), ModParamsMain()));
  std::vector<ModularInt::Int> expected = RemoveError<ModularInt>(
      plaintext_sub.InverseNtt(NttParamsMain(), ModParamsMain()),
      ModParamsMain()->modulus, kPlaintextModulus, ModParamsMain());

  EXPECT_EQ(decrypted, expected);
}

// Checks that applying two Galois keys with substitution powers a, b, resp,
// is equivalent to substituting X with X^(ab) on plaintext.
TEST_P(AuxModGaloisKeyTest, ComposedKeySwitchingDecrypts) {
  constexpr int k_galois_power0 = 5;  // The first power a
  constexpr int k_galois_power1 = 7;  // The second power b
  constexpr int k_substitution_power = k_galois_power0 * k_galois_power1;

  ASSERT_OK_AND_ASSIGN(auto key, SampleKey());
  ASSERT_OK_AND_ASSIGN(auto galois_key0,
                       AuxModGaloisKey::Create(key, prng_type_, k_galois_power0,
                                               ModParamsAux(), NttParamsAux()));
  ASSERT_OK_AND_ASSIGN(auto galois_key1,
                       AuxModGaloisKey::Create(key, prng_type_, k_galois_power1,
                                               ModParamsAux(), NttParamsAux()));
  ASSERT_EQ(galois_key0.SubstitutionPower(), k_galois_power0);
  ASSERT_EQ(galois_key1.SubstitutionPower(), k_galois_power1);

  // Encrypt a plaintext m(X).
  ASSERT_OK_AND_ASSIGN(auto plaintext_ntt, SamplePlaintextNtt());
  ASSERT_OK_AND_ASSIGN(auto ciphertext, Encrypt(key, plaintext_ntt));

  // Substitute X with X^a to get Enc(s(X^a), m(X^a)).
  ASSERT_OK_AND_ASSIGN(auto ciphertext_sub0,
                       ciphertext.Substitute(k_galois_power0, NttParamsMain()));
  ASSERT_EQ(ciphertext_sub0.PowerOfS(), k_galois_power0);
  // Apply the first Galois key to get Enc(s(X), m(X^a))
  ASSERT_OK_AND_ASSIGN(auto ciphertext_switched0,
                       galois_key0.ApplyTo(ciphertext_sub0));
  ASSERT_EQ(ciphertext_switched0.PowerOfS(), 1);
  ASSERT_EQ(ciphertext_switched0.Len(), 2);
  ASSERT_EQ(ciphertext_switched0.ModulusParams()->modulus,
            ModParamsMain()->modulus);

  // Substitute X with X^b to get Enc(s(X^b), m(X^(ab)))
  ASSERT_OK_AND_ASSIGN(
      auto ciphertext_sub1,
      ciphertext_switched0.Substitute(k_galois_power1, NttParamsMain()));
  ASSERT_EQ(ciphertext_sub1.PowerOfS(), k_galois_power1);
  // Apply the second Galois key to get Enc(s(X), m(X^(ab)))
  ASSERT_OK_AND_ASSIGN(auto ciphertext_switched1,
                       galois_key1.ApplyTo(ciphertext_sub1));
  ASSERT_EQ(ciphertext_switched1.PowerOfS(), 1);
  ASSERT_EQ(ciphertext_switched1.Len(), 2);
  ASSERT_EQ(ciphertext_switched1.ModulusParams()->modulus,
            ModParamsMain()->modulus);

  // Decrypt the final ciphertext.
  ASSERT_OK_AND_ASSIGN(std::vector<ModularInt::Int> decrypted,
                       Decrypt<ModularInt>(key, ciphertext_switched1));

  // Create the polynomial m(X^(ab)) we expect.
  ASSERT_OK_AND_ASSIGN(
      auto plaintext_sub,
      plaintext_ntt.Substitute(k_substitution_power, NttParamsMain(),
                               ModParamsMain()));
  std::vector<ModularInt::Int> expected = RemoveError<ModularInt>(
      plaintext_sub.InverseNtt(NttParamsMain(), ModParamsMain()),
      ModParamsMain()->modulus, kPlaintextModulus, ModParamsMain());

  EXPECT_EQ(decrypted, expected);
}

TEST_P(AuxModGaloisKeyTest, ErrorIfDeserializeWithNullModParamsMain) {
  ASSERT_OK_AND_ASSIGN(auto key, SampleKey());
  ASSERT_OK_AND_ASSIGN(
      auto galois_key,
      AuxModGaloisKey::Create(key, prng_type_, /*substitution_power=*/5,
                              ModParamsAux(), NttParamsAux()));

  // Serialize.
  ASSERT_OK_AND_ASSIGN(SerializedAuxModGaloisKey serialized,
                       galois_key.Serialize());
  // Deserialize with a null mod_params_main argument.
  EXPECT_THAT(AuxModGaloisKey::Deserialize(
                  serialized, /*mod_params_main=*/nullptr, NttParamsMain(),
                  ModParamsAux(), NttParamsAux(),
                  key.PlaintextModulus().ExportInt(ModParamsMain())),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("mod_params_main must not be null")));
}

TEST_P(AuxModGaloisKeyTest, ErrorIfDeserializeWithNullNttParamsMain) {
  ASSERT_OK_AND_ASSIGN(auto key, SampleKey());
  ASSERT_OK_AND_ASSIGN(
      auto galois_key,
      AuxModGaloisKey::Create(key, prng_type_, /*substitution_power=*/5,
                              ModParamsAux(), NttParamsAux()));

  // Serialize.
  ASSERT_OK_AND_ASSIGN(SerializedAuxModGaloisKey serialized,
                       galois_key.Serialize());
  // Deserialize with a null ntt_params_main argument.
  EXPECT_THAT(AuxModGaloisKey::Deserialize(
                  serialized, ModParamsMain(),
                  /*ntt_params_main=*/nullptr, ModParamsAux(), NttParamsAux(),
                  key.PlaintextModulus().ExportInt(ModParamsMain())),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("ntt_params_main must not be null")));
}

TEST_P(AuxModGaloisKeyTest, ErrorIfDeserializeWithNullModParamsAux) {
  ASSERT_OK_AND_ASSIGN(auto key, SampleKey());
  ASSERT_OK_AND_ASSIGN(
      auto galois_key,
      AuxModGaloisKey::Create(key, prng_type_, /*substitution_power=*/5,
                              ModParamsAux(), NttParamsAux()));

  // Serialize.
  ASSERT_OK_AND_ASSIGN(SerializedAuxModGaloisKey serialized,
                       galois_key.Serialize());
  // Deserialize with a null mod_params_main argument.
  EXPECT_THAT(AuxModGaloisKey::Deserialize(
                  serialized, ModParamsMain(), NttParamsMain(),
                  /*mod_params_aux=*/nullptr, NttParamsAux(),
                  key.PlaintextModulus().ExportInt(ModParamsMain())),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("mod_params_aux must not be null")));
}

TEST_P(AuxModGaloisKeyTest, ErrorIfDeserializeWithNullNttParamsAux) {
  ASSERT_OK_AND_ASSIGN(auto key, SampleKey());
  ASSERT_OK_AND_ASSIGN(
      auto galois_key,
      AuxModGaloisKey::Create(key, prng_type_, /*substitution_power=*/5,
                              ModParamsAux(), NttParamsAux()));

  // Serialize.
  ASSERT_OK_AND_ASSIGN(SerializedAuxModGaloisKey serialized,
                       galois_key.Serialize());
  // Deserialize with a null ntt_params_main argument.
  EXPECT_THAT(AuxModGaloisKey::Deserialize(
                  serialized, ModParamsMain(), NttParamsMain(), ModParamsAux(),
                  /*ntt_params_aux=*/nullptr,
                  key.PlaintextModulus().ExportInt(ModParamsMain())),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("ntt_params_aux must not be null")));
}

TEST_P(AuxModGaloisKeyTest, SerializeAuxModGaloisKey) {
  constexpr int k_power = 5;
  ASSERT_OK_AND_ASSIGN(auto key, SampleKey());
  ASSERT_OK_AND_ASSIGN(auto galois_key,
                       AuxModGaloisKey::Create(key, prng_type_,
                                               /*substitution_power=*/k_power,
                                               ModParamsAux(), NttParamsAux()));

  // Serialize and deserialize.
  ASSERT_OK_AND_ASSIGN(SerializedAuxModGaloisKey serialized,
                       galois_key.Serialize());
  ASSERT_OK_AND_ASSIGN(
      auto deserialized,
      AuxModGaloisKey::Deserialize(
          serialized, ModParamsMain(), NttParamsMain(), ModParamsAux(),
          NttParamsAux(), key.PlaintextModulus().ExportInt(ModParamsMain())));
  ASSERT_EQ(deserialized.SubstitutionPower(), k_power);

  // Encrypt a random plaintext.
  ASSERT_OK_AND_ASSIGN(Polynomial plaintext_ntt, SamplePlaintextNtt());
  ASSERT_OK_AND_ASSIGN(auto ciphertext, Encrypt(key, plaintext_ntt));
  ASSERT_OK_AND_ASSIGN(auto ciphertext_sub,
                       ciphertext.Substitute(k_power, NttParamsMain()));
  // Apply the deserialized galois key.
  ASSERT_OK_AND_ASSIGN(auto ciphertext_switched,
                       deserialized.ApplyTo(ciphertext_sub));
  // The switched ciphertext should be in standard form (c0, c1) with power
  // of the secret key being 1, and wrt the main modulus.
  ASSERT_EQ(ciphertext_switched.PowerOfS(), 1);
  ASSERT_EQ(ciphertext_switched.Len(), 2);
  ASSERT_EQ(ciphertext_switched.ModulusParams()->modulus,
            ModParamsMain()->modulus);

  // Decrypt the switched ciphertext.
  ASSERT_OK_AND_ASSIGN(std::vector<ModularInt::Int> decrypted,
                       Decrypt<ModularInt>(key, ciphertext_switched));

  // Create the polynomial we expect.
  ASSERT_OK_AND_ASSIGN(
      auto plaintext_sub,
      plaintext_ntt.Substitute(k_power, NttParamsMain(), ModParamsMain()));
  std::vector<ModularInt::Int> expected = RemoveError<ModularInt>(
      plaintext_sub.InverseNtt(NttParamsMain(), ModParamsMain()),
      ModParamsMain()->modulus, kPlaintextModulus, ModParamsMain());
  EXPECT_EQ(decrypted, expected);
}

INSTANTIATE_TEST_SUITE_P(ParameterizedTest, AuxModGaloisKeyTest,
                         ::testing::Values(PRNG_TYPE_CHACHA, PRNG_TYPE_HKDF));

}  // namespace
}  // namespace rlwe
