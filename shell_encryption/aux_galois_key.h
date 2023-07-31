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

#ifndef RLWE_AUX_GALOIS_KEY_H_
#define RLWE_AUX_GALOIS_KEY_H_

#include <memory>
#include <vector>

#include "shell_encryption/aux_relinearization_key.h"
#include "shell_encryption/montgomery.h"
#include "shell_encryption/ntt_parameters.h"
#include "shell_encryption/polynomial.h"
#include "shell_encryption/prng/prng.h"
#include "shell_encryption/statusor.h"
#include "shell_encryption/symmetric_encryption.h"

namespace rlwe {

// Implementation of Galois keys using the auxiliary modulus technique suggested
// in "Homomorphic Evaluation of the AES Circuit", by Craig Gentry, Shai Halevi,
// and Nigel Smart, https://eprint.iacr.org/2012/099.
//
// A Galois key is a special type of relinearization (aka key-switching) key
// that can transform a ciphertext encrypting m under a secret key (1, s(X^k))
// to a ciphertext encrypting the same plaintext m under the canonical secret
// key (1, s(X)), where k is a "substitution power". The auxiliary modulus p is
// used to reduce the size of error introduced during the key-switching process.
// For details please see document in `AuxModRelinearizationKey`.
template <typename ModularInt>
class AuxModGaloisKey {
 public:
  using ModularIntParams = typename ModularInt::Params;

  // AuxModGaloisKey is copyable and movable.
  AuxModGaloisKey(const AuxModGaloisKey&) = default;
  AuxModGaloisKey& operator=(const AuxModGaloisKey&) = default;
  AuxModGaloisKey(AuxModGaloisKey&&) = default;
  AuxModGaloisKey& operator=(AuxModGaloisKey&&) = default;
  ~AuxModGaloisKey() = default;

  // Generates a AuxModGaloisKey that can switch a ciphertext of 2 components
  // encrypted under a secret key (1, s(X^k)) to a ciphertext encrypted under
  // `secret_key` (1, s(X)), where k = `substitution_power`.
  static StatusOr<AuxModGaloisKey> Create(
      const SymmetricRlweKey<ModularInt>& secret_key, PrngType prng_type,
      int substitution_power, const ModularIntParams* mod_params_aux,
      const NttParameters<ModularInt>* ntt_params_aux) {
    if (mod_params_aux == nullptr) {
      return absl::InvalidArgumentError("mod_params_aux must not be null.");
    }
    if (ntt_params_aux == nullptr) {
      return absl::InvalidArgumentError("ntt_params_aux must not be null.");
    }
    // A Galois key has degree 1, corresponding to ciphertexts of 2 components.
    RLWE_ASSIGN_OR_RETURN(
        auto relinearization_key,
        AuxModRelinearizationKey<ModularInt>::Create(
            secret_key, prng_type, /*degree=*/1, mod_params_aux, ntt_params_aux,
            substitution_power));
    return AuxModGaloisKey(std::move(relinearization_key));
  }

  // Applies the Galois key to a ciphertext of 2 components encrypted under key
  // (1, s(X^k)) and returns the resulting ciphertext of 2 components encrypting
  // the same plaintext under the canonical secret key (1, s(X)), where k is
  // `SubstitutionPower`.
  StatusOr<SymmetricRlweCiphertext<ModularInt>> ApplyTo(
      const SymmetricRlweCiphertext<ModularInt>& ciphertext) const {
    if (ciphertext.PowerOfS() != SubstitutionPower()) {
      return absl::InvalidArgumentError(absl::StrCat(
          "Ciphertext PowerOfS: ", ciphertext.PowerOfS(),
          " doesn't match the key substitution power: ", SubstitutionPower()));
    }
    return relinearization_key_.ApplyTo(ciphertext);
  }

  // Returns a SerializedAuxModGaloisKey containing a representation of the
  // underlying relinearization key.
  StatusOr<SerializedAuxModGaloisKey> Serialize() const {
    SerializedAuxModGaloisKey output;
    RLWE_ASSIGN_OR_RETURN(*output.mutable_key(),
                          relinearization_key_.Serialize());
    return output;
  }

  // Returns an AuxModGaloisKey represented as in `serialized`, which should
  // contain the underlying relinearization key.
  static StatusOr<AuxModGaloisKey> Deserialize(
      const SerializedAuxModGaloisKey& serialized,
      const ModularIntParams* mod_params_main,
      const NttParameters<ModularInt>* ntt_params_main,
      const ModularIntParams* mod_params_aux,
      const NttParameters<ModularInt>* ntt_params_aux,
      typename ModularInt::Int t) {
    if (mod_params_main == nullptr) {
      return absl::InvalidArgumentError("mod_params_main must not be null.");
    }
    if (ntt_params_main == nullptr) {
      return absl::InvalidArgumentError("ntt_params_main must not be null.");
    }
    if (mod_params_aux == nullptr) {
      return absl::InvalidArgumentError("mod_params_aux must not be null.");
    }
    if (ntt_params_aux == nullptr) {
      return absl::InvalidArgumentError("ntt_params_aux must not be null.");
    }

    RLWE_ASSIGN_OR_RETURN(
        AuxModRelinearizationKey<ModularInt> key,
        AuxModRelinearizationKey<ModularInt>::Deserialize(
            serialized.key(), mod_params_main, ntt_params_main, mod_params_aux,
            ntt_params_aux, t));
    return AuxModGaloisKey(std::move(key));
  }

  // Accessors.
  int SubstitutionPower() const {
    return relinearization_key_.SubstitutionPower();
  }

 private:
  explicit AuxModGaloisKey(
      AuxModRelinearizationKey<ModularInt> relinearization_key)
      : relinearization_key_(std::move(relinearization_key)) {}

  // The underlying key-switching / relinearization key.
  AuxModRelinearizationKey<ModularInt> relinearization_key_;
};

}  // namespace rlwe

#endif  // RLWE_AUX_GALOIS_KEY_H_
