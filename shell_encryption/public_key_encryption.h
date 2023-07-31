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

#ifndef RLWE_PUBLIC_KEY_ENCRYPTION_H_
#define RLWE_PUBLIC_KEY_ENCRYPTION_H_

#include <utility>
#include <vector>

#include "absl/status/status.h"
#include "shell_encryption/error_params.h"
#include "shell_encryption/polynomial.h"
#include "shell_encryption/prng/prng.h"
#include "shell_encryption/serialization.pb.h"
#include "shell_encryption/status_macros.h"
#include "shell_encryption/symmetric_encryption.h"

namespace rlwe {

// This file implements public key encryption for the BGV scheme as specified
// in "Fully Homomorphic Encryption from Ring-LWE and Security for Key
// Dependent Messages" by Zvika Brakerski and Vinod Vaikuntanathan.
// http://www.wisdom.weizmann.ac.il/~zvikab/localpapers/IdealHom.pdf
//
// The scheme has CPA security under the hardness of the
// Ring-Learning with Errors problem (see reference above for details). We do
// not implement protections against timing attacks.

// Holds a public key pair (b = a*s + t*e, -a) as an symmetric encryption of 0.
//
// This class is thread-safe.
template <typename ModularInt>
class PublicRlweKey {
  using ModularIntParams = typename ModularInt::Params;

 public:
  // Allow copy and move, disallow copy-assign and move-assign.
  PublicRlweKey(const PublicRlweKey&) = default;
  PublicRlweKey& operator=(const PublicRlweKey&) = delete;
  PublicRlweKey(PublicRlweKey&&) = default;
  PublicRlweKey& operator=(PublicRlweKey&&) = delete;
  ~PublicRlweKey() = default;

  // Generate a public key (b = a*s + t*e, -a) derived from the given secret
  // key, where the randomness a is freshly sampled uniform over the key's
  // modulus, and the error term e has coefficients sampled from a centered
  // binomial distribution of the given variance.
  static rlwe::StatusOr<PublicRlweKey> Create(
      const SymmetricRlweKey<ModularInt>& secret_key, int variance,
      PrngType prng_type);

  // Encrypt the plaintext polynomial using the public key. Assume the public
  // key is (b, a), then the public key encryption of a plaintext polynomial m
  // is (v * b + t * e0 + m, v * a + t * e1), where v is the encryption
  // randomness with coefficients sampled from the error distribution of the
  // given variance, and e0, e1 are also sampled from the same error
  // distribution.
  //
  // This ciphertext (c0, c1) can be decrypted using the underlying secret key
  // (1, s), in the same way as a ciphertext encrypted using the secret key.
  // This is because the public key (b, a) is an symmetric encryption of 0 with
  // error te = <(1, s), (b, a)>, and so decrypting the ciphertext (c0, c1) we
  // get <(1, s), (c0, c1)> = t * (v * e + e0 + s * e1) + m, where the noise is
  // small and can be removed.
  //
  // The ciphertext type of public key encryption is same as in symmetric key
  // encryption, so we reuse the SymmetricRlweCiphertext type.
  rlwe::StatusOr<SymmetricRlweCiphertext<ModularInt>> Encrypt(
      const Polynomial<ModularInt>& plaintext, int variance,
      const ErrorParams<ModularInt>& error_params, SecurePrng* prng) const;

  // Return a SerializedPublicRlweKey containing a serialized representation of
  // the "b" polynomial in the public key and a PRNG seed used to sample the
  // random "a" polynomial in the public key.
  rlwe::StatusOr<SerializedPublicRlweKey> Serialize() const;

  // Return a PublicRlweKey represented as in `serialized`, which should contain
  // the "b" polynomial and a PRNG seed used to sample the "a" polynomial.
  // Crashes for non-valid input parameters.
  static rlwe::StatusOr<PublicRlweKey> Deserialize(
      const SerializedPublicRlweKey& serialized, ModularInt t_mod,
      const ModularIntParams* mod_params,
      const NttParameters<ModularInt>* ntt_params);

  // Accessors to the key components.
  const Polynomial<ModularInt>& GetA() const { return key_a_; }
  const Polynomial<ModularInt>& GetB() const { return key_b_; }
  // The number of coefficients in each key component.
  const unsigned int Len() const { return key_a_.Len(); }

  const ModularInt& PlaintextModulus() const { return t_mod_; }
  const NttParameters<ModularInt>& NttParams() const { return ntt_params_; }
  const ModularIntParams& ModulusParams() const { return mod_params_; }

 private:
  // Construct a public key from polynomials a and b
  PublicRlweKey(Polynomial<ModularInt> a, Polynomial<ModularInt> b,
                absl::string_view prng_seed, PrngType prng_type,
                ModularInt t_mod, const NttParameters<ModularInt>* ntt_params,
                const ModularIntParams* mod_params);

  // The two components in a public key pair (b = a*s + t*e, -a).
  Polynomial<ModularInt> key_a_;
  Polynomial<ModularInt> key_b_;

  // Prng seed used to sample the pseudorandom polynomial key_a_
  std::string prng_seed_;
  // Prng type
  PrngType prng_type_;

  // The plaintext modulus.
  ModularInt t_mod_;

  // The NTT and Montgomery parameters used to create the key polynomials.
  const NttParameters<ModularInt>& ntt_params_;
  const ModularIntParams& mod_params_;
};

}  // namespace rlwe

#endif  // RLWE_PUBLIC_KEY_ENCRYPTION_H_
