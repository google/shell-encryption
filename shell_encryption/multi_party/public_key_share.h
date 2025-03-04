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

#ifndef RLWE_MULTI_PARTY_PUBLIC_KEY_SHARE_H_
#define RLWE_MULTI_PARTY_PUBLIC_KEY_SHARE_H_

#include <vector>

#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "absl/types/span.h"
#include "shell_encryption/multi_party/public_parameter.h"
#include "shell_encryption/multi_party/secret_key_share.h"
#include "shell_encryption/prng/prng.h"
#include "shell_encryption/rns/rns_modulus.h"
#include "shell_encryption/rns/rns_polynomial.h"
#include "shell_encryption/status_macros.h"

namespace rlwe {
namespace multi_party {

// Public key share of the RLWE-based multi-party homomorphic encryption scheme.
// Every secret key share holder derives its own public key share, and by
// combining all public key shares we get the public key of the multi-party
// scheme.
//
// In the current protocol, a public key share is a polynomial b = -a * s + e,
// where `a` is the random "a" component of the public key and is given in the
// public parameter, and `e` is a fresh error polynomial.
template <typename ModularInt>
class PublicKeyShare {
 public:
  // Returns a public key share based on the given `secret_key_share`.
  static absl::StatusOr<PublicKeyShare> Create(
      const SecretKeyShare<ModularInt>* secret_key_share,
      const PublicParameter<ModularInt>* public_parameter, PrngType prng_type);

  // Creates a new public key share derived from `secret_key_share` and the "a"
  // component of the public key as given in `public_parameter`, and populates
  // `key_b` with the raw key share polynomial in NTT form. Randomness is
  // sampled using `prng`.
  // `key_error` is optional and can be null. When it is non-null, `key_error`
  // is populated with the error polynomial in NTT form such that
  //    key_b = -a(X) * secret_key_share + key_error \in Z[X]/(Q, X^N+1),
  // where a(X) is the "a" component of the public key.
  // `wrap_around` is optional and can be null. When it is non-null,
  // `wrap_around` is populated in coefficient form such that
  //    key_b = -a(X) * secret_key_share + key_error + wrap_around * (X^N + 1)
  // over the ring Z[X]/(Q), i.e. the RHS is computed without reduction modulo
  // X^N + 1, and the LHS is defined in the first equation above.
  //
  // The optional `key_error` and `wrap_around` are useful to generate ZK proofs
  // about valid public key shares.
  //
  // Note: When `wrap_around` is non-null, the parameters must be NTT-friendly
  //       wrt to the larger cyclotomic X^{2N} + 1, i.e. N is a power of 2 and
  //       the modulus Q must be a product of primes q such that 4N factors q-1.
  static absl::Status CreateExplicit(
      const RnsPolynomial<ModularInt>& secret_key_share,
      const PublicParameter<ModularInt>* public_parameter, SecurePrng* prng,
      RnsPolynomial<ModularInt>* key_b, RnsPolynomial<ModularInt>* key_error,
      RnsPolynomial<ModularInt>* wrap_around);

  static absl::StatusOr<PublicKeyShare> Deserialize(
      const SerializedPublicKeyShare& serialized,
      absl::Span<const PrimeModulus<ModularInt>* const> moduli);

  absl::StatusOr<SerializedPublicKeyShare> Serialize() const {
    SerializedPublicKeyShare serialized;
    RLWE_ASSIGN_OR_RETURN(*serialized.mutable_key_b(),
                          key_b_.Serialize(moduli_));
    return serialized;
  }

  // Accessor to the "b" component in a public key share.
  const RnsPolynomial<ModularInt>& ComponentB() const { return key_b_; }

 private:
  explicit PublicKeyShare(RnsPolynomial<ModularInt> key_b,
                          std::vector<const PrimeModulus<ModularInt>*> moduli)
      : key_b_(std::move(key_b)), moduli_(std::move(moduli)) {}

  // The "b" component of the public key share.
  RnsPolynomial<ModularInt> key_b_;

  // The prime moduli constituting the modulus of this ciphertext.
  std::vector<const PrimeModulus<ModularInt>*> moduli_;
};

}  // namespace multi_party
}  // namespace rlwe

#endif  // RLWE_MULTI_PARTY_PUBLIC_KEY_SHARE_H_
