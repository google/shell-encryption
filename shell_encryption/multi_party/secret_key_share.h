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

#ifndef RLWE_MULTI_PARTY_SECRET_KEY_SHARE_H_
#define RLWE_MULTI_PARTY_SECRET_KEY_SHARE_H_

#include <vector>

#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "absl/types/span.h"
#include "shell_encryption/prng/prng.h"
#include "shell_encryption/rns/rns_context.h"
#include "shell_encryption/rns/rns_modulus.h"
#include "shell_encryption/rns/rns_polynomial.h"
#include "shell_encryption/sampler/discrete_gaussian.h"

namespace rlwe {
namespace multi_party {

// Secret key share of the RLWE-based multi-party homomorphic encryption scheme.
// Every secret key share is held by a party in the protocol, and it is used to
// partially decrypt a ciphertext encrypted under the multi-party public key.
template <typename ModularInt>
class SecretKeyShare {
 public:
  using Integer = typename ModularInt::Int;

  // Samples a secret key share from the uniform ternary distribution, wrt the
  // RNS `moduli`.
  static absl::StatusOr<SecretKeyShare> Sample(
      const RnsContext<ModularInt>* rns_context, SecurePrng* prng);

  // Returns the partial decryption contribution from this secret key share
  // holder, for a ciphertext whose "a" component is given.
  absl::StatusOr<RnsPolynomial<ModularInt>> PartialDecrypt(
      const RnsPolynomial<ModularInt>& ciphertext_component_a, double s_flood,
      const DiscreteGaussianSampler<Integer>* dg_sampler,
      SecurePrng* prng) const;

  // Accessors
  int LogN() const { return key_.LogN(); }
  int NumCoeffs() const { return key_.NumCoeffs(); }
  int NumModuli() const { return moduli_.size(); }

  // Accessor for the prime moduli chain.
  absl::Span<const PrimeModulus<ModularInt>* const> Moduli() const {
    return moduli_;
  }

  // Accessor for the key polynomial
  const RnsPolynomial<ModularInt>& Key() const { return key_; }

  // For Rust interoperability. Defined in the wrapper library, not here.
  friend class SecretKeyShareRawFactory;

 private:
  explicit SecretKeyShare(RnsPolynomial<ModularInt> key,
                          std::vector<const PrimeModulus<ModularInt>*> moduli)
      : key_(std::move(key)), moduli_(std::move(moduli)) {}

  // The key polynomial
  RnsPolynomial<ModularInt> key_;

  // The prime moduli constituting the modulus of this ciphertext.
  std::vector<const PrimeModulus<ModularInt>*> moduli_;
};

}  // namespace multi_party
}  // namespace rlwe

#endif  // RLWE_MULTI_PARTY_SECRET_KEY_SHARE_H_
