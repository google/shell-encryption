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

#ifndef RLWE_MULTI_PARTY_PUBLIC_PARAMETER_H_
#define RLWE_MULTI_PARTY_PUBLIC_PARAMETER_H_

#include <memory>
#include <string>
#include <vector>

#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "absl/types/span.h"
#include "shell_encryption/multi_party/serialization.pb.h"
#include "shell_encryption/prng/prng.h"
#include "shell_encryption/rns/rns_context.h"
#include "shell_encryption/rns/rns_modulus.h"
#include "shell_encryption/rns/rns_polynomial.h"
#include "shell_encryption/status_macros.h"

namespace rlwe {
namespace multi_party {

// This class stores the public parameter of multi-party additive homomorphic
// encryption protocol, which is used by all parties to generate their public
// key shares.
template <typename ModularInt>
class PublicParameter {
 public:
  // Deterministic factory function to create a public parameter from a seed.
  static absl::StatusOr<std::unique_ptr<const PublicParameter>> CreateFromSeed(
      const RnsContext<ModularInt>* rns_context, int error_variance,
      std::string prng_seed, PrngType prng_type);

  // Factory function to create a fresh public parameter for the parties to
  // generate their public key shares.
  static absl::StatusOr<std::unique_ptr<const PublicParameter>> Create(
      const RnsContext<ModularInt>* rns_context, int error_variance,
      PrngType prng_type);

  static absl::StatusOr<std::unique_ptr<PublicParameter>> Deserialize(
      const SerializedPublicParameter& serialized,
      const RnsContext<ModularInt>* rns_context);

  absl::StatusOr<SerializedPublicParameter> Serialize() const {
    SerializedPublicParameter serialized;
    serialized.set_prng_seed(prng_seed_);
    serialized.set_prng_type(prng_type_);
    serialized.set_error_variance(error_variance_);
    return serialized;
  }

  // Accessors.
  int ErrorVariance() const { return error_variance_; }

  const RnsPolynomial<ModularInt>& PublicKeyComponentA() const {
    return key_a_;
  }

  absl::Span<const PrimeModulus<ModularInt>* const> Moduli() const {
    return moduli_;
  }

  int LogN() const { return key_a_.LogN(); }
  int NumCoeffs() const { return key_a_.NumCoeffs(); }
  int NumModuli() const { return moduli_.size(); }

 private:
  explicit PublicParameter(std::string prng_seed, PrngType prng_type,
                           int error_variance, RnsPolynomial<ModularInt> key_a,
                           std::vector<const PrimeModulus<ModularInt>*> moduli)
      : prng_seed_(std::move(prng_seed)),
        prng_type_(prng_type),
        error_variance_(error_variance),
        key_a_(std::move(key_a)),
        moduli_(std::move(moduli)) {}

  // PRNG seed and type for sampling the random polynomial `key_a_`.
  const std::string prng_seed_;
  const PrngType prng_type_;

  // The variance for generating the public key and for encrypting using the
  // public key.
  const int error_variance_;

  // The "a" component of the public key.
  const RnsPolynomial<ModularInt> key_a_;

  // The RNS moduli used by the public key.
  const std::vector<const PrimeModulus<ModularInt>*> moduli_;
};

}  // namespace multi_party
}  // namespace rlwe

#endif  // RLWE_MULTI_PARTY_PUBLIC_PARAMETER_H_
