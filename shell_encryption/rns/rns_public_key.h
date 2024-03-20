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

#ifndef RLWE_RNS_RNS_PUBLIC_KEY_H_
#define RLWE_RNS_RNS_PUBLIC_KEY_H_

#include <string>
#include <vector>

#include "absl/status/statusor.h"
#include "absl/strings/string_view.h"
#include "absl/types/span.h"
#include "shell_encryption/rns/rns_modulus.h"
#include "shell_encryption/rns/rns_polynomial.h"
#include "shell_encryption/rns/rns_secret_key.h"

namespace rlwe {

template <typename ModularInt>
class RnsRlwePublicKey {
 public:
  using Integer = typename ModularInt::Int;

  // Accessors.
  int LogN() const { return key_a_.LogN(); }

  int Level() const { return moduli_.size() - 1; }

  int NumCoeffs() const { return key_a_.NumCoeffs(); }

  int NumModuli() const { return moduli_.size(); }

  // Accessor for the prime moduli chain.
  absl::Span<const PrimeModulus<ModularInt>* const> Moduli() const {
    return moduli_;
  }

  // Accessors to the two components in a RLWE public key.
  const RnsPolynomial<ModularInt>& KeyA() const { return key_a_; }
  const RnsPolynomial<ModularInt>& KeyB() const { return key_b_; }

 protected:
  // Generate a public key (b = a*s + scale*e, -a) derived from the given secret
  // key, where the randomness a is freshly sampled uniformly random polynomial,
  // and the error term e has coefficients sampled from a centered binomial
  // distribution of the given variance. The scaling factor `error_scalar`
  // should be set to the plaintext modulus for BGV public key, and it should be
  // 1 in other schemes.
  static absl::StatusOr<RnsRlwePublicKey> Create(
      const RnsRlweSecretKey<ModularInt>& secret_key, int variance,
      PrngType prng_type, Integer error_scalar = 1);

  // Allow derived public key classes to access data members.
  const std::vector<const PrimeModulus<ModularInt>*>& moduli() const {
    return moduli_;
  }

  int variance() const { return variance_; }

 private:
  explicit RnsRlwePublicKey(RnsPolynomial<ModularInt> key_a,
                            RnsPolynomial<ModularInt> key_b,
                            std::vector<const PrimeModulus<ModularInt>*> moduli,
                            int variance, absl::string_view prng_seed,
                            PrngType prng_type)
      : key_a_(std::move(key_a)),
        key_b_(std::move(key_b)),
        moduli_(std::move(moduli)),
        variance_(variance),
        prng_seed_(std::string(prng_seed)),
        prng_type_(prng_type) {}

  // The two components (b, a) of a RLWE public key.
  RnsPolynomial<ModularInt> key_a_;
  RnsPolynomial<ModularInt> key_b_;

  // The RNS moduli used to construct this public key.
  std::vector<const PrimeModulus<ModularInt>*> moduli_;

  // The variance of the binomial distribution for the error terms in the public
  // key and the public key encryptions.
  int variance_;

  // Prng seed used to sample the pseudorandom polynomial key_a_
  std::string prng_seed_;
  // Prng type
  PrngType prng_type_;
};

}  // namespace rlwe

#endif  // RLWE_RNS_RNS_PUBLIC_KEY_H_
