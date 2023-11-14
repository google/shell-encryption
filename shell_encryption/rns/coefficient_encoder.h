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

#ifndef RLWE_RNS_COEFFICIENT_ENCODER_H_
#define RLWE_RNS_COEFFICIENT_ENCODER_H_

#include <vector>

#include "absl/status/statusor.h"
#include "absl/types/span.h"
#include "shell_encryption/rns/rns_context.h"
#include "shell_encryption/rns/rns_modulus.h"
#include "shell_encryption/rns/rns_polynomial.h"

namespace rlwe {

// This class specifies how to encode messages in Z_t^N as a polynomial
// in Z[X]/(Q, X^N+1), where t is the plaintext modulus and Q is the ciphertext
// modulus. This implements the coefficient encoding for both BGV and BFV.
template <typename ModularInt>
class CoefficientEncoder {
 public:
  using Integer = typename ModularInt::Int;
  using ModularIntParams = typename ModularInt::Params;

  // Returns a CoefficientEncoder that can encode n messages in Z_t, where
  // t is the plaintext modulus defined in `context`.
  static absl::StatusOr<CoefficientEncoder<ModularInt>> Create(
      const RnsContext<ModularInt>* context);

  // Returns a polynomial (in coefficient form) whose coefficients encode
  // `messages` in the least significant bits, using balanced representation.
  // For example, if messages[i] > plaintext_modulus / 2, then the corresponding
  // coefficient is Q - (plaintext_modulus - messages[i]), where Q is product of
  // `moduli`; otherwise the corresponding coefficient is messages[i] (mod Q).
  // Note that values in `messages` must all be smaller than plaintext_modulus.
  absl::StatusOr<RnsPolynomial<ModularInt>> EncodeBgv(
      absl::Span<const Integer> messages,
      absl::Span<const PrimeModulus<ModularInt>* const> moduli) const;

  // Returns the messages encoded in the least significant bits of coefficients
  // of `noisy_plaintext`.
  absl::StatusOr<std::vector<Integer>> DecodeBgv(
      RnsPolynomial<ModularInt> noisy_plaintext,
      absl::Span<const PrimeModulus<ModularInt>* const> moduli) const;

  // Returns a polynomial (in NTT form) whose coefficients encode `messages`.
  // When `is_scaled` is set, this encoding method can be used in BFV scheme to
  // encrypt `messages`. When `is_scaled` is set to false, this encoding method
  // produces a plaintext polynomial suitable for ciphertext-plaintext
  // multiplication in BFV.
  absl::StatusOr<RnsPolynomial<ModularInt>> EncodeBfv(
      absl::Span<const Integer> messages,
      absl::Span<const PrimeModulus<ModularInt>* const> moduli,
      bool is_scaled = true) const;

  // Returns the messages encoded in the most significant bits of coefficients
  // of `noisy_plaintext`.
  absl::StatusOr<std::vector<Integer>> DecodeBfv(
      RnsPolynomial<ModularInt> noisy_plaintext,
      absl::Span<const PrimeModulus<ModularInt>* const> moduli) const;

  // Accessors
  int LogN() const { return context_->LogN(); }

  // The plaintext modulus.
  Integer PlaintextModulus() const { return context_->PlaintextModulus(); }

  // The underlying RNS context.
  const RnsContext<ModularInt>* Context() const { return context_; }

 private:
  explicit CoefficientEncoder(const RnsContext<ModularInt>* context)
      : context_(context) {}

  // Doesn't own the RnsContext, which must live longer than this encoder.
  const RnsContext<ModularInt>* context_;
};

}  // namespace rlwe

#endif  // RLWE_RNS_COEFFICIENT_ENCODER_H_
