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

#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "absl/strings/str_cat.h"
#include "absl/types/span.h"
#include "shell_encryption/rns/crt_interpolation.h"
#include "shell_encryption/rns/error_correction.h"
#include "shell_encryption/rns/rns_context.h"
#include "shell_encryption/rns/rns_modulus.h"
#include "shell_encryption/rns/rns_polynomial.h"
#include "shell_encryption/status_macros.h"

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

  // Extended version of `EncodeBgv` above that allows the `messages` to be
  // a vector of big integers, where each coefficient of the returned polynomial
  // encodes a member of `messages`. All members of `messages` are assumed to be
  // smaller than `plaintext_modulus`.
  // Note that `plaintext_modulus` may be different from the plaintext modulus
  // defined in the RNS context.
  template <typename BigInteger>
  absl::StatusOr<RnsPolynomial<ModularInt>> EncodeBgv(
      absl::Span<const BigInteger> messages, BigInteger plaintext_modulus,
      absl::Span<const PrimeModulus<ModularInt>* const> moduli) const;

  // Returns the messages encoded in the least significant bits of coefficients
  // of `noisy_plaintext`.
  absl::StatusOr<std::vector<Integer>> DecodeBgv(
      RnsPolynomial<ModularInt> noisy_plaintext,
      absl::Span<const PrimeModulus<ModularInt>* const> moduli) const;

  // Alternative decoding that allows `plaintext_modulus` being larger than
  // an Integer, where the `plaintext` polynomial is defined wrt RNS `moduli`.
  // The parameter `modulus_hats` stores {Q/q_i}_i where {q_i}_i are the
  // RNS moduli, and `modulus_hat_invs` stores the CRT factors
  // {[(Q/q_i)^(-1) (mod q_i)]}_i.
  // Note that `plaintext_modulus` may be different from the plaintext modulus
  // defined in the RNS context.
  template <typename BigInteger>
  absl::StatusOr<std::vector<BigInteger>> DecodeBgv(
      RnsPolynomial<ModularInt> plaintext, BigInteger plaintext_modulus,
      absl::Span<const PrimeModulus<ModularInt>* const> moduli,
      absl::Span<const BigInteger> modulus_hats,
      absl::Span<const ModularInt> modulus_hat_invs) const;

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

template <typename ModularInt>
template <typename BigInteger>
absl::StatusOr<RnsPolynomial<ModularInt>>
CoefficientEncoder<ModularInt>::EncodeBgv(
    absl::Span<const BigInteger> messages, BigInteger plaintext_modulus,
    absl::Span<const PrimeModulus<ModularInt>* const> moduli) const {
  int num_coeffs = 1 << LogN();
  if (messages.size() > num_coeffs) {
    return absl::InvalidArgumentError(absl::StrCat(
        "`messages` can contain at most ", num_coeffs, " elements."));
  }
  if (moduli.empty()) {
    return absl::InvalidArgumentError("`moduli` cannot be empty.");
  }

  BigInteger plaintext_modulus_half = plaintext_modulus / 2;

  std::vector<std::vector<ModularInt>> coeff_vectors(moduli.size());
  for (int i = 0; i < moduli.size(); ++i) {
    const ModularIntParams* mod_params_qi = moduli[i]->ModParams();
    Integer qi = mod_params_qi->modulus;
    BigInteger qi_big = static_cast<BigInteger>(qi);

    // Convert message to balanced representation mod qi.
    coeff_vectors[i].reserve(num_coeffs);
    for (auto message : messages) {
      bool is_negative = message > plaintext_modulus_half;
      BigInteger message_qi_abs =
          (is_negative ? (plaintext_modulus - message) : message) % qi_big;
      RLWE_ASSIGN_OR_RETURN(
          auto coeff_mod_qi,
          ModularInt::ImportInt(is_negative
                                    ? qi - static_cast<Integer>(message_qi_abs)
                                    : static_cast<Integer>(message_qi_abs),
                                mod_params_qi));
      coeff_vectors[i].push_back(std::move(coeff_mod_qi));
    }
    // Pad the remaining coefficients with 0.
    coeff_vectors[i].resize(num_coeffs, ModularInt::ImportZero(mod_params_qi));
  }

  return RnsPolynomial<ModularInt>::Create(coeff_vectors, /*is_ntt=*/false);
}

template <typename ModularInt>
template <typename BigInteger>
absl::StatusOr<std::vector<BigInteger>>
CoefficientEncoder<ModularInt>::DecodeBgv(
    RnsPolynomial<ModularInt> noisy_plaintext, BigInteger plaintext_modulus,
    absl::Span<const PrimeModulus<ModularInt>* const> moduli,
    absl::Span<const BigInteger> modulus_hats,
    absl::Span<const ModularInt> modulus_hat_invs) const {
  // Ensure the two auxiliary arrays of constants are consistent with
  // the existing moduli.
  if (moduli.size() != noisy_plaintext.NumModuli() ||
      moduli.size() != modulus_hats.size() ||
      moduli.size() != modulus_hat_invs.size()) {
    return absl::InvalidArgumentError(
        "`moduli`, `modulus_hats` and `modulus_hat_invs` must have the same "
        "size.");
  }

  // CRT interpolation
  if (noisy_plaintext.IsNttForm()) {
    RLWE_RETURN_IF_ERROR(noisy_plaintext.ConvertToCoeffForm(moduli));
  }
  RLWE_ASSIGN_OR_RETURN(
      std::vector<BigInteger> noisy_coeffs,
      (CrtInterpolation<ModularInt, BigInteger>(
          noisy_plaintext.Coeffs(), moduli, absl::MakeSpan(modulus_hats),
          absl::MakeSpan(modulus_hat_invs))));
  // Compute the composite modulus Q
  BigInteger q{1};
  for (auto modulus : moduli) {
    q *= static_cast<BigInteger>(modulus->Modulus());
  }

  // Remove error and get coefficients of the plaintext polynomial (mod t).
  return RemoveErrorOnMsb<BigInteger>(std::move(noisy_coeffs), q,
                                      plaintext_modulus);
}

}  // namespace rlwe

#endif  // RLWE_RNS_COEFFICIENT_ENCODER_H_
