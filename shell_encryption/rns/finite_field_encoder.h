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

#ifndef RLWE_RNS_FINITE_FIELD_ENCODER_H_
#define RLWE_RNS_FINITE_FIELD_ENCODER_H_

#include <memory>
#include <utility>
#include <vector>

#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "absl/types/span.h"
#include "shell_encryption/ntt_parameters.h"
#include "shell_encryption/rns/coefficient_encoder.h"
#include "shell_encryption/rns/crt_interpolation.h"
#include "shell_encryption/rns/error_correction.h"
#include "shell_encryption/rns/rns_context.h"
#include "shell_encryption/rns/rns_modulus.h"
#include "shell_encryption/rns/rns_polynomial.h"
#include "shell_encryption/status_macros.h"

namespace rlwe {

// This class specifies how to encode messages in F_t^N as a polynomial
// in Z[X]/(Q, X^N+1), where t is the plaintext modulus and Q is the ciphertext
// modulus. The plaintext modulus t must be a prime such that t == 1 (mod 2*N).
// This class implements the so-called "slot encoding" for BGV and BFV.
template <typename ModularInt>
class FiniteFieldEncoder {
 public:
  using Integer = typename ModularInt::Int;
  using ModularIntParams = typename ModularInt::Params;
  using NttParams = NttParameters<ModularInt>;

  // Returns a FiniteFieldEncoder that can encode messages in a finite field
  // F_t, where t is the  plaintext modulus defined in `context`.
  static absl::StatusOr<FiniteFieldEncoder<ModularInt>> Create(
      const RnsContext<ModularInt>* context);

  // Returns a polynomial (in coefficient form) whose slots encode
  // `messages` as in the BGV scheme. `messages` are expected to be in the
  // range [0, t).
  absl::StatusOr<RnsPolynomial<ModularInt>> EncodeBgv(
      absl::Span<const Integer> messages,
      absl::Span<const PrimeModulus<ModularInt>* const> moduli) const;

  // Returns the messages that are encoded in the slots of `plaintext` as in the
  // BGV scheme. Messages are returned in the range [0, t).
  absl::StatusOr<std::vector<Integer>> DecodeBgv(
      RnsPolynomial<ModularInt> plaintext,
      absl::Span<const PrimeModulus<ModularInt>* const> moduli) const;

  // Alternative decoding that allows for additional parameter flexibility
  // at the expense of being slower.
  template <typename BigInteger>
  absl::StatusOr<std::vector<Integer>> DecodeBgvWithCrt(
      RnsPolynomial<ModularInt> plaintext,
      absl::Span<const PrimeModulus<ModularInt>* const> moduli,
      absl::Span<const BigInteger> modulus_hats,
      absl::Span<const ModularInt> modulus_hat_invs) const;

  // Returns the plaintext polynomial where `messages` are encoded in the slots
  // of `plaintext` as in the BFV scheme. When `is_scaled` is set, this encoding
  // method can be used in BFV scheme to encrypt `messages`. When `is_scaled` is
  // set to false, this encoding method produces a polynomial suitable for
  // ciphertext-plaintext multiplication in BFV.
  absl::StatusOr<RnsPolynomial<ModularInt>> EncodeBfv(
      absl::Span<const Integer> messages,
      absl::Span<const PrimeModulus<ModularInt>* const> moduli,
      bool is_scaled = true) const;

  // Returns the messages that are encoded in the slots of `plaintext` as in the
  // BGV scheme.
  absl::StatusOr<std::vector<Integer>> DecodeBfv(
      RnsPolynomial<ModularInt> plaintext,
      absl::Span<const PrimeModulus<ModularInt>* const> moduli) const;

  // Return a set of unsigned messages where all elements are in the range
  // [0, t) from a set of signed `messages` in range [-t/2, t/2).
  // All elements of `messages` are expected to be between [-t/2, t/2).
  template <typename SignedInteger>
  absl::StatusOr<std::vector<Integer>> WrapSigned(
      absl::Span<const SignedInteger> messages) const;

  // Return a set of signed messages where all elements are in the range
  // [-t/2, t/2) from a set of unsigned `messages` in range [0, t).
  // All elements of `messages` are expected to be between [0, t).
  template <typename SignedInteger>
  absl::StatusOr<std::vector<SignedInteger>> UnwrapToSigned(
      absl::Span<const Integer> messages) const;

  // Accessors
  int LogN() const { return coeff_encoder_->LogN(); }

  // The plaintext modulus.
  Integer PlaintextModulus() const { return context_->PlaintextModulus(); }

  // The underlying RNS context.
  const RnsContext<ModularInt>* Context() const { return context_; }

 private:
  explicit FiniteFieldEncoder(
      const RnsContext<ModularInt>* context,
      std::unique_ptr<const CoefficientEncoder<ModularInt>> coeff_encoder,
      std::vector<int> slot_indices)
      : context_(context),
        coeff_encoder_(std::move(coeff_encoder)),
        slot_indices_(std::move(slot_indices)) {}

  // Doesn't own the RnsContext, which must live longer than this encoder.
  const RnsContext<ModularInt>* context_;

  // The coefficient encoder used to remove error during decoding.
  std::unique_ptr<const CoefficientEncoder<ModularInt>> coeff_encoder_;

  // The permutation from messages indices to their slot indices.
  std::vector<int> slot_indices_;
};

// The standard `Decode` algorithm starts with a value
//
// (m + t*e) mod q_0, ..., (m + t*e) mod q_L,
//
// reduces this to modulus q_0 via modulus reduction, and proceeds. This assumes
// that the plaintext modulus t is smaller than the first level prime modulus
// q_0, but modulus reduction may introduce errors depending on the value of q_i
// mod t, and therefore a constraint on parameters so that q_1 mod t is small.
//
// This alternative decoder first CRT interpolates the RNS values above to
// recover (m + p*e) mod (\prod_i q_i), and then reduces this value mod t.
// E.g. it has greater parameter flexibility (at higher computational cost).
template <typename ModularInt>
template <typename BigInteger>
absl::StatusOr<std::vector<typename ModularInt::Int>>
FiniteFieldEncoder<ModularInt>::DecodeBgvWithCrt(
    RnsPolynomial<ModularInt> noisy_plaintext,
    absl::Span<const PrimeModulus<ModularInt>* const> moduli,
    absl::Span<const BigInteger> modulus_hats,
    absl::Span<const ModularInt> modulus_hat_invs) const {
  // Ensure the two auxiliary inputs are consistent with
  // the existing moduli.
  if (modulus_hats.size() != modulus_hat_invs.size() ||
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

  // Remove error and get the coefficients of the plaintext polynomial (mod t).
  const ModularIntParams* mod_params_t =
      context_->PlaintextModulusParams().ModParams();
  const NttParams* ntt_params_t =
      context_->PlaintextModulusParams().NttParams();
  BigInteger t = static_cast<BigInteger>(mod_params_t->modulus);
  RLWE_ASSIGN_OR_RETURN(
      std::vector<BigInteger> coeffs,
      RemoveErrorOnMsb<BigInteger>(std::move(noisy_coeffs), q, t));

  // Convert coefficients of the message polynomial to slot values.
  int num_coeffs = coeffs.size();
  std::vector<ModularInt> slots;
  slots.reserve(num_coeffs);
  for (auto coeff : coeffs) {
    RLWE_ASSIGN_OR_RETURN(
        ModularInt slot,
        ModularInt::ImportInt(static_cast<typename ModularInt::Int>(coeff),
                              mod_params_t));
    slots.push_back(std::move(slot));
  }
  RLWE_RETURN_IF_ERROR(
      ForwardNumberTheoreticTransform(slots, *ntt_params_t, *mod_params_t));

  // Lastly we move slots to their positions in the message vector.
  std::vector<typename ModularInt::Int> values(num_coeffs, 0);
  for (int i = 0; i < num_coeffs; ++i) {
    values[i] = slots[slot_indices_[i]].ExportInt(mod_params_t);
  }
  return values;
}

template <typename ModularInt>
template <typename SignedInteger>
absl::StatusOr<std::vector<typename ModularInt::Int>>
FiniteFieldEncoder<ModularInt>::WrapSigned(
    absl::Span<const SignedInteger> messages) const {
  const Integer t = context_->PlaintextModulus();

  std::vector<Integer> unsigned_messages(messages.size());
  for (int i = 0; i < messages.size(); ++i) {
    auto unsigned_message = messages[i];
    // this conversion assumes t is odd, enforced by RnsContext.
    unsigned_messages[i] = static_cast<Integer>(
        (unsigned_message < 0) ? unsigned_message + t : unsigned_message);
  }
  return unsigned_messages;
}

template <typename ModularInt>
template <typename SignedInteger>
absl::StatusOr<std::vector<SignedInteger>>
FiniteFieldEncoder<ModularInt>::UnwrapToSigned(
    absl::Span<const Integer> messages) const {
  const Integer t = context_->PlaintextModulus();
  const Integer t_half = t / 2;  // really is floor(t/2) because t is odd.

  std::vector<SignedInteger> signed_messages(messages.size());
  for (int i = 0; i < messages.size(); ++i) {
    auto signed_message = static_cast<SignedInteger>(messages[i]);
    // this conversion assumes t is odd, enforced by RnsContext.
    signed_messages[i] =
        (signed_message > t_half) ? signed_message - t : signed_message;
  }
  return signed_messages;
}

}  // namespace rlwe

#endif  // RLWE_RNS_FINITE_FIELD_ENCODER_H_
