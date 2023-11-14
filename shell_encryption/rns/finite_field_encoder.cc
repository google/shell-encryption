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

#include "shell_encryption/rns/finite_field_encoder.h"

#include <memory>
#include <utility>
#include <vector>

#include "absl/numeric/int128.h"
#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "absl/strings/str_cat.h"
#include "absl/types/span.h"
#include "shell_encryption/dft_transformations.h"
#include "shell_encryption/integral_types.h"
#include "shell_encryption/montgomery.h"
#include "shell_encryption/rns/coefficient_encoder.h"
#include "shell_encryption/rns/rns_context.h"
#include "shell_encryption/rns/rns_modulus.h"
#include "shell_encryption/rns/rns_polynomial.h"
#include "shell_encryption/status_macros.h"

namespace rlwe {

template <typename ModularInt>
absl::StatusOr<FiniteFieldEncoder<ModularInt>>
FiniteFieldEncoder<ModularInt>::Create(const RnsContext<ModularInt>* context) {
  if (context == nullptr) {
    return absl::InvalidArgumentError("`context` must not be null.");
  }
  // Montgomery and NTT parameters for the plaintext algebra should have been
  // initialized in `context`.
  const ModularIntParams* mod_params_t =
      context->PlaintextModulusParams().ModParams();
  const NttParams* ntt_params_t = context->PlaintextModulusParams().NttParams();
  if (mod_params_t == nullptr || ntt_params_t == nullptr) {
    return absl::InvalidArgumentError(
        "`context` must be created with finite field encoding support.");
  }

  // Coefficient encoder is used to remove error and extract the plaintext
  // polynomial.
  RLWE_ASSIGN_OR_RETURN(CoefficientEncoder<ModularInt> coeff_encoder,
                        CoefficientEncoder<ModularInt>::Create(context));

  // Compute the permutation that moves the encoding input values to their
  // corresponding slots. Since the encoding essentially computes inverse NTT,
  // the slot index of each input value is the bitreversed index of the output
  // of NTT.
  int n = 1 << context->LogN();  // number of slots.
  std::vector<int> slot_indices(n, 0);
  int curr_power = 1;
  for (int i = 0; i < n / 2; ++i) {
    slot_indices[i] = ntt_params_t->bitrevs[(curr_power - 1) / 2];
    int co_power = curr_power * (2 * n - 1) % (2 * n);
    slot_indices[i + n / 2] = ntt_params_t->bitrevs[(co_power - 1) / 2];
    curr_power = curr_power * 5 % (2 * n);
  }

  return FiniteFieldEncoder<ModularInt>(
      context,
      std::make_unique<const CoefficientEncoder<ModularInt>>(
          std::move(coeff_encoder)),
      std::move(slot_indices));
}

template <typename ModularInt>
absl::StatusOr<RnsPolynomial<ModularInt>>
FiniteFieldEncoder<ModularInt>::EncodeBgv(
    absl::Span<const Integer> messages,
    absl::Span<const PrimeModulus<ModularInt>* const> moduli) const {
  int num_coeffs = 1 << LogN();
  if (messages.size() > num_coeffs) {
    return absl::InvalidArgumentError(absl::StrCat(
        "`messages` can contain at most ", num_coeffs, " elements."));
  }
  if (moduli.empty()) {
    return absl::InvalidArgumentError("`moduli` cannot be empty.");
  }

  // First we move messages (mod t) to their slots in Z[X]/(t, X^N+1).
  const ModularIntParams* mod_params_t =
      context_->PlaintextModulusParams().ModParams();
  const NttParams* ntt_params_t =
      context_->PlaintextModulusParams().NttParams();
  std::vector<ModularInt> slots(num_coeffs,
                                ModularInt::ImportZero(mod_params_t));
  for (int i = 0; i < messages.size(); ++i) {
    RLWE_ASSIGN_OR_RETURN(slots[slot_indices_[i]],
                          ModularInt::ImportInt(messages[i], mod_params_t));
  }
  RLWE_RETURN_IF_ERROR(
      InverseNumberTheoreticTransform(slots, *ntt_params_t, *mod_params_t));

  // The plaintext polynomial is converted from (mod t) to (mod Q).
  return RnsPolynomial<ModularInt>::ConvertBalancedFromPolynomialCoeffs(
      slots, mod_params_t, moduli);
}

template <typename ModularInt>
absl::StatusOr<std::vector<typename ModularInt::Int>>
FiniteFieldEncoder<ModularInt>::DecodeBgv(
    RnsPolynomial<ModularInt> noisy_plaintext,
    absl::Span<const PrimeModulus<ModularInt>* const> moduli) const {
  // Remove error and get the coefficients of the plaintext polynomial (mod t).
  RLWE_ASSIGN_OR_RETURN(
      std::vector<Integer> coeffs,
      coeff_encoder_->DecodeBgv(std::move(noisy_plaintext), moduli));

  // Convert coefficients of the message polynomial to slot values.
  const ModularIntParams* mod_params_t =
      context_->PlaintextModulusParams().ModParams();
  const NttParams* ntt_params_t =
      context_->PlaintextModulusParams().NttParams();
  int num_coeffs = coeffs.size();
  std::vector<ModularInt> slots;
  slots.reserve(num_coeffs);
  for (auto coeff : coeffs) {
    RLWE_ASSIGN_OR_RETURN(ModularInt slot,
                          ModularInt::ImportInt(coeff, mod_params_t));
    slots.push_back(std::move(slot));
  }
  RLWE_RETURN_IF_ERROR(
      ForwardNumberTheoreticTransform(slots, *ntt_params_t, *mod_params_t));

  // Lastly we move slots to their positions in the message vector.
  std::vector<Integer> messages(num_coeffs, 0);
  for (int i = 0; i < num_coeffs; ++i) {
    messages[i] = slots[slot_indices_[i]].ExportInt(mod_params_t);
  }
  return messages;
}

template <typename ModularInt>
absl::StatusOr<RnsPolynomial<ModularInt>>
FiniteFieldEncoder<ModularInt>::EncodeBfv(
    absl::Span<const Integer> messages,
    absl::Span<const PrimeModulus<ModularInt>* const> moduli,
    bool is_scaled) const {
  int num_coeffs = 1 << LogN();
  if (messages.size() > num_coeffs) {
    return absl::InvalidArgumentError(absl::StrCat(
        "`messages` can contain at most ", num_coeffs, " elements."));
  }
  if (moduli.empty()) {
    return absl::InvalidArgumentError("`moduli` cannot be empty.");
  }

  // First we move messages (mod t) to their slots in Z[X]/(t, X^N+1).
  const ModularIntParams* mod_params_t =
      context_->PlaintextModulusParams().ModParams();
  const NttParams* ntt_params_t =
      context_->PlaintextModulusParams().NttParams();
  std::vector<ModularInt> slots(num_coeffs,
                                ModularInt::ImportZero(mod_params_t));
  for (int i = 0; i < messages.size(); ++i) {
    RLWE_ASSIGN_OR_RETURN(slots[slot_indices_[i]],
                          ModularInt::ImportInt(messages[i], mod_params_t));
  }
  RLWE_RETURN_IF_ERROR(
      InverseNumberTheoreticTransform(slots, *ntt_params_t, *mod_params_t));

  std::vector<Integer> slot_values;
  slot_values.reserve(slots.size());
  for (auto const& slot : slots) {
    slot_values.push_back(slot.ExportInt(mod_params_t));
  }
  return coeff_encoder_->EncodeBfv(slot_values, moduli, is_scaled);
}

template <typename ModularInt>
absl::StatusOr<std::vector<typename ModularInt::Int>>
FiniteFieldEncoder<ModularInt>::DecodeBfv(
    RnsPolynomial<ModularInt> noisy_plaintext,
    absl::Span<const PrimeModulus<ModularInt>* const> moduli) const {
  // Remove error and get the coefficients of the plaintext polynomial (mod t).
  RLWE_ASSIGN_OR_RETURN(
      std::vector<Integer> coeffs,
      coeff_encoder_->DecodeBfv(std::move(noisy_plaintext), moduli));

  // Convert coefficients of the message polynomial to slot values.
  const ModularIntParams* mod_params_t =
      context_->PlaintextModulusParams().ModParams();
  const NttParams* ntt_params_t =
      context_->PlaintextModulusParams().NttParams();
  int num_coeffs = coeffs.size();
  std::vector<ModularInt> slots;
  slots.reserve(num_coeffs);
  for (auto coeff : coeffs) {
    RLWE_ASSIGN_OR_RETURN(ModularInt slot,
                          ModularInt::ImportInt(coeff, mod_params_t));
    slots.push_back(std::move(slot));
  }
  RLWE_RETURN_IF_ERROR(
      ForwardNumberTheoreticTransform(slots, *ntt_params_t, *mod_params_t));

  // Lastly we move slots to their positions in the message vector.
  std::vector<Integer> messages(num_coeffs, 0);
  for (int i = 0; i < num_coeffs; ++i) {
    messages[i] = slots[slot_indices_[i]].ExportInt(mod_params_t);
  }
  return messages;
}

template class FiniteFieldEncoder<MontgomeryInt<Uint16>>;
template class FiniteFieldEncoder<MontgomeryInt<Uint32>>;
template class FiniteFieldEncoder<MontgomeryInt<Uint64>>;
template class FiniteFieldEncoder<MontgomeryInt<absl::uint128>>;
#ifdef ABSL_HAVE_INTRINSIC_INT128
template class FiniteFieldEncoder<MontgomeryInt<unsigned __int128>>;
#endif

}  // namespace rlwe
