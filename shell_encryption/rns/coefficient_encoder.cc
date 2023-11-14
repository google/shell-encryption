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

#include "shell_encryption/rns/coefficient_encoder.h"

#include <utility>
#include <vector>

#include "absl/numeric/int128.h"
#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "absl/strings/str_cat.h"
#include "absl/types/span.h"
#include "shell_encryption/integral_types.h"
#include "shell_encryption/modulus_conversion.h"
#include "shell_encryption/montgomery.h"
#include "shell_encryption/rns/error_correction.h"
#include "shell_encryption/rns/rns_context.h"
#include "shell_encryption/rns/rns_integer.h"
#include "shell_encryption/rns/rns_modulus.h"
#include "shell_encryption/rns/rns_polynomial.h"
#include "shell_encryption/status_macros.h"

namespace rlwe {

template <typename ModularInt>
absl::StatusOr<CoefficientEncoder<ModularInt>>
CoefficientEncoder<ModularInt>::Create(const RnsContext<ModularInt>* context) {
  if (context == nullptr) {
    return absl::InvalidArgumentError("`context` must not be null.");
  }
  return CoefficientEncoder<ModularInt>(context);
}

template <typename ModularInt>
absl::StatusOr<RnsPolynomial<ModularInt>>
CoefficientEncoder<ModularInt>::EncodeBgv(
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

  std::vector<std::vector<ModularInt>> coeff_vectors(moduli.size());
  for (int i = 0; i < moduli.size(); ++i) {
    const ModularIntParams* mod_params_qi = moduli[i]->ModParams();
    // Convert message to balanced representation mod qi.
    RLWE_ASSIGN_OR_RETURN(coeff_vectors[i],
                          ImportBalancedModularInt<ModularInt>(
                              messages, PlaintextModulus(), *mod_params_qi));
    // Pad the remaining coefficients with 0.
    coeff_vectors[i].resize(num_coeffs, ModularInt::ImportZero(mod_params_qi));
  }
  return RnsPolynomial<ModularInt>::Create(coeff_vectors, /*is_ntt=*/false);
}

template <typename ModularInt>
absl::StatusOr<std::vector<typename ModularInt::Int>>
CoefficientEncoder<ModularInt>::DecodeBgv(
    RnsPolynomial<ModularInt> noisy_plaintext,
    absl::Span<const PrimeModulus<ModularInt>* const> moduli) const {
  int num_moduli = noisy_plaintext.NumModuli();
  if (moduli.size() != num_moduli) {
    return absl::InvalidArgumentError(
        absl::StrCat("`moduli` must contain ", num_moduli, " RNS moduli."));
  }

  if (noisy_plaintext.IsNttForm()) {
    RLWE_RETURN_IF_ERROR(noisy_plaintext.ConvertToCoeffForm(moduli));
  }

  // First we reduce the modulus of the noisy plaintext to a single prime
  // modulus q_0.
  Integer t = PlaintextModulus();
  absl::Span<const RnsInt<ModularInt>> q_inv_mod_qs =
      context_->MainPrimeModulusInverseResidues();
  while (noisy_plaintext.NumModuli() > 1) {
    int level = noisy_plaintext.NumModuli() - 1;
    RLWE_RETURN_IF_ERROR(noisy_plaintext.ModReduceLsb(
        t, q_inv_mod_qs[level].Prefix(level), moduli.subspan(0, level + 1)));
  }

  // Next, extract the coefficients (mod q_0) and remove errors.
  int num_coeffs = noisy_plaintext.NumCoeffs();
  const ModularIntParams* mod_params_q0 = moduli[0]->ModParams();
  std::vector<Integer> noisy_coeffs(num_coeffs, Integer{0});
  for (int i = 0; i < num_coeffs; ++i) {
    noisy_coeffs[i] = noisy_plaintext.Coeffs()[0][i].ExportInt(mod_params_q0);
  }
  return RemoveErrorOnMsb(std::move(noisy_coeffs), mod_params_q0->modulus, t);
}

// This is an implementation of the improved BFV encoding algorithm from
// ``Revisiting Homomorphic Encryption Schemes for Finite Fields'' by Kim et.al,
// https://eprint.iacr.org/2021/204.
template <typename ModularInt>
absl::StatusOr<RnsPolynomial<ModularInt>>
CoefficientEncoder<ModularInt>::EncodeBfv(
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

  const ModularIntParams* mod_params_t =
      context_->PlaintextModulusParams().ModParams();
  if (mod_params_t == nullptr) {
    return absl::InvalidArgumentError(
        "RnsContext does not contain a valid plaintext modulus parameters.");
  }

  // First convert to Montgomery integers mod t.
  std::vector<ModularInt> scaled_messages_mod_t;
  scaled_messages_mod_t.reserve(messages.size());
  for (auto message : messages) {
    RLWE_ASSIGN_OR_RETURN(ModularInt message_mod_t,
                          ModularInt::ImportInt(message, mod_params_t));
    scaled_messages_mod_t.push_back(std::move(message_mod_t));
  }
  // Pad the remaining coefficients with 0.
  scaled_messages_mod_t.resize(num_coeffs,
                               ModularInt::ImportZero(mod_params_t));

  // If the encoded messages are to be encrypted, then we need to scale up the
  // polynomial coefficients in order to tolerate noises. In BFV we scale up
  // messages to ([-Q * messages] mod t) / t (mod Q).
  if (is_scaled) {
    RLWE_ASSIGN_OR_RETURN(
        ModularInt q_mod_t,
        ModularInt::ImportInt(context_->MainModulusPlaintextResidue(),
                              mod_params_t));
    // [-Q * messages] mod t.
    RLWE_RETURN_IF_ERROR(ModularInt::BatchMulInPlace(
        &scaled_messages_mod_t, q_mod_t.Negate(mod_params_t), mod_params_t));
  }

  // Convert the coefficients (mod t) to a polynomial (mod Q).
  RLWE_ASSIGN_OR_RETURN(RnsPolynomial<ModularInt> plaintext,
                        RnsPolynomial<ModularInt>::ConvertFromPolynomialCoeffs(
                            scaled_messages_mod_t, mod_params_t, moduli));

  if (is_scaled) {
    // ([-Q * messages] mod t) / t (mod Q).
    RLWE_RETURN_IF_ERROR(plaintext.MulInPlace(
        context_->PlaintextModulusInverseMainResidues(), moduli));
  }
  return plaintext;
}

template <typename ModularInt>
absl::StatusOr<std::vector<typename ModularInt::Int>>
CoefficientEncoder<ModularInt>::DecodeBfv(
    RnsPolynomial<ModularInt> noisy_plaintext,
    absl::Span<const PrimeModulus<ModularInt>* const> moduli) const {
  int num_moduli = noisy_plaintext.NumModuli();
  if (moduli.size() != num_moduli) {
    return absl::InvalidArgumentError(
        absl::StrCat("`moduli` must contain ", num_moduli, " RNS moduli."));
  }

  const PrimeModulus<ModularInt>& modulus_t =
      context_->PlaintextModulusParams();
  if (modulus_t.ModParams() == nullptr) {
    return absl::InvalidArgumentError(
        "RnsContext does not contain valid plaintext modulus parameters.");
  }

  if (noisy_plaintext.IsNttForm()) {
    RLWE_RETURN_IF_ERROR(noisy_plaintext.ConvertToCoeffForm(moduli));
  }

  // Compute round(t / Q * a(X)) mod t.
  int modulus_level = moduli.size() - 1;
  RLWE_ASSIGN_OR_RETURN(std::vector<ModularInt> q_hat_inv_mod_qs,
                        context_->MainPrimeModulusCrtFactors(modulus_level));
  std::vector<ModularInt> moduli_invs_mod_t(
      context_->MainPrimeModulusInversePlaintextResidues().begin(),
      context_->MainPrimeModulusInversePlaintextResidues().end());
  RnsInt<ModularInt> q_inv_mod_t{std::move(moduli_invs_mod_t)};

  RLWE_ASSIGN_OR_RETURN(
      RnsPolynomial<ModularInt> plaintext,
      noisy_plaintext.ScaleAndSwitchRnsBasis(
          moduli, {&modulus_t}, q_hat_inv_mod_qs, {q_inv_mod_t},
          context_->PlaintextModulusMainResidues()));

  int num_coeffs = plaintext.NumCoeffs();
  std::vector<Integer> coeffs;
  coeffs.reserve(num_coeffs);
  const ModularIntParams* mod_params_t = modulus_t.ModParams();
  for (int i = 0; i < num_coeffs; ++i) {
    coeffs.push_back(plaintext.Coeffs()[0][i].ExportInt(mod_params_t));
  }
  return coeffs;
}

template class CoefficientEncoder<MontgomeryInt<Uint16>>;
template class CoefficientEncoder<MontgomeryInt<Uint32>>;
template class CoefficientEncoder<MontgomeryInt<Uint64>>;
template class CoefficientEncoder<MontgomeryInt<absl::uint128>>;
#ifdef ABSL_HAVE_INTRINSIC_INT128
template class CoefficientEncoder<MontgomeryInt<unsigned __int128>>;
#endif

}  // namespace rlwe
