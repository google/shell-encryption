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

#include <vector>

#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "absl/strings/str_cat.h"
#include "shell_encryption/integral_types.h"
#include "shell_encryption/modulus_conversion.h"
#include "shell_encryption/montgomery.h"
#include "shell_encryption/rns/error_correction.h"
#include "shell_encryption/rns/rns_context.h"
#include "shell_encryption/rns/rns_integer.h"
#include "shell_encryption/rns/rns_modulus.h"

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
  if (noisy_plaintext.IsNttForm()) {
    RLWE_RETURN_IF_ERROR(
        noisy_plaintext.ConvertToCoeffForm(moduli.subspan(0, 1)));
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

template class CoefficientEncoder<MontgomeryInt<Uint16>>;
template class CoefficientEncoder<MontgomeryInt<Uint32>>;
template class CoefficientEncoder<MontgomeryInt<Uint64>>;
template class CoefficientEncoder<MontgomeryInt<absl::uint128>>;
#ifdef ABSL_HAVE_INTRINSIC_INT128
template class CoefficientEncoder<MontgomeryInt<unsigned __int128>>;
#endif

}  // namespace rlwe
