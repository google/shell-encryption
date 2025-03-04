/*
 * Copyright 2025 Google LLC.
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

#include "shell_encryption/rns/approximate_encoder.h"

#include <algorithm>
#include <cmath>
#include <complex>
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
#include "shell_encryption/ntt_parameters.h"
#include "shell_encryption/rns/crt_interpolation.h"
#include "shell_encryption/rns/rns_context.h"
#include "shell_encryption/rns/rns_integer.h"
#include "shell_encryption/rns/rns_modulus.h"
#include "shell_encryption/rns/rns_polynomial.h"
#include "shell_encryption/status_macros.h"

namespace rlwe {

namespace {

template <typename T>
static void BitReversal(const std::vector<unsigned int>& bitrevs,
                        std::vector<T>& item_to_reverse) {
  for (int i = 0; i < item_to_reverse.size(); ++i) {
    // Only swap in one direction - don't accidentally swap twice.
    unsigned int r = bitrevs[i];
    if (i < r) {
      std::swap(item_to_reverse[i], item_to_reverse[r]);
    }
  }
}

}  // namespace

template <typename ModularInt>
absl::StatusOr<ApproximateEncoder<ModularInt>>
ApproximateEncoder<ModularInt>::Create(const RnsContext<ModularInt>* context,
                                       double scaling_factor) {
  if (context == nullptr) {
    return absl::InvalidArgumentError("`context` must not be null.");
  }
  if (scaling_factor < 1) {
    return absl::InvalidArgumentError("`scaling_factor` must be at least 1.");
  }

  int log_n = context->LogN();
  int n = 1 << log_n;

  // Generate the primitive 2N'th roots exp(j * pi / N * I) and their inverse.
  std::vector<std::complex<double>> psis_bitrev(n, {0, 0});
  std::vector<std::complex<double>> psis_bitrev_inv(n, {0, 0});
  double theta = M_PI / n;
  for (int j = 0; j < n; ++j) {
    psis_bitrev[j].real(cos(theta * j));
    psis_bitrev[j].imag(sin(theta * j));
    psis_bitrev_inv[j].real(cos(theta * j));
    psis_bitrev_inv[j].imag(-sin(theta * j));
  }

  // Cooley-Tukey uses the primitive roots in the bit reversed order, and for
  // Gentleman-Sande we need the inverses of the reversed primitive roots.
  std::vector<unsigned int> psis_bitrev_indices = internal::BitrevArray(log_n);
  BitReversal(psis_bitrev_indices, psis_bitrev);
  std::reverse(psis_bitrev_inv.begin() + 1, psis_bitrev_inv.end());
  BitReversal(psis_bitrev_indices, psis_bitrev_inv);
  std::reverse(psis_bitrev_inv.begin(), psis_bitrev_inv.end());

  // Bit reversal indices for the slots.
  std::vector<unsigned int> slot_bitrev_indices =
      internal::BitrevArray(log_n - 1);

  return ApproximateEncoder<ModularInt>(
      context, scaling_factor, std::move(psis_bitrev),
      std::move(psis_bitrev_inv), std::move(slot_bitrev_indices));
}

template <typename ModularInt>
absl::StatusOr<std::vector<std::complex<double>>>
ApproximateEncoder<ModularInt>::InverseTransform(
    absl::Span<const std::complex<double>> slots) const {
  // Move the input slot values to coefficients for Gentleman-Sande: the slots
  // correspond to primitive roots of power [5**j % (2N)]_j, which are reordered
  // to [4j+1]_j to run Gentleman-Sande.
  int n = 1 << context_->LogN();
  int num_slots = n >> 1;
  std::vector<std::complex<double>> coeffs(num_slots, {0, 0});
  for (int j = 0, power = 1; j < slots.size();
       ++j, power = (power * 5) % (2 * n)) {
    int k = (power - 1) / 4;
    coeffs[k] = slots[j];
  }
  BitReversal(bitrevs_, coeffs);

  RLWE_RETURN_IF_ERROR(IterativeHalfGentlemanSande(coeffs, psis_bitrev_inv_));

  // Normalize the transformation.
  for (auto& coeff : coeffs) {
    coeff /= num_slots;
  }
  return coeffs;
}

template <typename ModularInt>
absl::StatusOr<std::vector<std::complex<double>>>
ApproximateEncoder<ModularInt>::ForwardTransform(
    std::vector<std::complex<double>> coeffs) const {
  int n = 1 << context_->LogN();
  int num_slots = n >> 1;

  RLWE_RETURN_IF_ERROR(IterativeHalfCooleyTukey(coeffs, psis_bitrev_));

  BitReversal(bitrevs_, coeffs);

  // Cooley-Tukey evaluates the input polynomial on primitive roots of power
  // [4j+1]_j, and we now move them to their slots which correspond to power
  // [5**j % 2N]_j.
  std::vector<std::complex<double>> slots;
  slots.reserve(num_slots);
  for (int j = 0, power = 1; j < num_slots;
       ++j, power = (power * 5) % (2 * n)) {
    int k = (power - 1) / 4;
    slots.push_back(coeffs[k]);
  }
  return slots;
}

template <typename ModularInt>
absl::StatusOr<RnsPolynomial<ModularInt>>
ApproximateEncoder<ModularInt>::EncodeCkks(
    absl::Span<const std::complex<double>> values,
    absl::Span<const PrimeModulus<ModularInt>* const> moduli) const {
  int n = 1 << context_->LogN();
  int num_slots = n >> 1;
  if (values.size() > num_slots) {
    return absl::InvalidArgumentError(absl::StrCat(
        "`values` cannot have more than ", num_slots, " elements."));
  }
  if (moduli.empty()) {
    return absl::InvalidArgumentError("`moduli` cannot be empty.");
  }

  RLWE_ASSIGN_OR_RETURN(std::vector<std::complex<double>> coeff_values,
                        InverseTransform(values));

  // Scale DFT^-1(coeff_values) and encode it to balanced representation.
  int encoding_threshold_bits = sizeof(BigInteger) * 8 - 1;
  double encoding_threshold = std::ldexp(1, encoding_threshold_bits);
  std::vector<BigInteger> coeffs(n);
  for (int i = 0; i < num_slots; ++i) {
    double re = std::round(coeff_values[i].real() * scaling_factor_);
    double im = std::round(coeff_values[i].imag() * scaling_factor_);
    if (abs(re) > encoding_threshold || abs(im) > encoding_threshold) {
      return absl::InvalidArgumentError("not enough precision.");
    }
    auto re_abs = static_cast<BigInteger>(abs(re));
    auto im_abs = static_cast<BigInteger>(abs(im));
    coeffs[i] = re < 0 ? -re_abs : re_abs;
    coeffs[i + num_slots] = im < 0 ? -im_abs : im_abs;
  }

  // Finally convert the mod-q0 representation to mod-Q representation.
  BigInteger encoding_modulus_half = BigInteger(1) << encoding_threshold_bits;
  std::vector<std::vector<ModularInt>> coeff_vectors(moduli.size());
  for (int i = 0; i < moduli.size(); ++i) {
    auto mod_params_qi = moduli[i]->ModParams();
    Integer qi = mod_params_qi->modulus;
    coeff_vectors[i].reserve(n);
    for (int j = 0; j < n; ++j) {
      if (coeffs[j] > encoding_modulus_half) {
        Integer coeff_abs_qi = static_cast<Integer>((-coeffs[j]) % qi);
        RLWE_ASSIGN_OR_RETURN(
            ModularInt coeff_mod_qi,
            ModularInt::ImportInt(qi - coeff_abs_qi, mod_params_qi));
        coeff_vectors[i].push_back(std::move(coeff_mod_qi));
      } else {
        Integer coeff_abs_qi = static_cast<Integer>(coeffs[j] % qi);
        RLWE_ASSIGN_OR_RETURN(
            ModularInt coeff_mod_qi,
            ModularInt::ImportInt(coeff_abs_qi, mod_params_qi));
        coeff_vectors[i].push_back(std::move(coeff_mod_qi));
      }
    }
  }

  RLWE_ASSIGN_OR_RETURN(
      auto poly, RnsPolynomial<ModularInt>::Create(std::move(coeff_vectors),
                                                   /*is_ntt=*/false));
  RLWE_RETURN_IF_ERROR(poly.ConvertToNttForm(moduli));
  return poly;
}

template <typename ModularInt>
absl::StatusOr<std::vector<std::complex<double>>>
ApproximateEncoder<ModularInt>::DecodeCkks(
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
  int num_decoding_moduli = 1;  // At least one prime modulus.
  int decoding_modulus_bits = moduli[0]->ModParams()->log_modulus;
  BigInteger decoding_modulus = static_cast<BigInteger>(moduli[0]->Modulus());
  for (int i = 1; i < num_moduli; ++i) {
    // CRT interpolation of {a_j}_0^i mod {q_j}_0^i can be carried out over
    // BigInteger if sum(a_j * (q_0 * .. * q_i), j = 0..i) fits in BigInteger.
    int qi_bits = moduli[i]->ModParams()->log_modulus;
    int crt_bits =
        2 * qi_bits + decoding_modulus_bits + std::ceil(std::log2(i + 1));
    if (crt_bits < sizeof(BigInteger) * 8) {
      num_decoding_moduli++;
      decoding_modulus_bits += qi_bits;
      decoding_modulus *= static_cast<BigInteger>(moduli[i]->Modulus());
    } else {
      break;
    }
  }
  absl::Span<const PrimeModulus<ModularInt>* const> decoding_moduli =
      moduli.subspan(0, num_decoding_moduli);

  // Modulus reduce until we can interpolate the plaintext polynomial and get
  // BigInteger coefficients.
  double scaling_final = scaling_factor_;
  for (int i = num_decoding_moduli; i < moduli.size(); ++i) {
    RLWE_RETURN_IF_ERROR(
        noisy_plaintext.MulInPlace(moduli[i]->Modulus(), moduli));
  }
  absl::Span<const RnsInt<ModularInt>> q_inv_mod_qs =
      context_->MainPrimeModulusInverseResidues();
  while (noisy_plaintext.NumModuli() > num_decoding_moduli) {
    int level = noisy_plaintext.NumModuli() - 1;
    RLWE_RETURN_IF_ERROR(noisy_plaintext.ModReduceMsb(
        q_inv_mod_qs[level].Prefix(level), moduli.subspan(0, level + 1)));
  }

  // Interpolate the RNS polynomial.
  std::vector<BigInteger> decoding_modulus_hats =
      RnsModulusComplements<ModularInt, BigInteger>(decoding_moduli);
  RLWE_ASSIGN_OR_RETURN(
      std::vector<ModularInt> decoding_modulus_hat_invs,
      this->context_->MainPrimeModulusCrtFactors(num_decoding_moduli - 1));
  RLWE_ASSIGN_OR_RETURN(std::vector<BigInteger> noisy_coeffs,
                        (CrtInterpolation<ModularInt, BigInteger>(
                            noisy_plaintext.Coeffs(), decoding_moduli,
                            decoding_modulus_hats, decoding_modulus_hat_invs)));

  // Rescale the integer coefficients and compute forward DFT.
  BigInteger decoding_modulus_half = decoding_modulus >> 1;
  int log_n = context_->LogN();
  int num_slots = 1 << (log_n - 1);
  std::vector<std::complex<double>> coeffs(num_slots, {0, 0});
  for (int i = 0; i < num_slots; ++i) {
    BigInteger coeff_re = noisy_coeffs[i];
    BigInteger coeff_im = noisy_coeffs[i + num_slots];
    if (coeff_re > decoding_modulus_half) {
      coeffs[i].real(-static_cast<double>(decoding_modulus - coeff_re));
    } else {
      coeffs[i].real(static_cast<double>(coeff_re));
    }
    if (coeff_im > decoding_modulus_half) {
      coeffs[i].imag(-static_cast<double>(decoding_modulus - coeff_im));
    } else {
      coeffs[i].imag(static_cast<double>(coeff_im));
    }
    coeffs[i] /= scaling_final;
  }
  return ForwardTransform(std::move(coeffs));
}

template class ApproximateEncoder<MontgomeryInt<Uint64>>;
template class ApproximateEncoder<MontgomeryInt<absl::uint128>>;
#ifdef ABSL_HAVE_INTRINSIC_INT128
template class ApproximateEncoder<MontgomeryInt<unsigned __int128>>;
#endif

}  // namespace rlwe
