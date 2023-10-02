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

#ifndef RLWE_RNS_CRT_INTERPOLATION_H_
#define RLWE_RNS_CRT_INTERPOLATION_H_

#include <vector>

#include "absl/numeric/int128.h"
#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "shell_encryption/integral_types.h"
#include "shell_encryption/montgomery.h"
#include "shell_encryption/rns/rns_modulus.h"

namespace rlwe {

namespace {

// Converts an `Integer` value to a `BigInteger` value. We assume `BigInteger`
// is as large as `Integer`.
// This is needed in case there is no conversion operator from `Integer` to
// `BigInteger`.
template <typename Integer, typename BigInteger>
static BigInteger ConvertToBigInteger(Integer value) {
  return static_cast<BigInteger>(value);
}

}  // namespace

// Given a list of prime modulus {q_0,..,q_L} as specified in `moduli`, returns
// the list {Q/q_0, .., Q/q_L}, where Q = q_0 * ... * q_L.
template <typename ModularInt, typename BigInteger>
std::vector<BigInteger> RnsModulusComplements(
    absl::Span<const PrimeModulus<ModularInt>* const> moduli) {
  std::vector<BigInteger> modulus_hats(moduli.size(), BigInteger{1});
  for (int i = 0; i < moduli.size(); ++i) {
    for (int j = 0; j < moduli.size(); ++j) {
      if (j != i) {
        modulus_hats[i] *=
            ConvertToBigInteger<typename ModularInt::Int, BigInteger>(
                moduli[j]->Modulus());
      }
    }
  }
  return modulus_hats;
}

// Returns the coefficients mod Q of the RNS polynomial wrt prime moduli q_0,..,
// q_L, where Q = q_0 * .. * q_L. The RNS polynomial's coefficients are given in
// `coeff_vectors`.
// The CRT interpolation of (a_0,..,a_L), a_i \in Z_{q_i}, is
//    a = sum(a_i * [(Q/q_i)^(-1) (mod q_i)] * (Q/q_i), i=0,..,L) mod Q,
// where the constants [(Q/q_i)^(-1) (mod q_i)] are given in `modulus_hat_invs`,
// and the constants Q/q_i are given in `modulus_hats`.
template <typename ModularInt, typename BigInteger>
absl::StatusOr<std::vector<BigInteger>> CrtInterpolation(
    const std::vector<std::vector<ModularInt>>& coeff_vectors,
    absl::Span<const PrimeModulus<ModularInt>* const> moduli,
    absl::Span<const BigInteger> modulus_hats,
    absl::Span<const ModularInt> modulus_hat_invs) {
  using Integer = typename ModularInt::Int;
  if (coeff_vectors.empty()) {
    return absl::InvalidArgumentError("`coeff_vectors` must not be empty.");
  }
  int num_moduli = coeff_vectors.size();
  if (moduli.size() != num_moduli) {
    return absl::InvalidArgumentError(
        absl::StrCat("`moduli` must contain ", num_moduli, " RNS moduli."));
  }
  if (modulus_hats.size() != num_moduli) {
    return absl::InvalidArgumentError(
        absl::StrCat("`modulus_hats` must contain ", num_moduli, " elements."));
  }
  if (modulus_hat_invs.size() != num_moduli) {
    return absl::InvalidArgumentError(absl::StrCat(
        "`modulus_hat_invs` must contain ", num_moduli, " elements."));
  }

  // Compute the composite modulus Q = q_0 * ... * q_L.
  BigInteger q{1};
  for (auto prime_modulus : moduli) {
    q *= ConvertToBigInteger<Integer, BigInteger>(prime_modulus->Modulus());
  }

  // Interpolate mod-Q.
  int num_coeffs = coeff_vectors[0].size();
  std::vector<BigInteger> coeffs_q(num_coeffs, BigInteger{0});
  for (int i = 0; i < num_moduli; ++i) {
    auto mod_params_qi = moduli[i]->ModParams();
    for (int j = 0; j < num_coeffs; ++j) {
      auto factor = ConvertToBigInteger<Integer, BigInteger>(
          modulus_hat_invs[i].ExportInt(mod_params_qi));
      factor *= modulus_hats[i];
      factor %= q;
      factor *= ConvertToBigInteger<Integer, BigInteger>(
          coeff_vectors[i][j].ExportInt(mod_params_qi));
      factor %= q;
      coeffs_q[j] += factor;
      coeffs_q[j] %= q;
    }
  }
  return coeffs_q;
}

}  // namespace rlwe

#endif  // RLWE_RNS_CRT_INTERPOLATION_H_
