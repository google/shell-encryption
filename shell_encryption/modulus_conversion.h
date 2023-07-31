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

#ifndef RLWE_MODULUS_CONVERSION_H_
#define RLWE_MODULUS_CONVERSION_H_

#include "absl/types/span.h"
#include "shell_encryption/montgomery.h"
#include "shell_encryption/ntt_parameters.h"
#include "shell_encryption/polynomial.h"
#include "shell_encryption/status_macros.h"

namespace rlwe {

// Given a vector `vec_q` of Montgomery ints modulo q, returns a vector of
// Montgomery ints modulo p representing the values in `vec_q`. We assume that
// `vec_q` represents integer values in the range [-t/2, t/2], for t < min(q,p).
template <typename ModularInt>
absl::StatusOr<std::vector<ModularInt>> ConvertModulusBalanced(
    const std::vector<ModularInt>& vec_q,
    const typename ModularInt::Params& mod_params_q,
    const typename ModularInt::Params& mod_params_p) {
  typename ModularInt::Int modulus_q = mod_params_q.modulus;
  typename ModularInt::Int modulus_p = mod_params_p.modulus;
  if (modulus_q == modulus_p) {
    return vec_q;  // trivial case, no conversion is needed.
  }

  typename ModularInt::Int modulus_q_by_2 = modulus_q >> 1;
  typename ModularInt::Int diff =
      modulus_q > modulus_p ? (modulus_q - modulus_p) : (modulus_p - modulus_q);
  RLWE_ASSIGN_OR_RETURN(ModularInt diff_p,
                        ModularInt::ImportInt(diff, &mod_params_p));

  std::vector<ModularInt> vec_p;
  vec_p.reserve(vec_q.size());
  if (modulus_p > modulus_q) {  // new modulus > old modulus
    for (int i = 0; i < vec_q.size(); ++i) {
      typename ModularInt::Int x = vec_q[i].ExportInt(&mod_params_q);
      RLWE_ASSIGN_OR_RETURN(ModularInt x_p,
                            ModularInt::ImportInt(x, &mod_params_p));
      if (x > modulus_q_by_2) {
        // When x > q/2, it represents a negative value x - q, so the mod-p
        // representation is (x - q) + p = x + (p - q) = x + diff (mod p).
        x_p.AddInPlace(diff_p, &mod_params_p);
      } else {
        // For x <= q/2, it represents a positive value x, and by our assumption
        // that x < min(q,p), its mod-p representation is also x, so we do
        // nothing. Since this function is only applied on ciphertext that are
        // pseudorandom, timing information won't leak any secret.
      }
      vec_p.push_back(std::move(x_p));
    }
  } else {  // new modulus < old modulus
    for (int i = 0; i < vec_q.size(); ++i) {
      typename ModularInt::Int x = vec_q[i].ExportInt(&mod_params_q);
      RLWE_ASSIGN_OR_RETURN(ModularInt x_p,
                            ModularInt::ImportInt(x, &mod_params_p));
      if (x > modulus_q_by_2) {
        // Similar to the above when x > q/2, x represents a negative value
        // x - q, but since p < q, its mod-p representation is now
        // (x - q) + p = x - (q - p) = x - diff (mod p).
        x_p.SubInPlace(diff_p, &mod_params_p);
      } else {
        // Do nothing as x represents a positive value in [0, p).
      }
      vec_p.push_back(std::move(x_p));
    }
  }
  return vec_p;
}

// Given a vector `vs` of unsigned integers representing values modulo t,
// returns a vector of Montgomery ints modulo p representing the same values as
// in `vs`, where t is given in `modulus_t`. We assume that modulo t values are
// in the range [-t/2, t/2].
// Note: This is implements the same algorithm as `ConvertModulusBalanced`
// except that `vs` are now unsigned integers.
template <typename ModularInt>
absl::StatusOr<std::vector<ModularInt>> ImportBalancedModularInt(
    absl::Span<const typename ModularInt::Int> vs,
    typename ModularInt::Int modulus_t,
    const typename ModularInt::Params& mod_params_p) {
  typename ModularInt::Int modulus_t_by_2 = modulus_t >> 1;
  std::vector<ModularInt> vs_p;
  vs_p.reserve(vs.size());
  for (auto v : vs) {
    typename ModularInt::Int vv = v % modulus_t;
    if (vv > modulus_t_by_2) {
      RLWE_ASSIGN_OR_RETURN(
          ModularInt v_mod_p,
          ModularInt::ImportInt(modulus_t - vv, &mod_params_p));
      v_mod_p.NegateInPlace(&mod_params_p);
      vs_p.push_back(std::move(v_mod_p));
    } else {
      RLWE_ASSIGN_OR_RETURN(ModularInt v_mod_p,
                            ModularInt::ImportInt(vv, &mod_params_p));
      vs_p.push_back(std::move(v_mod_p));
    }
  }
  return vs_p;
}

// Given a polynomial (in NTT form) `a_q` whose coefficients are represented by
// values in Z_q, returns the same polynomial (in NTT form) whose coefficients
// are represented by values in Z_p, for a different modulus p.
// Same as above, we assume that the integer values of coefficients of `a_q` are
// in the range [-t/2, t/2], for t < min(q,p).
template <typename ModularInt>
absl::StatusOr<rlwe::Polynomial<ModularInt>>
ConvertModulusBalancedOnNttPolynomial(
    const Polynomial<ModularInt>& a_q,
    const typename ModularInt::Params& mod_params_q,
    const NttParameters<ModularInt>& ntt_params_q,
    const typename ModularInt::Params& mod_params_p,
    const NttParameters<ModularInt>& ntt_params_p) {
  std::vector<ModularInt> coeffs_q =
      a_q.InverseNtt(&ntt_params_q, &mod_params_q);
  RLWE_ASSIGN_OR_RETURN(
      std::vector<ModularInt> coeffs_p,
      ConvertModulusBalanced<ModularInt>(coeffs_q, mod_params_q, mod_params_p));
  return Polynomial<ModularInt>::ConvertToNtt(std::move(coeffs_p),
                                              &ntt_params_p, &mod_params_p);
}

// Given a vector `vs_q` representing values in [0,q), returns values reduced
// modulo p.
template <typename ModularInt>
absl::StatusOr<std::vector<ModularInt>> ConvertModulus(
    const std::vector<ModularInt>& vs_q,
    const typename ModularInt::Params& mod_params_q,
    const typename ModularInt::Params& mod_params_p) {
  if (mod_params_q.modulus == mod_params_p.modulus) {
    return vs_q;  // trivial case, no conversion is needed.
  }
  std::vector<ModularInt> vs_p;
  vs_p.reserve(vs_q.size());
  for (auto const& v_mod_q : vs_q) {
    typename ModularInt::Int v = v_mod_q.ExportInt(&mod_params_q);
    RLWE_ASSIGN_OR_RETURN(ModularInt v_mod_p,
                          ModularInt::ImportInt(v, &mod_params_p));
    vs_p.push_back(std::move(v_mod_p));
  }
  return vs_p;
}

}  // namespace rlwe

#endif  // RLWE_MODULUS_CONVERSION_H_
