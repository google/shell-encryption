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

#include "shell_encryption/rns/rns_gadget.h"

#include <cstddef>
#include <vector>

#include "absl/numeric/int128.h"
#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "absl/strings/str_cat.h"
#include "absl/types/span.h"
#include "shell_encryption/gadget.h"
#include "shell_encryption/integral_types.h"
#include "shell_encryption/modulus_conversion.h"
#include "shell_encryption/montgomery.h"
#include "shell_encryption/rns/rns_integer.h"
#include "shell_encryption/rns/rns_modulus.h"
#include "shell_encryption/rns/rns_polynomial.h"
#include "shell_encryption/status_macros.h"

namespace rlwe {

// Create a gadget with CRT representation given the bases and the modulus Q
// specified in the RNS context. Note that a gadget is the concatenation of
// gadgets in the residual ring Z_{qi}:
//   g_crt = [ qi_hat * qi_hat_inv * g_i for i in 1..L ] mod Q,
// where Q = q1 * ... * qL, qi_hat = Q / q_i, qi_hat_inv = qi_hat^(-1) mod qi.
template <typename ModularInt>
absl::StatusOr<RnsGadget<ModularInt>> RnsGadget<ModularInt>::Create(
    int log_n, std::vector<size_t> log_bs,
    absl::Span<const ModularInt> q_hats,      // {qi_hat mod qi}_i
    absl::Span<const ModularInt> q_hat_invs,  // {qi_hat_inv}_i
    absl::Span<const PrimeModulus<ModularInt>* const> moduli) {
  if (log_n <= 0) {
    return absl::InvalidArgumentError("`log_n` must be positive.");
  }
  if (moduli.empty()) {
    return absl::InvalidArgumentError("`moduli` must not be empty.");
  }
  int num_moduli = moduli.size();  // number of moduli in RNS chain.
  if (log_bs.size() != num_moduli) {
    return absl::InvalidArgumentError(
        absl::StrCat("`log_bs` must contain ", num_moduli, " elements."));
  }
  if (q_hats.size() != num_moduli) {
    return absl::InvalidArgumentError(
        absl::StrCat("`q_hats` must contain ", num_moduli, " elements."));
  }
  if (q_hat_invs.size() != num_moduli) {
    return absl::InvalidArgumentError(
        absl::StrCat("`q_hat_invs` must contain ", num_moduli, " elements."));
  }

  int k = 0;  // dimension of the gadget.
  std::vector<Integer> bs(num_moduli, 0);
  std::vector<size_t> ks(num_moduli, 0);  // dimension of sub-gadget wrt q_i.
  for (int i = 0; i < num_moduli; ++i) {
    bs[i] = 1 << log_bs[i];
    ks[i] = GadgetSize<ModularInt>(log_bs[i], moduli[i]->ModParams());
    k += ks[i];
  }

  std::vector<ModularInt> zeros;
  zeros.reserve(num_moduli);
  for (auto modulus : moduli) {
    zeros.push_back(ModularInt::ImportZero(modulus->ModParams()));
  }

  std::vector<RnsInt<ModularInt>> gs;
  gs.reserve(k);
  for (int i = 0; i < num_moduli; ++i) {
    auto mod_params_qi = moduli[i]->ModParams();
    RLWE_ASSIGN_OR_RETURN(auto bi, ModularInt::ImportInt(bs[i], mod_params_qi));
    // gs[i'th part] = qi_hat * qi_hat_inv * [1, bi, ..., bi^(ki-1)] (mod qi),
    // and it is [0, 0, 0, ...0] mod qj for j != i.
    ModularInt current_entry_qi =
        q_hats[i].Mul(q_hat_invs[i], mod_params_qi);  // qi_hat * qi_hat_inv
    for (int j = 0; j < ks[i]; ++j) {
      RnsInt<ModularInt> entry{zeros};
      entry.zs[i] = current_entry_qi;  // set the component mod qi.
      gs.push_back(std::move(entry));

      // Update to qi_hat * qi_hat_inv * bi^j for next j.
      current_entry_qi.MulInPlace(bi, mod_params_qi);
    }
  }
  return RnsGadget<ModularInt>(std::move(log_bs), std::move(ks), std::move(gs));
}

// Gadget-decompose a RNS polynomial a into xs, such that a equals to the inner
// product <xs, g> mod Q. The vector xs has the form
//   xs = [ g_1^(-1)(a), ..., g_l^(-1)(a) ],
// where g_i is the sub-gadget for the sub-modulus qi
template <typename ModularInt>
absl::StatusOr<std::vector<RnsPolynomial<ModularInt>>>
RnsGadget<ModularInt>::Decompose(
    const RnsPolynomial<ModularInt>& a,
    absl::Span<const PrimeModulus<ModularInt>* const> moduli) const {
  if (a.IsNttForm()) {
    return absl::InvalidArgumentError("`a` must be in coefficient form");
  }
  int num_moduli = moduli.size();  // number of moduli in RNS chain
  if (a.NumModuli() != num_moduli) {
    return absl::InvalidArgumentError(absl::StrCat(
        "`a` must be defined with respect to ", num_moduli, " RNS moduli."));
  }

  int k = Dimension();
  std::vector<RnsPolynomial<ModularInt>> xs;
  xs.reserve(k);
  for (int i = 0; i < num_moduli; ++i) {
    // The ith block of xs, g_i^(-1)(a), can be computed from ai = a mod qi,
    // and then convert to the other moduli.
    auto mod_params_qi = moduli[i]->ModParams();
    RLWE_ASSIGN_OR_RETURN(
        std::vector<std::vector<ModularInt>> x_coeffs_qi,
        BaseDecompose<ModularInt>(a.Coeffs()[i], mod_params_qi, log_bs_[i],
                                  dims_[i]));
    for (int j = 0; j < dims_[i]; ++j) {
      std::vector<std::vector<ModularInt>> x_coeff_vectors(num_moduli);
      for (int k = 0; k < num_moduli; ++k) {
        if (k != i) {
          // Convert from mod qi to mod qk in balanced form
          auto mod_params_qk = moduli[k]->ModParams();
          RLWE_ASSIGN_OR_RETURN(
              x_coeff_vectors[k],
              ConvertModulusBalanced<ModularInt>(x_coeffs_qi[j], *mod_params_qi,
                                                 *mod_params_qk));
        }
      }
      x_coeff_vectors[i] = std::move(x_coeffs_qi[j]);
      RLWE_ASSIGN_OR_RETURN(
          auto x, RnsPolynomial<ModularInt>::Create(std::move(x_coeff_vectors),
                                                    /*is_ntt=*/false));
      xs.push_back(std::move(x));
    }
  }
  return xs;
}

template class RnsGadget<MontgomeryInt<Uint16>>;
template class RnsGadget<MontgomeryInt<Uint32>>;
template class RnsGadget<MontgomeryInt<Uint64>>;
template class RnsGadget<MontgomeryInt<absl::uint128>>;
#ifdef ABSL_HAVE_INTRINSIC_INT128
template class RnsGadget<MontgomeryInt<unsigned __int128>>;
#endif

}  // namespace rlwe
