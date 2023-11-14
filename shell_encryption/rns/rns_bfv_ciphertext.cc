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

#include "shell_encryption/rns/rns_bfv_ciphertext.h"

#include <utility>
#include <vector>

#include "absl/numeric/int128.h"
#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "shell_encryption/integral_types.h"
#include "shell_encryption/montgomery.h"
#include "shell_encryption/rns/rns_integer.h"
#include "shell_encryption/rns/rns_modulus.h"
#include "shell_encryption/rns/rns_polynomial.h"
#include "shell_encryption/status_macros.h"

namespace rlwe {

template <typename ModularInt>
absl::Status RnsBfvCiphertext<ModularInt>::AbsorbInPlace(
    const RnsPolynomial<ModularInt>& plaintext) {
  if (!plaintext.IsNttForm()) {
    return absl::InvalidArgumentError("`plaintext` must be in NTT form.");
  }

  typename ModularInt::Int t = context_->PlaintextModulus();
  std::vector<const PrimeModulus<ModularInt>*> aux_moduli =
      context_->AuxPrimeModuli();
  int level = this->Level();
  RLWE_ASSIGN_OR_RETURN(std::vector<ModularInt> q_hat_inv_mod_qs,
                        context_->MainPrimeModulusCrtFactors(level));
  RLWE_ASSIGN_OR_RETURN(std::vector<RnsInt<ModularInt>> q_hat_mod_ps,
                        context_->MainPrimeModulusComplementResidues(level));
  RLWE_ASSIGN_OR_RETURN(auto q_mod_ps,
                        context_->MainLeveledModulusAuxResidues(level));

  // Q/t * plaintext (mod P).
  RnsPolynomial<ModularInt> plaintext_coeff = plaintext;
  RLWE_RETURN_IF_ERROR(plaintext_coeff.ConvertToCoeffForm(this->Moduli()));
  RLWE_ASSIGN_OR_RETURN(RnsPolynomial<ModularInt> plaintext_aux,
                        plaintext_coeff.SwitchRnsBasis(
                            this->Moduli(), aux_moduli, q_hat_inv_mod_qs,
                            q_hat_mod_ps, q_mod_ps.RnsRep()));
  // Convert `plaintext_aux` to NTT form for multiplication.
  RLWE_RETURN_IF_ERROR(plaintext_aux.ConvertToNttForm(aux_moduli));

  for (RnsPolynomial<ModularInt>& ci : this->components()) {
    if (ci.IsNttForm()) {
      RLWE_RETURN_IF_ERROR(ci.ConvertToCoeffForm(this->Moduli()));
    }
    // P/t * ci (mod P).
    RLWE_ASSIGN_OR_RETURN(RnsPolynomial<ModularInt> ci_scaled_aux,
                          ci.ScaleAndSwitchRnsBasis(
                              this->Moduli(), aux_moduli, q_hat_inv_mod_qs,
                              context_->MainPrimeModulusInverseAuxResidues(),
                              context_->AuxModulusResidues()));
    // P/t * ci (mod Q).
    RLWE_ASSIGN_OR_RETURN(
        RnsPolynomial<ModularInt> ci_scaled_main,
        ci_scaled_aux.SwitchRnsBasis(
            aux_moduli, this->Moduli(), context_->AuxPrimeModulusCrtFactors(),
            context_->AuxPrimeModulusComplementResidues(),
            context_->AuxModulusResidues()));

    // PQ/t^2 * ci * plaintext (mod Q).
    RLWE_RETURN_IF_ERROR(ci_scaled_main.ConvertToNttForm(this->Moduli()));
    RLWE_RETURN_IF_ERROR(ci_scaled_main.MulInPlace(plaintext, this->Moduli()));
    // PQ/t^2 * ci * plaintext (mod P).
    RLWE_RETURN_IF_ERROR(ci_scaled_aux.ConvertToNttForm(aux_moduli));
    RLWE_RETURN_IF_ERROR(ci_scaled_aux.MulInPlace(plaintext_aux, aux_moduli));

    // Q/t * ci * plaintext (mod Q).
    RLWE_RETURN_IF_ERROR(ci_scaled_main.ScaleAndReduceRnsBasisInPlace(
        std::move(ci_scaled_aux), t, this->Moduli(), aux_moduli,
        context_->AuxPrimeModulusCrtFactors(),
        context_->AuxPrimeModulusComplementResidues(),
        context_->AuxModulusResidues(), context_->AuxModulusInverseResidues()));
    ci = ci_scaled_main;
  }

  this->error() = this->Error() * this->ErrorParams()->B_plaintext();
  return absl::OkStatus();
}

template <typename ModularInt>
absl::StatusOr<RnsBfvCiphertext<ModularInt>> RnsBfvCiphertext<ModularInt>::Mul(
    const RnsBfvCiphertext& that) const {
  if (this->Level() != that.Level()) {
    return absl::InvalidArgumentError("`that` has a mismatched level.");
  }
  if (this->PowerOfS() != that.PowerOfS()) {
    return absl::InvalidArgumentError(
        "Ciphertexts must be encrypted with the same key power.");
  }

  // We follow the improved multiplication algorithm from
  // ``Revisiting Homomorphic Encryption Schemes for Finite Fields''
  // by Kim et. al, https://eprint.iacr.org/2021/204.
  // We compute the tensor of two vectors of polynomials mod Q, where we take
  // the following steps to multiply two polynomials c(X) and c'(X) mod Q:
  // 1) Scale c (mod Q) to c_aux = P/t * c (mod P);
  // 2) Convert RNS basis to get c_main = P/t * c (mod Q);
  // 3) Convert RNS basis to get c_aux' = Q/t * c' (mod P);
  // 4) Multiply polynomials from the previous steps wrt Q and P, and we get
  //    product_main = QP/t^2 * c * c' (mod Q) and
  //    product_aux  = QP/t^2 * c * c' (mod P);
  // 5) Scales down the polynomial (product_main, product_aux) mod (Q*P) to
  //    get Q/t * c * c' (mod Q).

  typename ModularInt::Int t = context_->PlaintextModulus();
  std::vector<const PrimeModulus<ModularInt>*> aux_moduli =
      context_->AuxPrimeModuli();
  int level = this->Level();
  RLWE_ASSIGN_OR_RETURN(std::vector<ModularInt> q_hat_inv_mod_qs,
                        context_->MainPrimeModulusCrtFactors(level));
  RLWE_ASSIGN_OR_RETURN(std::vector<RnsInt<ModularInt>> q_hat_mod_ps,
                        context_->MainPrimeModulusComplementResidues(level));
  RLWE_ASSIGN_OR_RETURN(auto q_mod_ps,
                        context_->MainLeveledModulusAuxResidues(level));

  // Q/t * that.ci (mod P) for all i.
  std::vector<RnsPolynomial<ModularInt>> that_components_aux;
  that_components_aux.reserve(that.components().size());
  for (RnsPolynomial<ModularInt> ci : that.components()) {
    RLWE_RETURN_IF_ERROR(ci.ConvertToCoeffForm(this->Moduli()));
    RLWE_ASSIGN_OR_RETURN(
        RnsPolynomial<ModularInt> ci_aux,
        ci.SwitchRnsBasis(this->Moduli(), aux_moduli, q_hat_inv_mod_qs,
                          q_hat_mod_ps, q_mod_ps.RnsRep()));
    // Convert `ci_aux` to NTT form for multiplication.
    RLWE_RETURN_IF_ERROR(ci_aux.ConvertToNttForm(aux_moduli));
    that_components_aux.push_back(std::move(ci_aux));
  }

  // Create a vector of zero RNS polynomials
  RLWE_ASSIGN_OR_RETURN(auto zero,
                        RnsPolynomial<ModularInt>::CreateZero(
                            this->LogN(), this->Moduli(), /*is_ntt=*/true));
  RLWE_ASSIGN_OR_RETURN(auto zero_aux,
                        RnsPolynomial<ModularInt>::CreateZero(
                            this->LogN(), aux_moduli, /*is_ntt=*/true));
  int tensor_size = this->components().size() + that.components().size() - 1;
  std::vector<RnsPolynomial<ModularInt>> tensor(tensor_size, zero);
  std::vector<RnsPolynomial<ModularInt>> tensor_aux(tensor_size, zero_aux);

  // Compute the tensor product mod (Q * P).
  for (int i = 0; i < this->components().size(); ++i) {
    RnsPolynomial<ModularInt> ci = this->components()[i];
    RLWE_RETURN_IF_ERROR(ci.ConvertToCoeffForm(this->Moduli()));

    // P/t * ci (mod P).
    RLWE_ASSIGN_OR_RETURN(RnsPolynomial<ModularInt> ci_scaled_aux,
                          ci.ScaleAndSwitchRnsBasis(
                              this->Moduli(), aux_moduli, q_hat_inv_mod_qs,
                              context_->MainPrimeModulusInverseAuxResidues(),
                              context_->AuxModulusResidues()));
    // P/t * ci (mod Q).
    RLWE_ASSIGN_OR_RETURN(
        RnsPolynomial<ModularInt> ci_scaled_main,
        ci_scaled_aux.SwitchRnsBasis(
            aux_moduli, this->Moduli(), context_->AuxPrimeModulusCrtFactors(),
            context_->AuxPrimeModulusComplementResidues(),
            context_->AuxModulusResidues()));
    RLWE_RETURN_IF_ERROR(ci_scaled_main.ConvertToNttForm(this->Moduli()));
    RLWE_RETURN_IF_ERROR(ci_scaled_aux.ConvertToNttForm(aux_moduli));

    for (int j = 0; j < that.components().size(); ++j) {
      // PQ/t^2 * ci * that.cj (mod Q).
      RLWE_RETURN_IF_ERROR(tensor[i + j].FusedMulAddInPlace(
          ci_scaled_main, that.components()[j], this->Moduli()));
      // PQ/t^2 * ci * that.cj (mod P).
      RLWE_RETURN_IF_ERROR(tensor_aux[i + j].FusedMulAddInPlace(
          ci_scaled_aux, that_components_aux[j], aux_moduli));
    }
  }

  // Scale down tensor (mod Q*P) to Q/t * tensor (mod Q).
  for (int i = 0; i < tensor_size; ++i) {
    // Q/t * tensor[i] (mod Q).
    RLWE_RETURN_IF_ERROR(tensor[i].ScaleAndReduceRnsBasisInPlace(
        std::move(tensor_aux[i]), t, this->Moduli(), aux_moduli,
        context_->AuxPrimeModulusCrtFactors(),
        context_->AuxPrimeModulusComplementResidues(),
        context_->AuxModulusResidues(), context_->AuxModulusInverseResidues()));
  }

  double error_product = this->Error() * that.Error();
  return RnsBfvCiphertext<ModularInt>(std::move(tensor), this->moduli(),
                                      this->PowerOfS(), error_product,
                                      this->ErrorParams(), this->context_);
}

template class RnsBfvCiphertext<MontgomeryInt<Uint16>>;
template class RnsBfvCiphertext<MontgomeryInt<Uint32>>;
template class RnsBfvCiphertext<MontgomeryInt<Uint64>>;
template class RnsBfvCiphertext<MontgomeryInt<absl::uint128>>;
#ifdef ABSL_HAVE_INTRINSIC_INT128
template class RnsBfvCiphertext<MontgomeryInt<unsigned __int128>>;
#endif

}  // namespace rlwe
