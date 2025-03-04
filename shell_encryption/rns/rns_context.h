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

#ifndef RLWE_RNS_RNS_CONTEXT_H_
#define RLWE_RNS_RNS_CONTEXT_H_

#include <memory>
#include <utility>
#include <vector>

#include "absl/types/span.h"
#include "shell_encryption/montgomery.h"
#include "shell_encryption/ntt_parameters.h"
#include "shell_encryption/rns/rns_integer.h"
#include "shell_encryption/rns/rns_modulus.h"

namespace rlwe {

// The RNS variant of the RLWE schemes utilize arithmetic in Z[X]/(Q,X^N+1) for
// a composite modulus Q = q0 * .. * q_{L-1} for distinct primes q_i. We call
// q_i's the main prime moduli, and their (partial) products the main modulus at
// a certain level. By Chinese Remainder Theorem, computation modulo Q can be
// decomposed into computation in smaller rings modulo q_i to use efficient
// modular arithmetic over q_i. We will also use an auxiliary modulus P = p_1 *
// ... * p_K for distinct primes p_j in certain special operations such as
// key-switching. Conversions between different representations and moduli are
// often needed, and to speed up computation we pre-compute relevant constants.
//
// This class holds the parameters of the prime moduli q_i's and p_j's, as well
// as pre-computed constants relevant to them.
// Note that we require q_i == p_j == 1 (mod 2*N) to enable NTT for polynomials
// modular q_i and p_j.
// Furthermore, we assume that the plaintext modulus t is smaller than all q_i.
template <typename ModularInt>
class RnsContext {
 public:
  using ModularIntParams = typename ModularInt::Params;

  // RnsContext is movable and move-assignable, but not copyable nor
  // copy-assignable.
  RnsContext(const RnsContext<ModularInt>&) = delete;
  RnsContext& operator=(const RnsContext<ModularInt>&) = delete;
  RnsContext(RnsContext<ModularInt>&&) = default;
  RnsContext& operator=(RnsContext<ModularInt>&&) = default;

  ~RnsContext() = default;

  // Creates a RnsContext with main prime moduli `qs` and auxiliary
  // prime moduli `ps`, and with a plaintext modulus.
  static absl::StatusOr<RnsContext> Create(
      int log_n, absl::Span<const typename ModularInt::Int> qs,
      absl::Span<const typename ModularInt::Int> ps,
      typename ModularInt::Int plaintext_modulus);

  // Creates a RnsContext suitable for using the finite field encoder with BGV
  // scheme, with main prime moduli `qs`, auxiliary prime moduli `ps`, and a
  // plaintext modulus.
  static absl::StatusOr<RnsContext> CreateForBgvFiniteFieldEncoding(
      int log_n, absl::Span<const typename ModularInt::Int> qs,
      absl::Span<const typename ModularInt::Int> ps,
      typename ModularInt::Int plaintext_modulus);

  // Creates a RnsContext suitable for instantiating the BFV scheme, with main
  // prime moduli `qs`, auxiliary prime moduli `ps`, and a plaintext modulus.
  static absl::StatusOr<RnsContext> CreateForBfv(
      int log_n, absl::Span<const typename ModularInt::Int> qs,
      absl::Span<const typename ModularInt::Int> ps,
      typename ModularInt::Int plaintext_modulus);

  // Creates a RnsContext suitable for using the finite field encoder with BFV
  // scheme, with main prime moduli `qs`, auxiliary prime moduli `ps`, and a
  // plaintext modulus.
  static absl::StatusOr<RnsContext> CreateForBfvFiniteFieldEncoding(
      int log_n, absl::Span<const typename ModularInt::Int> qs,
      absl::Span<const typename ModularInt::Int> ps,
      typename ModularInt::Int plaintext_modulus);

  // Creates a RnsContext suitable for instantiating the CKKS scheme, with main
  // prime moduli `qs`, auxiliary prime moduli `ps`.
  static absl::StatusOr<RnsContext> CreateForCkks(
      int log_n, absl::Span<const typename ModularInt::Int> qs,
      absl::Span<const typename ModularInt::Int> ps);

  // Returns the log of the dimension N.
  int LogN() const { return log_n_; }

  // Returns the number of prime moduli constituting the main modulus.
  int NumMainPrimeModuli() const { return modulus_qs_.size(); }

  // Returns the number of prime moduli constituting the auxiliary modulus.
  int NumAuxPrimeModuli() const { return modulus_ps_.size(); }

  // Returns a vector of raw pointers to main prime moduli.
  std::vector<const PrimeModulus<ModularInt>*> MainPrimeModuli() const {
    std::vector<const PrimeModulus<ModularInt>*> moduli;
    std::transform(modulus_qs_.begin(), modulus_qs_.end(),
                   std::back_inserter(moduli),
                   [](auto& ptr) { return ptr.get(); });
    return moduli;
  }

  // Returns a vector of raw pointers to auxiliary prime moduli.
  std::vector<const PrimeModulus<ModularInt>*> AuxPrimeModuli() const {
    std::vector<const PrimeModulus<ModularInt>*> moduli;
    std::transform(modulus_ps_.begin(), modulus_ps_.end(),
                   std::back_inserter(moduli),
                   [](auto& ptr) { return ptr.get(); });
    return moduli;
  }

  // Returns the plaintext modulus t.
  typename ModularInt::Int PlaintextModulus() const {
    return plaintext_modulus_;
  }

  // Returns the parameters (Montgomery integer and NTT) for the plaintext
  // modulus t.
  // Note that the Montgomery parameters is defined only when this context is
  // created for BFV scheme, and the NTT parameters is defined only when this
  // is created for using finite field encoding.
  const PrimeModulus<ModularInt>& PlaintextModulusParams() const {
    return modulus_t_;
  }

  // The i'th entry is q_i (mod Q) = (q_i (mod q_j) for all j).
  absl::Span<const RnsInt<ModularInt>> MainPrimeModulusResidues() const {
    return prime_qs_;
  }

  // The i'th entry is q_i^(-1) (mod Q) = (q_i^(-1) (mod q_j) for all j).
  // Note that q_i does not have an inverse modulo q_i itself, so we use 1 as
  // a placeholder.
  absl::Span<const RnsInt<ModularInt>> MainPrimeModulusInverseResidues() const {
    return prime_q_invs_;
  }

  // The i'th entry is Q^(-1) (mod p_i) = (q_j^(-1) (mod p_i) for all j).
  absl::Span<const RnsInt<ModularInt>> MainPrimeModulusInverseAuxResidues()
      const {
    return prime_q_inv_mod_ps_;
  }

  // The l'th entry is [Q_l (mod p_j) = \prod_{i=0}^l q_i (mod p_j) for all j].
  absl::StatusOr<RnsInt<ModularInt>> MainLeveledModulusAuxResidues(
      int level) const {
    if (level < 0 || level >= modulus_qs_.size()) {
      return absl::InvalidArgumentError(absl::StrCat(
          "`level` must be non-negative and at most ", modulus_qs_.size() - 1));
    }
    return leveled_q_mod_ps_[level];
  }

  // The i'th entry is \hat{Q}_{l,i} = Q_l / q_i (mod q_i), where Q_l = q_0 * ..
  // * q_l for l = `level`.
  absl::StatusOr<std::vector<ModularInt>> MainPrimeModulusComplements(
      int level) const;

  // The i'th entry is \hat{Q}_{l,i}^(-1) = (Q_l / q_i)^(-1) (mod q_i), where
  // Q_l = q_0 * .. * q_l for l = `level`. These values are useful to compute
  // CRT interpolation and RNS basis conversion.
  absl::StatusOr<std::vector<ModularInt>> MainPrimeModulusCrtFactors(
      int level) const;

  // The i'th entry is \hat{Q}_{l,i} (mod P) = ((Q_l/q_i) (mod p_j) for all j).
  absl::StatusOr<std::vector<RnsInt<ModularInt>>>
  MainPrimeModulusComplementResidues(int level) const;

  // The j'th entry is \hat{p}_j^(-1) = (P / p_j)^(-1) (mod p_j), which is
  // useful to compute CRT interpolation.
  absl::Span<const ModularInt> AuxPrimeModulusCrtFactors() const {
    return prime_p_hat_invs_;
  }

  // The j'th entry is P/p_j (mod Q) = ((P / p_j) (mod q_i) for all i).
  absl::Span<const RnsInt<ModularInt>> AuxPrimeModulusComplementResidues()
      const {
    return prime_p_hat_mod_qs_;
  }

  // The i'th entry is P (mod q_i).
  absl::Span<const ModularInt> AuxModulusResidues() const { return p_mod_qs_; }

  // The i'th entry is P^(-1) (mod q_i).
  absl::Span<const ModularInt> AuxModulusInverseResidues() const {
    return p_inv_mod_qs_;
  }

  // Returns the residue of Q (mod t), where t is the plaintext modulus.
  typename ModularInt::Int MainModulusPlaintextResidue() const {
    return q_mod_t_;
  }

  // Returns the residue of P (mod t), where t is the plaintext modulus.
  typename ModularInt::Int AuxModulusPlaintextResidue() const {
    return p_mod_t_;
  }

  // The i'th entry is t (mod q_i).
  absl::Span<const ModularInt> PlaintextModulusMainResidues() const {
    return t_mod_qs_;
  }

  // The i'th entry is t^(-1) (mod q_i).
  absl::Span<const ModularInt> PlaintextModulusInverseMainResidues() const {
    return t_inv_mod_qs_;
  }

  // The j'th entry is t^(-1) (mod p_j).
  absl::Span<const ModularInt> PlaintextModulusInverseAuxResidues() const {
    return t_inv_mod_ps_;
  }

  // The i'th entry is q_i (mod t).
  absl::Span<const ModularInt> MainPrimeModulusPlaintextResidues() const {
    return qs_mod_t_;
  }

  // The i'th entry is q_i^(-1) (mod t).
  absl::Span<const ModularInt> MainPrimeModulusInversePlaintextResidues()
      const {
    return q_invs_mod_t_;
  }

  // The l'th entry is Q_l^(-1) (mod t) = \prod_{i=0}^l q_i^(-1) (mod t).
  absl::Span<const ModularInt> MainLeveledModulusInversePlaintextResidue()
      const {
    return ql_invs_mod_t_;
  }

 private:
  explicit RnsContext(
      int log_n,
      std::vector<std::unique_ptr<const PrimeModulus<ModularInt>>> modulus_qs,
      std::vector<std::unique_ptr<const PrimeModulus<ModularInt>>> modulus_ps)
      : log_n_(log_n),
        modulus_qs_(std::move(modulus_qs)),
        modulus_ps_(std::move(modulus_ps)) {}

  // Creates a RnsContext for common operations among RLWE-based schemes, with
  // main prime moduli `qs`, auxiliary prime moduli `ps`.
  static absl::StatusOr<RnsContext> CreateCommon(
      int log_n, absl::Span<const typename ModularInt::Int> qs,
      absl::Span<const typename ModularInt::Int> ps);

  // Computes the CRT constants relevant to main and auxiliary prime moduli.
  absl::Status GenerateCrtConstantsForMainModulus();
  absl::Status GenerateCrtConstantsForAuxModulus();
  absl::Status GeneratePlaintextModulusConstants();
  absl::Status GeneratePlaintextModulusConstantsForFiniteFieldEncoding();
  absl::Status GeneratePlaintextModulusConstantsForBfv();

  int log_n_;

  // Montgomery and NTT parameters for the main prime moduli q_i.
  std::vector<std::unique_ptr<const PrimeModulus<ModularInt>>> modulus_qs_;

  // Montgomery and NTT parameters for the auxiliary prime moduli p_j.
  std::vector<std::unique_ptr<const PrimeModulus<ModularInt>>> modulus_ps_;

  // Plaintext modulus t.
  typename ModularInt::Int plaintext_modulus_;

  // Montgomery and (optional) NTT parameters for the plaintext modulus t.
  PrimeModulus<ModularInt> modulus_t_;

  // The entry [i] is q_i (mod Q) = [q_i (mod q_j) for all j].
  std::vector<RnsInt<ModularInt>> prime_qs_;

  // The entry [i] is q_i^(-1) (mod Q) = [q_i^(-1) (mod q_j) for all j].
  std::vector<RnsInt<ModularInt>> prime_q_invs_;

  // The entry [i] is Q^(-1) (mod p_i), stored as [q_j^-1 (mod p_i) for all j].
  std::vector<RnsInt<ModularInt>> prime_q_inv_mod_ps_;

  // The entry [l] is [Q_l (mod p_j) = \prod_{i=0}^l q_i (mod p_j) for all j].
  std::vector<RnsInt<ModularInt>> leveled_q_mod_ps_;

  // The entry [i] is \hat{p}_i^(-1) = (P / p_i)^(-1) (mod p_i).
  std::vector<ModularInt> prime_p_hat_invs_;

  // The entry [i] is P / p_i (mod Q) = [(P / p_i) (mod q_j) for all j].
  std::vector<RnsInt<ModularInt>> prime_p_hat_mod_qs_;

  // The entry [i] is P (mod q_i).
  std::vector<ModularInt> p_mod_qs_;

  // The entry [i] is P^(-1) (mod q_i).
  std::vector<ModularInt> p_inv_mod_qs_;

  // The residue of Q (mod t).
  typename ModularInt::Int q_mod_t_;

  // The residue of P (mod t).
  typename ModularInt::Int p_mod_t_;

  // The entry [i] is t (mod q_i).
  std::vector<ModularInt> t_mod_qs_;

  // The entry [i] is t^(-1) (mod q_i).
  std::vector<ModularInt> t_inv_mod_qs_;

  // The entry [i] is t^(-1) (mod p_i).
  std::vector<ModularInt> t_inv_mod_ps_;

  // The entry [i] is q_i (mod t).
  std::vector<ModularInt> qs_mod_t_;

  // The entry [i] is q_i^(-1) (mod t).
  std::vector<ModularInt> q_invs_mod_t_;

  // The entry [i] is Q_l^(-1) (mod t) = \prod_{j=0}^l q_j^(-1) (mod t).
  std::vector<ModularInt> ql_invs_mod_t_;
};

}  // namespace rlwe

#endif  // RLWE_RNS_RNS_CONTEXT_H_
