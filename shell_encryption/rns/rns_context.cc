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

#include "shell_encryption/rns/rns_context.h"

#include "absl/memory/memory.h"

namespace rlwe {

template <typename ModularInt>
absl::StatusOr<RnsContext<ModularInt>> RnsContext<ModularInt>::Create(
    int log_n, absl::Span<const typename ModularInt::Int> qs,
    absl::Span<const typename ModularInt::Int> ps,
    typename ModularInt::Int plaintext_modulus) {
  if (log_n <= 0) {
    return absl::InvalidArgumentError("`log_n` must be positive.");
  }
  if (plaintext_modulus == 0) {
    return absl::InvalidArgumentError("`plaintext_modulus` must be positive.");
  }

  std::vector<std::unique_ptr<const PrimeModulus<ModularInt>>> modulus_qs;
  modulus_qs.reserve(qs.size());
  for (const auto& q : qs) {
    RLWE_ASSIGN_OR_RETURN(std::unique_ptr<const ModularIntParams> mod_params_q,
                          ModularInt::Params::Create(q));
    RLWE_ASSIGN_OR_RETURN(
        NttParameters<ModularInt> ntt_params_q,
        rlwe::InitializeNttParameters<ModularInt>(log_n, mod_params_q.get()));
    auto ntt_params_q_ptr = std::make_unique<const NttParameters<ModularInt>>(
        std::move(ntt_params_q));
    // Append a unique_ptr owning a PrimeModulus that includes both Montgomery
    // and NTT parameters about q.
    auto modulus_q = absl::WrapUnique(new PrimeModulus<ModularInt>{
        std::move(mod_params_q), std::move(ntt_params_q_ptr)});
    modulus_qs.push_back(std::move(modulus_q));
  }

  std::vector<std::unique_ptr<const PrimeModulus<ModularInt>>> modulus_ps;
  modulus_ps.reserve(ps.size());
  for (const auto& p : ps) {
    RLWE_ASSIGN_OR_RETURN(std::unique_ptr<const ModularIntParams> mod_params_p,
                          ModularInt::Params::Create(p));
    RLWE_ASSIGN_OR_RETURN(
        NttParameters<ModularInt> ntt_params_p,
        rlwe::InitializeNttParameters<ModularInt>(log_n, mod_params_p.get()));
    auto ntt_params_p_ptr = std::make_unique<const NttParameters<ModularInt>>(
        std::move(ntt_params_p));
    // Append a unique_ptr owning a PrimeModulus that includes both Montgomery
    // and NTT parameters about p.
    auto modulus_p = absl::WrapUnique(new PrimeModulus<ModularInt>{
        std::move(mod_params_p), std::move(ntt_params_p_ptr)});
    modulus_ps.push_back(std::move(modulus_p));
  }

  RnsContext<ModularInt> context(log_n, std::move(modulus_qs),
                                 std::move(modulus_ps), plaintext_modulus);
  RLWE_RETURN_IF_ERROR(context.GenerateCrtConstantsForMainModulus());
  RLWE_RETURN_IF_ERROR(context.GenerateCrtConstantsForAuxModulus());
  RLWE_RETURN_IF_ERROR(context.GeneratePlaintextModulusConstants());
  return context;
}

template <typename ModularInt>
absl::Status RnsContext<ModularInt>::GenerateCrtConstantsForMainModulus() {
  // Computes q_i and q_i^(-1) (mod q_j) for all main moduli q_i, q_j.
  prime_qs_.reserve(modulus_qs_.size());
  prime_q_invs_.reserve(modulus_qs_.size());
  for (int i = 0; i < modulus_qs_.size(); ++i) {
    std::vector<ModularInt> qi_mod_qs;      // hold (q_i mod q_j for all j)
    std::vector<ModularInt> qi_inv_mod_qs;  // hold (q_i^(-1) mod q_j for all j)
    qi_mod_qs.reserve(modulus_qs_.size());
    qi_inv_mod_qs.reserve(modulus_qs_.size());
    for (int j = 0; j < modulus_qs_.size(); ++j) {
      const ModularIntParams* params_qj = modulus_qs_[j]->ModParams();
      if (j == i) {
        qi_mod_qs.push_back(ModularInt::ImportZero(params_qj));
        // q_i does not have an inverse mod q_i itself, so we use 1 as a
        // placeholder.
        qi_inv_mod_qs.push_back(ModularInt::ImportOne(params_qj));
      } else {
        RLWE_ASSIGN_OR_RETURN(
            ModularInt qi_mod_qj,
            ModularInt::ImportInt(modulus_qs_[i]->Modulus(), params_qj));
        qi_inv_mod_qs.push_back(qi_mod_qj.MultiplicativeInverse(params_qj));
        qi_mod_qs.push_back(std::move(qi_mod_qj));
      }
    }
    prime_qs_.push_back({std::move(qi_mod_qs)});
    prime_q_invs_.push_back({std::move(qi_inv_mod_qs)});
  }

  return absl::OkStatus();
}

template <typename ModularInt>
absl::StatusOr<std::vector<ModularInt>>
RnsContext<ModularInt>::MainPrimeModulusComplements(int level) const {
  if (level < 0 || level >= modulus_qs_.size()) {
    return absl::InvalidArgumentError(absl::StrCat(
        "`level` must be non-negative and at most ", modulus_qs_.size() - 1));
  }

  // Computes the complement modulus values \hat{Q}_{l,i} = Q_l/q_i (mod q_i),
  // where Q_l = q_0 * .. * q_l.
  std::vector<ModularInt> prime_q_hats;
  prime_q_hats.reserve(level + 1);
  for (int i = 0; i <= level; ++i) {
    const ModularIntParams* params_qi = modulus_qs_[i]->ModParams();
    auto qi_hat = ModularInt::ImportOne(params_qi);
    for (int j = 0; j <= level; ++j) {
      if (j != i) {
        qi_hat.MulInPlace(prime_qs_[j].Component(i), params_qi);
      }
    }
    prime_q_hats.push_back(std::move(qi_hat));
  }
  return prime_q_hats;
}

template <typename ModularInt>
absl::StatusOr<std::vector<ModularInt>>
RnsContext<ModularInt>::MainPrimeModulusCrtFactors(int level) const {
  if (level < 0 || level >= modulus_qs_.size()) {
    return absl::InvalidArgumentError(absl::StrCat(
        "`level` must be non-negative and at most ", modulus_qs_.size() - 1));
  }

  // Computes the inverse of complement modulus values \hat{Q}_{l,i}^(-1) =
  // (Q_l/q_i)^(-1) (mod q_i), where Q_l = q_0 * .. * q_l.
  std::vector<ModularInt> prime_q_hat_invs;
  prime_q_hat_invs.reserve(level + 1);
  for (int i = 0; i <= level; ++i) {
    const ModularIntParams* params_qi = modulus_qs_[i]->ModParams();
    auto qi_hat_inv = ModularInt::ImportOne(params_qi);
    for (int j = 0; j <= level; ++j) {
      if (j != i) {
        qi_hat_inv.MulInPlace(prime_q_invs_[j].Component(i), params_qi);
      }
    }
    prime_q_hat_invs.push_back(std::move(qi_hat_inv));
  }
  return prime_q_hat_invs;
}

template <typename ModularInt>
absl::StatusOr<std::vector<RnsInt<ModularInt>>>
RnsContext<ModularInt>::MainPrimeModulusComplementResidues(int level) const {
  if (level < 0 || level >= modulus_qs_.size()) {
    return absl::InvalidArgumentError(absl::StrCat(
        "`level` must be non-negative and at most ", modulus_qs_.size() - 1));
  }

  // Compute Q_l / q_i (mod P) = [(Q_L / q_i) (mod p_j) for all j].
  std::vector<RnsInt<ModularInt>> prime_q_hat_mod_ps;
  prime_q_hat_mod_ps.reserve(level + 1);
  for (int i = 0; i <= level; ++i) {
    std::vector<ModularInt> qi_hat_mod_ps;
    qi_hat_mod_ps.reserve(modulus_ps_.size());
    for (int j = 0; j < modulus_ps_.size(); ++j) {
      const ModularIntParams* params_pj = modulus_ps_[j]->ModParams();
      ModularInt qi_hat_mod_pj = ModularInt::ImportOne(params_pj);
      for (int k = 0; k <= level; ++k) {
        if (k != i) {
          RLWE_ASSIGN_OR_RETURN(
              ModularInt qk_mod_pj,
              ModularInt::ImportInt(modulus_qs_[k]->Modulus(), params_pj));
          qi_hat_mod_pj.MulInPlace(qk_mod_pj, params_pj);
        }
      }
      qi_hat_mod_ps.push_back(std::move(qi_hat_mod_pj));
    }
    prime_q_hat_mod_ps.push_back({std::move(qi_hat_mod_ps)});
  }
  return prime_q_hat_mod_ps;
}

template <typename ModularInt>
absl::Status RnsContext<ModularInt>::GenerateCrtConstantsForAuxModulus() {
  // Computes \hat{p}_j^(-1) = (P / p_j)^(-1) (mod p_j).
  prime_p_hat_invs_.reserve(modulus_ps_.size());
  for (int j = 0; j < modulus_ps_.size(); ++j) {
    const ModularIntParams* params_pj = modulus_ps_[j]->ModParams();
    auto p_hat_inv_mod_pj = ModularInt::ImportOne(params_pj);
    for (int k = 0; k < modulus_ps_.size(); ++k) {
      if (k != j) {
        RLWE_ASSIGN_OR_RETURN(
            ModularInt pk_mod_pj,
            ModularInt::ImportInt(modulus_ps_[k]->Modulus(), params_pj));
        p_hat_inv_mod_pj.MulInPlace(pk_mod_pj.MultiplicativeInverse(params_pj),
                                    params_pj);
      }
    }
    prime_p_hat_invs_.push_back(std::move(p_hat_inv_mod_pj));
  }

  // Computes P (mod q_i) and P^(-1) (mod q_i).
  p_mod_qs_.reserve(modulus_qs_.size());
  p_inv_mod_qs_.reserve(modulus_qs_.size());
  for (int i = 0; i < modulus_qs_.size(); ++i) {
    const ModularIntParams* params_qi = modulus_qs_[i]->ModParams();
    auto p_mod_qi = ModularInt::ImportOne(params_qi);
    auto p_inv_mod_qi = ModularInt::ImportOne(params_qi);
    for (int j = 0; j < modulus_ps_.size(); ++j) {
      RLWE_ASSIGN_OR_RETURN(
          ModularInt pj_mod_qi,
          ModularInt::ImportInt(modulus_ps_[j]->Modulus(), params_qi));
      p_mod_qi.MulInPlace(pj_mod_qi, params_qi);
      p_inv_mod_qi.MulInPlace(pj_mod_qi.MultiplicativeInverse(params_qi),
                              params_qi);
    }
    p_mod_qs_.push_back(std::move(p_mod_qi));
    p_inv_mod_qs_.push_back(std::move(p_inv_mod_qi));
  }

  // Computes (P / p_j) (mod q_i).
  prime_p_hat_mod_qs_.reserve(modulus_ps_.size());
  for (int j = 0; j < modulus_ps_.size(); ++j) {
    std::vector<ModularInt> pj_hat_mod_qs;
    pj_hat_mod_qs.reserve(modulus_qs_.size());
    for (int i = 0; i < modulus_qs_.size(); ++i) {
      const ModularIntParams* params_qi = modulus_qs_[i]->ModParams();
      auto pj_hat_mod_qi = ModularInt::ImportOne(params_qi);
      for (int k = 0; k < modulus_ps_.size(); ++k) {
        if (k != j) {
          RLWE_ASSIGN_OR_RETURN(
              ModularInt pk_mod_qi,
              ModularInt::ImportInt(modulus_ps_[k]->Modulus(), params_qi));
          pj_hat_mod_qi.MulInPlace(pk_mod_qi, params_qi);
        }
      }
      pj_hat_mod_qs.push_back(std::move(pj_hat_mod_qi));
    }
    prime_p_hat_mod_qs_.push_back({std::move(pj_hat_mod_qs)});
  }
  return absl::OkStatus();
}

template <typename ModularInt>
absl::Status RnsContext<ModularInt>::GeneratePlaintextModulusConstants() {
  // Computes P (mod t) over the big integer type rather than Montgomery Int
  // to allow an even t.
  using BigInt = typename ModularInt::BigInt;
  BigInt p_rem{1};
  BigInt t_big = static_cast<BigInt>(plaintext_modulus_);
  for (auto const& modulus_p : modulus_ps_) {
    p_rem *= static_cast<BigInt>(modulus_p->Modulus());
    p_rem %= t_big;
  }
  p_mod_t_ = static_cast<typename ModularInt::Int>(p_rem);

  // Computes t^(-1) (mod q_i).
  t_inv_mod_qs_.reserve(modulus_qs_.size());
  for (int i = 0; i < modulus_qs_.size(); ++i) {
    const ModularIntParams* mod_params_qi = modulus_qs_[i]->ModParams();
    RLWE_ASSIGN_OR_RETURN(
        ModularInt t_mod_qi,
        ModularInt::ImportInt(plaintext_modulus_, mod_params_qi));
    t_inv_mod_qs_.push_back(t_mod_qi.MultiplicativeInverse(mod_params_qi));
  }

  // Computes t^(-1) (mod p_j).
  t_inv_mod_ps_.reserve(modulus_ps_.size());
  for (int j = 0; j < modulus_ps_.size(); ++j) {
    const ModularIntParams* mod_params_pj = modulus_ps_[j]->ModParams();
    RLWE_ASSIGN_OR_RETURN(
        ModularInt t_mod_pj,
        ModularInt::ImportInt(plaintext_modulus_, mod_params_pj));
    t_inv_mod_ps_.push_back(t_mod_pj.MultiplicativeInverse(mod_params_pj));
  }
  return absl::OkStatus();
}

template class RnsContext<MontgomeryInt<Uint16>>;
template class RnsContext<MontgomeryInt<Uint32>>;
template class RnsContext<MontgomeryInt<Uint64>>;
template class RnsContext<MontgomeryInt<absl::uint128>>;
#ifdef ABSL_HAVE_INTRINSIC_INT128
template class RnsContext<MontgomeryInt<unsigned __int128>>;
#endif

}  // namespace rlwe
