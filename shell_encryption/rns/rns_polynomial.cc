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

#include "shell_encryption/rns/rns_polynomial.h"

#include "shell_encryption/dft_transformations.h"
#include "shell_encryption/modulus_conversion.h"

namespace rlwe {

namespace {

// Returns true if n <= 0 or n is not a power of two.
inline bool IsNotPowerOfTwo(int n) {
  if (n <= 0) return true;
  return (n & (n - 1)) != 0;
}

}  // namespace

template <typename ModularInt>
absl::StatusOr<RnsPolynomial<ModularInt>> RnsPolynomial<ModularInt>::Create(
    std::vector<std::vector<ModularInt>> coeff_vectors, bool is_ntt) {
  if (coeff_vectors.empty()) {
    return absl::InvalidArgumentError("`coeff_vectors` cannot be empty.");
  }
  int num_coeffs = coeff_vectors[0].size();
  if (IsNotPowerOfTwo(num_coeffs)) {
    return absl::InvalidArgumentError(
        "`coeff_vectors` must contain vectors of length a power of two.");
  }
  for (auto const& coeff_vector : coeff_vectors) {
    if (coeff_vector.size() != num_coeffs) {
      return absl::InvalidArgumentError(
          "`coeff_vectors` must contain vectors of equal length.");
    }
  }
  return RnsPolynomial(log2(num_coeffs), std::move(coeff_vectors), is_ntt);
}

template <typename ModularInt>
absl::StatusOr<RnsPolynomial<ModularInt>> RnsPolynomial<ModularInt>::CreateZero(
    int log_n, absl::Span<const PrimeModulus<ModularInt>* const> moduli,
    bool is_ntt) {
  if (log_n <= 0) {
    return absl::InvalidArgumentError("`log_n` must be positive.");
  }
  if (moduli.empty()) {
    return absl::InvalidArgumentError("`moduli` must not be empty.");
  }
  std::vector<std::vector<ModularInt>> coeff_vectors;
  coeff_vectors.reserve(moduli.size());
  int num_coeffs = 1 << log_n;
  for (auto modulus : moduli) {
    // Create a zero polynomial mod the ith main modulus.
    std::vector<ModularInt> coeffs(
        num_coeffs, ModularInt::ImportZero(modulus->ModParams()));
    coeff_vectors.push_back(std::move(coeffs));
  }
  return RnsPolynomial<ModularInt>(log_n, std::move(coeff_vectors), is_ntt);
}

template <typename ModularInt>
absl::StatusOr<RnsPolynomial<ModularInt>> RnsPolynomial<ModularInt>::CreateOne(
    int log_n, absl::Span<const PrimeModulus<ModularInt>* const> moduli) {
  if (log_n <= 0) {
    return absl::InvalidArgumentError("`log_n` must be positive.");
  }
  if (moduli.empty()) {
    return absl::InvalidArgumentError("`moduli` must not be empty.");
  }
  std::vector<std::vector<ModularInt>> coeff_vectors;
  coeff_vectors.reserve(moduli.size());
  int num_coeffs = 1 << log_n;
  for (auto modulus : moduli) {
    // The NTT coefficients of the polynomial 1 are all 1 wrt q_i.
    std::vector<ModularInt> coeffs(num_coeffs,
                                   ModularInt::ImportOne(modulus->ModParams()));
    coeff_vectors.push_back(std::move(coeffs));
  }
  return RnsPolynomial<ModularInt>(log_n, std::move(coeff_vectors),
                                   /*is_ntt=*/true);
}

template <typename ModularInt>
absl::StatusOr<RnsPolynomial<ModularInt>>
RnsPolynomial<ModularInt>::ConvertFromPolynomialCoeffs(
    const std::vector<ModularInt>& coeffs_q,  // coefficients (mod q)
    const ModularIntParams* mod_params_q,
    absl::Span<const PrimeModulus<ModularInt>* const> moduli) {
  int num_coeffs = coeffs_q.size();
  if (IsNotPowerOfTwo(num_coeffs)) {
    return absl::InvalidArgumentError(
        "`coeffs_q` must have length a power of 2.");
  }
  if (mod_params_q == nullptr) {
    return absl::InvalidArgumentError("`mod_params_q` must not be null.");
  }
  if (moduli.empty()) {
    return absl::InvalidArgumentError("`moduli` must not be empty.");
  }
  std::vector<std::vector<ModularInt>> coeff_vectors;
  coeff_vectors.reserve(moduli.size());
  for (const auto* modulus : moduli) {
    RLWE_ASSIGN_OR_RETURN(std::vector<ModularInt> coeffs_qi,
                          ConvertModulus<ModularInt>(coeffs_q, *mod_params_q,
                                                     *modulus->ModParams()));
    RLWE_RETURN_IF_ERROR(ForwardNumberTheoreticTransform(
        coeffs_qi, *modulus->NttParams(), *modulus->ModParams()));
    coeff_vectors.push_back(std::move(coeffs_qi));
  }
  return RnsPolynomial<ModularInt>(log2(num_coeffs), std::move(coeff_vectors),
                                   /*is_ntt=*/true);
}

template <typename ModularInt>
absl::StatusOr<RnsPolynomial<ModularInt>>
RnsPolynomial<ModularInt>::ConvertBalancedFromPolynomialCoeffs(
    const std::vector<ModularInt>& coeffs_q,  // coefficients (mod q)
    const ModularIntParams* mod_params_q,
    absl::Span<const PrimeModulus<ModularInt>* const> moduli) {
  int num_coeffs = coeffs_q.size();
  if ((num_coeffs & (num_coeffs - 1)) != 0) {
    return absl::InvalidArgumentError(
        "`coeffs_q` must have length a power of 2.");
  }
  if (mod_params_q == nullptr) {
    return absl::InvalidArgumentError("`mod_params_q` must not be null.");
  }
  if (moduli.empty()) {
    return absl::InvalidArgumentError("`moduli` must not be empty.");
  }
  std::vector<std::vector<ModularInt>> coeff_vectors;
  coeff_vectors.reserve(moduli.size());
  for (int i = 0; i < moduli.size(); ++i) {
    RLWE_ASSIGN_OR_RETURN(
        std::vector<ModularInt> coeffs_qi,
        ConvertModulusBalanced<ModularInt>(coeffs_q, *mod_params_q,
                                           *moduli[i]->ModParams()));
    RLWE_RETURN_IF_ERROR(ForwardNumberTheoreticTransform(
        coeffs_qi, *moduli[i]->NttParams(), *moduli[i]->ModParams()));
    coeff_vectors.push_back(std::move(coeffs_qi));
  }
  return RnsPolynomial<ModularInt>(log2(num_coeffs), std::move(coeff_vectors),
                                   /*is_ntt=*/true);
}

template <typename ModularInt>
absl::Status RnsPolynomial<ModularInt>::ConvertToNttForm(
    absl::Span<const PrimeModulus<ModularInt>* const> moduli) {
  if (is_ntt_) {
    return absl::InvalidArgumentError("Polynomial already in NTT form");
  }
  int num_moduli = coeff_vectors_.size();
  if (moduli.size() != num_moduli) {
    return absl::InvalidArgumentError(
        absl::StrCat("`moduli` must contain ", num_moduli, " RNS moduli."));
  }
  for (int i = 0; i < num_moduli; ++i) {
    RLWE_RETURN_IF_ERROR(ForwardNumberTheoreticTransform(
        coeff_vectors_[i], *moduli[i]->NttParams(), *moduli[i]->ModParams()));
  }
  is_ntt_ = true;
  return absl::OkStatus();
}

template <typename ModularInt>
absl::Status RnsPolynomial<ModularInt>::ConvertToCoeffForm(
    absl::Span<const PrimeModulus<ModularInt>* const> moduli) {
  if (!is_ntt_) {
    return absl::InvalidArgumentError("Polynomial already in Coefficient form");
  }
  int num_moduli = coeff_vectors_.size();
  if (moduli.size() != num_moduli) {
    return absl::InvalidArgumentError(
        absl::StrCat("`moduli` must contain ", num_moduli, " RNS moduli."));
  }
  for (int i = 0; i < num_moduli; ++i) {
    RLWE_RETURN_IF_ERROR(InverseNumberTheoreticTransform(
        coeff_vectors_[i], *moduli[i]->NttParams(), *moduli[i]->ModParams()));
  }
  is_ntt_ = false;
  return absl::OkStatus();
}

template <typename ModularInt>
absl::StatusOr<RnsPolynomial<ModularInt>> RnsPolynomial<ModularInt>::Substitute(
    int power, absl::Span<const PrimeModulus<ModularInt>* const> moduli) const {
  int num_coeffs = 1 << log_n_;
  if (power < 0 || (power % 2) == 0 || power >= 2 * num_coeffs) {
    return absl::InvalidArgumentError(
        absl::StrCat("Substitution power must be a non-negative odd "
                     "integer less than 2*n."));
  }
  int num_moduli = coeff_vectors_.size();
  if (moduli.size() != num_moduli) {
    return absl::InvalidArgumentError(
        absl::StrCat("`moduli` must contain ", num_moduli, " RNS moduli."));
  }
  if (!is_ntt_) {
    return absl::InvalidArgumentError("RnsPolynomial must be in NTT form");
  }

  // The NTT representation of the polynomial modulo prime modulus q_i consists
  // of the evaluations of the polynomial at roots
  //    psi^brv[N/2], psi^brv[N/2+1], ..., psi^brv[N/2+N/2-1],
  //    psi^(N/2+brv[N/2+1]), ...,         psi^(N/2+brv[N/2+N/2-1]),
  // where psi is a primitive 2N-th root of unity in Z_{q_i}^*.
  // Thus, the NTT representation of f(X^power) is the evaluation of f(X) at
  // powers of psi^power.
  //
  // Get the index of the psi^power evaluation
  int psi_power_index = (power - 1) / 2;

  // Since substitution f(X) -> f(X^k) is an automorphism, we have that
  // f(X^k) (mod q_i) = f_i(X^k), where f_i(X) = f(X) (mod q_i).
  std::vector<std::vector<ModularInt>> subbed_coeff_vectors = coeff_vectors_;
  for (int i = 0; i < num_moduli; ++i) {
    const NttParams* ntt_params_qi = moduli[i]->NttParams();
    for (int j = 0; j < num_coeffs; j++) {
      int coeff_index = ntt_params_qi->bitrevs[psi_power_index];
      if (coeff_index >= num_coeffs) {
        return absl::InternalError(absl::StrFormat(
            "Index %d out-of-bounds in coeff_vectors_[%d] of size %d.",
            coeff_index, i, num_coeffs));
      }
      subbed_coeff_vectors[i][ntt_params_qi->bitrevs[j]] =
          coeff_vectors_[i][coeff_index];
      // Note that (psi^j)^power = psi^{(j * power) % (2 * N).
      // So each time the index j increases by 1, psi_power_index increases by
      // power mod N.
      psi_power_index = (psi_power_index + power) % num_coeffs;
    }
  }
  return RnsPolynomial<ModularInt>(log_n_, std::move(subbed_coeff_vectors),
                                   is_ntt_);
}

template <typename ModularInt>
absl::Status RnsPolynomial<ModularInt>::ModReduceLsb(
    typename ModularInt::Int t, const RnsInt<ModularInt>& ql_inv,
    absl::Span<const PrimeModulus<ModularInt>* const> moduli) {
  const int num_moduli = coeff_vectors_.size();
  if (moduli.size() != num_moduli) {
    return absl::InvalidArgumentError(
        absl::StrCat("`moduli` must contain ", num_moduli, " RNS moduli."));
  }
  if (num_moduli <= 1) {
    return absl::InvalidArgumentError(
        "ModReduceLsb cannot apply with only one prime modulus.");
  }
  const int num_moduli_reduced = num_moduli - 1;
  if (ql_inv.NumModuli() != num_moduli_reduced) {
    return absl::InvalidArgumentError(
        absl::StrCat("`ql_inv` must be defined with respect to ",
                     num_moduli_reduced, " moduli."));
  }

  // Detach the last sub-polynomial a_l = (a mod ql), and then compute a
  // correction term delta = t * [a_l * k (mod ql)] (mod Q/ql) for k and r such
  // that ql = k*t + r. The result of modulus reduction is (r * a + delta)/ql
  // (mod Q/ql), where r * a + delta is divisible by ql.
  RLWE_ASSIGN_OR_RETURN(std::vector<ModularInt> delta_ql_coeffs,
                        DetachLastCoeffVector());
  const ModularIntParams* mod_params_ql = moduli[num_moduli - 1]->ModParams();
  const NttParams* ntt_params_ql = moduli[num_moduli - 1]->NttParams();
  typename ModularInt::Int r = mod_params_ql->modulus % t;
  typename ModularInt::Int k = (mod_params_ql->modulus - r) / t;
  RLWE_ASSIGN_OR_RETURN(ModularInt k_mod_ql,
                        ModularInt::ImportInt(k, mod_params_ql));
  RLWE_RETURN_IF_ERROR(ModularInt::BatchMulInPlace(
      &delta_ql_coeffs, k_mod_ql, mod_params_ql));  // (a_l * k) (mod ql)

  // From now on we work on smaller RNS modulus Q/ql = q0 * .. * q{l-1}.
  absl::Span<const PrimeModulus<ModularInt>* const> moduli_reduced =
      moduli.first(num_moduli_reduced);
  RLWE_RETURN_IF_ERROR(MulInPlace(r, moduli_reduced));  // r * a (mod Q/ql).

  if (IsNttForm()) {
    RLWE_RETURN_IF_ERROR(InverseNumberTheoreticTransform(
        delta_ql_coeffs, *ntt_params_ql, *mod_params_ql));

    // Convert (a_l * k) from mod-ql to mod-(Q/ql).
    std::vector<std::vector<ModularInt>> delta_coeff_vectors;
    delta_coeff_vectors.reserve(num_moduli_reduced);
    for (int i = 0; i < num_moduli_reduced; ++i) {
      RLWE_ASSIGN_OR_RETURN(
          std::vector<ModularInt> delta_qi_coeffs,
          ConvertModulusBalanced(delta_ql_coeffs, *mod_params_ql,
                                 *moduli[i]->ModParams()));
      delta_coeff_vectors.push_back(std::move(delta_qi_coeffs));
    }
    RnsPolynomial<ModularInt> delta(log_n_, std::move(delta_coeff_vectors),
                                    /*is_ntt=*/false);
    RLWE_RETURN_IF_ERROR(delta.ConvertToNttForm(moduli_reduced));
    RLWE_RETURN_IF_ERROR(delta.MulInPlace(t, moduli_reduced));

    // Next we compute (r * a + delta) / ql.
    RLWE_RETURN_IF_ERROR(AddInPlace(delta, moduli_reduced));
    RLWE_RETURN_IF_ERROR(MulInPlace(ql_inv, moduli_reduced));
  } else {
    // Now this polynomial is in coeff form. In the following we compute
    // (r * a + delta) / ql wrt each prime moduli q0, .., q{l-1}.
    for (int i = 0; i < num_moduli_reduced; ++i) {
      const ModularIntParams* mod_params_qi = moduli[i]->ModParams();
      RLWE_ASSIGN_OR_RETURN(
          std::vector<ModularInt> delta_qi_coeffs,
          ConvertModulusBalanced<ModularInt>(delta_ql_coeffs, *mod_params_ql,
                                             *mod_params_qi));
      RLWE_ASSIGN_OR_RETURN(auto t_mod_qi,
                            ModularInt::ImportInt(t, mod_params_qi));
      RLWE_RETURN_IF_ERROR(ModularInt::BatchMulInPlace(
          &delta_qi_coeffs, t_mod_qi, mod_params_qi));
      RLWE_RETURN_IF_ERROR(ModularInt::BatchAddInPlace(
          &coeff_vectors_[i], delta_qi_coeffs, mod_params_qi));
      RLWE_RETURN_IF_ERROR(ModularInt::BatchMulInPlace(
          &coeff_vectors_[i], ql_inv.Component(i), mod_params_qi));
    }
  }
  return absl::OkStatus();
}

// Updates a RNS polynomial a (mod Q) to round(a / q_l) mod (Q/q_l), where Q =
// q_0 * .. * q_l. Assume a = b * q_l + a_l, where a_l = a mod q_l, then a - a_l
// is divisible by q_l and hence b = (a - a_l) / q_l (mod Q/q_l) = round(a/q_l).
template <typename ModularInt>
absl::Status RnsPolynomial<ModularInt>::ModReduceMsb(
    const RnsInt<ModularInt>& ql_inv,
    absl::Span<const PrimeModulus<ModularInt>* const> moduli) {
  const int num_moduli = coeff_vectors_.size();
  if (moduli.size() != num_moduli) {
    return absl::InvalidArgumentError(
        absl::StrCat("`moduli` must contain ", num_moduli, " RNS moduli."));
  }
  if (num_moduli <= 1) {
    return absl::InvalidArgumentError(
        "ModReduceMsb cannot apply with only one prime modulus.");
  }

  const int num_moduli_reduced = num_moduli - 1;
  if (ql_inv.NumModuli() != num_moduli_reduced) {
    return absl::InvalidArgumentError(
        absl::StrCat("`ql_inv` must be defined with respect to ",
                     num_moduli_reduced, " moduli."));
  }

  // Detach the last sub-polynomial [a]_ql = a (mod ql), and then compute a
  // correction term delta = a_l (mod Q/ql).
  RLWE_ASSIGN_OR_RETURN(std::vector<ModularInt> delta_ql_coeffs,
                        DetachLastCoeffVector());
  const ModularIntParams* mod_params_ql = moduli[num_moduli - 1]->ModParams();
  const NttParams* ntt_params_ql = moduli[num_moduli - 1]->NttParams();

  // From now on we work on smaller RNS modulus Q/ql = q0 * .. * q{l-1}.
  if (IsNttForm()) {
    RLWE_RETURN_IF_ERROR(InverseNumberTheoreticTransform(
        delta_ql_coeffs, *ntt_params_ql, *mod_params_ql));

    // Convert a_l from mod-ql to mod-(Q/ql).
    std::vector<std::vector<ModularInt>> delta_coeff_vectors;
    delta_coeff_vectors.reserve(num_moduli_reduced);
    for (int i = 0; i < num_moduli_reduced; ++i) {
      RLWE_ASSIGN_OR_RETURN(
          std::vector<ModularInt> delta_qi_coeffs,
          ConvertModulusBalanced(delta_ql_coeffs, *mod_params_ql,
                                 *moduli[i]->ModParams()));
      delta_coeff_vectors.push_back(std::move(delta_qi_coeffs));
    }
    RnsPolynomial<ModularInt> delta(log_n_, std::move(delta_coeff_vectors),
                                    /*is_ntt=*/false);
    absl::Span<const PrimeModulus<ModularInt>* const> moduli_reduced =
        moduli.first(num_moduli_reduced);
    RLWE_RETURN_IF_ERROR(delta.ConvertToNttForm(moduli_reduced));

    // We now have delta = [a (mod ql)] (mod Q/ql).
    RLWE_RETURN_IF_ERROR(SubInPlace(delta, moduli_reduced));   // (a - delta)
    RLWE_RETURN_IF_ERROR(MulInPlace(ql_inv, moduli_reduced));  // (a - delta)/ql
  } else {
    // Now this polynomial is in coeff form. In the following we compute
    // (a - delta) / ql wrt each prime moduli q0, .., q{l-1}.
    for (int i = 0; i < num_moduli_reduced; ++i) {
      const ModularIntParams* mod_params_qi = moduli[i]->ModParams();
      RLWE_ASSIGN_OR_RETURN(
          std::vector<ModularInt> delta_qi_coeffs,
          ConvertModulusBalanced<ModularInt>(delta_ql_coeffs, *mod_params_ql,
                                             *mod_params_qi));
      RLWE_RETURN_IF_ERROR(ModularInt::BatchSubInPlace(
          &coeff_vectors_[i], delta_qi_coeffs, mod_params_qi));
      RLWE_RETURN_IF_ERROR(ModularInt::BatchMulInPlace(
          &coeff_vectors_[i], ql_inv.Component(i), mod_params_qi));
    }
  }
  return absl::OkStatus();
}

// Switch to another CRT basis P = [p0, ..., p{k-1}] in NTT form
// Given (a_0,..,a_L) mod (q_0,..,q_L) which represents a \in [0,Q) for Q
// the product of q_0,..,q_L, by Chinese Reminder Theorem we have that
// \sum_i (a_i * [(Q/q_i)^(-1) mod q_i] * (Q/q_i)) = a (*).
// So, a can be represented by (b_1,..,b_K) wrt the CRT basis (p_1,..,p_K)
// where b_j = a mod p_j. The fast basis conversion is to avoid the expensive
// CRT interpolation process by computing b_j = a mod p_j. The result may not
// exactly represent a but a + e * Q, where |e| < L/2. Thus, if the modulus
// P = p_1 * .. * p_K is much larger than Q, then such error can be ignored
// and hence we have an approximate basis switching algorithm.
template <typename ModularInt>
absl::StatusOr<RnsPolynomial<ModularInt>>
RnsPolynomial<ModularInt>::ApproxSwitchRnsBasis(
    absl::Span<const PrimeModulus<ModularInt>* const> this_moduli,    // q_i's
    absl::Span<const PrimeModulus<ModularInt>* const> output_moduli,  // p_j's
    absl::Span<const ModularInt> prime_q_hat_inv_mod_qs,
    absl::Span<const RnsInt<ModularInt>> prime_q_hat_mod_ps,
    bool is_balanced_rep) const {
  const int num_this_moduli = coeff_vectors_.size();
  if (this_moduli.size() != num_this_moduli) {
    return absl::InvalidArgumentError(absl::StrCat(
        "`this_moduli` must contain ", num_this_moduli, " RNS moduli."));
  }
  if (output_moduli.empty()) {
    return absl::InvalidArgumentError("`output_moduli` must not be empty.");
  }
  if (prime_q_hat_inv_mod_qs.size() != num_this_moduli) {
    return absl::InvalidArgumentError(
        absl::StrCat("`prime_q_hat_inv_mod_qs` must contain ", num_this_moduli,
                     " elements."));
  }
  if (prime_q_hat_mod_ps.size() != num_this_moduli) {
    return absl::InvalidArgumentError(absl::StrCat(
        "`prime_q_hat_mod_ps` must contain ", num_this_moduli, " elements."));
  }

  std::vector<std::vector<ModularInt>> output_coeffs(output_moduli.size());
  for (int j = 0; j < output_moduli.size(); ++j) {
    output_coeffs[j].resize(
        NumCoeffs(), ModularInt::ImportZero(output_moduli[j]->ModParams()));
  }

  for (int i = 0; i < coeff_vectors_.size(); ++i) {
    const ModularIntParams* mod_params_qi = this_moduli[i]->ModParams();
    const NttParams* ntt_params_qi = this_moduli[i]->NttParams();
    // a_qi = a_i * (Q/q_i)^(-1) mod q_i.
    RLWE_ASSIGN_OR_RETURN(
        std::vector<ModularInt> a_coeffs_qi,
        ModularInt::BatchMul(coeff_vectors_[i], prime_q_hat_inv_mod_qs[i],
                             mod_params_qi));
    if (IsNttForm()) {
      RLWE_RETURN_IF_ERROR(InverseNumberTheoreticTransform(
          a_coeffs_qi, *ntt_params_qi, *mod_params_qi));
    }

    for (int j = 0; j < output_moduli.size(); ++j) {
      const ModularIntParams* mod_params_pj = output_moduli[j]->ModParams();
      const NttParams* ntt_params_pj = output_moduli[j]->NttParams();
      // Compute a_pj = ([a_i * (Q/q_i)^(-1) mod q_i] * Q/q_i) (mod p_j).
      RLWE_ASSIGN_OR_RETURN(
          std::vector<ModularInt> a_coeffs_pj,
          is_balanced_rep
              ? ConvertModulusBalanced(a_coeffs_qi, *mod_params_qi,
                                       *mod_params_pj)
              : ConvertModulus(a_coeffs_qi, *mod_params_qi, *mod_params_pj));
      RLWE_RETURN_IF_ERROR(ForwardNumberTheoreticTransform(
          a_coeffs_pj, *ntt_params_pj, *mod_params_pj));
      RLWE_RETURN_IF_ERROR(ModularInt::BatchMulInPlace(
          &a_coeffs_pj, prime_q_hat_mod_ps[i].Component(j), mod_params_pj));
      RLWE_RETURN_IF_ERROR(ModularInt::BatchAddInPlace(
          &output_coeffs[j], a_coeffs_pj, mod_params_pj));
    }
  }
  // Return the RNS polynomial in NTT form
  return RnsPolynomial<ModularInt>(log_n_, std::move(output_coeffs),
                                   /*is_ntt=*/true);
}

template <typename ModularInt>
absl::Status RnsPolynomial<ModularInt>::ApproxModReduceLsb(
    RnsPolynomial<ModularInt> polynomial_aux,  // mod P = mod (p_1 * .. p_K)
    absl::Span<const PrimeModulus<ModularInt>* const> this_moduli,  // q_i's
    absl::Span<const PrimeModulus<ModularInt>* const> aux_moduli,   // p_j's
    absl::Span<const ModularInt> p_inv_mod_qs,
    absl::Span<const ModularInt> prime_p_hat_inv_mod_ps,
    absl::Span<const RnsInt<ModularInt>> prime_p_hat_mod_qs,
    absl::Span<const ModularInt> t_inv_mod_ps, typename ModularInt::Int p_mod_t,
    typename ModularInt::Int t) {
  const int num_this_moduli = coeff_vectors_.size();
  if (this_moduli.size() != num_this_moduli) {
    return absl::InvalidArgumentError(absl::StrCat(
        "`this_moduli` must contain ", num_this_moduli, " RNS moduli."));
  }
  const int num_aux_moduli = polynomial_aux.coeff_vectors_.size();
  if (aux_moduli.size() != num_aux_moduli) {
    return absl::InvalidArgumentError(absl::StrCat(
        "`aux_moduli` must contain ", num_aux_moduli, " RNS moduli."));
  }
  if (p_inv_mod_qs.size() != num_this_moduli) {
    return absl::InvalidArgumentError(absl::StrCat(
        "`p_inv_mod_qs` must contain ", num_this_moduli, " elements."));
  }
  if (prime_p_hat_inv_mod_ps.size() != num_aux_moduli) {
    return absl::InvalidArgumentError(
        absl::StrCat("`prime_p_hat_inv_mod_ps` must contain ", num_aux_moduli,
                     " elements."));
  }
  if (prime_p_hat_mod_qs.size() != num_aux_moduli) {
    return absl::InvalidArgumentError(absl::StrCat(
        "`prime_p_hat_mod_qs` must contain ", num_aux_moduli, " elements."));
  }
  if (t_inv_mod_ps.size() != num_aux_moduli) {
    return absl::InvalidArgumentError(absl::StrCat(
        "`t_inv_mod_ps` must contain ", num_aux_moduli, " elements."));
  }

  // Compute delta = [k * polynomial_aux] mod P, where P = k * t + p_mod_t.
  // Since k = (P - p_mod_t) / t, we have -p_mod_t * t^(-1) = k (mod p_j).
  for (int j = 0; j < num_aux_moduli; ++j) {
    auto mod_params_pj = aux_moduli[j]->ModParams();
    RLWE_ASSIGN_OR_RETURN(ModularInt r,
                          ModularInt::ImportInt(p_mod_t, mod_params_pj));
    RLWE_RETURN_IF_ERROR(ModularInt::BatchMulInPlace(
        &polynomial_aux.coeff_vectors_[j], r, mod_params_pj));
    RLWE_RETURN_IF_ERROR(ModularInt::BatchMulInPlace(
        &polynomial_aux.coeff_vectors_[j],
        t_inv_mod_ps[j].Negate(mod_params_pj), mod_params_pj));
  }

  // Switch delta from (mod P) to (mod Q). The result is in NTT form.
  RLWE_ASSIGN_OR_RETURN(
      RnsPolynomial<ModularInt> delta,
      polynomial_aux.ApproxSwitchRnsBasis(
          aux_moduli, this_moduli, prime_p_hat_inv_mod_ps, prime_p_hat_mod_qs,
          /*is_balanced_rep=*/true));

  // delta = t * [k * polynomial_aux (mod P)] (mod Q).
  RLWE_RETURN_IF_ERROR(delta.MulInPlace(t, this_moduli));

  // Make sure this RNS polynomial is in NTT form.
  if (!IsNttForm()) {
    RLWE_RETURN_IF_ERROR(ConvertToNttForm(this_moduli));
  }

  // P^(-1) * (p_mod_t * this polynomial + delta) (mod Q).
  RLWE_RETURN_IF_ERROR(MulInPlace(p_mod_t, this_moduli));
  RLWE_RETURN_IF_ERROR(AddInPlace(delta, this_moduli));
  RLWE_RETURN_IF_ERROR(MulInPlace(p_inv_mod_qs, this_moduli));

  return absl::OkStatus();
}

template <typename ModularInt>
absl::Status RnsPolynomial<ModularInt>::ApproxModReduceMsb(
    RnsPolynomial<ModularInt> polynomial_aux,
    absl::Span<const PrimeModulus<ModularInt>* const> this_moduli,  // q_i's
    absl::Span<const PrimeModulus<ModularInt>* const> aux_moduli,   // p_j's
    absl::Span<const ModularInt> p_inv_mod_qs,
    absl::Span<const ModularInt> prime_p_hat_inv_mod_ps,
    absl::Span<const RnsInt<ModularInt>> prime_p_hat_mod_qs) {
  const int num_this_moduli = coeff_vectors_.size();
  if (this_moduli.size() != num_this_moduli) {
    return absl::InvalidArgumentError(absl::StrCat(
        "`this_moduli` must contain ", num_this_moduli, " RNS moduli."));
  }
  const int num_aux_moduli = polynomial_aux.coeff_vectors_.size();
  if (aux_moduli.size() != num_aux_moduli) {
    return absl::InvalidArgumentError(absl::StrCat(
        "`aux_moduli` must contain ", num_aux_moduli, " RNS moduli."));
  }
  if (p_inv_mod_qs.size() != num_this_moduli) {
    return absl::InvalidArgumentError(absl::StrCat(
        "`p_inv_mod_qs` must contain ", num_this_moduli, " elements."));
  }
  if (prime_p_hat_inv_mod_ps.size() != num_aux_moduli) {
    return absl::InvalidArgumentError(
        absl::StrCat("`prime_p_hat_inv_mod_ps` must contain ", num_aux_moduli,
                     " elements."));
  }
  if (prime_p_hat_mod_qs.size() != num_aux_moduli) {
    return absl::InvalidArgumentError(absl::StrCat(
        "`prime_p_hat_mod_qs` must contain ", num_aux_moduli, " elements."));
  }

  // delta = polynomial_aux (mod Q).
  RLWE_ASSIGN_OR_RETURN(
      RnsPolynomial<ModularInt> delta,
      polynomial_aux.ApproxSwitchRnsBasis(
          aux_moduli, this_moduli, prime_p_hat_inv_mod_ps, prime_p_hat_mod_qs,
          /*is_balanced_rep=*/true));

  // Make sure this RNS polynomial is in NTT form.
  if (!IsNttForm()) {
    RLWE_RETURN_IF_ERROR(ConvertToNttForm(this_moduli));
  }

  // P^(-1) * (this polynomial - delta) (mod Q).
  RLWE_RETURN_IF_ERROR(SubInPlace(delta, this_moduli));
  RLWE_RETURN_IF_ERROR(MulInPlace(p_inv_mod_qs, this_moduli));

  return absl::OkStatus();
}

template class RnsPolynomial<MontgomeryInt<Uint16>>;
template class RnsPolynomial<MontgomeryInt<Uint32>>;
template class RnsPolynomial<MontgomeryInt<Uint64>>;
template class RnsPolynomial<MontgomeryInt<absl::uint128>>;
#ifdef ABSL_HAVE_INTRINSIC_INT128
template class RnsPolynomial<MontgomeryInt<unsigned __int128>>;
#endif

}  // namespace rlwe
