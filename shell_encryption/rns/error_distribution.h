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

#ifndef RLWE_RNS_ERROR_DISTRIBUTION_H_
#define RLWE_RNS_ERROR_DISTRIBUTION_H_

#include <vector>

#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "absl/types/span.h"
#include "shell_encryption/prng/prng.h"
#include "shell_encryption/rns/rns_modulus.h"
#include "shell_encryption/rns/rns_polynomial.h"
#include "shell_encryption/sample_error.h"
#include "shell_encryption/sampler/discrete_gaussian.h"
#include "shell_encryption/sampler/uniform_ternary.h"
#include "shell_encryption/status_macros.h"

namespace rlwe {

// Samples from the error distribution, where the returned sample is represented
// as a RNS polynomial mod Q, whose prime decomposition is given in `moduli`.
// The error polynomial has independent coefficients sampled from a centered
// binomial distribution of the given variance, such that the absolute value of
// every coefficient is smaller than a single prime modulus.
template <typename ModularInt>
absl::StatusOr<RnsPolynomial<ModularInt>> SampleError(
    int log_n, int variance,
    absl::Span<const PrimeModulus<ModularInt>* const> moduli,
    SecurePrng* prng) {
  if (log_n <= 0) {
    return absl::InvalidArgumentError("`log_n` must be positive.");
  }
  if (variance <= 0) {
    return absl::InvalidArgumentError("`variance` must be positive.");
  }
  if (moduli.empty()) {
    return absl::InvalidArgumentError("`moduli` must not be empty.");
  }

  int num_coeffs = 1 << log_n;
  const typename ModularInt::Params* mod_params_q0 = moduli[0]->ModParams();

  // Sample coefficients of e (mod q0).
  RLWE_ASSIGN_OR_RETURN(std::vector<ModularInt> e_coeffs_q0,
                        SampleFromErrorDistribution<ModularInt>(
                            num_coeffs, variance, prng, mod_params_q0));
  // Convert e (mod q0) to e (mod Q).
  return RnsPolynomial<ModularInt>::ConvertBalancedFromPolynomialCoeffs(
      e_coeffs_q0, mod_params_q0, moduli);
}

// Samples a RNS polynomial with uniform ternary coefficients {-1, 0, 1} modulo
// the RNS modulus Q given in `moduli`.
template <typename ModularInt>
absl::StatusOr<RnsPolynomial<ModularInt>> SampleUniformTernary(
    int log_n, absl::Span<const PrimeModulus<ModularInt>* const> moduli,
    SecurePrng* prng) {
  if (log_n <= 0) {
    return absl::InvalidArgumentError("`log_n` must be positive.");
  }
  if (moduli.empty()) {
    return absl::InvalidArgumentError("`moduli` must not be empty.");
  }
  // Sample coefficients of v (mod q0), and then convert to v (mod Q).
  const typename ModularInt::Params* mod_params_q0 = moduli[0]->ModParams();
  RLWE_ASSIGN_OR_RETURN(
      std::vector<ModularInt> v_coeffs_q0,
      SampleFromUniformTernary<ModularInt>(1 << log_n, mod_params_q0, prng));
  return RnsPolynomial<ModularInt>::ConvertBalancedFromPolynomialCoeffs(
      v_coeffs_q0, mod_params_q0, moduli);
}

// Samples a RNS polynomial with independent discrete Gaussian coefficients
// with Gaussian parameter `s`. That is, each coefficient c has probability
// proportional to exp(-c^2 / (2 * s^2)).
// Note that coefficients are represented modulo Q as given in `moduli`, where
// negative c's are represented as Q - abs(c).
template <typename ModularInt>
absl::StatusOr<RnsPolynomial<ModularInt>> SampleDiscreteGaussian(
    int log_n, double s,
    absl::Span<const PrimeModulus<ModularInt>* const> moduli,
    const DiscreteGaussianSampler<typename ModularInt::Int>* dg_sampler,
    SecurePrng* prng) {
  using Integer = typename ModularInt::Int;
  if (log_n <= 0) {
    return absl::InvalidArgumentError("`log_n` must be positive.");
  }
  if (moduli.empty()) {
    return absl::InvalidArgumentError("`moduli` must not be empty.");
  }
  if (dg_sampler == nullptr) {
    return absl::InvalidArgumentError("`dg_sampler` must not be null.");
  }

  // Precompute the number of iterations needed to run the sampling algorithm.
  RLWE_ASSIGN_OR_RETURN(int num_iterations, dg_sampler->NumIterations(s));

  int num_coeffs = 1 << log_n;
  int num_moduli = moduli.size();
  std::vector<std::vector<ModularInt>> coeff_vectors(num_moduli);
  for (int j = 0; j < num_moduli; ++j) {
    coeff_vectors[j].reserve(num_coeffs);
  }
  for (int i = 0; i < num_coeffs; ++i) {
    RLWE_ASSIGN_OR_RETURN(Integer coeff, dg_sampler->SampleWithIterations(
                                             s, num_iterations, *prng));
    bool is_negative =
        coeff > DiscreteGaussianSampler<Integer>::kNegativeThreshold;

    for (int j = 0; j < num_moduli; ++j) {
      auto qj = moduli[j]->Modulus();
      Integer coeff_mod_qj =
          is_negative ? qj - (static_cast<Integer>(-coeff) % qj) : coeff;
      RLWE_ASSIGN_OR_RETURN(
          ModularInt coeff_qj,
          ModularInt::ImportInt(coeff_mod_qj, moduli[j]->ModParams()));
      coeff_vectors[j].push_back(std::move(coeff_qj));
    }
  }
  return RnsPolynomial<ModularInt>::Create(std::move(coeff_vectors),
                                           /*is_ntt=*/false);
}

}  // namespace rlwe

#endif  // RLWE_RNS_ERROR_DISTRIBUTION_H_
