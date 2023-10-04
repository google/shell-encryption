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

}  // namespace rlwe

#endif  // RLWE_RNS_ERROR_DISTRIBUTION_H_
