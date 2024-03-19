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

#ifndef RLWE_SAMPLER_UNIFORM_TERNARY_H_
#define RLWE_SAMPLER_UNIFORM_TERNARY_H_

#include <algorithm>
#include <vector>

#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "shell_encryption/integral_types.h"
#include "shell_encryption/prng/prng.h"
#include "shell_encryption/status_macros.h"

namespace rlwe {

// Samples a vector of `num_coeffs` coefficients from the uniform ternary
// distribution on {-1, 0, 1}, where coefficients are reduced modulo q as given
// in `mod_params`.
template <typename ModularInt>
static absl::StatusOr<std::vector<ModularInt>> SampleFromUniformTernary(
    int num_coeffs, const typename ModularInt::Params* mod_params,
    SecurePrng* prng) {
  using Integer = typename ModularInt::Int;
  if (num_coeffs <= 0) {
    return absl::InvalidArgumentError("`num_coeffs` must be positive.");
  }
  if (mod_params == nullptr) {
    return absl::InvalidArgumentError("`mod_params` must not be null.");
  }
  if (prng == nullptr) {
    return absl::InvalidArgumentError("`prng` must not be null.");
  }

  // Below implements a parallel rejection sampling algorithm as suggested in
  // "A constant-time sampler for close-to-uniform bitsliced ternary vectors"
  // by Pierre Karpman, https://hal.archives-ouvertes.fr/hal-03777885
  // An element from {-1, 0, 1} is represented using two bits: 0 as (0,0), 1 as
  // (1,0), and -1 as (1,1). The algorithm samples in batches two uniformly
  // random 8-bit integers, r0 and r1, and uses a mask to indicate if a certain
  // bit in r0 and r1 is the invalid representation (0,1) and needs re-sample.
  Integer plus{1};
  Integer minus = mod_params->modulus - 1;
  std::vector<ModularInt> coeffs;
  coeffs.reserve(num_coeffs);
  while (num_coeffs > 0) {
    int num_filled_coeffs = std::min(num_coeffs, 8);
    // The mask to indicate if a bit index requires re-sampling.
    Uint8 missing_bits =
        num_filled_coeffs < 8 ? (1 << num_filled_coeffs) - 1 : ~0;
    Uint8 encoding_bits0 = 0;
    Uint8 encoding_bits1 = 0;
    while (missing_bits != 0) {
      RLWE_ASSIGN_OR_RETURN(Uint8 rand_bits0, prng->Rand8());
      RLWE_ASSIGN_OR_RETURN(Uint8 rand_bits1, prng->Rand8());
      encoding_bits0 ^= (rand_bits0 & missing_bits);
      encoding_bits1 ^= (rand_bits1 & missing_bits);
      missing_bits = ~encoding_bits0 & encoding_bits1;
    }

    for (int i = 0; i < num_filled_coeffs; ++i) {
      Uint8 mask = 1 << i;
      bool bit0 = ((encoding_bits0 & mask) > 0);
      bool bit1 = ((encoding_bits1 & mask) > 0);
      Integer is_plus = -static_cast<Integer>(bit0 && !bit1);
      Integer is_minus = -static_cast<Integer>(bit0 && bit1);
      Integer value = (is_plus & plus) | (is_minus & minus);
      RLWE_ASSIGN_OR_RETURN(auto coeff,
                            ModularInt::ImportInt(value, mod_params));
      coeffs.push_back(std::move(coeff));
    }
    num_coeffs -= num_filled_coeffs;
  }
  return coeffs;
}

}  // namespace rlwe

#endif  // RLWE_SAMPLER_UNIFORM_TERNARY_H_
