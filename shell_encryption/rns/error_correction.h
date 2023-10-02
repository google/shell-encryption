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

#ifndef RLWE_RNS_ERROR_CORRECTION_H_
#define RLWE_RNS_ERROR_CORRECTION_H_

#include <utility>
#include <vector>

#include "absl/status/statusor.h"

namespace rlwe {

// Given `noisy_coeffs` representing integers m + t * e using mod-q values,
// where m in [0, t) and the error t * e are in [-q/2, q/2), returns the
// messages m. This implements the error correction for BGV where the errors
// reside in most significant bits of ciphertext modulus q.
//
// This is essentially the same algorithm as `rlwe::RemoveError()` defined in
// symmetric_encryption.h used by the non-RNS version of rlwe. The difference is
// that this version is parameterized by an integral type that support modulo
// reduction, so it is possible to compute CRT interpolation of noisy plaintext
// RNS polynomial based on small integer types and then error correct it using
// a bigger integer type.
template <typename Integer>
absl::StatusOr<std::vector<Integer>> RemoveErrorOnMsb(
    std::vector<Integer> noisy_coeffs, Integer q, Integer t) {
  Integer zero{0};
  if (t == zero) {
    return absl::InvalidArgumentError("`t` cannot be zero.");
  }

  Integer q_mod_t = q % t;
  Integer q_half = q >> 1;
  for (int i = 0; i < noisy_coeffs.size(); ++i) {
    if (noisy_coeffs[i] > q_half) {
      noisy_coeffs[i] -= q_mod_t;
    }
    noisy_coeffs[i] %= t;
  }
  return noisy_coeffs;
}

}  // namespace rlwe

#endif  // RLWE_RNS_ERROR_CORRECTION_H_
