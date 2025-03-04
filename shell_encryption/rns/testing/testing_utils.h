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

// This file contains commonly used utilities for testing the RNS extension.

#ifndef RLWE_RNS_TESTING_TESTING_UTILS_H_
#define RLWE_RNS_TESTING_TESTING_UTILS_H_

#include <complex>
#include <cstdint>
#include <vector>

#include "absl/log/check.h"
#include "absl/numeric/int128.h"
#include "absl/random/random.h"
#include "absl/types/span.h"
#include "shell_encryption/integral_types.h"
#include "shell_encryption/rns/rns_modulus.h"
#include "shell_encryption/rns/rns_polynomial.h"

namespace rlwe {
namespace testing {

// Returns a vector of `num_values` many random integers in [0, max_value).
// Assumes that `max_value` is smaller than 2^64.
template <typename Integer>
std::vector<Integer> SampleMessages(int num_values, Integer max_value) {
  absl::BitGen bitgen;
  std::vector<Integer> messages;
  for (int i = 0; i < num_values; ++i) {
    Integer message = static_cast<Integer>(
        absl::Uniform<uint64_t>(bitgen, 0, static_cast<uint64_t>(max_value)));
    messages.push_back(message);
  }
  return messages;
}

inline std::vector<std::complex<double>> SampleComplexValues(
    int num_values, uint64_t max_value, uint64_t precision) {
  absl::BitGen bitgen;
  std::vector<std::complex<double>> values;
  for (int i = 0; i < num_values; ++i) {
    uint64_t u = absl::Uniform<uint64_t>(bitgen, 0, max_value * 2 * precision);
    uint64_t v = absl::Uniform<uint64_t>(bitgen, 0, max_value * 2 * precision);
    int64_t real = static_cast<int64_t>(u) - max_value * precision;
    int64_t imag = static_cast<int64_t>(v) - max_value * precision;
    std::complex<double> value{static_cast<double>(real) / precision,
                               static_cast<double>(imag) / precision};
    values.push_back(value);
  }
  return values;
}

// Returns a RnsPolynomial in Z[X]/(Q, X^N+1) whose coefficients are random
// 0, 1, or -1 (mod Q).
template <typename ModularInt>
RnsPolynomial<ModularInt> SampleTernaryNoises(
    int log_n, absl::Span<const PrimeModulus<ModularInt>* const> moduli) {
  int num_coeffs = 1 << log_n;
  int num_moduli = moduli.size();
  absl::BitGen bitgen;
  std::vector<std::vector<ModularInt>> coeff_vectors(num_moduli);
  for (int i = 0; i < num_coeffs; ++i) {
    int8_t e = absl::Uniform<int8_t>(absl::IntervalClosed, bitgen, -1, 1);
    for (int j = 0; j < num_moduli; ++j) {
      auto mod_params_qj = moduli[j]->ModParams();
      auto e_mod_qj = ModularInt::ImportInt(
          e == -1 ? mod_params_qj->modulus : (e == 0 ? 0 : 1), mod_params_qj);
      CHECK(e_mod_qj.ok());
      coeff_vectors[j].push_back(std::move(e_mod_qj.value()));
    }
  }
  auto a = RnsPolynomial<ModularInt>::Create(coeff_vectors, /*is_ntt=*/false);
  CHECK(a.ok());
  return a.value();
}

}  // namespace testing
}  // namespace rlwe

#endif  // RLWE_RNS_TESTING_TESTING_UTILS_H_
