/*
 * Copyright 2025 Google LLC.
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

#ifndef RLWE_MULTI_PARTY_TESTING_PARAMETERS_H_
#define RLWE_MULTI_PARTY_TESTING_PARAMETERS_H_

#include <vector>

#include <gtest/gtest.h>
#include "shell_encryption/integral_types.h"
#include "shell_encryption/montgomery.h"

namespace rlwe {
namespace multi_party {
namespace testing {

using ModularInt32 = MontgomeryInt<Uint32>;
using ModularInt64 = MontgomeryInt<Uint64>;
using ModularIntTypesForMultiParty =
    ::testing::Types<rlwe::MontgomeryInt<Uint32>, rlwe::MontgomeryInt<Uint64>>;
using ModularIntTypesForNegativeTests = ::testing::Types<ModularInt64>;

// Aggregate type that holds the parameters defining multi-party BFV protocol.
template <typename ModularInt>
struct MpaheParameters {
  int log_n;
  std::vector<typename ModularInt::Int> qs;  // main prime moduli.
  std::vector<typename ModularInt::Int> ps;  // auxiliary prime moduli.
  typename ModularInt::Int t;                // plaintext modulus.
  int log_gadget_base;                       // Bit size of gadget base.
  double s_flood;  // Gaussian parameter for flooding noise.
};

// Returns the testing parameters for the underlying integer type.
template <typename ModularInt>
MpaheParameters<ModularInt> GetMultiPartyParameters();

// Returns the testing parameters when instantiating RNS using 32-bit integers.
// The prime moduli in `qs` must be NTT-friendly for the higher order cyclotomic
// X^{2N} + 1 to allow efficient computation of the "wrap around" polynomial for
// the public key share.
template <>
inline MpaheParameters<ModularInt32> GetMultiPartyParameters<ModularInt32>() {
  return
      // 60 bits main modulus.
      MpaheParameters<ModularInt32>{.log_n = 11,
                                    .qs = {1073692673, 1073643521},
                                    .ps = {},
                                    .t = 10001,
                                    .log_gadget_base = 5,
                                    .s_flood = 1.0e+10};
}

// Returns the testing parameters when instantiating RNS using 64-bit integers.
// Same as above, the prime moduli in `qs` must be NTT friendly for X^{2N} + 1.
template <>
inline MpaheParameters<ModularInt64> GetMultiPartyParameters<ModularInt64>() {
  return
      // 69 bits main modulus.
      MpaheParameters<ModularInt64>{.log_n = 12,
                                    .qs = {34359410689ULL, 34359214081ULL},
                                    .ps = {},
                                    .t = 54001,
                                    .log_gadget_base = 5,
                                    .s_flood = 4.25839e+13};
}

}  // namespace testing
}  // namespace multi_party
}  // namespace rlwe

#endif  // RLWE_MULTI_PARTY_TESTING_PARAMETERS_H_
