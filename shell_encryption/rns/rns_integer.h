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

#ifndef RLWE_RNS_RNS_INTEGER_H_
#define RLWE_RNS_RNS_INTEGER_H_

#include <vector>

#include "absl/types/span.h"

namespace rlwe {

// This class stores an integer z (mod Q) using RNS representation, where
// Q = q_0 * ... * q_{L-1} for distinct primes q_i. The RNS representation
// of z is the vector [z (mod q_i) for i = 0 .. L-1].
// We assume that L >= 1 and the RNS integer z is not vacuous.
template <typename ModularInt>
struct RnsInt {
  std::vector<ModularInt> zs;

  // The number of prime moduli for this RNS integer.
  int NumModuli() const { return zs.size(); }

  // Returns the RnsInt represented by the first `num_moduli` residue values.
  RnsInt Prefix(int num_moduli) const {
    std::vector<ModularInt> zs_prefix(zs.begin(), zs.begin() + num_moduli);
    return RnsInt{std::move(zs_prefix)};
  }

  // Returns the vector of residue values representing this RNS integer.
  absl::Span<const ModularInt> RnsRep() const { return zs; }

  // Returns the vector of the first `num_moduli` residue values.
  absl::Span<const ModularInt> RnsRepPrefix(int num_moduli) const {
    return absl::MakeSpan(zs).subspan(0, num_moduli);
  }

  // Accessor to the residue value at the given level index.
  const ModularInt& Component(int index) const { return zs[index]; }
};

}  // namespace rlwe

#endif  // RLWE_RNS_RNS_INTEGER_H_
