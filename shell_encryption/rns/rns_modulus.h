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

#ifndef RLWE_RNS_RNS_MODULUS_H_
#define RLWE_RNS_RNS_MODULUS_H_

#include "shell_encryption/montgomery.h"
#include "shell_encryption/ntt_parameters.h"

namespace rlwe {

// A RNS modulus Q = q_0 * ... * q_L is the product of prime moduli q_i's.
// By Chinese Remainder Theorem, we have a ring isomorphism between Z_Q and
// Z_{q_0} * ... * Z_{q_L} and so, addition and multiplication mod-Q can be done
// with respect to each q_i independently. This type is an aggregate that holds
// the Montgomery and NTT parameters for a prime modulus q_i that define the
// arithmetic mod q_i and mod (X^N+1) for N a power of two.
template <typename ModularInt>
struct PrimeModulus {
  using ModularIntParams = typename ModularInt::Params;

  typename ModularInt::Int Modulus() const { return mod_params->modulus; }
  const ModularIntParams* ModParams() const { return mod_params.get(); }
  const NttParameters<ModularInt>* NttParams() const {
    return ntt_params.get();
  }

  std::unique_ptr<const ModularIntParams> mod_params;
  std::unique_ptr<const NttParameters<ModularInt>> ntt_params;
};

}  // namespace rlwe

#endif  // RLWE_RNS_RNS_MODULUS_H_
