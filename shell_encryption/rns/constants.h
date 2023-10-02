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

#ifndef RLWE_RNS_CONSTANTS_H_
#define RLWE_RNS_CONSTANTS_H_

#include "absl/numeric/int128.h"
#include "shell_encryption/integral_types.h"

namespace rlwe {

// This file declares some commonly used constants for the RNS variant of the
// RLWE encryption scheme. We consider power-of-2 cyclotomic ring Z[X]/(X^N+1)
// and its residue ring Z[X]/(Q*P,X^N+1) for a main modulus Q = q_0 *..* q_{L-1}
// and an auxiliary modulus P = p_0 *..* p_{K-1}, where q_i and p_j are distinct
// primes satisfying q_i == p_j == 1 (mod 2*N) to enable radix-2 NTT.
// Furthermore, to enable modulus reduction for the BGV scheme with minimal
// noise growth, the plaintext modulus t satisfies that q_i == 1 (mod t).

// Prime moduli for a 60-bit modulus Q*P and N = 2^11
constexpr Uint32 k60BitModulusQ0 = 974849;  // 20-bits
constexpr Uint32 k60BitModulusQ1 = 765953;  // 20-bits
constexpr Uint32 k60BitModulusP0 = 557057;  // 20-bits
constexpr Uint32 k60BitModulusLogN = 11;    // N = 2^11
constexpr Uint32 k60BitModulusPlaintextModulus = 17;

// Prime moduli for a 59-bit modulus Q*P and N = 2^11.
constexpr Uint32 k59BitModulusQ0 = 536375297;   // 29-bits
constexpr Uint32 k59BitModulusP0 = 1073655809;  // 30-bits
constexpr Uint32 k59BitModulusLogN = 11;
constexpr Uint32 k59BitModulusPlaintextModulus = 17;

}  // namespace rlwe

#endif  // RLWE_RNS_CONSTANTS_H_
