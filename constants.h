/*
 * Copyright 2017 Google Inc.
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

#ifndef RLWE_CONSTANTS_H_
#define RLWE_CONSTANTS_H_

#include <cstdint>

namespace rlwe {

// Parameters from the New Hope key exchange protocol. Note that we are using
// these parameters for symmetric key encryption, not key exchange.
constexpr uint64_t kNewhopeModulus = 12289;
constexpr uint64_t kNewhopeLogR = 18;
constexpr uint64_t kNewhopeLogDegreeBound = 10;
constexpr uint64_t kNewhopeDegreeBound = 1 << kNewhopeLogDegreeBound;

// Montgomery parameters for a 60-bit modulus.
constexpr uint64_t kModulus60 = 332366567264636929;
constexpr uint64_t kInvModulus60 = 7124357790306815999;
constexpr uint16_t kLogR60 = 64;
constexpr uint64_t kInvR60 = 128364026370599367;

// RLWE parameters for a 30-bit modulus.
constexpr uint64_t kModulus30 = 536881153;
constexpr uint64_t kLogR30 = 31;
constexpr uint64_t kLogDegreeBound30 = 10;
constexpr uint64_t kDegreeBound30 = 1L << kLogDegreeBound30;

}  // namespace rlwe

#endif  // RLWE_CONSTANTS_H_
