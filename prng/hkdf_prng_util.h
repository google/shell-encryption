/*
 * Copyright 2021 Google LLC
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

// An implementation of a PRNG using a HMAC-based key derivation function.
//
// HMAC-based key derivation functions (HKDF, for short) consist of two
// important functions: extract and expand. Given an input key with sufficient
// entropy, the HKDF will extract the entropy into a more uniform, unbiased
// entropy. The HKDF can expand this entropy into many pseudorandom outputs.
// Therefore, the input key must have sufficient entropy to ensure the outputs
// are pseudorandom. The pseudorandom outputs of HKDF are deterministic for any
// fixed input key allowing replay of the pseudorandom outputs by multiple
// clients by sharing the input key. For more information about HKDFs, see [1]
// for an overview and [2] for a full description.
//
// [1] https://en.wikipedia.org/wiki/HKDF
// [2] https://tools.ietf.org/html/rfc5869

#ifndef RLWE_HKDF_PRNG_UTIL_H_
#define RLWE_HKDF_PRNG_UTIL_H_

#include "absl/status/status.h"
#include "absl/strings/string_view.h"
#include "integral_types.h"
#include "statusor.h"

namespace rlwe::internal {

const int kHkdfKeyBytesSize = 64;
const int kHkdfMaxOutputBytes = 255 * 32;

// Once pseudorandom output is exhausted, the salt is updated to construct
// new pseudorandom output.
absl::Status HkdfPrngResalt(absl::string_view key, int buffer_size,
                            int* salt_counter, int* position_in_buffer,
                            std::vector<Uint8>* buffer);

// Generates a secure key for instantiating an HKDF.
rlwe::StatusOr<std::string> HkdfPrngGenerateKey();

// Returns 8 bits of randomness.
//
// Fails on internal cryptographic errors.
rlwe::StatusOr<Uint8> HkdfPrngRand8(absl::string_view key,
                                    int* position_in_buffer, int* salt_counter,
                                    std::vector<Uint8>* buffer);

// Returns 64 bits of randomness.
//
// Fails on internal cryptographic errors.
rlwe::StatusOr<Uint64> HkdfPrngRand64(absl::string_view key,
                                      int* position_in_buffer,
                                      int* salt_counter,
                                      std::vector<Uint8>* buffer);

}  // namespace rlwe::internal

#endif  // RLWE_HKDF_PRNG_UTIL_H_
