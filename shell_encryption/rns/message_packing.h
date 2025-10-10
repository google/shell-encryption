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

#ifndef RLWE_RNS_MESSAGE_PACKING_H_
#define RLWE_RNS_MESSAGE_PACKING_H_

#include <algorithm>
#include <vector>

#include "absl/types/span.h"

namespace rlwe {

namespace internal {

// Returns ceil(x / y).
template <typename T>
T DivideAndRoundUp(const T& x, const T& y) {
  return (x + y - 1) / y;
}

}  // namespace internal

// Packs integers in [0, `packing_base`) from `messages` into a smaller number
// of larger integers in [0, `packing_base`^`packing_dimension`). Expects
// packing_base > 1, packing_dimension > 0, and
// packing_base^packing_dimension < std::numeric_limits<BigInteger>::max().
template <typename Integer, typename BigInteger>
std::vector<BigInteger> PackMessagesFlat(absl::Span<const Integer> messages,
                                         Integer packing_base,
                                         int packing_dimension) {
  int num_packed_coeffs =
      internal::DivideAndRoundUp<int>(messages.size(), packing_dimension);
  std::vector<BigInteger> packed_messages;
  packed_messages.reserve(num_packed_coeffs);

  // Each packed_message is a base-packing_base value with messages as digits,
  // i.e. packed_message = sum(messages[k * j + i] * B^i, i=0..k-1),
  // for k = packing_dimension.
  BigInteger packed_message = 0;
  BigInteger multiplier = 1;
  int packing_idx = 0;
  for (int i = 0; i < messages.size(); ++i) {
    Integer x = messages[i];
    packed_message += x * multiplier;
    multiplier *= packing_base;
    packing_idx++;

    // Packed enough messages on an integer
    if (packing_idx == packing_dimension || i == messages.size() - 1) {
      packed_messages.push_back(packed_message);
      packed_message = 0;
      multiplier = 1;
      packing_idx = 0;
    }
  }
  return packed_messages;
}

// Unpacks large integers from [0, `packing_base`^`packing_dimension`) into
// smaller integers, by decomposing them on base `packing_base`. If the number
// of messages is not a multiple of `packing_dimension`, we pad with zeros.
// Expects packing_base > 1, packing_dimension > 0, and
// packing_base^packing_dimension < std::numeric_limits<BigInteger>::max().
template <typename Integer, typename BigInteger>
std::vector<Integer> UnpackMessagesFlat(
    absl::Span<const BigInteger> packed_messages, Integer packing_base,
    int packing_dimension) {
  std::vector<Integer> unpacked_messages;
  unpacked_messages.reserve(packing_dimension * packed_messages.size());

  for (BigInteger packed_message : packed_messages) {
    // Decompose `packed_message` on base `packing_base`
    for (int k = 0; k < packing_dimension; ++k) {
      unpacked_messages.push_back(
          static_cast<Integer>(packed_message % packing_base));
      packed_message /= packing_base;
    }
  }
  return unpacked_messages;
}

// Packs integers in [0, `packing_base`) from `messages` into a smaller number
// of larger integers in [0, `packing_base`^`packing_dimension`). Groups the
// packed messages into vectors with at most `num_coeffs_per_packed_vector`
// elements. If the number of messages is not a multiple of `packing_dimension`,
// we pad with zeros. Expects packing_base > 1, packing_dimension > 0,
// num_coeffs_per_packed_vector > 0, and packing_base^packing_dimension <
// std::numeric_limits<BigInteger>::max().
template <typename Integer, typename BigInteger>
std::vector<std::vector<BigInteger>> PackMessages(
    absl::Span<const Integer> messages, Integer packing_base,
    int packing_dimension, int num_coeffs_per_packed_vector) {
  std::vector<BigInteger> packed_messages =
      PackMessagesFlat<Integer, BigInteger>(messages, packing_base,
                                            packing_dimension);
  int num_packed_coeffs = packed_messages.size();
  int num_packed_vectors = internal::DivideAndRoundUp<int>(
      num_packed_coeffs, num_coeffs_per_packed_vector);
  std::vector<std::vector<BigInteger>> all_packed_messages(num_packed_vectors);
  auto packed_messages_it = packed_messages.begin();
  for (int i = 0; i < num_packed_vectors;
       ++i, packed_messages_it += num_coeffs_per_packed_vector) {
    all_packed_messages[i].resize(num_coeffs_per_packed_vector, 0);
    std::copy_n(packed_messages_it,
                std::min(num_coeffs_per_packed_vector, num_packed_coeffs),
                all_packed_messages[i].begin());
    num_packed_coeffs -= num_coeffs_per_packed_vector;
  }
  return all_packed_messages;
}

// Unpacks large integers from [0, `packing_base`^`packing_dimension`) into
// smaller integers, by decomposing them on base `packing_base`. Expects
// packing_base > 1, packing_dimension > 0, packing_base^packing_dimension <
// std::numeric_limits<BigInteger>::max().
template <typename Integer, typename BigInteger>
std::vector<Integer> UnpackMessages(
    const std::vector<std::vector<BigInteger>>& all_packed_messages,
    Integer packing_base, int packing_dimension) {
  std::vector<BigInteger> packed_messages;
  for (const auto& packed_vector : all_packed_messages) {
    packed_messages.insert(packed_messages.end(), packed_vector.begin(),
                           packed_vector.end());
  }
  return UnpackMessagesFlat<Integer, BigInteger>(packed_messages, packing_base,
                                                 packing_dimension);
}

}  // namespace rlwe

#endif  // RLWE_RNS_MESSAGE_PACKING_H_
