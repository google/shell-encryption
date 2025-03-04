/*
 * Copyright 2024 Google LLC.
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

#include <vector>

namespace rlwe {

// Packs integers in [0, `packing_base`) from `messages` into a smaller number
// of larger integers in [0, `packing_base`^`num_packing`). Groups the packed
// messages into vectors with at most `num_coeffs` elements. If the number of
// messages is not a multiple of `num_packing`, we pad with zeros.
// Expects packing_base > 1, num_packing > 0, num_coeffs > 0,
// packing_base^num_packing < std::numeric_limits<BigInteger>::max().
template <typename Integer, typename BigInteger>
std::vector<std::vector<BigInteger>> PackMessages(
    const std::vector<Integer>& messages, Integer packing_base, int num_packing,
    int num_coeffs) {
  std::vector<std::vector<BigInteger>> all_packed_messages;
  int max_coeffs = messages.size() / num_packing + 1;
  int max_poly = max_coeffs / num_coeffs + 1;
  all_packed_messages.reserve(max_poly);

  std::vector<BigInteger> packed_messages;
  packed_messages.reserve(num_coeffs);
  Integer packed_message{0};
  int packing_idx = 0;
  for (int i = 0; i < messages.size(); ++i) {
    Integer x = messages[i];
    packed_message *= packing_base;
    packed_message += x;
    packing_idx++;

    // Pad with zeros if we reached the end of input
    if (i == messages.size() - 1) {
      while (packing_idx < num_packing) {
        packed_message *= packing_base;
        packing_idx++;
      }
    }

    // Packed enough messages on an integer
    if (packing_idx == num_packing) {
      packed_messages.push_back(packed_message);
      packed_message = 0;
      packing_idx = 0;
    }

    // Stored enough packed_messages in a vector, or reached the end of input
    if (packed_messages.size() == num_coeffs || i == messages.size() - 1) {
      all_packed_messages.push_back(std::move(packed_messages));
      packed_messages.clear();
    }
  }
  return all_packed_messages;
}

// Unpacks large integers from [0, `packing_base`^`num_packing`) into smaller
// integers, by decomposing them on base `packing_base`.
// Expects packing_base > 1, num_packing > 0, packing_base^num_packing <
// std::numeric_limits<BigInteger>::max().
template <typename Integer, typename BigInteger>
std::vector<Integer> UnpackMessages(
    const std::vector<std::vector<BigInteger>>& all_packed_messages,
    Integer packing_base, int num_packing) {
  std::vector<Integer> all_unpacked_messages;
  std::vector<Integer> unpacked_messages;
  unpacked_messages.reserve(num_packing);

  for (auto& packed_messages : all_packed_messages) {
    for (BigInteger packed_message : packed_messages) {
      // Decompose `packed_message` on base `packing_base`
      for (int k = 0; k < num_packing; ++k) {
        BigInteger quot = packed_message / packing_base;
        Integer rem = static_cast<Integer>(packed_message % packing_base);
        unpacked_messages.push_back(rem);
        packed_message = quot;
      }

      // Reverse order because we packed the first message in the most
      // significant position
      for (int k = num_packing - 1; k >= 0; --k) {
        all_unpacked_messages.push_back(unpacked_messages[k]);
      }
      unpacked_messages.clear();
    }
  }
  return all_unpacked_messages;
}

}  // namespace rlwe

#endif  // RLWE_RNS_MESSAGE_PACKING_H_
