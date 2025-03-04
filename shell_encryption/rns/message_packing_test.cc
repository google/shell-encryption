// Copyright 2024 Google LLC
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include "shell_encryption/rns/message_packing.h"

#include <cstdint>
#include <vector>

#include <gtest/gtest.h>
#include "absl/types/span.h"
#include "shell_encryption/int256.h"

namespace rlwe {
namespace {

using Integer = uint64_t;
using BigInteger = uint256;

TEST(PackingTest, PackU64Integers) {
  int num_packing = 3;
  int num_coeffs = 4;
  Integer input_domain = 10;

  std::vector<Integer> input_values;
  std::vector<std::vector<BigInteger>> packed_messages;
  std::vector<std::vector<BigInteger>> expected_packed_messages;

  input_values = {1, 2, 3, 4, 5, 6, 7, 8};
  packed_messages = PackMessages<Integer, BigInteger>(
      input_values, input_domain, num_packing, num_coeffs);
  expected_packed_messages = {{123, 456, 780}};
  EXPECT_EQ(packed_messages, expected_packed_messages);

  input_values = {1, 2, 3, 4, 5, 6, 7, 8, 9};
  packed_messages = PackMessages<Integer, BigInteger>(
      input_values, input_domain, num_packing, num_coeffs);
  expected_packed_messages = {{123, 456, 789}};
  EXPECT_EQ(packed_messages, expected_packed_messages);

  input_values = {};
  packed_messages = PackMessages<Integer, BigInteger>(
      input_values, input_domain, num_packing, num_coeffs);
  expected_packed_messages = {};
  EXPECT_EQ(packed_messages, expected_packed_messages);

  input_values = {0};
  packed_messages = PackMessages<Integer, BigInteger>(
      input_values, input_domain, num_packing, num_coeffs);
  expected_packed_messages = {{0}};
  EXPECT_EQ(packed_messages, expected_packed_messages);

  num_packing = 1;
  num_coeffs = 1;
  input_values = {1, 2, 3};
  packed_messages = PackMessages<Integer, BigInteger>(
      input_values, input_domain, num_packing, num_coeffs);
  expected_packed_messages = {{1}, {2}, {3}};
  EXPECT_EQ(packed_messages, expected_packed_messages);

  num_packing = 1;
  num_coeffs = 2;
  input_values = {1, 2, 3};
  packed_messages = PackMessages<Integer, BigInteger>(
      input_values, input_domain, num_packing, num_coeffs);
  expected_packed_messages = {{1, 2}, {3}};
  EXPECT_EQ(packed_messages, expected_packed_messages);
}

TEST(PackingTest, UnpackU64Integers) {
  int num_packing = 3;
  Integer input_domain = 10;

  std::vector<std::vector<BigInteger>> packed_messages;
  std::vector<Integer> unpacked_messages;
  std::vector<Integer> expected_unpacked_messages;

  packed_messages = {{401, 321}, {450}};
  unpacked_messages =
      UnpackMessages(packed_messages, input_domain, num_packing);
  expected_unpacked_messages = {4, 0, 1, 3, 2, 1, 4, 5, 0};
  EXPECT_EQ(unpacked_messages, expected_unpacked_messages);

  num_packing = 2;
  input_domain = 2;

  packed_messages = {{1, 2}, {0, 3}};
  unpacked_messages =
      UnpackMessages(packed_messages, input_domain, num_packing);
  expected_unpacked_messages = {0, 1, 1, 0, 0, 0, 1, 1};
  EXPECT_EQ(unpacked_messages, expected_unpacked_messages);
}

TEST(PackingTest, PackUnpackU64Integers) {
  int num_packing = 3;
  int num_coeffs = 4;
  Integer input_domain = 10;

  std::vector<Integer> input_values;
  std::vector<std::vector<BigInteger>> packed_messages;
  input_values = {1, 2, 3, 4, 5, 6, 7, 8, 9};
  packed_messages = PackMessages<Integer, BigInteger>(
      input_values, input_domain, num_packing, num_coeffs);
  std::vector<Integer> unpacked_values =
      UnpackMessages(packed_messages, input_domain, num_packing);
  EXPECT_EQ(unpacked_values, input_values);

  input_values = {1, 2, 3, 4, 5, 6, 7, 8};
  unpacked_values = UnpackMessages(packed_messages, input_domain, num_packing);
  EXPECT_NE(unpacked_values, input_values);
}

TEST(PackingTest, PackUnpackU256Integers) {
  constexpr int num_packing = 8;
  constexpr int num_coeffs = 50;
  constexpr int num_messages = 30;
  constexpr uint64_t input_domain = 8;

  std::vector<Integer> input_vec = {2, 2, 4, 7, 4, 0, 0, 1, 4, 3,
                                    5, 3, 4, 5, 4, 4, 4, 3, 3, 2,
                                    2, 0, 1, 3, 1, 5, 2, 2, 2, 0};

  std::vector<std::vector<BigInteger>> packed_messages =
      rlwe::PackMessages<Integer, BigInteger>(input_vec, input_domain,
                                              num_packing, num_coeffs);

  std::vector<Integer> unpacked_messages =
      rlwe::UnpackMessages(packed_messages, input_domain, num_packing);

  EXPECT_EQ(absl::MakeSpan(unpacked_messages).subspan(0, num_messages),
            absl::MakeSpan(input_vec).subspan(0, num_messages));
}

}  // namespace
}  // namespace rlwe
