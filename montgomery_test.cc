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
#include "montgomery.h"
#include <cstdint>
#include <random>
#include <vector>
#include <gtest/gtest.h>
#include "constants.h"
#include "serialization.pb.h"

namespace {

using uint_m = rlwe::MontgomeryInt;

// Access New Hope parameters.
rlwe::MontgomeryIntParams params14(rlwe::kNewhopeLogR, rlwe::kNewhopeModulus);

// Montgomery parameters for a 10-bit modulus.
constexpr uint64_t kModulus10 = 997;
constexpr uint64_t kLogR10 = 10;
constexpr uint64_t kInvModulus10 = 531;
rlwe::MontgomeryIntParams params10(kLogR10, kModulus10);

// Montgomery parameters for a 30-bit modulus.
rlwe::MontgomeryIntParams params30(rlwe::kLogR30, rlwe::kModulus30);

// Verifies that the MontgomeryIntParams code computes the inverse modulus.
TEST(MontgomeryTest, ParamsInvModulus) {
  EXPECT_EQ(
      1, params14.r * params14.inv_r - params14.modulus * params14.inv_modulus);
  EXPECT_EQ(
      1, params10.r * params10.inv_r - params10.modulus * params10.inv_modulus);
  EXPECT_EQ(
      1, params30.r * params30.inv_r - params30.modulus * params30.inv_modulus);
}

// Verifies that numbers can be imported and exported properly.
TEST(MontgomeryTest, ImportExportInt) {
  for (uint64_t i = 0; i < kModulus10; i++) {
    uint_m m = uint_m::ImportInt(&params10, i);
    uint64_t after = m.ExportInt();
    uint64_t expected = i % kModulus10;
    EXPECT_EQ(expected, after);
  }
}

// Verifies that numbers can be added correctly.
TEST(MontgomeryTest, AddSub) {
  // Test over a selection of the possible input space.
  unsigned int seed = 0;

  for (int i = 0; i < 10000; i++) {
    uint64_t a = rand_r(&seed) % rlwe::kNewhopeModulus;
    uint64_t b = rand_r(&seed) % rlwe::kNewhopeModulus;
    uint_m ma = uint_m::ImportInt(&params14, a);
    uint_m mb = uint_m::ImportInt(&params14, b);
    uint_m mc = ma + mb;
    uint64_t c = mc.ExportInt();

    uint64_t expected = (a + b) % rlwe::kNewhopeModulus;
    EXPECT_EQ(expected, c);

    uint_m md = ma - mb;
    uint64_t d = md.ExportInt();

    uint64_t expected2 =
        (a + rlwe::kNewhopeModulus - b) % rlwe::kNewhopeModulus;
    EXPECT_EQ(expected2, d);
  }
}

// Verifies that equality functions properly.
TEST(MontgomeryTest, Equality) {
  // Test over a selection of the possible input space.
  unsigned int seed = 0;

  for (int i = 0; i < 10000; i++) {
    uint64_t a = rand_r(&seed) % rlwe::kNewhopeModulus;
    uint64_t b = rand_r(&seed) % rlwe::kNewhopeModulus;
    while (b == a) {
      b = rand_r(&seed) % rlwe::kNewhopeModulus;
    }

    uint_m ma1 = uint_m::ImportInt(&params14, a);
    uint_m ma2 = uint_m::ImportInt(&params14, a);
    uint_m mb1 = uint_m::ImportInt(&params14, b);
    uint_m mb2 = uint_m::ImportInt(&params14, b);

    EXPECT_TRUE(ma1 == ma2);
    EXPECT_TRUE(ma2 == ma1);
    EXPECT_FALSE(ma1 != ma2);
    EXPECT_FALSE(ma2 != ma1);

    EXPECT_TRUE(mb1 == mb2);
    EXPECT_TRUE(mb2 == mb1);
    EXPECT_FALSE(mb1 != mb2);
    EXPECT_FALSE(mb2 != mb1);

    EXPECT_TRUE(ma1 != mb1);
    EXPECT_TRUE(mb1 != ma1);
    EXPECT_FALSE(ma1 == mb1);
    EXPECT_FALSE(mb1 == ma1);
  }
}

// Verifies that numbers can be negated correctly.
TEST(MontgomeryTest, Negate) {
  for (uint64_t i = 0; i < 4 * rlwe::kNewhopeModulus; i++) {
    uint_m mi = uint_m::ImportInt(&params14, i);

    EXPECT_EQ(0, (mi + (-mi)).ExportInt());
  }
}

// Verifies that repeated addition works properly.
TEST(MontgomeryTest, AddRepeatedly) {
  // Test over a selection of the possible input space.
  unsigned int seed = 0;

  for (int i = 0; i < 10000; i++) {
    uint64_t sum = 0;
    uint64_t diff = 0;
    uint_m mont_sum = uint_m::ImportInt(&params14, 0);
    uint_m mont_diff = uint_m::ImportInt(&params14, 0);

    for (int j = 0; j < 1000; j++) {
      uint64_t a = rand_r(&seed) % rlwe::kNewhopeModulus;

      sum = (sum + a) % rlwe::kNewhopeModulus;
      mont_sum = mont_sum + uint_m::ImportInt(&params14, a);

      diff = (diff + rlwe::kNewhopeModulus - a) % rlwe::kNewhopeModulus;
      mont_diff = mont_diff - uint_m::ImportInt(&params14, a);
    }

    EXPECT_EQ(sum, mont_sum.ExportInt());
    EXPECT_EQ(diff, mont_diff.ExportInt());
  }
}

// Verifies that numbers can be multiplied correctly.
TEST(MontgomeryTest, Multiply) {
  unsigned int seed = 0;

  // Test over many random values.
  for (int i = 0; i < 10000; i++) {
    uint64_t a = rand_r(&seed) % rlwe::kNewhopeModulus;
    uint64_t b = rand_r(&seed) % rlwe::kNewhopeModulus;
    uint_m ma = uint_m::ImportInt(&params14, a);
    uint_m mb = uint_m::ImportInt(&params14, b);
    uint_m mc = ma * mb;
    uint64_t c = mc.ExportInt();

    uint64_t expected = (a * b) % rlwe::kNewhopeModulus;
    EXPECT_EQ(expected, c);
  }
}

// Verifies that repeated addition works properly.
TEST(MontgomeryTest, MultiplyRepeatedly) {
  // Test over a selection of the possible input space.
  unsigned int seed = 0;

  for (int i = 0; i < 10000; i++) {
    uint64_t prod = 1;
    uint_m mont_prod = uint_m::ImportInt(&params14, 1);

    for (int j = 0; j < 1000; j++) {
      uint64_t a = rand_r(&seed) % rlwe::kNewhopeModulus;

      prod = (prod * a) % rlwe::kNewhopeModulus;
      mont_prod = mont_prod * uint_m::ImportInt(&params14, a);
    }

    EXPECT_EQ(prod, mont_prod.ExportInt());
  }
}

// Test the entire space for a small modulus.
TEST(MontgomeryTest, SmallModulus) {
  for (uint64_t a = 0; a < 2 * kModulus10; a++) {
    for (uint64_t b = 0; b < 2 * kModulus10; b++) {
      uint_m ma = uint_m::ImportInt(&params10, a);
      uint_m mb = uint_m::ImportInt(&params10, b);
      uint_m mc = ma + mb;

      // Equality.
      if (a % kModulus10 == b % kModulus10) {
        EXPECT_TRUE(ma == mb);
        EXPECT_FALSE(ma != mb);
      } else {
        EXPECT_TRUE(ma != mb);
        EXPECT_FALSE(ma == mb);
      }

      // Addition.
      EXPECT_EQ((a + b) % kModulus10, mc.ExportInt());
      EXPECT_EQ((a + b) % kModulus10, mc.ExportInt());

      // Negation.
      EXPECT_EQ((kModulus10 * 2 - a) % kModulus10, (-ma).ExportInt());
      EXPECT_EQ((kModulus10 * 2 - b) % kModulus10, (-mb).ExportInt());
      EXPECT_EQ((kModulus10 * 4 - a - b) % kModulus10, (-mc).ExportInt());

      // Subtraction.
      EXPECT_EQ((kModulus10 * 2 - a + b) % kModulus10, (mb - ma).ExportInt());
      EXPECT_EQ((kModulus10 * 2 - b + a) % kModulus10, (ma - mb).ExportInt());

      // Multiplication.
      EXPECT_EQ((a * b) % kModulus10, (ma * mb).ExportInt());
      EXPECT_EQ((a * b) % kModulus10, (mb * ma).ExportInt());
    }
  }
}

TEST(MontgomeryTest, ModExp) {
  for (uint64_t base = 0; base < 2 * kModulus10; base++) {
    uint64_t expected = 1;
    for (uint64_t exp = 0; exp < 1000; exp++) {
      uint_m base_m = uint_m::ImportInt(&params10, base);
      uint_m actual_m = base_m ^ exp;
      uint64_t actual = actual_m.ExportInt();
      ASSERT_EQ(actual, expected);

      expected *= base;
      expected %= kModulus10;
    }
  }
}

TEST(MontgomeryTest, Inverse) {
  for (uint64_t i = 1; i < kModulus10; i++) {
    uint_m i_m = uint_m::ImportInt(&params10, i);
    uint_m inv = i_m.MultiplicativeInverse();
    ASSERT_EQ((i_m * inv).ExportInt(), 1);
  }
}

TEST(MontgomeryTest, Serialization) {
  // Try all possible values in kModulus10.
  // Try these values in both kModulus14 and kModulus10.
  for (uint64_t i = 0; i < kModulus10; i++) {
    // Serialize and ensure the byte length is as expected.
    uint_m i10 = uint_m::ImportInt(&params10, i);
    uint_m i14 = uint_m::ImportInt(&params14, i);

    rlwe::SerializedModularInt serialized10 = i10.Serialize();
    rlwe::SerializedModularInt serialized14 = i14.Serialize();
    std::string payload10 = serialized10.payload();
    std::string payload14 = serialized14.payload();

    EXPECT_EQ(payload10.length(), (params10.log_modulus + 7) / 8);
    EXPECT_EQ(payload14.length(), (params14.log_modulus + 7) / 8);

    // Ensure that deserialization works properly.
    EXPECT_EQ(uint_m::Deserialize(&params10, serialized10), i10);
    EXPECT_EQ(uint_m::Deserialize(&params14, serialized14), i14);

    // Ensure that that any bit beyond bit the serialized bit length can be
    // wiped out without issue. That is, ensure that the bit size is accurate.
    payload10[1] &= (1 << 3) - 1;
    serialized10.set_payload(payload10);
    EXPECT_EQ(uint_m::Deserialize(&params10, serialized10), i10);

    payload14[1] &= (1 << 7) - 1;
    serialized14.set_payload(payload14);
    EXPECT_EQ(uint_m::Deserialize(&params14, serialized14), i14);
  }
}

}  // namespace
