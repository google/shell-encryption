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

#include "ntt_parameters.h"
#include <cstdint>
#include <vector>
#include <gtest/gtest.h>
#include "constants.h"
#include "montgomery.h"

using uint_m = rlwe::MontgomeryInt;

namespace {

static const uint64_t kQ = rlwe::kNewhopeModulus;  // The modulus.
static const uint64_t kLogN = rlwe::kNewhopeLogDegreeBound;
static const uint64_t kN = rlwe::kNewhopeDegreeBound;
static const rlwe::MontgomeryIntParams params14(rlwe::kNewhopeLogR,
                                                rlwe::kNewhopeModulus);

TEST(NttParametersTest, PrimitiveNthRootOfUnity) {
  unsigned int log_ns[] = {2u, 4u, 6u, 8u, 11u};
  unsigned int len = 5;

  for (unsigned int i = 0; i < len; i++) {
    uint_m w =
        rlwe::internal::PrimitiveNthRootOfUnity<uint_m>(&params14, log_ns[i]);
    unsigned int n = 1 << log_ns[i];

    // Ensure it is really a n-th root of unity.
    EXPECT_EQ(w ^ n, uint_m::ImportInt(&params14, 1)) << "hi";

    // Ensure it is really a primitive n-th root of unity.
    EXPECT_NE(w ^ (n / 2), uint_m::ImportInt(&params14, 1)) << "hi2";
  }
}

TEST(NttParametersTest, NttPsis) {
  // The values of psi should be the powers of the primitive 2n-th root of
  // unity.
  // Obtain the psis.
  std::vector<uint_m> psis = rlwe::NttPsis<uint_m>(&params14, kLogN);

  // Verify that that the 0th entry is 1.
  uint_m one = uint_m::ImportInt(&params14, 1);
  EXPECT_EQ(one, psis[0]);

  // Verify that the 1th entry is a primitive 2n-th root of unity.
  EXPECT_EQ(one, psis[1] ^ (kN << 1));
  EXPECT_NE(one, psis[1] ^ kN);

  // Verify that each subsequent entry is the appropriate power of the 1th
  // entry.
  for (unsigned int i = 2; i < kN; i++) {
    EXPECT_EQ(psis[i], psis[1] ^ i);
  }
}

TEST(NttParametersTest, NttPsisInv) {
  // Obtain the psis.
  std::vector<uint_m> psis = rlwe::NttPsis<uint_m>(&params14, kLogN);
  std::vector<uint_m> psis_inv = rlwe::NttPsisInv<uint_m>(&params14, kLogN);
  uint_m n = uint_m::ImportInt(&params14, 1 << kLogN);

  for (unsigned int i = 0; i < 1 << kLogN; i++) {
    EXPECT_EQ(1, (psis[i] * psis_inv[i] * n).ExportInt());
  }
}

TEST(NttParametersTest, NttOmegas) {
  // Obtain the omegas.
  std::vector<std::vector<uint_m>> omegas =
      rlwe::NttOmegas<uint_m>(&params14, kLogN);

  for (unsigned int i = 0; i < kLogN; i++) {
    // Ensure that the 0th entry of row i is 1.
    uint_m one = uint_m::ImportInt(&params14, 1);
    EXPECT_EQ(one, omegas[i][0]);

    // Ensure that the 1th entry of row i is the 2^(i+1)th primitive root
    // of unity.
    uint_m w = omegas[i][1];
    EXPECT_EQ(one, w ^ (1 << (i + 1)));
    EXPECT_NE(one, w ^ (1 << i));

    // Ensure the remaining entries of each row are powers of w.
    for (unsigned int j = 2; j < kN; j++) {
      EXPECT_EQ(w ^ j, omegas[i][j]);
    }
  }
}

TEST(NttParametersTest, NttOmegasInv) {
  // Obtain the omegas.
  std::vector<std::vector<uint_m>> omegas =
      rlwe::NttOmegas<uint_m>(&params14, kLogN);
  std::vector<std::vector<uint_m>> omegas_inv =
      rlwe::NttOmegasInv<uint_m>(&params14, kLogN);

  for (unsigned int i = 0; i < kLogN; i++) {
    for (unsigned int j = 0; j < kN; j++) {
      EXPECT_EQ(1, (omegas[i][j] * omegas_inv[i][j]).ExportInt());
    }
  }
}

TEST(NttParametersTest, Bitrev) {
  for (unsigned int log_N = 2; log_N < 11; log_N++) {
    unsigned int N = 1 << log_N;
    std::vector<unsigned int> mapping = rlwe::BitrevArray(log_N);

    // Visit each entry of the array.
    for (unsigned int i = 0; i < N; i++) {
      for (unsigned int j = 0; j < log_N; j++) {
        // Ensure bit j of i is equal to bit (log_N - j) of mapping[i].
        uint64_t mask1 = 1 << j;
        uint64_t mask2 = 1 << (log_N - j - 1);
        EXPECT_EQ((i & mask1) == 0, (mapping[i] & mask2) == 0);
      }
    }
  }
}

}  // namespace
