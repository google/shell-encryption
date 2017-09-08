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

#include "ntt_polynomial.h"
#include <cmath>
#include <random>
#include <vector>
#include <gtest/gtest.h>
#include "constants.h"
#include "montgomery.h"
#include "ntt_parameters.h"
#include "serialization.pb.h"
#include "testing/polynomial_test_util.h"

namespace {

// Useful typedefs.
using uint_m = rlwe::MontgomeryInt;
rlwe::MontgomeryIntParams params14(rlwe::kNewhopeLogR, rlwe::kNewhopeModulus);
using PolynomialTestUtil =
    rlwe::testing::PolynomialTestUtil<rlwe::kNewhopeModulus>;
using NttPolynomial = rlwe::NttPolynomial<uint_m>;

unsigned int seed = 0;

// Test fixture to take care of messy setup.
class NttPolynomialTest : public ::testing::Test {
 protected:
  static std::vector<uint_m> ConvertToMontgomery(
      std::vector<uint_m::Int> coeffs) {
    std::vector<uint_m> output(coeffs.size(), uint_m::ImportInt(&params14, 0));
    for (unsigned int i = 0; i < output.size(); i++) {
      output[i] = uint_m::ImportInt(&params14, coeffs[i]);
    }
    return output;
  }

  static std::vector<uint_m::Int> ConvertFromMontgomery(
      std::vector<uint_m> coeffs) {
    std::vector<uint_m::Int> output(coeffs.size());
    for (unsigned int i = 0; i < coeffs.size(); i++) {
      output[i] = coeffs[i].ExportInt();
    }

    return output;
  }

  void SetUp() override { srand(0); }

  void SetParams(int n, int log_n) {
    std::vector<uint_m::Int> p_coeffs(n);
    std::vector<uint_m::Int> q_coeffs(n);

    // Create some random polynomials. Ensure that they are different.
    for (int j = 0; j < n; j++) {
      p_coeffs[j] = rand_r(&seed) % rlwe::kNewhopeModulus;
      q_coeffs[j] = rand_r(&seed) % rlwe::kNewhopeModulus;
    }

    // Ensure the polynomials are different.
    uint_m::Int rand_index = rand_r(&seed) % n;
    p_coeffs[rand_index] = (q_coeffs[rand_index] + 1) % rlwe::kNewhopeModulus;

    p_.reset(new PolynomialTestUtil(p_coeffs));
    q_.reset(new PolynomialTestUtil(q_coeffs));

    // Acquire all of the NTT parameters.
    params_ = rlwe::InitializeNttParameters<uint_m>(&params14, log_n);

    // Put p and q in the NTT domain.
    ntt_p_ =
        NttPolynomial::ConvertToNtt(ConvertToMontgomery(p_coeffs), params_);
    ntt_q_ =
        NttPolynomial::ConvertToNtt(ConvertToMontgomery(q_coeffs), params_);
  }

  rlwe::NttParameters<uint_m> params_;
  std::unique_ptr<PolynomialTestUtil> p_;
  std::unique_ptr<PolynomialTestUtil> q_;
  NttPolynomial ntt_p_;
  NttPolynomial ntt_q_;
};

// Ensure that a default NTT polynomial is invalid.
TEST_F(NttPolynomialTest, DefaultIsInvalid) {
  EXPECT_FALSE(NttPolynomial().IsValid());
}

// Ensure that a polynomial converted to NTT form can be converted back.
TEST_F(NttPolynomialTest, Symmetry) {
  for (int i = 2; i < 11; i++) {
    for (int k = 0; k < 100; k++) {
      SetParams(1 << i, i);
      EXPECT_TRUE(ntt_p_.IsValid());
      PolynomialTestUtil p_prime(
          ConvertFromMontgomery(ntt_p_.InverseNtt(params_)));
      EXPECT_EQ(*p_, p_prime);
    }
  }
}

// Ensure that equality holds properly.
TEST_F(NttPolynomialTest, Equality) {
  for (int i = 2; i < 11; i++) {
    for (int k = 0; k < 100; k++) {
      SetParams(1 << i, i);
      NttPolynomial ntt_p_cpy = ntt_p_;
      NttPolynomial ntt_q_cpy = ntt_q_;

      EXPECT_TRUE(ntt_p_ == ntt_p_cpy);
      EXPECT_TRUE(ntt_q_ == ntt_q_cpy);
      EXPECT_FALSE(ntt_p_ != ntt_p_cpy);
      EXPECT_FALSE(ntt_q_ != ntt_q_cpy);

      EXPECT_TRUE(ntt_p_ != ntt_q_);
      EXPECT_TRUE(ntt_q_ != ntt_p_);
      EXPECT_FALSE(ntt_p_ == ntt_q_);
      EXPECT_FALSE(ntt_q_ == ntt_p_);
    }
  }
}

// Ensure that a polynomial whose size is not a power of two gets rejected.
TEST_F(NttPolynomialTest, NotPowerOfTwo) {
  for (int i = 2; i < 11; i++) {
    // j is any value that isn't a power of 2.
    for (int j = 1 + (1 << (i - 1)); j < (1 << i); j++) {
      SetParams(j, i);
      EXPECT_FALSE(ntt_p_.IsValid());
    }
  }
}

// Ensure that adding or multiplying two polynomials of different lengths gets
// rejected.
TEST_F(NttPolynomialTest, BinopOfDifferentLengths) {
  for (int i = 2; i < 11; i++) {
    for (int j = 2; j < 11; j++) {
      if (j == i) {
        continue;
      }

      int bigger = std::max(i, j);
      SetParams(1 << bigger, bigger);

      std::vector<uint_m::Int> x(1 << i);
      std::vector<uint_m::Int> y(1 << j);

      params_.bitrevs = rlwe::BitrevArray(i);
      auto ntt_x = NttPolynomial::ConvertToNtt(ConvertToMontgomery(x), params_);
      params_.bitrevs = rlwe::BitrevArray(j);
      auto ntt_y = NttPolynomial::ConvertToNtt(ConvertToMontgomery(y), params_);

      EXPECT_TRUE(ntt_x.IsValid());
      EXPECT_TRUE(ntt_y.IsValid());

      EXPECT_FALSE((ntt_x * ntt_y).IsValid());
      EXPECT_FALSE((ntt_y * ntt_x).IsValid());
      EXPECT_FALSE((ntt_x + ntt_y).IsValid());
      EXPECT_FALSE((ntt_y + ntt_x).IsValid());
    }
  }
}

// Test that the convolution property holds. Let p, q be polynomials.
// Polynomial multiplication of p and q =
// NTT_INV(the coordinate-wise product of NTT(p) and NTT(q))
TEST_F(NttPolynomialTest, Multiply) {
  for (int i = 2; i < 11; i++) {
    for (int k = 0; k < 100; k++) {
      SetParams(1 << i, i);

      EXPECT_TRUE(ntt_p_.IsValid());
      EXPECT_TRUE(ntt_q_.IsValid());

      NttPolynomial ntt_res1 = ntt_p_ * ntt_q_;
      NttPolynomial ntt_res2 = ntt_q_ * ntt_p_;

      EXPECT_TRUE(ntt_res1.IsValid());
      EXPECT_TRUE(ntt_res2.IsValid());

      PolynomialTestUtil res1(
          ConvertFromMontgomery(ntt_res1.InverseNtt(params_)));
      PolynomialTestUtil res2(
          ConvertFromMontgomery(ntt_res2.InverseNtt(params_)));

      PolynomialTestUtil expected = ((*p_) * (*q_)).RemainderModXNPlus1(1 << i);

      EXPECT_EQ(res1, expected);
      EXPECT_EQ(res2, expected);
      EXPECT_EQ(res1, res2);
    }
  }
}

// Test scalar multiplication.
TEST_F(NttPolynomialTest, ScalarMultiply) {
  uint_m scalar =
      uint_m::ImportInt(&params14, rand_r(&seed) % rlwe::kNewhopeModulus);
  for (int i = 2; i < 11; i++) {
    for (int k = 0; k < 100; k++) {
      SetParams(1 << i, i);

      EXPECT_TRUE(ntt_p_.IsValid());

      NttPolynomial ntt_res = ntt_p_ * scalar;

      EXPECT_TRUE(ntt_res.IsValid());

      PolynomialTestUtil res(
          ConvertFromMontgomery(ntt_res.InverseNtt(params_)));

      PolynomialTestUtil expected = ((*p_) * scalar.ExportInt());

      EXPECT_EQ(res, expected);
    }
  }
}

// Test that p + (-p) = 0.
TEST_F(NttPolynomialTest, Negate) {
  for (int i = 2; i < 11; i++) {
    for (int k = 0; k < 100; k++) {
      SetParams(1 << i, i);

      // An NTT polynomial of all zeros.
      NttPolynomial zeros_ntt = NttPolynomial::ConvertToNtt(
          std::vector<uint_m>(1 << i, uint_m::ImportInt(&params14, 0)),
          params_);

      EXPECT_EQ(zeros_ntt, ntt_p_ + (-ntt_p_));
      EXPECT_EQ(zeros_ntt, -ntt_p_ + ntt_p_);
    }
  }
}

// Test that p + q = NTT_INV(NTT(p) + NTT(q)).
TEST_F(NttPolynomialTest, Add) {
  for (int i = 2; i < 11; i++) {
    for (int k = 0; k < 100; k++) {
      SetParams(1 << i, i);

      EXPECT_TRUE(ntt_p_.IsValid());
      EXPECT_TRUE(ntt_q_.IsValid());

      NttPolynomial ntt_res1 = ntt_p_ + ntt_q_;
      NttPolynomial ntt_res2 = ntt_q_ + ntt_p_;

      EXPECT_TRUE(ntt_res1.IsValid());
      EXPECT_TRUE(ntt_res2.IsValid());

      PolynomialTestUtil res1(
          ConvertFromMontgomery(ntt_res1.InverseNtt(params_)));
      PolynomialTestUtil res2(
          ConvertFromMontgomery(ntt_res2.InverseNtt(params_)));

      PolynomialTestUtil expected = (*p_) + (*q_);

      EXPECT_EQ(res1, expected);
      EXPECT_EQ(res2, expected);
      EXPECT_EQ(res1, res2);
    }
  }
}

TEST_F(NttPolynomialTest, Serialize) {
  for (int i = 2; i < 11; i++) {
    for (int j = 0; j < 200; j++) {
      SetParams(1 << i, i);
      EXPECT_EQ(ntt_p_,
                NttPolynomial::Deserialize(&params14, ntt_p_.Serialize()));
      EXPECT_EQ(ntt_q_,
                NttPolynomial::Deserialize(&params14, ntt_q_.Serialize()));
    }
  }
}

}  // namespace
