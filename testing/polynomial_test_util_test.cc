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

#include "testing/polynomial_test_util.h"
#include <random>
#include <vector>
#include <gtest/gtest.h>
#include "constants.h"

namespace {

using Polynomial = rlwe::testing::PolynomialTestUtil<rlwe::kNewhopeModulus>;

TEST(PolynomialTestUtilTest, Len) {
  unsigned int seed = 0;

  for (int i = 1; i < 20; i++) {
    std::vector<uint64_t> v(i, 0);
    for (int j = 0; j < i; j++) {
      v[j] = rand_r(&seed) % 1000 + 1;
    }

    Polynomial p(v);
    EXPECT_EQ(i, p.Len());
  }
}

TEST(PolynomialTestUtilTest, Equality) {
  unsigned int seed = 0;

  for (int i = 1; i < 20; i++) {
    std::vector<uint64_t> v(i), w(i);
    for (int j = 0; j < i; j++) {
      v[j] = rand_r(&seed) % 1000 + 1;
      w[j] = rand_r(&seed) % 1000 + 1;
    }

    // Ensure v is really different than w.
    int k = rand_r(&seed) % i;
    v[k] = w[k] + 1;

    Polynomial p(v), q(v), r(w);

    // Ensure the equality relations hold.
    EXPECT_EQ(p, q);
    EXPECT_NE(p, r);
    EXPECT_NE(q, r);

    // Try flipping one value at a time.
    for (int j = 0; j < i; j++) {
      v[j]++;
      Polynomial s(v);
      EXPECT_NE(p, s);
      v[j]--;
    }
  }
}

TEST(PolynomialTestUtilTest, Degree) {
  unsigned int seed = 0;

  for (int i = 1; i < 20; i++) {
    // Slowly increase the degree of the polynomial.
    std::vector<uint64_t> v(i, 0);
    EXPECT_EQ(0, Polynomial(v).Degree());
    for (int j = 0; j < i; j++) {
      v[j] = rand_r(&seed) % 1000 + 1;
      EXPECT_EQ(j, Polynomial(v).Degree());
    }
  }
}

TEST(PolynomialTestUtilTest, PointwiseAdd) {
  unsigned int seed = 0;

  for (int i = 1; i < 20; i++) {
    std::vector<uint64_t> v(i), w(i), v_plus_w(i);

    for (int j = 0; j < i; j++) {
      v[j] = rand_r(&seed) % 1000 + 1;
      w[j] = rand_r(&seed) % 1000 + 1;
      v_plus_w[j] = (v[j] + w[j]) % rlwe::kNewhopeModulus;
    }

    Polynomial p(v), q(w), p_plus_q(v_plus_w);

    Polynomial r = p + q;
    Polynomial s = q + p;

    EXPECT_EQ(p_plus_q, r);
    EXPECT_EQ(p_plus_q, s);
    EXPECT_EQ(r, s);
  }
}

TEST(PolynomialTestUtilTest, PointwiseAddWithOverflow) {
  unsigned int seed = 0;

  for (int i = 1; i < 20; i++) {
    std::vector<uint64_t> v(i), w(i), v_plus_w(i);

    for (int j = 0; j < i; j++) {
      v[j] = rlwe::kNewhopeModulus - (rand_r(&seed) % 1000 + 1);
      w[j] = rlwe::kNewhopeModulus - (rand_r(&seed) % 1000 + 1);
      v_plus_w[j] = (v[j] + w[j]) % rlwe::kNewhopeModulus;
    }

    Polynomial p(v), q(w), p_plus_q(v_plus_w);

    Polynomial r = p + q;
    Polynomial s = q + p;

    EXPECT_EQ(p_plus_q, r);
    EXPECT_EQ(p_plus_q, s);
    EXPECT_EQ(r, s);
  }
}

TEST(PolynomialTestUtilTest, ScalarMultiply) {
  unsigned int seed = 0;

  for (int i = 1; i < 20; i++) {
    std::vector<uint64_t> v(i), v_times_k(i);

    int k = rand_r(&seed) % 1000 + 1;

    for (int j = 0; j < i; j++) {
      v[j] = rand_r(&seed) % 1000 + 1;
      v_times_k[j] = (v[j] * k) % rlwe::kNewhopeModulus;
    }

    Polynomial p(v), p_times_k(v_times_k);
    Polynomial r = p * k;

    EXPECT_EQ(p_times_k, r);
  }
}

TEST(PolynomialTestUtilTest, ScalarMultiplyOverflow) {
  unsigned int seed = 0;

  for (int i = 1; i < 20; i++) {
    std::vector<uint64_t> v(i), v_times_k(i);

    int k = rand_r(&seed) % 1000 + 2;

    for (int j = 0; j < i; j++) {
      v[j] = rlwe::kNewhopeModulus - (rand_r(&seed) % 1000 + 1);
      v_times_k[j] = (v[j] * k) % rlwe::kNewhopeModulus;
    }

    Polynomial p(v), p_times_k(v_times_k);
    Polynomial r = p * k;

    EXPECT_EQ(p_times_k, r);
  }
}

TEST(PolynomialTestUtilTest, Multiply) {
  unsigned int seed = 0;

  for (int i = 2; i < 20; i++) {
    std::vector<uint64_t> v(i), w(i);

    for (int j = 0; j < i; j++) {
      v[j] = rand_r(&seed) % 1000 + 1;
      w[j] = rand_r(&seed) % 1000 + 1;
    }

    Polynomial p(v), q(w);
    Polynomial r = p * q;
    Polynomial s = q * p;

    EXPECT_GE(r.Degree(), 2 * (i - 1));
    EXPECT_GE(r.Len(), 2 * i - 1);
    EXPECT_GE(s.Degree(), 2 * (i - 1));
    EXPECT_GE(s.Len(), 2 * i - 1);

    std::vector<uint64_t> result(2 * i, 0);

    for (int j = 0; j < i; j++) {
      for (int k = 0; k < i; k++) {
        uint64_t product = (v[j] * w[k]) % rlwe::kNewhopeModulus;

        // Handle modular arithmetic to avoid getting negative numbers.
        result[j + k] += product;
        result[j + k] %= rlwe::kNewhopeModulus;
      }
    }
    EXPECT_EQ(Polynomial(result), r);
    EXPECT_EQ(Polynomial(result), s);
    EXPECT_EQ(r, s);
  }
}

TEST(PolynomialTestUtilTest, MultiplyOverflow) {
  unsigned int seed = 0;

  for (int i = 2; i < 20; i++) {
    std::vector<uint64_t> v(i), w(i);

    for (int j = 0; j < i; j++) {
      v[j] = rlwe::kNewhopeModulus - (rand_r(&seed) % 1000 + 1);
      w[j] = rlwe::kNewhopeModulus - (rand_r(&seed) % 1000 + 1);
    }

    Polynomial p(v), q(w);
    Polynomial r = p * q;
    Polynomial s = q * p;

    EXPECT_GE(r.Degree(), 2 * (i - 1));
    EXPECT_GE(r.Len(), 2 * i - 1);
    EXPECT_GE(s.Degree(), 2 * (i - 1));
    EXPECT_GE(s.Len(), 2 * i - 1);

    std::vector<uint64_t> result(2 * i, 0);

    for (int j = 0; j < i; j++) {
      for (int k = 0; k < i; k++) {
        uint64_t product = (v[j] * w[k]) % rlwe::kNewhopeModulus;

        // Handle modular arithmetic to avoid getting negative numbers.
        result[j + k] += product;
        result[j + k] %= rlwe::kNewhopeModulus;
      }
    }
    EXPECT_EQ(Polynomial(result), r);
    EXPECT_EQ(Polynomial(result), s);
    EXPECT_EQ(r, s);
  }
}

TEST(PolynomialTestUtilTest, RemainderModXNPlus1) {
  unsigned int seed = 0;

  for (int i = 2; i < 10; i++) {
    // Create a random polynomial.
    std::vector<uint64_t> v(i);

    for (int j = 0; j < i; j++) {
      v[j] = rand_r(&seed) % 1000 + 1;
    }

    // Test for how underflow is handled by setting a very big first coeff.
    Polynomial p(v);
    v[i - 1] = rlwe::kNewhopeModulus - (rand_r(&seed) % 1000 + 1);
    Polynomial q(v);

    // Pick a value of N to generate our modulus X^N + 1.
    for (int N = 1; N < i * 2; N++) {
      std::vector<uint64_t> w(N + 1);
      w[0] = 1;
      w[N] = 1;
      Polynomial modulus(w);

      Polynomial divisible_p = p * modulus;
      Polynomial divisible_q = q * modulus;
      Polynomial zero(std::vector<uint64_t>(i, 0));

      // Each divisible should currently be cleanly divisible by the modulus.
      EXPECT_EQ(zero, divisible_p.RemainderModXNPlus1(N));
      EXPECT_EQ(zero, divisible_q.RemainderModXNPlus1(N));

      // Pick a remainder smaller than N.
      for (int k = 1; k < N; k++) {
        std::vector<uint64_t> x(k);

        for (int l = 0; l < k; l++) {
          x[l] = rand_r(&seed) % 1000 + 1;
        }

        Polynomial r(x);

        Polynomial indivisible_p = divisible_p + r;
        Polynomial indivisible_q = divisible_q + r;

        EXPECT_EQ(r, indivisible_p.RemainderModXNPlus1(N));
        EXPECT_EQ(r, indivisible_q.RemainderModXNPlus1(N));
      }
    }
  }
}

}  // namespace
