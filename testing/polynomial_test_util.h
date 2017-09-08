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

#ifndef RLWE_TESTING_POLYNOMIAL_TEST_UTIL_H_
#define RLWE_TESTING_POLYNOMIAL_TEST_UTIL_H_

#include <cmath>
#include <cstdint>
#include <string>
#include <vector>

namespace rlwe {
namespace testing {

// A standard polynomial with coefficients modulo kModulus.
template <uint64_t kModulus>
class PolynomialTestUtil {
 public:
  // Copy constructor.
  PolynomialTestUtil(const PolynomialTestUtil& that) = default;

  // Constructor. The polynomial is initialized to the values of a vector.
  explicit PolynomialTestUtil(const std::vector<uint64_t>& coeffs) {
    coeffs_ = coeffs;
  }

  // Accessor for length.
  int Len() const { return coeffs_.size(); }

  // Compute the degree.
  int Degree() const {
    for (int i = Len() - 1; i >= 0; i--) {
      if (coeffs_[i] != 0) {
        return i;
      }
    }

    return 0;
  }

  // Equality.
  bool operator==(const PolynomialTestUtil& that) const {
    if (Degree() != that.Degree()) {
      return false;
    }

    for (int i = 0; i <= Degree(); i++) {
      if (coeffs_[i] != that.coeffs_[i]) {
        return false;
      }
    }

    return true;
  }

  bool operator!=(const PolynomialTestUtil& that) const {
    return !(*this == that);
  }

  // Polynomial addition.
  PolynomialTestUtil operator+(const PolynomialTestUtil& that) const {
    // Find the polynomial with the bigger degree.
    const PolynomialTestUtil* bigger = this;
    const PolynomialTestUtil* smaller = &that;
    if (Degree() < that.Degree()) {
      bigger = &that;
      smaller = this;
    }

    // Copy the polynomial with the bigger degree and add the smaller one to it.
    PolynomialTestUtil out(*bigger);
    for (int i = 0; i <= smaller->Degree(); i++) {
      out.coeffs_[i] = (out.coeffs_[i] + smaller->coeffs_[i]) % kModulus;
    }

    return out;
  }

  // Polynomial scalar multiplication.
  PolynomialTestUtil operator*(uint64_t c) const {
    PolynomialTestUtil out(*this);
    for (int i = 0; i < Len(); i++) {
      out.coeffs_[i] = (out.coeffs_[i] * c) % kModulus;
    }
    return out;
  }

  // Polynomial multiplication.
  PolynomialTestUtil operator*(const PolynomialTestUtil& that) const {
    std::vector<uint64_t> out(Degree() + that.Degree() + 1);

    for (int i = 0; i <= Degree(); i++) {
      for (int j = 0; j <= that.Degree(); j++) {
        out[i + j] += coeffs_[i] * that.coeffs_[j] % kModulus;
        out[i + j] %= kModulus;
      }
    }

    return PolynomialTestUtil(out);
  }

  // Remainder when divided by the polynomial (x^N + 1).
  PolynomialTestUtil RemainderModXNPlus1(int N) const {
    PolynomialTestUtil out(*this);

    while (out.Degree() >= N) {
      // Determine the difference in degrees.
      int diff = out.Degree() - N;

      // Determine the number of times X^N + 1 goes into the polynomial.
      int q = out.coeffs_[out.Degree()];

      // Zero out the highest degree term of out.
      out.coeffs_[out.Degree()] = 0;

      // Subtract out q * x^diff.
      if (out.coeffs_[diff] < q) {
        out.coeffs_[diff] += kModulus - q;
      } else {
        out.coeffs_[diff] -= q;
      }
    }

    return out;
  }

 private:
  std::vector<uint64_t> coeffs_;
};

}  // namespace testing
}  // namespace rlwe

#endif  // RLWE_TESTING_POLYNOMIAL_TEST_UTIL_H_
