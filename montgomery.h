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

// Defines types that are necessary for Ring-Learning with Errors (rlwe)
// encryption.

#ifndef RLWE_MONTGOMERY_H_
#define RLWE_MONTGOMERY_H_

#include <cmath>
#include <cstdint>
#include <tuple>
#include <vector>
#include "serialization.pb.h"
#include "utils.h"

namespace rlwe {

// The parameters necessary for a Montgomery integer.
struct MontgomeryIntParams {
  MontgomeryIntParams(uint64_t log_r, uint64_t modulus) {
    this->log_r = log_r;
    this->r = 1L << log_r;  // 1L to handle moduli > 31 bits.

    this->modulus = modulus;
    this->log_modulus = ceil(log2(modulus));
    this->r_mod_modulus = r % modulus;

    std::tuple<uint64_t, uint64_t> result = Inverses(log_r, modulus);
    this->inv_r = std::get<0>(result);
    this->inv_modulus = std::get<1>(result);

    // Select the barrett parameters. 1L to handle moduli > 31 bits.
    this->barrett_k = log_modulus + 1;
    this->barrett_numerator = (1L << barrett_k) / modulus;
  }

  // The binary log of the value of R to be used in Montgomery representation.
  // This class assumes R = 2^kLogR (i.e., a 1 followed by kLogR-1 zeros.)
  uint64_t log_r;

  // The concrete value of R.
  uint64_t r;

  // The modulus over which these modular operations are being performed.
  uint64_t modulus;

  // The number of bits in the modulus.
  unsigned int log_modulus;

  // The value of R taken modulo the modulus.
  uint64_t r_mod_modulus;

  // A value that, together with log_R and Modulus, satisfies the
  // formula:
  //     2^kLogR * inv_r - modulus * inv_modulus = 1
  // It is not literally the multiplicative inverse of Modulus.
  uint64_t inv_modulus;

  // The inverse of R needed to translate from Montgomery to normal rep.
  uint64_t inv_r;

  // The exponent k used in the denominator of the Barrett reduction. 2^k must
  // be greater than the modulus N. The greater the value of k, the larger the
  // range of values on which the Barrett reduction will work. k is set to 1
  // greater than the number of bits in N to ensure 2^k > N.
  uint64_t barrett_k;

  // The numerator used in the Barrett reduction.
  uint64_t barrett_numerator;

 private:
  // From Hacker's Delight.
  // Computes the Montgomery inverse coefficients for r and modulus using
  // the Extended Euclidean Algorithm.
  //
  // modulus must be odd.
  // Returns a tuple of (inv_r, inv_modulus) such that:
  //     r * inv_r - modulus * inv_modulus = 1
  std::tuple<int64_t, int64_t> Inverses(uint64_t log_r, uint64_t modulus) {
    // Invariants:
    //   1) sum = x * 2r - y * modulus.
    //   2) sum is always a power of 2.
    //   3) modulus is odd.
    //   4) y is always even.
    uint64_t log_sum = log_r + 1;
    uint64_t r = 1L << log_r;  // 1L to handle moduli > 31 bits.
    uint64_t x = 1;
    uint64_t y = 0;
    log_r -= 1;

    // Repeatedly divide both sides by 2 until sum == 1.
    while (log_sum > 0) {
      // Divide the left side by 2.
      log_sum -= 1;

      // Ensure that x is even.
      if (x % 2 != 0) {
        // If x is odd, make x even by adding modulus to x and changing the
        // value of y accordingly (y remains even).
        //
        //     sum = (x + modulus) * 2r - y * modulus - 2 * r * modulus
        //     sum = (x + modulus) * 2r - (y + 2r) * modulus
        //
        // We can then divide the new values of x and y by 2 safely.
        x += modulus;
        y += 2 * r;
      }

      // Divide the right side by 2.
      x >>= 1;
      y >>= 1;
    }

    return std::make_tuple(x << 1, y);
  }
};

// Stores an integer in Montgomery representation. The goal of this
// class is to provide a static type that differentiates between integers
// in Montgomery representation and standard integers. Once a Montgomery
// integer is created, it has the * and + operators necessary to be treated
// as another integer.
class MontgomeryInt {
 public:
  // Expose the underlying integer type.
  using Int = uint64_t;

  // Expose the parameter type.
  using Params = MontgomeryIntParams;

  // Static factory that converts a non-Montgomery representation integer
  // into a Montgomery representation integer. Does not take ownership of
  // params.
  // i.e., import "a".
  static MontgomeryInt ImportInt(const Params* params, uint64_t n) {
    return MontgomeryInt(
        params,
        ((n % params->modulus) * params->r_mod_modulus) % params->modulus);
  }

  // Import a random integer. Does not take ownership of params.
  static MontgomeryInt ImportRandom(const Params* params,
                                    utils::SecurePrng* rand) {
    return MontgomeryInt(params, rand->Rand64() % params->modulus);
  }

  // No default constructor.
  MontgomeryInt() = delete;

  // Default copy constructor.
  MontgomeryInt(const MontgomeryInt& that) = default;
  MontgomeryInt& operator=(const MontgomeryInt& that) = default;

  // Convert a Montgomery representation integer back to a normal integer.
  // i.e., export "a".
  uint64_t ExportInt() const {
    return (n_ * params_->inv_r) % params_->modulus;
  }

  SerializedModularInt Serialize() const {
    unsigned int byte_length = (params_->log_modulus + 7) / 8;
    std::string output(byte_length, 0);

    uint64_t remaining_bits = n_;
    for (unsigned int i = 0; i < byte_length; i++) {
      // Extract the low-order 8 bits of remaining_bits.
      output[i] = remaining_bits & 0xFF;

      // Shift right.
      remaining_bits >>= 8;
    }

    SerializedModularInt s;
    s.set_payload(output);

    return s;
  }

  // Does not take ownership of params.
  static MontgomeryInt Deserialize(const Params* params,
                                   const SerializedModularInt& serialized) {
    MontgomeryInt output(params, 0);
    std::string payload = serialized.payload();

    for (int i = payload.length() - 1; i >= 0; i--) {
      output.n_ <<= 8;
      output.n_ |= payload[i] & 0xFF;
    }

    return output;
  }

  // Perform a multiply followed by a modular reduction. Produces an output
  // in the range [0, N).
  // Taken from Hacker's Delight chapter on Montgomery multiplication.
  MontgomeryInt operator*(MontgomeryInt that) {
    // The product of the two numbers (a and b). Will equal a * R * b * R.
    uint64_t t = n_ * that.n_;

    // A bit mask used for computing (mod R). Creates all 1's in the log_R-1
    // bits of the integer. Computing x & mask is the same as x % R.
    uint64_t mask = params_->r - 1;

    // Perform the reduction, creating a * b * R (mod N).
    uint64_t u = t + params_->modulus * ((t * params_->inv_modulus) & mask);
    u >>= params_->log_r;

    // The steps above produce an integer that is in the range [0, 2N).
    if (u > params_->modulus) {
      return MontgomeryInt(params_, u - params_->modulus);
    }
    return MontgomeryInt(params_, u);
  }

  // Montgomery addition. Performs a Barrett reduction.
  MontgomeryInt operator+(const MontgomeryInt& that) const {
    uint64_t sum = n_ + that.n_;
    uint64_t u = (sum * params_->barrett_numerator) >> params_->barrett_k;
    uint64_t result = sum - u * params_->modulus;

    // u is in the range [0, 2n).
    if (result >= params_->modulus) {
      return MontgomeryInt(params_, result - params_->modulus);
    }

    return MontgomeryInt(params_, result);
  }

  // Negation.
  MontgomeryInt operator-() const {
    return MontgomeryInt(params_, params_->modulus - n_);
  }
  MontgomeryInt operator-(const MontgomeryInt& that) const {
    return *this + (-that);
  }

  // Equality.
  bool operator==(const MontgomeryInt& that) const { return (n_ == that.n_); }
  bool operator!=(const MontgomeryInt& that) const { return !(*this == that); }

  // Modular exponentiation.
  MontgomeryInt operator^(Int exponent) const {
    MontgomeryInt result = ImportInt(params_, 1);
    MontgomeryInt base = *this;

    // Uses the bits of the exponent to gradually compute the result.
    // When bit k of the exponent is 1, the result is multiplied by
    // base^{2^k}.
    while (exponent > 0) {
      // If the current bit (bit k) is 1, multiply base^{2^k} into the result.
      if (exponent % 2 == 1) {
        result = result * base;
      }

      // Update base from base^{2^k} to base^{2^{k+1}}.
      base = base * base;
      exponent >>= 1;
    }

    return result;
  }

  MontgomeryInt MultiplicativeInverse() const {
    return (*this) ^ (params_->modulus - 2);
  }

 private:
  explicit MontgomeryInt(const Params* params, uint64_t n)
      : n_(n), params_(params) {}
  uint64_t n_;
  const Params* params_;
};

}  // namespace rlwe

#endif  // RLWE_MONTGOMERY_H_
