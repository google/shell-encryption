/*
 * Copyright 2023 Google LLC.
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

#ifndef RLWE_RNS_ERROR_PARAMS_H_
#define RLWE_RNS_ERROR_PARAMS_H_

#include <cmath>
#include <numeric>
#include <utility>
#include <vector>

#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "absl/types/span.h"
#include "shell_encryption/rns/rns_modulus.h"

namespace rlwe {

// This class stores the error constants and contains functions to compute error
// bounds for various homomorphic operations.
//
// We consider RLWE homomorphic encryption schemes instantiated using the ring
// R_{Q*P} = Z[X]/(Q*P, X^N+1), where Q = q_0 * .. * q_{L-1} is the main modulus
// P = p_1 * .. * p_K is the auxiliary modulus, and N is a power of two.
// In addition, we assume plaintext modulus t is smaller than all q_i's, and we
// assume both the secret and fresh error are polynomials whose coefficients are
// independently sampled from a centered binomial distribution.
//
template <typename ModularInt>
class RnsErrorParams {
 public:
  // Creates an error parameter instance, where `sigma` is the standard
  // deviation of the error distribution.
  static absl::StatusOr<RnsErrorParams> Create(
      int log_n, absl::Span<const PrimeModulus<ModularInt>* const> main_moduli,
      absl::Span<const PrimeModulus<ModularInt>* const> aux_moduli,
      int log_plaintext_modulus, double sigma) {
    if (log_n <= 0) {
      return absl::InvalidArgumentError("`log_n` must be positive.");
    }
    if (main_moduli.empty()) {
      return absl::InvalidArgumentError("`main_moduli` must not be empty.");
    }
    // We allow the auxiliary moduli to be empty.
    if (log_plaintext_modulus <= 0) {
      return absl::InvalidArgumentError(
          "`log_plaintext_modulus` must be positive.");
    }
    if (sigma <= 0) {
      return absl::InvalidArgumentError("`sigma` must be positive.");
    }
    std::vector<int> log_main_prime_moduli;
    for (const auto& modulus : main_moduli) {
      log_main_prime_moduli.push_back(
          std::ceil(std::log2(static_cast<double>(modulus->Modulus()))));
    }
    std::vector<int> log_aux_prime_moduli;
    for (const auto& modulus : aux_moduli) {
      log_aux_prime_moduli.push_back(
          std::ceil(std::log2(static_cast<double>(modulus->Modulus()))));
    }
    return RnsErrorParams(log_n, std::move(log_main_prime_moduli),
                          std::move(log_aux_prime_moduli),
                          log_plaintext_modulus, sigma);
  }

  // A polynomial consisting of error terms is added to the ciphertext during
  // key switching. The noise of a ciphertext increases additively by the size
  // of the polynomial, which depends on the decomposition base of the gadget
  // used in the key-switching key and the number of ciphertext components
  // applied on.
  double BoundOnGadgetBasedKeySwitching(int num_components, int log_gadget_base,
                                        int num_digits) const {
    int gadget_base = 1 << log_gadget_base;
    int dimension = 1 << log_n_;
    return (8.0 / sqrt(3.0)) * ExportDoubleT() * num_digits * sigma_ *
           dimension * gadget_base * num_components;
  }

  // Accessors for frequently used error bound constants.
  double B_plaintext() const { return b_plaintext_; }
  double B_secretkey_encryption() const { return b_secretkey_encryption_; }
  double B_publickey_encryption() const { return b_publickey_encryption_; }
  double B_scale() const { return b_scale_; }

 private:
  // Constructor to set up the params.
  RnsErrorParams(int log_n, std::vector<int> log_main_prime_moduli,
                 std::vector<int> log_aux_prime_moduli,
                 int log_plaintext_modulus, double sigma)
      : log_n_(log_n),
        log_main_prime_moduli_(std::move(log_main_prime_moduli)),
        log_aux_prime_moduli_(std::move(log_aux_prime_moduli)),
        log_plaintext_modulus_(log_plaintext_modulus),
        sigma_(sigma) {
    log_main_modulus_ = std::accumulate(log_main_prime_moduli_.begin(),
                                        log_main_prime_moduli_.end(), 0);
    log_aux_modulus_ = std::accumulate(log_aux_prime_moduli_.begin(),
                                       log_aux_prime_moduli_.end(), 0);

    // Compute error constants.
    int dimension = 1 << log_n;
    b_plaintext_ = B_plaintext(dimension);
    b_secretkey_encryption_ = B_secretkey_encryption(dimension, sigma);
    b_publickey_encryption_ = B_publickey_encryption(dimension, sigma);
    b_scale_ = B_scale(dimension);
  }

  // This represents the "size" of an NTT coefficient of a randomly sampled
  // plaintext polynomial. The ciphertext error grows multiplicatively by this
  // constant under an absorb. Assume the plaintext polynomial has coefficients
  // chosen uniformly at random from the range [0, t), where t is the plaintext
  // modulus. Then the variance of a coefficient is  V = t ^ 2 / 12. In the NTT
  // domain, the variance is (dimension * t ^ 2 / 12).
  double B_plaintext(int dimension) const {
    // Return 6 * sqrt(V) where V is the variance of a coefficient in the NTT
    // domain.
    return ExportDoubleT() * sqrt(3.0 * dimension);
  }

  // This represents the "size" of a freshly encrypted ciphertext with a secret
  // key and error sampled from a centered binomial distribution with the
  // specified variance. The error and message have size |m + et|. Like
  // B_plaintext, the variance of a coefficient of m is V = t ^ 2 / 12, and the
  // variance of a coefficient of e is sigma ^ 2. In the NTT domain we can bound
  // the coefficient's variance by (dimension * (t ^ 2 / 12 + t * sigma)). The
  // bound 6 * sqrt(V) is then t * sqrt(dimension) * (sqrt(3.0) + 6.0 * sigma).
  double B_secretkey_encryption(int dimension, double sigma) const {
    return ExportDoubleT() * sqrt(dimension) * (sqrt(3.0) + 6.0 * sigma);
  }

  // This represents the "size" of a freshly encrypted ciphertext using a public
  // key, where the public key's error term, the public-key encryption's random
  // element and error terms are all sampled from a centered binomial
  // distribution with the specified standard deviation `sigma`. The error in a
  // fresh public-key encryption is t * (v * e + e' + s * e''), where s, v, e,
  // e', e'' are all sampled from the same error distribution of variance
  // sigma^2. In the NTT domain, the norm of this error term is bounded by t *
  // (72 * N * sigma^2 + 6 * sqrt(N) * sigma). Then adding the bound on the
  // message t * sqrt(3 * N) and we get the bound on the error and message.
  double B_publickey_encryption(int dimension, double sigma) const {
    return ExportDoubleT() * (sqrt(dimension) * (6 * sigma + sqrt(3)) +
                              72 * dimension * sigma * sigma);
  }

  // When modulus switching a ciphertext from a modulus q to a smaller modulus
  // p, the polynomial is scaled by (p / q) and a small rounding polynomial is
  // added so that the result is the closest integer polynomial with c' = c mod
  // t. The rounding polynomial's size contributes additively to the ciphertext
  // error, and its size is given by this constant.
  double B_scale(int dimension) const {
    return ExportDoubleT() *
           (sqrt(3.0 * dimension) + 8.0 * dimension * sqrt(1 / 3.0));
  }

  // Returns ceil(plaintext modulus).
  double ExportDoubleT() const { return std::pow(log_plaintext_modulus_, 2); }

  int log_n_;
  std::vector<int> log_main_prime_moduli_;
  std::vector<int> log_aux_prime_moduli_;
  int log_plaintext_modulus_;  // ceil of log2(plaintext modulus)
  double sigma_;

  int log_main_modulus_;
  int log_aux_modulus_;

  double b_plaintext_;
  double b_secretkey_encryption_;
  double b_publickey_encryption_;
  double b_scale_;
};

}  //  namespace rlwe

#endif  // RLWE_RNS_ERROR_PARAMS_H_
