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

#ifndef RLWE_RNS_APPROXIMATE_ENCODER_H_
#define RLWE_RNS_APPROXIMATE_ENCODER_H_

#include <complex>
#include <utility>
#include <vector>

#include "absl/numeric/int128.h"
#include "absl/status/statusor.h"
#include "absl/types/span.h"
#include "shell_encryption/int256.h"
#include "shell_encryption/integral_types.h"
#include "shell_encryption/rns/rns_context.h"
#include "shell_encryption/rns/rns_modulus.h"
#include "shell_encryption/rns/rns_polynomial.h"

namespace rlwe {

namespace {

// To capture the big integer type that can hold coefficient values of
// CKKS plaintext polynomials during encoding and decoding.
template <typename T>
struct EncodingModulus;

template <>
struct EncodingModulus<Uint64> {
  using value_type = uint256;
};

template <>
struct EncodingModulus<Uint128> {
  using value_type = uint256;
};

#ifdef ABSL_HAVE_INTRINSIC_INT128
template <>
struct EncodingModulus<absl::uint128> {
  using value_type = uint256;
};
#endif

}  // namespace

// This class specifies how to encode a complex vector in CC^{N/2} as an integer
// polynomial with representation in Z[X]/(Q, X^N+1), where Q is the ciphertext
// modulus.
//
// In CKKS a complex vector z in CC^{N/2} is identified as the image of
// a real polynomial f(X) under mapping \pi * \tau, where \tau is the canonical
// embedding map which takes a polynomial f(X) to its evaluations {a(omega_j)}_j
// under all primitive 2N'th root of unity, i.e. all roots of X^N+1, and \pi
// maps pairs of conjugate complex numbers (a + bI, a - bI) to a + bI. To gain
// enough precision and represent f(X) as a discrete element, we further scale
// up f(X) and round its coefficients to integers. So, specifically,
// - Encode(z) = round(scaling_factor * \tau^-1(\pi^-1(z)));
// - Decode(f) = 1/scaling_factor * \pi(\tau(f)).
// Note that Decode is not exactly the inverse of Encode, and precision loss
// will happen due to rounding error.
//
// When combining the encoding scheme described above with RLWE encryption to
// build CKKS scheme, secret-dependent information may be leaked in decryption
// results even when considering just passive adversary, as described in
// https://eprint.iacr.org/2020/1533.pdf. In this library we add Gaussian noise
// as described in https://eprint.iacr.org/2022/816.
template <typename ModularInt>
class ApproximateEncoder {
 public:
  using Integer = typename ModularInt::Int;
  using BigInteger = typename EncodingModulus<Integer>::value_type;

  // Returns an ApproximateEncoder supporting CKKS encoding with a given scaling
  // factor.
  static absl::StatusOr<ApproximateEncoder<ModularInt>> Create(
      const RnsContext<ModularInt>* context, double scaling_factor);

  // Encodes a vector of complex numbers into a polynomial f(X) mod Q, where
  // the modulus Q is given in the RNS `moduli`.
  absl::StatusOr<RnsPolynomial<ModularInt>> EncodeCkks(
      absl::Span<const std::complex<double>> values,
      absl::Span<const PrimeModulus<ModularInt>* const> moduli) const;

  // Decodes a polynomial f(X) mod Q and returns the vector of complex numbers
  // encoded by f(X), where Q is given in the RNS `moduli`.
  absl::StatusOr<std::vector<std::complex<double>>> DecodeCkks(
      RnsPolynomial<ModularInt> noisy_plaintext,
      absl::Span<const PrimeModulus<ModularInt>* const> moduli) const;

 private:
  explicit ApproximateEncoder(const RnsContext<ModularInt>* context,
                              double scaling_factor,
                              std::vector<std::complex<double>> psis_bitrev,
                              std::vector<std::complex<double>> psis_bitrev_inv,
                              std::vector<unsigned int> bitrevs)
      : context_(context),
        scaling_factor_(scaling_factor),
        psis_bitrev_(std::move(psis_bitrev)),
        psis_bitrev_inv_(std::move(psis_bitrev_inv)),
        bitrevs_(std::move(bitrevs)) {}

  // Computes the forward DFT \pi * \tau used in CKKS decoding.
  absl::StatusOr<std::vector<std::complex<double>>> ForwardTransform(
      std::vector<std::complex<double>> coeffs) const;

  // Computes the inverse DFT \tau^-1 * \pi^-1 used in CKKS encoding.
  absl::StatusOr<std::vector<std::complex<double>>> InverseTransform(
      absl::Span<const std::complex<double>> slots) const;

  // Doesn't own the RnsContext, which must live longer than this encoder.
  const RnsContext<ModularInt>* context_;

  // The scaling factor.
  double scaling_factor_;

  // The primitive complex roots of unity, in bit reversed order.
  std::vector<std::complex<double>> psis_bitrev_;
  std::vector<std::complex<double>> psis_bitrev_inv_;

  // Bit reversal indices for a vector of length N/2.
  std::vector<unsigned int> bitrevs_;
};

}  // namespace rlwe

#endif  // RLWE_RNS_APPROXIMATE_ENCODER_H_
