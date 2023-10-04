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

#ifndef RLWE_RNS_RNS_SECRET_KEY_H_
#define RLWE_RNS_RNS_SECRET_KEY_H_

#include <algorithm>
#include <utility>
#include <vector>

#include "absl/types/span.h"
#include "shell_encryption/prng/prng.h"
#include "shell_encryption/rns/coefficient_encoder.h"
#include "shell_encryption/rns/error_distribution.h"
#include "shell_encryption/rns/rns_bgv_ciphertext.h"
#include "shell_encryption/rns/rns_ciphertext.h"
#include "shell_encryption/rns/rns_error_params.h"
#include "shell_encryption/rns/rns_integer.h"
#include "shell_encryption/rns/rns_polynomial.h"
#include "shell_encryption/status_macros.h"

namespace rlwe {

// This class holds a RLWE secret key represented by a RNS polynomial.
// A canonical RLWE secret key is (1, s) \in R^2, for R = Z[X]/(Q, X^N+1),
// where N is a power of 2 and Q is a product of prime moduli. The secret
// polynomial is thus s.
template <typename ModularInt>
class RnsRlweSecretKey {
 public:
  using ModularIntParams = typename ModularInt::Params;

  // Samples a secret key from the error distribution, and returns its
  // RNS representation wrt `moduli`.
  static absl::StatusOr<RnsRlweSecretKey> Sample(
      int log_n, int variance,
      std::vector<const PrimeModulus<ModularInt>*> moduli, SecurePrng* prng);

  // Returns a ciphertext that encrypts `messages` under this secret key as in
  // the BGV scheme, where `messages` are encoded using the given encoder, the
  // encryption randomness is supplied using `prng`, and the error parameters
  // are given in `error_params`.
  // Note that the encoder type is a template parameter, and by default we use
  // `CoefficientEncoder` to use messages as coefficients of the plaintext
  // polynomial.
  template <typename Encoder = CoefficientEncoder<ModularInt>>
  absl::StatusOr<RnsBgvCiphertext<ModularInt>> EncryptBgv(
      absl::Span<const typename ModularInt::Int> messages,
      const Encoder* encoder, const RnsErrorParams<ModularInt>* error_params,
      SecurePrng* prng) const;

  // Decrypts a BGV ciphertext using this secret key. The returned result is the
  // underlying messages decoded using `encoder`.
  template <typename Encoder = CoefficientEncoder<ModularInt>>
  absl::StatusOr<std::vector<typename ModularInt::Int>> DecryptBgv(
      const RnsBgvCiphertext<ModularInt>& ciphertext,
      const Encoder* encoder) const;

  // Alternative decoding that allows for additional parameter flexabillity
  // at the expense of being slower.
  template <typename Encoder = CoefficientEncoder<ModularInt>,
            typename BigInteger>
  absl::StatusOr<std::vector<typename ModularInt::Int>> DecryptBgvWithCrt(
      const RnsRlweCiphertext<ModularInt>& ciphertext, const Encoder* encoder,
      absl::Span<const BigInteger> modulus_hats,
      absl::Span<const ModularInt> modulus_hat_invs) const;

  // Reduces the modulus of the secret polynomial from Q = q0 * .. * ql to Q/ql.
  absl::Status ModReduce();

  // Accessors
  int LogN() const { return key_.LogN(); }
  int Level() const { return moduli_.size() - 1; }

  int NumCoeffs() const { return key_.NumCoeffs(); }
  int NumModuli() const { return moduli_.size(); }

  // Accessor for the prime moduli chain.
  absl::Span<const PrimeModulus<ModularInt>* const> Moduli() const {
    return moduli_;
  }

  // Accessor for the key polynomial
  const RnsPolynomial<ModularInt>& Key() const { return key_; }

  int Variance() const { return variance_; }

 private:
  explicit RnsRlweSecretKey(RnsPolynomial<ModularInt> key,
                            std::vector<const PrimeModulus<ModularInt>*> moduli,
                            int variance)
      : key_(std::move(key)), moduli_(std::move(moduli)), variance_(variance) {}

  // The "raw" decryption of RLWE scheme, which computes the inner product of
  // the canonical secret key (1, s) and a ciphertext (c0, c1) without removing
  // error nor decoding the result. Returns INVALID_ARGUMENT_ERROR if the input
  // `ciphertext` has degree other than 1 (i.e. its length is not 2).
  absl::StatusOr<RnsPolynomial<ModularInt>> RawDecrypt(
      const RnsRlweCiphertext<ModularInt>& ciphertext) const;

  // The key polynomial
  RnsPolynomial<ModularInt> key_;

  // The prime moduli constituting the modulus of this ciphertext.
  std::vector<const PrimeModulus<ModularInt>*> moduli_;

  // The variance of the binomial distribution of the key and error
  int variance_;
};

template <typename ModularInt>
template <typename Encoder>
absl::StatusOr<RnsBgvCiphertext<ModularInt>>
RnsRlweSecretKey<ModularInt>::EncryptBgv(
    absl::Span<const typename ModularInt::Int> messages, const Encoder* encoder,
    const RnsErrorParams<ModularInt>* error_params, SecurePrng* prng) const {
  if (encoder == nullptr) {
    return absl::InvalidArgumentError("`encoder` must not be null.");
  }
  if (error_params == nullptr) {
    return absl::InvalidArgumentError("`error_params` must not be null.");
  }
  if (prng == nullptr) {
    return absl::InvalidArgumentError("`prng` must not be null.");
  }

  // Encode messages into a plaintext polynomial.
  RLWE_ASSIGN_OR_RETURN(RnsPolynomial<ModularInt> plaintext,
                        encoder->EncodeBgv(messages, moduli_));
  if (!plaintext.IsNttForm()) {
    RLWE_RETURN_IF_ERROR(plaintext.ConvertToNttForm(moduli_));
  }

  // Sample a from the uniform distribution.
  RLWE_ASSIGN_OR_RETURN(
      auto a, RnsPolynomial<ModularInt>::SampleUniform(LogN(), prng, moduli_));

  // Sample the error term e (mod Q) from the error distribution.
  RLWE_ASSIGN_OR_RETURN(
      RnsPolynomial<ModularInt> c0,
      SampleError<ModularInt>(LogN(), variance_, moduli_, prng));

  // c0 = e * t (mod Q).
  RLWE_RETURN_IF_ERROR(c0.MulInPlace(encoder->PlaintextModulus(), moduli_));

  // c0 = e * t + m (mod Q).
  RLWE_RETURN_IF_ERROR(c0.AddInPlace(plaintext, moduli_));

  // c0 = e * t + m + a * key (mod Q).
  RLWE_RETURN_IF_ERROR(c0.FusedMulAddInPlace(a, key_, moduli_));

  // c1 = -a (mod Q).
  RLWE_RETURN_IF_ERROR(a.NegateInPlace(moduli_));

  return RnsBgvCiphertext<ModularInt>(
      {std::move(c0), std::move(a)}, moduli_,
      /*power=*/1, error_params->B_secretkey_encryption(), error_params);
}

template <typename ModularInt>
template <typename Encoder>
absl::StatusOr<std::vector<typename ModularInt::Int>>
RnsRlweSecretKey<ModularInt>::DecryptBgv(
    const RnsBgvCiphertext<ModularInt>& ciphertext,
    const Encoder* encoder) const {
  if (ciphertext.PowerOfS() != 1) {
    return absl::InvalidArgumentError(
        "Cannot decrypt `ciphertext` with power of s not equal to 1.");
  }
  if (ciphertext.LogN() != LogN()) {
    return absl::InvalidArgumentError(
        "Cannot decrypt `ciphertext` with a mismatching polynomial degree.");
  }
  if (ciphertext.NumModuli() != NumModuli()) {
    return absl::InvalidArgumentError(
        "Cannot decrypt `ciphertext` with a mismatching number"
        " of prime moduli.");
  }
  if (encoder == nullptr) {
    return absl::InvalidArgumentError("`encoder` must not be null.");
  }

  // Do the raw decryption to get the inner product of ciphertext and secret
  // key.
  RLWE_ASSIGN_OR_RETURN(RnsPolynomial<ModularInt> noisy_plaintext,
                        RawDecrypt(ciphertext));

  // Decode the noisy plaintext polynomial.
  return encoder->DecodeBgv(noisy_plaintext, moduli_);
}

// Decrypt with the alternative method `DecodeBgvWithCrt` in
// `FiniteFieldDecoder`.
template <typename ModularInt>
template <typename Encoder, typename BigInteger>
absl::StatusOr<std::vector<typename ModularInt::Int>>
RnsRlweSecretKey<ModularInt>::DecryptBgvWithCrt(
    const RnsRlweCiphertext<ModularInt>& ciphertext, const Encoder* encoder,
    absl::Span<const BigInteger> modulus_hats,
    absl::Span<const ModularInt> modulus_hat_invs) const {
  if (ciphertext.PowerOfS() != 1) {
    return absl::InvalidArgumentError(
        "Cannot decrypt `ciphertext` with power of s not equal to 1.");
  }
  if (ciphertext.LogN() != LogN()) {
    return absl::InvalidArgumentError(
        "Cannot decrypt `ciphertext` with a mismatching polynomial degree.");
  }
  if (ciphertext.NumModuli() != NumModuli() ||
      modulus_hats.size() != modulus_hat_invs.size() ||
      modulus_hats.size() != NumModuli()) {
    return absl::InvalidArgumentError(
        "Cannot decrypt `ciphertext` with a mismatching number"
        " of prime moduli.");
  }
  if (encoder == nullptr) {
    return absl::InvalidArgumentError("`encoder` must not be null.");
  }

  // Do the raw decryption to get the inner product of ciphertext and secret
  // key.
  RLWE_ASSIGN_OR_RETURN(RnsPolynomial<ModularInt> noisy_plaintext,
                        RawDecrypt(ciphertext));

  // Decode the noisy plaintext polynomial.
  return encoder->template DecodeBgvWithCrt<BigInteger>(
      std::move(noisy_plaintext), moduli_, modulus_hats, modulus_hat_invs);
}

}  // namespace rlwe

#endif  // RLWE_RNS_RNS_SECRET_KEY_H_
