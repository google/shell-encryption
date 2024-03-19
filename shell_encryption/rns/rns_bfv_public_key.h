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

#ifndef RLWE_RNS_RNS_BFV_PUBLIC_KEY_H_
#define RLWE_RNS_RNS_BFV_PUBLIC_KEY_H_

#include <vector>

#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "absl/strings/string_view.h"
#include "absl/types/span.h"
#include "shell_encryption/prng/prng.h"
#include "shell_encryption/rns/coefficient_encoder.h"
#include "shell_encryption/rns/rns_bfv_ciphertext.h"
#include "shell_encryption/rns/rns_error_params.h"
#include "shell_encryption/rns/rns_modulus.h"
#include "shell_encryption/rns/rns_polynomial.h"
#include "shell_encryption/rns/rns_public_key.h"
#include "shell_encryption/rns/rns_secret_key.h"
#include "shell_encryption/status_macros.h"

namespace rlwe {

template <typename ModularInt>
class RnsBfvPublicKey : public RnsRlwePublicKey<ModularInt> {
 public:
  using Integer = typename ModularInt::Int;

  // Allow copy and move, disallow copy-assign and move-assign.
  RnsBfvPublicKey(const RnsBfvPublicKey&) = default;
  RnsBfvPublicKey& operator=(const RnsBfvPublicKey&) = delete;
  RnsBfvPublicKey(RnsBfvPublicKey&&) = default;
  RnsBfvPublicKey& operator=(RnsBfvPublicKey&&) = delete;
  ~RnsBfvPublicKey() = default;

  // Generate a public key (b = a*s + t*e, -a) derived from the given secret
  // key, where the randomness a is freshly sampled uniform over the key's
  // modulus, and the error term e has coefficients sampled from a centered
  // binomial distribution of the given variance.
  static absl::StatusOr<RnsBfvPublicKey> Create(
      const RnsRlweSecretKey<ModularInt>& secret_key, int variance,
      PrngType prng_type) {
    RLWE_ASSIGN_OR_RETURN(
        RnsRlwePublicKey<ModularInt> public_key,
        RnsRlwePublicKey<ModularInt>::Create(secret_key, variance, prng_type));
    return RnsBfvPublicKey<ModularInt>(std::move(public_key));
  }

  // Returns a ciphertext that encrypts `messages` under this public key, where
  // `messages` are encoded using the given encoder, the encryption noises and
  // randomness have the given variance and are sampled using `prng`, and the
  // error parameters are given in `error_params`.
  // Note that the encoder type is a template parameter, and by default we use
  // `CoefficientEncoder` to use messages as coefficients of the plaintext
  // polynomial.
  template <typename Encoder = CoefficientEncoder<ModularInt>>
  absl::StatusOr<RnsBfvCiphertext<ModularInt>> Encrypt(
      absl::Span<const typename ModularInt::Int> messages,
      const Encoder* encoder, const RnsErrorParams<ModularInt>* error_params,
      SecurePrng* prng) const;

 private:
  explicit RnsBfvPublicKey(RnsRlwePublicKey<ModularInt> public_key)
      : RnsRlwePublicKey<ModularInt>(std::move(public_key)) {}
};

template <typename ModularInt>
template <typename Encoder>
absl::StatusOr<RnsBfvCiphertext<ModularInt>>
RnsBfvPublicKey<ModularInt>::Encrypt(
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
                        encoder->EncodeBfv(messages, this->Moduli()));

  if (!plaintext.IsNttForm()) {
    RLWE_RETURN_IF_ERROR(plaintext.ConvertToNttForm(this->Moduli()));
  }

  // Sample encryption randomness r.
  int log_n = this->LogN();
  RLWE_ASSIGN_OR_RETURN(
      RnsPolynomial<ModularInt> r,
      SampleError<ModularInt>(log_n, this->variance(), this->Moduli(), prng));

  // c0 = b * r + e' + Encode(messages).
  RLWE_ASSIGN_OR_RETURN(
      RnsPolynomial<ModularInt> c0,
      SampleError<ModularInt>(log_n, this->variance(), this->Moduli(), prng));
  RLWE_RETURN_IF_ERROR(c0.FusedMulAddInPlace(this->KeyB(), r, this->Moduli()));
  RLWE_RETURN_IF_ERROR(c0.AddInPlace(plaintext, this->Moduli()));

  // c1 = a * r + e''.
  RLWE_ASSIGN_OR_RETURN(
      RnsPolynomial<ModularInt> c1,
      SampleError<ModularInt>(log_n, this->variance(), this->Moduli(), prng));
  RLWE_RETURN_IF_ERROR(c1.FusedMulAddInPlace(this->KeyA(), r, this->Moduli()));

  return RnsBfvCiphertext<ModularInt>(
      {std::move(c0), std::move(c1)}, this->moduli(),
      /*power_of_s=*/1, error_params->B_publickey_encryption(), error_params,
      encoder->Context());
}

}  // namespace rlwe

#endif  // RLWE_RNS_RNS_BFV_PUBLIC_KEY_H_
