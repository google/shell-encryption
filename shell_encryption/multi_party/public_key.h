/*
 * Copyright 2025 Google LLC.
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

#ifndef RLWE_MULTI_PARTY_PUBLIC_KEY_H_
#define RLWE_MULTI_PARTY_PUBLIC_KEY_H_

#include <vector>

#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "absl/types/span.h"
#include "shell_encryption/multi_party/public_key_share.h"
#include "shell_encryption/multi_party/public_parameter.h"
#include "shell_encryption/prng/prng.h"
#include "shell_encryption/rns/coefficient_encoder.h"
#include "shell_encryption/rns/error_distribution.h"
#include "shell_encryption/rns/rns_bfv_ciphertext.h"
#include "shell_encryption/rns/rns_error_params.h"
#include "shell_encryption/rns/rns_modulus.h"
#include "shell_encryption/rns/rns_polynomial.h"
#include "shell_encryption/status_macros.h"

namespace rlwe {
namespace multi_party {

// Public key of the RLWE-based multi-party homomorphic encryption scheme,
// which is the sum of public key shares from all parties in the protocol.
template <typename ModularInt>
class PublicKey {
 public:
  using Integer = typename ModularInt::Int;

  // Returns a public key given the public parameter and all parties' public key
  // shares.
  static absl::StatusOr<PublicKey> Create(
      const PublicParameter<ModularInt>* public_parameter,
      absl::Span<const PublicKeyShare<ModularInt>> public_key_shares);

  static absl::StatusOr<PublicKey> Deserialize(
      const SerializedPublicKey& serialized,
      const PublicParameter<ModularInt>* public_parameter);

  absl::StatusOr<SerializedPublicKey> Serialize() const {
    SerializedPublicKey serialized;
    RLWE_ASSIGN_OR_RETURN(*serialized.mutable_key_b(),
                          key_b_.Serialize(public_parameter_->Moduli()));
    return serialized;
  }

  // For Rust interoperability. Defined in the wrapper library, not here.
  friend class PublicKeyRawFactory;

  // Returns a ciphertext that encrypts `messages` under this public key, where
  // `messages` are encoded using the given encoder, the encryption noises and
  // randomness have the given variance and are sampled using `prng`, and the
  // error parameters are given in `error_params`.
  // Note that the encoder type is a template parameter, and by default we use
  // `CoefficientEncoder` to use messages as coefficients of the plaintext
  // polynomial.
  template <typename Encoder = CoefficientEncoder<ModularInt>>
  absl::StatusOr<RnsBfvCiphertext<ModularInt>> Encrypt(
      absl::Span<const Integer> messages, const Encoder* encoder,
      const RnsErrorParams<ModularInt>* error_params, SecurePrng* prng) const;

  // Encodes and encrypts `messages` under this public key. Stores the raw
  // components of the ciphertext in `ciphertext_component_b` (a.k.a. ct0) and
  // `ciphertext_component_a` (a.k.a. ct1). `ciphertext_secret_r` and
  // `ciphertext_error_e` are optional and can be nullptr. If they are not
  // nullptr, they will be populated with the secret randomness and the error.
  template <typename Encoder = CoefficientEncoder<ModularInt>>
  absl::Status EncryptExplicit(
      absl::Span<const typename ModularInt::Int> messages,
      const Encoder* encoder, const RnsErrorParams<ModularInt>* error_params,
      SecurePrng* prng, RnsPolynomial<ModularInt>* ciphertext_component_b,
      RnsPolynomial<ModularInt>* ciphertext_component_a,
      RnsPolynomial<ModularInt>* ciphertext_secret_r,
      RnsPolynomial<ModularInt>* ciphertext_error_e) const;

  // Accessor to the "b" component in a public key.
  const RnsPolynomial<ModularInt>& ComponentB() const { return key_b_; }

 private:
  explicit PublicKey(const PublicParameter<ModularInt>* public_parameter,
                     RnsPolynomial<ModularInt> key_b)
      : public_parameter_(public_parameter), key_b_(std::move(key_b)) {}

  const PublicParameter<ModularInt>* public_parameter_;

  // The "b" component of a RLWE public key.
  RnsPolynomial<ModularInt> key_b_;
};

template <typename ModularInt>
template <typename Encoder>
absl::StatusOr<RnsBfvCiphertext<ModularInt>> PublicKey<ModularInt>::Encrypt(
    absl::Span<const typename ModularInt::Int> messages, const Encoder* encoder,
    const RnsErrorParams<ModularInt>* error_params, SecurePrng* prng) const {
  // Compute ciphertext components, discard error e'' and secret r.
  auto ciphertext_component_b =
      RnsPolynomial<ModularInt>::CreateEmpty();  // ct0
  auto ciphertext_component_a =
      RnsPolynomial<ModularInt>::CreateEmpty();  // ct1
  RLWE_RETURN_IF_ERROR(EncryptExplicit(
      messages, encoder, error_params, prng, &ciphertext_component_b,
      &ciphertext_component_a, /*ciphertext_secret_r=*/nullptr,
      /*ciphertext_error_e=*/nullptr));

  auto moduli = public_parameter_->Moduli();
  std::vector<const PrimeModulus<ModularInt>*> moduli_vector;
  moduli_vector.insert(moduli_vector.begin(), moduli.begin(), moduli.end());

  return RnsBfvCiphertext<ModularInt>(
      {std::move(ciphertext_component_b), std::move(ciphertext_component_a)},
      std::move(moduli_vector),
      /*power_of_s=*/1, error_params->B_publickey_encryption(), error_params,
      encoder->Context());
}

template <typename ModularInt>
template <typename Encoder>
absl::Status PublicKey<ModularInt>::EncryptExplicit(
    absl::Span<const typename ModularInt::Int> messages, const Encoder* encoder,
    const RnsErrorParams<ModularInt>* error_params, SecurePrng* prng,
    RnsPolynomial<ModularInt>* ciphertext_component_b,
    RnsPolynomial<ModularInt>* ciphertext_component_a,
    RnsPolynomial<ModularInt>* ciphertext_secret_r,
    RnsPolynomial<ModularInt>* ciphertext_error_e) const {
  if (encoder == nullptr) {
    return absl::InvalidArgumentError("`encoder` must not be null.");
  }
  if (error_params == nullptr) {
    return absl::InvalidArgumentError("`error_params` must not be null.");
  }
  if (prng == nullptr) {
    return absl::InvalidArgumentError("`prng` must not be null.");
  }
  if (ciphertext_component_a == nullptr) {
    return absl::InvalidArgumentError(
        "`ciphertext_component_a` must not be null.");
  }
  if (ciphertext_component_b == nullptr) {
    return absl::InvalidArgumentError(
        "`ciphertext_component_b` must not be null.");
  }

  // Encode message into a plaintext polynomial.
  auto moduli = public_parameter_->Moduli();
  Integer t = encoder->PlaintextModulus();
  RLWE_ASSIGN_OR_RETURN(RnsPolynomial<ModularInt> plaintext,
                        encoder->EncodeBfv(messages, moduli));
  if (!plaintext.IsNttForm()) {
    RLWE_RETURN_IF_ERROR(plaintext.ConvertToNttForm(moduli));
  }

  // Sample encryption randomness r and optionally keep a copy of it.
  int log_n = public_parameter_->LogN();
  RLWE_ASSIGN_OR_RETURN(RnsPolynomial<ModularInt> r,
                        SampleUniformTernary<ModularInt>(log_n, moduli, prng));
  if (ciphertext_secret_r != nullptr) {
    *ciphertext_secret_r = r;
  }

  // c0 = b * r + e' + Encode(message).
  RLWE_ASSIGN_OR_RETURN(
      *ciphertext_component_b,
      SampleError<ModularInt>(log_n, public_parameter_->ErrorVariance(), moduli,
                              prng));
  RLWE_RETURN_IF_ERROR(
      ciphertext_component_b->FusedMulAddInPlace(key_b_, r, moduli));
  RLWE_RETURN_IF_ERROR(ciphertext_component_b->AddInPlace(plaintext, moduli));

  // Sample e'' and optionally keep a copy of it.
  RLWE_ASSIGN_OR_RETURN(
      *ciphertext_component_a,
      SampleError<ModularInt>(log_n, public_parameter_->ErrorVariance(), moduli,
                              prng));
  if (ciphertext_error_e != nullptr) {
    *ciphertext_error_e = *ciphertext_component_a;
  }

  // c1 = a * r + e''.
  RLWE_RETURN_IF_ERROR(ciphertext_component_a->FusedMulAddInPlace(
      public_parameter_->PublicKeyComponentA(), r, moduli));

  return absl::OkStatus();
}

}  // namespace multi_party
}  // namespace rlwe

#endif  // RLWE_MULTI_PARTY_PUBLIC_KEY_H_
