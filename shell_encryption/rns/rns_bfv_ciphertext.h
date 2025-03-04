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

#ifndef RLWE_RNS_RNS_BFV_CIPHERTEXT_H_
#define RLWE_RNS_RNS_BFV_CIPHERTEXT_H_

#include <memory>
#include <utility>
#include <vector>

#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "absl/types/span.h"
#include "shell_encryption/rns/rns_ciphertext.h"
#include "shell_encryption/rns/rns_context.h"
#include "shell_encryption/rns/rns_error_params.h"
#include "shell_encryption/rns/rns_integer.h"
#include "shell_encryption/rns/rns_modulus.h"
#include "shell_encryption/rns/rns_polynomial.h"
#include "shell_encryption/status_macros.h"

namespace rlwe {

// This class stores a BFV ciphertext in the ring R_Q = Z[X]/(Q, X^N+1), where
// polynomials are represented in RNS form.
template <typename ModularInt>
class RnsBfvCiphertext : public RnsRlweCiphertext<ModularInt> {
 public:
  using ModularIntParams = typename ModularInt::Params;

  // Constructs a BFV ciphertext that consists of the polynomials in
  // `components` mod Q where Q is specified in `moduli`. The ciphertext
  // encrypts a plaintext under the secret key (1, s(X^j), .., s^d(X^j)), where
  // d = len(components) - 1 and j = `power_of_s`.
  explicit RnsBfvCiphertext(std::vector<RnsPolynomial<ModularInt>> components,
                            std::vector<const PrimeModulus<ModularInt>*> moduli,
                            int power_of_s, double error,
                            const RnsErrorParams<ModularInt>* error_params,
                            const RnsContext<ModularInt>* context)
      : RnsRlweCiphertext<ModularInt>(std::move(components), std::move(moduli),
                                      power_of_s, error, error_params),
        context_(context) {}

  // Construct a BFV ciphertext from a plain RLWE ciphertext.
  explicit RnsBfvCiphertext(RnsRlweCiphertext<ModularInt> ciphertext)
      : RnsRlweCiphertext<ModularInt>(std::move(ciphertext)) {}

  static RnsBfvCiphertext CreateZero(
      std::vector<const PrimeModulus<ModularInt>*> moduli,
      const RnsErrorParams<ModularInt>* error_params,
      const RnsContext<ModularInt>* context) {
    return RnsBfvCiphertext({}, moduli, /*power_of_s=*/1, 0, error_params,
                            context);
  }

  // Returns the homomorphic negation of this ciphertext.
  absl::StatusOr<RnsBfvCiphertext> Negate() const {
    RnsBfvCiphertext out = *this;
    RLWE_RETURN_IF_ERROR(out.NegateInPlace());
    return out;
  }

  // Returns the homomorphic addition of this ciphertext with `that`.
  // It is assumed that the two ciphertexts have the same modulus and the same
  // plaintext modulus.
  absl::StatusOr<RnsBfvCiphertext> operator+(
      const RnsBfvCiphertext& that) const {
    RnsBfvCiphertext out = *this;
    RLWE_RETURN_IF_ERROR(out.AddInPlace(that));
    return out;
  }

  // Returns the homomorphic addition of this ciphertext with `plaintext`.
  // It is assumed that `plaintext` is a polynomial in the plaintext space.
  absl::StatusOr<RnsBfvCiphertext> operator+(
      const RnsPolynomial<ModularInt>& plaintext) const {
    RnsBfvCiphertext out = *this;
    RLWE_RETURN_IF_ERROR(out.AddInPlace(plaintext));
    return out;
  }

  // Returns the homomorphic subtraction of this ciphertext by `that`.
  // It is assumed that the two ciphertexts have the same modulus and the same
  // plaintext modulus.
  absl::StatusOr<RnsBfvCiphertext> operator-(
      const RnsBfvCiphertext& that) const {
    RnsBfvCiphertext out = *this;
    RLWE_RETURN_IF_ERROR(out.SubInPlace(that));
    return out;
  }

  // Returns the homomorphic subtraction of this ciphertext by `plaintext`.
  // It is assumed that `plaintext` is a polynomial in the plaintext space.
  absl::StatusOr<RnsBfvCiphertext> operator-(
      const RnsPolynomial<ModularInt>& plaintext) const {
    RnsBfvCiphertext out = *this;
    RLWE_RETURN_IF_ERROR(out.SubInPlace(plaintext));
    return out;
  }

  absl::StatusOr<RnsBfvCiphertext> SubWithoutPad(
      const RnsBfvCiphertext& that) const {
    RnsBfvCiphertext out = *this;
    RLWE_RETURN_IF_ERROR(out.SubInPlaceWithoutPad(that));
    return out;
  }

  // Returns the homomorphic absorbtion of this ciphertext and a plaintext.
  // It is assumed that the plaintext has coefficients modulo the plaintext
  // modulus of this ciphertext.
  absl::StatusOr<RnsBfvCiphertext> operator*(
      const RnsPolynomial<ModularInt>& plaintext) const {
    RnsBfvCiphertext out = *this;
    RLWE_RETURN_IF_ERROR(out.AbsorbInPlace(plaintext));
    return out;
  }

  // Returns the homomorphic absorbtion of this ciphertext and a plaintext
  // scalar given its RNS representation `scalar_mod_qs`.
  // It is assumed that the scalar is wrt the same RNS moduli as this ciphertext
  // and that it represents a value bounded by the plaintext modulus.
  absl::StatusOr<RnsBfvCiphertext> operator*(
      absl::Span<const ModularInt> scalar_mod_qs) const {
    RnsBfvCiphertext out = *this;
    RLWE_RETURN_IF_ERROR(out.AbsorbInPlace(scalar_mod_qs));
    return out;
  }

  // Returns the homomorphic absorption of this ciphertext and `scalar`.
  // The scalar is assumed to be modulo the plaintext modulus.
  absl::StatusOr<RnsBfvCiphertext> operator*(
      typename ModularInt::Int scalar) const {
    RnsBfvCiphertext out = *this;
    RLWE_RETURN_IF_ERROR(out.AbsorbInPlace(scalar));
    return out;
  }

  // Homomorphic multiplication. This function returns the tensor product
  // of the two ciphertext vectors, and the resulting ciphertext has length
  // m + n - 1 where the two input ciphertexts have length m and n.
  // It is assumed that the two ciphertexts have the same modulus.
  absl::StatusOr<RnsBfvCiphertext> operator*(
      const RnsBfvCiphertext& that) const {
    return Mul(that);
  }

  absl::StatusOr<RnsBfvCiphertext> Mul(const RnsBfvCiphertext& that) const;

  // Make all variants of `AbsorbInPlace` from `RnsRlweCiphertext` visible.
  using RnsRlweCiphertext<ModularInt>::AbsorbInPlace;

  // Homomorphically multiply `plaintext` to this ciphertext, assuming
  // `plaintext` is a polynomial that encodes the plaintext messages, i.e.
  // the messages are round(t / Q * `plaintext`) (mod t).
  absl::Status AbsorbInPlace(const RnsPolynomial<ModularInt>& plaintext);

  // Homomorphically multiply `plaintext` to this ciphertext, assuming
  // `plaintext` is a polynomial that encodes the plaintext messages without
  // scaling, i.e. the messages are `plaintext` (mod t).
  absl::Status AbsorbInPlaceSimple(const RnsPolynomial<ModularInt>& plaintext);

  absl::StatusOr<RnsBfvCiphertext> AbsorbSimple(
      const RnsPolynomial<ModularInt>& plaintext) const {
    RnsBfvCiphertext out = *this;
    RLWE_RETURN_IF_ERROR(out.AbsorbInPlaceSimple(plaintext));
    return out;
  }

  // Performs modulus reduction on BFV ciphertext. If `t` is the plaintext
  // modulus, and q_L is the last prime modulus in the moduli chain, then the
  // modulus reduced ciphertext is [round((q_L/Q)*c_i) mod (Q/q_L) : i = 0..d]
  // where Q = q_0 * .. * q_L is the product of all prime moduli, d is the
  // degree of this ciphertext [c_0,..,c_d].
  absl::Status ModReduce() {
    if (this->moduli().size() <= 1) {
      return absl::FailedPreconditionError(
          "Cannot perform ModReduce with insufficient number of prime moduli.");
    }

    int level = this->Level();
    RnsInt<ModularInt> ql_inv =
        context_->MainPrimeModulusInverseResidues()[level];
    for (RnsPolynomial<ModularInt>& c : this->components()) {
      RLWE_RETURN_IF_ERROR(
          c.ModReduceMsb(ql_inv.Prefix(level), this->moduli()));
    }
    this->moduli().pop_back();
    return absl::OkStatus();
  }

  // Returns the ciphertext where all polynomials are substituted with the
  // given power, ie [c0(X^j), c1(X^j), ...]
  absl::StatusOr<RnsBfvCiphertext> Substitute(int substitution_power) const {
    std::vector<RnsPolynomial<ModularInt>> subbed_components;
    subbed_components.reserve(this->components().size());
    for (const RnsPolynomial<ModularInt>& c : this->components()) {
      RLWE_ASSIGN_OR_RETURN(auto subbed_c,
                            c.Substitute(substitution_power, this->moduli()));
      subbed_components.push_back(std::move(subbed_c));
    }
    int power_of_s = (this->PowerOfS() * substitution_power) %
                     (2 * this->components()[0].NumCoeffs());
    return RnsBfvCiphertext(subbed_components, this->moduli(), power_of_s,
                            this->Error(), this->ErrorParams(), context_);
  }

  absl::StatusOr<RnsBfvCiphertext> SubstituteWithoutPad(
      int substitution_power) const {
    std::vector<RnsPolynomial<ModularInt>> subbed_components;
    subbed_components.reserve(this->components().size());

    RLWE_ASSIGN_OR_RETURN(
        auto subbed_c,
        this->components()[0].Substitute(substitution_power, this->moduli()));
    subbed_components.push_back(std::move(subbed_c));

    int power_of_s = (this->PowerOfS() * substitution_power) %
                     (2 * this->components()[0].NumCoeffs());
    return RnsBfvCiphertext(subbed_components, this->moduli(), power_of_s,
                            this->Error(), this->ErrorParams(), context_);
  }

  // Accessors.
  const RnsContext<ModularInt>* Context() const { return context_; }

 private:
  // The RNS context which contains pre-computed modulus constants.
  const RnsContext<ModularInt>* context_;
};

}  // namespace rlwe

#endif  // RLWE_RNS_RNS_BFV_CIPHERTEXT_H_
