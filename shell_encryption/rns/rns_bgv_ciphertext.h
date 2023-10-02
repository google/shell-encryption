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

#ifndef RLWE_RNS_RNS_BGV_CIPHERTEXT_H_
#define RLWE_RNS_RNS_BGV_CIPHERTEXT_H_

#include <memory>
#include <utility>
#include <vector>

#include "absl/strings/str_cat.h"
#include "absl/types/span.h"
#include "shell_encryption/montgomery.h"
#include "shell_encryption/rns/rns_ciphertext.h"
#include "shell_encryption/rns/rns_error_params.h"
#include "shell_encryption/rns/rns_polynomial.h"
#include "shell_encryption/status_macros.h"

namespace rlwe {

// This class stores a BGV ciphertext in the ring R_Q = Z[X]/(Q, X^N+1), where
// polynomials are represented in RNS form.
template <typename ModularInt>
class RnsBgvCiphertext : public RnsRlweCiphertext<ModularInt> {
 public:
  using ModularIntParams = typename ModularInt::Params;

  // Constructs a BGV ciphertext that consists of the polynomials in
  // `components` mod Q where Q is specified in `moduli`. The ciphertext
  // encrypts a plaintext under the secret key (1, s(X^j), .., s^d(X^j)), where
  // d = len(components) - 1 and j = `power_of_s`.
  explicit RnsBgvCiphertext(std::vector<RnsPolynomial<ModularInt>> components,
                            std::vector<const PrimeModulus<ModularInt>*> moduli,
                            int power_of_s, double error,
                            const RnsErrorParams<ModularInt>* error_params)
      : RnsRlweCiphertext<ModularInt>(std::move(components), std::move(moduli),
                                      power_of_s, error, error_params) {}

  // Construct a BGV ciphertext from a plain RLWE ciphertext.
  explicit RnsBgvCiphertext(RnsRlweCiphertext<ModularInt> ciphertext)
      : RnsRlweCiphertext<ModularInt>(std::move(ciphertext)) {}

  // Returns the homomorphic addition of this ciphertext with `that`.
  // It is assumed that the two ciphertexts have the same modulus and the same
  // plaintext modulus.
  absl::StatusOr<RnsBgvCiphertext> operator+(
      const RnsBgvCiphertext& that) const {
    RnsBgvCiphertext out = *this;
    RLWE_RETURN_IF_ERROR(out.AddInPlace(that));
    return out;
  }

  // Returns the homomorphic subtraction of this ciphertext by `that`.
  // It is assumed that the two ciphertexts have the same modulus and the same
  // plaintext modulus.
  absl::StatusOr<RnsBgvCiphertext> operator-(
      const RnsBgvCiphertext& that) const {
    RnsBgvCiphertext out = *this;
    RLWE_RETURN_IF_ERROR(out.SubInPlace(that));
    return out;
  }

  // Returns the homomorphic absorbtion of this ciphertext and a plaintext.
  // It is assumed that the plaintext has coefficients modulo the plaintext
  // modulus of this ciphertext.
  absl::StatusOr<RnsBgvCiphertext> operator*(
      const RnsPolynomial<ModularInt>& plaintext) const {
    RnsBgvCiphertext out = *this;
    RLWE_RETURN_IF_ERROR(out.AbsorbInPlace(plaintext));
    return out;
  }

  // Returns the homomorphic absorbtion of this ciphertext and a plaintext
  // scalar given its RNS representation `scalar_mod_qs`.
  // It is assumed that the scalar is wrt the same RNS moduli as this ciphertext
  // and that it represents a value bounded by the plaintext modulus.
  absl::StatusOr<RnsBgvCiphertext> operator*(
      absl::Span<const ModularInt> scalar_mod_qs) const {
    RnsBgvCiphertext out = *this;
    RLWE_RETURN_IF_ERROR(out.AbsorbInPlace(scalar_mod_qs));
    return out;
  }

  // Returns the homomorphic absorbtion of this ciphertext and `scalar`.
  // The scalar is assumed to be modulo the plaintext modulus.
  absl::StatusOr<RnsBgvCiphertext> operator*(
      typename ModularInt::Int scalar) const {
    RnsBgvCiphertext out = *this;
    RLWE_RETURN_IF_ERROR(out.AbsorbInPlace(scalar));
    return out;
  }

  // Homomorphic multiplication. This function returns the tensor product
  // of the two ciphertext vectors, and the resulting ciphertext has length
  // m + n - 1 where the two input ciphertexts have length m and n.
  // It is assumed that the two ciphertexts have the same modulus.
  absl::StatusOr<RnsBgvCiphertext> operator*(
      const RnsBgvCiphertext& that) const {
    RLWE_ASSIGN_OR_RETURN(RnsRlweCiphertext<ModularInt> product,
                          this->Mul(that));
    return RnsBgvCiphertext<ModularInt>(std::move(product));
  }

  // Performs modulus reduction on BGV ciphertext. If `t` is the plaintext
  // modulus, and q_L is the last prime modulus in the moduli chain, then the
  // modulus reduced ciphertext is [round((q_L/Q)*c_i) mod (Q/q_L) : i = 0..d]
  // where Q = q_0 * .. * q_L is the product of all prime moduli, d is the
  // degree of this ciphertext [c_0,..,c_d].
  absl::Status ModReduce(typename ModularInt::Int t,
                         const RnsInt<ModularInt>& ql_inv) {
    if (this->moduli().size() <= 1) {
      return absl::FailedPreconditionError(
          "Cannot perform ModReduce with insufficient number of prime moduli.");
    }
    for (RnsPolynomial<ModularInt>& c : this->components()) {
      RLWE_RETURN_IF_ERROR(c.ModReduceLsb(t, ql_inv, this->moduli()));
    }
    this->moduli().pop_back();
    return absl::OkStatus();
  }

  // Returns the ciphertext where all polynomials are substituted with the
  // given power, ie [c0(X^j), c1(X^j), ...]
  absl::StatusOr<RnsBgvCiphertext> Substitute(int substitution_power) const {
    std::vector<RnsPolynomial<ModularInt>> subbed_components;
    subbed_components.reserve(this->components().size());
    for (const RnsPolynomial<ModularInt>& c : this->components()) {
      RLWE_ASSIGN_OR_RETURN(auto subbed_c,
                            c.Substitute(substitution_power, this->moduli()));
      subbed_components.push_back(std::move(subbed_c));
    }
    int power_of_s = (this->PowerOfS() * substitution_power) %
                     (2 * this->components()[0].NumCoeffs());

    return RnsBgvCiphertext(subbed_components, this->moduli(), power_of_s,
                            this->Error(), this->ErrorParams());
  }
};

}  // namespace rlwe

#endif  // RLWE_RNS_RNS_BGV_CIPHERTEXT_H_
