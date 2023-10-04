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

#ifndef RLWE_RNS_RNS_CIPHERTEXT_H_
#define RLWE_RNS_RNS_CIPHERTEXT_H_

#include <memory>
#include <utility>
#include <vector>

#include "absl/strings/str_cat.h"
#include "absl/types/span.h"
#include "shell_encryption/montgomery.h"
#include "shell_encryption/rns/rns_error_params.h"
#include "shell_encryption/rns/rns_polynomial.h"
#include "shell_encryption/status_macros.h"

namespace rlwe {

// This class stores a RLWE ciphertext in the ring R_Q = Z[X]/(Q, X^N+1), where
// polynomials are represented in RNS form. This is the base RLWE ciphertext
// class which defines the basic data structure and operations that are common
// in all RLWE-based schemes.
//
// A RLWE ciphertext is a tuple of R_Q polynomials (c_0(X^j), c_1(X^j), ...,
// c_d(X^j)), where d is the degree of the ciphertext and j = `power_of_s` is
// the substitution power of the underlying secret polynomial s. Such a
// ciphertext encrypts a plaintext polynomial under the secret key (1, s(X^j),
// .., s^d(X^j)).
template <typename ModularInt>
class RnsRlweCiphertext {
 public:
  using ModularIntParams = typename ModularInt::Params;

  // Constructs a ciphertext that consists of the polynomials in `components`
  // mod Q where Q is specified in `moduli`. The ciphertext encrypts a plaintext
  // under the secret key (1, s(X^j), .., s^d(X^j)), where d = len(components) -
  // 1 and j = `power_of_s`.
  explicit RnsRlweCiphertext(
      std::vector<RnsPolynomial<ModularInt>> components,
      std::vector<const PrimeModulus<ModularInt>*> moduli, int power_of_s,
      double error, const RnsErrorParams<ModularInt>* error_params)
      : components_(std::move(components)),
        moduli_(std::move(moduli)),
        error_params_(error_params),
        power_of_s_(power_of_s),
        error_(error) {}

  static RnsRlweCiphertext CreateZero(
      std::vector<const PrimeModulus<ModularInt>*> moduli,
      const RnsErrorParams<ModularInt>* error_params) {
    return RnsRlweCiphertext({}, moduli, /*power_of_s=*/1, 0, error_params);
  }

  // Homomorphically add `that` to this ciphertext.
  // It is assumed that the two ciphertexts have the same modulus and the same
  // plaintext modulus.
  absl::Status AddInPlace(const RnsRlweCiphertext& that) {
    if (Degree() != that.Degree()) {
      return absl::InvalidArgumentError("`that` has a mismatched degree.");
    }
    if (Level() != that.Level()) {
      return absl::InvalidArgumentError("`that` has a mismatched level.");
    }
    if (PowerOfS() != that.PowerOfS()) {
      return absl::InvalidArgumentError(
          "`that` is encrypted with a different key power.");
    }
    for (int i = 0; i < components_.size(); ++i) {
      RLWE_RETURN_IF_ERROR(
          components_[i].AddInPlace(that.components_[i], moduli_));
    }
    error_ += that.error_;
    return absl::OkStatus();
  }

  // Homomorphically subtract `that` from this ciphertext.
  // It is assumed that the two ciphertexts have the same modulus and the same
  // plaintext modulus.
  absl::Status SubInPlace(const RnsRlweCiphertext& that) {
    if (Degree() != that.Degree()) {
      return absl::InvalidArgumentError("`that` has a mismatched degree.");
    }
    if (Level() != that.Level()) {
      return absl::InvalidArgumentError("`that` has a mismatched level.");
    }
    if (PowerOfS() != that.PowerOfS()) {
      return absl::InvalidArgumentError(
          "`that` is encrypted with a different key power.");
    }
    for (int i = 0; i < components_.size(); ++i) {
      RLWE_RETURN_IF_ERROR(
          components_[i].SubInPlace(that.components_[i], moduli_));
    }
    error_ += that.error_;
    return absl::OkStatus();
  }

  // Homomorphically multiply `plaintext` to this ciphertext.
  // It is assumed that the plaintext has coefficients modulo the plaintext
  // modulus of this ciphertext.
  absl::Status AbsorbInPlace(const RnsPolynomial<ModularInt>& plaintext) {
    for (auto& c : components_) {
      RLWE_RETURN_IF_ERROR(c.MulInPlace(plaintext, moduli_));
    }
    error_ *= error_params_->B_plaintext();
    return absl::OkStatus();
  }

  // Homomorphic Absorb of a monomial stored in a Polynomial. This function is
  // exactly like the homomorphic absorb, except that since the plaintext is a
  // monomial, the error estimate in the ciphertext is _not_ modified.
  absl::Status AbsorbMonomialInPlace(
      const RnsPolynomial<ModularInt>& monomial) {
    for (auto& c : components_) {
      RLWE_RETURN_IF_ERROR(c.MulInPlace(monomial, moduli_));
    }
    return absl::OkStatus();
  }

  // Homomorphically multiply the plaintext scalar to this ciphertext, where
  // the scalar is given in RNS representation `scalar_mod_qs`.
  // It is assumed that the scalar is wrt the same RNS moduli as this ciphertext
  // and that it represents a value bounded by the plaintext modulus.
  absl::StatusOr<RnsRlweCiphertext> Absorb(
      absl::Span<const ModularInt> scalar_mod_qs) const {
    RnsRlweCiphertext out = *this;
    RLWE_RETURN_IF_ERROR(out.AbsorbInPlace(scalar_mod_qs));
    return out;
  }

  absl::Status AbsorbInPlace(absl::Span<const ModularInt> scalar_mod_qs) {
    if (scalar_mod_qs.size() != moduli_.size()) {
      return absl::InvalidArgumentError(
          absl::StrCat("`scalar_mod_qs` must contain ", moduli_.size(),
                       " modular integers."));
    }
    for (auto& c : components_) {
      RLWE_RETURN_IF_ERROR(c.MulInPlace(scalar_mod_qs, moduli_));
    }
    error_ *= error_params_->B_plaintext();
    return absl::OkStatus();
  }

  // Homomorphically multiply `scalar` to of this ciphertext.
  // The scalar is assumed to be modulo the plaintext modulus.
  absl::StatusOr<RnsRlweCiphertext> Absorb(
      typename ModularInt::Int scalar) const {
    RnsRlweCiphertext out = *this;
    RLWE_RETURN_IF_ERROR(out.AbsorbInPlace(scalar));
    return out;
  }

  absl::Status AbsorbInPlace(typename ModularInt::Int scalar) {
    for (auto& c : components_) {
      RLWE_RETURN_IF_ERROR(c.MulInPlace(scalar, moduli_));
    }
    // The scalar is represented as a value in [-scalar/2, scalar/2].
    error_ *= static_cast<double>(scalar) / 2;
    return absl::OkStatus();
  }

  // Returns the ciphertext of the fused operation this + `a` * `b`.
  absl::Status FusedAbsorbAddInPlace(const RnsRlweCiphertext& a,
                                     const RnsPolynomial<ModularInt>& b) {
    if (Degree() != a.Degree()) {
      return absl::InvalidArgumentError("`a` has a mismatched degree.");
    }
    if (Level() != a.Level()) {
      return absl::InvalidArgumentError("`a` has a mismatched level.");
    }
    if (PowerOfS() != a.PowerOfS()) {
      return absl::InvalidArgumentError(
          "Ciphertexts must be encrypted with the same key power.");
    }

    // Compute components_[i] += a.components_[i] * b.
    for (int i = 0; i < components_.size(); i++) {
      RLWE_RETURN_IF_ERROR(
          components_[i].FusedMulAddInPlace(a.components_[i], b, moduli_));
    }

    // Update the error
    error_ += a.error_ * error_params_->B_plaintext();
    return absl::OkStatus();
  }

  // Homomorphic multiplication. This function returns the tensor product
  // of the two ciphertext vectors, and the resulting ciphertext has length
  // m + n - 1 where the two input ciphertexts have length m and n.
  // It is assumed that thw two ciphertexts have the same modulus.
  absl::StatusOr<RnsRlweCiphertext> Mul(const RnsRlweCiphertext& that) const {
    if (Level() != that.Level()) {
      return absl::InvalidArgumentError("`that` has a mismatched level.");
    }
    if (PowerOfS() != that.PowerOfS()) {
      return absl::InvalidArgumentError(
          "Ciphertexts must be encrypted with the same key power.");
    }

    // Create a vector of zero RNS polynomials
    RLWE_ASSIGN_OR_RETURN(auto zero, RnsPolynomial<ModularInt>::CreateZero(
                                         LogN(), moduli_, /*is_ntt=*/true));
    std::vector<RnsPolynomial<ModularInt>> result(
        components_.size() + that.components_.size() - 1, zero);
    // Compute the tensor product
    for (int i = 0; i < components_.size(); ++i) {
      for (int j = 0; j < that.components_.size(); ++j) {
        RLWE_RETURN_IF_ERROR(result[i + j].FusedMulAddInPlace(
            components_[i], that.components_[j], moduli_));
      }
    }
    return RnsRlweCiphertext(std::move(result), moduli_, power_of_s_,
                             error_ * that.error_, error_params_);
  }

  // Converts all components of this ciphertext to NTT form.
  absl::Status ConvertToNttForm() {
    for (auto& ci : components_) {
      RLWE_RETURN_IF_ERROR(ci.ConvertToNttForm(moduli_));
    }
    return absl::OkStatus();
  }

  // Converts all components of this ciphertext to Coefficient form.
  absl::Status ConvertToCoeffForm() {
    for (auto& ci : components_) {
      RLWE_RETURN_IF_ERROR(ci.ConvertToCoeffForm(moduli_));
    }
    return absl::OkStatus();
  }

  // Returns true if the modulus of this ciphertext equals to the modulus of
  // `that`. We consider equality of prime RNS moduli with order, that is, the
  // two ciphertexts have the same modulus if their prime RNS moduli are equal
  // and are in the same order.
  bool IsModulusSameAs(const RnsRlweCiphertext<ModularInt>& that) const {
    if (moduli_.size() != that.moduli_.size()) {
      return false;
    }
    for (int i = 0; i < moduli_.size(); ++i) {
      if (moduli_[i]->ModParams()->modulus !=
          that.moduli_[i]->ModParams()->modulus) {
        return false;
      }
    }
    return true;
  }

  ////////////////////////////////////////////////////////////////////////////////
  // Accessors
  ////////////////////////////////////////////////////////////////////////////////

  // Length of the ciphertext vector
  int Len() const { return components_.size(); }

  // The degree of this ciphertext is the number of member polynomials - 1.
  int Degree() const { return components_.size() - 1; }

  // The level of this ciphertext is the number of prime moduli - 1.
  int Level() const { return moduli_.size() - 1; }

  // Accessor for the ciphertext polynomial of the given index
  absl::StatusOr<RnsPolynomial<ModularInt>> Component(int index) const {
    if (index < 0 || index >= static_cast<int>(components_.size())) {
      return absl::InvalidArgumentError("Index out of range.");
    }
    return components_[index];
  }

  int LogN() const { return components_[0].LogN(); }

  int NumCoeffs() const { return components_[0].NumCoeffs(); }

  int NumModuli() const { return moduli_.size(); }

  absl::Span<const PrimeModulus<ModularInt>* const> Moduli() const {
    return moduli_;
  }

  const RnsErrorParams<ModularInt>* ErrorParams() const {
    return error_params_;
  }

  int PowerOfS() const { return power_of_s_; }

  double Error() const { return error_; }

 protected:
  // Allow derived ciphertext classes to access / mutate data members.
  const std::vector<RnsPolynomial<ModularInt>>& components() const {
    return components_;
  }
  std::vector<RnsPolynomial<ModularInt>>& components() { return components_; }

  const std::vector<const PrimeModulus<ModularInt>*>& moduli() const {
    return moduli_;
  }
  std::vector<const PrimeModulus<ModularInt>*>& moduli() { return moduli_; }

  double& error() { return error_; }

 private:
  // The ciphertext polynomials
  std::vector<RnsPolynomial<ModularInt>> components_;

  // The prime moduli constituting the modulus of this ciphertext.
  std::vector<const PrimeModulus<ModularInt>*> moduli_;

  // Error parameters.
  const RnsErrorParams<ModularInt>* error_params_;

  // The power a in s(x^a) that the ciphertext can be decrypted with.
  int power_of_s_;

  // A heuristic on the error of the ciphertext.
  double error_;
};

}  // namespace rlwe

#endif  // RLWE_RNS_RNS_CIPHERTEXT_H_
