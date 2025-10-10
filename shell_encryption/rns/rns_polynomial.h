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

#ifndef RLWE_RNS_RNS_POLYNOMIAL_H_
#define RLWE_RNS_RNS_POLYNOMIAL_H_

#include <utility>
#include <vector>

#include "absl/strings/str_cat.h"
#include "absl/types/span.h"
#include "shell_encryption/montgomery.h"
#include "shell_encryption/ntt_parameters.h"
#include "shell_encryption/polynomial.h"
#include "shell_encryption/rns/rns_integer.h"
#include "shell_encryption/rns/rns_modulus.h"
#include "shell_encryption/rns/serialization.pb.h"
#include "shell_encryption/status_macros.h"

namespace rlwe {

// This class stores a polynomial in Z[X]/(Q, X^N+1) using RNS representation,
// where N is a power of two and Q = q_0 * .. * q_{L-1} for distinct primes
// q_0,..,q_{L-1} such that 2N | q_i - 1 and q_i can fit into the MontgomeryInt
// type `ModularInt`. The RNS representation of a polynomial a(X) in this ring
// is a vector of polynomials a_i(X) = a(X) mod q_i, for i = 0,..,L-1.
// RnsPolynomial stores the coefficient vectors of these sub-polynomials a_i(X).
//
// By default, a RNS polynomial is stored in the NTT form, where the coefficient
// vectors are evaluations [a_i(\omega_j), j = 1..N] for \omega_j the primitive
// 2N's root of unity mod q_i. It can be converted to the coefficient form, in
// which case the coefficient vectors contain coefficients of a_i(X) wrt the
// power basis {1, X, .., X^{N-1}}. Note that certain operations require the
// polynomial be in NTT form, e.g. power substitution.
template <typename ModularInt>
class RnsPolynomial {
 public:
  using ModularIntParams = typename ModularInt::Params;
  using NttParams = NttParameters<ModularInt>;

  // Creates a polynomial whose RNS representation is given in `coeff_vectors`.
  static absl::StatusOr<RnsPolynomial> Create(
      std::vector<std::vector<ModularInt>> coeff_vectors, bool is_ntt = true);

  // Creates an empty polynomial with empty coefficient vectors. This empty
  // polynomial should only be used for placeholder purposes (e.g. initialize a
  // return argument) and is not meant to be operated on. Note: this empty
  // polynomial has `log_n_` = 0 and thus `NumCoeffs()` = 1.
  static RnsPolynomial CreateEmpty();

  // Creates a polynomial 0 wrt the given moduli representing q_i's.
  static absl::StatusOr<RnsPolynomial> CreateZero(
      int log_n, absl::Span<const PrimeModulus<ModularInt>* const> moduli,
      bool is_ntt = true);

  // Creates a polynomial 1 wrt the given moduli, in NTT form.
  static absl::StatusOr<RnsPolynomial> CreateOne(
      int log_n, absl::Span<const PrimeModulus<ModularInt>* const> moduli);

  // Returns a randomly sampled RNS polynomial wrt the given the prime moduli.
  template <typename Prng = rlwe::SecurePrng>
  static absl::StatusOr<RnsPolynomial> SampleUniform(
      int log_n, Prng* prng,
      absl::Span<const PrimeModulus<ModularInt>* const> moduli);

  // Returns a RnsPolynomial representing a polynomial in Z[X]/(X^N+1), whose
  // coefficients are given in `coeffs_q` in balanced form modulo-q.
  // The returned RNS polynomial is in NTT form wrt `moduli`.
  // Note: A value a (mod q) is in balanced form if it represents an integer
  //       b in [-q/2, q/2) such that b == a (mod q). In particular, b may
  //       take negative values.
  static absl::StatusOr<RnsPolynomial> ConvertBalancedFromPolynomialCoeffs(
      const std::vector<ModularInt>& coeffs_q,  // coefficients (mod q)
      const ModularIntParams* mod_params_q,
      absl::Span<const PrimeModulus<ModularInt>* const> moduli);

  // Returns a RnsPolynomial whose coefficients are given in `coeffs_q` modulo
  // a prime modulus q. The returned RNS polynomial is in NTT form wrt `moduli`.
  // Note: We do not consider balanced representation in this function, and
  //       in particular a modular value a (mod q) represents an integer b in
  //       [0, q) such that b == a (mod q).
  static absl::StatusOr<RnsPolynomial> ConvertFromPolynomialCoeffs(
      const std::vector<ModularInt>& coeffs_q,  // in NTT form, (mod q)
      const ModularIntParams* mod_params_q,
      absl::Span<const PrimeModulus<ModularInt>* const> moduli);

  static absl::StatusOr<RnsPolynomial> Deserialize(
      const SerializedRnsPolynomial& serialized,
      absl::Span<const PrimeModulus<ModularInt>* const> moduli) {
    int log_n = serialized.log_n();
    if (log_n <= 0) {
      return absl::InvalidArgumentError("`log_n` must be positive.");
    }
    int num_coeff_vectors = serialized.coeff_vectors_size();
    if (num_coeff_vectors != moduli.size()) {
      return absl::InvalidArgumentError(absl::StrCat(
          "Number of serialized coefficient vectors must be ", moduli.size()));
    }
    int num_coeffs = 1 << log_n;
    std::vector<std::vector<ModularInt>> coeff_vectors;
    coeff_vectors.reserve(num_coeff_vectors);
    for (int i = 0; i < num_coeff_vectors; ++i) {
      RLWE_ASSIGN_OR_RETURN(
          std::vector<ModularInt> coeffs_qi,
          ModularInt::DeserializeVector(num_coeffs, serialized.coeff_vectors(i),
                                        moduli[i]->ModParams()));
      coeff_vectors.push_back(std::move(coeffs_qi));
    }
    return Create(std::move(coeff_vectors), serialized.is_ntt());
  }

  absl::StatusOr<SerializedRnsPolynomial> Serialize(
      absl::Span<const PrimeModulus<ModularInt>* const> moduli) const {
    if (moduli.size() != coeff_vectors_.size()) {
      return absl::InvalidArgumentError(absl::StrCat(
          "`moduli` must contain ", coeff_vectors_.size(), " RNS moduli."));
    }

    SerializedRnsPolynomial output;
    output.set_log_n(log_n_);
    output.set_is_ntt(is_ntt_);
    for (int i = 0; i < coeff_vectors_.size(); ++i) {
      RLWE_ASSIGN_OR_RETURN(*(output.add_coeff_vectors()),
                            ModularInt::SerializeVector(
                                coeff_vectors_[i], moduli[i]->ModParams()));
    }
    return output;
  }

  // Converts this RNS polynomial from Coefficient form to NTT form.
  absl::Status ConvertToNttForm(
      absl::Span<const PrimeModulus<ModularInt>* const> moduli);

  // Converts this RNS polynomial from NTT form to Coefficient form.
  absl::Status ConvertToCoeffForm(
      absl::Span<const PrimeModulus<ModularInt>* const> moduli);

  // Returns the negation of this polynomial.
  absl::StatusOr<RnsPolynomial> Negate(
      absl::Span<const PrimeModulus<ModularInt>* const> moduli) const {
    RnsPolynomial output = *this;
    output.NegateInPlace(moduli);
    return output;
  }

  // Negates this polynomial in-place and returns the result.
  absl::Status NegateInPlace(
      absl::Span<const PrimeModulus<ModularInt>* const> moduli) {
    int num_moduli = coeff_vectors_.size();
    if (moduli.size() != num_moduli) {
      return absl::InvalidArgumentError(
          absl::StrCat("`moduli` must contain ", num_moduli, " RNS moduli."));
    }
    for (int i = 0; i < num_moduli; ++i) {
      const auto& mod_params = moduli[i]->ModParams();
      for (ModularInt& coeff : coeff_vectors_[i]) {
        coeff.NegateInPlace(mod_params);
      }
    }
    return absl::OkStatus();
  }

  // Returns the polynomial `this` + `that`.
  absl::StatusOr<RnsPolynomial> Add(
      const RnsPolynomial& that,
      absl::Span<const PrimeModulus<ModularInt>* const> moduli) const {
    RnsPolynomial output = *this;
    RLWE_RETURN_IF_ERROR(output.AddInPlace(that, moduli));
    return output;
  }

  // Returns the polynomial `this` - `that`.
  absl::StatusOr<RnsPolynomial> Sub(
      const RnsPolynomial& that,
      absl::Span<const PrimeModulus<ModularInt>* const> moduli) const {
    RnsPolynomial output = *this;
    RLWE_RETURN_IF_ERROR(output.SubInPlace(that, moduli));
    return output;
  }

  // Adds `that` to `this` in-place.
  absl::Status AddInPlace(
      const RnsPolynomial& that,
      absl::Span<const PrimeModulus<ModularInt>* const> moduli) {
    int num_moduli = coeff_vectors_.size();
    if (that.coeff_vectors_.size() != num_moduli) {
      return absl::InvalidArgumentError(absl::StrCat(
          "`that` must contain ", num_moduli, " coefficient vectors."));
    }
    if (moduli.size() != num_moduli) {
      return absl::InvalidArgumentError(
          absl::StrCat("`moduli` must contain ", num_moduli, " RNS moduli."));
    }
    for (int i = 0; i < num_moduli; ++i) {
      RLWE_RETURN_IF_ERROR(ModularInt::BatchAddInPlace(
          &coeff_vectors_[i], that.coeff_vectors_[i], moduli[i]->ModParams()));
    }
    return absl::OkStatus();
  }

  // Substracts `that` from `this` in-place.
  absl::Status SubInPlace(
      const RnsPolynomial& that,
      absl::Span<const PrimeModulus<ModularInt>* const> moduli) {
    int num_moduli = coeff_vectors_.size();
    if (that.coeff_vectors_.size() != num_moduli) {
      return absl::InvalidArgumentError(absl::StrCat(
          "`that` must contain ", num_moduli, " coefficient vectors."));
    }
    if (moduli.size() != num_moduli) {
      return absl::InvalidArgumentError(
          absl::StrCat("`moduli` must contain ", num_moduli, " RNS moduli."));
    }
    for (int i = 0; i < num_moduli; ++i) {
      RLWE_RETURN_IF_ERROR(ModularInt::BatchSubInPlace(
          &coeff_vectors_[i], that.coeff_vectors_[i], moduli[i]->ModParams()));
    }
    return absl::OkStatus();
  }

  // Returns the product of this polynomial and an integer scalar.
  absl::StatusOr<RnsPolynomial> Mul(
      typename ModularInt::Int scalar,
      absl::Span<const PrimeModulus<ModularInt>* const> moduli) const {
    RnsPolynomial output = *this;
    RLWE_RETURN_IF_ERROR(output.MulInPlace(scalar, moduli));
    return output;
  }

  // Multiplies this polynomial by an integer scalar.
  absl::Status MulInPlace(
      typename ModularInt::Int scalar,
      absl::Span<const PrimeModulus<ModularInt>* const> moduli) {
    int num_moduli = coeff_vectors_.size();
    if (moduli.size() != num_moduli) {
      return absl::InvalidArgumentError(
          absl::StrCat("`moduli` must contain ", num_moduli, " RNS moduli."));
    }
    for (int i = 0; i < num_moduli; ++i) {
      const ModularIntParams* mod_params_qi = moduli[i]->ModParams();
      RLWE_ASSIGN_OR_RETURN(auto scalar_mod_qi,
                            ModularInt::ImportInt(scalar, mod_params_qi));
      RLWE_RETURN_IF_ERROR(ModularInt::BatchMulInPlace(
          &coeff_vectors_[i], scalar_mod_qi, mod_params_qi));
    }
    return absl::OkStatus();
  }

  // Multiplies this polynomial by a scalar that is given in RNS form.
  absl::Status MulInPlace(
      const RnsInt<ModularInt>& scalar,
      absl::Span<const PrimeModulus<ModularInt>* const> moduli) {
    return MulInPlace(scalar.RnsRep(), moduli);
  }

  // Multiplies this polynomial by a scalar that is given in RNS form.
  absl::Status MulInPlace(
      absl::Span<const ModularInt> scalar_mod_qs,
      absl::Span<const PrimeModulus<ModularInt>* const> moduli) {
    int num_moduli = coeff_vectors_.size();
    if (scalar_mod_qs.size() != num_moduli) {
      return absl::InvalidArgumentError(absl::StrCat(
          "`scalar_mod_qs` must contain ", num_moduli, " modular integers."));
    }
    if (moduli.size() != num_moduli) {
      return absl::InvalidArgumentError(
          absl::StrCat("`moduli` must contain ", num_moduli, " RNS moduli."));
    }
    for (int i = 0; i < moduli.size(); ++i) {
      RLWE_RETURN_IF_ERROR(ModularInt::BatchMulInPlace(
          &coeff_vectors_[i], scalar_mod_qs[i], moduli[i]->ModParams()));
    }
    return absl::OkStatus();
  }

  // Returns the product of this polynomial and `that`.
  // Both `this` and `that` must be in NTT form.
  absl::StatusOr<RnsPolynomial> Mul(
      const RnsPolynomial& that,
      absl::Span<const PrimeModulus<ModularInt>* const> moduli) const {
    RnsPolynomial output = *this;
    RLWE_RETURN_IF_ERROR(output.MulInPlace(that, moduli));
    return output;
  }

  // Multiplies this polynomial by `that`.
  // Both `this` and `that` must be in NTT form.
  absl::Status MulInPlace(
      const RnsPolynomial& that,
      absl::Span<const PrimeModulus<ModularInt>* const> moduli) {
    int num_moduli = coeff_vectors_.size();
    if (that.coeff_vectors_.size() != num_moduli) {
      return absl::InvalidArgumentError(absl::StrCat(
          "`that` must contain ", num_moduli, " coefficient vectors."));
    }
    if (moduli.size() != num_moduli) {
      return absl::InvalidArgumentError(
          absl::StrCat("`moduli` must contain ", num_moduli, " RNS moduli."));
    }
    if (!is_ntt_) {
      return absl::InvalidArgumentError(
          "RNS polynomial `this` must be in NTT form.");
    }
    if (!that.is_ntt_) {
      return absl::InvalidArgumentError(
          "RNS polynomial `that` must be in NTT form.");
    }

    for (int i = 0; i < num_moduli; ++i) {
      RLWE_RETURN_IF_ERROR(ModularInt::BatchMulInPlace(
          &coeff_vectors_[i], that.coeff_vectors_[i], moduli[i]->ModParams()));
    }
    return absl::OkStatus();
  }

  // Adds the polynomial product a * b to this polynomial.
  // Polynomials `this`, `a`, and `b` must be all in NTT form.
  absl::Status FusedMulAddInPlace(
      const RnsPolynomial& a, const RnsPolynomial& b,
      absl::Span<const PrimeModulus<ModularInt>* const> moduli) {
    int num_moduli = coeff_vectors_.size();
    if (a.coeff_vectors_.size() != num_moduli) {
      return absl::InvalidArgumentError(absl::StrCat(
          "`a` must contain ", num_moduli, " coefficient vectors."));
    }
    if (b.coeff_vectors_.size() != num_moduli) {
      return absl::InvalidArgumentError(absl::StrCat(
          "`b` must contain ", num_moduli, " coefficient vectors."));
    }
    if (moduli.size() != num_moduli) {
      return absl::InvalidArgumentError(
          absl::StrCat("`moduli` must contain ", num_moduli, " RNS moduli."));
    }
    if (!is_ntt_) {
      return absl::InvalidArgumentError(
          "RNS polynomial `this` must be in NTT form.");
    }
    if (!a.is_ntt_) {
      return absl::InvalidArgumentError(
          "RNS polynomial `a` must be in NTT form.");
    }
    if (!b.is_ntt_) {
      return absl::InvalidArgumentError(
          "RNS polynomial `b` must be in NTT form.");
    }

    for (int i = 0; i < num_moduli; ++i) {
      RLWE_RETURN_IF_ERROR(ModularInt::BatchFusedMulAddInPlace(
          &coeff_vectors_[i], a.coeff_vectors_[i], b.coeff_vectors_[i],
          moduli[i]->ModParams()));
    }
    return absl::OkStatus();
  }

  // Returns the polynomial a(X^power) by substituting X with X^power in this
  // polynomial a.
  // Note: 1) `power` must be an odd non-negative integer less than 2*N;
  //       2) this polynomial must be in NTT form.
  absl::StatusOr<RnsPolynomial> Substitute(
      int power,
      absl::Span<const PrimeModulus<ModularInt>* const> moduli) const;

  // Converts this polynomial from mod Q to mod (Q/q_l), where Q = q_0 ... q_l.
  // This version of modulus reduction is compatible with BGV scheme where the
  // plaintext polynomial is embedded at LSB of ciphertext modulus, and modulus
  // reduction makes sure that the plaintext is kept at the LSB after reducing
  // to modulus Q/q_l, i.e. the resulting polynomial a' (mod Q/q_l) is such that
  // a' % t == a % t for plaintext modulus t.
  //
  // `ql_inv` stores the RNS integer q_l (mod Q/q_l) wrt reduced RNS basis
  // (q_0, .., q_{l-1}), and `moduli` stores the RNS prime moduli q_0, .., q_l
  // for this RNS polynomial.
  absl::Status ModReduceLsb(
      typename ModularInt::Int t, const RnsInt<ModularInt>& ql_inv,
      absl::Span<const PrimeModulus<ModularInt>* const> moduli);

  // Converts this polynomial from mod Q to mod (Q/q_l), where Q = q_0 ... q_l
  // as specified in `moduli`.
  // This version of modulus reduction updates a RNS polynomial a (mod Q) to
  // round(a / q_l (mod Q/q_l)), and it is compatible with BFV and CKKS schemes
  // where the plaintext polynomial resides at MSB of ciphertext modulus. Still,
  // modulus reduction keeps the plaintext at the MSB of the reduced modulus
  // Q/q_l, and it only introduces error in the low order bits.
  absl::Status ModReduceMsb(
      const RnsInt<ModularInt>& ql_inv,
      absl::Span<const PrimeModulus<ModularInt>* const> moduli);

  // Returns `this` (mod P), where this polynomial is modulo Q = q_0 ... q_L.
  // This polynomial must be in coefficient form, and the returned polynomial is
  // in coefficient form too.
  // The parameters are:
  // `this_moduli` specifies the prime moduli q_i,
  // `output_moduli` specifies the prime moduli p_j of P,
  // `prime_q_hat_inv_mod_qs` = {(Q/q_i)^(-1) mod q_i}_i,
  // `prime_q_hat_mod_ps` = {Q / q_i mod p_j}_j,
  // `q_mod_ps` = {Q mod p_j}_j.
  absl::StatusOr<RnsPolynomial<ModularInt>> SwitchRnsBasis(
      absl::Span<const PrimeModulus<ModularInt>* const> this_moduli,    // q_i's
      absl::Span<const PrimeModulus<ModularInt>* const> output_moduli,  // p_j's
      absl::Span<const ModularInt> prime_q_hat_inv_mod_qs,
      absl::Span<const RnsInt<ModularInt>> prime_q_hat_mod_ps,
      absl::Span<const ModularInt> q_mod_ps) const;

  // Returns round(P / Q * `this`) mod P, where this polynomial is modulo Q =
  // q_0 * ... * q_L, and the output RNS modulus is P = p_1 * ... * p_K.
  // This polynomial must be in coefficient form, and the returned polynomial is
  // in coefficient form too.
  // The parameters are:
  // `this_moduli` specifies the prime moduli q_i,
  // `output_moduli` specifies the prime moduli p_j of P,
  // `prime_q_hat_inv_mod_qs` = {(Q/q_i)^(-1) mod q_i}_i,
  // `prime_q_inv_mod_ps` = {q_i^(-1) mod p_j}_j,
  // `p_mod_qs` = {P mod q_i}_i.
  absl::StatusOr<RnsPolynomial<ModularInt>> ScaleAndSwitchRnsBasis(
      absl::Span<const PrimeModulus<ModularInt>* const> this_moduli,    // q_i's
      absl::Span<const PrimeModulus<ModularInt>* const> output_moduli,  // p_j's
      absl::Span<const ModularInt> prime_q_hat_inv_mod_qs,
      absl::Span<const RnsInt<ModularInt>> prime_q_inv_mod_ps,
      absl::Span<const ModularInt> p_mod_qs) const;

  // Extends this polynomial from mod Q to mod (Q*P), where Q = q_0 ... q_l
  // is specified in `this_moduli`, and P = p_1 ... p_k is specified in
  // `aux_moduli`.
  absl::Status ExtendRnsBasisInPlace(
      absl::Span<const PrimeModulus<ModularInt>* const> this_moduli,  // q_i's
      absl::Span<const PrimeModulus<ModularInt>* const> aux_moduli,   // p_j's
      absl::Span<const ModularInt> prime_q_hat_inv_mod_qs,
      absl::Span<const RnsInt<ModularInt>> prime_q_hat_mod_ps,
      absl::Span<const ModularInt> q_mod_ps) {
    RLWE_ASSIGN_OR_RETURN(
        RnsPolynomial<ModularInt> a_mod_aux,
        SwitchRnsBasis(this_moduli, aux_moduli, prime_q_hat_inv_mod_qs,
                       prime_q_hat_mod_ps, q_mod_ps));
    coeff_vectors_.reserve(coeff_vectors_.size() +
                           a_mod_aux.coeff_vectors_.size());
    for (auto& coeffs_mod_pj : a_mod_aux.coeff_vectors_) {
      coeff_vectors_.push_back(std::move(coeffs_mod_pj));
    }
    return absl::OkStatus();
  }

  // Scales the polynomial (this, polynomial_aux) mod (Q, P) up by `t` and then
  // scales down by P, storing the resulting polynomial (mod Q) in `this`.
  // This function computes round(t / P * polynomial) mod Q, where polynomial is
  // modulo Q * P, represented by `this` (mod Q) and `polynomial_aux` (mod P).
  absl::Status ScaleAndReduceRnsBasisInPlace(
      RnsPolynomial<ModularInt> polynomial_aux, typename ModularInt::Int t,
      absl::Span<const PrimeModulus<ModularInt>* const> this_moduli,  // q_i's
      absl::Span<const PrimeModulus<ModularInt>* const> aux_moduli,   // p_j's
      absl::Span<const ModularInt> prime_p_hat_inv_mod_ps,
      absl::Span<const RnsInt<ModularInt>> prime_p_hat_mod_qs,
      absl::Span<const ModularInt> p_mod_qs,
      absl::Span<const ModularInt> p_inv_mod_qs) {
    // t * `this` (mod Q).
    RLWE_RETURN_IF_ERROR(MulInPlace(t, this_moduli));
    // t * `polynomial_aux` (mod P).
    RLWE_RETURN_IF_ERROR(polynomial_aux.MulInPlace(t, aux_moduli));

    // Next we scale down the polynomial (this, polynomial_aux) mod (Q*P) by P,
    // and take the result mod Q, which can be computed as
    // (polynomial - polynomial_aux) / P (mod Q).
    if (polynomial_aux.IsNttForm()) {
      RLWE_RETURN_IF_ERROR(polynomial_aux.ConvertToCoeffForm(aux_moduli));
    }
    // Compute the delta term = `polynomial_aux` (mod Q).
    RLWE_ASSIGN_OR_RETURN(RnsPolynomial<ModularInt> delta,
                          polynomial_aux.SwitchRnsBasis(
                              aux_moduli, this_moduli, prime_p_hat_inv_mod_ps,
                              prime_p_hat_mod_qs, p_mod_qs));
    if (IsNttForm()) {
      RLWE_RETURN_IF_ERROR(delta.ConvertToNttForm(this_moduli));
    }
    // (polynomial - polynomial_aux) (mod Q).
    RLWE_RETURN_IF_ERROR(SubInPlace(delta, this_moduli));
    // (polynomial - polynomial_aux) / P (mod Q).
    RLWE_RETURN_IF_ERROR(MulInPlace(p_inv_mod_qs, this_moduli));
    return absl::OkStatus();
  }

  // Returns the polynomial b(X) wrt RNS modulus P = [p_0, ..., p_{K-1}] that
  // approximately represents the same polynomial over the integers as this
  // polynomial mod Q. The returned polynomial b(X) is in coefficient form.
  //
  // This is an implementation of the fast approximate basis conversion
  // algorithm from
  // ``A Full RNS Variant of FV like Somewhat Homomorphic Encryption Schemes''
  // by Bajard et al, https://eprint.iacr.org/2016/510
  // Note that the conversion introduces an approximation error, namely for
  // a (mod Q) under the RNS basis Q = [q_0, ..., q_L], then the conversion
  // outputs b = (a + u * Q) mod P, and hence the error is (u * Q) mod P, where
  // u \in [0, L). So, for P larger than Q, the result is approximately correct
  // and can be fully corrected by taking modulo Q.
  // The parameter `is_balanced_rep` indicates whether mod-q values represent
  // integers in [-q/2, q/2), which is useful when plaintext values reside in
  // lower order bits of modulus.
  absl::StatusOr<RnsPolynomial<ModularInt>> ApproxSwitchRnsBasis(
      absl::Span<const PrimeModulus<ModularInt>* const> this_moduli,    // q_i's
      absl::Span<const PrimeModulus<ModularInt>* const> output_moduli,  // p_j's
      absl::Span<const ModularInt> prime_q_hat_inv_mod_qs,
      absl::Span<const RnsInt<ModularInt>> prime_q_hat_mod_ps,
      bool is_balanced_rep = false) const;

  // Approximately scales down the polynomial c(X) (mod Q*P), represented by
  // (this (mod Q), polynomial_aux (mod P)) wrt extended RNS basis Q \cup P, to
  // mod-Q and updates this polynomial using the result in NTT form.
  // The composite modulus Q is the product of prime moduli in `this_moduli`,
  // and P is the product of prime moduli in `aux_moduli`.
  //
  // This version of approximate modulus reduction is compatible with BGV where
  // plaintext resides in the LSB of modulus Q and error is scaled by plaintext
  // modulus t. It essentially computes (r * this + delta)/P (mod Q), where
  // P = k * t + r for residue r = P (mod t), and delta is approximately
  // converted from t * [k * polynomial_aux (mod P)] to mod Q.
  // If c = m + t * e (mod Q*P) with this = c mod Q and polynomial_aux = c mod P
  // == this + u * Q (mod P) for small u, then the approximate reduction will
  // produce m + t * round(e' * r / P) (mod Q). So the error size is reduced by
  // a factor roughly P / r.
  absl::Status ApproxModReduceLsb(
      RnsPolynomial<ModularInt> polynomial_aux,
      absl::Span<const PrimeModulus<ModularInt>* const> this_moduli,  // q_i's
      absl::Span<const PrimeModulus<ModularInt>* const> aux_moduli,   // p_j's
      absl::Span<const ModularInt> p_inv_mod_qs,
      absl::Span<const ModularInt> prime_p_hat_inv_mod_ps,
      absl::Span<const RnsInt<ModularInt>> prime_p_hat_mod_qs,
      absl::Span<const ModularInt> t_inv_mod_ps,
      typename ModularInt::Int p_mod_t, typename ModularInt::Int t);

  // Approximately scales down the polynomial c(X) (mod Q*P), represented by
  // (this (mod Q), polynomial_aux (mod P)) wrt extended RNS basis Q \cup P, to
  // mod-Q and updates this polynomial with the result in NTT form.
  // The composite modulus Q is the product of prime moduli in `this_moduli` and
  // P is the product of prime moduli in `aux_moduli`.
  //
  // This version of approximate modulus reduction is compatible with BFV and
  // CKKS schemes where plaintext resides in the MSB of modulus Q. It
  // approximately computes round(c / P) (mod Q), by computing (this - delta) /
  // P (mod Q), where delta is approximately converted from polynomial_aux from
  // mod-P to mod-Q. If c = e + P * T * m (mod Q*P) for some scaling factor T,
  // such that this = c mod Q and polynomial_aux = c mod P == this + u * Q (mod
  // P) for small u, then the result will be approximately round(c / P) (mod Q),
  // i.e. e' + T * m (mod Q) with e' roughly equal to e / P.
  absl::Status ApproxModReduceMsb(
      RnsPolynomial<ModularInt> polynomial_aux,
      absl::Span<const PrimeModulus<ModularInt>* const> this_moduli,  // q_i's
      absl::Span<const PrimeModulus<ModularInt>* const> aux_moduli,   // p_j's
      absl::Span<const ModularInt> p_inv_mod_qs,
      absl::Span<const ModularInt> prime_p_hat_inv_mod_ps,
      absl::Span<const RnsInt<ModularInt>> prime_p_hat_mod_qs);

  // Appends the coefficient vector of a new sub-polynomial to the RNS
  // coefficient vectors of this polynomial.
  void AppendCoeffVector(std::vector<ModularInt> coeffs) {
    coeff_vectors_.push_back(std::move(coeffs));
  }

  // Removes and returns the coefficients of the last sub-polynomial.
  absl::StatusOr<std::vector<ModularInt>> DetachLastCoeffVector() {
    if (coeff_vectors_.empty()) {
      return absl::InvalidArgumentError("RnsPolynomial is empty.");
    }
    std::vector<ModularInt> last_coeffs = coeff_vectors_.back();
    coeff_vectors_.pop_back();
    return last_coeffs;
  }

  // Replaces a sub-polynomial at the given index.
  absl::Status ReplaceSubPolynomialAt(int index,
                                      const Polynomial<ModularInt>& a_qi) {
    if (index < 0 || index >= coeff_vectors_.size()) {
      return absl::InvalidArgumentError(
          absl::StrCat("Index ", index, " is out of bound."));
    }
    coeff_vectors_[index] = a_qi.Coeffs();
    return absl::OkStatus();
  }

  // Comparison with another RnsPolynomial.
  bool operator==(const RnsPolynomial& that) const {
    return std::tie(log_n_, is_ntt_, coeff_vectors_) ==
           std::tie(that.log_n_, that.is_ntt_, that.coeff_vectors_);
  }
  bool operator!=(const RnsPolynomial& that) const { return !(*this == that); }

  // Accessors.
  // Log of the polynomial ring dimension N.
  int LogN() const { return log_n_; }

  // Number of coefficients in the polynomial (mod X^N+1).
  int NumCoeffs() const { return 1 << log_n_; }

  // Number of prime moduli comprising the RNS modulus chain.
  int NumModuli() const { return coeff_vectors_.size(); }

  // Returns the entire collection of sub-polynomial coefficients.
  const std::vector<std::vector<ModularInt>>& Coeffs() const {
    return coeff_vectors_;
  }

  // Returns whether this polynomial is in the NTT form.
  bool IsNttForm() const { return is_ntt_; }

  // Clones this polynomial.
  RnsPolynomial Clone() const {
    return RnsPolynomial(log_n_, coeff_vectors_, is_ntt_);
  }

 private:
  explicit RnsPolynomial(int log_n,
                         std::vector<std::vector<ModularInt>> coeff_vectors,
                         bool is_ntt)
      : log_n_(log_n),
        coeff_vectors_(std::move(coeff_vectors)),
        is_ntt_(is_ntt) {}

  // Log of ring degree N.
  int log_n_;

  // Coefficients of the polynomial modulo prime moduli.
  // Each vector corresponds to a prime modulus in moduli_ in the same order.
  std::vector<std::vector<ModularInt>> coeff_vectors_;

  // Is the polynomial in the NTT form?
  bool is_ntt_;
};

template <typename ModularInt>
template <typename Prng>
absl::StatusOr<RnsPolynomial<ModularInt>>
RnsPolynomial<ModularInt>::SampleUniform(
    int log_n, Prng* prng,
    absl::Span<const PrimeModulus<ModularInt>* const> moduli) {
  if (log_n <= 0) {
    return absl::InvalidArgumentError("`log_n` must be positive.");
  }
  if (prng == nullptr) {
    return absl::InvalidArgumentError("`prng` must not be null.");
  }
  if (moduli.empty()) {
    return absl::InvalidArgumentError("`moduli` must not be empty.");
  }

  int num_coeffs = 1 << log_n;
  std::vector<std::vector<ModularInt>> coeff_vectors(moduli.size());
  for (int i = 0; i < moduli.size(); ++i) {
    coeff_vectors[i].reserve(num_coeffs);
    const typename ModularInt::Params* mod_params_qi = moduli[i]->ModParams();
    for (int j = 0; j < num_coeffs; ++j) {
      RLWE_ASSIGN_OR_RETURN(ModularInt coeff,
                            ModularInt::ImportRandom(prng, mod_params_qi));
      coeff_vectors[i].push_back(std::move(coeff));
    }
  }
  return RnsPolynomial<ModularInt>(log_n, std::move(coeff_vectors),
                                   /*is_ntt=*/true);
}

}  // namespace rlwe

#endif  // RLWE_RNS_RNS_POLYNOMIAL_H_
