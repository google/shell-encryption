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

#ifndef RLWE_RNS_RNS_GADGET_H_
#define RLWE_RNS_RNS_GADGET_H_

#include <algorithm>
#include <utility>
#include <vector>

#include "absl/algorithm/container.h"
#include "absl/status/status.h"
#include "absl/types/span.h"
#include "shell_encryption/rns/rns_integer.h"
#include "shell_encryption/rns/rns_modulus.h"
#include "shell_encryption/rns/rns_polynomial.h"
#include "shell_encryption/status_macros.h"

namespace rlwe {

// The Gadget class represents a gadget vector g with respect to the RNS modulus
// Q = q_0 * ... * q_L, for representing mod-Q numbers using "digits" of quality
// bounded by certain bases B = {b_0, ..., b_L}. The basic operation supported
// by a gadget is decomposition, g^-1(x mod Q) = (x_0 mod q_0, ..., x_L mod q_L)
// such that the inner product between the decomposition and the gadaget vector
// satisfies that <g^-1(x mod Q), g> mod Q == x. Each component x_i of the
// decomposition is a vector whose entries are bounded by b_i. The length of the
// decomposition which is also the dimension of the gadget vector g, is the sum
// of the lengths of x_i. The implementation assumes all b_i are powers of 2.
// Note that a gadget matrix can be efficiently computed as I \tensor g.
//
// This class implements the so-called digit decomposition gadget with CRT
// representation. A "gadget" in lattice-based cryptography in general is a way
// to encode an element with large norm into a higher-dimension vector of
// "smaller" elements, which has benefits for controlling noise growth. For a
// formal exposition see https://eprint.iacr.org/2018/946.
//
template <typename ModularInt>
class RnsGadget {
  using Integer = typename ModularInt::Int;
  using ModularIntParams = typename ModularInt::Params;

 public:
  // Create a gadget with CRT representation given the log of bases b_i and the
  // RNS moduli. Note that a gadget is the concatenation of gadgets in the
  // residual ring Z_{qi}:
  //   g_crt = [ qi_hat * qi_hat_inv * g_i for i in 1..L ] mod Q,
  // where Q = q1 * ... * qL, qi_hat = Q / q_i, qi_hat_inv = qi_hat^(-1) mod qi.
  static absl::StatusOr<RnsGadget> Create(
      int log_n, std::vector<size_t> log_bs,
      absl::Span<const ModularInt> q_hats,      // {qi_hat mod qi}_i
      absl::Span<const ModularInt> q_hat_invs,  // {qi_hat_inv}_i
      absl::Span<const PrimeModulus<ModularInt>* const> moduli);

  // Returns the gadget decomposition of a RNS polynomial.
  // Note that the input polynomial `a` must be in coefficient form, and the
  // returned polynomial is in the coefficient form.
  absl::StatusOr<std::vector<RnsPolynomial<ModularInt>>> Decompose(
      const RnsPolynomial<ModularInt>& a,
      absl::Span<const PrimeModulus<ModularInt>* const> moduli) const;

  // Returns the gadget decomposition of a vector of RNS polynomials.
  // The result is a vector containing decompositions of every polynomial in
  // the input vector.
  // Note that the input polynomials in `as` must all be in coefficient form.
  absl::StatusOr<std::vector<std::vector<RnsPolynomial<ModularInt>>>> Decompose(
      absl::Span<const RnsPolynomial<ModularInt>> as,
      absl::Span<const PrimeModulus<ModularInt>* const> moduli) const {
    std::vector<std::vector<RnsPolynomial<ModularInt>>> as_digits;
    as_digits.reserve(as.size());
    for (const auto& a : as) {
      RLWE_ASSIGN_OR_RETURN(std::vector<RnsPolynomial<ModularInt>> a_digits,
                            Decompose(a, moduli));
      as_digits.push_back(std::move(a_digits));
    }
    return as_digits;
  }

  int LogGadgetBase() const {
    auto max = std::max_element(log_bs_.begin(), log_bs_.end());
    return *max;
  }

  // The dimension of the gadget vector.
  int Dimension() const { return absl::c_accumulate(dims_, size_t{0}); }

  const RnsInt<ModularInt>& Component(int index) const { return gs_[index]; }

 private:
  explicit RnsGadget(std::vector<size_t> log_bs, std::vector<size_t> dimensions,
                     std::vector<RnsInt<ModularInt>> gs)
      : log_bs_(std::move(log_bs)),
        dims_(std::move(dimensions)),
        gs_(std::move(gs)) {}

  // The base must be small enough to be represented by the residual int.
  std::vector<size_t> log_bs_;

  // The sub-dimensions, which are the lengths of components of a gadget vector
  // wrt RNS prime moduli.
  std::vector<size_t> dims_;

  // The gadget vector {g_i mod q_i}_i, stored as RNS integers.
  std::vector<RnsInt<ModularInt>> gs_;
};

}  // namespace rlwe

#endif  // RLWE_RNS_RNS_GADGET_H_
