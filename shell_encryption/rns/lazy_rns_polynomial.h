/*
 * Copyright 2024 Google LLC.
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

#ifndef RLWE_RNS_LAZY_RNS_POLYNOMIAL_H_
#define RLWE_RNS_LAZY_RNS_POLYNOMIAL_H_

#include <utility>
#include <vector>

#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "absl/strings/str_cat.h"
#include "absl/types/span.h"
#include "hwy/aligned_allocator.h"
#include "shell_encryption/rns/rns_modulus.h"
#include "shell_encryption/rns/rns_polynomial.h"
#include "shell_encryption/status_macros.h"

namespace rlwe {

// This class defines a `lazy` polynomial in the RNS system, which stores its
// RNS coefficients as BigInt, and only performs modulo reduction when they
// reach the maximal value a BigInt can hold.
//
// Currently LazyRnsPolynomial can only be operated on by the fused-multiply-add
// operation.
template <typename ModularInt>
class LazyRnsPolynomial {
 public:
  using Integer = typename ModularInt::Int;
  using BigInt = typename ModularInt::BigInt;

  // Returns a lazy polynomial representing the product of `a` and `b`.
  static absl::StatusOr<LazyRnsPolynomial> Create(
      const RnsPolynomial<ModularInt>& a, const RnsPolynomial<ModularInt>& b,
      absl::Span<const PrimeModulus<ModularInt>* const> moduli) {
    if (!a.IsNttForm() || !b.IsNttForm()) {
      return absl::InvalidArgumentError(
          "Polynomials `a` and `b` must be in the NTT form.");
    }
    int num_coeffs = a.NumCoeffs();
    if (b.NumCoeffs() != num_coeffs) {
      return absl::InvalidArgumentError(
          "Polynomials `a` and `b` must have the same number of coefficients.");
    }
    int num_moduli = moduli.size();
    if (a.NumModuli() != num_moduli || b.NumModuli() != num_moduli) {
      return absl::InvalidArgumentError(
          "Polynomials `a` and `b` must be defined wrt `moduli`");
    }

    // Get the maximal number of FMA operations can be made wrt the RNS moduli.
    RLWE_ASSIGN_OR_RETURN(int log_maximum_level,
                          ComputeLogMaximalLevel(moduli));

    // Let's multiply the CRT coefficients.
    std::vector<hwy::AlignedVector<BigInt>> coeff_vectors(num_moduli);
    for (int i = 0; i < num_moduli; ++i) {
      coeff_vectors[i].reserve(num_coeffs);
      const auto& a_coeffs = a.Coeffs()[i];
      const auto& b_coeffs = b.Coeffs()[i];
      for (int j = 0; j < num_coeffs; ++j) {
        coeff_vectors[i].push_back(
            static_cast<BigInt>(a_coeffs[j].GetMontgomeryRepresentation()) *
            b_coeffs[j].GetMontgomeryRepresentation());
      }
    }

    return LazyRnsPolynomial(std::move(coeff_vectors),
                             /*current_level=*/static_cast<BigInt>(1),
                             static_cast<BigInt>(1) << log_maximum_level);
  }

  // Returns a lazy polynomial representing the zero polynomial.
  static absl::StatusOr<LazyRnsPolynomial> CreateZero(
      int log_n, absl::Span<const PrimeModulus<ModularInt>* const> moduli) {
    if (log_n <= 0) {
      return absl::InvalidArgumentError("`log_n` must be positive.");
    }
    if (moduli.empty()) {
      return absl::InvalidArgumentError("`moduli` must not be empty.");
    }
    int num_moduli = moduli.size();
    int num_coeffs = 1 << log_n;
    RLWE_ASSIGN_OR_RETURN(int log_maximum_level,
                          ComputeLogMaximalLevel(moduli));
    std::vector<hwy::AlignedVector<BigInt>> coeff_vectors(num_moduli);
    for (int i = 0; i < num_moduli; ++i) {
      coeff_vectors[i].resize(num_coeffs, 0);
    }
    return LazyRnsPolynomial(std::move(coeff_vectors),
                             /*current_level=*/static_cast<BigInt>(0),
                             static_cast<BigInt>(1) << log_maximum_level);
  }

  // Returns this polynomial reduced wrt `moduli`, as a RnsPolynomial in NTT
  // form.
  absl::StatusOr<RnsPolynomial<ModularInt>> Export(
      absl::Span<const PrimeModulus<ModularInt>* const> moduli) {
    int num_moduli = moduli.size();
    if (num_moduli != coeff_vectors_.size()) {
      return absl::InvalidArgumentError(
          "`moduli` does not contain enough RNS moduli.");
    }
    int num_coeffs = coeff_vectors_[0].size();
    std::vector<std::vector<ModularInt>> mod_coeff_vectors(num_moduli);
    for (int i = 0; i < num_moduli; ++i) {
      mod_coeff_vectors.reserve(num_coeffs);
      auto mod_params_qi = moduli[i]->ModParams();
      for (auto const& coeff : coeff_vectors_[i]) {
        mod_coeff_vectors[i].emplace_back(mod_params_qi->ExportInt(
            mod_params_qi->BarrettReduceBigInt(coeff)));
      }
    }
    return RnsPolynomial<ModularInt>::Create(std::move(mod_coeff_vectors),
                                             /*is_ntt=*/true);
  }

  // Adds the product `a`* `b` (mod `moduli`) to this polynomial.
  absl::Status FusedMulAddInPlace(
      const RnsPolynomial<ModularInt>& a, const RnsPolynomial<ModularInt>& b,
      absl::Span<const PrimeModulus<ModularInt>* const> moduli);

 private:
  explicit LazyRnsPolynomial(
      std::vector<hwy::AlignedVector<BigInt>> coeff_vectors,
      BigInt current_level, BigInt maximum_level)
      : coeff_vectors_(std::move(coeff_vectors)),
        current_level_(current_level),
        maximum_level_(maximum_level) {}

  static absl::StatusOr<int> ComputeLogMaximalLevel(
      absl::Span<const PrimeModulus<ModularInt>* const> moduli) {
    // Let's compute the maximal level: after multiplication of two numbers
    // the coefficients of a lazy polynomial will be
    // < modulus^2 < 2^(2 * log_modulus). LazyPolynomial can only be added
    // to vectors with coefficients < 2^(2 * log_modulus), resulting to a vector
    // with coefficients 2^(2 * log_modulus + 1), and so on. The maximum level,
    // at which we need to reduce the coefficients of the lazy polynomial, is
    // 2^(log_modulus+maximum_level) = 2^(bitsize(BigInt) - 1).
    int log_maximum_level = 8 * sizeof(Integer) - 2;
    for (int i = 0; i < moduli.size(); ++i) {
      auto mod_params_qi = moduli[i]->ModParams();
      int log_maximum_level_qi =
          mod_params_qi->bitsize_bigint - 1 - 2 * mod_params_qi->log_modulus;
      if (log_maximum_level_qi < 1) {
        return absl::InvalidArgumentError(absl::StrCat(
            "The RNS moduli do not allow for lazy polynomials, as the "
            "logarithm of the maximum level, ",
            log_maximum_level_qi, " is stricly smaller than 1."));
      }
      if (log_maximum_level_qi < log_maximum_level) {
        log_maximum_level = log_maximum_level_qi;
      }
    }
    return log_maximum_level;
  }

  // Reduces the coefficients using the provided modulus parameters.
  void Refresh(absl::Span<const PrimeModulus<ModularInt>* const> moduli) {
    for (int i = 0; i < coeff_vectors_.size(); ++i) {
      auto mod_params_qi = moduli[i]->ModParams();
      for (auto& coeff : coeff_vectors_[i]) {
        coeff = mod_params_qi->BarrettReduceBigInt(coeff);
      }
    }
    current_level_ = static_cast<BigInt>(1);
  }

  absl::Status CheckFusedMulAddInPlaceParameters(
      const RnsPolynomial<ModularInt>& a, const RnsPolynomial<ModularInt>& b,
      absl::Span<const PrimeModulus<ModularInt>* const> moduli);

  // Coefficients of the polynomial modulo prime moduli.
  // Each vector corresponds to a prime modulus in moduli_ in the same order.
  std::vector<hwy::AlignedVector<BigInt>> coeff_vectors_;

  BigInt current_level_;
  BigInt maximum_level_;
};

}  // namespace rlwe

#endif  // RLWE_RNS_LAZY_RNS_POLYNOMIAL_H_
