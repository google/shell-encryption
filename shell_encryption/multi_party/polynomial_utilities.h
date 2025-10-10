#ifndef RLWE_MULTI_PARTY_POLYNOMIAL_UTILITIES_H_
#define RLWE_MULTI_PARTY_POLYNOMIAL_UTILITIES_H_

#include <memory>
#include <vector>

#include "absl/status/statusor.h"
#include "absl/types/span.h"
#include "shell_encryption/rns/rns_context.h"
#include "shell_encryption/rns/rns_modulus.h"
#include "shell_encryption/rns/rns_polynomial.h"
#include "shell_encryption/status_macros.h"

namespace rlwe {
namespace multi_party {
namespace rlwe_internal {

// Returns the RNS constants for the ambient ring Z_Q[X]/(X^{2N} + 1), where
// `moduli` contains RNS moduli parameters of the ring Z_Q[X]/(X^N + 1).
template <typename ModularInt>
absl::StatusOr<std::unique_ptr<const RnsContext<ModularInt>>>
CreateAmbientRnsContext(
    int log_n, absl::Span<const PrimeModulus<ModularInt>* const> moduli) {
  std::vector<typename ModularInt::Int> qs;
  for (auto modulus : moduli) {
    qs.push_back(modulus->Modulus());
  }
  // We don't really need a plaintext modulus for the ambient ring; use 2 as
  // it is compatible with any RNS moduli.
  RLWE_ASSIGN_OR_RETURN(
      RnsContext<ModularInt> rns_context_ambient,
      RnsContext<ModularInt>::Create(log_n + 1, qs, /*ps=*/{},
                                     /*plaintext_modulus=*/2));
  return std::make_unique<const RnsContext<ModularInt>>(
      std::move(rns_context_ambient));
}

// Given polynomial a(X) \in Z_Q[X]/(X^N + 1), returns the polynomial a'(X) \in
// Z_Q[X]/(X^{2N} + 1) whose coefficients of the lower order degrees are those
// of `a`, i.e. a'[i] = a[i] for all 0 <= i < N.
template <typename ModularInt>
absl::StatusOr<RnsPolynomial<ModularInt>> LiftRnsPolynomial(
    RnsPolynomial<ModularInt> a,
    absl::Span<const PrimeModulus<ModularInt>* const> moduli) {
  if (a.IsNttForm()) {
    RLWE_RETURN_IF_ERROR(a.ConvertToCoeffForm(moduli));
  }
  int m = 1 << (a.LogN() + 1);
  std::vector<std::vector<ModularInt>> a_ambient_coeffs = a.Coeffs();
  for (int i = 0; i < a.NumModuli(); ++i) {
    auto mod_params_qi = moduli[i]->ModParams();
    auto zero_mod_qi = ModularInt::ImportZero(mod_params_qi);
    a_ambient_coeffs[i].resize(m, zero_mod_qi);
  }
  return RnsPolynomial<ModularInt>::Create(std::move(a_ambient_coeffs),
                                           /*is_ntt=*/false);
}

// Given polynomials a(X), b(X), c(X) \in Z_Q[X]/(X^N + 1) such that c = a * b,
// returns the quotient polynomial v(X) \in Z_Q[X] such that
// c(X) = a(X) * b(X) + v(X) * (X^N + 1), where `a`, `b`, and `c` are
// considered as polynomials in Z_Q[X].
// Note that the RNS modulus Q of a, b, c must be a product of NTT-friendly
// primes wrt the higher order cyclotomic X^{2N} + 1.
template <typename ModularInt>
absl::StatusOr<RnsPolynomial<ModularInt>> QuotientOf(
    const RnsPolynomial<ModularInt>& a, const RnsPolynomial<ModularInt>& b,
    const RnsPolynomial<ModularInt>& c,
    absl::Span<const PrimeModulus<ModularInt>* const> moduli) {
  // We compute v(X) in the ambient ring R' = Z_Q[X]/(X^{2N} + 1), which is
  // large enough such that, when both a and b are lifted to R', a(X) * b(X)
  // has degree < 2N and so a(X) * b(X) does not wrap around modulo X^{2N} + 1.
  int log_n = a.LogN();
  RLWE_ASSIGN_OR_RETURN(auto ambient_rns_context,
                        CreateAmbientRnsContext(log_n, moduli));
  auto ambient_moduli = ambient_rns_context->MainPrimeModuli();
  RLWE_ASSIGN_OR_RETURN(auto amb_ab, LiftRnsPolynomial(a, moduli));
  RLWE_RETURN_IF_ERROR(amb_ab.ConvertToNttForm(ambient_moduli));
  RLWE_ASSIGN_OR_RETURN(auto amb_b, LiftRnsPolynomial(b, moduli));
  RLWE_RETURN_IF_ERROR(amb_b.ConvertToNttForm(ambient_moduli));
  RLWE_ASSIGN_OR_RETURN(auto amb_c, LiftRnsPolynomial(c, moduli));
  RLWE_RETURN_IF_ERROR(amb_c.ConvertToNttForm(ambient_moduli));

  // a(X) * b(X) \in R'
  RLWE_RETURN_IF_ERROR(amb_ab.MulInPlace(amb_b, ambient_moduli));
  // c(X) - a(X) * b(X) \in R'
  RLWE_ASSIGN_OR_RETURN(auto amb_v, amb_c.Sub(amb_ab, ambient_moduli));

  // (X^N + 1)^-1 mod (Q, X^2N + 1) = (Q-1)/2 * X^N + (Q+1)/2
  int n = 1 << log_n;
  int num_moduli = moduli.size();
  std::vector<std::vector<ModularInt>> f_inv_coeffs(num_moduli);
  for (int i = 0; i < num_moduli; ++i) {
    auto mod_params_qi = moduli[i]->ModParams();
    auto qi = moduli[i]->Modulus();
    f_inv_coeffs[i].resize(2 * n, ModularInt::ImportZero(mod_params_qi));
    RLWE_ASSIGN_OR_RETURN(f_inv_coeffs[i][0],
                          ModularInt::ImportInt((qi + 1) / 2, mod_params_qi));
    RLWE_ASSIGN_OR_RETURN(f_inv_coeffs[i][n],
                          ModularInt::ImportInt((qi - 1) / 2, mod_params_qi));
  }
  RLWE_ASSIGN_OR_RETURN(
      auto f_inv, RnsPolynomial<ModularInt>::Create(std::move(f_inv_coeffs),
                                                    /*is_ntt=*/false));
  RLWE_RETURN_IF_ERROR(f_inv.ConvertToNttForm(ambient_moduli));

  // v = (c - a * b) * (X^N + 1)^-1 mod (Q, X^2N + 1)
  RLWE_RETURN_IF_ERROR(amb_v.MulInPlace(f_inv, ambient_moduli));
  RLWE_RETURN_IF_ERROR(amb_v.ConvertToCoeffForm(ambient_moduli));

  // The high order coefficients of v(X) should all be 0.
  std::vector<std::vector<ModularInt>> v_coeffs = amb_v.Coeffs();
  for (int i = 0; i < num_moduli; ++i) {
    v_coeffs[i].resize(n, ModularInt::ImportZero(moduli[i]->ModParams()));
  }
  return RnsPolynomial<ModularInt>::Create(std::move(v_coeffs),
                                           /*is_ntt=*/false);
}
}  // namespace rlwe_internal

}  // namespace multi_party
}  // namespace rlwe
#endif  // RLWE_MULTI_PARTY_POLYNOMIAL_UTILITIES_H_
