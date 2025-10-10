// Copyright 2025 Google LLC
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include "shell_encryption/multi_party/secret_key_share.h"

#include <utility>
#include <vector>

#include "absl/numeric/int128.h"
#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "absl/types/span.h"
#include "shell_encryption/integral_types.h"
#include "shell_encryption/montgomery.h"
#include "shell_encryption/multi_party/polynomial_utilities.h"
#include "shell_encryption/prng/prng.h"
#include "shell_encryption/rns/error_distribution.h"
#include "shell_encryption/rns/rns_context.h"
#include "shell_encryption/rns/rns_modulus.h"
#include "shell_encryption/rns/rns_polynomial.h"
#include "shell_encryption/sampler/discrete_gaussian.h"
#include "shell_encryption/status_macros.h"

namespace rlwe {
namespace multi_party {

template <typename ModularInt>
absl::StatusOr<SecretKeyShare<ModularInt>> SecretKeyShare<ModularInt>::Sample(
    const RnsContext<ModularInt>* rns_context, SecurePrng* prng) {
  if (rns_context == nullptr) {
    return absl::InvalidArgumentError("`rns_context` must not be null.");
  }
  if (prng == nullptr) {
    return absl::InvalidArgumentError("`prng` must not be null.");
  }

  // Sample the secret s with uniform ternary coefficients.
  int log_n = rns_context->LogN();
  std::vector<const PrimeModulus<ModularInt>*> moduli =
      rns_context->MainPrimeModuli();
  RLWE_ASSIGN_OR_RETURN(RnsPolynomial<ModularInt> s,
                        SampleUniformTernary<ModularInt>(log_n, moduli, prng));
  return SecretKeyShare(std::move(s), std::move(moduli));
}

template <typename ModularInt>
absl::StatusOr<RnsPolynomial<ModularInt>>
SecretKeyShare<ModularInt>::PartialDecrypt(
    const RnsPolynomial<ModularInt>& ciphertext_component_a, double s_flood,
    const DiscreteGaussianSampler<Integer>* dg_sampler, SecurePrng* prng,
    RnsPolynomial<ModularInt>* error_flood,
    RnsPolynomial<ModularInt>* wrap_around) const {
  if (dg_sampler == nullptr) {
    return absl::InvalidArgumentError("`dg_sampler` must not be null.");
  }
  if (prng == nullptr) {
    return absl::InvalidArgumentError("`prng` must not be null.");
  }

  // Sample e_flood.
  RLWE_ASSIGN_OR_RETURN(RnsPolynomial<ModularInt> d,
                        SampleDiscreteGaussian<ModularInt>(
                            LogN(), s_flood, moduli_, dg_sampler, prng));
  if (error_flood != nullptr) {
    *error_flood = d;
  }

  if (!d.IsNttForm()) {
    RLWE_RETURN_IF_ERROR(d.ConvertToNttForm(moduli_));
  }

  RnsPolynomial<ModularInt> c1 = ciphertext_component_a;
  if (!c1.IsNttForm()) {
    RLWE_RETURN_IF_ERROR(c1.ConvertToNttForm(moduli_));
  }

  // c1s = c1 * s.
  RLWE_ASSIGN_OR_RETURN(RnsPolynomial<ModularInt> c1s, c1.Mul(key_, moduli_));

  // d = c1 * s + e_flood.
  RLWE_RETURN_IF_ERROR(d.AddInPlace(c1s, moduli_));

  // Optionally save wraparound such that
  // d = c1 * s + e_flood + wraparound * (X^N + 1) over Z_Q[X].
  if (wrap_around != nullptr) {
    RLWE_ASSIGN_OR_RETURN(
        *wrap_around,
        rlwe_internal::QuotientOf(c1, key_, c1s, absl::MakeSpan(moduli_)));
  }

  return d;
}

template class SecretKeyShare<MontgomeryInt<Uint32>>;
template class SecretKeyShare<MontgomeryInt<Uint64>>;
template class SecretKeyShare<MontgomeryInt<absl::uint128>>;
#ifdef ABSL_HAVE_INTRINSIC_INT128
template class SecretKeyShare<MontgomeryInt<unsigned __int128>>;
#endif

}  // namespace multi_party
}  // namespace rlwe
