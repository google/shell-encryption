// Copyright 2023 Google LLC
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
#include "shell_encryption/integral_types.h"
#include "shell_encryption/montgomery.h"
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
    const DiscreteGaussianSampler<Integer>* dg_sampler,
    SecurePrng* prng) const {
  if (dg_sampler == nullptr) {
    return absl::InvalidArgumentError("`dg_sampler` must not be null.");
  }
  if (prng == nullptr) {
    return absl::InvalidArgumentError("`prng` must not be null.");
  }

  // Compute partial decryption d = s * ciphertext_component_a + e_flood.
  RnsPolynomial<ModularInt> d = ciphertext_component_a;
  if (!d.IsNttForm()) {
    RLWE_RETURN_IF_ERROR(d.ConvertToNttForm(moduli_));
  }
  RLWE_RETURN_IF_ERROR(d.MulInPlace(key_, moduli_));

  RLWE_ASSIGN_OR_RETURN(RnsPolynomial<ModularInt> e_flood,
                        SampleDiscreteGaussian<ModularInt>(
                            LogN(), s_flood, moduli_, dg_sampler, prng));
  if (!e_flood.IsNttForm()) {
    RLWE_RETURN_IF_ERROR(e_flood.ConvertToNttForm(moduli_));
  }
  RLWE_RETURN_IF_ERROR(d.AddInPlace(e_flood, moduli_));
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
