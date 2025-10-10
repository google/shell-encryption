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

#include "shell_encryption/multi_party/public_key_share.h"

#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "absl/numeric/int128.h"
#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "absl/types/span.h"
#include "shell_encryption/integral_types.h"
#include "shell_encryption/montgomery.h"
#include "shell_encryption/multi_party/public_parameter.h"
#include "shell_encryption/multi_party/secret_key_share.h"
#include "shell_encryption/multi_party/polynomial_utilities.h"
#include "shell_encryption/prng/prng.h"
#include "shell_encryption/prng/single_thread_chacha_prng.h"
#include "shell_encryption/prng/single_thread_hkdf_prng.h"
#include "shell_encryption/rns/error_distribution.h"
#include "shell_encryption/rns/rns_context.h"
#include "shell_encryption/rns/rns_modulus.h"
#include "shell_encryption/rns/rns_polynomial.h"
#include "shell_encryption/status_macros.h"

namespace rlwe {
namespace multi_party {

namespace {



}  // namespace

template <typename ModularInt>
absl::StatusOr<PublicKeyShare<ModularInt>> PublicKeyShare<ModularInt>::Create(
    const SecretKeyShare<ModularInt>* secret_key_share,
    const PublicParameter<ModularInt>* public_parameter, PrngType prng_type) {
  if (secret_key_share == nullptr) {
    return absl::InvalidArgumentError("`secret_key_share` must not be null.");
  }
  if (public_parameter == nullptr) {
    return absl::InvalidArgumentError("`public_parameter` must not be null.");
  }

  // Create a PRNG to sample the error term.
  std::unique_ptr<SecurePrng> prng_error;
  std::string prng_error_seed;
  if (prng_type == PRNG_TYPE_HKDF) {
    RLWE_ASSIGN_OR_RETURN(prng_error_seed,
                          SingleThreadHkdfPrng::GenerateSeed());
    RLWE_ASSIGN_OR_RETURN(prng_error,
                          SingleThreadHkdfPrng::Create(prng_error_seed));
  } else if (prng_type == PRNG_TYPE_CHACHA) {
    RLWE_ASSIGN_OR_RETURN(prng_error_seed,
                          SingleThreadChaChaPrng::GenerateSeed());
    RLWE_ASSIGN_OR_RETURN(prng_error,
                          SingleThreadChaChaPrng::Create(prng_error_seed));
  } else {
    return absl::InvalidArgumentError("`prng_type` not specified correctly.");
  }

  // Initialize return value
  auto moduli = public_parameter->Moduli();
  auto key_b = RnsPolynomial<ModularInt>::CreateEmpty();
  auto status = CreateExplicit(secret_key_share->Key(), public_parameter,
                               prng_error.get(), &key_b, /*key_error=*/nullptr,
                               /*wrap_around=*/nullptr);
  if (!status.ok()) {
    return status;
  }

  std::vector<const PrimeModulus<ModularInt>*> moduli_vector;
  moduli_vector.insert(moduli_vector.begin(), moduli.begin(), moduli.end());

  return PublicKeyShare<ModularInt>(std::move(key_b), std::move(moduli_vector));
}

template <typename ModularInt>
absl::Status PublicKeyShare<ModularInt>::CreateExplicit(
    const RnsPolynomial<ModularInt>& secret_key_share,
    const PublicParameter<ModularInt>* public_parameter, SecurePrng* prng,
    RnsPolynomial<ModularInt>* key_b, RnsPolynomial<ModularInt>* key_error,
    RnsPolynomial<ModularInt>* wrap_around) {
  if (public_parameter == nullptr) {
    return absl::InvalidArgumentError("`public_parameter` must not be null.");
  }
  if (prng == nullptr) {
    return absl::InvalidArgumentError("`prng` must not be null.");
  }
  if (key_b == nullptr) {
    return absl::InvalidArgumentError("`key_b` must not be null.");
  }

  auto moduli = public_parameter->Moduli();
  int log_n = public_parameter->LogN();

  // Sample the error polynomial e and optionally save it.
  RLWE_ASSIGN_OR_RETURN(
      RnsPolynomial<ModularInt> b,
      SampleError<ModularInt>(log_n, public_parameter->ErrorVariance(), moduli,
                              prng));
  if (key_error != nullptr) {
    *key_error = b;
  }

  // u = -a where a is the public random polynomial.
  RLWE_ASSIGN_OR_RETURN(RnsPolynomial<ModularInt> u,
                        public_parameter->PublicKeyComponentA().Negate(moduli));
  // us = -a * s.
  RLWE_ASSIGN_OR_RETURN(RnsPolynomial<ModularInt> us,
                        u.Mul(secret_key_share, moduli));

  // b = -a * s + e.
  RLWE_RETURN_IF_ERROR(b.AddInPlace(us, moduli));
  *key_b = b;

  // Optionally save wrap_around such that
  // b = -a * s + e + wrap_around * (X^N + 1) over Z_Q[X].
  if (wrap_around != nullptr) {
    RLWE_ASSIGN_OR_RETURN(*wrap_around,
      rlwe_internal::QuotientOf(u, secret_key_share, us, moduli));
  }

  return absl::OkStatus();
}

template <typename ModularInt>
absl::StatusOr<PublicKeyShare<ModularInt>>
PublicKeyShare<ModularInt>::Deserialize(
    const SerializedPublicKeyShare& serialized,
    absl::Span<const PrimeModulus<ModularInt>* const> moduli) {
  RLWE_ASSIGN_OR_RETURN(
      RnsPolynomial<ModularInt> key_b,
      RnsPolynomial<ModularInt>::Deserialize(serialized.key_b(), moduli));
  std::vector<const PrimeModulus<ModularInt>*> moduli_vector;
  moduli_vector.insert(moduli_vector.begin(), moduli.begin(), moduli.end());
  return PublicKeyShare<ModularInt>(std::move(key_b), std::move(moduli_vector));
}

template class PublicKeyShare<MontgomeryInt<Uint32>>;
template class PublicKeyShare<MontgomeryInt<Uint64>>;
template class PublicKeyShare<MontgomeryInt<absl::uint128>>;
#ifdef ABSL_HAVE_INTRINSIC_INT128
template class PublicKeyShare<MontgomeryInt<unsigned __int128>>;
#endif

}  // namespace multi_party
}  // namespace rlwe
