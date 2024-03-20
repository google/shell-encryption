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

#include "shell_encryption/rns/rns_public_key.h"

#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "absl/numeric/int128.h"
#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "shell_encryption/integral_types.h"
#include "shell_encryption/montgomery.h"
#include "shell_encryption/prng/prng.h"
#include "shell_encryption/prng/single_thread_chacha_prng.h"
#include "shell_encryption/prng/single_thread_hkdf_prng.h"
#include "shell_encryption/rns/error_distribution.h"
#include "shell_encryption/rns/rns_modulus.h"
#include "shell_encryption/rns/rns_polynomial.h"
#include "shell_encryption/rns/rns_secret_key.h"
#include "shell_encryption/status_macros.h"

namespace rlwe {

template <typename ModularInt>
absl::StatusOr<RnsRlwePublicKey<ModularInt>>
RnsRlwePublicKey<ModularInt>::Create(
    const RnsRlweSecretKey<ModularInt>& secret_key, int variance,
    PrngType prng_type, Integer error_scalar) {
  if (variance <= 0) {
    return absl::InvalidArgumentError("`variance` must be positive.");
  }

  // Create two PRNGs: `prng_pad` is used to sample the random polynomial "a"
  // and its seed is stored such that it can be serialized. The other PRNG
  // `prng_error` is used to sample the error term and it is not stored.
  std::unique_ptr<SecurePrng> prng_pad, prng_error;
  std::string prng_pad_seed, prng_error_seed;
  if (prng_type == PRNG_TYPE_HKDF) {
    RLWE_ASSIGN_OR_RETURN(prng_pad_seed, SingleThreadHkdfPrng::GenerateSeed());
    RLWE_ASSIGN_OR_RETURN(prng_pad,
                          SingleThreadHkdfPrng::Create(prng_pad_seed));
    RLWE_ASSIGN_OR_RETURN(prng_error_seed,
                          SingleThreadHkdfPrng::GenerateSeed());
    RLWE_ASSIGN_OR_RETURN(prng_error,
                          SingleThreadHkdfPrng::Create(prng_error_seed));
  } else if (prng_type == PRNG_TYPE_CHACHA) {
    RLWE_ASSIGN_OR_RETURN(prng_pad_seed,
                          SingleThreadChaChaPrng::GenerateSeed());
    RLWE_ASSIGN_OR_RETURN(prng_pad,
                          SingleThreadChaChaPrng::Create(prng_pad_seed));
    RLWE_ASSIGN_OR_RETURN(prng_error_seed,
                          SingleThreadChaChaPrng::GenerateSeed());
    RLWE_ASSIGN_OR_RETURN(prng_error,
                          SingleThreadChaChaPrng::Create(prng_error_seed));
  } else {
    return absl::InvalidArgumentError("PrngType not specified correctly.");
  }

  // b = u * s + t * e for a uniformly random u.
  auto moduli = secret_key.Moduli();
  int log_n = secret_key.LogN();
  RLWE_ASSIGN_OR_RETURN(auto u, RnsPolynomial<ModularInt>::SampleUniform(
                                    log_n, prng_pad.get(), moduli));
  RLWE_ASSIGN_OR_RETURN(
      RnsPolynomial<ModularInt> b,
      SampleError<ModularInt>(log_n, variance, moduli, prng_error.get()));
  RLWE_RETURN_IF_ERROR(b.MulInPlace(error_scalar, moduli));
  RLWE_RETURN_IF_ERROR(b.FusedMulAddInPlace(u, secret_key.Key(), moduli));

  // a = -u
  RLWE_RETURN_IF_ERROR(u.NegateInPlace(moduli));

  // Store the RNS moduli.
  std::vector<const PrimeModulus<ModularInt>*> moduli_vector;
  moduli_vector.insert(moduli_vector.begin(), moduli.begin(), moduli.end());

  return RnsRlwePublicKey<ModularInt>(std::move(u), std::move(b),
                                      std::move(moduli_vector), variance,
                                      prng_pad_seed, prng_type);
}

template class RnsRlwePublicKey<MontgomeryInt<Uint16>>;
template class RnsRlwePublicKey<MontgomeryInt<Uint32>>;
template class RnsRlwePublicKey<MontgomeryInt<Uint64>>;
template class RnsRlwePublicKey<MontgomeryInt<absl::uint128>>;
#ifdef ABSL_HAVE_INTRINSIC_INT128
template class RnsRlwePublicKey<MontgomeryInt<unsigned __int128>>;
#endif

}  // namespace rlwe
