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

#include "shell_encryption/multi_party/public_parameter.h"

#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "absl/memory/memory.h"
#include "absl/numeric/int128.h"
#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "shell_encryption/integral_types.h"
#include "shell_encryption/montgomery.h"
#include "shell_encryption/multi_party/serialization.pb.h"
#include "shell_encryption/prng/prng.h"
#include "shell_encryption/prng/single_thread_chacha_prng.h"
#include "shell_encryption/prng/single_thread_hkdf_prng.h"
#include "shell_encryption/rns/rns_context.h"
#include "shell_encryption/rns/rns_modulus.h"
#include "shell_encryption/rns/rns_polynomial.h"
#include "shell_encryption/status_macros.h"

namespace rlwe {
namespace multi_party {

template <typename ModularInt>
absl::StatusOr<std::unique_ptr<const PublicParameter<ModularInt>>>
PublicParameter<ModularInt>::CreateFromSeed(
    const RnsContext<ModularInt>* rns_context, int error_variance,
    std::string prng_seed, PrngType prng_type) {
  if (rns_context == nullptr) {
    return absl::InvalidArgumentError("`rns_context` must not be null.");
  }
  if (error_variance <= 0) {
    return absl::InvalidArgumentError("`error_variance` must be positive.");
  }

  std::unique_ptr<SecurePrng> prng;
  if (prng_type == PRNG_TYPE_HKDF) {
    RLWE_ASSIGN_OR_RETURN(prng, SingleThreadHkdfPrng::Create(prng_seed));

  } else if (prng_type == PRNG_TYPE_CHACHA) {
    RLWE_ASSIGN_OR_RETURN(prng, SingleThreadChaChaPrng::Create(prng_seed));
  } else {
    return absl::InvalidArgumentError("`prng_type` not specified correctly.");
  }

  int log_n = rns_context->LogN();
  std::vector<const PrimeModulus<ModularInt>*> moduli =
      rns_context->MainPrimeModuli();
  RLWE_ASSIGN_OR_RETURN(auto key_a, RnsPolynomial<ModularInt>::SampleUniform(
                                        log_n, prng.get(), moduli));
  return absl::WrapUnique(new PublicParameter(prng_seed, prng_type,
                                              error_variance, std::move(key_a),
                                              std::move(moduli)));
}

template <typename ModularInt>
absl::StatusOr<std::unique_ptr<const PublicParameter<ModularInt>>>
PublicParameter<ModularInt>::Create(const RnsContext<ModularInt>* rns_context,
                                    int error_variance, PrngType prng_type) {
  // Generate a PRNG seed for sampling the public random polynomial `a`.
  std::string prng_seed;
  if (prng_type == PRNG_TYPE_HKDF) {
    RLWE_ASSIGN_OR_RETURN(prng_seed, SingleThreadHkdfPrng::GenerateSeed());

  } else if (prng_type == PRNG_TYPE_CHACHA) {
    RLWE_ASSIGN_OR_RETURN(prng_seed, SingleThreadChaChaPrng::GenerateSeed());
  } else {
    return absl::InvalidArgumentError("`prng_type` not specified correctly.");
  }

  return CreateFromSeed(rns_context, error_variance, std::move(prng_seed),
                        prng_type);
}

template <typename ModularInt>
absl::StatusOr<std::unique_ptr<PublicParameter<ModularInt>>>
PublicParameter<ModularInt>::Deserialize(
    const SerializedPublicParameter& serialized,
    const RnsContext<ModularInt>* rns_context) {
  if (rns_context == nullptr) {
    return absl::InvalidArgumentError("`rns_context` must not be null.");
  }
  if (serialized.error_variance() <= 0) {
    return absl::InvalidArgumentError("`error_variance` must be positive.");
  }

  // Generate a PRNG seed for sampling the public random polynomial `a`.
  std::unique_ptr<SecurePrng> prng;
  if (serialized.prng_type() == PRNG_TYPE_HKDF) {
    RLWE_ASSIGN_OR_RETURN(prng,
                          SingleThreadHkdfPrng::Create(serialized.prng_seed()));

  } else if (serialized.prng_type() == PRNG_TYPE_CHACHA) {
    RLWE_ASSIGN_OR_RETURN(
        prng, SingleThreadChaChaPrng::Create(serialized.prng_seed()));
  } else {
    return absl::InvalidArgumentError("Invalid `prng_type` deserialized.");
  }

  int log_n = rns_context->LogN();
  std::vector<const PrimeModulus<ModularInt>*> moduli =
      rns_context->MainPrimeModuli();
  RLWE_ASSIGN_OR_RETURN(auto key_a, RnsPolynomial<ModularInt>::SampleUniform(
                                        log_n, prng.get(), moduli));
  return absl::WrapUnique(new PublicParameter(
      serialized.prng_seed(), serialized.prng_type(),
      serialized.error_variance(), std::move(key_a), std::move(moduli)));
}

template class PublicParameter<MontgomeryInt<Uint32>>;
template class PublicParameter<MontgomeryInt<Uint64>>;
template class PublicParameter<MontgomeryInt<absl::uint128>>;
#ifdef ABSL_HAVE_INTRINSIC_INT128
template class PublicParameter<MontgomeryInt<unsigned __int128>>;
#endif

}  // namespace multi_party
}  // namespace rlwe
