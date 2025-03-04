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

#include "shell_encryption/multi_party/public_key.h"

#include <utility>

#include "absl/numeric/int128.h"
#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "absl/types/span.h"
#include "shell_encryption/integral_types.h"
#include "shell_encryption/montgomery.h"
#include "shell_encryption/multi_party/public_key_share.h"
#include "shell_encryption/multi_party/public_parameter.h"
#include "shell_encryption/rns/rns_polynomial.h"
#include "shell_encryption/status_macros.h"

namespace rlwe {
namespace multi_party {

template <typename ModularInt>
absl::StatusOr<PublicKey<ModularInt>> PublicKey<ModularInt>::Create(
    const PublicParameter<ModularInt>* public_parameter,
    absl::Span<const PublicKeyShare<ModularInt>> public_key_shares) {
  if (public_parameter == nullptr) {
    return absl::InvalidArgumentError("`public_parameter` must not be null.");
  }
  if (public_key_shares.empty()) {
    return absl::InvalidArgumentError("`public_key_shares` must not be empty.");
  }

  RLWE_ASSIGN_OR_RETURN(RnsPolynomial<ModularInt> key_b,
                        RnsPolynomial<ModularInt>::CreateZero(
                            public_parameter->LogN(),
                            public_parameter->Moduli(), /*is_ntt=*/true));
  for (auto const& public_key_share : public_key_shares) {
    RLWE_RETURN_IF_ERROR(key_b.AddInPlace(public_key_share.ComponentB(),
                                          public_parameter->Moduli()));
  }
  return PublicKey<ModularInt>(public_parameter, std::move(key_b));
}

template <typename ModularInt>
absl::StatusOr<PublicKey<ModularInt>> PublicKey<ModularInt>::Deserialize(
    const SerializedPublicKey& serialized,
    const PublicParameter<ModularInt>* public_parameter) {
  if (public_parameter == nullptr) {
    return absl::InvalidArgumentError("`public_parameter` must not be null.");
  }

  RLWE_ASSIGN_OR_RETURN(RnsPolynomial<ModularInt> key_b,
                        RnsPolynomial<ModularInt>::Deserialize(
                            serialized.key_b(), public_parameter->Moduli()));
  return PublicKey<ModularInt>(public_parameter, std::move(key_b));
}

template class PublicKey<MontgomeryInt<Uint32>>;
template class PublicKey<MontgomeryInt<Uint64>>;
template class PublicKey<MontgomeryInt<absl::uint128>>;
#ifdef ABSL_HAVE_INTRINSIC_INT128
template class PublicKey<MontgomeryInt<unsigned __int128>>;
#endif

}  // namespace multi_party
}  // namespace rlwe
