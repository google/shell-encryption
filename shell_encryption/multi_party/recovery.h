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

#ifndef RLWE_MULTI_PARTY_RECOVERY_H_
#define RLWE_MULTI_PARTY_RECOVERY_H_

#include <vector>

#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "absl/types/span.h"
#include "shell_encryption/multi_party/public_parameter.h"
#include "shell_encryption/rns/coefficient_encoder.h"
#include "shell_encryption/rns/rns_polynomial.h"
#include "shell_encryption/status_macros.h"

namespace rlwe {
namespace multi_party {

template <typename ModularInt,
          typename Encoder = CoefficientEncoder<ModularInt>>
absl::StatusOr<std::vector<typename ModularInt::Int>> RecoverMessages(
    absl::Span<const RnsPolynomial<ModularInt>> partial_decryptions,
    const RnsPolynomial<ModularInt>& ciphertext_component_b,
    const PublicParameter<ModularInt>& public_parameter,
    const Encoder* encoder) {
  if (partial_decryptions.empty()) {
    return absl::InvalidArgumentError(
        "`partial_decryptions` must not be empty.");
  }

  auto moduli = public_parameter.Moduli();
  RLWE_ASSIGN_OR_RETURN(RnsPolynomial<ModularInt> noisy_plaintext,
                        RnsPolynomial<ModularInt>::CreateZero(
                            public_parameter.LogN(), moduli, /*is_ntt=*/true));

  for (auto const& partial_decryption : partial_decryptions) {
    RLWE_RETURN_IF_ERROR(
        noisy_plaintext.AddInPlace(partial_decryption, moduli));
  }

  return RecoverMessagesFromSum(noisy_plaintext, ciphertext_component_b,
                                public_parameter, encoder);
}

template <typename ModularInt,
          typename Encoder = CoefficientEncoder<ModularInt>>
absl::StatusOr<std::vector<typename ModularInt::Int>> RecoverMessagesFromSum(
    const RnsPolynomial<ModularInt>& sum_partial_decryptions,
    const RnsPolynomial<ModularInt>& ciphertext_component_b,
    const PublicParameter<ModularInt>& public_parameter,
    const Encoder* encoder) {
  if (encoder == nullptr) {
    return absl::InvalidArgumentError("`encoder` must not be null.");
  }

  auto moduli = public_parameter.Moduli();

  RLWE_ASSIGN_OR_RETURN(
      RnsPolynomial<ModularInt> noisy_plaintext,
      sum_partial_decryptions.Add(ciphertext_component_b, moduli));

  // Decode the noisy plaintext polynomial.
  return encoder->DecodeBfv(std::move(noisy_plaintext), moduli);
}

}  // namespace multi_party
}  // namespace rlwe

#endif  // RLWE_MULTI_PARTY_RECOVERY_H_
