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

#include "shell_encryption/rns/rns_secret_key.h"

#include <utility>
#include <vector>

#include "absl/numeric/int128.h"
#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "shell_encryption/integral_types.h"
#include "shell_encryption/montgomery.h"
#include "shell_encryption/prng/prng.h"
#include "shell_encryption/rns/error_distribution.h"
#include "shell_encryption/rns/rns_ciphertext.h"
#include "shell_encryption/rns/rns_modulus.h"
#include "shell_encryption/rns/rns_polynomial.h"
#include "shell_encryption/status_macros.h"

namespace rlwe {

template <typename ModularInt>
absl::StatusOr<RnsRlweSecretKey<ModularInt>>
RnsRlweSecretKey<ModularInt>::Sample(
    int log_n, int variance,
    std::vector<const PrimeModulus<ModularInt>*> moduli, SecurePrng* prng) {
  if (log_n <= 0) {
    return absl::InvalidArgumentError("`log_n` must be positive.");
  }
  if (variance <= 0) {
    return absl::InvalidArgumentError("`variance` must be positive.");
  }
  if (moduli.empty()) {
    return absl::InvalidArgumentError(
        "`moduli` must contain at least one element.");
  }
  if (prng == nullptr) {
    return absl::InvalidArgumentError("`prng` must not be null.");
  }

  // Sample the secret s from the error distribution.
  RLWE_ASSIGN_OR_RETURN(RnsPolynomial<ModularInt> s,
                        SampleError<ModularInt>(log_n, variance, moduli, prng));
  return RnsRlweSecretKey(std::move(s), std::move(moduli), variance);
}

template <typename ModularInt>
absl::StatusOr<RnsPolynomial<ModularInt>>
RnsRlweSecretKey<ModularInt>::RawDecrypt(
    const RnsRlweCiphertext<ModularInt>& ciphertext) const {
  RnsPolynomial<ModularInt> s_power = key_;
  RLWE_ASSIGN_OR_RETURN(RnsPolynomial<ModularInt> output,
                        RnsPolynomial<ModularInt>::CreateZero(LogN(), moduli_));
  int ciphertext_len = ciphertext.Len();
  for (int i = 0; i < ciphertext_len; ++i) {
    // Get the i-th component
    RLWE_ASSIGN_OR_RETURN(RnsPolynomial<ModularInt> ci,
                          ciphertext.Component(i));

    // Compute the next power of the secret polynomial s.
    if (i > 1) {
      RLWE_RETURN_IF_ERROR(s_power.MulInPlace(key_, moduli_));
    }
    // Compute c[i] * s^i.
    if (i > 0) {
      RLWE_RETURN_IF_ERROR(ci.MulInPlace(s_power, moduli_));
    }
    // Add c[i] * s^i to the result.
    RLWE_RETURN_IF_ERROR(output.AddInPlace(ci, moduli_));
  }
  return output;
}

template <typename ModularInt>
absl::Status RnsRlweSecretKey<ModularInt>::ModReduce() {
  if (moduli_.size() <= 1) {
    return absl::FailedPreconditionError(
        "Cannot perform ModReduce with insufficient number of prime moduli.");
  }
  RLWE_ASSIGN_OR_RETURN(auto k_ql, key_.DetachLastCoeffVector());
  moduli_.pop_back();
  return absl::OkStatus();
}

template class RnsRlweSecretKey<MontgomeryInt<Uint16>>;
template class RnsRlweSecretKey<MontgomeryInt<Uint32>>;
template class RnsRlweSecretKey<MontgomeryInt<Uint64>>;
template class RnsRlweSecretKey<MontgomeryInt<absl::uint128>>;
#ifdef ABSL_HAVE_INTRINSIC_INT128
template class RnsRlweSecretKey<MontgomeryInt<unsigned __int128>>;
#endif

}  // namespace rlwe
