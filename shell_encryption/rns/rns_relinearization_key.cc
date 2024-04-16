// Copyright 2024 Google LLC
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

#include "shell_encryption/rns/rns_relinearization_key.h"

#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "absl/log/check.h"
#include "absl/numeric/int128.h"
#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "absl/strings/str_cat.h"
#include "absl/strings/string_view.h"
#include "absl/types/span.h"
#include "shell_encryption/integral_types.h"
#include "shell_encryption/montgomery.h"
#include "shell_encryption/prng/prng.h"
#include "shell_encryption/prng/single_thread_chacha_prng.h"
#include "shell_encryption/prng/single_thread_hkdf_prng.h"
#include "shell_encryption/rns/error_distribution.h"
#include "shell_encryption/rns/rns_ciphertext.h"
#include "shell_encryption/rns/rns_gadget.h"
#include "shell_encryption/rns/rns_modulus.h"
#include "shell_encryption/rns/rns_polynomial.h"
#include "shell_encryption/rns/rns_secret_key.h"
#include "shell_encryption/status_macros.h"

namespace rlwe {

template <typename ModularInt>
absl::StatusOr<RnsRelinKey<ModularInt>> RnsRelinKey<ModularInt>::Create(
    const RnsRlweSecretKey<ModularInt>& secret_key, int degree, int variance,
    const RnsGadget<ModularInt>* gadget, PrngType prng_type,
    Integer error_scalar) {
  if (degree <= 1) {
    return absl::InvalidArgumentError("`degree` must be at least 2.");
  }
  if (gadget == nullptr) {
    return absl::InvalidArgumentError("`gadget` must not be null.");
  }

  // Sample a PRNG seed for sampling the random polynomials in `key_as`.
  std::string prng_pad_seed;
  if (prng_type == PRNG_TYPE_HKDF) {
    RLWE_ASSIGN_OR_RETURN(prng_pad_seed, SingleThreadHkdfPrng::GenerateSeed());
  } else if (prng_type == PRNG_TYPE_CHACHA) {
    RLWE_ASSIGN_OR_RETURN(prng_pad_seed,
                          SingleThreadChaChaPrng::GenerateSeed());
  } else {
    return absl::InvalidArgumentError("PrngType not specified correctly.");
  }

  // The relinearization key is a k*(l-1) x 2 matrix [[b2, a2], ..., [bl, al]],
  // where each block [bi, ai] consists of k RLWE encryptions of the target
  // secret g * s^i under the canonical secret key s(X). In the following, we
  // first sample the random polynomials ai's, and then generate bi's.
  RLWE_ASSIGN_OR_RETURN(
      std::vector<RnsPolynomial<ModularInt>> key_as,
      SampleRandomPad(gadget->Dimension(), degree, secret_key.LogN(),
                      secret_key.Moduli(), prng_pad_seed, prng_type));

  return CreateWithRandomPad(std::move(key_as), secret_key, degree, variance,
                             gadget, prng_pad_seed, prng_type, error_scalar);
}

template <typename ModularInt>
absl::StatusOr<std::vector<RnsPolynomial<ModularInt>>>
RnsRelinKey<ModularInt>::SampleRandomPad(
    int dimension, int degree, int log_n,
    absl::Span<const PrimeModulus<ModularInt>* const> moduli,
    absl::string_view prng_seed, PrngType prng_type) {
  if (dimension <= 0) {
    return absl::InvalidArgumentError("`dimension` must be positive.");
  }
  if (degree <= 1) {
    return absl::InvalidArgumentError("`degree` must be at least 2.");
  }
  if (log_n <= 0) {
    return absl::InvalidArgumentError("`log_n` must be positive.");
  }
  // Create a PRNG for sampling the random polynomials in `key_as`.
  std::unique_ptr<SecurePrng> prng_pad;
  if (prng_type == PRNG_TYPE_HKDF) {
    RLWE_ASSIGN_OR_RETURN(prng_pad, SingleThreadHkdfPrng::Create(prng_seed));
  } else if (prng_type == PRNG_TYPE_CHACHA) {
    RLWE_ASSIGN_OR_RETURN(prng_pad, SingleThreadChaChaPrng::Create(prng_seed));
  } else {
    return absl::InvalidArgumentError("PrngType not specified correctly.");
  }

  std::vector<RnsPolynomial<ModularInt>> key_as;
  key_as.reserve(dimension * (degree - 1));
  for (int i = 0; i < dimension * (degree - 1); ++i) {
    RLWE_ASSIGN_OR_RETURN(auto a, RnsPolynomial<ModularInt>::SampleUniform(
                                      log_n, prng_pad.get(), moduli));
    key_as.push_back(std::move(a));
  }
  return key_as;
}

template <typename ModularInt>
absl::StatusOr<RnsRelinKey<ModularInt>>
RnsRelinKey<ModularInt>::CreateWithRandomPad(
    std::vector<RnsPolynomial<ModularInt>> pads,
    const RnsRlweSecretKey<ModularInt>& secret_key, int degree, int variance,
    const RnsGadget<ModularInt>* gadget, absl::string_view prng_pad_seed,
    PrngType prng_type, Integer error_scalar) {
  if (variance <= 0) {
    return absl::InvalidArgumentError("`variance` must be positive.");
  }

  // Create the PRNGs for sampling the encryption randomness.
  std::unique_ptr<SecurePrng> prng_encryption;
  std::string prng_encryption_seed;
  if (prng_type == PRNG_TYPE_HKDF) {
    RLWE_ASSIGN_OR_RETURN(prng_encryption_seed,
                          SingleThreadHkdfPrng::GenerateSeed());
    RLWE_ASSIGN_OR_RETURN(prng_encryption,
                          SingleThreadHkdfPrng::Create(prng_encryption_seed));
  } else {
    RLWE_ASSIGN_OR_RETURN(prng_encryption_seed,
                          SingleThreadChaChaPrng::GenerateSeed());
    RLWE_ASSIGN_OR_RETURN(prng_encryption,
                          SingleThreadChaChaPrng::Create(prng_encryption_seed));
  }

  const RnsPolynomial<ModularInt>& target_key = secret_key.Key();
  RnsPolynomial<ModularInt> secret = target_key;

  // The relinearization key is a k*(l-1) x 2 matrix [[b2, a2], ..., [bl, al]],
  // where ai consists of k uniformly random polynomials mod `moduli` and bi =
  // -ai * s + t * ei + gadget * s^i, for t = `error_scalar`.
  int log_n = secret_key.LogN();
  int k = gadget->Dimension();
  std::vector<RnsPolynomial<ModularInt>> key_bs;
  key_bs.reserve(k);
  int index = 0;  // for polynomials in `pads`.
  auto moduli = secret_key.Moduli();
  for (int i = 2; i <= degree; ++i) {
    RLWE_RETURN_IF_ERROR(secret.MulInPlace(target_key, moduli));  // s^i
    for (int j = 0; j < k; ++j) {
      // a = -u
      RLWE_ASSIGN_OR_RETURN(RnsPolynomial<ModularInt> u,
                            pads[index++].Negate(moduli));

      RnsPolynomial<ModularInt> z = secret;
      RLWE_RETURN_IF_ERROR(z.MulInPlace(gadget->Component(j), moduli));

      RLWE_ASSIGN_OR_RETURN(
          RnsPolynomial<ModularInt> b,
          SampleError<ModularInt>(log_n, variance, moduli,
                                  prng_encryption.get()));       // b = e
      RLWE_RETURN_IF_ERROR(b.MulInPlace(error_scalar, moduli));  // b = t * e
      RLWE_RETURN_IF_ERROR(b.AddInPlace(z, moduli));  // b = t * e + g[j] * s'
      RLWE_RETURN_IF_ERROR(b.FusedMulAddInPlace(u, target_key, moduli));

      key_bs.push_back(std::move(b));
    }
  }

  // Store the RNS moduli.
  std::vector<const PrimeModulus<ModularInt>*> moduli_vector;
  moduli_vector.insert(moduli_vector.begin(), moduli.begin(), moduli.end());

  return RnsRelinKey(/*key_as=*/std::move(pads), std::move(key_bs), gadget,
                     degree, std::move(moduli_vector), prng_pad_seed,
                     prng_type);
}

template <typename ModularInt>
absl::StatusOr<std::vector<RnsPolynomial<ModularInt>>>
RnsRelinKey<ModularInt>::ApplyToRlweCiphertext(
    const RnsRlweCiphertext<ModularInt>& ciphertext) const {
  if (ciphertext.Degree() > degree_) {
    return absl::InvalidArgumentError(
        absl::StrCat("`ciphertext` degree is larger than degree of this "
                     "relinearization key, ",
                     degree_, "."));
  }
  if (ciphertext.NumModuli() != moduli_.size()) {
    return absl::InvalidArgumentError(
        "`ciphertext` does not have a matching RNS moduli set.");
  }
  if (ciphertext.PowerOfS() != 1) {
    return absl::InvalidArgumentError(
        "Relinearization key can only apply to a ciphertext of power 1.");
  }

  // Apply the relinearization key with blocks [b2, a2], ..., [bl, al] to a
  // degree-l ciphertext (c0, ..., cl) to get a new ciphertext (c0', c1') =
  // (c0, c1) + sum(g^-1(ci) * [bi, ai], i = 2..l).
  int k = gadget_->Dimension();
  int l = ciphertext.Degree();
  RLWE_ASSIGN_OR_RETURN(RnsPolynomial<ModularInt> c0_new,
                        ciphertext.Component(0));
  RLWE_ASSIGN_OR_RETURN(RnsPolynomial<ModularInt> c1_new,
                        ciphertext.Component(1));
  for (int i = 2; i <= l; ++i) {
    RLWE_ASSIGN_OR_RETURN(RnsPolynomial<ModularInt> ci,
                          ciphertext.Component(i));
    if (ci.IsNttForm()) {
      RLWE_RETURN_IF_ERROR(ci.ConvertToCoeffForm(moduli_));
    }
    RLWE_ASSIGN_OR_RETURN(std::vector<RnsPolynomial<ModularInt>> ci_digits,
                          gadget_->Decompose(ci, moduli_));
    for (int j = 0; j < k; ++j) {
      RLWE_RETURN_IF_ERROR(ci_digits[j].ConvertToNttForm(moduli_));
      int index = (i - 2) * l + j;
      RLWE_RETURN_IF_ERROR(
          c0_new.FusedMulAddInPlace(ci_digits[j], key_bs_[index], moduli_));
      RLWE_RETURN_IF_ERROR(
          c1_new.FusedMulAddInPlace(ci_digits[j], key_as_[index], moduli_));
    }
  }

  return std::vector<RnsPolynomial<ModularInt>>{std::move(c0_new),
                                                std::move(c1_new)};
}

template class RnsRelinKey<MontgomeryInt<Uint16>>;
template class RnsRelinKey<MontgomeryInt<Uint32>>;
template class RnsRelinKey<MontgomeryInt<Uint64>>;
template class RnsRelinKey<MontgomeryInt<absl::uint128>>;
#ifdef ABSL_HAVE_INTRINSIC_INT128
template class RnsRelinKey<MontgomeryInt<unsigned __int128>>;
#endif

}  // namespace rlwe
