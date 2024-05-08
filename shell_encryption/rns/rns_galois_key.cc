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

#include "shell_encryption/rns/rns_galois_key.h"

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
#include "shell_encryption/rns/lazy_rns_polynomial.h"
#include "shell_encryption/rns/rns_bfv_ciphertext.h"
#include "shell_encryption/rns/rns_bgv_ciphertext.h"
#include "shell_encryption/rns/rns_ciphertext.h"
#include "shell_encryption/rns/rns_gadget.h"
#include "shell_encryption/rns/rns_modulus.h"
#include "shell_encryption/rns/rns_polynomial.h"
#include "shell_encryption/rns/rns_secret_key.h"
#include "shell_encryption/status_macros.h"

namespace rlwe {

template <typename ModularInt>
absl::StatusOr<RnsGaloisKey<ModularInt>> RnsGaloisKey<ModularInt>::Create(
    const RnsRlweSecretKey<ModularInt>& secret_key, int power, int variance,
    const RnsGadget<ModularInt>* gadget, PrngType prng_type,
    Integer error_scalar) {
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

  // The Galois key is a k x 2 matrix [[b1, a1], ..., [bk, ak]], where each pair
  // (bi, ai) is an RLWE encryption of the target secret gi * s(X^power) under
  // the canonical secret key s(X). In the following, we first sample the random
  // polynomials ai's, and then create bi's and hence the Galois key.
  RLWE_ASSIGN_OR_RETURN(
      std::vector<RnsPolynomial<ModularInt>> key_as,
      SampleRandomPad(gadget->Dimension(), secret_key.LogN(),
                      secret_key.Moduli(), prng_pad_seed, prng_type));

  return CreateWithRandomPad(std::move(key_as), secret_key, power, variance,
                             gadget, prng_pad_seed, prng_type, error_scalar);
}

template <typename ModularInt>
absl::StatusOr<std::vector<RnsPolynomial<ModularInt>>>
RnsGaloisKey<ModularInt>::SampleRandomPad(
    int dimension, int log_n,
    absl::Span<const PrimeModulus<ModularInt>* const> moduli,
    absl::string_view prng_seed, PrngType prng_type) {
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
  key_as.reserve(dimension);
  for (int i = 0; i < dimension; ++i) {
    RLWE_ASSIGN_OR_RETURN(auto a, RnsPolynomial<ModularInt>::SampleUniform(
                                      log_n, prng_pad.get(), moduli));
    key_as.push_back(std::move(a));
  }
  return key_as;
}

template <typename ModularInt>
absl::StatusOr<RnsGaloisKey<ModularInt>>
RnsGaloisKey<ModularInt>::CreateWithRandomPad(
    std::vector<RnsPolynomial<ModularInt>> pads,
    const RnsRlweSecretKey<ModularInt>& secret_key, int power, int variance,
    const RnsGadget<ModularInt>* gadget, absl::string_view prng_pad_seed,
    PrngType prng_type, Integer error_scalar) {
  int log_n = secret_key.LogN();
  int num_coeffs = 1 << log_n;
  if (power < 0 || (power % 2) == 0 || power >= 2 * num_coeffs) {
    return absl::InvalidArgumentError(
        absl::StrCat("`power` must be a non-negative odd "
                     "integer less than 2*n."));
  }
  if (variance <= 0) {
    return absl::InvalidArgumentError("`variance` must be positive.");
  }
  if (gadget == nullptr) {
    return absl::InvalidArgumentError("`gadget` must not be null.");
  }

  // Create the PRNGs for sampling the encryption randomness.
  std::unique_ptr<SecurePrng> prng_encryption;
  std::string prng_encryption_seed;
  if (prng_type == PRNG_TYPE_HKDF) {
    RLWE_ASSIGN_OR_RETURN(prng_encryption_seed,
                          SingleThreadHkdfPrng::GenerateSeed());
    RLWE_ASSIGN_OR_RETURN(prng_encryption,
                          SingleThreadHkdfPrng::Create(prng_encryption_seed));
  } else if (prng_type == PRNG_TYPE_CHACHA) {
    RLWE_ASSIGN_OR_RETURN(prng_encryption_seed,
                          SingleThreadChaChaPrng::GenerateSeed());
    RLWE_ASSIGN_OR_RETURN(prng_encryption,
                          SingleThreadChaChaPrng::Create(prng_encryption_seed));
  } else {
    return absl::InvalidArgumentError("PrngType not specified correctly.");
  }

  // The source key s(X^power), which is the secret polynomial under which the
  // source ciphertexts are encrypted.
  auto moduli = secret_key.Moduli();
  const RnsPolynomial<ModularInt>& target_key = secret_key.Key();
  RLWE_ASSIGN_OR_RETURN(auto source_key, target_key.Substitute(power, moduli));

  // The galois key is a k * 2 matrix [[b1, a1], ..., [bk, ak]], where
  // ai is uniformly random and bi = -ai * s + t * ei + gi * s(X^power), for
  // gi = i'th component of the gadget and t = `error_scalar`.
  int k = gadget->Dimension();
  std::vector<RnsPolynomial<ModularInt>> key_bs;
  key_bs.reserve(k);
  for (int i = 0; i < k; ++i) {
    RnsPolynomial<ModularInt> secret = source_key;
    RLWE_RETURN_IF_ERROR(secret.MulInPlace(gadget->Component(i), moduli));

    // a = -u
    RLWE_ASSIGN_OR_RETURN(RnsPolynomial<ModularInt> u, pads[i].Negate(moduli));

    RLWE_ASSIGN_OR_RETURN(
        RnsPolynomial<ModularInt> b,
        SampleError<ModularInt>(log_n, variance, moduli,
                                prng_encryption.get()));       // b = e
    RLWE_RETURN_IF_ERROR(b.MulInPlace(error_scalar, moduli));  // b = t * e
    RLWE_RETURN_IF_ERROR(b.AddInPlace(secret, moduli));  // b = t * e + gi * s'
    RLWE_RETURN_IF_ERROR(b.FusedMulAddInPlace(u, target_key, moduli));

    key_bs.push_back(std::move(b));
  }

  // Store the RNS moduli.
  std::vector<const PrimeModulus<ModularInt>*> moduli_vector;
  moduli_vector.insert(moduli_vector.begin(), moduli.begin(), moduli.end());
  return RnsGaloisKey(/*key_as=*/std::move(pads), std::move(key_bs), gadget,
                      power, std::move(moduli_vector), prng_pad_seed,
                      prng_type);
}

template <typename ModularInt>
absl::StatusOr<RnsGaloisKey<ModularInt>>
RnsGaloisKey<ModularInt>::CreateFromKeyComponents(
    std::vector<RnsPolynomial<ModularInt>> key_as,
    std::vector<RnsPolynomial<ModularInt>> key_bs, int power,
    const RnsGadget<ModularInt>* gadget,
    std::vector<const PrimeModulus<ModularInt>*> moduli,
    absl::string_view prng_seed, PrngType prng_type) {
  if (gadget == nullptr) {
    return absl::InvalidArgumentError("`gadget` must not be null.");
  }
  if (key_as.size() != gadget->Dimension() ||
      key_bs.size() != gadget->Dimension()) {
    return absl::InvalidArgumentError(
        "`key_as` and `key_bs` must have the same length as the gadget "
        "dimension.");
  }
  return RnsGaloisKey(std::move(key_as), std::move(key_bs), gadget, power,
                      std::move(moduli), prng_seed, prng_type);
}

template <typename ModularInt>
absl::StatusOr<RnsGaloisKey<ModularInt>> RnsGaloisKey<ModularInt>::Deserialize(
    const SerializedRnsGaloisKey& serialized,
    const RnsGadget<ModularInt>* gadget,
    std::vector<const PrimeModulus<ModularInt>*> moduli) {
  if (gadget == nullptr) {
    return absl::InvalidArgumentError("`gadget` must not be null.");
  }

  int dimension = serialized.key_bs_size();
  if (dimension <= 0) {
    return absl::InvalidArgumentError("`key_bs` must not be empty.");
  }

  std::vector<RnsPolynomial<ModularInt>> key_bs;
  key_bs.reserve(dimension);
  for (int i = 0; i < dimension; ++i) {
    RLWE_ASSIGN_OR_RETURN(
        RnsPolynomial<ModularInt> key_b,
        RnsPolynomial<ModularInt>::Deserialize(serialized.key_bs(i), moduli));
    key_bs.push_back(std::move(key_b));
  }

  // Sample the random polynomials in `key_as` using the given PRNG seed.
  RLWE_ASSIGN_OR_RETURN(
      std::vector<RnsPolynomial<ModularInt>> key_as,
      SampleRandomPad(dimension, key_bs[0].LogN(), moduli,
                      serialized.prng_seed(), serialized.prng_type()));

  return RnsGaloisKey(std::move(key_as), std::move(key_bs), gadget,
                      serialized.power(), std::move(moduli),
                      serialized.prng_seed(), serialized.prng_type());
}

template <typename ModularInt>
absl::StatusOr<std::vector<RnsPolynomial<ModularInt>>>
RnsGaloisKey<ModularInt>::ApplyToRlweCiphertext(
    const RnsRlweCiphertext<ModularInt>& ciphertext) const {
  if (ciphertext.PowerOfS() != power_) {
    return absl::InvalidArgumentError(
        "`ciphertext` does not have a matching substitution power.");
  }
  if (ciphertext.NumModuli() != moduli_.size()) {
    return absl::InvalidArgumentError(
        "`ciphertext` does not have a matching RNS moduli set.");
  }
  if (ciphertext.Degree() != 1) {
    return absl::InvalidArgumentError(
        "Galois key can only apply to a ciphertext of degree 1.");
  }

  // Apply the Galois key (key_bs, key_as) to a ciphertext (c0, c1) to get
  // a new ciphertext (c0', c1') = (c0, 0) + g^-1(c1) * (key_bs, key_as).
  RLWE_ASSIGN_OR_RETURN(RnsPolynomial<ModularInt> c1, ciphertext.Component(1));
  if (c1.IsNttForm()) {
    RLWE_RETURN_IF_ERROR(c1.ConvertToCoeffForm(moduli_));
  }
  // Gadget decomposition g^-1(c1).
  RLWE_ASSIGN_OR_RETURN(std::vector<RnsPolynomial<ModularInt>> c1_digits,
                        gadget_->Decompose(c1, moduli_));

  int k = gadget_->Dimension();
  RLWE_ASSIGN_OR_RETURN(RnsPolynomial<ModularInt> c0_new,
                        RnsPolynomial<ModularInt>::CreateZero(
                            c1.LogN(), moduli_, /*is_ntt=*/true));
  for (int i = 0; i < k; ++i) {
    RLWE_RETURN_IF_ERROR(c1_digits[i].ConvertToNttForm(moduli_));
    RLWE_RETURN_IF_ERROR(
        c0_new.FusedMulAddInPlace(c1_digits[i], key_bs_[i], moduli_));
  }

  RLWE_ASSIGN_OR_RETURN(RnsPolynomial<ModularInt> c1_new,
                        RnsPolynomial<ModularInt>::CreateZero(
                            c1.LogN(), moduli_, /*is_ntt=*/true));
  for (int i = 0; i < k; ++i) {
    RLWE_RETURN_IF_ERROR(
        c1_new.FusedMulAddInPlace(c1_digits[i], key_as_[i], moduli_));
  }
  RLWE_ASSIGN_OR_RETURN(RnsPolynomial<ModularInt> c0, ciphertext.Component(0));
  RLWE_RETURN_IF_ERROR(c0_new.AddInPlace(c0, moduli_));

  return std::vector<RnsPolynomial<ModularInt>>{std::move(c0_new),
                                                std::move(c1_new)};
}

template <typename ModularInt>
absl::StatusOr<RnsPolynomial<ModularInt>>
RnsGaloisKey<ModularInt>::PrecomputeRandomPad(
    const std::vector<RnsPolynomial<ModularInt>>& ciphertext_pad_digits) const {
  int k = gadget_->Dimension();
  if (ciphertext_pad_digits.size() != k) {
    return absl::InvalidArgumentError(
        "`ciphertext_pad_digits` has incorrect size");
  }

  int log_n = key_as_[0].LogN();
  RLWE_ASSIGN_OR_RETURN(auto c1_new, RnsPolynomial<ModularInt>::CreateZero(
                                         log_n, moduli_, /*is_ntt=*/true));
  for (int i = 0; i < ciphertext_pad_digits.size(); ++i) {
    RLWE_RETURN_IF_ERROR(c1_new.FusedMulAddInPlace(ciphertext_pad_digits[i],
                                                   key_as_[i], moduli_));
  }
  return c1_new;
}

template <typename ModularInt>
absl::StatusOr<RnsPolynomial<ModularInt>>
RnsGaloisKey<ModularInt>::ApplyToRlweCiphertextWithRandomPad(
    const RnsRlweCiphertext<ModularInt>& ciphertext,
    const std::vector<RnsPolynomial<ModularInt>>& ciphertext_pad_digits,
    const RnsPolynomial<ModularInt>& ciphertext_pad_new) const {
  if (ciphertext.PowerOfS() != power_) {
    return absl::InvalidArgumentError(
        "`ciphertext` does not have a matching substitution power.");
  }
  if (ciphertext.NumModuli() != moduli_.size()) {
    return absl::InvalidArgumentError(
        "`ciphertext` does not have a matching RNS moduli set.");
  }
  if (ciphertext.Degree() != 1) {
    return absl::InvalidArgumentError(
        "Galois key can only apply to a ciphertext of degree 1.");
  }
  int k = gadget_->Dimension();
  if (ciphertext_pad_digits.size() != k) {
    return absl::InvalidArgumentError(
        "`ciphertext_pad_digits` has incorrect size");
  }

  RLWE_ASSIGN_OR_RETURN(
      auto lazy_c0_new,
      LazyRnsPolynomial<ModularInt>::CreateZero(key_bs_[0].LogN(), moduli_));
  for (int i = 0; i < k; ++i) {
    CHECK(ciphertext_pad_digits[i].IsNttForm())
        << "`ciphertext_pad_digits` must be in NTT form";
    RLWE_RETURN_IF_ERROR(lazy_c0_new.FusedMulAddInPlace(
        ciphertext_pad_digits[i], key_bs_[i], moduli_));
  }
  RLWE_ASSIGN_OR_RETURN(RnsPolynomial<ModularInt> c0_new,
                        lazy_c0_new.Export(moduli_));
  RLWE_ASSIGN_OR_RETURN(RnsPolynomial<ModularInt> c0, ciphertext.Component(0));
  RLWE_RETURN_IF_ERROR(c0_new.AddInPlace(c0, moduli_));
  return c0_new;
}

template class RnsGaloisKey<MontgomeryInt<Uint16>>;
template class RnsGaloisKey<MontgomeryInt<Uint32>>;
template class RnsGaloisKey<MontgomeryInt<Uint64>>;
template class RnsGaloisKey<MontgomeryInt<absl::uint128>>;
#ifdef ABSL_HAVE_INTRINSIC_INT128
template class RnsGaloisKey<MontgomeryInt<unsigned __int128>>;
#endif

}  // namespace rlwe
