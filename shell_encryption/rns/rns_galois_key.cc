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

#include "absl/numeric/int128.h"
#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "absl/strings/str_cat.h"
#include "shell_encryption/integral_types.h"
#include "shell_encryption/montgomery.h"
#include "shell_encryption/prng/prng.h"
#include "shell_encryption/prng/single_thread_chacha_prng.h"
#include "shell_encryption/prng/single_thread_hkdf_prng.h"
#include "shell_encryption/rns/error_distribution.h"
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

  // The source key s(X^power), which is the secret polynomial under which the
  // source ciphertexts are encrypted.
  auto moduli = secret_key.Moduli();
  const RnsPolynomial<ModularInt>& target_key = secret_key.Key();
  RLWE_ASSIGN_OR_RETURN(auto source_key, target_key.Substitute(power, moduli));

  // Create the PRNGs for sampling the random polynomials in `key_as` and for
  // sampling the encryption randomness.
  std::unique_ptr<SecurePrng> prng_pad, prng_encryption;
  std::string prng_pad_seed, prng_encryption_seed;
  if (prng_type == PRNG_TYPE_HKDF) {
    RLWE_ASSIGN_OR_RETURN(prng_pad_seed, SingleThreadHkdfPrng::GenerateSeed());
    RLWE_ASSIGN_OR_RETURN(prng_pad,
                          SingleThreadHkdfPrng::Create(prng_pad_seed));
    RLWE_ASSIGN_OR_RETURN(prng_encryption_seed,
                          SingleThreadHkdfPrng::GenerateSeed());
    RLWE_ASSIGN_OR_RETURN(prng_encryption,
                          SingleThreadHkdfPrng::Create(prng_encryption_seed));
  } else if (prng_type == PRNG_TYPE_CHACHA) {
    RLWE_ASSIGN_OR_RETURN(prng_pad_seed,
                          SingleThreadChaChaPrng::GenerateSeed());
    RLWE_ASSIGN_OR_RETURN(prng_pad,
                          SingleThreadChaChaPrng::Create(prng_pad_seed));
    RLWE_ASSIGN_OR_RETURN(prng_encryption_seed,
                          SingleThreadChaChaPrng::GenerateSeed());
    RLWE_ASSIGN_OR_RETURN(prng_encryption,
                          SingleThreadChaChaPrng::Create(prng_encryption_seed));
  } else {
    return absl::InvalidArgumentError("PrngType not specified correctly.");
  }

  // The Galois key is a k x 2 matrix [[b1, a1], ..., [bk, ak]], where
  // ai is uniformly random and bi = -ai * s + t * ei + gi * s(X^power), for
  // gi = i'th component of the gadget and t = `error_scalar`. That is,
  // conceptually, (bi, ai) is an encryption of the source secret s(X^power)
  // under the target secret s(X). Note that t should be the plaintext modulus
  // for BGV, and it should be 1 for other schemes.
  int k = gadget->Dimension();
  std::vector<RnsPolynomial<ModularInt>> key_as;
  std::vector<RnsPolynomial<ModularInt>> key_bs;
  key_as.reserve(k);
  key_bs.reserve(k);
  for (int i = 0; i < k; ++i) {
    // Sample a uniformly random polynomial u.
    RLWE_ASSIGN_OR_RETURN(auto u, RnsPolynomial<ModularInt>::SampleUniform(
                                      log_n, prng_pad.get(), moduli));

    // gi * s(X^power)
    RnsPolynomial<ModularInt> secret = source_key;
    RLWE_RETURN_IF_ERROR(secret.MulInPlace(gadget->Component(i), moduli));

    // bi = u * s(X) + t * ei + gi * s(X^power)
    RLWE_ASSIGN_OR_RETURN(RnsPolynomial<ModularInt> b,
                          SampleError<ModularInt>(log_n, variance, moduli,
                                                  prng_encryption.get()));
    RLWE_RETURN_IF_ERROR(b.MulInPlace(error_scalar, moduli));
    RLWE_RETURN_IF_ERROR(b.AddInPlace(secret, moduli));
    RLWE_RETURN_IF_ERROR(b.FusedMulAddInPlace(u, target_key, moduli));

    RLWE_RETURN_IF_ERROR(u.NegateInPlace(moduli));  // a = -u
    key_as.push_back(std::move(u));
    key_bs.push_back(std::move(b));
  }

  // Store the RNS moduli.
  std::vector<const PrimeModulus<ModularInt>*> moduli_vector;
  moduli_vector.insert(moduli_vector.begin(), moduli.begin(), moduli.end());

  return RnsGaloisKey(std::move(key_as), std::move(key_bs), gadget, power,
                      std::move(moduli_vector), prng_pad_seed, prng_type);
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
absl::StatusOr<RnsBgvCiphertext<ModularInt>> RnsGaloisKey<ModularInt>::ApplyTo(
    const RnsBgvCiphertext<ModularInt>& ciphertext) const {
  // First, key-switch the ciphertext as a generic RLWE ciphertext.
  RLWE_ASSIGN_OR_RETURN(std::vector<RnsPolynomial<ModularInt>> components_new,
                        ApplyToRlweCiphertext(ciphertext));

  // Compute the error bound in the new ciphertext.
  double error =
      ciphertext.Error() +
      ciphertext.ErrorParams()->BoundOnGadgetBasedKeySwitching(
          /*num_components=*/1, gadget_->LogGadgetBase(), gadget_->Dimension());

  return RnsBgvCiphertext<ModularInt>(std::move(components_new), moduli_,
                                      /*power=*/1, error,
                                      ciphertext.ErrorParams());
}

template <typename ModularInt>
absl::StatusOr<RnsBfvCiphertext<ModularInt>> RnsGaloisKey<ModularInt>::ApplyTo(
    const RnsBfvCiphertext<ModularInt>& ciphertext) const {
  // First, key-switch the ciphertext as a generic RLWE ciphertext.
  RLWE_ASSIGN_OR_RETURN(std::vector<RnsPolynomial<ModularInt>> components_new,
                        ApplyToRlweCiphertext(ciphertext));

  // Compute the error bound in the new ciphertext.
  double error =
      ciphertext.Error() +
      ciphertext.ErrorParams()->BoundOnGadgetBasedKeySwitching(
          /*num_components=*/1, gadget_->LogGadgetBase(), gadget_->Dimension());

  return RnsBfvCiphertext<ModularInt>(
      std::move(components_new), moduli_,
      /*power=*/1, error, ciphertext.ErrorParams(), ciphertext.Context());
}

template class RnsGaloisKey<MontgomeryInt<Uint16>>;
template class RnsGaloisKey<MontgomeryInt<Uint32>>;
template class RnsGaloisKey<MontgomeryInt<Uint64>>;
template class RnsGaloisKey<MontgomeryInt<absl::uint128>>;
#ifdef ABSL_HAVE_INTRINSIC_INT128
template class RnsGaloisKey<MontgomeryInt<unsigned __int128>>;
#endif

}  // namespace rlwe
