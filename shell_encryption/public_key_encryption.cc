/*
 * Copyright 2022 Google LLC.
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

#include "shell_encryption/public_key_encryption.h"

#include <memory>

#include "absl/numeric/int128.h"
#include "absl/status/status.h"
#include "shell_encryption/montgomery.h"
#include "shell_encryption/prng/prng.h"
#include "shell_encryption/prng/single_thread_chacha_prng.h"
#include "shell_encryption/prng/single_thread_hkdf_prng.h"
#include "shell_encryption/sample_error.h"
#include "shell_encryption/serialization.pb.h"
#include "shell_encryption/status_macros.h"
#include "shell_encryption/statusor.h"

namespace rlwe {

template <typename ModularInt>
PublicRlweKey<ModularInt>::PublicRlweKey(
    Polynomial<ModularInt> a, Polynomial<ModularInt> b,
    absl::string_view prng_seed, PrngType prng_type, ModularInt t_mod,
    const NttParameters<ModularInt>* ntt_params,
    const ModularIntParams* mod_params)
    : key_a_(std::move(a)),
      key_b_(std::move(b)),
      prng_seed_(prng_seed),
      prng_type_(prng_type),
      t_mod_(t_mod),
      ntt_params_(*ntt_params),
      mod_params_(*mod_params) {}

template <typename ModularInt>
rlwe::StatusOr<PublicRlweKey<ModularInt>> PublicRlweKey<ModularInt>::Create(
    const SymmetricRlweKey<ModularInt>& secret_key, int variance,
    PrngType prng_type) {
  if (variance <= 0) {
    return absl::InvalidArgumentError(
        absl::StrCat("The variance, ", variance, ", must be positive."));
  }

  // Use the same parameters as in the secret key
  const typename ModularInt::Params* mod_params = secret_key.ModulusParams();
  const NttParameters<ModularInt>* ntt_params = secret_key.NttParams();
  if (mod_params == nullptr) {
    return absl::InvalidArgumentError(
        "Symmetric key must have a valid ModulusParams.");
  }
  if (ntt_params == nullptr) {
    return absl::InvalidArgumentError(
        "Symmetric key must have a valid NttParams.");
  }
  // Use the same plaintext modulus as in the secret key
  ModularInt t_mod = secret_key.PlaintextModulus();
  unsigned int num_coeffs = secret_key.Len();

  // Create two PRNGs: `prng` is used to sample the random polynomial "a" in the
  // public key and its seed is stored such that it can be serialized. The other
  // PRNG `prng_error` is used ro sample the error term and it is not stored.
  std::unique_ptr<SecurePrng> prng, prng_error;
  std::string prng_seed, prng_error_seed;
  if (prng_type == PRNG_TYPE_HKDF) {
    RLWE_ASSIGN_OR_RETURN(prng_seed, SingleThreadHkdfPrng::GenerateSeed());
    RLWE_ASSIGN_OR_RETURN(prng, SingleThreadHkdfPrng::Create(prng_seed));
    RLWE_ASSIGN_OR_RETURN(prng_error_seed,
                          SingleThreadHkdfPrng::GenerateSeed());
    RLWE_ASSIGN_OR_RETURN(prng_error,
                          SingleThreadHkdfPrng::Create(prng_error_seed));
  } else if (prng_type == PRNG_TYPE_CHACHA) {
    RLWE_ASSIGN_OR_RETURN(prng_seed, SingleThreadChaChaPrng::GenerateSeed());
    RLWE_ASSIGN_OR_RETURN(prng, SingleThreadChaChaPrng::Create(prng_seed));
    RLWE_ASSIGN_OR_RETURN(prng_error_seed,
                          SingleThreadChaChaPrng::GenerateSeed());
    RLWE_ASSIGN_OR_RETURN(prng_error,
                          SingleThreadChaChaPrng::Create(prng_error_seed));
  } else {
    return absl::InvalidArgumentError("PrngType not specified correctly.");
  }

  // Sample a from the uniform distribution.
  RLWE_ASSIGN_OR_RETURN(
      Polynomial<ModularInt> a,
      SamplePolynomialFromPrng<ModularInt>(num_coeffs, prng.get(), mod_params));

  // Sample the coefficients for the error term e
  RLWE_ASSIGN_OR_RETURN(
      std::vector<ModularInt> e_coeffs,
      SampleFromErrorDistribution<ModularInt>(num_coeffs, variance,
                                              prng_error.get(), mod_params));
  // b = e in NTT format
  Polynomial<ModularInt> b = Polynomial<ModularInt>::ConvertToNtt(
      std::move(e_coeffs), ntt_params, mod_params);
  // b = e * t
  RLWE_RETURN_IF_ERROR(b.MulInPlace(t_mod, mod_params));
  // b = e * t + a * s
  RLWE_RETURN_IF_ERROR(b.FusedMulAddInPlace(a, secret_key.Key(), mod_params));

  // Construct the public key (b, -a)
  return PublicRlweKey<ModularInt>(a.Negate(mod_params), b, prng_seed,
                                   prng_type, t_mod, ntt_params, mod_params);
}

template <typename ModularInt>
rlwe::StatusOr<SymmetricRlweCiphertext<ModularInt>>
PublicRlweKey<ModularInt>::Encrypt(const Polynomial<ModularInt>& plaintext,
                                   int variance,
                                   const ErrorParams<ModularInt>& error_params,
                                   SecurePrng* prng) const {
  if (variance <= 0) {
    return absl::InvalidArgumentError(
        absl::StrCat("The variance, ", variance, ", must be positive."));
  }
  if (prng == nullptr) {
    return absl::InvalidArgumentError("prng must not be null");
  }

  // Sample the encryption randomness
  RLWE_ASSIGN_OR_RETURN(std::vector<ModularInt> v_coeffs,
                        SampleFromErrorDistribution<ModularInt>(
                            Len(), variance, prng, &mod_params_));
  Polynomial<ModularInt> v = Polynomial<ModularInt>::ConvertToNtt(
      std::move(v_coeffs), &ntt_params_, &mod_params_);

  // Sample the error terms e0, e1
  RLWE_ASSIGN_OR_RETURN(std::vector<ModularInt> e0_coeffs,
                        SampleFromErrorDistribution<ModularInt>(
                            Len(), variance, prng, &mod_params_));
  RLWE_ASSIGN_OR_RETURN(std::vector<ModularInt> e1_coeffs,
                        SampleFromErrorDistribution<ModularInt>(
                            Len(), variance, prng, &mod_params_));

  // Generate the ciphertext polynomials (c0, c1)
  Polynomial<ModularInt> c0 = Polynomial<ModularInt>::ConvertToNtt(
      std::move(e0_coeffs), &ntt_params_, &mod_params_);  // c0 = e0
  RLWE_RETURN_IF_ERROR(c0.MulInPlace(t_mod_,
                                     &mod_params_));  // c0 = e0 * t
  RLWE_RETURN_IF_ERROR(
      c0.FusedMulAddInPlace(key_b_, v, &mod_params_));  // c0 = e0 * t + b * v
  RLWE_RETURN_IF_ERROR(
      c0.AddInPlace(plaintext, &mod_params_));  // c0 = e0 * t + b * v + m

  Polynomial<ModularInt> c1 = Polynomial<ModularInt>::ConvertToNtt(
      std::move(e1_coeffs), &ntt_params_, &mod_params_);  // c1 = e1
  RLWE_RETURN_IF_ERROR(c1.MulInPlace(t_mod_,
                                     &mod_params_));  // c1 = e1 * t
  RLWE_RETURN_IF_ERROR(
      c1.FusedMulAddInPlace(key_a_, v, &mod_params_));  // c1 = e1 * t + a * v

  RLWE_ASSIGN_OR_RETURN(double error,
                        error_params.B_publickey_encryption(Len(), variance));
  return SymmetricRlweCiphertext<ModularInt>(
      std::vector<Polynomial<ModularInt>>{std::move(c0), std::move(c1)}, 1,
      std::move(error), &mod_params_, &error_params);
}

template <typename ModularInt>
rlwe::StatusOr<SerializedPublicRlweKey> PublicRlweKey<ModularInt>::Serialize()
    const {
  SerializedPublicRlweKey output;
  output.set_prng_seed(prng_seed_);
  output.set_prng_type(prng_type_);
  RLWE_ASSIGN_OR_RETURN(*output.mutable_b(), key_b_.Serialize(&mod_params_));
  return output;
}

template <typename ModularInt>
rlwe::StatusOr<PublicRlweKey<ModularInt>>
PublicRlweKey<ModularInt>::Deserialize(
    const SerializedPublicRlweKey& serialized, ModularInt t_mod,
    const ModularIntParams* mod_params,
    const NttParameters<ModularInt>* ntt_params) {
  if (mod_params == nullptr) {
    return absl::InvalidArgumentError("mod_params must not be null.");
  }
  if (ntt_params == nullptr) {
    return absl::InvalidArgumentError("ntt_params must not be null.");
  }

  // Create prng based on seed and type, which is to sample the polynomial "a"
  std::unique_ptr<SecurePrng> prng;
  if (serialized.prng_type() == PRNG_TYPE_HKDF) {
    RLWE_ASSIGN_OR_RETURN(prng,
                          SingleThreadHkdfPrng::Create(serialized.prng_seed()));
  } else if (serialized.prng_type() == PRNG_TYPE_CHACHA) {
    RLWE_ASSIGN_OR_RETURN(
        prng, SingleThreadChaChaPrng::Create(serialized.prng_seed()));
  } else {
    return absl::InvalidArgumentError("Invalid PRNG type is specified.");
  }

  // Generate "a" using the deserialized prng.
  RLWE_ASSIGN_OR_RETURN(Polynomial<ModularInt> a,
                        SamplePolynomialFromPrng<ModularInt>(
                            ntt_params->number_coeffs, prng.get(), mod_params));

  // Deserialize the polynomial "b"
  RLWE_ASSIGN_OR_RETURN(
      Polynomial<ModularInt> b,
      Polynomial<ModularInt>::Deserialize(serialized.b(), mod_params));

  // Construct the public key (b, -a)
  return PublicRlweKey<ModularInt>(
      a.Negate(mod_params), b, serialized.prng_seed(), serialized.prng_type(),
      t_mod, ntt_params, mod_params);
}

// Instantiations of PublicRlweKey with specific MontgomeryInt classes.
// If any new types are added, montgomery.h should be updated accordingly (such
// as ensuring BigInt is correctly specialized, etc.).
template class PublicRlweKey<MontgomeryInt<Uint16>>;
template class PublicRlweKey<MontgomeryInt<Uint32>>;
template class PublicRlweKey<MontgomeryInt<Uint64>>;
template class PublicRlweKey<MontgomeryInt<absl::uint128>>;
#ifdef ABSL_HAVE_INTRINSIC_INT128
template class PublicRlweKey<MontgomeryInt<unsigned __int128>>;
#endif

}  // namespace rlwe
