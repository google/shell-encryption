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

#include "shell_encryption/aux_relinearization_key.h"

#include <memory>

#include "absl/numeric/int128.h"
#include "shell_encryption/modulus_conversion.h"
#include "shell_encryption/montgomery.h"
#include "shell_encryption/prng/prng.h"
#include "shell_encryption/prng/single_thread_chacha_prng.h"
#include "shell_encryption/prng/single_thread_hkdf_prng.h"
#include "shell_encryption/serialization.pb.h"
#include "shell_encryption/status_macros.h"
#include "shell_encryption/statusor.h"

namespace rlwe {

template <typename ModularInt>
rlwe::StatusOr<typename AuxModRelinearizationKey<ModularInt>::KeyComponent>
AuxModRelinearizationKey<ModularInt>::KeyComponent::Create(
    const rlwe::Polynomial<ModularInt>& s_target_main,
    const rlwe::Polynomial<ModularInt>& s_target_aux,
    const rlwe::Polynomial<ModularInt>& s_source_main,
    const typename ModularInt::Params& mod_params_main,
    const NttParameters<ModularInt>& ntt_params_main,
    const typename ModularInt::Params& mod_params_aux,
    const NttParameters<ModularInt>& ntt_params_aux,
    const ModularInt& p_mod_main, const ModularInt& t_mod_main,
    const ModularInt& t_mod_aux, int num_coeffs, Uint64 variance,
    SecurePrng& prng, SecurePrng& prng_error) {
  // Sample random a for both (mod q) and (mod p)
  RLWE_ASSIGN_OR_RETURN(Polynomial<ModularInt> a_main,
                        rlwe::SamplePolynomialFromPrng<ModularInt>(
                            num_coeffs, &prng, &mod_params_main));
  RLWE_ASSIGN_OR_RETURN(Polynomial<ModularInt> a_aux,
                        rlwe::SamplePolynomialFromPrng<ModularInt>(
                            num_coeffs, &prng, &mod_params_aux));

  // Sample error term e (mod q) and convert it to mod-p
  RLWE_ASSIGN_OR_RETURN(
      std::vector<ModularInt> e_coeffs_main,
      rlwe::SampleFromErrorDistribution<ModularInt>(
          num_coeffs, variance, &prng_error, &mod_params_main));
  RLWE_ASSIGN_OR_RETURN(
      std::vector<ModularInt> e_coeffs_aux,
      ConvertModulusBalanced(e_coeffs_main, mod_params_main, mod_params_aux));

  // b = t * e (mod q) and (mod p)
  auto b_main = Polynomial<ModularInt>::ConvertToNtt(
      std::move(e_coeffs_main), &ntt_params_main, &mod_params_main);
  auto b_aux = Polynomial<ModularInt>::ConvertToNtt(
      std::move(e_coeffs_aux), &ntt_params_aux, &mod_params_aux);
  RLWE_RETURN_IF_ERROR(b_main.MulInPlace(t_mod_main, &mod_params_main));
  RLWE_RETURN_IF_ERROR(b_aux.MulInPlace(t_mod_aux, &mod_params_aux));

  // b = a * s' + t * e (mod q) and (mod p)
  RLWE_RETURN_IF_ERROR(
      b_main.FusedMulAddInPlace(a_main, s_target_main, &mod_params_main));
  RLWE_RETURN_IF_ERROR(
      b_aux.FusedMulAddInPlace(a_aux, s_target_aux, &mod_params_aux));

  // b = p * s^j + a * s' + t * e (mod q)
  RLWE_RETURN_IF_ERROR(b_main.AddInPlace(s_source_main, &mod_params_main));

  // Negate the random part "a" since we want each pair (b, a) to be
  // an encryption (under the target key s') of the source key power s^j
  a_main.NegateInPlace(&mod_params_main);
  a_aux.NegateInPlace(&mod_params_aux);

  return KeyComponent{std::move(b_main), std::move(a_main), std::move(b_aux),
                      std::move(a_aux)};
}

template <typename ModularInt>
rlwe::StatusOr<AuxModRelinearizationKey<ModularInt>>
AuxModRelinearizationKey<ModularInt>::Create(
    const SymmetricRlweKey<ModularInt>& secret_key, PrngType prng_type,
    int degree, const ModularIntParams* mod_params_aux,
    const NttParameters<ModularInt>* ntt_params_aux, int substitution_power) {
  if (degree <= 0) {
    return absl::InvalidArgumentError(
        absl::StrCat("degree: ", degree, " must be positive."));
  }
  if (mod_params_aux == nullptr) {
    return absl::InvalidArgumentError("mod_params_aux must not be null.");
  }
  if (ntt_params_aux == nullptr) {
    return absl::InvalidArgumentError("ntt_params_aux must not be null.");
  }
  if (substitution_power <= 0) {
    return absl::InvalidArgumentError(absl::StrCat(
        "substitution_power: ", substitution_power, " must be positive."));
  }

  // Generate the PRNGs used to sample the uniform `a` part and the error terms
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
    return absl::InvalidArgumentError(absl::StrCat(
        "PrngType not specified correctly: ", PrngType_Name(prng_type)));
  }

  // The main modulus q and its parameters
  const ModularIntParams* mod_params_main = secret_key.ModulusParams();
  const NttParameters<ModularInt>* ntt_params_main = secret_key.NttParams();
  if (mod_params_main == nullptr) {
    return absl::InvalidArgumentError(
        "`secret_key` must not have a null ModulusParams");
  }
  if (ntt_params_main == nullptr) {
    return absl::InvalidArgumentError(
        "`secret_key` must not have a null NttParams");
  }

  // Compute p (mod q), where p = auxiliary modulus
  RLWE_ASSIGN_OR_RETURN(
      ModularInt p_mod_main,
      ModularInt::ImportInt(mod_params_aux->modulus, mod_params_main));
  // The plaintext modulus t (mod q) and (mod p)
  ModularInt t_mod_main = secret_key.PlaintextModulus();
  typename ModularInt::Int t = t_mod_main.ExportInt(mod_params_main);
  RLWE_ASSIGN_OR_RETURN(ModularInt t_mod_aux,
                        ModularInt::ImportInt(t, mod_params_aux));

  // The target secret key s' (mod q) and (mod p)
  const Polynomial<ModularInt>& ss_main = secret_key.Key();
  RLWE_ASSIGN_OR_RETURN(const Polynomial<ModularInt> ss_aux,
                        ConvertModulusBalancedOnNttPolynomial(
                            ss_main, *mod_params_main, *ntt_params_main,
                            *mod_params_aux, *ntt_params_aux));

  // The source secret key s (mod q). No need to compute s (mod p) since
  // it vanishes when scaling up by p.
  RLWE_ASSIGN_OR_RETURN(
      Polynomial<ModularInt> s_base_main,
      ss_main.Substitute(substitution_power, ntt_params_main, mod_params_main));
  // p * s (mod q)
  RLWE_ASSIGN_OR_RETURN(Polynomial<ModularInt> s_power_main,
                        s_base_main.Mul(p_mod_main, mod_params_main));

  // We store `degree` many key components that correspond to powers s^1,...,s^k
  // of the source secret key, unless substitution_power is 1, in which case we
  // ignore the first power s^1 since it is the same as the target key s'.
  const bool has_identical_base_key = substitution_power == 1;
  const int key_length = has_identical_base_key ? (degree - 1) : degree;
  std::vector<KeyComponent> keys;
  keys.reserve(key_length);
  for (int i = 0; i < key_length; ++i) {
    if (has_identical_base_key || i > 0) {
      // p * s^j, where j = has_identical_base_key ? i+1 : i
      RLWE_RETURN_IF_ERROR(
          s_power_main.MulInPlace(s_base_main, mod_params_main));
    }
    RLWE_ASSIGN_OR_RETURN(
        auto key_component,
        KeyComponent::Create(
            ss_main, ss_aux, s_power_main, *mod_params_main, *ntt_params_main,
            *mod_params_aux, *ntt_params_aux, p_mod_main, t_mod_main, t_mod_aux,
            secret_key.Len(), secret_key.Variance(), *prng, *prng_error));
    keys.push_back(std::move(key_component));
  }
  return AuxModRelinearizationKey(
      mod_params_main, ntt_params_main, mod_params_aux, ntt_params_aux,
      prng_seed, prng_type, std::move(keys), t, substitution_power);
}

template <typename ModularInt>
rlwe::StatusOr<rlwe::SymmetricRlweCiphertext<ModularInt>>
AuxModRelinearizationKey<ModularInt>::ApplyTo(
    const SymmetricRlweCiphertext<ModularInt>& ciphertext) const {
  if (ciphertext.ModulusParams()->modulus != mod_params_main_->modulus) {
    return absl::InvalidArgumentError(
        "Ciphertext modulus does not match key-switching key modulus");
  }
  // If this key is generated for substitution_power == 1, then we store l - 2
  // key components where l is the ciphertext length; otherwise we store l - 1
  // key components.
  const bool has_identical_base_key = substitution_power_ == 1;
  const int expected_num_components =
      has_identical_base_key ? keys_.size() + 2 : keys_.size() + 1;
  if (ciphertext.Len() != expected_num_components) {
    return absl::InvalidArgumentError(
        "Ciphertext has incompatible number of components");
  }
  if (ciphertext.ErrorParams() == nullptr) {
    return absl::InvalidArgumentError(
        "Ciphertext must not have a null ErrorParams");
  }

  RLWE_ASSIGN_OR_RETURN(
      ModularInt p_mod_q,
      ModularInt::ImportInt(mod_params_aux_->modulus, mod_params_main_));

  // Initialize the result ciphertext (c0', c1') (mod q).
  // We always have c0' = p * c0. For c1', if substitution_power == 1, then
  // we know s' == s and so we set c0' = p * c1 as an optimization; otherwise,
  // we have c0' = 0.
  RLWE_ASSIGN_OR_RETURN(Polynomial<ModularInt> c0, ciphertext.Component(0));
  int n = ntt_params_main_->number_coeffs;
  Polynomial<ModularInt> c1(n, mod_params_main_);
  if (has_identical_base_key) {
    RLWE_ASSIGN_OR_RETURN(Polynomial<ModularInt> c, ciphertext.Component(1));
    RLWE_RETURN_IF_ERROR(c1.AddInPlace(c, mod_params_main_));
  }
  RLWE_RETURN_IF_ERROR(c0.MulInPlace(p_mod_q, mod_params_main_));
  RLWE_RETURN_IF_ERROR(c1.MulInPlace(p_mod_q, mod_params_main_));

  // Initialize (c0', c1') (mod p) to (0, 0) for the auxiliary modulus p.
  // No need to convert from (c0', c1') (mod q) since c0' and c1' are both
  // divisible by p now.
  Polynomial<ModularInt> c0_aux(n, mod_params_aux_);
  Polynomial<ModularInt> c1_aux(n, mod_params_aux_);

  // Add <key, (c[first],..,c[degree])> to (c0', c1') for both mod-q and mod-p
  int first_key_power = has_identical_base_key ? 2 : 1;
  for (int i = first_key_power; i < ciphertext.Len(); ++i) {
    int key_index = i - first_key_power;
    RLWE_ASSIGN_OR_RETURN(Polynomial<ModularInt> ci, ciphertext.Component(i));
    // c0' (mod q) = p * c0 + sum(key[idx].b * ci)
    RLWE_RETURN_IF_ERROR(
        c0.FusedMulAddInPlace(keys_[key_index].b_main, ci, mod_params_main_));
    // c1' (mod q) = p * c1 + sum(key[idx].a * ci)
    RLWE_RETURN_IF_ERROR(
        c1.FusedMulAddInPlace(keys_[key_index].a_main, ci, mod_params_main_));

    RLWE_ASSIGN_OR_RETURN(Polynomial<ModularInt> ci_aux,
                          ConvertModulusBalancedOnNttPolynomial(
                              ci, *mod_params_main_, *ntt_params_main_,
                              *mod_params_aux_, *ntt_params_aux_));
    // c0' (mod p) = sum(key[idx].b * ci)
    RLWE_RETURN_IF_ERROR(c0_aux.FusedMulAddInPlace(keys_[key_index].b_aux,
                                                   ci_aux, mod_params_aux_));
    // c1' (mod p) = sum(key[idx].a * ci)
    RLWE_RETURN_IF_ERROR(c1_aux.FusedMulAddInPlace(keys_[key_index].a_aux,
                                                   ci_aux, mod_params_aux_));
  }

  // Convert the ciphertext (c0, c1) from mod-q*p to mod-q.
  RLWE_RETURN_IF_ERROR(ConvertToMainModulus(p_mod_q, std::move(c0_aux),
                                            std::move(c1_aux), c0, c1));

  // Update error in the new ciphertext
  double error = ciphertext.Error() +
                 ciphertext.ErrorParams()->B_aux_mod_relinearize(
                     ciphertext.Len() - first_key_power, *mod_params_aux_);

  return SymmetricRlweCiphertext<ModularInt>(
      {std::move(c0), std::move(c1)},
      /*power_of_s=*/1, error, mod_params_main_, ciphertext.ErrorParams());
}

template <typename ModularInt>
absl::Status AuxModRelinearizationKey<ModularInt>::ConvertToMainModulus(
    const ModularInt& p_mod_main, rlwe::Polynomial<ModularInt> c0_aux,
    rlwe::Polynomial<ModularInt> c1_aux, rlwe::Polynomial<ModularInt>& c0_main,
    rlwe::Polynomial<ModularInt>& c1_main) const {
  // Multiply [c0']_p and [c1']_p by t^(-1) (mod p) before converting to mod-p
  RLWE_ASSIGN_OR_RETURN(auto t_mod_p,
                        ModularInt::ImportInt(t_, mod_params_aux_));
  ModularInt t_inv_mod_p = t_mod_p.MultiplicativeInverse(mod_params_aux_);
  RLWE_RETURN_IF_ERROR(c0_aux.MulInPlace(t_inv_mod_p, mod_params_aux_));
  RLWE_RETURN_IF_ERROR(c1_aux.MulInPlace(t_inv_mod_p, mod_params_aux_));

  // Convert [c0']_p and [c1']_p to mod-q
  RLWE_ASSIGN_OR_RETURN(Polynomial<ModularInt> c0_aux_modq,
                        ConvertModulusBalancedOnNttPolynomial(
                            c0_aux, *mod_params_aux_, *ntt_params_aux_,
                            *mod_params_main_, *ntt_params_main_));
  RLWE_ASSIGN_OR_RETURN(Polynomial<ModularInt> c1_aux_modq,
                        ConvertModulusBalancedOnNttPolynomial(
                            c1_aux, *mod_params_aux_, *ntt_params_aux_,
                            *mod_params_main_, *ntt_params_main_));

  // Multiply by t (mod q) after conversion
  RLWE_ASSIGN_OR_RETURN(auto t_mod_q,
                        ModularInt::ImportInt(t_, mod_params_main_));
  RLWE_RETURN_IF_ERROR(c0_aux_modq.MulInPlace(t_mod_q, mod_params_main_));
  RLWE_RETURN_IF_ERROR(c1_aux_modq.MulInPlace(t_mod_q, mod_params_main_));

  // Convert (c0', c1') from mod-q*p to mod-q.
  // [c0']_q = p^(-1) * (c0 - [c0']_p) mod q
  ModularInt p_inv_mod_q = p_mod_main.MultiplicativeInverse(mod_params_main_);
  RLWE_RETURN_IF_ERROR(c0_main.SubInPlace(c0_aux_modq, mod_params_main_));
  RLWE_RETURN_IF_ERROR(c0_main.MulInPlace(p_inv_mod_q, mod_params_main_));

  // [c1']_q = p^(-1) * (c1 - [c1']_p) mod q
  RLWE_RETURN_IF_ERROR(c1_main.SubInPlace(c1_aux_modq, mod_params_main_));
  RLWE_RETURN_IF_ERROR(c1_main.MulInPlace(p_inv_mod_q, mod_params_main_));
  return absl::OkStatus();
}

template <typename ModularInt>
rlwe::StatusOr<SerializedAuxModRelinearizationKey>
AuxModRelinearizationKey<ModularInt>::Serialize() const {
  SerializedAuxModRelinearizationKey output;
  output.set_num_components(keys_.size());
  output.set_prng_seed(prng_seed_);
  output.set_prng_type(prng_type_);
  output.set_power_of_s(substitution_power_);
  // Only serialize the "b" part of the key matrix.
  for (const auto& [b_main, unused_a_main, b_aux, unused_a_aux] : keys_) {
    RLWE_ASSIGN_OR_RETURN(*output.add_b(), b_main.Serialize(mod_params_main_));
    RLWE_ASSIGN_OR_RETURN(*output.add_b(), b_aux.Serialize(mod_params_aux_));
  }
  return output;
}

template <typename ModularInt>
rlwe::StatusOr<AuxModRelinearizationKey<ModularInt>>
AuxModRelinearizationKey<ModularInt>::Deserialize(
    const SerializedAuxModRelinearizationKey& serialized,
    const ModularIntParams* mod_params_main,  // q
    const NttParameters<ModularInt>* ntt_params_main,
    const ModularIntParams* mod_params_aux,  // p
    const NttParameters<ModularInt>* ntt_params_aux,
    typename ModularInt::Int t) {
  // Verifies that the number of components in `serialized` is expected.
  // There should be `num_components` * 2 polynomials in `serialized`, which
  // correspond to the "b" row the key matrix wrt both mod-q and mod-p.
  if (serialized.num_components() <= 0) {
    return absl::InvalidArgumentError(absl::StrCat("The number of components, ",
                                                   serialized.num_components(),
                                                   ", must be positive."));
  } else if (serialized.b_size() != (serialized.num_components() * 2)) {
    return absl::InvalidArgumentError(absl::StrCat(
        "The length of serialized, ", serialized.b_size(), ", ",
        "must be 2 * num_components, ", serialized.num_components(), "."));
  }

  if (mod_params_main == nullptr) {
    return absl::InvalidArgumentError("mod_params_main must not be null.");
  }
  if (ntt_params_main == nullptr) {
    return absl::InvalidArgumentError("ntt_params_main must not be null.");
  }
  if (mod_params_aux == nullptr) {
    return absl::InvalidArgumentError("mod_params_aux must not be null.");
  }
  if (ntt_params_aux == nullptr) {
    return absl::InvalidArgumentError("ntt_params_aux must not be null.");
  }

  // Create prng for sampling the random "a" part based on seed and type.
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

  // Generate the "a" and "b" parts wrt both mod-q and mod-p.
  std::vector<SerializedNttPolynomial> serialized_b(serialized.b().begin(),
                                                    serialized.b().end());
  int key_length = serialized.num_components();
  std::vector<KeyComponent> keys;
  keys.reserve(key_length);
  for (int i = 0; i < key_length; i += 2) {
    RLWE_ASSIGN_OR_RETURN(auto b_main, Polynomial<ModularInt>::Deserialize(
                                           serialized_b[i], mod_params_main));
    RLWE_ASSIGN_OR_RETURN(auto b_aux, Polynomial<ModularInt>::Deserialize(
                                          serialized_b[i + 1], mod_params_aux));

    RLWE_ASSIGN_OR_RETURN(auto a_main, SamplePolynomialFromPrng<ModularInt>(
                                           ntt_params_main->number_coeffs,
                                           prng.get(), mod_params_main));
    a_main.NegateInPlace(mod_params_main);
    RLWE_ASSIGN_OR_RETURN(auto a_aux, SamplePolynomialFromPrng<ModularInt>(
                                          ntt_params_aux->number_coeffs,
                                          prng.get(), mod_params_aux));
    a_aux.NegateInPlace(mod_params_aux);

    keys.push_back(KeyComponent{std::move(b_main), std::move(a_main),
                                std::move(b_aux), std::move(a_aux)});
  }

  return AuxModRelinearizationKey(
      mod_params_main, ntt_params_main, mod_params_aux, ntt_params_aux,
      serialized.prng_seed(), serialized.prng_type(), std::move(keys), t,
      serialized.power_of_s());
}

template class AuxModRelinearizationKey<MontgomeryInt<Uint16>>;
template class AuxModRelinearizationKey<MontgomeryInt<Uint32>>;
template class AuxModRelinearizationKey<MontgomeryInt<Uint64>>;
template class AuxModRelinearizationKey<MontgomeryInt<absl::uint128>>;
#ifdef ABSL_HAVE_INTRINSIC_INT128
template class AuxModRelinearizationKey<MontgomeryInt<unsigned __int128>>;
#endif

}  // namespace rlwe
