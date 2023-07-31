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

#ifndef RLWE_AUX_RELINEARIZATION_KEY_H_
#define RLWE_AUX_RELINEARIZATION_KEY_H_

#include <memory>
#include <vector>

#include "shell_encryption/montgomery.h"
#include "shell_encryption/ntt_parameters.h"
#include "shell_encryption/polynomial.h"
#include "shell_encryption/prng/prng.h"
#include "shell_encryption/statusor.h"
#include "shell_encryption/symmetric_encryption.h"

namespace rlwe {

// Implementation of the key-switching techniques using an auxiliary
// modulus as suggested in "Homomorphic Evaluation of the AES Circuit" by
// Craig Gentry, Shai Halevi, and Nigel Smart, https://eprint.iacr.org/2012/099
//
// Assume the ciphertexts are polynomials in the ring R_q = Z[X]/(X^N + 1, q).
// This key-switching key uses an auxiliary modulus p that is co-prime to the
// main modulus q to reduce the error size introduced during the key-switching
// operations. More specifically, to switch from a cipertext ct = (c0,..,ck) of
// degree k wrt the secret key sk = (1,s,..,s^k) to a ciphertext ct' = (c0',c1')
// of degree 1 wrt the secret key sk' = (1,s'), we encrypt the source secret key
// sk under the target secret key sk' to generate a key-switching key rk over
// R_{q*p}:
//
//      (b_0, ..., b_k)
// rk = (-------------) (mod q*p) \in R_{q*p}^{2 \times k+1}
//      (a_0, ..., a_k)
// where for each i:
// - a_i is sampled uniformly random over R_{q*p} and
// - b_i = - a_i * s' + t * e_i + p * s^i (mod q*p), where t is the plaintext
//   modulus.
//
// To apply a key-switching key rk on ciphertext ct_q = (c0_q,..,ck_q) (mod q),
// we first convert it to ct_p = (c0_p,..,ck_p) (mod p); (ct_q, ct_p) represents
// the ciphertext ct over R_{q*p} due to Chinese Reminder Theorem.
// We then compute the matrix-vector product ct' = p^(-1) * rk * ct and convert
// it back to the mod-q ring. The resulting ciphertext ct' encrypts the same
// plaintext as ct but under a different secret key sk'.
//
// Note that the first component of sk and sk' are both 1, so we don't need to
// store the first column in rk. Furthermore, when s' = s, we can also ignore
// the second column of sk, which is the case of the ``relinearization'' key
// that can bring a ciphertext of length k + 1 due to multiplications to a
// ciphertext of standard length 2.
//
// Note that the matrix-vector product rk * ct (mod q*p) is almost an encryption
// of the plaintext as in ct but with a large error sum(e_i * ci_q). In order to
// reduce the error size we then multiply it with p^(-1). So we must choose the
// auxiliary modulus that is sufficiently large. In general p can be roughly the
// size of sqrt(N)*q.
//
// In addition, we now should consider LWE parameters (N,q*p,sigma) for concrete
// security, as this relinearization key is RLWE encryptions wrt modulus q*p.
template <typename ModularInt>
class AuxModRelinearizationKey {
 public:
  using ModularIntParams = typename ModularInt::Params;

  // AuxModRelinearizationKey is copyable and movable.
  AuxModRelinearizationKey(const AuxModRelinearizationKey&) = default;
  AuxModRelinearizationKey& operator=(const AuxModRelinearizationKey&) =
      default;
  AuxModRelinearizationKey(AuxModRelinearizationKey&&) = default;
  AuxModRelinearizationKey& operator=(AuxModRelinearizationKey&&) = default;
  ~AuxModRelinearizationKey() = default;

  // Generates a AuxModRelinearizationKey that can switch a ciphertext encrypted
  // under a secret key sk to a ciphertext encrypted under `secret_key` (1, s').
  //
  // We consider two cases: 1) `substitution_power` = 1, in which case we assume
  // sk = (1, s, ..., s^k) for k = `degree`; 2) `substitution_power` != 1, in
  // which case we assume sk = (1, s(X^substitution_power)) and the generated
  // key is a Galois key.
  static StatusOr<AuxModRelinearizationKey> Create(
      const SymmetricRlweKey<ModularInt>& secret_key, PrngType prng_type,
      int degree, const ModularIntParams* mod_params_aux,
      const NttParameters<ModularInt>* ntt_params_aux,
      int substitution_power = 1);

  // Applies the relinearization key to a ciphertext.
  StatusOr<SymmetricRlweCiphertext<ModularInt>> ApplyTo(
      const SymmetricRlweCiphertext<ModularInt>& ciphertext) const;

  // Returns a SerializedAuxModRelinearizationKey containing a flattened
  // representation of the SerializedNttPolynomials in the key, the
  // substitution_power, and the number of components the key is comprised of.
  StatusOr<SerializedAuxModRelinearizationKey> Serialize() const;

  // Returns an AuxModRelinearizationKey represented as in `serialized`, which
  // should contain `num_components` * 2 polynomials in serialized form for the
  // "b" part of the key wrt the main and the auxiliary moduli,
  // and also the PRNG seed and type used to generate the "a" part polynomials.
  static StatusOr<AuxModRelinearizationKey> Deserialize(
      const SerializedAuxModRelinearizationKey& serialized,
      const ModularIntParams* mod_params_main,
      const NttParameters<ModularInt>* ntt_params_main,
      const ModularIntParams* mod_params_aux,
      const NttParameters<ModularInt>* ntt_params_aux,
      typename ModularInt::Int t);

  // Accessors.
  // Returns the "degree" of a ciphertext ct that this relinearization key can
  // be applied to, which is the highest power k of its underlying secret s,
  // i.e. ct can be decrypted as <ct, (1, s, ..., s^k)>.
  int Degree() const {
    if (substitution_power_ == 1) {
      return keys_.size() + 1;
    }
    return keys_.size();
  }

  int SubstitutionPower() const { return substitution_power_; }

 private:
  // This struct holds a component of a relinearization key, i.e. a column of
  // the relinearization key matrix, with respect to both the main and the
  // auxiliary moduli.
  struct KeyComponent {
    Polynomial<ModularInt> b_main;
    Polynomial<ModularInt> a_main;
    Polynomial<ModularInt> b_aux;
    Polynomial<ModularInt> a_aux;

    // Returns the relinearization key component (i.e. a column in the key
    // matrix) corresponding to the given source key power s in
    // `s_source_main`, wrt both the main modulus (q) and the
    // auxiliary modulus (p). A key component is an encryption of p*s under the
    // target secret key s' given in `s_target_main` and `s_target_aux`
    // wrt both moduli, and so it is represented by a quadruple (b_main,
    // a_main, b_aux, a_aux). The parameter `p_mod_main` stores p
    // (mod q), `t_mod_main` and `t_mod_aux` store t (mod q) and (mod p),
    // where t is the plaintext modulus, `prng` and `prng_error` are PRNGs for
    // sampling the "a" and the error terms, resp. Since p*s (mod p) == 0, we do
    // not need the source key s (mod p).
    static StatusOr<KeyComponent> Create(
        const Polynomial<ModularInt>& s_target_main,
        const Polynomial<ModularInt>& s_target_aux,
        const Polynomial<ModularInt>& s_source_main,
        const typename ModularInt::Params& mod_params_main,
        const NttParameters<ModularInt>& ntt_params_main,
        const typename ModularInt::Params& mod_params_aux,
        const NttParameters<ModularInt>& ntt_params_aux,
        const ModularInt& p_mod_main, const ModularInt& t_mod_main,
        const ModularInt& t_mod_aux, int num_coeffs, Uint64 variance,
        SecurePrng& prng, SecurePrng& prng_error);
  };

  // Constructor.
  explicit AuxModRelinearizationKey(
      const ModularIntParams* mod_params_main,
      const NttParameters<ModularInt>* ntt_params_main,
      const ModularIntParams* mod_params_aux,
      const NttParameters<ModularInt>* ntt_params_aux,
      absl::string_view prng_seed, PrngType prng_type,
      std::vector<KeyComponent> keys, typename ModularInt::Int t,
      int substitution_power)
      : mod_params_main_(mod_params_main),
        ntt_params_main_(ntt_params_main),
        mod_params_aux_(mod_params_aux),
        ntt_params_aux_(ntt_params_aux),
        prng_seed_(prng_seed),
        prng_type_(prng_type),
        keys_(std::move(keys)),
        t_(t),
        substitution_power_(substitution_power) {}

  // Converts a ciphertext (c0, c1) wrt modulus q*p with auxiliary modulus p,
  // represented by the four input polynomials `c0_main`, `c1_main`,
  // `c0_aux`, and `c1_aux`, to a ciphertext encrypting the same plaintext under
  // the main modulus q. The result ciphertext is modified in-place over
  // the input `c0_main` and `c1_main`.
  absl::Status ConvertToMainModulus(const ModularInt& p_mod_main,
                                    Polynomial<ModularInt> c0_aux,
                                    Polynomial<ModularInt> c1_aux,
                                    Polynomial<ModularInt>& c0_main,
                                    Polynomial<ModularInt>& c1_main) const;

  // Modulus and NTT parameters for the main modulus.
  const ModularIntParams* mod_params_main_;
  const NttParameters<ModularInt>* ntt_params_main_;

  // Modulus and NTT parameters for the auxiliary modulus.
  const ModularIntParams* mod_params_aux_;
  const NttParameters<ModularInt>* ntt_params_aux_;

  // Prng seed and type for sampling the random parts key_a_main_ and
  // key_a_aux_.
  std::string prng_seed_;
  PrngType prng_type_;

  // The columns of the relinearization key matrix, where each column consists
  // of two entries represented by four polynomials (wrt both the main and the
  // auxiliary modulus) b_mod_main, a_mod_main, b_mod_aux, a_mod_aux. There are
  // "degree" many columns if substitution_power != 1, and degree - 1 otherwise.
  std::vector<KeyComponent> keys_;

  // Plaintext modulus.
  typename ModularInt::Int t_;

  // Substitution power.
  int substitution_power_;
};

extern template class AuxModRelinearizationKey<MontgomeryInt<Uint16>>;
extern template class AuxModRelinearizationKey<MontgomeryInt<Uint32>>;
extern template class AuxModRelinearizationKey<MontgomeryInt<Uint64>>;
extern template class AuxModRelinearizationKey<MontgomeryInt<absl::uint128>>;
#ifdef ABSL_HAVE_INTRINSIC_INT128
extern template class AuxModRelinearizationKey<
    MontgomeryInt<unsigned __int128>>;
#endif

}  // namespace rlwe

#endif  // RLWE_AUX_RELINEARIZATION_KEY_H_
