/*
 * Copyright 2024 Google LLC.
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

#ifndef RLWE_RNS_RNS_RELINEARIZATION_KEY_H_
#define RLWE_RNS_RNS_RELINEARIZATION_KEY_H_

#include <algorithm>
#include <string>
#include <vector>

#include "absl/status/statusor.h"
#include "absl/strings/string_view.h"
#include "absl/types/span.h"
#include "shell_encryption/rns/rns_bfv_ciphertext.h"
#include "shell_encryption/rns/rns_bgv_ciphertext.h"
#include "shell_encryption/rns/rns_ciphertext.h"
#include "shell_encryption/rns/rns_gadget.h"
#include "shell_encryption/rns/rns_modulus.h"
#include "shell_encryption/rns/rns_polynomial.h"
#include "shell_encryption/rns/rns_secret_key.h"
#include "shell_encryption/rns/serialization.pb.h"
#include "shell_encryption/status_macros.h"

namespace rlwe {

// This class implements the relinearization key for RNS variants of RLWE
// homomorphic encryption schemes.
//
// The current RnsRelinKey class implements the gadget-based relinearization key
// in the power-of-2 cyclotomic ring R = Z[X]/(X^N + 1) with RNS modulus Q. A
// gadget-based relinearization key is a k*(l-1)-by-2 matrix (gk_a, gk_b):
//   gk_a = -u \in R_Q^(k * (l-1)),
//   gk_b = [u * s + t * e + s^i * g, for 2 <= i <= l],
// where u consists of independent and uniformly random polynomials, e is a
// vector of error polynomials, and g is the gadget vector of dimension k.
// In other words, the relinearization key consists of encryptions of s^2,..,
// s^l under the secret key s. This relinearization key can "key switch" a
// ciphertext of degree up to l to a canonical ciphertext (of degree 1).
template <typename ModularInt>
class RnsRelinKey {
  using Integer = typename ModularInt::Int;
  using ModularIntParams = typename ModularInt::Params;

 public:
  // Generates a relinearization key for working with BGV ciphertexts whose
  // degree is at most `degree`.
  static absl::StatusOr<RnsRelinKey> CreateForBgv(
      const RnsRlweSecretKey<ModularInt>& secret_key, int degree, int variance,
      const RnsGadget<ModularInt>* gadget, Integer plaintext_modulus,
      PrngType prng_type) {
    return Create(secret_key, degree, variance, gadget, prng_type,
                  /*error_scalar=*/plaintext_modulus);
  }

  // Generates a relinearization key for working with BFV ciphertexts whose
  // degree is at most `degree`.
  static absl::StatusOr<RnsRelinKey> CreateForBfv(
      const RnsRlweSecretKey<ModularInt>& secret_key, int degree, int variance,
      const RnsGadget<ModularInt>* gadget, PrngType prng_type) {
    return Create(secret_key, degree, variance, gadget, prng_type);
  }

  // Samples the random polynomial pads "gk_as" used in a relinerization key.
  static absl::StatusOr<std::vector<RnsPolynomial<ModularInt>>> SampleRandomPad(
      int dimension, int degree, int log_n,
      absl::Span<const PrimeModulus<ModularInt>* const> moduli,
      absl::string_view prng_seed, PrngType prng_type);

  // Applies the relinearization key to a BGV ciphertext.
  absl::StatusOr<RnsBgvCiphertext<ModularInt>> ApplyTo(
      const RnsBgvCiphertext<ModularInt>& ciphertext) const {
    // First, key-switch the ciphertext as a generic RLWE ciphertext.
    RLWE_ASSIGN_OR_RETURN(std::vector<RnsPolynomial<ModularInt>> components_new,
                          ApplyToRlweCiphertext(ciphertext));
    // Compute the error bound in the new ciphertext.
    double error = ErrorAfterApplication(ciphertext);
    return RnsBgvCiphertext<ModularInt>(std::move(components_new), moduli_,
                                        /*power=*/1, error,
                                        ciphertext.ErrorParams());
  }

  // Applies the relinearization key to a BFV ciphertext
  absl::StatusOr<RnsBfvCiphertext<ModularInt>> ApplyTo(
      const RnsBfvCiphertext<ModularInt>& ciphertext) const {
    // First, key-switch the ciphertext as a generic RLWE ciphertext.
    RLWE_ASSIGN_OR_RETURN(std::vector<RnsPolynomial<ModularInt>> components_new,
                          ApplyToRlweCiphertext(ciphertext));
    // Compute the error bound in the new ciphertext.
    double error = ErrorAfterApplication(ciphertext);
    return RnsBfvCiphertext<ModularInt>(
        std::move(components_new), moduli_,
        /*power=*/1, error, ciphertext.ErrorParams(), ciphertext.Context());
  }

  // Accessors to the key components.
  const std::vector<RnsPolynomial<ModularInt>>& GetKeyA() const {
    return key_as_;
  }
  const std::vector<RnsPolynomial<ModularInt>>& GetKeyB() const {
    return key_bs_;
  }

  const RnsGadget<ModularInt>* Gadget() { return gadget_; }

  int Dimension() const { return gadget_->Dimension(); }

  int Degree() const { return degree_; }

 private:
  // Factory function that samples a relinearization key for different RLWE
  // schemes. In particular, `error_scalar` should be set to the plaintext
  // modulus for Relinearization key in BGV, and it should be set to 1
  // otherwise.
  static absl::StatusOr<RnsRelinKey> Create(
      const RnsRlweSecretKey<ModularInt>& secret_key, int degree, int variance,
      const RnsGadget<ModularInt>* gadget, PrngType prng_type,
      Integer error_scalar = 1);

  // Factory function that samples a relinearization key derived from
  // `secret_key` for the given degree, using the given random polynomials
  // `pads` as the "a" part of the relinearization key.
  static absl::StatusOr<RnsRelinKey> CreateWithRandomPad(
      std::vector<RnsPolynomial<ModularInt>> pads,
      const RnsRlweSecretKey<ModularInt>& secret_key, int degree, int variance,
      const RnsGadget<ModularInt>* gadget, absl::string_view prng_pad_seed,
      PrngType prng_type, Integer error_scalar = 1);

  explicit RnsRelinKey(std::vector<RnsPolynomial<ModularInt>> key_as,
                       std::vector<RnsPolynomial<ModularInt>> key_bs,
                       const RnsGadget<ModularInt>* gadget, int degree,
                       std::vector<const PrimeModulus<ModularInt>*> moduli,
                       absl::string_view prng_seed, PrngType prng_type)
      : key_as_(std::move(key_as)),
        key_bs_(std::move(key_bs)),
        gadget_(gadget),
        degree_(degree),
        moduli_(std::move(moduli)),
        prng_seed_(std::string(prng_seed)),
        prng_type_(prng_type) {}

  // Applies the relinearization key to a generic RLWE ciphertext, and returns
  // the component polynomials of the resulting ciphertext.
  absl::StatusOr<std::vector<RnsPolynomial<ModularInt>>> ApplyToRlweCiphertext(
      const RnsRlweCiphertext<ModularInt>& ciphertext) const;

  double ErrorAfterApplication(
      const RnsRlweCiphertext<ModularInt>& ciphertext) const {
    return ciphertext.Error() +
           ciphertext.ErrorParams()->BoundOnGadgetBasedKeySwitching(
               /*num_components=*/ciphertext.Degree() - 1,
               gadget_->LogGadgetBase(), gadget_->Dimension());
  }

  // The twp columns of the key matrix.
  std::vector<RnsPolynomial<ModularInt>> key_as_;
  std::vector<RnsPolynomial<ModularInt>> key_bs_;

  // The gadget for creating this Relinearization key; does not take ownership.
  const RnsGadget<ModularInt>* gadget_;

  // The max degree of the source secret key encrypted under this relin key.
  int degree_;

  // The RNS moduli used to construct this Relinearization key.
  std::vector<const PrimeModulus<ModularInt>*> moduli_;

  // PRNG seed for sampling the random polynomials in `key_as_`.
  std::string prng_seed_;

  // PRNG type for sampling the random polynomials in `key_as_`.
  PrngType prng_type_;
};

}  // namespace rlwe

#endif  // RLWE_RNS_RNS_RELINEARIZATION_KEY_H_
