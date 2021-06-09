/*
 * Copyright 2021 Google LLC
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef PRIVACY_PRIVATE_RETRIEVAL_EXPAND_H_
#define PRIVACY_PRIVATE_RETRIEVAL_EXPAND_H_

#include "absl/memory/memory.h"
#include "galois_key.h"
#include "integral_types.h"
#include "ntt_parameters.h"
#include "polynomial.h"
#include "status_macros.h"
#include "statusor.h"
#include "symmetric_encryption.h"

namespace rlwe {

namespace internal {

// Homomorphic Absorb of a monomial stored in a Polynomial. This function is
// exactly like the homomorphic absorb in the SymmetricRlweCiphertext
// class, except that since the plaintext is a monomial, the error estimate in
// the ciphertext is _not_ modified.
template <typename ModularInt>
absl::Status AbsorbNttMonomialInPlace(
    SymmetricRlweCiphertext<ModularInt>* ciphertext,
    const Polynomial<ModularInt>& monomial_ntt,
    const typename ModularInt::Params* modulus_params) {
  double error_estimate = ciphertext->Error();
  RLWE_RETURN_IF_ERROR(ciphertext->AbsorbInPlace(monomial_ntt));
  ciphertext->SetError(error_estimate);  // Reset the error estimate.
  return absl::OkStatus();
}

}  // namespace internal

// Creates an ObliviousExpander object used to perform the ObliviousExpand
// protocol. This objects holds information about the available GaloisKeys and
// uses this to specialize the way ciphertexts are relinearized after
// Substitutions are performed in Expand.
//
// The types of implementations available are:
//   - DefaultObliviousExpander: no GaloisKeys are supplied, and ObliviousExpand
//   is a no-op.
//   - GaloisKeysObliviousExpander: GaloisKeys corresponding to each level of
//   ObliviousExpand are supplied.
//   - GaloisGeneratorObliviousExpander: a single GaloisKey with a substitution
//   power
//     corresponding to a generator for the Galois group is supplied.
template <typename ModularInt>
class ObliviousExpander {
 public:
  ObliviousExpander(const int log_t,
                    const typename ModularInt::Params* modulus_params,
                    const NttParameters<ModularInt>* ntt_params)
      : modulus_params_(modulus_params), ntt_params_(ntt_params) {}

  virtual ~ObliviousExpander() = default;

  // Oblivious Expansion of n ciphertexts. Given a vector of ciphertexts where
  // the i'th ciphertext encrypts x ^ j and the rest encrypt 0,  and a set of k
  // GaloisKeys, this produces an output_length vector of m ciphertexts such
  // that the ((i * (2 ^ k)) + j)-th one is an encryption of one, and the rest
  // are encryptions of zero. This is the inverse function of pir_client's
  // BuildRequest.
  // Note that a normalization step is performed in pir_client's BuildRequest.
  // Parameters:
  //   ciphertexts   - The i'th ciphertext encrypts x ^ j and the rest encrypt
  //                   0, normalized by 2^-log_compression_factor.
  //                   j must be less than 2 ^ k.
  //   output_length - The output length m must be less than (n * (2 ^ k)).
  // Outputs:
  //    A vector of m ciphertexts such that the ((i * (2 ^ k)) + j)-th one
  //    encrypts 1, and the rest encrypt 0.
  // Failure Cases:
  //    - the output_length is out of range.
  virtual StatusOr<std::vector<SymmetricRlweCiphertext<ModularInt>>>
  ObliviousExpand(std::vector<SymmetricRlweCiphertext<ModularInt>> ciphertexts,
                  int levels_of_expand, int output_length) {
    // If the number of ciphertexts is already output_length, return
    // ciphertexts.
    if (ciphertexts.size() == output_length) {
      return ciphertexts;
    }

    // Check that the output length is within range.
    int max_length = ciphertexts.size() * (1 << levels_of_expand);
    if (output_length > max_length) {
      return absl::InvalidArgumentError(
          absl::StrCat("Output length must be at most ", max_length));
    }

    std::vector<SymmetricRlweCiphertext<ModularInt>> result;
    result.reserve(max_length);

    for (const auto& ciphertext : ciphertexts) {
      RLWE_ASSIGN_OR_RETURN(
          auto expanded,
          ExpandCiphertext(std::move(ciphertext), levels_of_expand));
      result.insert(result.end(), std::make_move_iterator(expanded.begin()),
                    std::make_move_iterator(expanded.end()));
    }

    return std::vector<SymmetricRlweCiphertext<ModularInt>>(
        std::make_move_iterator(result.begin()),
        std::make_move_iterator(result.begin() + output_length));
  }

  // Computes the normalizer user in ObliviousExpand with
  // log_compression_factor = k.
  // This is equal to 2^(-k) mod ((2 ^ log_t) + 1).
  static Uint64 ComputeNormalizer(int k, int log_t) {
    // Let k = p*log_t + r. Then the inverse of (2 ^ k) mod ((2 ^ log_t) + 1) is
    // (2 ^ (log_t - r)) * (2 ^ ((-p - 1) * log_t)).
    int log_t_minus_r = log_t - (k % log_t);
    // Note that (2 ^ ((-p - 1) * log_t)) = ((-1) ^ (p + 1)) mod ((2 ^ log_t)
    // + 1).
    if (int p = k / log_t; p % 2 == 0) {
      // Return - 2 ^ (log_t - r)
      return ((static_cast<Uint64>(1) << log_t) + 1) -
             (static_cast<Uint64>(1) << log_t_minus_r);
    } else {
      // Otherwise, Return 2 ^ (log_t - r)
      return (static_cast<Uint64>(1) << log_t_minus_r);
    }
  }

 protected:
  // Modulus params.
  const typename ModularInt::Params* modulus_params_;

  // Ntt params.
  const NttParameters<ModularInt>* ntt_params_;

  // Computes a substitution and relinearization with the available galois keys.
  // The substitution power must be of the form (N/(2 ^ j) + 1).
  virtual StatusOr<SymmetricRlweCiphertext<ModularInt>> SubstitutionWithKeys(
      const SymmetricRlweCiphertext<ModularInt>& ciphertext,
      int substitution_power) = 0;

  // Expands a single ciphertext.
  StatusOr<std::vector<SymmetricRlweCiphertext<ModularInt>>> ExpandCiphertext(
      SymmetricRlweCiphertext<ModularInt> ciphertext, int levels_of_expand) {
    int output_size = 1 << levels_of_expand;
    RLWE_ASSIGN_OR_RETURN(auto comp0, ciphertext.Component(0));
    int num_coeffs = comp0.Len();

    std::vector<SymmetricRlweCiphertext<ModularInt>> result(
        output_size, SymmetricRlweCiphertext<ModularInt>(
                         modulus_params_, ciphertext.ErrorParams()));
    result[0] = std::move(ciphertext);

    ModularInt zero = ModularInt::ImportZero(modulus_params_);
    ModularInt minus_one = ModularInt::ImportOne(modulus_params_);
    minus_one.NegateInPlace(modulus_params_);

    // Each iteration of the loop doubles the number of ciphertexts.
    for (int j = 0; j < levels_of_expand; j++) {
      // Computes the monomial x^{-(2 ^ j)} to_absorb as -x^{N - (2 ^ j)} since
      // our polynomials are modulo (x ^ N + 1).
      std::vector<ModularInt> monomial(num_coeffs, zero);
      monomial[num_coeffs - (1 << j)] = minus_one;
      auto to_absorb = Polynomial<ModularInt>::ConvertToNtt(
          monomial, ntt_params_, modulus_params_);
      int substitution_power = (num_coeffs >> j) + 1;

      for (int k = 0; k < (1 << j); k++) {
        // Let c0 = result[k], and c1 = c0 * to_absorb.
        // This inner loop computes result[k] and result[k + 2 ^ j] as:
        //    result[k] = c0 + Sub(c0, substitution_power)
        //    result[k + 2 ^ j] = c1 + Sub(c1, substitution_power)
        //                      = (c0 - Sub(c0, substitution_power)) * to_absorb
        // since Sub(c1, substitution_power)
        //       = Sub(c0, substitution_power) * Sub(to_absorb,
        //                     substitution_power) = substituted * (-to_absorb).

        // Computes Sub(c0, substitution_power).
        RLWE_ASSIGN_OR_RETURN(
            auto sub, SubstitutionWithKeys(result[k], substitution_power));

        // result[k + 2 ^ j] = (c0 - Sub(c0, substitution_power)) * to_absorb
        RLWE_ASSIGN_OR_RETURN(result[k + (1 << j)], result[k] - sub);
        RLWE_RETURN_IF_ERROR(internal::AbsorbNttMonomialInPlace(
            &result[k + (1 << j)], to_absorb, modulus_params_));
        // result[k] = c0 + Sub(c0, substitution_power)
        RLWE_RETURN_IF_ERROR(result[k].AddInPlace(sub));
      }
    }

    return result;
  }
};

// ObliviousExpand using each a list of Galois keys corresponding to the
// substitutions needed at each level of expand.
template <typename ModularInt>
class GaloisKeysObliviousExpander : public ObliviousExpander<ModularInt> {
 public:
  // The galois keys must be a vector such that the i'th GaloisKey relinearizes
  // the substitution power (N / (2 ^ i)) + 1.
  static StatusOr<std::unique_ptr<GaloisKeysObliviousExpander>> Create(
      std::vector<GaloisKey<ModularInt>> galois_keys, const int log_t,
      const typename ModularInt::Params* modulus_params,
      const NttParameters<ModularInt>* ntt_params) {
    // Expect that the galois keys have the correct substitution powers.
    int num_coeffs = ntt_params->number_coeffs;
    std::map<int, int> key_index_for_substitution;

    for (int i = 0; i < galois_keys.size(); i++) {
      int substitution_power = (num_coeffs >> i) + 1;
      if (galois_keys[i].SubstitutionPower() != substitution_power) {
        return absl::InvalidArgumentError(
            "Substitution power of the i'th GaloisKey must be (N / (2 ^ i)) + "
            "1");
      }
      key_index_for_substitution[substitution_power] = i;
    }

    return absl::make_unique<GaloisKeysObliviousExpander>(
        std::move(galois_keys), std::move(key_index_for_substitution), log_t,
        modulus_params, ntt_params);
  }

  GaloisKeysObliviousExpander(std::vector<GaloisKey<ModularInt>> galois_keys,
                              std::map<int, int> key_index_for_substitution,
                              const int log_t,
                              const typename ModularInt::Params* modulus_params,
                              const NttParameters<ModularInt>* ntt_params)
      : ObliviousExpander<ModularInt>(log_t, modulus_params, ntt_params),
        key_index_for_substitution_(std::move(key_index_for_substitution)),
        galois_keys_(std::move(galois_keys)) {}

 protected:
  // If a galois key with a matching substitution_power is found, the
  // substitution is applied directly and the ciphertext is relinearized with
  // the corresponding galois key.
  StatusOr<SymmetricRlweCiphertext<ModularInt>> SubstitutionWithKeys(
      const SymmetricRlweCiphertext<ModularInt>& ciphertext,
      int substitution_power) override {
    // Use GaloisKey matching substitution power.
    RLWE_ASSIGN_OR_RETURN(
        auto sub, ciphertext.Substitute(substitution_power, this->ntt_params_));
    return galois_keys_[key_index_for_substitution_.find(substitution_power)
                            ->second]
        .ApplyTo(sub);
  }

 private:
  // Mapping between substitution power corresponding index in galois keys.
  const std::map<int, int> key_index_for_substitution_;

  // The galois keys.
  const std::vector<GaloisKey<ModularInt>> galois_keys_;
};

// ObliviousExpander using substitutions by a generator for the galois group.
// This is not compatible with using a levels_of_expand equal to the log of the
// number of coefficients, and will raise an InvalidArgumentError during
// expansion in that case.
template <typename ModularInt>
class GaloisGeneratorObliviousExpander : public ObliviousExpander<ModularInt> {
 public:
  static StatusOr<std::unique_ptr<GaloisGeneratorObliviousExpander>> Create(
      GaloisKey<ModularInt> galois_generator, const int log_t,
      const typename ModularInt::Params* modulus_params,
      const NttParameters<ModularInt>* ntt_params) {
    // Verify that the galois generator matches.
    if (galois_generator.SubstitutionPower() != kSubstitutionGenerator) {
      return absl::InvalidArgumentError(
          absl::StrCat("Substitution power of the generator GaloisKey must be ",
                       kSubstitutionGenerator));
    }

    return absl::make_unique<GaloisGeneratorObliviousExpander>(
        std::move(galois_generator), log_t, modulus_params, ntt_params);
  }

  GaloisGeneratorObliviousExpander(
      GaloisKey<ModularInt> galois_generator, const int log_t,
      const typename ModularInt::Params* modulus_params,
      const NttParameters<ModularInt>* ntt_params)
      : ObliviousExpander<ModularInt>(log_t, modulus_params, ntt_params),
        generator_(std::move(galois_generator)),
        number_coefficients_(ntt_params->number_coeffs) {}

  static const int GetSubstitutionGenerator() { return kSubstitutionGenerator; }

 protected:
  // The substitution generator is applied repeatedly to obtain the
  // substitution_power and the galois key generator is used to relinearize
  // after each step.
  StatusOr<SymmetricRlweCiphertext<ModularInt>> SubstitutionWithKeys(
      const SymmetricRlweCiphertext<ModularInt>& ciphertext,
      int substitution_power) override {
    // Use a generator by applying the substitution num_substitution_k times and
    // relinearizing with GaloisKey[k].
    auto result = ciphertext;
    RLWE_ASSIGN_OR_RETURN(auto comp0, ciphertext.Component(0));
    RLWE_ASSIGN_OR_RETURN(
        int num_substitutions,
        ComputeNumberSubstitutionsWithGenerator(substitution_power));
    for (int i = 0; i < num_substitutions; i++) {
      RLWE_ASSIGN_OR_RETURN(auto sub, result.Substitute(kSubstitutionGenerator,
                                                        this->ntt_params_));
      RLWE_ASSIGN_OR_RETURN(result, generator_.ApplyTo(sub));
    }
    return result;
  }

 private:
  // For any polynomial modulus f(x) = x ^ (N/(2 ^ j) + 1) for 1 <= 2 ^ j <=
  // N/4, the substitution x->x^5 can be applied repeatedly to get all possible
  // substitution powers.
  static const int kSubstitutionGenerator = 5;

  // On input substitution_power, this function outputs b such that the
  // substitution x->x^substitution_power can be computed by applying
  // (x->x^kSubstitutionGenerator) b times successively.
  rlwe::StatusOr<int> ComputeNumberSubstitutionsWithGenerator(
      int substitution_power) {
    if (substitution_power < 5 ||
        substitution_power >= 2 * number_coefficients_ ||
        internal::CountOnes64(static_cast<uint64_t>(substitution_power - 1)) !=
            1) {
      return absl::InvalidArgumentError(
          absl::StrCat("The substitution power, ", substitution_power,
                       " must be an integer of the form "
                       "(2N/2^j) + 1 >= 5, where N=",
                       number_coefficients_,
                       ". This may be due to using levels_of_expand = ",
                       ceil(log2(number_coefficients_)),
                       " when generating the ObliviousExpander."));
    }

    int power = 1;  // This invariant represents kSubstitutionGenerator^b % 2N.
    for (int b = 1; b < 2 * number_coefficients_; b++) {
      power = (power * kSubstitutionGenerator) % (2 * number_coefficients_);
      if (power == substitution_power) {
        return b;
      }
    }

    // This error should not be raised. If it is raised, it means that a cycle
    // has been found, and therefore x->x^5 cannot be used repeatedly.
    return absl::InvalidArgumentError(
        absl::StrCat("The substitution power, ", substitution_power,
                     " cannot be obtained in less than ",
                     2 * number_coefficients_, " iterations."));
  }

  // The galois key corresponding to a generator for the galois group.
  const GaloisKey<ModularInt> generator_;
  int number_coefficients_;
};

// Does not expand or manipulate ciphertexts.
template <typename ModularInt>
class DefaultObliviousExpander : public ObliviousExpander<ModularInt> {
 public:
  static StatusOr<std::unique_ptr<DefaultObliviousExpander>> Create(
      const int log_t, const typename ModularInt::Params* modulus_params,
      const NttParameters<ModularInt>* ntt_params) {
    return absl::make_unique<DefaultObliviousExpander>(log_t, modulus_params,
                                                       ntt_params);
  }

  DefaultObliviousExpander(const int log_t,
                           const typename ModularInt::Params* modulus_params,
                           const NttParameters<ModularInt>* ntt_params)
      : ObliviousExpander<ModularInt>(log_t, modulus_params, ntt_params) {}

  StatusOr<std::vector<SymmetricRlweCiphertext<ModularInt>>> ObliviousExpand(
      std::vector<SymmetricRlweCiphertext<ModularInt>> ciphertexts,
      int levels_of_expand, int output_length) override {
    // No expansion, so ciphertexts must match the output length.
    if (ciphertexts.size() != output_length) {
      return absl::InvalidArgumentError(
          absl::StrCat("Output length must be ", ciphertexts.size()));
    }

    return std::move(ciphertexts);
  }

  StatusOr<SymmetricRlweCiphertext<ModularInt>> SubstitutionWithKeys(
      const SymmetricRlweCiphertext<ModularInt>& ciphertext,
      int substitution_power) override {
    return ciphertext;
  }
};

// Returns a compressed vector with length (total_size / compression_factor)
// such that, for every index in the vector indices, the entry at
// (index / compression_factor) is a Polynomial containing the monomial
// x^{index % compression_factor}. This is the compression of
// a vector of length total_size with 1's at the positions specified by indices,
// and the rest 0.
// Returns an error if any index is greater than total_size, or
// compression_factor is greater than the number of coefficients of the
// Polynomial.
template <typename ModularInt>
StatusOr<std::vector<Polynomial<ModularInt>>> MakeCompressedVector(
    int total_size, const std::vector<int>& indices, int log_compression_factor,
    const typename ModularInt::Params* modulus_params,
    const NttParameters<ModularInt>* ntt_params) {
  int compression_factor = 1 << log_compression_factor;
  // Check error conditions.
  if (compression_factor > ntt_params->number_coeffs) {
    return absl::InvalidArgumentError(absl::StrCat(
        "Compression factor must be less than ", ntt_params->number_coeffs));
  }

  // Initialize the coefficients of the response to 0.
  std::vector<ModularInt> zero_vec(ntt_params->number_coeffs,
                                   ModularInt::ImportZero(modulus_params));
  int size_compression =
      ceil(total_size / static_cast<double>(compression_factor));
  std::vector<std::vector<ModularInt>> vector_coefficients(size_compression,
                                                           zero_vec);

  for (const int& index : indices) {
    if (index < 0 || index >= total_size) {
      return absl::InvalidArgumentError("Index out of range for total size.");
    }

    // Set the correct coefficient to 1.
    int polynomial_index = index >> log_compression_factor;
    vector_coefficients[polynomial_index][index % compression_factor] =
        ModularInt::ImportOne(modulus_params);
  }

  std::vector<Polynomial<ModularInt>> result;
  result.reserve(size_compression);

  for (auto& v : vector_coefficients) {
    result.push_back(Polynomial<ModularInt>::ConvertToNtt(
        std::move(v), ntt_params, modulus_params));
  }

  return result;
}

}  // namespace rlwe

#endif  // PRIVACY_PRIVATE_RETRIEVAL_EXPAND_H_
