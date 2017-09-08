/*
 * Copyright 2017 Google Inc.
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

#ifndef RLWE_SYMMETRIC_ENCRYPTION_H_
#define RLWE_SYMMETRIC_ENCRYPTION_H_

#include <cstdint>
#include <vector>
#include "ntt_polynomial.h"
#include "serialization.pb.h"
#include "utils.h"

namespace rlwe {

// This file implements the somewhat homomorphic symmetric-key encryption scheme
// from "Fully Homomorphic Encryption from Ring-LWE and Security for Key
// Dependent Messages" by Zvika Brakerski and Vinod Vaikuntanathan. This
// encryption scheme uses Ring Learning with Errors (RLWE).
//
// The encryption scheme in this file is not fully homomorphic. It does not
// implement any sort of bootstrapping.

// Forward references.
template <typename ModularInt>
class SymmetricRlweCiphertext;

template <typename ModularInt>
class SymmetricRlweKey;

template <typename ModularInt>
SymmetricRlweCiphertext<ModularInt> Encrypt(
    const SymmetricRlweKey<ModularInt>& key,
    const NttPolynomial<ModularInt>& plaintext);

template <typename ModularInt>
std::vector<typename ModularInt::Int> Decrypt(
    const SymmetricRlweKey<ModularInt>& key,
    const SymmetricRlweCiphertext<ModularInt>& ciphertext);

// Represents a ciphertext encrypted using a symmetric-key version of the ring
// learning-with-errors (RLWE) encryption scheme. See the comments that follow
// throughout this file for full details on the particular encryption scheme.
//
// This implementation supports the following homomorphic operations:
//  - Homomorphic addition.
//  - Scalar multiplication by a polynomial (absorbtion)
//  - Homomorphic multiplication.
//
// This implementation is only "somewhat homomorphic," not fully homomorphic.
// There is no bootstrapping, so a limited number of homomorphic operations can
// be performed before so much error accumulates that decryption is impossible.
//
// Each ciphertext comprises a vector of polynomials <c0, ..., cN>. Initially,
// a ciphertext comprises a pair <c0, c1>. Homomorphic multiplications cause
// the vector to grow longer.
template <typename ModularInt>
class SymmetricRlweCiphertext {
 public:
  // Default and copy constructors.
  explicit SymmetricRlweCiphertext(const typename ModularInt::Params* params)
      : modulus_params_(params) {}
  SymmetricRlweCiphertext(const SymmetricRlweCiphertext& that) = default;

  // Create a ciphertext by supplying the vector of components.
  explicit SymmetricRlweCiphertext(
      const std::vector<NttPolynomial<ModularInt>>& c,
      const typename ModularInt::Params* params)
      : c_(c), modulus_params_(params) {}

  // Homomorphic addition: add the polynomials representing the ciphertexts
  // component-wise. The example below demonstrates why this procedure works
  // properly in the two-component case. The quantities a, s, m, t, and e are
  // introduced during encryption and are explained in the SymmetricRlweKey
  // class.
  //
  //   (a1 * s + m1 + t * e1, -a1)
  // + (a2 * s + m2 + t * e2, -a2)
  // ------------------------------
  //   ((a1 + a2) * s + (m1 + m2) + t * (e1 + e2), -(a1 + a2))
  //
  // Substitute (a1 + a2) = a3, (e1 + e2) = e3:
  //
  //   (a3 * s + (m1 + m2) + t * e3, -a3)
  //
  // This result is a valid ciphertext where the value of a has changed, the
  // error has increased, and the encoded plaintext contains the sum of the
  // plaintexts that were encoded in the original two ciphertexts.
  SymmetricRlweCiphertext operator+(const SymmetricRlweCiphertext& that) const {
    const SymmetricRlweCiphertext* longer = this;
    const SymmetricRlweCiphertext* shorter = &that;
    if (c_.size() < that.c_.size()) {
      std::swap(longer, shorter);
    }

    std::vector<NttPolynomial<ModularInt>> result(longer->c_);

    for (int i = 0; i < shorter->c_.size(); i++) {
      result[i] = result[i] + shorter->c_[i];
    }

    return SymmetricRlweCiphertext(result, modulus_params_);
  }

  // Homomorphic absorbtion. Multiplies the current ciphertext {m1}_s (plaintext
  // m1 encrypted  with symmetric key s) by a plaintext m2, resulting in a
  // ciphertext {m1 * m2}_s that stores m1 * m2 encrypted with symmetric key s.
  //
  // DO NOT CONFUSE THIS OPERATION WITH HOMOMORPHIC MULTIPLICATION.
  //
  // To perform this operation, multiply the each component of the
  // ciphertext by the plaintext polynomial. The example below demonstrates why
  // this procedure works properly in the two-component case. The quantities a,
  // s, m, t, and e are introduced during encryption and are explained in the
  // Encrypt() function later in this file.
  //
  //    (a1 * s + m1 + t * e1, -a1) * p
  //  = (a1 * s * p + m1 * p + t * e1 * p)
  //
  // Substitute (a1 * p) = a2 and (e1 * p) = e2:
  //
  //    (a2 * s + m1 * p + t * e2)
  //
  // This result is a valid ciphertext where the value of a has changed, the
  // error has increased, and the encoded plaintext contains the product of
  // m1 and p.
  //
  // A few more details about the multiplication that takes place:
  //
  // The value stored in the resulting ciphertext is (m1 * m2) (mod 2^N + 1)
  // (mod t), where N is the number of coefficients in s (or m1 or m2, since
  // the all have the same number of coefficients). In other words, the
  // result is the remainder of (m1 * m2) mod the polynomial (2^N + 1) with
  // each of the coefficients the ntaken mod t. Any coefficient between 0 and
  // modulus / 2 is treated as a positive number for the purposes of the final
  // (mod t); any coefficient between modulus/2 and modulus is treated as
  // a negative number for the purposes of the final (mod t).
  SymmetricRlweCiphertext operator*(
      const NttPolynomial<ModularInt>& that) const {
    std::vector<NttPolynomial<ModularInt>> result(c_.size());

    for (int i = 0; i < c_.size(); i++) {
      result[i] = c_[i] * that;
    }

    return SymmetricRlweCiphertext(result, modulus_params_);
  }

  // Homomorphic multiply. Given two ciphertexts {m1}_s, {m2}_s containing
  // messages m1 and m2 encrypted with the same secret key s, return the
  // ciphertext {m1 * m2}_s containing the product of the messages.
  //
  // To perform this operation, treat the two ciphertext vectors as polynomials
  // and perform a polynomial multiplication:
  //
  //   <c0, c1> * <c0', c1'> = <c0 * c0, c0 * c1 + c1 * c0, c1 * c1>
  //
  // If the two ciphertext vectors are of length m and n, the resulting
  // ciphertext is of length m + n - 1.
  //
  // The details of the multiplication that takes place between m1 and m2 are
  // the same as in the homomorphic absorb operation above (the other overload
  // of the * operator).
  SymmetricRlweCiphertext operator*(const SymmetricRlweCiphertext& that) {
    std::vector<NttPolynomial<ModularInt>> result(
        c_.size() + that.c_.size() - 1,
        NttPolynomial<ModularInt>(c_[0].Len(), modulus_params_));
    for (int i = 0; i < c_.size(); i++) {
      for (int j = 0; j < that.c_.size(); j++) {
        result[i + j] = result[i + j] + c_[i] * that.c_[j];
      }
    }

    return SymmetricRlweCiphertext(result, modulus_params_);
  }

  // Convert this ciphertext from (mod p) to (mod q).
  // Assumes that ModularInt::Int and ModularIntQ::Int are the same type.
  //
  // The current modulus (mod t) must be equal to modulus q (mod t).
  // This will always be true. For NTT to work properly, any modulus must be
  // of the form 2N + 1, where N is a power of 2. Likewise, the implementation
  // requires that t is a power of 2. This means that, for any modulus q and
  // modulus t allowed by the RLWE implementation, q % t == 1.
  template <typename ModularIntQ>
  SymmetricRlweCiphertext<ModularIntQ> SwitchModulus(
      const NttParameters<ModularInt>& ntt_params_p,
      const typename ModularIntQ::Params* modulus_params_q,
      const NttParameters<ModularIntQ>& ntt_params_q,
      typename ModularInt::Int t) {
    typename ModularInt::Int p = modulus_params_->modulus;
    typename ModularInt::Int q = modulus_params_q->modulus;

    // Configuration error.
    if (p % t != q % t) {
      std::cerr << "p % t != q % t\n";
      abort();
    }

    SymmetricRlweCiphertext<ModularIntQ> output(modulus_params_q);

    for (const NttPolynomial<ModularInt>& c : c_) {
      // Extract each component of the ciphertext from NTT form.
      std::vector<ModularInt> coeffs_p = c.InverseNtt(ntt_params_p);
      std::vector<ModularIntQ> coeffs_q;

      // Convert each coefficient of the polynomial from (mod p) to (mod q)
      for (const ModularInt& coeff_p : coeffs_p) {
        typename ModularInt::Int int_p = coeff_p.ExportInt();

        // Scale the integer.
        typename ModularInt::Int int_q = int_p * q / p;

        // Ensure that int_p = int_q mod t by changing int_q as little as
        // possible.
        typename ModularInt::Int int_p_mod_t = int_p % t;
        typename ModularInt::Int int_q_mod_t = int_q % t;
        typename ModularInt::Int adjustment_up, adjustment_down;

        // Determine whether to adjust int_q up or down to make sure int_q =
        // int_p (mod t).
        adjustment_up = int_p_mod_t - int_q_mod_t;
        adjustment_down = t + int_q_mod_t - int_p_mod_t;
        if (int_p_mod_t < int_q_mod_t) {
          adjustment_up += t;
          adjustment_down -= t;
        }

        if (adjustment_up > adjustment_down) {
          // Adjust up.
          coeffs_q.push_back(
              ModularIntQ::ImportInt(modulus_params_q, int_q) +
              ModularIntQ::ImportInt(modulus_params_q, adjustment_up));
        } else {
          // Adjust down.
          coeffs_q.push_back(
              ModularIntQ::ImportInt(modulus_params_q, int_q) +
              ModularIntQ::ImportInt(modulus_params_q, q - adjustment_down));
        }
      }

      // Convert back to NTT.
      output.c_.push_back(
          NttPolynomial<ModularIntQ>::ConvertToNtt(coeffs_q, ntt_params_q));
    }

    return output;
  }

  SerializedSymmetricRlweCiphertext Serialize() const {
    SerializedSymmetricRlweCiphertext output;

    for (const NttPolynomial<ModularInt>& c : c_) {
      *output.add_c() = c.Serialize();
    }

    return output;
  }

  static SymmetricRlweCiphertext Deserialize(
      const typename ModularInt::Params* modulus_params,
      const SerializedSymmetricRlweCiphertext& serialized) {
    SymmetricRlweCiphertext output(modulus_params);

    for (int i = 0; i < serialized.c_size(); i++) {
      output.c_.push_back(NttPolynomial<ModularInt>::Deserialize(
          modulus_params, serialized.c(i)));
    }

    return output;
  }

 private:
  // The ciphertext.
  std::vector<NttPolynomial<ModularInt>> c_;

  // ModularInt parameters.
  const typename ModularInt::Params* modulus_params_;

  friend SymmetricRlweCiphertext<ModularInt> Encrypt<ModularInt>(
      const SymmetricRlweKey<ModularInt>& key,
      const NttPolynomial<ModularInt>& plaintext);
  friend std::vector<typename ModularInt::Int> Decrypt<ModularInt>(
      const SymmetricRlweKey<ModularInt>& key,
      const SymmetricRlweCiphertext<ModularInt>& ciphertext);

  // Make this class a friend of any version of this class, no matter the
  // template.
  template <typename Q>
  friend class SymmetricRlweCiphertext;
};

template <typename ModularInt>
class SymmetricRlweKey {
 public:
  // Disallow copy constructor.
  SymmetricRlweKey(SymmetricRlweKey&& that) = default;

  // Static factory that samples a key from the error distribution. The
  // polynomial representing the key must have a number of coefficients that is
  // a power of two, which is enforced by the first argument.
  //
  // Does not take ownership of rand.
  static SymmetricRlweKey Sample(
      typename ModularInt::Int log_num_coeffs,
      typename ModularInt::Int variance, typename ModularInt::Int log_t,
      utils::SecurePrng* rand,
      const typename ModularInt::Params* modulus_params,
      const NttParameters<ModularInt>& ntt_params) {
    NttPolynomial<ModularInt> key = NttPolynomial<ModularInt>::ConvertToNtt(
        SampleFromErrorDistribution(modulus_params, 1 << log_num_coeffs,
                                    variance, rand),
        ntt_params);
    return SymmetricRlweKey(key, variance, log_t, modulus_params, ntt_params,
                            rand);
  }

  // Generate a copy of this key in modulus q.
  //
  // The current modulus (mod t) must be equal to modulus q (mod t). This
  // property is implicitly enforced by the design of the code as described
  // by the corresponding comment on SymmetricRlweKey::SwitchModulus. This
  // property is also dynamically enforced.
  template <typename ModularIntQ>
  SymmetricRlweKey<ModularIntQ> SwitchModulus(
      const typename ModularIntQ::Params* modulus_params_q,
      const NttParameters<ModularIntQ>& ntt_params_q) const {
    // Configuration failure.
    typename ModularInt::Int t = 1 << log_t_;
    if (modulus_params_->modulus % t != modulus_params_q->modulus % t) {
      std::cerr << "p % t != q % t\n";
      abort();
    }

    typename ModularIntQ::Int p_mod_q =
        modulus_params_->modulus % modulus_params_q->modulus;
    std::vector<ModularInt> coeffs_p = key_.InverseNtt(ntt_params_);
    std::vector<ModularIntQ> coeffs_q;

    // Convert each coefficient of the polynomial from (mod p) to (mod q)
    for (const ModularInt& coeff_p : coeffs_p) {
      // Ensure that negative numbers (mod p) are translated into negative
      // numbers (mod q).
      typename ModularInt::Int int_p = coeff_p.ExportInt();
      if (int_p > modulus_params_->modulus / 2) {
        int_p -= p_mod_q;
      }

      coeffs_q.push_back(ModularIntQ::ImportInt(modulus_params_q, int_p));
    }

    // Convert back to NTT.
    auto key_q =
        NttPolynomial<ModularIntQ>::ConvertToNtt(coeffs_q, ntt_params_q);

    return SymmetricRlweKey<ModularIntQ>(key_q, variance_, log_t_,
                                         modulus_params_q, ntt_params_q, rand_);
  }

  // Accessors.
  unsigned int Len() const { return key_.Len(); }
  NttParameters<ModularInt> NttParams() const { return ntt_params_; }
  const typename ModularInt::Params* ModulusParams() const {
    return modulus_params_;
  }
  unsigned int BitsPerCoeff() const { return log_t_; }

 private:
  // The contents of the key itself.
  NttPolynomial<ModularInt> key_;

  // A secure random number generator. Does not have ownership.
  utils::SecurePrng* rand_;

  // The variance of the binomial distribution from which the key and error are
  // drawn.
  typename ModularInt::Int variance_;

  // The maximum size of any one coefficient of the polynomial representing a
  // plaintext message.
  typename ModularInt::Int log_t_;
  ModularInt t_mod_;

  // NTT parameters.
  NttParameters<ModularInt> ntt_params_;

  // ModularInt parameters.
  const typename ModularInt::Params* modulus_params_;

  // A constructor. Does not take ownership of rand.
  SymmetricRlweKey(const NttPolynomial<ModularInt>& key,
                   typename ModularInt::Int variance,
                   typename ModularInt::Int log_t,
                   const typename ModularInt::Params* modulus_params,
                   const NttParameters<ModularInt>& ntt_params,
                   utils::SecurePrng* rand)
      : key_(key),
        rand_(rand),
        variance_(variance),
        log_t_(log_t),
        t_mod_(ModularInt::ImportInt(modulus_params, 1 << log_t)),
        ntt_params_(ntt_params),
        modulus_params_(modulus_params) {}

  // Samples a vector of coefficients from the centered binomial distribution
  // with the specified variance. The RLWE proofs rely on
  // sampling keys and error values from a discrete Gaussian distribution, but
  // the NewHope paper indicates that a centered binomial distribution is
  // indistinguishable and is far easier to sample.
  //
  // All values sampled are multiplied by scalar.
  static std::vector<ModularInt> SampleFromErrorDistribution(
      const typename ModularInt::Params* modulus_params,
      unsigned int num_coeffs, typename ModularInt::Int variance,
      utils::SecurePrng* rand) {
    typename ModularInt::Int k = variance << 1;
    std::vector<ModularInt> coeffs(num_coeffs,
                                   ModularInt::ImportInt(modulus_params, 0));

    for (int i = 0; i < num_coeffs; i++) {
      // Sample from the centered binomial distribution. To do so, we sample
      // k pairs of bits (a, b), where k = 2 * variance. The sample of the
      // binomial distribution is the sum of the differences between each pair
      // of bits.
      //
      // Concretely, when variance = 2, k = 4 and our sample is:
      //   (a0 - b0) + (a1 - b1) + (a2 - b2) + (a3 - b3)
      for (int j = 0; j < k; j++) {
        ModularInt b0 =
            ModularInt::ImportInt(modulus_params, rand->Rand8() & 0x01);
        ModularInt b1 =
            ModularInt::ImportInt(modulus_params, rand->Rand8() & 0x01);

        coeffs[i] = coeffs[i] + (b0 - b1);
      }
    }

    return coeffs;
  }

  friend SymmetricRlweCiphertext<ModularInt> Encrypt<ModularInt>(
      const SymmetricRlweKey<ModularInt>& key,
      const NttPolynomial<ModularInt>& plaintext);
  friend std::vector<typename ModularInt::Int> Decrypt<ModularInt>(
      const SymmetricRlweKey<ModularInt>& key,
      const SymmetricRlweCiphertext<ModularInt>& ciphertext);

  // Make this class a friend of any version of this class, no matter the
  // template.
  template <typename Q>
  friend class SymmetricRlweKey;
};

// Encrypts the plaintext using ring learning-with-errors (RLWE) encryption.
//
// The scheme works as follows:
//   KeyGen(n, modulus q, error distr):
//     Sample a degree (n-1) polynomial whose coefficients are drawn from the
//     error distribution (mod q). This is our secret key. Call it s.
//
//   Encrypt(secret key s, plaintext m, modulus q, modulus t, error distr):
//     1) Sample a degree (n-1) polynomial whose coefficients are drawn
//        uniformly from any integer (mod q). Call this polynomial a.
//     2) Sample a degree (n-1) polynomial whose coefficients are drawn from
//        the error distribution (mod q). Call this polynomial e.
//     3) Our secret key s and plaintext m are both degree (n-1) polynomials.
//        For decryption to work, each coefficient of m must be < t.
//        Compute (a * s + t * e + m) (mod x^n + 1). Call this polynomial b.
//     4) The ciphertext is the pair (b, -a). We refer to the pair of
//        polynomials representing a ciphertext as (c0, c1) =
//        (a * s + m + e * t, -a).
//
//    Decrypt(secret key s, ciphertext (b, -a), modulus t):
//      // Decryption when the ciphertext has two components.
//      Compute and return (b - as) (mod t). Doing out the algebra:
//          b - as (mod t)
//        = as + te + m - as (mod t)
//        = te + m (mod t)
//        = m
//      Quoting the paper, "the condition for correct decryption is that the
//      L_infinity norm of the polynomial [te + m] is smaller than q/2." In
//      other words, the largest of the values te + m (recall that e is
//      sampled from a distribution) cannot exceed q/2.
//
//   When the ciphertext has more than two components <c0, c1, ..., cN>,
//   it can be decrypted by taking the dot product with the vector
//   <s^0, s^1, ..., s^N> containing powers of the secret key:
//       te + m = <c0, 1, ..., cN> dot <s^0, s^1, ..., s^N>
//              = c0 * s^0 + c1 * s^1 + ... + cN * s^N
//
// Note that the Encrypt() function takes the original plaintext as
// an NttPolynomial<ModularInt>, while the corresponding Decrypt() method
// returns a std::vector<uint64_t>. The two values will be the same once
// the original plaintext is converted out of NTT and Montgomery form.
//   - The Encrypt() function takes an NTT polynomial so that, if the same
//     plaintext is to be encrypted repeatedly, the NTT conversion only needs
//     to be performed once by the caller.
//   - The Decrypt() function returns a vector of integers because the final
//     (mod t) step requires taking the polynomial (te + m) out of NTT and
//     Montgomery form.
// It would be straightforward to write a wrapper of Encrypt() that takes
// a vector of integers as input, thereby making the plaintext types of the
// Encrypt() and Decrypt() functions symmetric.
template <typename ModularInt>
SymmetricRlweCiphertext<ModularInt> Encrypt(
    const SymmetricRlweKey<ModularInt>& key,
    const NttPolynomial<ModularInt>& plaintext) {
  // Sample the error term from the error distribution.
  unsigned int num_coeffs = key.key_.Len();
  std::vector<ModularInt> e_coeffs = key.SampleFromErrorDistribution(
      key.modulus_params_, num_coeffs, key.variance_, key.rand_);

  // Sample a from the uniform distribution.
  std::vector<ModularInt> a_coeffs(
      num_coeffs, ModularInt::ImportInt(key.modulus_params_, 0));
  for (int i = 0; i < num_coeffs; i++) {
    a_coeffs[i] = ModularInt::ImportRandom(key.modulus_params_, key.rand_);
  }

  // Create c0.
  auto a = NttPolynomial<ModularInt>::ConvertToNtt(a_coeffs, key.ntt_params_);
  auto e = NttPolynomial<ModularInt>::ConvertToNtt(e_coeffs, key.ntt_params_);
  NttPolynomial<ModularInt> c0 = (a * key.key_) + (e * key.t_mod_) + plaintext;

  // Compute c1 = -a and return the ciphertext.
  return SymmetricRlweCiphertext<ModularInt>(
      std::vector<NttPolynomial<ModularInt>>{c0, -a}, key.modulus_params_);
}

// Takes as input the result of decrypting a RLWE plaintext that still contains
// the error. Concretely, it contains m + e * t (mod q). This function
// eliminates the error and returns the message. For reasons described below,
// this operation is more complicated than a simple (mod t).
//
// The error is drawn from a binomial distribution centered at zero and
// multiplied by t, meaning error values are either positive or negative
// multiples of t. Since each coefficient of the plaintext is smaller than
// t, some coefficients of the quantity m + e * t (which is all that's
// left in the vector error_and_message) could be negative. We are using
// modular arithmetic, so negative values become large positive values.
//
// Unfortunately, these negative values caues the naive error elimination
// strategy to fail. In theory we could take (m + e * t) mod t to
// eliminate the error portion and extract the message. However, consider
// a case where the error is negative. Suppose that t=2, m=1, and e=-1
// with a modulus q=7:
//
//    m +  e * t (mod q) =
//    1 + -1 * 2 (mod 7) =
//            -1 (mod 7) =
//             6 (mod 7)
//
// When we take 6 (mod t) = 6 (mod 2), we get 0, which is not the original
// bit of m. To avoid this problem, we treat negative values as negative
// values, not as their equivalents mod q.
//
// We consider (m + e * t) to be negative whenever it is between q/2
// and q. Recall that, if |m + e * t| is greater than q/2, decryption
// fails.
//
// When the quantity (m + e * t) (mod q) represents a negative number
// mod q, we can re-create its non-modular negative form by computing
// ((m + e * t) - q). We can then take this value mod t to extract the
// correct answer.
//
// 1. (m + e * t (mod q)) =                    // in the range [q/2, q)
// 2. (m + e * t - q)     =                    // in the range [-q/2, 0)
// 3. m (mod t) + e * t (mod t) - q (mod t) =  // taken (mod t)
// 4. m - (q (mod t))
//
// If we subtract q at step 2, we return negative numbers to their
// original form. Since we are going to perform a (mod t) operation
// anyway, we can subtract q (mod t) at step 2 to get the same result.
// Subtracting q (mod t) instead ensures that the quantity at step 2
// does not become negative, which is convenient because we are using
// an unsigned integer type.
//
// Concluding the example from before with the fix:
//
//    m +  e * t (mod q) - q (mod t) =
//    1 + -1 * 2 (mod 7) - 7 (mod 2) =
//            -1 (mod 7) - 7 (mod 2) = 6 - 1 = 5
//
// 5 (mod t) = 1, which is the original message.
template <typename ModularInt>
std::vector<typename ModularInt::Int> RemoveError(
    std::vector<ModularInt> error_and_message, typename ModularInt::Int q,
    typename ModularInt::Int t) {
  typename ModularInt::Int q_mod_t = q % t;
  std::vector<typename ModularInt::Int> plaintext(error_and_message.size());

  for (int i = 0; i < error_and_message.size(); i++) {
    plaintext[i] = error_and_message[i].ExportInt();

    if (plaintext[i] > q / 2) {
      plaintext[i] -= q_mod_t;
    }

    plaintext[i] %= t;
  }

  return plaintext;
}

template <typename ModularInt>
std::vector<typename ModularInt::Int> Decrypt(
    const SymmetricRlweKey<ModularInt>& key,
    const SymmetricRlweCiphertext<ModularInt>& ciphertext) {
  // Extract the error and message. To do so, take the dot product of the
  // ciphertext vector <c0, c1, ..., cN> and the vector of the powers of
  // the key <s^0, s^1, ..., s^N>.

  // Accumulator variables.
  NttPolynomial<ModularInt> error_and_message_ntt(key.key_.Len(),
                                                  key.modulus_params_);
  NttPolynomial<ModularInt> key_powers = key.key_;
  unsigned int ciphertext_len = ciphertext.c_.size();

  for (int i = 0; i < ciphertext_len; i++) {
    // Extract component i.
    NttPolynomial<ModularInt> ci = ciphertext.c_[i];

    // Lazily increase the exponent of the key.
    if (i > 1) {
      key_powers = key_powers * key.key_;
    }

    // Beyond c0, multiply the exponentiated key in.
    if (i > 0) {
      ci = ci * key_powers;
    }

    error_and_message_ntt = error_and_message_ntt + ci;
  }

  // Invert the NTT process.
  std::vector<ModularInt> error_and_message =
      error_and_message_ntt.InverseNtt(key.ntt_params_);

  // Extract the message.
  return RemoveError<ModularInt>(error_and_message,
                                 key.modulus_params_->modulus, 1 << key.log_t_);
}

}  // namespace rlwe

#endif  // RLWE_SYMMETRIC_ENCRYPTION_H_
