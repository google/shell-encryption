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

#ifndef RLWE_NTT_POLYNOMIAL_H_
#define RLWE_NTT_POLYNOMIAL_H_

#include <cmath>
#include <vector>
#include "ntt_parameters.h"
#include "serialization.pb.h"

namespace rlwe {

// A polynomial in NTT form. The length of the polynomial must be a power of 2.
template <typename ModularInt>
class NttPolynomial {
 public:
  // Default constructor.
  NttPolynomial() = default;

  // Copy constructor.
  NttPolynomial(const NttPolynomial& p) = default;
  NttPolynomial& operator=(const NttPolynomial& that) = default;

  // Create an empty polynomial of the specified length. The length must be
  // a power of 2.
  explicit NttPolynomial(int len, const typename ModularInt::Params* params)
      : NttPolynomial(
            std::vector<ModularInt>(len, ModularInt::ImportInt(params, 0))) {}

  // This is an implementation of the iterative FFT from CLRS 2nd ed. Ch30.3.
  // with the addition of the psi-trick as described in (1) below.
  //
  // This factory method takes as an argument a vector of coefficients for a
  // polynomial of degree n-1 (where n is a power of 2) with coefficients mod
  // modulus. This polynomial is converted into NTT form. Entry i of this
  // vector represents the coefficient for the x^i term of the polynomial it
  // represents.
  //
  // All polynomial arithmetic performed while in NTT form is mod (x^n + 1),
  // with the coefficients operated on mod modulus.
  //
  // The function also takes as arguments three tables of parameters that are
  // necessary to perform the NTT transformation. These parameters are
  // described in detail in the inline comments in the body of the function.
  static NttPolynomial ConvertToNtt(
      const std::vector<ModularInt>& poly_coeffs,
      const NttParameters<ModularInt>& ntt_params) {
    // Check to ensure that the coefficient vector is of the correct length.
    int len = poly_coeffs.size();
    if (len <= 0 || (len & (len - 1)) != 0) {
      // An error value.
      return NttPolynomial();
    }

    NttPolynomial output(poly_coeffs);

    // For the purpose of the exposition that follows, let D be the degree of
    // the polynomial (i.e., the number of coefficients this object stores).
    // Implementation-wise, D is len.

    // 1) Multiply by powers of psi. Psi is the (2*D)-th primitive root of
    //    unity. We multiply element i of the polynomial by psi^i. According to
    //    "Software Speed Records for Lattice-Based Encryption," multiplying
    //    by these values before doing the NTT conversion:
    //       "avoids the doubling of the input length of the NTT and also
    //        gives us a modular reduction by (x^n + 1) for free."
    //
    //    For the purposes of this function, we assume that psis
    //    is a one-dimensional array of the powers of psi for each corresponding
    //    value in coeffs_, i.e., psis[i] = psi^i
    output.InPlaceCoordinatewiseMultiplyHelper(ntt_params.psis);

    // 2) Use the bitrev array to put the elements of the polynomial in the
    //    proper order for performing NTT.
    //    Each element of the bitrev array is the bitwise reverse of the
    //    corresponding index. For example, in a 16-item array:
    //        bitrev[0]  = 0  (0000 in binary -> 0000 reversed = 0)
    //        bitrev[4]  = 2  (0100 in binary -> 0010 reversed = 2)
    //        bitrev[14] = 11 (1101 in binary -> 1011 reversed = 11)
    //    Note that these mappings are symmetric, meaning that each operation of
    //    this transformation can be implemented as a swap.
    //
    //    This bit reversal puts the elements of the array in the proper order
    //    for an iterative, bottom-up implementation of NTT. In other words,
    //    this permutation puts the array in the same order as the leaves of
    //    the recursion tree for the recursive version of the NTT algorithm.
    BitrevHelper(ntt_params.bitrevs, &output.coeffs_);

    // 3) The NTT transformation itself. This procedure exactly follows the
    //    algorithm described in section of 30.3 of CLRS 2nd edition.
    //
    //    Each item in the array omegas is a power of
    //    a primitive root of unity (mod N). Specifically, omegas[i][j] is
    //    the (2^(i+1))-th primitive root of unity to the j-th power.
    output.NttTransformHelper(ntt_params.omegas);

    return output;
  }

  // Perform the inverse NTT transform, returning the coefficients of a normal
  // polynomial corresponding to the coefficients of the NttPolynomial.
  //
  // The inv_omegas and inv_psis tables contain the inverses of the values found
  // in the omegas and psis tables of the ConvertToNtt method. Each entry of
  // the inv_psis table has also been multiplied by N inverse. The bitrev
  // table is the same as in the ConvertToNtt function.
  std::vector<ModularInt> InverseNtt(
      const NttParameters<ModularInt>& ntt_params) const {
    NttPolynomial cpy = NttPolynomial(*this);

    // Run the bit reversal to prepare for the inverse NTT transform.
    BitrevHelper(ntt_params.bitrevs, &cpy.coeffs_);

    // Run the NTT transform with the inverses of the omegas.
    cpy.NttTransformHelper(ntt_params.omegas_inv);

    // CLRS specifies that, after performing the inverse NTT transformation,
    // we need to divide each entry by N. Since each entry of the inv_psis
    // table has been multiplied by N inverse, we accomplish this multiplication
    // in the next step.

    // Coordinate-wise multiply by the inverses of the psis.
    cpy.InPlaceCoordinatewiseMultiplyHelper(ntt_params.psis_inv);

    return cpy.coeffs_;
  }

  // Specifies whether the NttPolynomial is valid.
  bool IsValid() const { return !coeffs_.empty(); }

  // Scalar multiply.
  NttPolynomial operator*(const ModularInt& scalar) const {
    NttPolynomial output = *this;

    for (int i = 0; i < Len(); i++) {
      output.coeffs_[i] = output.coeffs_[i] * scalar;
    }

    return output;
  }

  // Coordinate-wise multiplication.
  NttPolynomial operator*(const NttPolynomial& that) const {
    // If this operation is invalid, return an invalid result.
    if (!IsValid() || !that.IsValid() || Len() != that.Len()) {
      return NttPolynomial();
    }

    // Create the output polynomial.
    NttPolynomial output = *this;

    // Perform the multiplication.
    for (int i = 0; i < Len(); i++) {
      output.coeffs_[i] = output.coeffs_[i] * that.coeffs_[i];
    }

    return output;
  }

  // Negation.
  NttPolynomial operator-() const {
    NttPolynomial output = *this;

    for (int i = 0; i < Len(); i++) {
      output.coeffs_[i] = -output.coeffs_[i];
    }

    return output;
  }

  // Coordinate-wise addition.
  NttPolynomial operator+(const NttPolynomial& that) const {
    // If this operation is invalid, return an invalid result.
    if (!IsValid() || !that.IsValid() || Len() != that.Len()) {
      return NttPolynomial();
    }

    // Create the output polynomial.
    NttPolynomial output = *this;

    // Perform the addition.
    for (int i = 0; i < Len(); i++) {
      output.coeffs_[i] = output.coeffs_[i] + that.coeffs_[i];
    }

    return output;
  }

  // Boolean comparison.
  bool operator==(const NttPolynomial& that) const {
    if (Len() != that.Len()) {
      return false;
    }

    for (int i = 0; i < Len(); i++) {
      if (coeffs_[i] != that.coeffs_[i]) {
        return false;
      }
    }

    return true;
  }
  bool operator!=(const NttPolynomial& that) const { return !(*this == that); }

  int Len() const { return coeffs_.size(); }

  SerializedNttPolynomial Serialize() const {
    SerializedNttPolynomial output;

    for (const ModularInt& coeff : coeffs_) {
      *output.add_coeffs() = coeff.Serialize();
    }

    return output;
  }

  static NttPolynomial Deserialize(
      const typename ModularInt::Params* modular_params,
      const SerializedNttPolynomial& serialized) {
    NttPolynomial output(serialized.coeffs_size(), modular_params);

    for (int i = 0; i < serialized.coeffs_size(); i++) {
      output.coeffs_[i] =
          ModularInt::Deserialize(modular_params, serialized.coeffs(i));
    }

    return output;
  }

 private:
  // Instance variables.
  size_t log_len_;
  std::vector<ModularInt> coeffs_;

  // Basic constructor.
  explicit NttPolynomial(const std::vector<ModularInt>& poly_coeffs)
      : log_len_(log2(poly_coeffs.size())),
        coeffs_(std::vector<ModularInt>(poly_coeffs)) {}

  // Helper function: Perform the bit-reversal operation in-place on coeffs_.
  static void BitrevHelper(const std::vector<unsigned int>& bitrevs,
                           std::vector<ModularInt>* item_to_reverse) {
    using std::swap;
    for (int i = 0; i < item_to_reverse->size(); i++) {
      unsigned int r = bitrevs[i];

      // Only swap in one direction - don't accidentally swap twice.
      if (i < r) {
        swap((*item_to_reverse)[i], (*item_to_reverse)[r]);
      }
    }
  }

  // Helper function: Multiply each element in coeffs_ by the corresponding
  // element in a and store the result in coeffs_.
  void InPlaceCoordinatewiseMultiplyHelper(const std::vector<ModularInt>& a) {
    for (int i = 0; i < Len(); i++) {
      coeffs_[i] = coeffs_[i] * a[i];
    }
  }

  // Helper function: Perform the NTT transform itself.
  void NttTransformHelper(const std::vector<std::vector<ModularInt>>& omegas) {
    for (int i = 1; i <= log_len_; i++) {
      unsigned int half_m = 1 << (i - 1);
      unsigned int m = half_m << 1;

      for (int k = 0; k < Len(); k += m) {
        int power_of_root_of_unity = 0;

        for (int j = 0; j < half_m; j++) {
          ModularInt w = omegas[i - 1][power_of_root_of_unity];

          // The butterfly operation.
          ModularInt t = w * coeffs_[k + j + half_m];
          ModularInt u = coeffs_[k + j];
          coeffs_[k + j] = u + t;
          coeffs_[k + j + half_m] = u - t;

          power_of_root_of_unity++;
        }
      }
    }
  }
};

}  // namespace rlwe

#endif  // RLWE_NTT_POLYNOMIAL_H_
