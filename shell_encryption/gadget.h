#ifndef RLWE_GADGET_H_
#define RLWE_GADGET_H_

#include <cstddef>
#include <vector>

#include "shell_encryption/status_macros.h"
#include "shell_encryption/statusor.h"

// Common methods for working with "gadgets". These are ways of representing
// large norm (but small dimension) objects with small norm (but large
// dimension) objects. For example, one can take an integer |x| <= q, and write
//
// x = \sum_{i = 0}^{\log_2 q} x_i 2^i
//
// where the x_i are small. One can of course generalize this to digits B != 2.
// The other main type of gadget is taking |x| <= \prod_i p_i for coprime values
// p_i, and mapping
//
// x -> (x mod p_0, x mod p_1, ..., x mod p_k)
//
// For a formal introduction to gadgets, see https://eprint.iacr.org/2018/946.

namespace rlwe {

// Method to compute the number of digits needed to represent integers mod
// q in base T.
template <typename ModularInt>
inline int GadgetSize(int log_base,
                      const typename ModularInt::Params* mod_params) {
  return (mod_params->log_modulus + (log_base - 1)) / log_base;
}

// Return the vector of base-B decomposition of each x mod q, where q is the
// modulus defined in mod_params, ie return a vector of [v_0, ..., v_{k-1}]
// such that sum(v_{j} * B^j) = x[i] mod q, and v_j \in [0, B)
template <typename ModularInt>
StatusOr<std::vector<std::vector<ModularInt>>> BaseDecompose(
    const std::vector<ModularInt>& coeffs,
    const typename ModularInt::Params* mod_params, const size_t log_base,
    int dimension) {
  // Determine the dimension, which is ceil(log_base(q))
  std::vector<typename ModularInt::Int> curr_digits(coeffs.size(), 0);
  std::transform(
      coeffs.begin(), coeffs.end(), curr_digits.begin(),
      [mod_params](ModularInt x) { return x.ExportInt(mod_params); });

  // Compute the mask to extract the log_base least significant
  // bits
  typename ModularInt::Int mask =
      (static_cast<typename ModularInt::Int>(1) << log_base) - 1;
  std::vector<std::vector<ModularInt>> result(dimension);
  for (int i = 0; i < dimension; i++) {
    result[i].reserve(curr_digits.size());
    for (int j = 0; j < curr_digits.size(); ++j) {
      RLWE_ASSIGN_OR_RETURN(
          auto coefficient_part,
          ModularInt::ImportInt((curr_digits[j] & mask), mod_params));
      result[i].push_back(std::move(coefficient_part));
      curr_digits[j] >>= log_base;
    }
  }
  return result;
}

}  // namespace rlwe

#endif  // RLWE_GADGET_H_
