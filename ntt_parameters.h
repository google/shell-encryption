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
#ifndef RLWE_NTT_PARAMETERS_H_
#define RLWE_NTT_PARAMETERS_H_

#include <algorithm>
#include <vector>

namespace rlwe {
namespace internal {

// Fill row with every power in {0, 1, ..., n-1} (mod modulus) of base .
template <typename ModularInt>
std::vector<ModularInt> FillWithEveryPower(
    const typename ModularInt::Params* params, ModularInt base,
    unsigned int n) {
  std::vector<ModularInt> row(n, ModularInt::ImportInt(params, 0));
  for (int i = 0; i < n; i++) {
    row[i] = base ^ i;
  }

  return row;
}

template <typename ModularInt>
ModularInt PrimitiveNthRootOfUnity(const typename ModularInt::Params* params,
                                   unsigned int log_n) {
  unsigned int n = 1 << log_n;
  unsigned int half_n = n >> 1;

  // The value k is a power such that any number raised to it will be a
  // n-th root of unity. (It will not necessarily be a *primitive* root of
  // unity, however).
  typename ModularInt::Int k = (params->modulus - 1) / n;

  // Test each number t to see whether t^k is a primitive 2n-th root
  // of unity - that t^{nk} is a root of unity but t^{(n/2)k} is not.
  ModularInt one = ModularInt::ImportInt(params, 1);
  for (typename ModularInt::Int t = 2; t < params->modulus; t++) {
    // Produce a candidate root of unity.
    ModularInt candidate = ModularInt::ImportInt(params, t) ^ k;

    // Check whether candidate^n = 1. If not, it is a primitive root of unity.
    if ((candidate ^ half_n) != one) {
      return candidate;
    }
  }

  // Failure state. The above loop should always return successfully assuming
  // the parameters were set properly.
  abort();
}

}  // namespace internal

// Creates the vector of psis. When performing the NTT transformation on a
// degree n-1 polynomial, psi is the primitive 2n-th root of unity. This
// function produces a vector of length n containing every power of psi between
// 0 and n-1. In other words, the vector produced by this function is:
//
//     {psi^0, psi^1, psi^2, ..., psi^(n-1)}
//
// As the first step of the NTT transformation,
// each element in the vector is multiplied by the corresponding element of the
// of the polynomial (i.e., psi^i * the coefficient for the x^i term).
//
// Each item of the vector is in modular integer representation.
template <typename ModularInt>
std::vector<ModularInt> NttPsis(const typename ModularInt::Params* params,
                                unsigned int log_n) {
  // Obtain a primitive 2n-th root of unity (hence log_n + 1).
  ModularInt w_2n =
      internal::PrimitiveNthRootOfUnity<ModularInt>(params, log_n + 1);
  unsigned int n = 1 << log_n;

  return internal::FillWithEveryPower<ModularInt>(params, w_2n, n);
}

// Creates a vector of the powers of psi^(-1), which is used in the inverse
// NTT transformation to invert the use of the psis from NttPsis. Each of these
// values is multiplied by the multiplicative inverse of n, which is necessary
// as part of the inverse NTT transformation.
template <typename ModularInt>
std::vector<ModularInt> NttPsisInv(const typename ModularInt::Params* params,
                                   unsigned int log_n) {
  unsigned int n = 1 << log_n;
  auto n_inv = ModularInt::ImportInt(params, n).MultiplicativeInverse();
  std::vector<ModularInt> row = NttPsis<ModularInt>(params, log_n);

  // Reverse the items at indices 1 through (n - 1). Multiplying index i
  // of the reversed row by index i of the original row will yield w_2n^n = -1.
  // (The exception is w_2n^0 = 1, which is already its own inverse.)
  std::reverse(row.begin() + 1, row.end());

  // Finally, multiply each of the items at indices 1 to (n-1) by -1. Multiply
  // every entry by n_inv.
  row[0] = row[0] * n_inv;
  for (int i = 1; i < n; i++) {
    row[i] = -row[i] * n_inv;
  }

  return row;
}

// Creates the table of omegas. Let omega_t be the primitive t-th root of unity
// for a given value of t. The j-th nested vector of this function's output
// table contains the powers from 0 to n-1 of (2^(j+1))-th primitive root of
// unity. In other words, the item at index (j, k) of the nested vectors is
// the (2^(j+1))-th primitive root of unity to the k^th power.
//
// Each item of the table is in modular integer representation.
//
// This table provides the omega values that are necessary to run the forward
// NTT transformation.
template <typename ModularInt>
std::vector<std::vector<ModularInt>> NttOmegas(
    const typename ModularInt::Params* params, unsigned int log_n) {
  unsigned int n = 1 << log_n;

  // Nested vectors that store a table of the primitive roots of unity and
  // their powers. Specifically, output[i][j] is the (2^(i+1))-th primitive root
  // of unity to the j-th power.
  std::vector<std::vector<ModularInt>> output(
      log_n, std::vector<ModularInt>(n, ModularInt::ImportInt(params, 0)));

  for (int row = log_n; row > 0; row--) {
    ModularInt w = internal::PrimitiveNthRootOfUnity<ModularInt>(params, row);
    output[row - 1] = internal::FillWithEveryPower<ModularInt>(params, w, n);
  }

  return output;
}

// Generates the omegas for the inverse  NTT transformation.
// Each row of the table will be the powers of the inverse
// of the corresponding primitive root of unity, producing the table necessary
// for the inverse NTT transformation.
//
// See NttOmegas for the structure of this table: index (i, j) of the output
// of this function is the multiplicative inverse of index (i, j) of the output
// of NttOmegas.
template <typename ModularInt>
std::vector<std::vector<ModularInt>> NttOmegasInv(
    const typename ModularInt::Params* params, unsigned int log_n) {
  // Retrieve the table for the forward transformation.
  std::vector<std::vector<ModularInt>> table =
      NttOmegas<ModularInt>(params, log_n);

  // The table for the reverse transformation is entries 1 through (n-1) in
  // reverse. The first entry in the table is 1, which is its own inverse.
  for (std::vector<ModularInt>& row : table) {
    std::reverse(row.begin() + 1, row.end());
  }

  return table;
}

// Creates a vector containing the indices necessary to perform the NTT bit
// reversal operation. Index i of the returned vector contains an integer with
// the rightmost log_n bits of i reversed.
std::vector<unsigned int> BitrevArray(unsigned int log_n);

// A struct that stores a package of NTT Parameters
template <typename ModularInt>
struct NttParameters {
  std::vector<ModularInt> psis;
  std::vector<ModularInt> psis_inv;
  std::vector<std::vector<ModularInt>> omegas;
  std::vector<std::vector<ModularInt>> omegas_inv;
  std::vector<unsigned int> bitrevs;
};

// A convenient function that sets up all NTT parameters at once.
// Does not take ownership of params.
template <typename ModularInt>
NttParameters<ModularInt> InitializeNttParameters(
    const typename ModularInt::Params* params, int log_n) {
  NttParameters<ModularInt> output;

  output.psis = NttPsis<ModularInt>(params, log_n);
  output.psis_inv = NttPsisInv<ModularInt>(params, log_n);
  output.omegas = NttOmegas<ModularInt>(params, log_n);
  output.omegas_inv = NttOmegasInv<ModularInt>(params, log_n);
  output.bitrevs = BitrevArray(log_n);

  return output;
}

}  // namespace rlwe

#endif  // RLWE_NTT_PARAMETERS_H_
