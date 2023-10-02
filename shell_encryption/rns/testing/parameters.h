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

#ifndef RLWE_RNS_TESTING_PARAMETERS_H_
#define RLWE_RNS_TESTING_PARAMETERS_H_

#include <gmock/gmock.h>
#include "shell_encryption/integral_types.h"
#include "shell_encryption/montgomery.h"

namespace rlwe {
namespace testing {

// Useful typedefs.
using ModularInt16 = MontgomeryInt<Uint16>;
using ModularInt32 = MontgomeryInt<Uint32>;
using ModularInt64 = MontgomeryInt<Uint64>;
using ModularInt128 = MontgomeryInt<Uint128>;

// Constants used for test purpose only.
constexpr Uint16 kModulus29Q0 = 12289;  // 14-bits
constexpr Uint16 kModulus29P0 = 15361;  // 14-bits
constexpr int kLogN29 = 9;
constexpr unsigned int kPlaintextModulus29 = 8;  // 12289 % 8 == 15361 % 8

// Plaintext modulus for other sets of prime moduli.
constexpr unsigned int kPlaintextModulus = (1 << 3) + 1;

// The following sets of moduli are chosen to be 1 modulo kPlaintextModulus
// such that they introduce minimal errors during modulus switching.
constexpr Uint32 kModulus60Q0 = 995329;   // 20-bits
constexpr Uint32 kModulus60Q1 = 921601;   // 20-bits
constexpr Uint32 kModulus60P0 = 1032193;  // 20-bits
constexpr int kLogN60 = 10;

constexpr Uint64 kModulus120Q0 = 1073442817ULL;  // 30-bits
constexpr Uint64 kModulus120Q1 = 1073184769ULL;  // 30-bits
constexpr Uint64 kModulus120Q2 = 1072705537ULL;  // 30-bits
constexpr Uint64 kModulus120P0 = 1073479681ULL;  // 30-bits
constexpr int kLogN120 = 11;

constexpr Uint128 kModulus300Q0 = static_cast<Uint128>(
    absl::MakeUint128(63, 18446744073708380161ULL));  // 70-bits
constexpr Uint128 kModulus300Q1 =
    static_cast<Uint128>(1099510824961ULL);  // 40-bits
constexpr Uint128 kModulus300Q2 =
    static_cast<Uint128>(1099510456321ULL);  // 40-bits
constexpr Uint128 kModulus300P0 =
    static_cast<Uint128>(1099508981761ULL);  // 40-bits
constexpr Uint128 kModulus300P1 =
    static_cast<Uint128>(1099508760577ULL);  // 40-bits
constexpr Uint128 kModulus300P2 =
    static_cast<Uint128>(1099508391937ULL);  // 40-bits
constexpr int kLogN300 = 12;

// Aggregate type that holds the parameters defining a RNS context.
template <typename ModularInt>
struct RnsParameters {
  int log_n;
  std::vector<typename ModularInt::Int> qs;  // main prime moduli.
  std::vector<typename ModularInt::Int> ps;  // auxiliary prime moduli.
  typename ModularInt::Int t;                // plaintext modulus.
};

// Returns the testing parameters for the underlying integer type.
template <typename ModularInt>
RnsParameters<ModularInt> GetRnsParameters();

template <>
inline RnsParameters<ModularInt16> GetRnsParameters<ModularInt16>() {
  return RnsParameters<ModularInt16>{
      .log_n = kLogN29, .qs = {kModulus29Q0}, .ps = {kModulus29P0}, .t = 8};
}

template <>
inline RnsParameters<ModularInt32> GetRnsParameters<ModularInt32>() {
  return RnsParameters<ModularInt32>{.log_n = kLogN60,
                                     .qs = {kModulus60Q0, kModulus60Q1},
                                     .ps = {kModulus60P0},
                                     .t = kPlaintextModulus};
}

template <>
inline RnsParameters<ModularInt64> GetRnsParameters<ModularInt64>() {
  return RnsParameters<ModularInt64>{
      .log_n = kLogN120,
      .qs = {kModulus120Q0, kModulus120Q1, kModulus120Q2},
      .ps = {kModulus120P0},
      .t = kPlaintextModulus};
}

template <>
inline RnsParameters<ModularInt128> GetRnsParameters<ModularInt128>() {
  return RnsParameters<ModularInt128>{
      .log_n = kLogN300,
      .qs = {kModulus300Q0, kModulus300Q1, kModulus300Q2},
      .ps = {kModulus300P0, kModulus300P1, kModulus300P2},
      .t = kPlaintextModulus};
}

#ifdef ABSL_HAVE_INTRINSIC_INT128
template <>
RnsParameters<MontgomeryInt<absl::uint128>> inline GetRnsParameters<
    MontgomeryInt<absl::uint128>>() {
  return RnsParameters<MontgomeryInt<absl::uint128>>{
      .log_n = kLogN300,
      .qs = {kModulus300Q0, kModulus300Q1, kModulus300Q2},
      .ps = {kModulus300P0, kModulus300P1, kModulus300P2},
      .t = kPlaintextModulus};
}
#endif

// We use the following structure (specialied for each integral type T)
// to define the product Q of prime moduli of integer type T defined above,
// and to define the integral type that is big enough to hold values mod Q.
template <typename T>
struct CompositeModulus;

// Specialization for uint16, uint32, uint64, and uint128.
template <>
struct CompositeModulus<Uint16> {
  using value_type = Uint32;
  static value_type Value() { return kModulus29Q0; }
};
template <>
struct CompositeModulus<Uint32> {
  using value_type = Uint64;
  static value_type Value() {
    return static_cast<value_type>(kModulus60Q0) *
           static_cast<value_type>(kModulus60Q1);
  }
};
template <>
struct CompositeModulus<Uint64> {
  using value_type = absl::uint128;
  static value_type Value() {
    return static_cast<value_type>(kModulus120Q0) *
           static_cast<value_type>(kModulus120Q1) *
           static_cast<value_type>(kModulus120Q2);
  }
};
template <>
struct CompositeModulus<absl::uint128> {
  using value_type = uint256;
  static value_type Value() {
    return static_cast<value_type>(kModulus300Q0) *
           static_cast<value_type>(kModulus300Q1) *
           static_cast<value_type>(kModulus300Q2);
  }
};
#ifdef ABSL_HAVE_INTRINSIC_INT128
template <>
struct CompositeModulus<unsigned __int128> {
  using value_type = uint256;
  static value_type Value() {
    return static_cast<value_type>(kModulus300Q0) *
           static_cast<value_type>(kModulus300Q1) *
           static_cast<value_type>(kModulus300Q2);
  }
};
#endif

////////////////////////////////////////////////////////////////////////////////
// Finite field encoding
////////////////////////////////////////////////////////////////////////////////

// Below we define some RNS parameters suitable for encoding finite field values
// into plaintext slots. We currently only support finite fields of prime order,
// i.e. Z_t for prime t, and plaintext encoding requires that t == 1 (mod 2*N).
template <typename ModularInt>
std::vector<RnsParameters<ModularInt>> GetRnsParametersForFiniteFieldEncoding();

// The following parameters are chosen such that q_i == 1 (mod 2*N) and
// q_i (mod t) = q_j (mod t) for all q_i, q_j, and furthermore t == 1 (mod 2*N).
// The test parameters are grouped by the underlying integer types, and then we
// have variants such as single main prime modulus and multiple prime moduli.
template <>
inline std::vector<RnsParameters<ModularInt32>>
GetRnsParametersForFiniteFieldEncoding<ModularInt32>() {
  return {
      // Two main prime moduli: q0 = 29 bits, q1 = 30 bits.
      RnsParameters<ModularInt32>{
          .log_n = 10, .qs = {335552513, 754993153}, .ps = {}, .t = 40961},
  };
}

template <>
inline std::vector<RnsParameters<ModularInt64>>
GetRnsParametersForFiniteFieldEncoding<ModularInt64>() {
  return {
      // A single main prime modulus: q0 = 40 bits.
      RnsParameters<ModularInt64>{
          .log_n = 11, .qs = {1095746727937ULL}, .ps = {}, .t = 40961},

      // Two main prime moduli: q0 = 33 bits, q1 = 33 bits.
      RnsParameters<ModularInt64>{.log_n = 11,
                                  .qs = {8556589057ULL, 8388812801ULL},
                                  .ps = {},
                                  .t = 40961},
  };
}

template <>
inline std::vector<RnsParameters<ModularInt128>>
GetRnsParametersForFiniteFieldEncoding<ModularInt128>() {
  return {
      // A single main prime modulus: q0 = 40 bits.
      RnsParameters<ModularInt128>{
          .log_n = 11,
          .qs = {static_cast<Uint128>(1125889168998401ULL)},
          .ps = {},
          .t = 65537},

      // Three main prime moduli: q0 = 35 bits, q1 = 35 bits, q2 = 35 bits.
      RnsParameters<ModularInt128>{.log_n = 12,
                                   .qs = {static_cast<Uint128>(28521963521ULL),
                                          static_cast<Uint128>(27179753473ULL),
                                          static_cast<Uint128>(26173095937ULL)},
                                   .ps = {},
                                   .t = 40961},
  };
}

#ifdef ABSL_HAVE_INTRINSIC_INT128
template <>
inline std::vector<RnsParameters<MontgomeryInt<absl::uint128>>>
GetRnsParametersForFiniteFieldEncoding<MontgomeryInt<absl::uint128>>() {
  return {
      RnsParameters<MontgomeryInt<absl::uint128>>{
          .log_n = 11,
          .qs = {static_cast<Uint128>(1125899644440577ULL)},
          .ps = {},
          .t = 40961},

      RnsParameters<MontgomeryInt<absl::uint128>>{
          .log_n = 12,
          .qs = {static_cast<Uint128>(28521963521ULL),
                 static_cast<Uint128>(27179753473ULL),
                 static_cast<Uint128>(26173095937ULL)},
          .ps = {},
          .t = 40961},
  };
}
#endif

// The smallest integer type that support finite field encoding is uint32_t.
#ifdef ABSL_HAVE_INTRINSIC_INT128
typedef ::testing::Types<
    rlwe::MontgomeryInt<Uint32>, rlwe::MontgomeryInt<Uint64>,
    rlwe::MontgomeryInt<absl::uint128>, rlwe::MontgomeryInt<unsigned __int128>>
    ModularIntTypesForFiniteFieldEncoding;
#else
typedef ::testing::Types<rlwe::MontgomeryInt<Uint32>,
                         rlwe::MontgomeryInt<Uint64>,
                         rlwe::MontgomeryInt<absl::uint128>>
    ModularIntTypesForFiniteFieldEncoding;
#endif

}  // namespace testing
}  // namespace rlwe

#endif  // RLWE_RNS_TESTING_PARAMETERS_H_
