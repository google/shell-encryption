// Copyright 2024 Google LLC
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include "shell_encryption/rns/rns_polynomial_hwy.h"

#include <cstdint>
#include <vector>

#include "absl/numeric/int128.h"
#include "absl/types/span.h"
#include "hwy/detect_targets.h"
#include "shell_encryption/integral_types.h"
#include "shell_encryption/montgomery.h"

// Highway implementations.
// clang-format off
#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "shell_encryption/rns/rns_polynomial_hwy.cc"
#include "hwy/foreach_target.h"  // IWYU pragma: keep
// clang-format on

// Must come after foreach_target.h to avoid redefinition errors.
#include "hwy/aligned_allocator.h"
#include "hwy/highway.h"

HWY_BEFORE_NAMESPACE();
namespace rlwe::internal {
namespace HWY_NAMESPACE {

namespace hn = hwy::HWY_NAMESPACE;

#if HWY_TARGET == HWY_SCALAR

template <typename Integer>
void BatchFusedMulAddMontgomeryRepHwy(
    absl::Span<const MontgomeryInt<Integer>> a,
    absl::Span<const MontgomeryInt<Integer>> b,
    hwy::AlignedVector<typename BigInt<Integer>::value_type>& output) {
  BatchFusedMulAddMontgomeryRepNoHwy(a, b, output);
}

template <typename Integer>
void BatchFusedMulSumAddMontgomeryRepHwy(
    absl::Span<const MontgomeryInt<Integer>> a,
    absl::Span<const MontgomeryInt<Integer>> b,
    absl::Span<const MontgomeryInt<Integer>> c,
    hwy::AlignedVector<typename BigInt<Integer>::value_type>& output) {
  BatchFusedMulSumAddMontgomeryRepNoHwy(a, b, c, output);
}

template <typename Integer>
void BatchFusedMulDifferenceAddMontgomeryRepHwy(
    absl::Span<const MontgomeryInt<Integer>> a,
    absl::Span<const MontgomeryInt<Integer>> b,
    absl::Span<const MontgomeryInt<Integer>> c, Integer q,
    hwy::AlignedVector<typename BigInt<Integer>::value_type>& output) {
  BatchFusedMulDifferenceAddMontgomeryRepNoHwy(a, b, c, q, output);
}

#else

// General version that should work for `Integer` smaller than 64-bit unsigned
// types.
template <typename Integer>
void BatchFusedMulAddMontgomeryRepHwy(
    absl::Span<const MontgomeryInt<Integer>> a,
    absl::Span<const MontgomeryInt<Integer>> b,
    hwy::AlignedVector<typename BigInt<Integer>::value_type>& output) {
  using BigInteger = typename BigInt<Integer>::value_type;
  using D = hn::ScalableTag<Integer>;
  using Wide = hwy::MakeWide<hn::TFromD<D>>;
  const D d;                          // For lanes containing Integers.
  const hn::Repartition<Wide, D> d2;  // For double width products of Integers.

  const int N = hn::Lanes(d);  // This is guaranteed to be a power of two.
  const int64_t num_coeffs = a.size();
  int64_t i = 0;
  for (; i + N <= num_coeffs; i += N) {
    BigInteger* output0_ptr = &output[i];
    BigInteger* output1_ptr = &output[i] + N / 2;
    const Integer* a_ptr = reinterpret_cast<const Integer*>(&a[i]);
    const Integer* b_ptr = reinterpret_cast<const Integer*>(&b[i]);
    auto output0 = hn::Load(d2, output0_ptr);
    auto output1 = hn::Load(d2, output1_ptr);
    auto a_vec = hn::LoadU(d, a_ptr);  // a[i..i+N]
    auto b_vec = hn::LoadU(d, b_ptr);  // b[i..i+N]
    // Compute a[even_idx] * b[even_idx].
    auto mul_even = hn::MulEven(a_vec, b_vec);
    // Compute a[odd_idx] * b[odd_idx].
    auto mul_odd = hn::MulOdd(a_vec, b_vec);
    // Merge the lower blocks into mul0, and upper blocks into mul1.
    auto mul0 = hn::InterleaveWholeLower(d2, mul_even, mul_odd);
    auto mul1 = hn::InterleaveWholeUpper(d2, mul_even, mul_odd);
    // Add the lower and higher blocks of products to the output vector.
    hn::Store(hn::Add(output0, mul0), d2, output0_ptr);
    hn::Store(hn::Add(output1, mul1), d2, output1_ptr);
  }

  // For the remaining elements in `a` and `b`.
  for (; i < num_coeffs; ++i) {
    output[i] += static_cast<BigInteger>(a[i].GetMontgomeryRepresentation()) *
                 b[i].GetMontgomeryRepresentation();
  }
}

using ModularInt64 = MontgomeryInt<Uint64>;
using BigInt64 = ModularInt64::BigInt;

// Specialized version for 64-bit integers, as multiplication needs more careful
// treatment.
template <>
void BatchFusedMulAddMontgomeryRepHwy(absl::Span<const ModularInt64> a,
                                      absl::Span<const ModularInt64> b,
                                      hwy::AlignedVector<BigInt64>& output) {
  using D = hn::ScalableTag<Uint64>;
  const D d;
  const int N = hn::Lanes(d);
  const int64_t num_coeffs = a.size();

  int64_t i = 0;
  for (; i + N <= num_coeffs; i += N) {
    const Uint64* a_ptr = reinterpret_cast<const Uint64*>(&a[i]);
    const Uint64* b_ptr = reinterpret_cast<const Uint64*>(&b[i]);
    auto a_vec = hn::LoadU(d, a_ptr);
    auto b_vec = hn::LoadU(d, b_ptr);
    auto mul_lo = hn::Mul(a_vec, b_vec);
    auto mul_hi = hn::MulHigh(a_vec, b_vec);

    Uint64* output_ptr = reinterpret_cast<Uint64*>(&output[i]);
    hn::Vec<decltype(d)> out_lo, out_hi;
    hn::LoadInterleaved2(d, output_ptr, out_lo, out_hi);

    auto new_lo = hn::Add(out_lo, mul_lo);
    auto carry = hn::IfThenElseZero(hn::Lt(new_lo, out_lo), hn::Set(d, 1));
    auto new_hi = hn::Add(hn::Add(out_hi, mul_hi), carry);

    hn::StoreInterleaved2(new_lo, new_hi, d, output_ptr);
  }

  // Handle the remaining elements in the input vectors.
  for (; i < num_coeffs; ++i) {
    output[i] += static_cast<BigInt64>(a[i].GetMontgomeryRepresentation()) *
                 b[i].GetMontgomeryRepresentation();
  }
}

////////////////////////////////////////////////////////////////////////////////
// FMSumAdd
////////////////////////////////////////////////////////////////////////////////

template <typename Integer>
void BatchFusedMulSumAddMontgomeryRepHwy(
    absl::Span<const MontgomeryInt<Integer>> a,
    absl::Span<const MontgomeryInt<Integer>> b,
    absl::Span<const MontgomeryInt<Integer>> c,
    hwy::AlignedVector<typename BigInt<Integer>::value_type>& output) {
  using BigInteger = typename BigInt<Integer>::value_type;
  using D = hn::ScalableTag<Integer>;
  using Wide = hwy::MakeWide<hn::TFromD<D>>;
  const D d;                          // For lanes containing Integers.
  const hn::Repartition<Wide, D> d2;  // For double width products of Integers.

  const int N = hn::Lanes(d);  // This is guaranteed to be a power of two.
  const int64_t num_coeffs = a.size();
  int64_t i = 0;
  for (; i + N <= num_coeffs; i += N) {
    BigInteger* output0_ptr = &output[i];
    BigInteger* output1_ptr = &output[i] + N / 2;
    const Integer* a_ptr = reinterpret_cast<const Integer*>(&a[i]);
    const Integer* b_ptr = reinterpret_cast<const Integer*>(&b[i]);
    const Integer* c_ptr = reinterpret_cast<const Integer*>(&c[i]);
    auto output0 = hn::Load(d2, output0_ptr);
    auto output1 = hn::Load(d2, output1_ptr);
    auto a_vec = hn::LoadU(d, a_ptr);  // a[i..i+N]
    auto b_vec = hn::LoadU(d, b_ptr);  // b[i..i+N]
    auto c_vec = hn::LoadU(d, c_ptr);  // b[i..i+N]

    // Compute a[i] + b[i].
    auto ab_vec = hn::Add(a_vec, b_vec);
    // Compute ab[even_idx] * c[even_idx].
    auto mul_even = hn::MulEven(ab_vec, c_vec);
    // Compute ab[odd_idx] * c[odd_idx].
    auto mul_odd = hn::MulOdd(ab_vec, c_vec);

    // Merge the lower blocks into mul0, and upper blocks into mul1.
    auto mul0 = hn::InterleaveWholeLower(d2, mul_even, mul_odd);
    auto mul1 = hn::InterleaveWholeUpper(d2, mul_even, mul_odd);
    // Add the lower and higher blocks of products to the output vector.
    hn::Store(hn::Add(output0, mul0), d2, output0_ptr);
    hn::Store(hn::Add(output1, mul1), d2, output1_ptr);
  }

  // For the remaining elements in `a`, `b`, and `c`.
  for (; i < num_coeffs; ++i) {
    output[i] += static_cast<BigInteger>(a[i].GetMontgomeryRepresentation() +
                                         b[i].GetMontgomeryRepresentation()) *
                 c[i].GetMontgomeryRepresentation();
  }
}

using ModularInt64 = MontgomeryInt<Uint64>;
using BigInt64 = ModularInt64::BigInt;

// Specialized version for 64-bit integers, as multiplication needs more careful
// treatment.
template <>
void BatchFusedMulSumAddMontgomeryRepHwy(absl::Span<const ModularInt64> a,
                                         absl::Span<const ModularInt64> b,
                                         absl::Span<const ModularInt64> c,
                                         hwy::AlignedVector<BigInt64>& output) {
  using D = hn::ScalableTag<Uint64>;
  const D d;
  const int N = hn::Lanes(d);
  const int64_t num_coeffs = a.size();

  int64_t i = 0;

  for (; i + N <= num_coeffs; i += N) {
    const Uint64* a_ptr = reinterpret_cast<const Uint64*>(&a[i]);
    const Uint64* b_ptr = reinterpret_cast<const Uint64*>(&b[i]);
    const Uint64* c_ptr = reinterpret_cast<const Uint64*>(&c[i]);
    auto a_vec = hn::LoadU(d, a_ptr);
    auto b_vec = hn::LoadU(d, b_ptr);

    // Compute (a[i] + b[i]) * c[i], assuming no overflow.
    auto ab_vec = hn::Add(a_vec, b_vec);
    auto c_vec = hn::LoadU(d, c_ptr);
    auto mul_lo = hn::Mul(ab_vec, c_vec);
    auto mul_hi = hn::MulHigh(ab_vec, c_vec);

    // Load the lower and upper uint64_t parts of uint128 values from `output`.
    Uint64* output_ptr = reinterpret_cast<Uint64*>(&output[i]);
    hn::Vec<decltype(d)> output_lo, output_hi;
    hn::LoadInterleaved2(d, output_ptr, output_lo, output_hi);

    auto output_lo_new = hn::Add(output_lo, mul_lo);

    auto carry =
        hn::IfThenElseZero(hn::Lt(output_lo_new, output_lo), hn::Set(d, 1));
    auto output_hi_new = hn::Add(output_hi, mul_hi);
    output_hi_new = hn::Add(output_hi_new, carry);

    hn::StoreInterleaved2(output_lo_new, output_hi_new, d, output_ptr);
  }

  // Handle the remaining elements in the input vectors.
  for (; i < num_coeffs; ++i) {
    output[i] += static_cast<BigInt64>(a[i].GetMontgomeryRepresentation() +
                                       b[i].GetMontgomeryRepresentation()) *
                 c[i].GetMontgomeryRepresentation();
  }
}

////////////////////////////////////////////////////////////////////////////////
// FMDifferenceAdd
////////////////////////////////////////////////////////////////////////////////

template <typename Integer>
void BatchFusedMulDifferenceAddMontgomeryRepHwy(
    absl::Span<const MontgomeryInt<Integer>> a,
    absl::Span<const MontgomeryInt<Integer>> b,
    absl::Span<const MontgomeryInt<Integer>> c, Integer q,
    hwy::AlignedVector<typename BigInt<Integer>::value_type>& output) {
  using BigInteger = typename BigInt<Integer>::value_type;
  using D = hn::ScalableTag<Integer>;
  using Wide = hwy::MakeWide<hn::TFromD<D>>;
  const D d;                          // For lanes containing Integers.
  const hn::Repartition<Wide, D> d2;  // For double width products of Integers.

  const int N = hn::Lanes(d);  // This is guaranteed to be a power of two.
  const int64_t num_coeffs = a.size();
  auto q_vec = hn::Set(d, q);
  int64_t i = 0;
  for (; i + N <= num_coeffs; i += N) {
    BigInteger* output0_ptr = &output[i];
    BigInteger* output1_ptr = &output[i] + N / 2;
    const Integer* a_ptr = reinterpret_cast<const Integer*>(&a[i]);
    const Integer* b_ptr = reinterpret_cast<const Integer*>(&b[i]);
    const Integer* c_ptr = reinterpret_cast<const Integer*>(&c[i]);
    auto output0 = hn::Load(d2, output0_ptr);
    auto output1 = hn::Load(d2, output1_ptr);
    auto a_vec = hn::LoadU(d, a_ptr);  // a[i..i+N]
    auto b_vec = hn::LoadU(d, b_ptr);  // b[i..i+N]
    auto c_vec = hn::LoadU(d, c_ptr);  // b[i..i+N]

    // Compute a[i] - b[i] (mod q).
    auto underflow_mask = hn::Lt(a_vec, b_vec);
    auto aa_vec = hn::IfThenElse(underflow_mask, hn::Add(a_vec, q_vec), a_vec);
    auto ab_vec = hn::Sub(aa_vec, b_vec);

    // Compute ab[even_idx] * c[even_idx].
    auto mul_even = hn::MulEven(ab_vec, c_vec);
    // Compute ab[odd_idx] * c[odd_idx].
    auto mul_odd = hn::MulOdd(ab_vec, c_vec);

    // Merge the lower blocks into mul0, and upper blocks into mul1.
    auto mul0 = hn::InterleaveWholeLower(d2, mul_even, mul_odd);
    auto mul1 = hn::InterleaveWholeUpper(d2, mul_even, mul_odd);
    // Add the lower and higher blocks of products to the output vector.
    hn::Store(hn::Add(output0, mul0), d2, output0_ptr);
    hn::Store(hn::Add(output1, mul1), d2, output1_ptr);
  }

  // For the remaining elements in `a`, `b`, and `c`.
  for (; i < num_coeffs; ++i) {
    output[i] +=
        static_cast<BigInteger>(a[i].GetMontgomeryRepresentation() + q -
                                b[i].GetMontgomeryRepresentation()) *
        c[i].GetMontgomeryRepresentation();
  }
}

using ModularInt64 = MontgomeryInt<Uint64>;
using BigInt64 = ModularInt64::BigInt;

// Specialized version for 64-bit integers, as multiplication needs more careful
// treatment.
template <>
void BatchFusedMulDifferenceAddMontgomeryRepHwy(
    absl::Span<const ModularInt64> a, absl::Span<const ModularInt64> b,
    absl::Span<const ModularInt64> c, Uint64 q,
    hwy::AlignedVector<BigInt64>& output) {
  using D = hn::ScalableTag<Uint64>;
  const D d;
  const int N = hn::Lanes(d);
  const int64_t num_coeffs = a.size();

  auto q_vec = hn::Set(d, q);

  int64_t i = 0;

  for (; i + N <= num_coeffs; i += N) {
    const Uint64* a_ptr = reinterpret_cast<const Uint64*>(&a[i]);
    const Uint64* b_ptr = reinterpret_cast<const Uint64*>(&b[i]);
    const Uint64* c_ptr = reinterpret_cast<const Uint64*>(&c[i]);
    auto a_vec = hn::LoadU(d, a_ptr);
    auto b_vec = hn::LoadU(d, b_ptr);

    // Compute a[i] - b[i] (mod q).
    auto underflow_mask = hn::Lt(a_vec, b_vec);
    auto aa_vec = hn::IfThenElse(underflow_mask, hn::Add(a_vec, q_vec), a_vec);
    auto ab_vec = hn::Sub(aa_vec, b_vec);

    // Compute (a[i] - b[i]) (mod q) * c[i]
    auto c_vec = hn::LoadU(d, c_ptr);
    auto mul_lo = hn::Mul(ab_vec, c_vec);
    auto mul_hi = hn::MulHigh(ab_vec, c_vec);

    // Load the lower and upper uint64_t parts of uint128 values from `output`.
    Uint64* output_ptr = reinterpret_cast<Uint64*>(&output[i]);
    hn::Vec<decltype(d)> output_lo, output_hi;
    hn::LoadInterleaved2(d, output_ptr, output_lo, output_hi);

    auto output_lo_new = hn::Add(output_lo, mul_lo);

    auto carry =
        hn::IfThenElseZero(hn::Lt(output_lo_new, output_lo), hn::Set(d, 1));
    auto output_hi_new = hn::Add(output_hi, mul_hi);
    output_hi_new = hn::Add(output_hi_new, carry);

    hn::StoreInterleaved2(output_lo_new, output_hi_new, d, output_ptr);
  }

  // Handle the remaining elements in the input vectors.
  for (; i < num_coeffs; ++i) {
    output[i] += static_cast<BigInt64>(a[i].GetMontgomeryRepresentation() + q -
                                       b[i].GetMontgomeryRepresentation()) *
                 c[i].GetMontgomeryRepresentation();
  }
}

#endif  // HWY_TARGET == HWY_SCALAR

}  // namespace HWY_NAMESPACE
}  // namespace rlwe::internal
HWY_AFTER_NAMESPACE();

#if HWY_ONCE || HWY_IDE
namespace rlwe::internal {

template <typename Integer>
void BatchFusedMulAddMontgomeryRepNoHwy(
    absl::Span<const MontgomeryInt<Integer>> a,
    absl::Span<const MontgomeryInt<Integer>> b,
    hwy::AlignedVector<typename BigInt<Integer>::value_type>& output) {
  using BigInteger = typename BigInt<Integer>::value_type;
  for (int j = 0; j < a.size(); ++j) {
    output[j] += static_cast<BigInteger>(a[j].GetMontgomeryRepresentation()) *
                 b[j].GetMontgomeryRepresentation();
  }
}

// For now we only instantiate the Uint32 and the specialized Uint64 versions,
// as they are the most common integer types used in RNS RLWE schemes.
HWY_EXPORT_T(BatchFusedMulAddMontgomeryRepHwy32,
             BatchFusedMulAddMontgomeryRepHwy<Uint32>);
HWY_EXPORT_T(BatchFusedMulAddMontgomeryRepHwy64,
             BatchFusedMulAddMontgomeryRepHwy<Uint64>);

template <typename T>
void BatchFusedMulAddMontgomeryRep(
    absl::Span<const MontgomeryInt<T>> a, absl::Span<const MontgomeryInt<T>> b,
    hwy::AlignedVector<typename BigInt<T>::value_type>& output) {
  BatchFusedMulAddMontgomeryRepNoHwy(a, b, output);
}

// Specialized instantiation to use the highway version of FMA computation.
template <>
void BatchFusedMulAddMontgomeryRep(
    absl::Span<const MontgomeryInt<Uint32>> a,
    absl::Span<const MontgomeryInt<Uint32>> b,
    hwy::AlignedVector<BigInt<Uint32>::value_type>& output) {
  HWY_DYNAMIC_DISPATCH_T(BatchFusedMulAddMontgomeryRepHwy32)(a, b, output);
}

// Specialized instantiation to use the highway version of FMA computation.
template <>
void BatchFusedMulAddMontgomeryRep(
    absl::Span<const MontgomeryInt<Uint64>> a,
    absl::Span<const MontgomeryInt<Uint64>> b,
    hwy::AlignedVector<BigInt<Uint64>::value_type>& output) {
  HWY_DYNAMIC_DISPATCH_T(BatchFusedMulAddMontgomeryRepHwy64)(a, b, output);
}

////////////////////////////////////////////////////////////////////////////////
// FMSumAdd
////////////////////////////////////////////////////////////////////////////////

template <typename Integer>
void BatchFusedMulSumAddMontgomeryRepNoHwy(
    absl::Span<const MontgomeryInt<Integer>> a,
    absl::Span<const MontgomeryInt<Integer>> b,
    absl::Span<const MontgomeryInt<Integer>> c,
    hwy::AlignedVector<typename BigInt<Integer>::value_type>& output) {
  using BigInteger = typename BigInt<Integer>::value_type;
  for (int j = 0; j < a.size(); ++j) {
    output[j] += static_cast<BigInteger>(a[j].GetMontgomeryRepresentation() +
                                         b[j].GetMontgomeryRepresentation()) *
                 c[j].GetMontgomeryRepresentation();
  }
}

// For now we only instantiate the Uint32 and the specialized Uint64 versions,
// as they are the most common integer types used in RNS RLWE schemes.
HWY_EXPORT_T(BatchFusedMulSumAddMontgomeryRepHwy32,
             BatchFusedMulSumAddMontgomeryRepHwy<Uint32>);
HWY_EXPORT_T(BatchFusedMulSumAddMontgomeryRepHwy64,
             BatchFusedMulSumAddMontgomeryRepHwy<Uint64>);

template <typename T>
void BatchFusedMulSumAddMontgomeryRep(
    absl::Span<const MontgomeryInt<T>> a, absl::Span<const MontgomeryInt<T>> b,
    absl::Span<const MontgomeryInt<T>> c,
    hwy::AlignedVector<typename BigInt<T>::value_type>& output) {
  BatchFusedMulSumAddMontgomeryRepNoHwy(a, b, c, output);
}

// Specialized instantiation to use the highway version of FMA computation.
template <>
void BatchFusedMulSumAddMontgomeryRep(
    absl::Span<const MontgomeryInt<Uint32>> a,
    absl::Span<const MontgomeryInt<Uint32>> b,
    absl::Span<const MontgomeryInt<Uint32>> c,
    hwy::AlignedVector<BigInt<Uint32>::value_type>& output) {
  HWY_DYNAMIC_DISPATCH_T(BatchFusedMulSumAddMontgomeryRepHwy32)(a, b, c,
                                                                output);
}

// Specialized instantiation to use the highway version of FMA computation.
template <>
void BatchFusedMulSumAddMontgomeryRep(
    absl::Span<const MontgomeryInt<Uint64>> a,
    absl::Span<const MontgomeryInt<Uint64>> b,
    absl::Span<const MontgomeryInt<Uint64>> c,
    hwy::AlignedVector<BigInt<Uint64>::value_type>& output) {
  HWY_DYNAMIC_DISPATCH_T(BatchFusedMulSumAddMontgomeryRepHwy64)(a, b, c,
                                                                output);
}

////////////////////////////////////////////////////////////////////////////////
// FMDifferenceAdd
////////////////////////////////////////////////////////////////////////////////

template <typename Integer>
void BatchFusedMulDifferenceAddMontgomeryRepNoHwy(
    absl::Span<const MontgomeryInt<Integer>> a,
    absl::Span<const MontgomeryInt<Integer>> b,
    absl::Span<const MontgomeryInt<Integer>> c, Integer q,
    hwy::AlignedVector<typename BigInt<Integer>::value_type>& output) {
  using BigInteger = typename BigInt<Integer>::value_type;
  for (int j = 0; j < a.size(); ++j) {
    output[j] +=
        static_cast<BigInteger>(a[j].GetMontgomeryRepresentation() + q -
                                b[j].GetMontgomeryRepresentation()) *
        c[j].GetMontgomeryRepresentation();
  }
}

// For now we only instantiate the Uint32 and the specialized Uint64 versions,
// as they are the most common integer types used in RNS RLWE schemes.
HWY_EXPORT_T(BatchFusedMulDifferenceAddMontgomeryRepHwy32,
             BatchFusedMulDifferenceAddMontgomeryRepHwy<Uint32>);
HWY_EXPORT_T(BatchFusedMulDifferenceAddMontgomeryRepHwy64,
             BatchFusedMulDifferenceAddMontgomeryRepHwy<Uint64>);

template <typename T>
void BatchFusedMulDifferenceAddMontgomeryRep(
    absl::Span<const MontgomeryInt<T>> a, absl::Span<const MontgomeryInt<T>> b,
    absl::Span<const MontgomeryInt<T>> c, T q,
    hwy::AlignedVector<typename BigInt<T>::value_type>& output) {
  BatchFusedMulDifferenceAddMontgomeryRepNoHwy(a, b, c, q, output);
}

// Specialized instantiation to use the highway version of FMA computation.
template <>
void BatchFusedMulDifferenceAddMontgomeryRep(
    absl::Span<const MontgomeryInt<Uint32>> a,
    absl::Span<const MontgomeryInt<Uint32>> b,
    absl::Span<const MontgomeryInt<Uint32>> c, Uint32 q,
    hwy::AlignedVector<BigInt<Uint32>::value_type>& output) {
  HWY_DYNAMIC_DISPATCH_T(BatchFusedMulDifferenceAddMontgomeryRepHwy32)(
      a, b, c, q, output);
}

// Specialized instantiation to use the highway version of FMA computation.
template <>
void BatchFusedMulDifferenceAddMontgomeryRep(
    absl::Span<const MontgomeryInt<Uint64>> a,
    absl::Span<const MontgomeryInt<Uint64>> b,
    absl::Span<const MontgomeryInt<Uint64>> c, Uint64 q,
    hwy::AlignedVector<BigInt<Uint64>::value_type>& output) {
  HWY_DYNAMIC_DISPATCH_T(BatchFusedMulDifferenceAddMontgomeryRepHwy64)(
      a, b, c, q, output);
}

}  // namespace rlwe::internal
#endif  // HWY_ONCE || HWY_IDE
