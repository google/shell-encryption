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

  // Generate the masks on the even lanes, which correspond to the lower 64 bits
  // of BigInt64 (unsigned 128-bit int) values in the output vector.
  auto mask_lo = hn::Eq(hn::And(hn::Iota(d, 0), hn::Set(d, 1)), hn::Zero(d));
  // A highway vector whose odd lanes (higher 64 bits of the output vector) are
  // all 1's and even lanes are all 0's.
  auto ones = hn::Slide1Up(d, hn::IfThenElseZero(mask_lo, hn::Set(d, 1)));

  int64_t i = 0;

  for (; i + N <= num_coeffs; i += N) {
    // Load the lower N/2 and upper N/2 values from `output`.
    Uint64* output0_ptr = reinterpret_cast<Uint64*>(&output[i]);
    Uint64* output1_ptr = reinterpret_cast<Uint64*>(&output[i]) + N;
    auto output0 = hn::Load(d, output0_ptr);
    auto output1 = hn::Load(d, output1_ptr);

    const Uint64* a_ptr = reinterpret_cast<const Uint64*>(&a[i]);
    const Uint64* b_ptr = reinterpret_cast<const Uint64*>(&b[i]);
    auto a_vec = hn::LoadU(d, a_ptr);
    auto b_vec = hn::LoadU(d, b_ptr);

    // `hn::MulEven` and `hn::MulOdd` on 64-bit inputs produce a vector of
    // 64-bit values with even lanes storing the lower half of product and odd
    // lanes storing the upper half. That is,
    //   mul_even = [ab[0].lo, ab[0].hi, ab[2].lo, ab[2].hi, ...]
    //   mul_odd  = [ab[1].lo, ab[1].hi, ab[3].lo, ab[3].hi, ...]
    // So we first merge the lower and upper halves, and then permute within
    // each block of four 64-bit values, and get
    //   mul0 = [ab[0].lo, ab[0].hi, ab[1].lo, ab[1].hi, ...]
    //   mul1 = [ab[N/2].lo, ab[N/2].hi, ab[N/2+1].lo, ab[N/2+1].hi, ...]
    auto mul_even = hn::MulEven(a_vec, b_vec);
    auto mul_odd = hn::MulOdd(a_vec, b_vec);
    auto mul_interleave0 = hn::InterleaveWholeLower(d, mul_even, mul_odd);
    auto mul_interleave1 = hn::InterleaveWholeUpper(d, mul_even, mul_odd);
    auto mul0 = hn::Per4LaneBlockShuffle<3, 1, 2, 0>(mul_interleave0);
    auto mul1 = hn::Per4LaneBlockShuffle<3, 1, 2, 0>(mul_interleave1);

    // Add the products to the output lanes. A carry bit occurs when the new
    // lower 64-bit value becomes smaller, and we add 1 to the upper 64-bit lane
    // if that happens.
    auto output0_new = hn::Add(output0, mul0);
    auto output1_new = hn::Add(output1, mul1);
    auto mask_carry0 = hn::And(hn::Lt(output0_new, output0), mask_lo);
    auto mask_carry1 = hn::And(hn::Lt(output1_new, output1), mask_lo);
    mask_carry0 = hn::SlideMask1Up(d, mask_carry0);  // move the mask to hi64
    mask_carry1 = hn::SlideMask1Up(d, mask_carry1);  // move the mask to hi64
    output0 = hn::MaskedAddOr(output0_new, mask_carry0, output0_new, ones);
    output1 = hn::MaskedAddOr(output1_new, mask_carry1, output1_new, ones);

    hn::Store(output0, d, output0_ptr);
    hn::Store(output1, d, output1_ptr);
  }

  // Handle the remaining elements in the input vectors.
  for (; i < num_coeffs; ++i) {
    output[i] += static_cast<BigInt64>(a[i].GetMontgomeryRepresentation()) *
                 b[i].GetMontgomeryRepresentation();
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

}  // namespace rlwe::internal
#endif  // HWY_ONCE || HWY_IDE
