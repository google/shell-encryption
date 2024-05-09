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

#ifndef RLWE_RNS_RNS_POLYNOMIAL_HWY_H_
#define RLWE_RNS_RNS_POLYNOMIAL_HWY_H_

#include "absl/types/span.h"
#include "hwy/aligned_allocator.h"
#include "shell_encryption/montgomery.h"

namespace rlwe {
namespace internal {

// Given vectors `a` and `b`, computes component-wise a[i] * b[i] and add them
// to output[i]. Modular multiplications are computed over the Montgomery form
// without reduction (and store the result in double width BigInt type), and
// are implemented using SIMD instructions via the highway library.
// Note: We assume that `a` and `b` have the same size, and `output` has at
// least the same number of elements.
template <typename Integer>
void BatchFusedMulAddMontgomeryRep(
    absl::Span<const MontgomeryInt<Integer>> a,
    absl::Span<const MontgomeryInt<Integer>> b,
    hwy::AlignedVector<typename BigInt<Integer>::value_type>& output);

// Given vectors `a` and `b`, computes component-wise a[i] * b[i] and add them
// to output[i]. This version does not use highway SIMD intrinsics.
// Note: We assume that `a` and `b` have the same size, and `output` has at
// least the same number of elements.
template <typename Integer>
void BatchFusedMulAddMontgomeryRepNoHwy(
    absl::Span<const MontgomeryInt<Integer>> a,
    absl::Span<const MontgomeryInt<Integer>> b,
    hwy::AlignedVector<typename BigInt<Integer>::value_type>& output);

}  // namespace internal
}  // namespace rlwe

#endif  // RLWE_RNS_RNS_POLYNOMIAL_HWY_H_
