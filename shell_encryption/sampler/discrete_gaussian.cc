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

#include "shell_encryption/sampler/discrete_gaussian.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <memory>
#include <tuple>
#include <utility>
#include <vector>

#include "absl/memory/memory.h"
#include "absl/numeric/int128.h"
#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "absl/strings/str_cat.h"
#include "shell_encryption/integral_types.h"
#include "shell_encryption/prng/prng.h"
#include "shell_encryption/status_macros.h"

namespace rlwe {

namespace {

constexpr int kPrecision = 8 * sizeof(Uint64) - 1;

// Returns the index in a monotonically increasing vector `cdt` such that
// cdt[index - 1] < u <= cdt[index].
inline size_t FindInCdt(const std::vector<Uint64>& cdt, Uint64 u) {
  size_t index = 0;
  for (int i = 0; i < cdt.size(); ++i) {
    size_t increment = (cdt[i] < u);  // Don't increment once found u.
    index += increment;
  }
  return index;
}

}  // namespace

template <typename Integer>
absl::StatusOr<std::unique_ptr<DiscreteGaussianSampler<Integer>>>
DiscreteGaussianSampler<Integer>::Create(double s) {
  // Approximate sum(exp(-x^2 / (2 * s^2)), for x over the integers).
  double mass = sqrt(2 * M_PI) * s;

  // Initialize the CDF table for the base sampler. We cut off the base
  // distribution at the following `bound` as the probability mass beyond it is
  // negligible.
  int cut_off_bound = kTailBoundMultiplier * std::ceil(s);
  std::vector<double> cdf;
  cdf.reserve(2 * cut_off_bound + 1);
  for (int x = -cut_off_bound; x <= cut_off_bound; ++x) {
    double y = static_cast<double>(x) / s;
    double p = std::exp(-(y * y) / 2);  // Gaussian function at x.
    cdf.push_back(p / mass);
  }
  for (int i = 0; i < 2 * cut_off_bound; ++i) {
    cdf[i + 1] += cdf[i];
  }
  // Scale up the CDF values by 2^kPrecision and store the integer values.
  std::vector<Uint64> cdt;
  cdt.reserve(cdf.size());
  for (auto const& p : cdf) {
    cdt.push_back(static_cast<Uint64>(p * std::exp2(kPrecision)));
  }
  return absl::WrapUnique(
      new DiscreteGaussianSampler<Integer>(s, std::move(cdt)));
}

template <typename Integer>
absl::StatusOr<std::unique_ptr<DiscreteGaussianSampler<Integer>>>
DiscreteGaussianSampler<Integer>::CreateGeneric(double s_base) {
  if (s_base < std::sqrt(2) * kSmoothingParameter) {
    return absl::InvalidArgumentError(
        "`s_base` must be at least sqrt(2) times the smoothing parameter.");
  }

  // Create a base sampler.
  return Create(s_base);
}

template <typename Integer>
absl::StatusOr<Integer> DiscreteGaussianSampler<Integer>::SampleWithIterations(
    double s, int num_iterations, SecurePrng& prng) const {
  if (s < s_base_) {
    return absl::InvalidArgumentError(
        absl::StrCat("`s` must be at least the base s parameter ", s_base_));
  }
  if (num_iterations < 0) {
    return absl::InvalidArgumentError("`num_iterations` must be non-negative.");
  }

  // Use the base sampler if `num_iterations` == 0.
  if (num_iterations == 0) {
    return Sample(prng);
  }

  // Implements the arbitrary standard deviation discrete Gaussian sampling in
  // Algorithm 2 of https://eprint.iacr.org/2017/259. `num_iterations` is the
  // largest "i" such that s_i < s, where s_i is defined in Lemma 5.1.
  Integer x1, x2;
  double s_i;
  RLWE_ASSIGN_OR_RETURN(std::tie(x1, s_i),
                        SampleIIterative(prng, num_iterations));
  RLWE_ASSIGN_OR_RETURN(std::tie(x2, s_i),
                        SampleIIterative(prng, num_iterations));
  // We have t > 1 and z > 1 as s > s_i by assumption.
  double t = s / s_i;
  double z = std::ceil((std::sqrt(2 * t * t - 1) + 1) / 2);
  Integer y1 = static_cast<Integer>(z) * x1;
  Integer y2 = static_cast<Integer>(z - 1) * x2;
  Integer sample = y1 + y2;
  return sample;
}

template <typename Integer>
absl::StatusOr<Integer> DiscreteGaussianSampler<Integer>::Sample(
    SecurePrng& prng) const {
  // Sample 64 uniformly random bits. We use the last bit as the sign, and the
  // first 63 random bits to look up the scaled CDF table.
  constexpr Uint64 kMask = (1ULL << kPrecision) - 1;
  RLWE_ASSIGN_OR_RETURN(const Uint64 u, prng.Rand64());
  Integer index = static_cast<Integer>(FindInCdt(cdt_, u & kMask));
  Integer center_index = static_cast<Integer>(cdt_.size() / 2);
  Integer sample = index - center_index;
  return sample;
}

template <typename Integer>
absl::StatusOr<int> DiscreteGaussianSampler<Integer>::NumIterations(
    double s) const {
  if (s < s_base_) {
    return absl::InvalidArgumentError(
        absl::StrCat("`s` must be at least the base s parameter ", s_base_));
  }

  // Find the largest i such that s_i < s, where
  // z_i = floor(s_{i-1} / sqrt(2) * kSmoothParameter),
  // s_i = sqrt(z_i^2 + max((z_i - 1)^2, 1)) * s_{i-1},
  // s_0 = s_base.
  double denum = std::sqrt(2) * kSmoothingParameter;
  double s_curr = s_base_;
  int i = 0;
  for (; s_curr <= s; ++i) {
    double z = std::floor(s_curr / denum);
    double t0 = std::max((z - 1) * (z - 1), 1.0);
    double t1 = z * z;
    s_curr *= std::sqrt(t0 + t1);
  }
  return i - 1;
}

template <typename Integer>
absl::StatusOr<std::pair<Integer, double>>
DiscreteGaussianSampler<Integer>::SampleIIterative(SecurePrng& prng,
                                                   int i) const {
  if (i < 0) {
    return absl::InvalidArgumentError("`i` cannot be negative.");
  }

  // Base case: sample from the base distribution.
  if (i == 0) {
    RLWE_ASSIGN_OR_RETURN(Integer x, Sample(prng));
    return std::make_pair(x, s_base_);
  }

  // Iterative implementation of SampleI of https://eprint.iacr.org/2017/259,
  // which returns a sample from the discrete Gaussian over integers, with
  // center zero and standard deviation s_i, where
  // s_i = sqrt(z_i^2 + max((z_i - 1)^2, 1)) * s_{i-1}, and
  // z_i = floor(s_{i-1} / sqrt(2) * kSmoothingParameter).
  // The algorithm SampleI(i) is recursively defined as:
  // 1. x_1 <- SampleI(i - 1);
  // 2. x_2 <- SampleI(i - 1);
  // 3. return z_i * x_1 + max(1, z_i - 1) * x2.
  int num_samples = 1 << i;
  std::vector<Integer> samples;
  samples.reserve(num_samples);
  for (int j = 0; j < num_samples; ++j) {
    RLWE_ASSIGN_OR_RETURN(auto sample, Sample(prng));
    samples.push_back(std::move(sample));
  }

  double s = s_base_;  // The value of s_lvl.
  for (int lvl = 0; lvl < i; ++lvl) {
    // Compute z_lvl and s_lvl from values in the previous iteration.
    double z = std::floor(s / (std::sqrt(2) * kSmoothingParameter));
    double t = std::max((z - 1) * (z - 1), 1.0) + z * z;
    s *= std::sqrt(t);
    Integer w = static_cast<Integer>(z);
    Integer v = static_cast<Integer>(std::max(1.0, z - 1));

    // samples[k * gap] stores the samples in the previous iteration, which
    // are used to update the current iteration stored in samples[2 * k * gap].
    int gap = 1 << lvl;
    for (int j = 0; j < num_samples; j += 2 * gap) {
      Integer x1 = samples[j];
      Integer x2 = samples[j + gap];
      Integer y1 = w * x1;
      Integer y2 = v * x2;
      samples[j] = y1 + y2;
    }
  }
  // The root of the recursion tree has the desired sample.
  return std::make_pair(samples[0], s);
}

template class DiscreteGaussianSampler<Uint16>;
template class DiscreteGaussianSampler<Uint32>;
template class DiscreteGaussianSampler<Uint64>;
template class DiscreteGaussianSampler<absl::uint128>;
#ifdef ABSL_HAVE_INTRINSIC_INT128
template class DiscreteGaussianSampler<unsigned __int128>;
#endif

}  // namespace rlwe
