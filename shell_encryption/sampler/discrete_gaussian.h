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

#ifndef RLWE_SAMPLER_DISCRETE_GAUSSIAN_H_
#define RLWE_SAMPLER_DISCRETE_GAUSSIAN_H_

#include <limits>
#include <memory>
#include <tuple>
#include <utility>
#include <vector>

#include "absl/status/statusor.h"
#include "shell_encryption/integral_types.h"
#include "shell_encryption/prng/prng.h"
#include "shell_encryption/status_macros.h"

namespace rlwe {

// This class implements a sampler for discrete Gaussian distributions on the
// integers, with center 0 and arbitrary standard deviation `s`, as specified in
// ``Gaussian Sampling over the Integers: Efficient, Generic, Constant-Time'' by
// Daniele Micciancio and Michael Walter, https://eprint.iacr.org/2017/259.
// Specifically, it samples from the distribution DG_s over the integers where
// the probability of any integer x is proportional to exp(-x^2 / (2 * s^2)).
//
// The sampler is built on a base discrete Gaussian sampler with a fixed, small,
// Gaussian parameter `s_base`, which is implemented using the inverse CDT
// method. The main sampler recursively combines base samples, resulting in the
// desired Gaussian parameter following a convolution theorem.
//
// Note that the base sampler's standard deviation `s_base` must be at least
// sqrt(2) * \eta, where \eta is the smoothing parameter of the integers whose
// concrete value is defined in `kSmoothParameter`. One can instantiate this
// sampler using a larger `s_base`, which can reduce the time to sample from
// a target discrete Gaussian distribution but at the cost of increased memory
// footprint to store the CDT of the base distribution.
//
// The template parameter `Integer` is assumed to be an unsigned integer type,
// with enough width to represent the range [+/- kTailBoundMultiplier * s].
template <typename Integer>
class DiscreteGaussianSampler {
 public:
  // The smoothing parameter of the integers, with an error parameter epsilon <=
  // 2^{-128}. The value is computed w/ the bound sqrt(ln(2 + 2/epsilon) / PI).
  static constexpr double kSmoothingParameter = 5.3349782;  // for eps = 2^-128

  // The tail bound cut-off multiplier such that the probability of a sample of
  // DG_s being outside of [+/- `kTailBoundMultiplier` * s] is negligible.
  static constexpr int kTailBoundMultiplier = 8;

  // The threshold such that Integer values larger than it represents negative
  // numbers.
  static constexpr Integer kNegativeThreshold =
      std::numeric_limits<Integer>::max() >> 1;

  // Factory function that returns a discrete Gaussian sampler with the given
  // standard deviation `s_base` for the base sampler.
  static absl::StatusOr<std::unique_ptr<DiscreteGaussianSampler>> Create(
      double s_base);

  // Returns a sample from the discrete Gaussian distribution DG_s with center 0
  // and standard deviation `s`, where the randomness is drawn using the secure
  // pseudorandom generator `prng`. The return value represents negative number
  // if it is larger than `kNegativeThreshold`.
  // Note: `s` must be at least `s_base`.
  absl::StatusOr<Integer> Sample(double s, SecurePrng& prng) const {
    RLWE_ASSIGN_OR_RETURN(int num_iterations, NumIterations(s));
    return SampleWithIterations(s, num_iterations, prng);
  }

  // Returns a sample from the discrete Gaussian distribution DG_s with center 0
  // and standard deviation `s`, with the given number of iterations in the
  // recursive convolution computation.
  //
  // This function should be used to generate a large number of samples with
  // the same `s`, such that the parameter `num_iterations` can be precomputed
  // by `NumIterations(s)` only once.
  // As in `Sample`, the randomness is supplied from the secure pseudorandom
  // generator `prng`, and the return value represents negative number if it is
  // larger than `kNegativeThreshold`.
  // Note: `s` must be at least `s_base`.
  absl::StatusOr<Integer> SampleWithIterations(double s, int num_iterations,
                                               SecurePrng& prng) const;

  // Returns the number of iterations needed to combine samples from the base
  // distribution to obtain a sample with the standard deviation `s`.
  absl::StatusOr<int> NumIterations(double s) const;

 private:
  explicit DiscreteGaussianSampler(double s_base, std::vector<Uint64> cdt)
      : s_base_(s_base), cdt_(std::move(cdt)) {}

  // Returns a sample from the base distribution. The return value represents a
  // negative number if it is larger than `kNegativeThreshold`.
  absl::StatusOr<Integer> SampleBase(SecurePrng& prng) const;

  // Returns a sample from running i iterations of the arbitrary deviation
  // discrete Gaussian sampling algorithm SampleI.
  absl::StatusOr<std::pair<Integer, double>> SampleIIterative(SecurePrng& prng,
                                                              int i) const;

  double s_base_;

  std::vector<Uint64> cdt_;
};

}  // namespace rlwe

#endif  // RLWE_SAMPLER_DISCRETE_GAUSSIAN_H_
