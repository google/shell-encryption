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

#include "shell_encryption/rns/approximate_encoder.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <memory>
#include <vector>

#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include "absl/log/check.h"
#include "absl/status/status.h"
#include "absl/types/span.h"
#include "shell_encryption/rns/rns_context.h"
#include "shell_encryption/rns/rns_modulus.h"
#include "shell_encryption/rns/rns_polynomial.h"
#include "shell_encryption/rns/testing/parameters.h"
#include "shell_encryption/rns/testing/testing_utils.h"
#include "shell_encryption/testing/status_matchers.h"
#include "shell_encryption/testing/status_testing.h"

namespace rlwe {
namespace {

using ::rlwe::testing::StatusIs;
using ::testing::HasSubstr;

template <typename ModularInt>
class ApproximateEncoderTest : public ::testing::Test {
 protected:
  using ModularIntParams = typename ModularInt::Params;
  using Integer = typename ModularInt::Int;

  void SetUpRnsParameters(
      const testing::RnsCkksParameters<ModularInt>& params) {
    // Create a RNS context suitable for CKKS.
    auto rns_context = RnsContext<ModularInt>::CreateForCkks(
        params.log_n, params.qs, params.ps);
    CHECK(rns_context.ok());
    rns_context_ = std::make_unique<const RnsContext<ModularInt>>(
        std::move(rns_context.value()));
    moduli_ = rns_context_->MainPrimeModuli();
  }

  void SetUpDefaultRnsParameters() {
    auto all_params =
        testing::GetCkksParametersForApproximateEncoding<ModularInt>();
    CHECK(!all_params.empty());
    SetUpRnsParameters(all_params[0]);
  }

  // Computes the additive precision loss between `actual` and `orig`.
  // We slack a bit and use distance on the real and the imaginary parts as the
  // modulus of a complex number.
  static double GetPrecision(const std::vector<std::complex<double>>& actual,
                             const std::vector<std::complex<double>>& orig) {
    CHECK_EQ(actual.size(), orig.size());
    double precision = 0;
    for (int i = 0; i < actual.size(); ++i) {
      std::complex<double> diff = orig[i] - actual[i];
      if (diff.real() > precision) {
        precision = diff.real();
      }
      if (diff.imag() > precision) {
        precision = diff.imag();
      }
    }
    return precision;
  }

  // Computes the minimal precision expected from Decode(Encode(messages)),
  // where precision is measured as l_inf norm of messages - decoded.
  // For random messages, we expect the encoding error (due to rounding the
  // scaled up messages) is about sqrt(N) * scaling_factor, and we add some
  // buffer to account for minor variations and additional small errors.
  static double GetPrecisionBound(
      const testing::RnsCkksParameters<ModularInt>& params) {
    return 1 / params.scaling_factor * ldexp(1, params.log_n / 2 + 2);
  }

  std::unique_ptr<const RnsContext<ModularInt>> rns_context_;
  std::vector<const PrimeModulus<ModularInt>*> moduli_;
};

TYPED_TEST_SUITE(ApproximateEncoderTest,
                 testing::ModularIntTypesForApproximateEncoding);

TYPED_TEST(ApproximateEncoderTest, CreateFailsIfNullContext) {
  EXPECT_THAT(ApproximateEncoder<TypeParam>::Create(/*context=*/nullptr,
                                                    /*scaling_factor=*/2.0),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`context` must not be null")));
}

TYPED_TEST(ApproximateEncoderTest, CreateFailsIfScalingFactorIsTooSmall) {
  this->SetUpDefaultRnsParameters();
  EXPECT_THAT(ApproximateEncoder<TypeParam>::Create(this->rns_context_.get(),
                                                    /*scaling_factor=*/-1.0),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`scaling_factor` must be at least 1")));
  EXPECT_THAT(ApproximateEncoder<TypeParam>::Create(this->rns_context_.get(),
                                                    /*scaling_factor=*/0),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`scaling_factor` must be at least 1")));
}

TYPED_TEST(ApproximateEncoderTest, CreateSucceeds) {
  this->SetUpDefaultRnsParameters();
  ASSERT_OK(ApproximateEncoder<TypeParam>::Create(this->rns_context_.get(),
                                                  /*scaling_factor=*/2.0));
}

TYPED_TEST(ApproximateEncoderTest, EncodeCkksFailsIfTooManyValues) {
  this->SetUpDefaultRnsParameters();
  ASSERT_OK_AND_ASSIGN(auto encoder, ApproximateEncoder<TypeParam>::Create(
                                         this->rns_context_.get(),
                                         /*scaling_factor=*/2.0));
  int num_slots = 1 << (this->rns_context_->LogN() - 1);
  std::vector<std::complex<double>> too_long_input(num_slots + 1, {0, 0});
  EXPECT_THAT(encoder.EncodeCkks(too_long_input, this->moduli_),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`values` cannot have more than")));
}

TYPED_TEST(ApproximateEncoderTest, EncodeCkksFailsIfEmptyModuli) {
  this->SetUpDefaultRnsParameters();
  ASSERT_OK_AND_ASSIGN(auto encoder, ApproximateEncoder<TypeParam>::Create(
                                         this->rns_context_.get(),
                                         /*scaling_factor=*/2.0));
  std::vector<std::complex<double>> values = {{0, 0}};
  EXPECT_THAT(encoder.EncodeCkks(values, /*moduli=*/{}),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`moduli` cannot be empty")));
}

TYPED_TEST(ApproximateEncoderTest, DecodeCkksFailsIfMismatchModuli) {
  // Find the first set of parameters with more than one prime modulus.
  for (auto params :
       testing::GetCkksParametersForApproximateEncoding<TypeParam>()) {
    if (params.qs.size() >= 2) {
      this->SetUpRnsParameters(params);
      break;
    }
  }
  ASSERT_TRUE(this->rns_context_ != nullptr);

  ASSERT_OK_AND_ASSIGN(auto encoder, ApproximateEncoder<TypeParam>::Create(
                                         this->rns_context_.get(),
                                         /*scaling_factor=*/2.0));

  int log_n = this->rns_context_->LogN();
  auto small_moduli = absl::MakeSpan(this->moduli_).subspan(0, 1);
  ASSERT_OK_AND_ASSIGN(
      auto poly_with_small_modulus,
      RnsPolynomial<TypeParam>::CreateZero(log_n, small_moduli));

  EXPECT_THAT(encoder.DecodeCkks(poly_with_small_modulus, this->moduli_),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`moduli` must contain")));
}

TYPED_TEST(ApproximateEncoderTest, EncodeCkksDecodes) {
  for (auto const& params :
       testing::GetCkksParametersForApproximateEncoding<TypeParam>()) {
    this->SetUpRnsParameters(params);
    ASSERT_OK_AND_ASSIGN(auto encoder,
                         ApproximateEncoder<TypeParam>::Create(
                             this->rns_context_.get(), params.scaling_factor));

    int num_slots = 1 << (this->rns_context_->LogN() - 1);
    std::vector<std::complex<double>> messages = testing::SampleComplexValues(
        num_slots, /*max_value=*/1, /*precision=*/1 << 10);

    // Encode the messages.
    ASSERT_OK_AND_ASSIGN(auto poly,
                         encoder.EncodeCkks(messages, this->moduli_));

    // Now decode the RNS polynomial.
    ASSERT_OK_AND_ASSIGN(auto decoded, encoder.DecodeCkks(poly, this->moduli_));
    double precision_bound = this->GetPrecisionBound(params);
    double precision_actual = this->GetPrecision(decoded, messages);
    EXPECT_LE(precision_actual, precision_bound);
  }
}

TYPED_TEST(ApproximateEncoderTest, EncodeCkksIsAdditiveHomomorphic) {
  for (auto const& params :
       testing::GetCkksParametersForApproximateEncoding<TypeParam>()) {
    this->SetUpRnsParameters(params);
    ASSERT_OK_AND_ASSIGN(auto encoder,
                         ApproximateEncoder<TypeParam>::Create(
                             this->rns_context_.get(), params.scaling_factor));

    // Sample two random vectors, and add the polynomials encoding them.
    int num_slots = 1 << (this->rns_context_->LogN() - 1);
    std::vector<std::complex<double>> messages0 = testing::SampleComplexValues(
        num_slots, /*max_value=*/1, /*precision=*/1 << 10);
    std::vector<std::complex<double>> messages1 = testing::SampleComplexValues(
        num_slots, /*max_value=*/1, /*precision=*/1 << 10);
    ASSERT_OK_AND_ASSIGN(auto poly0,
                         encoder.EncodeCkks(messages0, this->moduli_));
    ASSERT_OK_AND_ASSIGN(auto poly1,
                         encoder.EncodeCkks(messages1, this->moduli_));

    ASSERT_OK(poly0.AddInPlace(poly1, this->moduli_));

    // Now decode the polynomial which encodes the sum of two vectors.
    ASSERT_OK_AND_ASSIGN(auto decoded,
                         encoder.DecodeCkks(poly0, this->moduli_));
    std::vector<std::complex<double>> expected;
    for (int i = 0; i < num_slots; ++i) {
      expected.push_back(messages0[i] + messages1[i]);
    }
    double precision_bound = this->GetPrecisionBound(params);
    double precision = this->GetPrecision(decoded, expected);
    EXPECT_LE(precision, precision_bound);
  }
}

TYPED_TEST(ApproximateEncoderTest, EncodeCkksIsMultiplicativeHomomorphic) {
  for (auto const& params :
       testing::GetCkksParametersForApproximateEncoding<TypeParam>()) {
    this->SetUpRnsParameters(params);
    ASSERT_OK_AND_ASSIGN(auto encoder,
                         ApproximateEncoder<TypeParam>::Create(
                             this->rns_context_.get(), params.scaling_factor));

    // Sample random complex vectors and encode them.
    int num_slots = 1 << (this->rns_context_->LogN() - 1);
    std::vector<std::complex<double>> messages0 = testing::SampleComplexValues(
        num_slots, /*max_value=*/1, /*precision=*/1 << 10);
    std::vector<std::complex<double>> messages1 = testing::SampleComplexValues(
        num_slots, /*max_value=*/1, /*precision=*/1 << 10);
    ASSERT_OK_AND_ASSIGN(auto poly0,
                         encoder.EncodeCkks(messages0, this->moduli_));
    ASSERT_OK_AND_ASSIGN(auto poly1,
                         encoder.EncodeCkks(messages1, this->moduli_));

    // Multiply the two plaintext polynomials.
    ASSERT_OK_AND_ASSIGN(auto product, poly0.Mul(poly1, this->moduli_));

    // Now decode the product polynomial.
    ASSERT_OK_AND_ASSIGN(auto decoded,
                         encoder.DecodeCkks(product, this->moduli_));
    for (auto& d : decoded) {
      d /= params.scaling_factor;  // rescale the product values
    }
    std::vector<std::complex<double>> expected;
    for (int i = 0; i < num_slots; ++i) {
      expected.push_back(messages0[i] * messages1[i]);
    }
    double precision_bound = this->GetPrecisionBound(params);
    double precision = this->GetPrecision(decoded, expected);
    EXPECT_LE(precision, precision_bound);
  }
}

TYPED_TEST(ApproximateEncoderTest, SlotRotation) {
  for (auto const& params :
       testing::GetCkksParametersForApproximateEncoding<TypeParam>()) {
    this->SetUpRnsParameters(params);
    ASSERT_OK_AND_ASSIGN(auto encoder,
                         ApproximateEncoder<TypeParam>::Create(
                             this->rns_context_.get(), params.scaling_factor));
    int num_slots = 1 << (this->rns_context_->LogN() - 1);
    std::vector<std::complex<double>> messages = testing::SampleComplexValues(
        num_slots, /*max_value=*/1, /*precision=*/1 << 10);
    ASSERT_OK_AND_ASSIGN(auto poly,
                         encoder.EncodeCkks(messages, this->moduli_));

    constexpr int k_power = 5;  // rotate by one slot
    ASSERT_OK_AND_ASSIGN(auto rotated_poly,
                         poly.Substitute(k_power, this->moduli_));

    // Now decode the polynomial that encodes the rotated vector.
    ASSERT_OK_AND_ASSIGN(auto decoded,
                         encoder.DecodeCkks(rotated_poly, this->moduli_));
    std::vector<std::complex<double>> expected(num_slots, {0, 0});
    for (int i = 0; i < num_slots; ++i) {
      expected[i] = messages[(i + 1) % num_slots];
    }
    double precision_bound = this->GetPrecisionBound(params);
    double precision = this->GetPrecision(decoded, expected);
    EXPECT_LE(precision, precision_bound);
  }
}

TYPED_TEST(ApproximateEncoderTest, SlotConjugate) {
  for (auto const& params :
       testing::GetCkksParametersForApproximateEncoding<TypeParam>()) {
    this->SetUpRnsParameters(params);
    ASSERT_OK_AND_ASSIGN(auto encoder,
                         ApproximateEncoder<TypeParam>::Create(
                             this->rns_context_.get(), params.scaling_factor));
    int n = 1 << this->rns_context_->LogN();
    int num_slots = n >> 1;
    std::vector<std::complex<double>> messages = testing::SampleComplexValues(
        num_slots, /*max_value=*/1, /*precision=*/1 << 10);
    ASSERT_OK_AND_ASSIGN(auto poly,
                         encoder.EncodeCkks(messages, this->moduli_));

    // The polynomial poly(1/X) encodes the conjugate of messages.
    int power = 2 * n - 1;
    ASSERT_OK_AND_ASSIGN(auto poly_conjugate,
                         poly.Substitute(power, this->moduli_));
    ASSERT_OK_AND_ASSIGN(auto decoded,
                         encoder.DecodeCkks(poly_conjugate, this->moduli_));
    std::vector<std::complex<double>> expected(num_slots, {0, 0});
    for (int i = 0; i < num_slots; ++i) {
      expected[i].real(messages[i].real());
      expected[i].imag(-messages[i].imag());
    }
    double precision_bound = this->GetPrecisionBound(params);
    double precision = this->GetPrecision(decoded, expected);
    EXPECT_LE(precision, precision_bound);
  }
}

}  // namespace
}  // namespace rlwe
