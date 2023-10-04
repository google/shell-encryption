// Copyright 2023 Google LLC
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

#include "shell_encryption/rns/rns_error_params.h"

#include <cmath>
#include <string>
#include <vector>

#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include "absl/log/check.h"
#include "absl/random/random.h"
#include "absl/status/status.h"
#include "absl/types/span.h"
#include "shell_encryption/integral_types.h"
#include "shell_encryption/prng/single_thread_hkdf_prng.h"
#include "shell_encryption/rns/coefficient_encoder.h"
#include "shell_encryption/rns/crt_interpolation.h"
#include "shell_encryption/rns/error_distribution.h"
#include "shell_encryption/rns/rns_context.h"
#include "shell_encryption/rns/rns_modulus.h"
#include "shell_encryption/rns/rns_polynomial.h"
#include "shell_encryption/rns/testing/parameters.h"
#include "shell_encryption/rns/testing/testing_utils.h"
#include "shell_encryption/testing/parameters.h"
#include "shell_encryption/testing/status_matchers.h"
#include "shell_encryption/testing/status_testing.h"

namespace rlwe {
namespace {

using ::rlwe::testing::StatusIs;
using ::testing::HasSubstr;
using Prng = SingleThreadHkdfPrng;

// Number of samples used to compute the actual variance.
constexpr int kSamples = 50;
constexpr int kVariance = 8;

// This test checks the heuristics used to compute the error bounds are valid.
template <typename ModularInt>
class RnsErrorParamsTest : public ::testing::Test {
  using Integer = typename ModularInt::Int;
  using BigInteger = typename testing::CompositeModulus<Integer>::value_type;

 public:
  void SetUp() override {
    testing::RnsParameters<ModularInt> rns_params =
        testing::GetRnsParameters<ModularInt>();
    auto rns_context = RnsContext<ModularInt>::Create(
        rns_params.log_n, rns_params.qs, rns_params.ps, rns_params.t);
    CHECK(rns_context.ok());
    rns_context_ = std::make_unique<const RnsContext<ModularInt>>(
        std::move(rns_context.value()));
    main_moduli_ = rns_context_->MainPrimeModuli();
    aux_moduli_ = rns_context_->AuxPrimeModuli();

    // Precompute the constants for CRT interpolation over the main modulus.
    main_modulus_hats_ =
        RnsModulusComplements<ModularInt, BigInteger>(main_moduli_);

    // Compute the main modulus.
    main_modulus_ = 1;
    for (auto qi : rns_params.qs) {
      main_modulus_ *= static_cast<BigInteger>(qi);
    }

    // Cache the bit size of the plaintext modulus.
    log_t_ = log2(static_cast<double>(rns_context_->PlaintextModulus()));
  }

  // Computes the l-infinity norm of the coefficients of a RNS polynomial.
  // Note: We assume the polynomial is in coefficient form, using balanced
  // representation mod Q.
  BigInteger LInfinityNorm(
      const RnsPolynomial<ModularInt>& poly,
      absl::Span<const PrimeModulus<ModularInt>* const> moduli) const {
    // Compute the CRT interpolation of RNS coefficients.
    std::vector<ModularInt> modulus_hat_invs =
        rns_context_->MainPrimeModulusCrtFactors(moduli.size() - 1).value();
    std::vector<BigInteger> coeffs =
        CrtInterpolation<ModularInt, BigInteger>(
            poly.Coeffs(), moduli, main_modulus_hats_, modulus_hat_invs)
            .value();

    BigInteger main_modulus_half = main_modulus_ >> 1;
    BigInteger max_coeff_abs{0};
    for (const auto& coeff : coeffs) {
      BigInteger coeff_abs =
          coeff > main_modulus_half ? (main_modulus_ - coeff) : coeff;
      if (coeff_abs > max_coeff_abs) {
        max_coeff_abs = coeff_abs;
      }
    }
    return max_coeff_abs;
  }

  std::unique_ptr<const RnsContext<ModularInt>> rns_context_;
  std::vector<const PrimeModulus<ModularInt>*> main_moduli_;
  std::vector<const PrimeModulus<ModularInt>*> aux_moduli_;
  std::vector<BigInteger> main_modulus_hats_;
  BigInteger main_modulus_;
  int log_t_;
};
TYPED_TEST_SUITE(RnsErrorParamsTest, rlwe::testing::ModularIntTypes);

TYPED_TEST(RnsErrorParamsTest, CreateFailsIfLogNIsNotPositive) {
  EXPECT_THAT(RnsErrorParams<TypeParam>::Create(
                  /*log_n=*/0, this->main_moduli_, this->aux_moduli_,
                  this->log_t_, sqrt(kVariance)),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`log_n` must be positive")));
  EXPECT_THAT(RnsErrorParams<TypeParam>::Create(
                  /*log_n=*/-1, this->main_moduli_, this->aux_moduli_,
                  this->log_t_, sqrt(kVariance)),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`log_n` must be positive")));
}

TYPED_TEST(RnsErrorParamsTest, CreateFailsIfMainModuliIsEmpty) {
  int log_n = this->rns_context_->LogN();
  EXPECT_THAT(RnsErrorParams<TypeParam>::Create(log_n, /*main_moduli=*/{},
                                                this->aux_moduli_, this->log_t_,
                                                sqrt(kVariance)),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`main_moduli` must not be empty")));
}

TYPED_TEST(RnsErrorParamsTest, CreateFailsIfLogPlaintextModulusIsNotPositive) {
  int log_n = this->rns_context_->LogN();
  EXPECT_THAT(RnsErrorParams<TypeParam>::Create(
                  log_n, this->main_moduli_, this->aux_moduli_,
                  /*log_plainext_modulus=*/0, sqrt(kVariance)),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`log_plaintext_modulus` must be positive")));
  EXPECT_THAT(RnsErrorParams<TypeParam>::Create(
                  log_n, this->main_moduli_, this->aux_moduli_,
                  /*log_plainext_modulus=*/-1, sqrt(kVariance)),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`log_plaintext_modulus` must be positive")));
}

TYPED_TEST(RnsErrorParamsTest, CreateFailsIfSigmaIsNotPositive) {
  int log_n = this->rns_context_->LogN();
  EXPECT_THAT(RnsErrorParams<TypeParam>::Create(log_n, this->main_moduli_,
                                                this->aux_moduli_, this->log_t_,
                                                /*sigma=*/0),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`sigma` must be positive")));
  EXPECT_THAT(RnsErrorParams<TypeParam>::Create(log_n, this->main_moduli_,
                                                this->aux_moduli_, this->log_t_,
                                                /*sigma=*/-1),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`sigma` must be positive")));
}

TYPED_TEST(RnsErrorParamsTest, PlaintextBoundIsValid) {
  using Integer = typename TypeParam::Int;
  int log_n = this->rns_context_->LogN();
  ASSERT_OK_AND_ASSIGN(auto error_params,
                       RnsErrorParams<TypeParam>::Create(
                           log_n, this->main_moduli_, this->aux_moduli_,
                           this->log_t_, sqrt(kVariance)));

  // Create a coefficient encoder to convert messages to plaintext polynomials.
  ASSERT_OK_AND_ASSIGN(auto encoder, CoefficientEncoder<TypeParam>::Create(
                                         this->rns_context_.get()));

  // Sample uniform random plaintext polynomials and expect their l_inf norm is
  // bounded by b_plaintext.
  Integer t = this->rns_context_->PlaintextModulus();
  for (int i = 0; i < kSamples; ++i) {
    ASSERT_OK_AND_ASSIGN(std::string prng_seed,
                         rlwe::SingleThreadHkdfPrng::GenerateSeed());
    ASSERT_OK_AND_ASSIGN(auto prng,
                         rlwe::SingleThreadHkdfPrng::Create(prng_seed));

    std::vector<Integer> messages = testing::SampleMessages(1 << log_n, t);
    ASSERT_OK_AND_ASSIGN(RnsPolynomial<TypeParam> a,
                         encoder.EncodeBgv(messages, this->main_moduli_));
    // Polynomial should be in coefficient form to measure its norm.
    ASSERT_FALSE(a.IsNttForm());

    double norm =
        static_cast<double>(this->LInfinityNorm(a, this->main_moduli_));
    EXPECT_LT(norm, error_params.B_plaintext());
  }
}

TYPED_TEST(RnsErrorParamsTest, SecretKeyEncryptionBoundIsValid) {
  int log_n = this->rns_context_->LogN();
  ASSERT_OK_AND_ASSIGN(auto error_params,
                       RnsErrorParams<TypeParam>::Create(
                           log_n, this->main_moduli_, this->aux_moduli_,
                           this->log_t_, sqrt(kVariance)));

  // Sample fresh and scaled error polynomials and expect their l_inf norm is
  // bounded by b_secretkey_encryption.
  for (int i = 0; i < kSamples; ++i) {
    ASSERT_OK_AND_ASSIGN(std::string prng_seed,
                         rlwe::SingleThreadHkdfPrng::GenerateSeed());
    ASSERT_OK_AND_ASSIGN(auto prng,
                         rlwe::SingleThreadHkdfPrng::Create(prng_seed));

    ASSERT_OK_AND_ASSIGN(
        RnsPolynomial<TypeParam> e,
        (SampleError<TypeParam>(log_n, kVariance, this->main_moduli_,
                                prng.get())));
    ASSERT_OK(e.MulInPlace(this->rns_context_->PlaintextModulus(),
                           this->main_moduli_));
    ASSERT_OK(e.ConvertToCoeffForm(this->main_moduli_));

    double norm =
        static_cast<double>(this->LInfinityNorm(e, this->main_moduli_));
    EXPECT_LT(norm, error_params.B_secretkey_encryption());
  }
}

}  // namespace
}  // namespace rlwe
