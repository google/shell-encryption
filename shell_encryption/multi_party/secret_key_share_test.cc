// Copyright 2025 Google LLC
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

#include "shell_encryption/multi_party/secret_key_share.h"

#include <memory>
#include <vector>

#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include "absl/log/check.h"
#include "absl/status/status.h"
#include "shell_encryption/multi_party/testing/parameters.h"
#include "shell_encryption/rns/rns_context.h"
#include "shell_encryption/rns/rns_modulus.h"
#include "shell_encryption/rns/rns_polynomial.h"
#include "shell_encryption/sampler/discrete_gaussian.h"
#include "shell_encryption/testing/status_matchers.h"
#include "shell_encryption/testing/status_testing.h"
#include "shell_encryption/testing/testing_prng.h"

namespace rlwe {
namespace multi_party {
namespace {

using Prng = ::rlwe::testing::TestingPrng;
using ::rlwe::testing::StatusIs;
using ::testing::HasSubstr;
using testing::ModularInt32;

constexpr double kSBase = 10;
constexpr double kSFlood = 10;

template <typename ModularInt>
class SecretKeyShareTest : public ::testing::Test {
 protected:
  void SetUp() override {
    testing::MpaheParameters<ModularInt> params =
        testing::GetMultiPartyParameters<ModularInt>();

    auto rns_context = RnsContext<ModularInt>::CreateForBfv(
        params.log_n, params.qs, params.ps, params.t);
    CHECK(rns_context.ok());
    rns_context_ = std::make_unique<const RnsContext<ModularInt>>(
        std::move(rns_context.value()));
    moduli_ = rns_context_->MainPrimeModuli();

    prng_ = std::make_unique<Prng>(0);
  }

  std::unique_ptr<const RnsContext<ModularInt>> rns_context_;
  std::vector<const PrimeModulus<ModularInt>*> moduli_;
  std::unique_ptr<Prng> prng_;
};

TYPED_TEST_SUITE(SecretKeyShareTest, testing::ModularIntTypesForMultiParty);

class SecretKeyShareNegativeTest : public SecretKeyShareTest<ModularInt32> {};

TEST_F(SecretKeyShareNegativeTest, SampleFailsIfRnsContextIsNull) {
  EXPECT_THAT(SecretKeyShare<ModularInt32>::Sample(/*rns_context=*/nullptr,
                                                   this->prng_.get()),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`rns_context` must not be null")));
}

TEST_F(SecretKeyShareNegativeTest, SampleFailsIfPrngIsNull) {
  EXPECT_THAT(SecretKeyShare<ModularInt32>::Sample(this->rns_context_.get(),
                                                   /*prng=*/nullptr),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`prng` must not be null")));
}

TEST_F(SecretKeyShareNegativeTest, PartialDecryptFailsIfGaussianSamplerIsNull) {
  ASSERT_OK_AND_ASSIGN(auto secret_key_share,
                       SecretKeyShare<ModularInt32>::Sample(
                           this->rns_context_.get(), this->prng_.get()));
  ASSERT_OK_AND_ASSIGN(auto ct, RnsPolynomial<ModularInt32>::CreateZero(
                                    this->rns_context_->LogN(), this->moduli_));
  EXPECT_THAT(secret_key_share.PartialDecrypt(
                  ct, kSFlood, /*dg_sampler=*/nullptr, this->prng_.get()),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`dg_sampler` must not be null")));
}

TEST_F(SecretKeyShareNegativeTest, PartialDecryptFailsIfPrngIsNull) {
  ASSERT_OK_AND_ASSIGN(auto secret_key_share,
                       SecretKeyShare<ModularInt32>::Sample(
                           this->rns_context_.get(), this->prng_.get()));
  ASSERT_OK_AND_ASSIGN(auto ct, RnsPolynomial<ModularInt32>::CreateZero(
                                    this->rns_context_->LogN(), this->moduli_));
  ASSERT_OK_AND_ASSIGN(
      auto dg_sampler,
      DiscreteGaussianSampler<ModularInt32::Int>::Create(kSBase));
  EXPECT_THAT(secret_key_share.PartialDecrypt(ct, kSFlood, dg_sampler.get(),
                                              /*prng=*/nullptr),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`prng` must not be null")));
}

TYPED_TEST(SecretKeyShareTest, SampleSuccess) {
  ASSERT_OK_AND_ASSIGN(auto secret_key_share,
                       SecretKeyShare<TypeParam>::Sample(
                           this->rns_context_.get(), this->prng_.get()));
  EXPECT_EQ(secret_key_share.LogN(), this->rns_context_->LogN());
  EXPECT_EQ(secret_key_share.NumCoeffs(), 1 << this->rns_context_->LogN());
  EXPECT_EQ(secret_key_share.NumModuli(),
            this->rns_context_->NumMainPrimeModuli());

  const RnsPolynomial<TypeParam>& s = secret_key_share.Key();
  EXPECT_EQ(s.NumCoeffs(), 1 << this->rns_context_->LogN());
  EXPECT_EQ(s.NumModuli(), this->rns_context_->NumMainPrimeModuli());
}

}  // namespace
}  // namespace multi_party
}  // namespace rlwe
