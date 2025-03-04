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

#include "shell_encryption/multi_party/public_parameter.h"

#include <memory>
#include <utility>
#include <vector>

#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include "absl/log/check.h"
#include "absl/status/status.h"
#include "shell_encryption/multi_party/testing/parameters.h"
#include "shell_encryption/rns/rns_context.h"
#include "shell_encryption/rns/rns_modulus.h"
#include "shell_encryption/rns/rns_polynomial.h"
#include "shell_encryption/testing/status_matchers.h"
#include "shell_encryption/testing/status_testing.h"

namespace rlwe {
namespace multi_party {
namespace {

using ::rlwe::testing::StatusIs;
using ::testing::HasSubstr;
using testing::ModularInt32;

constexpr int kVariance = 8;
constexpr PrngType kPrngType = PRNG_TYPE_HKDF;

template <typename ModularInt>
class PublicParameterTest : public ::testing::Test {
 protected:
  void SetUp() override {
    testing::MpaheParameters<ModularInt> params =
        testing::GetMultiPartyParameters<ModularInt>();

    // Create the RNS context.
    auto rns_context = RnsContext<ModularInt>::CreateForBfv(
        params.log_n, params.qs, params.ps, params.t);
    CHECK(rns_context.ok());
    rns_context_ = std::make_unique<const RnsContext<ModularInt>>(
        std::move(rns_context.value()));
    moduli_ = rns_context_->MainPrimeModuli();
  }

  std::unique_ptr<const RnsContext<ModularInt>> rns_context_;
  std::vector<const PrimeModulus<ModularInt>*> moduli_;
};

TYPED_TEST_SUITE(PublicParameterTest, testing::ModularIntTypesForMultiParty);

class PublicParameterNegativeTest : public PublicParameterTest<ModularInt32> {};

TEST_F(PublicParameterNegativeTest, CreateFailsIfRnsContextIsNull) {
  EXPECT_THAT(PublicParameter<ModularInt32>::Create(/*rns_context=*/nullptr,
                                                    kVariance, kPrngType),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`rns_context` must not be null")));
}

TEST_F(PublicParameterNegativeTest, CreateFailsIfVarianceIsNotPositive) {
  EXPECT_THAT(PublicParameter<ModularInt32>::Create(this->rns_context_.get(),
                                                    /*variance=*/-1, kPrngType),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`error_variance` must be positive")));
  EXPECT_THAT(PublicParameter<ModularInt32>::Create(this->rns_context_.get(),
                                                    /*variance=*/0, kPrngType),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`error_variance` must be positive")));
}

TEST_F(PublicParameterNegativeTest, CreateFailsIfInvalidPrngType) {
  EXPECT_THAT(PublicParameter<ModularInt32>::Create(
                  this->rns_context_.get(), kVariance, PRNG_TYPE_INVALID),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`prng_type` not specified correctly")));
}

TEST_F(PublicParameterNegativeTest, DeserializeFailsIfRnsContextIsNull) {
  SerializedPublicParameter serialized;
  EXPECT_THAT(
      PublicParameter<ModularInt32>::Deserialize(serialized,
                                                 /*rns_context=*/nullptr),
      StatusIs(absl::StatusCode::kInvalidArgument,
               HasSubstr("`rns_context` must not be null")));
}

TEST_F(PublicParameterNegativeTest, DeserializeFailsIfVarianceIsNotPositive) {
  SerializedPublicParameter serialized;
  serialized.set_error_variance(0);
  EXPECT_THAT(PublicParameter<ModularInt32>::Deserialize(
                  serialized, this->rns_context_.get()),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`error_variance` must be positive")));
}

TEST_F(PublicParameterNegativeTest, DeserializeFailsIfInvalidPrngType) {
  SerializedPublicParameter serialized;
  serialized.set_error_variance(kVariance);
  // Default prng_type is PRNG_TYPE_INVALID
  ASSERT_EQ(serialized.prng_type(), PRNG_TYPE_INVALID);
  EXPECT_THAT(PublicParameter<ModularInt32>::Deserialize(
                  serialized, this->rns_context_.get()),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("Invalid `prng_type`")));
}

TYPED_TEST(PublicParameterTest, Create) {
  ASSERT_OK_AND_ASSIGN(auto public_parameter,
                       PublicParameter<TypeParam>::Create(
                           this->rns_context_.get(), kVariance, kPrngType));
  EXPECT_EQ(public_parameter->ErrorVariance(), kVariance);
  EXPECT_EQ(public_parameter->LogN(), this->rns_context_->LogN());
  EXPECT_EQ(public_parameter->NumCoeffs(), 1 << this->rns_context_->LogN());
  EXPECT_EQ(public_parameter->NumModuli(),
            this->rns_context_->NumMainPrimeModuli());

  const RnsPolynomial<TypeParam>& a = public_parameter->PublicKeyComponentA();
  EXPECT_EQ(a.NumCoeffs(), 1 << this->rns_context_->LogN());
  EXPECT_EQ(a.NumModuli(), this->rns_context_->NumMainPrimeModuli());
}

TYPED_TEST(PublicParameterTest, SerializeDeserialize) {
  ASSERT_OK_AND_ASSIGN(auto public_parameter,
                       PublicParameter<TypeParam>::Create(
                           this->rns_context_.get(), kVariance, kPrngType));
  ASSERT_OK_AND_ASSIGN(auto serialized, public_parameter->Serialize());
  ASSERT_OK_AND_ASSIGN(auto deserialized,
                       PublicParameter<TypeParam>::Deserialize(
                           serialized, this->rns_context_.get()));
  EXPECT_EQ(deserialized->LogN(), public_parameter->LogN());
  EXPECT_EQ(deserialized->NumCoeffs(), public_parameter->NumCoeffs());
  EXPECT_EQ(deserialized->NumModuli(), public_parameter->NumModuli());
  EXPECT_EQ(deserialized->ErrorVariance(), public_parameter->ErrorVariance());
  EXPECT_EQ(deserialized->PublicKeyComponentA(),
            public_parameter->PublicKeyComponentA());
}

}  // namespace
}  // namespace multi_party
}  // namespace rlwe
