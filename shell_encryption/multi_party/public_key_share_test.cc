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

#include "shell_encryption/multi_party/public_key_share.h"

#include <memory>
#include <vector>

#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include "absl/log/check.h"
#include "absl/status/status.h"
#include "shell_encryption/multi_party/public_parameter.h"
#include "shell_encryption/multi_party/secret_key_share.h"
#include "shell_encryption/multi_party/testing/parameters.h"
#include "shell_encryption/prng/single_thread_hkdf_prng.h"
#include "shell_encryption/rns/rns_context.h"
#include "shell_encryption/rns/rns_modulus.h"
#include "shell_encryption/rns/rns_polynomial.h"
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

constexpr int kVariance = 8;
constexpr PrngType kPrngType = PRNG_TYPE_HKDF;

template <typename ModularInt>
class PublicKeyShareTest : public ::testing::Test {
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

    // Create a PRNG for testing only.
    prng_ = std::make_unique<Prng>(0);
  }

  void SetupPublicParameter() {
    auto public_parameter = PublicParameter<ModularInt>::Create(
        rns_context_.get(), kVariance, kPrngType);
    CHECK(public_parameter.ok());
    public_parameter_ = std::move(public_parameter.value());
  }

  void SetupSecretKeyShare() {
    auto secret_key_share =
        SecretKeyShare<ModularInt>::Sample(rns_context_.get(), prng_.get());
    CHECK(secret_key_share.ok());
    secret_key_share_ = std::make_unique<const SecretKeyShare<ModularInt>>(
        std::move(secret_key_share.value()));
  }

  // Compute constants for the ambient ring R' = Z[X]/(Q, X^{2N} + 1).
  // The ring R' is used to verify the "wrap around" polynomial such that
  //    pk(X) = -a(X) * s(X) + e(X) + wrap_around(X) * (X^N + 1)
  // holds over the mod-Q polynomial ring Z[X]/(Q), where a(X), s(X), e(X), and
  // pk(X) are all polynomials in the ring R = Z[X]/(Q, X^N + 1). Since a(X) and
  // s(X) both have degree at most N-1, their product over Z[X]/(Q) has degree
  // 2N-2 and hence the same result can be computed efficiently in R'.
  void SetUpAmbientRnsContext() {
    testing::MpaheParameters<ModularInt> params =
        testing::GetMultiPartyParameters<ModularInt>();
    auto context = RnsContext<ModularInt>::CreateForBfv(
                       params.log_n + 2, params.qs, params.ps, params.t)
                       .value();
    ambient_rns_context_ =
        std::make_unique<const RnsContext<ModularInt>>(std::move(context));
    ambient_moduli_ = ambient_rns_context_->MainPrimeModuli();
  }

  // Returns the polynomial `a` lifted to the ambient ring R'.
  RnsPolynomial<ModularInt> LiftToAmbientRing(RnsPolynomial<ModularInt> a) {
    if (a.IsNttForm()) {
      auto status = a.ConvertToCoeffForm(moduli_);
      CHECK(status.ok());
    }
    int m = 1 << ambient_rns_context_->LogN();
    std::vector<std::vector<ModularInt>> a_ambient_coeffs = a.Coeffs();
    for (int i = 0; i < a.NumModuli(); ++i) {
      auto mod_params_qi = moduli_[i]->ModParams();
      auto zero_mod_qi = ModularInt::ImportZero(mod_params_qi);
      a_ambient_coeffs[i].resize(m, zero_mod_qi);
    }
    auto amb_a = RnsPolynomial<ModularInt>::Create(std::move(a_ambient_coeffs),
                                                   /*is_ntt=*/false);
    CHECK(amb_a.ok());
    auto status = amb_a->ConvertToNttForm(ambient_moduli_);
    CHECK(status.ok());
    return amb_a.value();
  }

  std::unique_ptr<const RnsContext<ModularInt>> rns_context_;
  std::vector<const PrimeModulus<ModularInt>*> moduli_;
  std::unique_ptr<Prng> prng_;
  std::unique_ptr<const PublicParameter<ModularInt>> public_parameter_;
  std::unique_ptr<const SecretKeyShare<ModularInt>> secret_key_share_;
  std::unique_ptr<const RnsContext<ModularInt>> ambient_rns_context_;
  std::vector<const PrimeModulus<ModularInt>*> ambient_moduli_;
};

TYPED_TEST_SUITE(PublicKeyShareTest, testing::ModularIntTypesForMultiParty);

class PublicKeyShareNegativeTest : public PublicKeyShareTest<ModularInt32> {};

TEST_F(PublicKeyShareNegativeTest, CreateFailsIfSecretKeyShareIsNull) {
  SetupPublicParameter();
  EXPECT_THAT(PublicKeyShare<ModularInt32>::Create(
                  /*secret_key_share=*/nullptr, this->public_parameter_.get(),
                  kPrngType),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`secret_key_share` must not be null")));
}

TEST_F(PublicKeyShareNegativeTest, CreateFailsIfPublicParameterIsNull) {
  SetupSecretKeyShare();
  EXPECT_THAT(PublicKeyShare<ModularInt32>::Create(secret_key_share_.get(),
                                                   /*public_parameter=*/nullptr,
                                                   kPrngType),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("`public_parameter` must not be null")));
}

TEST_F(PublicKeyShareNegativeTest, CreateFailsIfInvalidPrngType) {
  SetupPublicParameter();
  SetupSecretKeyShare();
  EXPECT_THAT(
      PublicKeyShare<ModularInt32>::Create(
          secret_key_share_.get(), public_parameter_.get(), PRNG_TYPE_INVALID),
      StatusIs(absl::StatusCode::kInvalidArgument,
               HasSubstr("`prng_type` not specified correctly")));
}

TYPED_TEST(PublicKeyShareTest, Create) {
  this->SetupPublicParameter();
  this->SetupSecretKeyShare();
  ASSERT_OK_AND_ASSIGN(PublicKeyShare<TypeParam> public_key_share,
                       PublicKeyShare<TypeParam>::Create(
                           this->secret_key_share_.get(),
                           this->public_parameter_.get(), kPrngType));
  const RnsPolynomial<TypeParam>& b = public_key_share.ComponentB();
  EXPECT_EQ(b.NumCoeffs(), 1 << this->rns_context_->LogN());
  EXPECT_EQ(b.NumModuli(), this->rns_context_->NumMainPrimeModuli());
  EXPECT_TRUE(b.IsNttForm());
}

TYPED_TEST(PublicKeyShareTest, CreateExplicit) {
  this->SetupPublicParameter();
  this->SetupSecretKeyShare();
  this->SetUpAmbientRnsContext();

  ASSERT_OK_AND_ASSIGN(auto prng_error_seed,
                       SingleThreadHkdfPrng::GenerateSeed());
  ASSERT_OK_AND_ASSIGN(auto prng_error,
                       SingleThreadHkdfPrng::Create(prng_error_seed));

  auto b = RnsPolynomial<TypeParam>::CreateEmpty();
  auto e = RnsPolynomial<TypeParam>::CreateEmpty();
  auto wrap_around = RnsPolynomial<TypeParam>::CreateEmpty();

  // The public key share is a polynomial `b` such that b = -a * s + e \in R
  // for R = Z[X]/(Q, X^N+1), where `a` is given by public_parameter and `s` is
  // the secret key share polynomial. We also capture the polynomial wrap_around
  // such that b = -a * s + e + wrap_around * (X^N + 1) over Z[X]/(Q).
  ASSERT_OK(PublicKeyShare<TypeParam>::CreateExplicit(
      this->secret_key_share_->Key(), this->public_parameter_.get(),
      prng_error.get(), &b, &e, &wrap_around));

  EXPECT_EQ(b.NumCoeffs(), 1 << this->rns_context_->LogN());
  EXPECT_EQ(b.NumModuli(), this->rns_context_->NumMainPrimeModuli());
  EXPECT_TRUE(b.IsNttForm());

  // Error is indeed populated.
  EXPECT_EQ(e.NumCoeffs(), 1 << this->rns_context_->LogN());
  EXPECT_EQ(e.NumModuli(), this->rns_context_->NumMainPrimeModuli());
  EXPECT_TRUE(e.IsNttForm());

  // wrap_around is correctly computed: it should hold over Z[X]/(Q) that
  // b = -a * s + e + wrap_around * (X^N + 1). We verify this using the larger
  // quotient ring R' = Z[X]/(Q, X^{2N} + 1), by lifting all polynomials to R'
  // and check equality there. This is faithful as a * s has degree <= 2N-2 and
  // so RHS of the equation has the same value in Z[X]/(Q) and in R'.
  ASSERT_OK_AND_ASSIGN(
      auto u,
      this->public_parameter_->PublicKeyComponentA().Negate(this->moduli_));
  auto amb_u = this->LiftToAmbientRing(u);
  auto amb_s = this->LiftToAmbientRing(this->secret_key_share_->Key());
  auto amb_e = this->LiftToAmbientRing(e);
  auto amb_wrap_around = this->LiftToAmbientRing(wrap_around);
  auto amb_b = this->LiftToAmbientRing(b);

  int n = 1 << this->rns_context_->LogN();          // Degree of R.
  int m = 1 << this->ambient_rns_context_->LogN();  // Degree of R'.

  // Generate f(X) = X^N + 1 as a polynomial in R'.
  std::vector<std::vector<TypeParam>> f_coeffs(this->ambient_moduli_.size());
  for (int i = 0; i < this->ambient_moduli_.size(); ++i) {
    auto mod_params_qi = this->ambient_moduli_[i]->ModParams();
    f_coeffs[i].resize(m, TypeParam::ImportZero(mod_params_qi));
    f_coeffs[i][0] = TypeParam::ImportOne(mod_params_qi);
    f_coeffs[i][n] = TypeParam::ImportOne(mod_params_qi);
  }
  ASSERT_OK_AND_ASSIGN(auto amb_f, RnsPolynomial<TypeParam>::Create(
                                       std::move(f_coeffs), /*is_ntt=*/false));
  ASSERT_OK(amb_f.ConvertToNttForm(this->ambient_moduli_));

  // RHS = -a * s + e + wrap_around * (X^N + 1) \in R'.
  ASSERT_OK_AND_ASSIGN(auto rhs, amb_u.Mul(amb_s, this->ambient_moduli_));
  ASSERT_OK(rhs.AddInPlace(amb_e, this->ambient_moduli_));
  ASSERT_OK(
      rhs.FusedMulAddInPlace(amb_wrap_around, amb_f, this->ambient_moduli_));
  ASSERT_OK(rhs.ConvertToCoeffForm(this->ambient_moduli_));
  // RHS = b \in R'.
  ASSERT_OK(amb_b.ConvertToCoeffForm(this->ambient_moduli_));
  // Verify LHS == RHS.
  ASSERT_EQ(amb_b, rhs);
}

TYPED_TEST(PublicKeyShareTest, SerializeDeserialize) {
  this->SetupPublicParameter();
  this->SetupSecretKeyShare();
  ASSERT_OK_AND_ASSIGN(PublicKeyShare<TypeParam> public_key_share,
                       PublicKeyShare<TypeParam>::Create(
                           this->secret_key_share_.get(),
                           this->public_parameter_.get(), kPrngType));
  ASSERT_OK_AND_ASSIGN(SerializedPublicKeyShare serialized,
                       public_key_share.Serialize());

  auto moduli = this->public_parameter_->Moduli();
  ASSERT_OK_AND_ASSIGN(
      PublicKeyShare<TypeParam> deserialized,
      PublicKeyShare<TypeParam>::Deserialize(serialized, moduli));
  EXPECT_EQ(public_key_share.ComponentB(), deserialized.ComponentB());
}

}  // namespace
}  // namespace multi_party
}  // namespace rlwe
