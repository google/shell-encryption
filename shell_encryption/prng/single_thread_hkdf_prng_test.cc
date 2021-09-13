/*
 * Copyright 2021 Google LLC.
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

// This file includes specific unit tests files to test populating the buffer
// every internal::kHkdfMaxOutputBytes bytes.

#include "shell_encryption/prng/single_thread_hkdf_prng.h"

#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include "absl/strings/str_cat.h"
#include "shell_encryption/testing/status_matchers.h"
#include "shell_encryption/testing/status_testing.h"

namespace rlwe {

namespace {

using ::testing::Eq;
using ::testing::Not;

class SingleThreadHkdfPrngTest : public ::testing::Test {
 protected:
  void SetUp() override {
    ASSERT_OK_AND_ASSIGN(seed_, SingleThreadHkdfPrng::GenerateSeed());
    ASSERT_OK_AND_ASSIGN(prng_, SingleThreadHkdfPrng::Create(seed_));
  }

  std::string seed_;
  std::unique_ptr<SingleThreadHkdfPrng> prng_;
};

TEST_F(SingleThreadHkdfPrngTest, TestRand8BeforeAndAfterResalting) {
  // Two random 8 bit strings have 1/256 probability of being equal. Instead,
  // we check that a sequence of 8 strings from the HKDF before and after
  // salting are not all equal.
  std::vector<std::vector<Uint8>> before_resalt(internal::kHkdfMaxOutputBytes);
  for (int i = 0; i < internal::kHkdfMaxOutputBytes / 8; ++i) {
    SCOPED_TRACE(absl::StrCat("Iteration ", i, "."));
    std::vector<Uint8> rand8s;
    for (int j = 0; j < 8; ++j) {
      ASSERT_OK_AND_ASSIGN(auto elt, prng_->Rand8());
      rand8s.push_back(elt);
    }
    before_resalt.push_back(rand8s);
  }
  for (int i = 0; i < internal::kHkdfMaxOutputBytes / 8; ++i) {
    SCOPED_TRACE(absl::StrCat("Iteration ", i, "."));
    std::vector<Uint8> rand8s;
    for (int j = 0; j < 8; ++j) {
      ASSERT_OK_AND_ASSIGN(auto elt, prng_->Rand8());
      rand8s.push_back(elt);
    }
    EXPECT_THAT(rand8s, Not(Eq(before_resalt[i])));
  }
}

TEST_F(SingleThreadHkdfPrngTest, TestRand64BeforeAndAfterResalting) {
  std::vector<Uint64> before_resalt(internal::kHkdfMaxOutputBytes / 8);
  for (int i = 0; i < internal::kHkdfMaxOutputBytes / 8; ++i) {
    SCOPED_TRACE(absl::StrCat("Iteration ", i, "."));
    ASSERT_OK_AND_ASSIGN(before_resalt[i], prng_->Rand64());
  }
  for (int i = 0; i < internal::kHkdfMaxOutputBytes / 8; ++i) {
    SCOPED_TRACE(absl::StrCat("Iteration ", i, "."));
    ASSERT_OK_AND_ASSIGN(auto r64, prng_->Rand64());
    EXPECT_NE(before_resalt[i], r64);
  }
}

TEST_F(SingleThreadHkdfPrngTest, TestRand64WithResaltingInBetween) {
  for (int i = 0; i < internal::kHkdfMaxOutputBytes - 1; ++i) {
    SCOPED_TRACE(absl::StrCat("Iteration ", i, "."));
    EXPECT_OK(prng_->Rand8());
  }
  EXPECT_OK(prng_->Rand64());
}

}  // namespace
}  // namespace rlwe
