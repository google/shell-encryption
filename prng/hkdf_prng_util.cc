// Copyright 2021 Google LLC
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

#include "prng/hkdf_prng_util.h"

#include <vector>

#include "absl/memory/memory.h"
#include "absl/strings/str_cat.h"
#include "status_macros.h"
#include "tink/subtle/common_enums.h"
#include "tink/subtle/hkdf.h"
#include "tink/subtle/random.h"

namespace rlwe::internal {

absl::Status HkdfPrngResalt(absl::string_view key, int buffer_size,
                            int* salt_counter, int* position_in_buffer,
                            std::vector<Uint8>* buffer) {
  std::string salt = absl::StrCat("salt", *salt_counter);
  RLWE_ASSIGN_OR_RETURN(
      auto buf, crypto::tink::subtle::Hkdf::ComputeHkdf(
                    crypto::tink::subtle::SHA256, key, salt, "", buffer_size));
  buffer->assign(buf.begin(), buf.end());
  ++(*salt_counter);
  *position_in_buffer = 0;

  return absl::OkStatus();
}

rlwe::StatusOr<std::string> HkdfPrngGenerateKey() {
  return crypto::tink::subtle::Random::GetRandomBytes(kHkdfKeyBytesSize);
}

rlwe::StatusOr<Uint8> HkdfPrngRand8(absl::string_view key,
                                    int* position_in_buffer, int* salt_counter,
                                    std::vector<Uint8>* buffer) {
  Uint8 rand;
  if (*position_in_buffer >= buffer->size()) {
    RLWE_RETURN_IF_ERROR(HkdfPrngResalt(key, kHkdfMaxOutputBytes, salt_counter,
                                        position_in_buffer, buffer));
  }
  rand = buffer->at(*position_in_buffer);
  ++(*position_in_buffer);
  return rand;
}

rlwe::StatusOr<Uint64> HkdfPrngRand64(absl::string_view key,
                                      int* position_in_buffer,
                                      int* salt_counter,
                                      std::vector<Uint8>* buffer) {
  Uint64 rand64 = 0;
  for (int i = 0; i < 8; ++i) {
    RLWE_ASSIGN_OR_RETURN(Uint8 rand8, HkdfPrngRand8(key, position_in_buffer,
                                                     salt_counter, buffer));
    rand64 += Uint64{rand8} << (8 * i);
  }
  return rand64;
}

}  // namespace rlwe::internal
