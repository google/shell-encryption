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

#include "shell_encryption/rns/rns_ciphertext.h"

#include "absl/numeric/int128.h"
#include "shell_encryption/integral_types.h"
#include "shell_encryption/montgomery.h"

namespace rlwe {

template class RnsRlweCiphertext<MontgomeryInt<Uint16>>;
template class RnsRlweCiphertext<MontgomeryInt<Uint32>>;
template class RnsRlweCiphertext<MontgomeryInt<Uint64>>;
template class RnsRlweCiphertext<MontgomeryInt<absl::uint128>>;
#ifdef ABSL_HAVE_INTRINSIC_INT128
template class RnsRlweCiphertext<MontgomeryInt<unsigned __int128>>;
#endif

}  // namespace rlwe
