# Copyright 2023 Google LLC.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     https://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

package(default_visibility = ["//visibility:public"])

cc_library(
    name = "parameters",
    testonly = 1,
    hdrs = ["parameters.h"],
    deps = [
        "//shell_encryption:integral_types",
        "//shell_encryption:montgomery",
        "@com_github_google_googletest//:gtest",
        "@com_google_absl//absl/numeric:int128",
    ],
)

cc_library(
    name = "testing_utils",
    testonly = 1,
    hdrs = ["testing_utils.h"],
    deps = [
        "//shell_encryption:integral_types",
        "//shell_encryption/rns:rns_modulus",
        "//shell_encryption/rns:rns_polynomial",
        "@com_google_absl//absl/log:check",
        "@com_google_absl//absl/numeric:int128",
        "@com_google_absl//absl/random",
        "@com_google_absl//absl/types:span",
    ],
)
