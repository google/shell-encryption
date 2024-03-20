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

# Error distribution samplers

package(default_visibility = ["//visibility:public"])

# Discrete Gaussian sampler implementations for lattice crypto.
cc_library(
    name = "discrete_gaussian",
    srcs = ["discrete_gaussian.cc"],
    hdrs = ["discrete_gaussian.h"],
    deps = [
        "//shell_encryption:integral_types",
        "//shell_encryption:statusor_fork",
        "//shell_encryption/prng",
        "@com_google_absl//absl/memory",
        "@com_google_absl//absl/numeric:int128",
        "@com_google_absl//absl/status",
        "@com_google_absl//absl/status:statusor",
        "@com_google_absl//absl/strings",
    ],
)

cc_test(
    name = "discrete_gaussian_test",
    srcs = ["discrete_gaussian_test.cc"],
    deps = [
        ":discrete_gaussian",
        "//shell_encryption:integral_types",
        "//shell_encryption/testing:status_is_fork",
        "//shell_encryption/testing:status_testing",
        "//shell_encryption/testing:testing_prng",
        "@com_github_google_googletest//:gtest_main",
        "@com_google_absl//absl/numeric:int128",
        "@com_google_absl//absl/status",
    ],
)

# Sampling a uniform ternary vector.
cc_library(
    name = "uniform_ternary",
    hdrs = ["uniform_ternary.h"],
    deps = [
        "//shell_encryption:integral_types",
        "//shell_encryption:statusor_fork",
        "//shell_encryption/prng",
        "@com_google_absl//absl/status",
        "@com_google_absl//absl/status:statusor",
    ],
)

cc_test(
    name = "uniform_ternary_test",
    srcs = ["uniform_ternary_test.cc"],
    deps = [
        ":uniform_ternary",
        "//shell_encryption:context",
        "//shell_encryption:integral_types",
        "//shell_encryption:montgomery",
        "//shell_encryption/testing:parameters",
        "//shell_encryption/testing:status_is_fork",
        "//shell_encryption/testing:status_testing",
        "//shell_encryption/testing:testing_prng",
        "@com_github_google_googletest//:gtest_main",
        "@com_google_absl//absl/status",
    ],
)
