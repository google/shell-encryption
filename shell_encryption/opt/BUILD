# Copyright 2021 Google LLC.
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

licenses(["notice"])

cc_library(
    name = "constant_polynomial",
    hdrs = ["constant_polynomial.h"],
    deps = [
        "//shell_encryption:statusor_fork",
        "@com_google_absl//absl/status",
    ],
)

cc_test(
    name = "constant_polynomial_test",
    srcs = ["constant_polynomial_test.cc"],
    deps = [
        ":constant_polynomial",
        "//shell_encryption:context",
        "//shell_encryption:polynomial",
        "//shell_encryption/testing:parameters",
        "//shell_encryption/testing:status_is_fork",
        "//shell_encryption/testing:status_testing",
        "@com_github_google_googletest//:gtest_main",
    ],
)

cc_library(
    name = "lazy_polynomial",
    hdrs = ["lazy_polynomial.h"],
    deps = [
        "//shell_encryption:statusor_fork",
        "@com_google_absl//absl/status",
        "@com_google_absl//absl/strings",
    ],
)

cc_test(
    name = "lazy_polynomial_test",
    srcs = ["lazy_polynomial_test.cc"],
    deps = [
        ":lazy_polynomial",
        "//shell_encryption:context",
        "//shell_encryption:montgomery",
        "//shell_encryption:polynomial",
        "//shell_encryption:statusor_fork",
        "//shell_encryption/prng:single_thread_hkdf_prng",
        "//shell_encryption/testing:parameters",
        "//shell_encryption/testing:status_is_fork",
        "//shell_encryption/testing:status_testing",
        "@com_github_google_googletest//:gtest_main",
        "@com_google_absl//absl/status",
    ],
)
