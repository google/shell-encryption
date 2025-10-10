# Copyright 2017 Google LLC.
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

load("@rules_cc//cc:defs.bzl", "cc_library")
load("@rules_proto//proto:defs.bzl", "proto_library")
load("//third_party/bazel_rules/rules_cc/cc:cc_library.bzl", "cc_library")
load("//third_party/bazel_rules/rules_cc/cc:cc_test.bzl", "cc_test")

package(default_visibility = ["//visibility:public"])

licenses(["notice"])

# Utilities.

cc_library(
    name = "parameters",
    testonly = 1,
    hdrs = ["parameters.h"],
    deps = [
        "//shell_encryption:constants",
        "//shell_encryption:context",
        "//shell_encryption:integral_types",
        "//shell_encryption:montgomery",
        "@com_github_google_googletest//:gtest",
        "@com_google_absl//absl/numeric:int128",
    ],
)

cc_library(
    name = "testing_prng",
    testonly = 1,
    hdrs = ["testing_prng.h"],
    deps = [
        "//shell_encryption:integral_types",
        "//shell_encryption/prng",
    ],
)

cc_library(
    name = "status_testing",
    testonly = 1,
    hdrs = ["status_testing.h"],
    deps = ["//shell_encryption:statusor_fork"],
)

cc_library(
    name = "testing_utils",
    testonly = 1,
    hdrs = ["testing_utils.h"],
    deps = [
        ":testing_prng",
        "//shell_encryption:constants",
        "//shell_encryption:montgomery",
        "//shell_encryption:polynomial",
        "//shell_encryption:statusor_fork",
        "//shell_encryption:symmetric_encryption",
        "//shell_encryption/prng:single_thread_chacha_prng",
        "//shell_encryption/prng:single_thread_hkdf_prng",
    ],
)

cc_library(
    name = "status_is_fork",
    testonly = 1,
    hdrs = [
        "protobuf_matchers.h",
        "status_matchers.h",
    ],
    deps = [
        "//shell_encryption:statusor_fork",
        "@com_github_google_googletest//:gtest",
        "@com_google_absl//absl/status",
        "@com_google_protobuf//:differencer",
        "@com_google_protobuf//:protobuf",
    ],
)

cc_library(
    name = "matchers",
    testonly = 1,
    hdrs = [
        "protobuf_matchers.h",
        "status_matchers.h",
    ],
    deps = [
        "//shell_encryption:statusor_fork",
        "@com_github_google_googletest//:gtest",
        "@com_google_absl//absl/status",
        "@com_google_protobuf//:differencer",
        "@com_google_protobuf//:protobuf",
    ],
)

cc_test(
    name = "protobuf_matchers_test",
    srcs = ["protobuf_matchers_test.cc"],
    deps = [
        ":coefficient_polynomial_cc_proto",
        ":matchers",
        ":status_testing",
        "@com_github_google_googletest//:gtest_main",
    ],
)

# Polynomials represented by their coefficient vector.

cc_library(
    name = "coefficient_polynomial",
    testonly = 1,
    hdrs = ["coefficient_polynomial.h"],
    deps = [
        ":coefficient_polynomial_cc_proto",
        "//shell_encryption:statusor_fork",
        "@com_github_google_glog//:glog",
        "@com_google_absl//absl/strings",
    ],
)

cc_test(
    name = "coefficient_polynomial_test",
    srcs = [
        "coefficient_polynomial_test.cc",
    ],
    deps = [
        ":coefficient_polynomial",
        ":status_is_fork",
        ":status_testing",
        ":testing_prng",
        "//shell_encryption:constants",
        "//shell_encryption:montgomery",
        "@com_github_google_googletest//:gtest_main",
    ],
)

cc_library(
    name = "coefficient_polynomial_ciphertext",
    testonly = 1,
    hdrs = ["coefficient_polynomial_ciphertext.h"],
    deps = [
        ":coefficient_polynomial",
        "//shell_encryption:error_params",
        "//shell_encryption:polynomial",
        "//shell_encryption:statusor_fork",
        "//shell_encryption:symmetric_encryption",
    ],
)

cc_test(
    name = "coefficient_polynomial_ciphertext_test",
    srcs = [
        "coefficient_polynomial_ciphertext_test.cc",
    ],
    deps = [
        ":coefficient_polynomial",
        ":coefficient_polynomial_ciphertext",
        ":status_is_fork",
        ":status_testing",
        ":testing_prng",
        ":testing_utils",
        "//shell_encryption:constants",
        "//shell_encryption:montgomery",
        "//shell_encryption:ntt_parameters",
        "//shell_encryption:polynomial",
        "//shell_encryption:symmetric_encryption",
        "@com_github_google_googletest//:gtest_main",
    ],
)

proto_library(
    name = "coefficient_polynomial_proto",
    srcs = ["coefficient_polynomial.proto"],
)

cc_proto_library(
    name = "coefficient_polynomial_cc_proto",
    deps = [":coefficient_polynomial_proto"],
)
