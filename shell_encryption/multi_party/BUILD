# Copyright 2025 Google LLC.
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

# Multi-party additive homomorphic encryption protocol

load("@rules_cc//cc:defs.bzl", "cc_library")
load("@rules_proto//proto:defs.bzl", "proto_library")

package(default_visibility = ["//visibility:public"])

proto_library(
    name = "serialization_proto",
    srcs = ["serialization.proto"],
    deps = [
        "//shell_encryption:serialization_proto",
        "//shell_encryption/rns:serialization_proto",
    ],
)

cc_proto_library(
    name = "serialization_cc_proto",
    deps = [":serialization_proto"],
)

# Public randomness.
cc_library(
    name = "public_parameter",
    srcs = ["public_parameter.cc"],
    hdrs = ["public_parameter.h"],
    deps = [
        ":serialization_cc_proto",
        "//shell_encryption:integral_types",
        "//shell_encryption:montgomery",
        "//shell_encryption:statusor_fork",
        "//shell_encryption/prng",
        "//shell_encryption/prng:single_thread_chacha_prng",
        "//shell_encryption/prng:single_thread_hkdf_prng",
        "//shell_encryption/rns:rns_context",
        "//shell_encryption/rns:rns_modulus",
        "//shell_encryption/rns:rns_polynomial",
        "@com_google_absl//absl/memory",
        "@com_google_absl//absl/numeric:int128",
        "@com_google_absl//absl/status",
        "@com_google_absl//absl/status:statusor",
        "@com_google_absl//absl/types:span",
    ],
)

cc_test(
    name = "public_parameter_test",
    size = "small",
    srcs = ["public_parameter_test.cc"],
    deps = [
        ":public_parameter",
        "//shell_encryption/multi_party/testing:parameters",
        "//shell_encryption/rns:rns_context",
        "//shell_encryption/rns:rns_modulus",
        "//shell_encryption/rns:rns_polynomial",
        "//shell_encryption/testing:matchers",
        "//shell_encryption/testing:status_testing",
        "@com_github_google_googletest//:gtest_main",
        "@com_google_absl//absl/log:check",
        "@com_google_absl//absl/status",
    ],
)

# Secret key shares.
cc_library(
    name = "secret_key_share",
    srcs = ["secret_key_share.cc"],
    hdrs = ["secret_key_share.h"],
    deps = [
        "//shell_encryption:integral_types",
        "//shell_encryption:montgomery",
        "//shell_encryption:statusor_fork",
        "//shell_encryption/prng",
        "//shell_encryption/rns:error_distribution",
        "//shell_encryption/rns:rns_context",
        "//shell_encryption/rns:rns_modulus",
        "//shell_encryption/rns:rns_polynomial",
        "//shell_encryption/sampler:discrete_gaussian",
        "@com_google_absl//absl/numeric:int128",
        "@com_google_absl//absl/status",
        "@com_google_absl//absl/status:statusor",
        "@com_google_absl//absl/types:span",
    ],
)

cc_test(
    name = "secret_key_share_test",
    size = "small",
    srcs = ["secret_key_share_test.cc"],
    deps = [
        ":secret_key_share",
        "//shell_encryption/multi_party/testing:parameters",
        "//shell_encryption/rns:rns_context",
        "//shell_encryption/rns:rns_modulus",
        "//shell_encryption/rns:rns_polynomial",
        "//shell_encryption/sampler:discrete_gaussian",
        "//shell_encryption/testing:matchers",
        "//shell_encryption/testing:status_testing",
        "//shell_encryption/testing:testing_prng",
        "@com_github_google_googletest//:gtest_main",
        "@com_google_absl//absl/log:check",
        "@com_google_absl//absl/status",
    ],
)

# Public key shares.
cc_library(
    name = "public_key_share",
    srcs = ["public_key_share.cc"],
    hdrs = ["public_key_share.h"],
    deps = [
        ":public_parameter",
        ":secret_key_share",
        "//shell_encryption:integral_types",
        "//shell_encryption:montgomery",
        "//shell_encryption:statusor_fork",
        "//shell_encryption/prng",
        "//shell_encryption/prng:single_thread_chacha_prng",
        "//shell_encryption/prng:single_thread_hkdf_prng",
        "//shell_encryption/rns:error_distribution",
        "//shell_encryption/rns:rns_context",
        "//shell_encryption/rns:rns_modulus",
        "//shell_encryption/rns:rns_polynomial",
        "@com_google_absl//absl/numeric:int128",
        "@com_google_absl//absl/status",
        "@com_google_absl//absl/status:statusor",
        "@com_google_absl//absl/types:span",
    ],
)

cc_test(
    name = "public_key_share_test",
    size = "small",
    srcs = ["public_key_share_test.cc"],
    deps = [
        ":public_key_share",
        ":public_parameter",
        ":secret_key_share",
        "//shell_encryption/multi_party/testing:parameters",
        "//shell_encryption/prng:single_thread_hkdf_prng",
        "//shell_encryption/rns:rns_context",
        "//shell_encryption/rns:rns_modulus",
        "//shell_encryption/rns:rns_polynomial",
        "//shell_encryption/testing:matchers",
        "//shell_encryption/testing:status_testing",
        "//shell_encryption/testing:testing_prng",
        "@com_github_google_googletest//:gtest_main",
        "@com_google_absl//absl/log:check",
        "@com_google_absl//absl/status",
    ],
)

# Public key.
cc_library(
    name = "public_key",
    srcs = ["public_key.cc"],
    hdrs = ["public_key.h"],
    deps = [
        ":public_key_share",
        ":public_parameter",
        "//shell_encryption:integral_types",
        "//shell_encryption:montgomery",
        "//shell_encryption:statusor_fork",
        "//shell_encryption/prng",
        "//shell_encryption/rns:coefficient_encoder",
        "//shell_encryption/rns:error_distribution",
        "//shell_encryption/rns:rns_bfv_ciphertext",
        "//shell_encryption/rns:rns_error_params",
        "//shell_encryption/rns:rns_modulus",
        "//shell_encryption/rns:rns_polynomial",
        "@com_google_absl//absl/numeric:int128",
        "@com_google_absl//absl/status",
        "@com_google_absl//absl/status:statusor",
        "@com_google_absl//absl/types:span",
    ],
)

cc_test(
    name = "public_key_test",
    size = "small",
    srcs = ["public_key_test.cc"],
    deps = [
        ":public_key",
        ":public_key_share",
        ":public_parameter",
        ":secret_key_share",
        "//shell_encryption/multi_party/testing:parameters",
        "//shell_encryption/rns:coefficient_encoder",
        "//shell_encryption/rns:rns_bfv_ciphertext",
        "//shell_encryption/rns:rns_context",
        "//shell_encryption/rns:rns_error_params",
        "//shell_encryption/rns:rns_modulus",
        "//shell_encryption/rns:rns_polynomial",
        "//shell_encryption/rns/testing:testing_utils",
        "//shell_encryption/testing:matchers",
        "//shell_encryption/testing:status_testing",
        "//shell_encryption/testing:testing_prng",
        "@com_github_google_googletest//:gtest_main",
        "@com_google_absl//absl/log:check",
        "@com_google_absl//absl/status",
        "@com_google_absl//absl/strings",
    ],
)

# Message recovery from partial decryptions.
cc_library(
    name = "recovery",
    hdrs = ["recovery.h"],
    deps = [
        ":public_parameter",
        "//shell_encryption:statusor_fork",
        "//shell_encryption/rns:coefficient_encoder",
        "//shell_encryption/rns:rns_polynomial",
        "@com_google_absl//absl/status",
        "@com_google_absl//absl/status:statusor",
        "@com_google_absl//absl/types:span",
    ],
)

cc_test(
    name = "recovery_test",
    size = "small",
    srcs = ["recovery_test.cc"],
    deps = [
        ":public_key",
        ":public_key_share",
        ":public_parameter",
        ":recovery",
        ":secret_key_share",
        "//shell_encryption/multi_party/testing:parameters",
        "//shell_encryption/rns:coefficient_encoder",
        "//shell_encryption/rns:rns_bfv_ciphertext",
        "//shell_encryption/rns:rns_context",
        "//shell_encryption/rns:rns_error_params",
        "//shell_encryption/rns:rns_modulus",
        "//shell_encryption/rns:rns_polynomial",
        "//shell_encryption/rns/testing:testing_utils",
        "//shell_encryption/sampler:discrete_gaussian",
        "//shell_encryption/testing:matchers",
        "//shell_encryption/testing:status_testing",
        "//shell_encryption/testing:testing_prng",
        "@com_github_google_googletest//:gtest_main",
        "@com_google_absl//absl/log:check",
        "@com_google_absl//absl/status",
        "@com_google_absl//absl/strings",
    ],
)
