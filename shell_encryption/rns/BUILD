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

# RNS extension of the rlwe library

load("@rules_cc//cc:defs.bzl", "cc_library")
load("@rules_proto//proto:defs.bzl", "proto_library")

package(default_visibility = ["//visibility:public"])

proto_library(
    name = "serialization_proto",
    srcs = ["serialization.proto"],
    deps = ["//shell_encryption:serialization_proto"],
)

cc_proto_library(
    name = "serialization_cc_proto",
    deps = [":serialization_proto"],
)

# RNS integer.
cc_library(
    name = "rns_integer",
    hdrs = ["rns_integer.h"],
    deps = [
        "@com_google_absl//absl/types:span",
    ],
)

# RNS modulus and its parameters.
cc_library(
    name = "rns_modulus",
    hdrs = ["rns_modulus.h"],
    deps = [
        "//shell_encryption:montgomery",
        "//shell_encryption:ntt_parameters",
        "@com_google_absl//absl/types:span",
    ],
)

# RNS context.
cc_library(
    name = "rns_context",
    srcs = ["rns_context.cc"],
    hdrs = ["rns_context.h"],
    deps = [
        ":rns_integer",
        ":rns_modulus",
        "//shell_encryption:montgomery",
        "//shell_encryption:ntt_parameters",
        "@com_google_absl//absl/memory",
        "@com_google_absl//absl/types:span",
    ],
)

cc_test(
    name = "rns_context_test",
    size = "small",
    srcs = ["rns_context_test.cc"],
    deps = [
        ":rns_context",
        "//shell_encryption:integral_types",
        "//shell_encryption:montgomery",
        "//shell_encryption:ntt_parameters",
        "//shell_encryption:statusor_fork",
        "//shell_encryption/rns/testing:parameters",
        "//shell_encryption/rns/testing:testing_utils",
        "//shell_encryption/testing:matchers",
        "//shell_encryption/testing:parameters",
        "//shell_encryption/testing:status_testing",
        "@com_github_google_googletest//:gtest_main",
        "@com_google_absl//absl/strings",
    ],
)

# RNS polynomial.
cc_library(
    name = "rns_polynomial",
    srcs = ["rns_polynomial.cc"],
    hdrs = ["rns_polynomial.h"],
    deps = [
        ":rns_integer",
        ":rns_modulus",
        ":serialization_cc_proto",
        "//shell_encryption:dft_transformations",
        "//shell_encryption:modulus_conversion",
        "//shell_encryption:montgomery",
        "//shell_encryption:ntt_parameters",
        "//shell_encryption:polynomial",
        "//shell_encryption:statusor_fork",
        "@com_google_absl//absl/strings",
        "@com_google_absl//absl/types:span",
    ],
)

cc_test(
    name = "rns_polynomial_test",
    size = "medium",
    srcs = ["rns_polynomial_test.cc"],
    deps = [
        ":crt_interpolation",
        ":rns_context",
        ":rns_modulus",
        ":rns_polynomial",
        ":serialization_cc_proto",
        "//shell_encryption:integral_types",
        "//shell_encryption:montgomery",
        "//shell_encryption:polynomial",
        "//shell_encryption:statusor_fork",
        "//shell_encryption/rns/testing:parameters",
        "//shell_encryption/rns/testing:testing_utils",
        "//shell_encryption/testing:matchers",
        "//shell_encryption/testing:parameters",
        "//shell_encryption/testing:status_testing",
        "//shell_encryption/testing:testing_prng",
        "@com_github_google_googletest//:gtest_main",
        "@com_google_absl//absl/log:check",
        "@com_google_absl//absl/random",
        "@com_google_absl//absl/status",
        "@com_google_absl//absl/strings",
    ],
)

# CRT interpolation utilities
cc_library(
    name = "crt_interpolation",
    hdrs = ["crt_interpolation.h"],
    deps = [
        ":rns_modulus",
        "//shell_encryption:integral_types",
        "//shell_encryption:montgomery",
        "@com_google_absl//absl/numeric:int128",
        "@com_google_absl//absl/status",
        "@com_google_absl//absl/status:statusor",
    ],
)

cc_test(
    name = "crt_interpolation_test",
    srcs = ["crt_interpolation_test.cc"],
    deps = [
        ":crt_interpolation",
        ":rns_context",
        "//shell_encryption:int256",
        "//shell_encryption:integral_types",
        "//shell_encryption:montgomery",
        "//shell_encryption:statusor_fork",
        "//shell_encryption/rns/testing:parameters",
        "//shell_encryption/testing:matchers",
        "//shell_encryption/testing:parameters",
        "//shell_encryption/testing:status_testing",
        "@com_github_google_googletest//:gtest_main",
        "@com_google_absl//absl/numeric:int128",
        "@com_google_absl//absl/random",
        "@com_google_absl//absl/status",
    ],
)

# Error removal utilities.
cc_library(
    name = "error_correction",
    hdrs = ["error_correction.h"],
    deps = [
        "@com_google_absl//absl/status:statusor",
    ],
)

# Encode to and decode from RNS polynomial coefficients.
cc_library(
    name = "coefficient_encoder",
    srcs = ["coefficient_encoder.cc"],
    hdrs = ["coefficient_encoder.h"],
    deps = [
        ":error_correction",
        ":rns_context",
        ":rns_integer",
        ":rns_modulus",
        ":rns_polynomial",
        "//shell_encryption:integral_types",
        "//shell_encryption:modulus_conversion",
        "//shell_encryption:montgomery",
        "//shell_encryption:statusor_fork",
        "@com_google_absl//absl/numeric:int128",
        "@com_google_absl//absl/status",
        "@com_google_absl//absl/status:statusor",
        "@com_google_absl//absl/strings",
        "@com_google_absl//absl/types:span",
    ],
)

cc_test(
    name = "coefficient_encoder_test",
    srcs = ["coefficient_encoder_test.cc"],
    deps = [
        ":coefficient_encoder",
        ":rns_context",
        ":rns_integer",
        ":rns_modulus",
        ":rns_polynomial",
        "//shell_encryption:montgomery",
        "//shell_encryption/rns/testing:parameters",
        "//shell_encryption/rns/testing:testing_utils",
        "//shell_encryption/testing:matchers",
        "//shell_encryption/testing:parameters",
        "//shell_encryption/testing:status_testing",
        "@com_github_google_googletest//:gtest_main",
        "@com_google_absl//absl/log:check",
        "@com_google_absl//absl/status",
        "@com_google_absl//absl/types:span",
    ],
)

# Encode to and decode from plaintext slots for finite field values.
cc_library(
    name = "finite_field_encoder",
    srcs = ["finite_field_encoder.cc"],
    hdrs = ["finite_field_encoder.h"],
    deps = [
        ":coefficient_encoder",
        ":crt_interpolation",
        ":error_correction",
        ":rns_context",
        ":rns_modulus",
        ":rns_polynomial",
        "//shell_encryption:dft_transformations",
        "//shell_encryption:integral_types",
        "//shell_encryption:montgomery",
        "//shell_encryption:ntt_parameters",
        "//shell_encryption:statusor_fork",
        "@com_google_absl//absl/numeric:int128",
        "@com_google_absl//absl/status",
        "@com_google_absl//absl/status:statusor",
        "@com_google_absl//absl/strings",
        "@com_google_absl//absl/types:span",
    ],
)

cc_test(
    name = "finite_field_encoder_test",
    srcs = ["finite_field_encoder_test.cc"],
    deps = [
        ":crt_interpolation",
        ":finite_field_encoder",
        ":rns_context",
        ":rns_integer",
        ":rns_modulus",
        ":rns_polynomial",
        "//shell_encryption:montgomery",
        "//shell_encryption/rns/testing:parameters",
        "//shell_encryption/rns/testing:testing_utils",
        "//shell_encryption/testing:matchers",
        "//shell_encryption/testing:status_testing",
        "@com_github_google_googletest//:gtest_main",
        "@com_google_absl//absl/log:check",
        "@com_google_absl//absl/random",
        "@com_google_absl//absl/status",
        "@com_google_absl//absl/types:span",
    ],
)

# RNS error parameters
cc_library(
    name = "rns_error_params",
    hdrs = ["rns_error_params.h"],
    deps = [
        ":rns_modulus",
        "@com_google_absl//absl/status",
        "@com_google_absl//absl/status:statusor",
        "@com_google_absl//absl/types:span",
    ],
)

cc_test(
    name = "rns_error_params_test",
    srcs = ["rns_error_params_test.cc"],
    deps = [
        ":coefficient_encoder",
        ":crt_interpolation",
        ":error_distribution",
        ":rns_context",
        ":rns_error_params",
        ":rns_modulus",
        ":rns_polynomial",
        "//shell_encryption:integral_types",
        "//shell_encryption/prng:single_thread_hkdf_prng",
        "//shell_encryption/rns/testing:parameters",
        "//shell_encryption/rns/testing:testing_utils",
        "//shell_encryption/testing:matchers",
        "//shell_encryption/testing:parameters",
        "//shell_encryption/testing:status_testing",
        "@com_github_google_googletest//:gtest_main",
        "@com_google_absl//absl/log:check",
        "@com_google_absl//absl/random",
        "@com_google_absl//absl/status",
        "@com_google_absl//absl/types:span",
    ],
)

# General RLWE ciphertext in RNS form.
cc_library(
    name = "rns_ciphertext",
    srcs = ["rns_ciphertext.cc"],
    hdrs = ["rns_ciphertext.h"],
    deps = [
        ":lazy_rns_polynomial",
        ":rns_error_params",
        ":rns_polynomial",
        ":serialization_cc_proto",
        "//shell_encryption:integral_types",
        "//shell_encryption:montgomery",
        "//shell_encryption:statusor_fork",
        "@com_google_absl//absl/numeric:int128",
        "@com_google_absl//absl/strings",
        "@com_google_absl//absl/types:span",
    ],
)

cc_test(
    name = "rns_ciphertext_test",
    srcs = ["rns_ciphertext_test.cc"],
    deps = [
        ":coefficient_encoder",
        ":finite_field_encoder",
        ":rns_ciphertext",
        ":rns_context",
        ":rns_error_params",
        ":rns_integer",
        ":rns_modulus",
        ":rns_polynomial",
        ":rns_secret_key",
        ":serialization_cc_proto",
        "//shell_encryption:integral_types",
        "//shell_encryption/prng:single_thread_hkdf_prng",
        "//shell_encryption/rns/testing:parameters",
        "//shell_encryption/testing:matchers",
        "//shell_encryption/testing:parameters",
        "//shell_encryption/testing:status_testing",
        "@com_github_google_googletest//:gtest_main",
        "@com_google_absl//absl/log:check",
        "@com_google_absl//absl/random",
        "@com_google_absl//absl/status",
        "@com_google_absl//absl/strings",
        "@com_google_absl//absl/types:span",
    ],
)

# BGV ciphertext in RNS form.
cc_library(
    name = "rns_bgv_ciphertext",
    srcs = ["rns_bgv_ciphertext.cc"],
    hdrs = ["rns_bgv_ciphertext.h"],
    deps = [
        ":rns_ciphertext",
        ":rns_error_params",
        ":rns_polynomial",
        "//shell_encryption:integral_types",
        "//shell_encryption:montgomery",
        "//shell_encryption:statusor_fork",
        "@com_google_absl//absl/numeric:int128",
        "@com_google_absl//absl/strings",
        "@com_google_absl//absl/types:span",
    ],
)

cc_test(
    name = "rns_bgv_ciphertext_test",
    srcs = ["rns_bgv_ciphertext_test.cc"],
    deps = [
        ":coefficient_encoder",
        ":finite_field_encoder",
        ":rns_bgv_ciphertext",
        ":rns_context",
        ":rns_error_params",
        ":rns_integer",
        ":rns_modulus",
        ":rns_polynomial",
        ":rns_secret_key",
        "//shell_encryption/prng:single_thread_hkdf_prng",
        "//shell_encryption/rns/testing:parameters",
        "//shell_encryption/rns/testing:testing_utils",
        "//shell_encryption/testing:matchers",
        "//shell_encryption/testing:parameters",
        "//shell_encryption/testing:status_testing",
        "@com_github_google_googletest//:gtest_main",
        "@com_google_absl//absl/log:check",
        "@com_google_absl//absl/status",
        "@com_google_absl//absl/strings",
        "@com_google_absl//absl/types:span",
    ],
)

# BFV ciphertext in RNS form.
cc_library(
    name = "rns_bfv_ciphertext",
    srcs = ["rns_bfv_ciphertext.cc"],
    hdrs = ["rns_bfv_ciphertext.h"],
    deps = [
        ":rns_ciphertext",
        ":rns_context",
        ":rns_error_params",
        ":rns_integer",
        ":rns_modulus",
        ":rns_polynomial",
        "//shell_encryption:integral_types",
        "//shell_encryption:montgomery",
        "//shell_encryption:statusor_fork",
        "@com_google_absl//absl/numeric:int128",
        "@com_google_absl//absl/status",
        "@com_google_absl//absl/status:statusor",
        "@com_google_absl//absl/types:span",
    ],
)

cc_test(
    name = "rns_bfv_ciphertext_test",
    srcs = ["rns_bfv_ciphertext_test.cc"],
    deps = [
        ":coefficient_encoder",
        ":finite_field_encoder",
        ":rns_bfv_ciphertext",
        ":rns_context",
        ":rns_error_params",
        ":rns_integer",
        ":rns_modulus",
        ":rns_polynomial",
        ":rns_secret_key",
        "//shell_encryption/prng:single_thread_hkdf_prng",
        "//shell_encryption/rns/testing:parameters",
        "//shell_encryption/rns/testing:testing_utils",
        "//shell_encryption/testing:matchers",
        "//shell_encryption/testing:parameters",
        "//shell_encryption/testing:status_testing",
        "@com_github_google_googletest//:gtest_main",
        "@com_google_absl//absl/log:check",
        "@com_google_absl//absl/status",
        "@com_google_absl//absl/strings",
        "@com_google_absl//absl/types:span",
    ],
)

# Error distribution with RNS polynomial representation.
cc_library(
    name = "error_distribution",
    hdrs = ["error_distribution.h"],
    deps = [
        ":rns_modulus",
        ":rns_polynomial",
        "//shell_encryption:sample_error",
        "//shell_encryption:statusor_fork",
        "//shell_encryption/prng",
        "//shell_encryption/sampler:discrete_gaussian",
        "//shell_encryption/sampler:uniform_ternary",
        "@com_google_absl//absl/status",
        "@com_google_absl//absl/status:statusor",
        "@com_google_absl//absl/types:span",
    ],
)

cc_test(
    name = "error_distribution_test",
    size = "small",
    srcs = ["error_distribution_test.cc"],
    deps = [
        ":error_distribution",
        ":rns_context",
        ":rns_modulus",
        ":rns_polynomial",
        "//shell_encryption/prng:single_thread_hkdf_prng",
        "//shell_encryption/rns/testing:parameters",
        "//shell_encryption/sampler:discrete_gaussian",
        "//shell_encryption/testing:matchers",
        "//shell_encryption/testing:parameters",
        "//shell_encryption/testing:status_testing",
        "@com_github_google_googletest//:gtest_main",
        "@com_google_absl//absl/log:check",
        "@com_google_absl//absl/status",
        "@com_google_absl//absl/strings",
    ],
)

# RNS secret key.
cc_library(
    name = "rns_secret_key",
    srcs = ["rns_secret_key.cc"],
    hdrs = ["rns_secret_key.h"],
    deps = [
        ":coefficient_encoder",
        ":error_distribution",
        ":rns_bfv_ciphertext",
        ":rns_bgv_ciphertext",
        ":rns_ciphertext",
        ":rns_error_params",
        ":rns_modulus",
        ":rns_polynomial",
        "//shell_encryption:integral_types",
        "//shell_encryption:montgomery",
        "//shell_encryption:statusor_fork",
        "//shell_encryption/prng",
        "@com_google_absl//absl/numeric:int128",
        "@com_google_absl//absl/status",
        "@com_google_absl//absl/status:statusor",
        "@com_google_absl//absl/types:span",
    ],
)

cc_test(
    name = "rns_secret_key_test",
    size = "small",
    srcs = ["rns_secret_key_test.cc"],
    deps = [
        ":coefficient_encoder",
        ":finite_field_encoder",
        ":rns_bfv_ciphertext",
        ":rns_bgv_ciphertext",
        ":rns_context",
        ":rns_error_params",
        ":rns_modulus",
        ":rns_polynomial",
        ":rns_secret_key",
        "//shell_encryption:integral_types",
        "//shell_encryption/prng:single_thread_hkdf_prng",
        "//shell_encryption/rns/testing:parameters",
        "//shell_encryption/rns/testing:testing_utils",
        "//shell_encryption/testing:matchers",
        "//shell_encryption/testing:parameters",
        "//shell_encryption/testing:status_testing",
        "@com_github_google_googletest//:gtest_main",
        "@com_google_absl//absl/log:check",
        "@com_google_absl//absl/random",
        "@com_google_absl//absl/status",
        "@com_google_absl//absl/strings",
        "@com_google_absl//absl/types:span",
    ],
)

# RNS public key.
cc_library(
    name = "rns_public_key",
    srcs = ["rns_public_key.cc"],
    hdrs = ["rns_public_key.h"],
    deps = [
        ":error_distribution",
        ":rns_modulus",
        ":rns_polynomial",
        ":rns_secret_key",
        "//shell_encryption:integral_types",
        "//shell_encryption:montgomery",
        "//shell_encryption:statusor_fork",
        "//shell_encryption/prng",
        "//shell_encryption/prng:single_thread_chacha_prng",
        "//shell_encryption/prng:single_thread_hkdf_prng",
        "@com_google_absl//absl/numeric:int128",
        "@com_google_absl//absl/status",
        "@com_google_absl//absl/status:statusor",
        "@com_google_absl//absl/strings",
        "@com_google_absl//absl/types:span",
    ],
)

# BGV variant of public key.
cc_library(
    name = "rns_bgv_public_key",
    srcs = ["rns_bgv_public_key.cc"],
    hdrs = ["rns_bgv_public_key.h"],
    deps = [
        ":coefficient_encoder",
        ":rns_bgv_ciphertext",
        ":rns_error_params",
        ":rns_modulus",
        ":rns_polynomial",
        ":rns_public_key",
        ":rns_secret_key",
        "//shell_encryption:integral_types",
        "//shell_encryption:montgomery",
        "//shell_encryption:statusor_fork",
        "//shell_encryption/prng",
        "@com_google_absl//absl/numeric:int128",
        "@com_google_absl//absl/status",
        "@com_google_absl//absl/status:statusor",
        "@com_google_absl//absl/strings",
        "@com_google_absl//absl/types:span",
    ],
)

cc_test(
    name = "rns_bgv_public_key_test",
    size = "small",
    srcs = ["rns_bgv_public_key_test.cc"],
    deps = [
        ":coefficient_encoder",
        ":rns_bgv_ciphertext",
        ":rns_bgv_public_key",
        ":rns_context",
        ":rns_error_params",
        ":rns_modulus",
        ":rns_polynomial",
        ":rns_secret_key",
        "//shell_encryption:integral_types",
        "//shell_encryption/prng:single_thread_hkdf_prng",
        "//shell_encryption/rns/testing:parameters",
        "//shell_encryption/rns/testing:testing_utils",
        "//shell_encryption/testing:matchers",
        "//shell_encryption/testing:parameters",
        "//shell_encryption/testing:status_testing",
        "@com_github_google_googletest//:gtest_main",
        "@com_google_absl//absl/log:check",
        "@com_google_absl//absl/random",
        "@com_google_absl//absl/status",
        "@com_google_absl//absl/strings",
        "@com_google_absl//absl/types:span",
    ],
)

# BFV variant of public key.
cc_library(
    name = "rns_bfv_public_key",
    srcs = ["rns_bfv_public_key.cc"],
    hdrs = ["rns_bfv_public_key.h"],
    deps = [
        ":coefficient_encoder",
        ":rns_bfv_ciphertext",
        ":rns_error_params",
        ":rns_modulus",
        ":rns_polynomial",
        ":rns_public_key",
        ":rns_secret_key",
        "//shell_encryption:integral_types",
        "//shell_encryption:montgomery",
        "//shell_encryption:statusor_fork",
        "//shell_encryption/prng",
        "@com_google_absl//absl/numeric:int128",
        "@com_google_absl//absl/status",
        "@com_google_absl//absl/status:statusor",
        "@com_google_absl//absl/strings",
        "@com_google_absl//absl/types:span",
    ],
)

cc_test(
    name = "rns_bfv_public_key_test",
    size = "small",
    srcs = ["rns_bfv_public_key_test.cc"],
    deps = [
        ":coefficient_encoder",
        ":rns_bfv_ciphertext",
        ":rns_bfv_public_key",
        ":rns_context",
        ":rns_error_params",
        ":rns_modulus",
        ":rns_polynomial",
        ":rns_secret_key",
        "//shell_encryption:integral_types",
        "//shell_encryption/prng:single_thread_hkdf_prng",
        "//shell_encryption/rns/testing:parameters",
        "//shell_encryption/rns/testing:testing_utils",
        "//shell_encryption/testing:matchers",
        "//shell_encryption/testing:parameters",
        "//shell_encryption/testing:status_testing",
        "@com_github_google_googletest//:gtest_main",
        "@com_google_absl//absl/log:check",
        "@com_google_absl//absl/random",
        "@com_google_absl//absl/status",
        "@com_google_absl//absl/strings",
        "@com_google_absl//absl/types:span",
    ],
)

# RNS gadget.
cc_library(
    name = "rns_gadget",
    srcs = ["rns_gadget.cc"],
    hdrs = ["rns_gadget.h"],
    deps = [
        ":rns_integer",
        ":rns_modulus",
        ":rns_polynomial",
        "//shell_encryption:gadget",
        "//shell_encryption:integral_types",
        "//shell_encryption:modulus_conversion",
        "//shell_encryption:montgomery",
        "//shell_encryption:statusor_fork",
        "@com_google_absl//absl/algorithm:container",
        "@com_google_absl//absl/numeric:int128",
        "@com_google_absl//absl/status",
        "@com_google_absl//absl/status:statusor",
        "@com_google_absl//absl/strings",
        "@com_google_absl//absl/types:span",
    ],
)

cc_test(
    name = "rns_gadget_test",
    srcs = ["rns_gadget_test.cc"],
    deps = [
        ":rns_context",
        ":rns_gadget",
        ":rns_modulus",
        ":rns_polynomial",
        "//shell_encryption/rns/testing:parameters",
        "//shell_encryption/testing:matchers",
        "//shell_encryption/testing:parameters",
        "//shell_encryption/testing:status_testing",
        "//shell_encryption/testing:testing_prng",
        "@com_github_google_googletest//:gtest_main",
        "@com_google_absl//absl/status",
    ],
)

# Gadget-based Galois key.
cc_library(
    name = "rns_galois_key",
    srcs = ["rns_galois_key.cc"],
    hdrs = ["rns_galois_key.h"],
    deps = [
        ":error_distribution",
        ":lazy_rns_polynomial",
        ":rns_bfv_ciphertext",
        ":rns_bgv_ciphertext",
        ":rns_ciphertext",
        ":rns_gadget",
        ":rns_modulus",
        ":rns_polynomial",
        ":rns_secret_key",
        ":serialization_cc_proto",
        "//shell_encryption:integral_types",
        "//shell_encryption:montgomery",
        "//shell_encryption:statusor_fork",
        "//shell_encryption/prng",
        "//shell_encryption/prng:single_thread_chacha_prng",
        "//shell_encryption/prng:single_thread_hkdf_prng",
        "@com_google_absl//absl/log:check",
        "@com_google_absl//absl/numeric:int128",
        "@com_google_absl//absl/status",
        "@com_google_absl//absl/status:statusor",
        "@com_google_absl//absl/strings",
        "@com_google_absl//absl/types:span",
    ],
)

cc_test(
    name = "rns_galois_key_test",
    srcs = ["rns_galois_key_test.cc"],
    deps = [
        ":finite_field_encoder",
        ":rns_bfv_ciphertext",
        ":rns_bgv_ciphertext",
        ":rns_context",
        ":rns_error_params",
        ":rns_gadget",
        ":rns_galois_key",
        ":rns_modulus",
        ":rns_polynomial",
        ":rns_secret_key",
        "//shell_encryption/rns/testing:parameters",
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

# Gadget-based relinearization key.
cc_library(
    name = "rns_relinearization_key",
    srcs = ["rns_relinearization_key.cc"],
    hdrs = ["rns_relinearization_key.h"],
    deps = [
        ":error_distribution",
        ":rns_bfv_ciphertext",
        ":rns_bgv_ciphertext",
        ":rns_ciphertext",
        ":rns_gadget",
        ":rns_modulus",
        ":rns_polynomial",
        ":rns_secret_key",
        ":serialization_cc_proto",
        "//shell_encryption:integral_types",
        "//shell_encryption:montgomery",
        "//shell_encryption:statusor_fork",
        "//shell_encryption/prng",
        "//shell_encryption/prng:single_thread_chacha_prng",
        "//shell_encryption/prng:single_thread_hkdf_prng",
        "@com_google_absl//absl/log:check",
        "@com_google_absl//absl/numeric:int128",
        "@com_google_absl//absl/status",
        "@com_google_absl//absl/status:statusor",
        "@com_google_absl//absl/strings",
        "@com_google_absl//absl/types:span",
    ],
)

cc_test(
    name = "rns_relinearization_key_test",
    srcs = ["rns_relinearization_key_test.cc"],
    deps = [
        ":finite_field_encoder",
        ":rns_bfv_ciphertext",
        ":rns_bgv_ciphertext",
        ":rns_context",
        ":rns_error_params",
        ":rns_gadget",
        ":rns_modulus",
        ":rns_polynomial",
        ":rns_relinearization_key",
        ":rns_secret_key",
        "//shell_encryption/rns/testing:parameters",
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

# Highway-based polynomial multiplication.
cc_library(
    name = "rns_polynomial_hwy",
    srcs = ["rns_polynomial_hwy.cc"],
    hdrs = ["rns_polynomial_hwy.h"],
    deps = [
        "//shell_encryption:integral_types",
        "//shell_encryption:montgomery",
        "@com_github_google_highway//:hwy",
        "@com_google_absl//absl/base:core_headers",
        "@com_google_absl//absl/numeric:int128",
        "@com_google_absl//absl/status",
        "@com_google_absl//absl/status:statusor",
        "@com_google_absl//absl/types:span",
    ],
)

# Lazy RNS polynomial.
cc_library(
    name = "lazy_rns_polynomial",
    srcs = ["lazy_rns_polynomial.cc"],
    hdrs = ["lazy_rns_polynomial.h"],
    deps = [
        ":rns_modulus",
        ":rns_polynomial",
        ":rns_polynomial_hwy",
        "//shell_encryption:integral_types",
        "//shell_encryption:montgomery",
        "//shell_encryption:statusor_fork",
        "@com_github_google_highway//:hwy",
        "@com_google_absl//absl/numeric:int128",
        "@com_google_absl//absl/status",
        "@com_google_absl//absl/status:statusor",
        "@com_google_absl//absl/strings",
        "@com_google_absl//absl/types:span",
    ],
)

cc_test(
    name = "lazy_rns_polynomial_test",
    srcs = ["lazy_rns_polynomial_test.cc"],
    deps = [
        ":lazy_rns_polynomial",
        ":rns_context",
        ":rns_modulus",
        ":rns_polynomial",
        "//shell_encryption:integral_types",
        "//shell_encryption:montgomery",
        "//shell_encryption:polynomial",
        "//shell_encryption:statusor_fork",
        "//shell_encryption/rns/testing:parameters",
        "//shell_encryption/rns/testing:testing_utils",
        "//shell_encryption/testing:matchers",
        "//shell_encryption/testing:parameters",
        "//shell_encryption/testing:status_testing",
        "//shell_encryption/testing:testing_prng",
        "@com_github_google_googletest//:gtest_main",
        "@com_google_absl//absl/log:check",
        "@com_google_absl//absl/status",
        "@com_google_absl//absl/strings",
    ],
)
