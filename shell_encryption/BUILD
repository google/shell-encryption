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
load("@rules_license//rules:license.bzl", "license")
load("@rules_proto//proto:defs.bzl", "proto_library")

package(default_visibility = ["//visibility:public"])

licenses(["notice"])

exports_files(["LICENSE"])

# Protos.
proto_library(
    name = "serialization_proto",
    srcs = ["serialization.proto"],
)

cc_proto_library(
    name = "serialization_cc_proto",
    deps = [":serialization_proto"],
)

# Context.

cc_library(
    name = "context",
    hdrs = ["context.h"],
    deps = [
        ":error_params",
        ":ntt_parameters",
        ":statusor_fork",
        "@com_google_absl//absl/memory",
        "@com_google_absl//absl/status",
    ],
)

cc_test(
    name = "context_test",
    size = "small",
    srcs = ["context_test.cc"],
    deps = [
        ":constants",
        ":context",
        ":integral_types",
        ":montgomery",
        ":statusor_fork",
        "//shell_encryption/testing:parameters",
        "//shell_encryption/testing:status_testing",
        "@com_github_google_googletest//:gtest_main",
        "@com_google_absl//absl/numeric:int128",
    ],
)

# Integral types for integers.

cc_library(
    name = "integral_types",
    hdrs = ["integral_types.h"],
    deps = [
        "@com_google_absl//absl/numeric:int128",
    ],
)

# Constants.

cc_library(
    name = "constants",
    hdrs = ["constants.h"],
    deps = [
        ":integral_types",
        "@com_google_absl//absl/numeric:int128",
    ],
)

# Utilities.

cc_library(
    name = "bits_util",
    hdrs = ["bits_util.h"],
    deps = [
        ":integral_types",
        "@com_google_absl//absl/numeric:int128",
    ],
)

cc_test(
    name = "bits_util_test",
    size = "small",
    srcs = ["bits_util_test.cc"],
    deps = [
        ":bits_util",
        "@com_github_google_googletest//:gtest_main",
        "@com_google_absl//absl/numeric:int128",
    ],
)

cc_library(
    name = "sample_error",
    hdrs = ["sample_error.h"],
    deps = [
        ":bits_util",
        ":constants",
        ":error_params",
        ":statusor_fork",
        "//shell_encryption/prng",
        "@com_google_absl//absl/strings",
    ],
)

cc_test(
    name = "sample_error_test",
    size = "small",
    srcs = [
        "sample_error_test.cc",
    ],
    deps = [
        ":context",
        ":montgomery",
        ":sample_error",
        ":symmetric_encryption",
        "//shell_encryption/testing:parameters",
        "//shell_encryption/testing:status_is_fork",
        "//shell_encryption/testing:status_testing",
        "//shell_encryption/testing:testing_prng",
        "@com_github_google_googletest//:gtest_main",
    ],
)

cc_library(
    name = "statusor_fork",
    hdrs = [
        "status_macros.h",
        "statusor.h",
    ],
    deps = [
        "@com_github_google_glog//:glog",
        "@com_google_absl//absl/base:core_headers",
        "@com_google_absl//absl/status",
        "@com_google_absl//absl/status:statusor",
        "@com_google_absl//absl/types:optional",
    ],
)

cc_test(
    name = "status_macros_test",
    size = "small",
    srcs = ["status_macros_test.cc"],
    deps = [
        ":statusor_fork",
        "@com_github_google_googletest//:gtest_main",
        "@com_google_absl//absl/status",
    ],
)

# Montgomery integers.

cc_library(
    name = "montgomery",
    srcs = ["montgomery.cc"],
    hdrs = ["montgomery.h"],
    deps = [
        ":bits_util",
        ":constants",
        ":int256",
        ":serialization_cc_proto",
        ":statusor_fork",
        ":transcription",
        "//shell_encryption/prng",
        "@com_github_google_glog//:glog",
        "@com_google_absl//absl/numeric:int128",
        "@com_google_absl//absl/status:statusor",
        "@com_google_absl//absl/strings",
    ],
)

cc_test(
    name = "montgomery_test",
    size = "small",
    srcs = [
        "montgomery_test.cc",
    ],
    deps = [
        ":constants",
        ":montgomery",
        ":serialization_cc_proto",
        ":statusor_fork",
        "//shell_encryption/testing:parameters",
        "//shell_encryption/testing:status_is_fork",
        "//shell_encryption/testing:status_testing",
        "//shell_encryption/testing:testing_prng",
        "@com_github_google_googletest//:gtest_main",
        "@com_google_absl//absl/numeric:int128",
    ],
)

# NTT parameters.

cc_library(
    name = "ntt_parameters",
    srcs = ["ntt_parameters.cc"],
    hdrs = ["ntt_parameters.h"],
    deps = [
        ":constants",
        ":statusor_fork",
        "@com_google_absl//absl/memory",
        "@com_google_absl//absl/strings",
    ],
)

cc_test(
    name = "ntt_parameters_test",
    size = "small",
    srcs = [
        "ntt_parameters_test.cc",
    ],
    deps = [
        ":constants",
        ":montgomery",
        ":ntt_parameters",
        ":statusor_fork",
        "//shell_encryption/testing:parameters",
        "//shell_encryption/testing:status_is_fork",
        "//shell_encryption/testing:status_testing",
        "@com_github_google_googletest//:gtest_main",
        "@com_google_absl//absl/numeric:int128",
    ],
)

# Discrete Fourier transformations.

cc_library(
    name = "dft_transformations",
    srcs = ["dft_transformations.cc"],
    hdrs = ["dft_transformations.h"],
    deps = [
        ":ntt_parameters",
        ":statusor_fork",
        "@com_google_absl//absl/status",
    ],
)

cc_test(
    name = "dft_transformations_test",
    size = "small",
    srcs = [
        "dft_transformations_test.cc",
    ],
    deps = [
        ":context",
        ":dft_transformations",
        ":statusor_fork",
        "//shell_encryption/testing:parameters",
        "//shell_encryption/testing:status_is_fork",
        "//shell_encryption/testing:status_testing",
        "//shell_encryption/testing:testing_prng",
        "@com_github_google_googletest//:gtest_main",
        "@com_google_absl//absl/status",
        "@com_google_absl//absl/strings",
    ],
)

cc_library(
    name = "dft_transformations_hwy",
    srcs = ["dft_transformations_hwy.cc"],
    hdrs = [
        "dft_transformations_hwy.h",
    ],
    deps = [
        ":dft_transformations",
        ":montgomery",
        ":ntt_parameters",
        ":statusor_fork",
        "@com_github_google_highway//:hwy",
        "@com_google_absl//absl/numeric:int128",
        "@com_google_absl//absl/status",
    ],
)

cc_test(
    name = "dft_transformations_hwy_test",
    size = "small",
    srcs = [
        "dft_transformations_hwy_test.cc",
    ],
    deps = [
        ":context",
        ":dft_transformations",
        ":dft_transformations_hwy",
        ":statusor_fork",
        "//shell_encryption/testing:parameters",
        "//shell_encryption/testing:status_is_fork",
        "//shell_encryption/testing:status_testing",
        "//shell_encryption/testing:testing_prng",
        "@com_github_google_googletest//:gtest_main",
        "@com_github_google_highway//:hwy",
        "@com_google_absl//absl/status",
        "@com_google_absl//absl/strings",
    ],
)

# NTT polynomial.

cc_library(
    name = "polynomial",
    hdrs = ["polynomial.h"],
    deps = [
        ":constants",
        ":dft_transformations",
        ":dft_transformations_hwy",
        ":ntt_parameters",
        ":serialization_cc_proto",
        ":statusor_fork",
        "//shell_encryption/opt:constant_polynomial",
        "//shell_encryption/prng",
        "@com_google_absl//absl/status",
        "@com_google_absl//absl/strings",
        "@com_google_absl//absl/strings:str_format",
    ],
)

cc_test(
    name = "polynomial_test",
    size = "small",
    timeout = "moderate",
    srcs = [
        "polynomial_test.cc",
    ],
    deps = [
        ":constants",
        ":montgomery",
        ":ntt_parameters",
        ":polynomial",
        ":serialization_cc_proto",
        ":statusor_fork",
        "//shell_encryption/prng:integral_prng_testing_types",
        "//shell_encryption/testing:coefficient_polynomial",
        "//shell_encryption/testing:status_is_fork",
        "//shell_encryption/testing:status_testing",
        "//shell_encryption/testing:testing_prng",
        "@com_github_google_googletest//:gtest_main",
    ],
)

# Error parameters.

cc_library(
    name = "error_params",
    hdrs = ["error_params.h"],
    deps = [
        ":montgomery",
        ":ntt_parameters",
        ":statusor_fork",
    ],
)

cc_test(
    name = "error_params_test",
    size = "small",
    srcs = [
        "error_params_test.cc",
    ],
    deps = [
        ":constants",
        ":context",
        ":error_params",
        ":montgomery",
        ":ntt_parameters",
        ":statusor_fork",
        ":symmetric_encryption",
        "//shell_encryption/testing:parameters",
        "//shell_encryption/testing:status_is_fork",
        "//shell_encryption/testing:status_testing",
        "//shell_encryption/testing:testing_prng",
        "//shell_encryption/testing:testing_utils",
        "@com_github_google_googletest//:gtest_main",
    ],
)

# Symmetric Encryption.

cc_library(
    name = "symmetric_encryption",
    hdrs = ["symmetric_encryption.h"],
    deps = [
        ":error_params",
        ":polynomial",
        ":sample_error",
        ":serialization_cc_proto",
        ":statusor_fork",
        "//shell_encryption/opt:lazy_polynomial",
        "//shell_encryption/prng",
        "@com_google_absl//absl/status",
    ],
)

cc_test(
    name = "symmetric_encryption_test",
    size = "medium",
    srcs = [
        "symmetric_encryption_test.cc",
    ],
    deps = [
        ":constants",
        ":context",
        ":montgomery",
        ":ntt_parameters",
        ":polynomial",
        ":serialization_cc_proto",
        ":statusor_fork",
        ":symmetric_encryption",
        "//shell_encryption/testing:parameters",
        "//shell_encryption/testing:status_is_fork",
        "//shell_encryption/testing:status_testing",
        "//shell_encryption/testing:testing_prng",
        "//shell_encryption/testing:testing_utils",
        "@com_github_google_googletest//:gtest_main",
        "@com_google_absl//absl/status",
    ],
)

cc_library(
    name = "symmetric_encryption_with_prng",
    hdrs = ["symmetric_encryption_with_prng.h"],
    deps = [
        ":polynomial",
        ":statusor_fork",
        ":symmetric_encryption",
        "//shell_encryption/prng",
    ],
)

cc_test(
    name = "symmetric_encryption_with_prng_test",
    size = "small",
    srcs = [
        "symmetric_encryption_with_prng_test.cc",
    ],
    deps = [
        ":context",
        ":montgomery",
        ":ntt_parameters",
        ":polynomial",
        ":statusor_fork",
        ":symmetric_encryption_with_prng",
        "//shell_encryption/testing:parameters",
        "//shell_encryption/testing:status_is_fork",
        "//shell_encryption/testing:status_testing",
        "//shell_encryption/testing:testing_utils",
        "@com_github_google_googletest//:gtest_main",
    ],
)

# Public Key Encryption

cc_library(
    name = "public_key_encryption",
    srcs = ["public_key_encryption.cc"],
    hdrs = ["public_key_encryption.h"],
    deps = [
        ":error_params",
        ":montgomery",
        ":polynomial",
        ":sample_error",
        ":serialization_cc_proto",
        ":statusor_fork",
        ":symmetric_encryption",
        "//shell_encryption/prng",
        "//shell_encryption/prng:single_thread_chacha_prng",
        "//shell_encryption/prng:single_thread_hkdf_prng",
        "@com_google_absl//absl/numeric:int128",
        "@com_google_absl//absl/status",
    ],
)

cc_test(
    name = "public_key_encryption_test",
    size = "medium",
    srcs = [
        "public_key_encryption_test.cc",
    ],
    deps = [
        ":constants",
        ":context",
        ":montgomery",
        ":ntt_parameters",
        ":polynomial",
        ":public_key_encryption",
        ":serialization_cc_proto",
        ":statusor_fork",
        ":symmetric_encryption",
        "//shell_encryption/testing:parameters",
        "//shell_encryption/testing:status_is_fork",
        "//shell_encryption/testing:status_testing",
        "//shell_encryption/testing:testing_prng",
        "//shell_encryption/testing:testing_utils",
        "@com_github_google_googletest//:gtest_main",
        "@com_google_absl//absl/status",
    ],
)

cc_library(
    name = "gadget",
    hdrs = ["gadget.h"],
    deps = [
        ":statusor_fork",
        "@com_google_absl//absl/status",
    ],
)

# Relinearization Key

cc_library(
    name = "relinearization_key",
    srcs = ["relinearization_key.cc"],
    hdrs = ["relinearization_key.h"],
    deps = [
        ":bits_util",
        ":gadget",
        ":montgomery",
        ":sample_error",
        ":serialization_cc_proto",
        ":statusor_fork",
        ":symmetric_encryption",
        ":symmetric_encryption_with_prng",
        "//shell_encryption/prng",
        "//shell_encryption/prng:single_thread_chacha_prng",
        "//shell_encryption/prng:single_thread_hkdf_prng",
        "@com_google_absl//absl/memory",
        "@com_google_absl//absl/numeric:int128",
        "@com_google_absl//absl/status",
    ],
)

cc_test(
    name = "relinearization_key_test",
    size = "small",
    srcs = ["relinearization_key_test.cc"],
    deps = [
        ":constants",
        ":gadget",
        ":montgomery",
        ":ntt_parameters",
        ":polynomial",
        ":relinearization_key",
        ":serialization_cc_proto",
        ":statusor_fork",
        ":symmetric_encryption",
        "//shell_encryption/prng:chacha_prng",
        "//shell_encryption/prng:hkdf_prng",
        "//shell_encryption/prng:single_thread_chacha_prng",
        "//shell_encryption/prng:single_thread_hkdf_prng",
        "//shell_encryption/testing:status_is_fork",
        "//shell_encryption/testing:status_testing",
        "//shell_encryption/testing:testing_prng",
        "//shell_encryption/testing:testing_utils",
        "@com_github_google_googletest//:gtest_main",
        "@com_google_absl//absl/status",
    ],
)

# Galois Key.

cc_library(
    name = "galois_key",
    hdrs = ["galois_key.h"],
    deps = [
        ":relinearization_key",
        ":statusor_fork",
        "@com_google_absl//absl/strings",
    ],
)

cc_test(
    name = "galois_key_test",
    size = "small",
    srcs = ["galois_key_test.cc"],
    deps = [
        ":constants",
        ":galois_key",
        ":montgomery",
        ":ntt_parameters",
        ":polynomial",
        ":statusor_fork",
        ":symmetric_encryption",
        "//shell_encryption/testing:status_is_fork",
        "//shell_encryption/testing:status_testing",
        "//shell_encryption/testing:testing_prng",
        "//shell_encryption/testing:testing_utils",
        "@com_github_google_googletest//:gtest_main",
    ],
)

# Utilities to convert between different modulus representations.
cc_library(
    name = "modulus_conversion",
    hdrs = ["modulus_conversion.h"],
    deps = [
        ":montgomery",
        ":ntt_parameters",
        ":polynomial",
        ":statusor_fork",
        "@com_google_absl//absl/types:span",
    ],
)

# Relinearization key using the auxiliary modulus technique
cc_library(
    name = "aux_relinearization_key",
    srcs = ["aux_relinearization_key.cc"],
    hdrs = ["aux_relinearization_key.h"],
    deps = [
        ":modulus_conversion",
        ":montgomery",
        ":ntt_parameters",
        ":polynomial",
        ":serialization_cc_proto",
        ":statusor_fork",
        ":symmetric_encryption",
        "//shell_encryption/prng",
        "//shell_encryption/prng:single_thread_chacha_prng",
        "//shell_encryption/prng:single_thread_hkdf_prng",
        "@com_google_absl//absl/memory",
        "@com_google_absl//absl/numeric:int128",
        "@com_google_absl//absl/status",
    ],
)

cc_test(
    name = "aux_relinearization_key_test",
    size = "small",
    srcs = ["aux_relinearization_key_test.cc"],
    deps = [
        ":aux_relinearization_key",
        ":constants",
        ":integral_types",
        ":montgomery",
        ":ntt_parameters",
        ":polynomial",
        ":serialization_cc_proto",
        ":statusor_fork",
        ":symmetric_encryption",
        "//shell_encryption/prng:chacha_prng",
        "//shell_encryption/prng:hkdf_prng",
        "//shell_encryption/prng:single_thread_chacha_prng",
        "//shell_encryption/prng:single_thread_hkdf_prng",
        "//shell_encryption/testing:status_is_fork",
        "//shell_encryption/testing:status_testing",
        "//shell_encryption/testing:testing_prng",
        "//shell_encryption/testing:testing_utils",
        "@com_github_google_googletest//:gtest_main",
        "@com_google_absl//absl/status",
    ],
)

# Galois key using the auxiliary modulus technique
cc_library(
    name = "aux_galois_key",
    hdrs = ["aux_galois_key.h"],
    deps = [
        ":aux_relinearization_key",
        ":montgomery",
        ":ntt_parameters",
        ":polynomial",
        ":serialization_cc_proto",
        ":statusor_fork",
        ":symmetric_encryption",
        "//shell_encryption/prng",
    ],
)

cc_test(
    name = "aux_galois_key_test",
    size = "small",
    srcs = ["aux_galois_key_test.cc"],
    deps = [
        ":aux_galois_key",
        ":constants",
        ":montgomery",
        ":ntt_parameters",
        ":polynomial",
        ":serialization_cc_proto",
        ":statusor_fork",
        ":symmetric_encryption",
        "//shell_encryption/testing:status_is_fork",
        "//shell_encryption/testing:status_testing",
        "//shell_encryption/testing:testing_utils",
        "@com_github_google_googletest//:gtest_main",
        "@com_google_absl//absl/status",
    ],
)

# int256 Type

cc_library(
    name = "int256",
    srcs = ["int256.cc"],
    hdrs = ["int256.h"],
    deps = [
        ":integral_types",
        "@com_github_google_glog//:glog",
        "@com_google_absl//absl/numeric:int128",
    ],
)

cc_test(
    name = "int256_test",
    size = "small",
    srcs = ["int256_test.cc"],
    deps = [
        ":int256",
        "@com_github_google_glog//:glog",
        "@com_github_google_googletest//:gtest_main",
        "@com_google_absl//absl/container:fixed_array",
        "@com_google_absl//absl/numeric:int128",
    ],
)

# Transcription

cc_library(
    name = "transcription",
    hdrs = ["transcription.h"],
    deps = [
        ":statusor_fork",
        "@com_google_absl//absl/strings",
    ],
)

cc_test(
    name = "transcription_test",
    size = "small",
    timeout = "moderate",
    srcs = ["transcription_test.cc"],
    deps = [
        ":integral_types",
        ":statusor_fork",
        ":transcription",
        "//shell_encryption/testing:status_is_fork",
        "//shell_encryption/testing:status_testing",
        "@com_github_google_googletest//:gtest_main",
        "@com_google_absl//absl/numeric:int128",
    ],
)

# Oblivious expand utilities

cc_library(
    name = "oblivious_expand",
    hdrs = ["oblivious_expand.h"],
    deps = [
        ":galois_key",
        ":integral_types",
        ":ntt_parameters",
        ":polynomial",
        ":statusor_fork",
        ":symmetric_encryption",
        "@com_google_absl//absl/memory",
        "@com_google_absl//absl/status",
    ],
)

cc_test(
    name = "oblivious_expand_test",
    size = "small",
    srcs = ["oblivious_expand_test.cc"],
    deps = [
        ":constants",
        ":galois_key",
        ":montgomery",
        ":ntt_parameters",
        ":oblivious_expand",
        ":symmetric_encryption",
        "//shell_encryption/testing:status_is_fork",
        "//shell_encryption/testing:status_testing",
        "//shell_encryption/testing:testing_utils",
        "@com_github_google_googletest//:gtest_main",
        "@com_google_absl//absl/memory",
        "@com_google_absl//absl/status",
    ],
)
