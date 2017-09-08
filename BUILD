# Copyright 2017 Google Inc.
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

licenses(["notice"])  # Apache 2.0

exports_files(["LICENSE"])

# Protos.

proto_library(
    name = "serialization_proto",
    srcs = ["serialization.proto"],
)

cc_proto_library(
    name = "serialization",
    deps = [":serialization_proto"],
)

# Constants.

cc_library(
    name = "constants",
    hdrs = ["constants.h"],
)

# Utilities.

cc_library(
    name = "utils",
    hdrs = ["utils.h"],
)

# Montgomery integers.

cc_library(
    name = "montgomery",
    hdrs = ["montgomery.h"],
    deps = [
        ":serialization",
        ":utils",
    ],
)

cc_test(
    name = "montgomery_test",
    srcs = [
        "montgomery_test.cc",
    ],
    deps = [
        ":constants",
        ":montgomery",
        ":serialization",
        "@com_google_googletest//:gtest_main",
    ],
)

# NTT parameters.

cc_library(
    name = "ntt_parameters",
    srcs = ["ntt_parameters.cc"],
    hdrs = ["ntt_parameters.h"],
)

cc_test(
    name = "ntt_parameters_test",
    srcs = [
        "ntt_parameters_test.cc",
    ],
    deps = [
        ":constants",
        ":montgomery",
        ":ntt_parameters",
        "@com_google_googletest//:gtest_main",
    ],
)

# NTT polynomial.

cc_library(
    name = "ntt_polynomial",
    hdrs = ["ntt_polynomial.h"],
    deps = [
        ":ntt_parameters",
        ":serialization",
    ],
)

cc_test(
    name = "ntt_polynomial_test",
    srcs = [
        "ntt_polynomial_test.cc",
    ],
    deps = [
        ":constants",
        ":montgomery",
        ":ntt_parameters",
        ":ntt_polynomial",
        ":serialization",
        "//testing:polynomial_test_util",
        "@com_google_googletest//:gtest_main",
    ],
)

# Encryption.

cc_library(
    name = "symmetric_encryption",
    hdrs = ["symmetric_encryption.h"],
    deps = [
        ":ntt_polynomial",
        ":serialization",
        ":utils",
    ],
)

cc_test(
    name = "symmetric_encryption_test",
    srcs = [
        "symmetric_encryption_test.cc",
    ],
    deps = [
        ":constants",
        ":montgomery",
        ":ntt_parameters",
        ":ntt_polynomial",
        ":serialization",
        ":symmetric_encryption",
        ":utils",
        "//testing:utils",
        "@com_google_googletest//:gtest_main",
    ],
)
