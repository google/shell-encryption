load("@bazel_tools//tools/build_defs/repo:git.bzl", "git_repository")
load("@bazel_tools//tools/build_defs/repo:http.bzl", "http_archive")

# Tink
http_archive(
    name = "com_github_tink_crypto_tink_cc",
    urls = ["https://github.com/tink-crypto/tink-cc/archive/refs/tags/v2.0.0.zip"],
    strip_prefix = "tink-cc-2.0.0",
    sha256 = "103ddfce800e77f3b3b6b2c808a8611bc734b31ddb12fbcfd8bebc1b96a7e963",
)

# Load Tink dependencies.

load("@com_github_tink_crypto_tink_cc//:tink_cc_deps.bzl", "tink_cc_deps")

tink_cc_deps()

load("@com_github_tink_crypto_tink_cc//:tink_cc_deps_init.bzl", "tink_cc_deps_init")

tink_cc_deps_init()

# rules_cc defines rules for generating C++ code from Protocol Buffers.
http_archive(
    name = "rules_cc",
    sha256 = "2037875b9a4456dce4a79d112a8ae885bbc4aad968e6587dca6e64f3a0900cdf",
    strip_prefix = "rules_cc-0.0.9",
    urls = [
        "https://github.com/bazelbuild/rules_cc/releases/download/0.0.9/rules_cc-0.0.9.tar.gz",
    ],
)

load("@rules_cc//cc:repositories.bzl", "rules_cc_dependencies")
rules_cc_dependencies()

# rules_proto defines abstract rules for building Protocol Buffers.
# https://github.com/bazelbuild/rules_proto
http_archive(
    name = "rules_proto",
    sha256 = "602e7161d9195e50246177e7c55b2f39950a9cf7366f74ed5f22fd45750cd208",
    strip_prefix = "rules_proto-97d8af4dc474595af3900dd85cb3a29ad28cc313",
    urls = [
        "https://mirror.bazel.build/github.com/bazelbuild/rules_proto/archive/97d8af4dc474595af3900dd85cb3a29ad28cc313.tar.gz",
        "https://github.com/bazelbuild/rules_proto/archive/97d8af4dc474595af3900dd85cb3a29ad28cc313.tar.gz",
    ],
)
load("@rules_proto//proto:repositories.bzl", "rules_proto_dependencies", "rules_proto_toolchains")
rules_proto_dependencies()
rules_proto_toolchains()

# Install gtest.
http_archive(
  name = "com_github_google_googletest",
  urls = ["https://github.com/google/googletest/archive/5ab508a01f9eb089207ee87fd547d290da39d015.zip"],
  strip_prefix = "googletest-5ab508a01f9eb089207ee87fd547d290da39d015",
)

# abseil-cpp
http_archive(
    name = "com_google_absl",
    sha256 = "3c743204df78366ad2eaf236d6631d83f6bc928d1705dd0000b872e53b73dc6a",
    strip_prefix = "abseil-cpp-20240116.1",
    urls = [
         "https://github.com/abseil/abseil-cpp/releases/download/20240116.1/abseil-cpp-20240116.1.tar.gz",
    ],
)

# BoringSSL
git_repository(
    name = "boringssl",
    commit = "67ffb9606462a1897d3a5edf5c06d329878ba600",  # https://boringssl.googlesource.com/boringssl/+/refs/heads/master-with-bazel
    remote = "https://boringssl.googlesource.com/boringssl",
    shallow_since = "1585767053 +0000"
)

# Logging
http_archive(
    name = "com_github_google_glog",
    urls = ["https://github.com/google/glog/archive/96a2f23dca4cc7180821ca5f32e526314395d26a.zip"],
    strip_prefix = "glog-96a2f23dca4cc7180821ca5f32e526314395d26a",
    sha256 = "6281aa4eeecb9e932d7091f99872e7b26fa6aacece49c15ce5b14af2b7ec050f",
)

# gflags, needed for glog
http_archive(
    name = "com_github_gflags_gflags",
    urls = ["https://github.com/gflags/gflags/archive/v2.2.2.tar.gz"],
    sha256 = "34af2f15cf7367513b352bdcd2493ab14ce43692d2dcd9dfc499492966c64dcf",
    strip_prefix = "gflags-2.2.2",
)
