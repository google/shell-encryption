load("@bazel_tools//tools/build_defs/repo:git.bzl", "git_repository")
load("@bazel_tools//tools/build_defs/repo:http.bzl", "http_archive")

# rules_cc defines rules for generating C++ code from Protocol Buffers.
http_archive(
    name = "rules_cc",
    sha256 = "35f2fb4ea0b3e61ad64a369de284e4fbbdcdba71836a5555abb5e194cf119509",
    strip_prefix = "rules_cc-624b5d59dfb45672d4239422fa1e3de1822ee110",
    urls = [
        "https://mirror.bazel.build/github.com/bazelbuild/rules_cc/archive/624b5d59dfb45672d4239422fa1e3de1822ee110.tar.gz",
        "https://github.com/bazelbuild/rules_cc/archive/624b5d59dfb45672d4239422fa1e3de1822ee110.tar.gz",
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
git_repository(
    name = "com_github_google_googletest",
    commit = "703bd9caab50b139428cea1aaff9974ebee5742e",  # tag = "release-1.10.0"
    remote = "https://github.com/google/googletest.git",
    shallow_since = "1570114335 -0400",
)

# abseil-cpp
git_repository(
    name = "com_github_abseil_abseil-cpp",
    commit = "df3ea785d8c30a9503321a3d35ee7d35808f190d",  # tag = "20200225.1"
    remote = "https://github.com/abseil/abseil-cpp.git",
    shallow_since = "1583355457 -0500",
)

# BoringSSL
git_repository(
    name = "boringssl",
    commit = "67ffb9606462a1897d3a5edf5c06d329878ba600",  # https://boringssl.googlesource.com/boringssl/+/refs/heads/master-with-bazel
    remote = "https://boringssl.googlesource.com/boringssl",
    shallow_since = "1585767053 +0000"
)

# Logging
git_repository(
    name = "com_github_google_glog",
    commit = "96a2f23dca4cc7180821ca5f32e526314395d26a",  # tag = "v0.4.0"
    remote = "https://github.com/google/glog.git",
    shallow_since = "1553223106 +0900",
)

# gflags, needed for glog
git_repository(
    name = "com_github_gflags_gflags",
    commit = "e171aa2d15ed9eb17054558e0b3a6a413bb01067",  # tag = "v2.2.2"
    remote = "https://github.com/gflags/gflags.git",
    shallow_since = "1541971260 +0000",
)
