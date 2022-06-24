from conans import ConanFile


class ReaDDyTests(ConanFile):
    options = {}
    name = "ReaDDyTests"
    version = "0.1"
    requires = (
        "catch2/3.0.1"
    )
    generators = "cmake", "gcc", "txt", "cmake_find_package"
