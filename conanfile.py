from conans import ConanFile, CMake


class ReaDDyTests(ConanFile):
    options = {}
    name = "ReaDDyTests"
    version = "0.1"
    requires = (
        "catch2/3.0.1"
    )
    generators = "cmake", "gcc", "txt", "cmake_find_package"

    def build(self):
        cmake = CMake(self)
        if self.settings.os == 'Macos':
            cmake.definitions["CONAN_CMAKE_CXX_STANDARD"] = "17"
        cmake.build()
