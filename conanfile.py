from conans import ConanFile, CMake


class ReaDDyTests(ConanFile):
    options = {}
    name = "ReaDDyTests"
    version = "0.1"
    requires = (
        "catch2/3.0.1"
    )
    generators = "cmake", "gcc", "txt", "cmake_find_package"

    def configure(self):
        self.settings.compiler.cppstd = 17

    def build(self):
        cmake = CMake(self)
        cmake.definitions["CONAN_CMAKE_CXX_STANDARD"] = "17"
        cmake.configure()
        cmake.build()
        cmake.install()
