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
        cmake.definitions['CATCH_CONFIG_NO_CPP17_UNCAUGHT_EXCEPTIONS'] = ''
        if self.settings.os == 'Macos':
            cmake.definitions["CONAN_LIBCXX"] = ""
        cmake.build()
