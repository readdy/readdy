from conans import ConanFile, CMake


class ReaDDyTests(ConanFile):
    options = {}
    name = "ReaDDyTests"
    version = "0.1"
    requires = (
        "nlohmann_json/3.10.3",
        "spdlog/1.10.0",
        "fmt/8.1.1"
    )
    generators = "cmake", "gcc", "txt", "cmake_find_package"

    def configure(self):
        self.options["spdlog"].header_only = True
        self.options["fmt"].header_only = True
        super().configure()
