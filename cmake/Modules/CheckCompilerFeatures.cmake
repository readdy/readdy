include(CheckCXXSourceCompiles)

CHECK_CXX_SOURCE_COMPILES("
    #include <memory>
    int main() {
        std::make_unique<int>();
    }" COMPILER_SUPPORTS_MAKE_UNIQUE)