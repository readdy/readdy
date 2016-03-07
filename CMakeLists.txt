#-----------------------------------
# ReaDDy2
#-----------------------------------
CMAKE_MINIMUM_REQUIRED(VERSION 3.2.0)

PROJECT(ReaDDy2_Main C CXX)

# for custom cmake modules
SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules")

## ensure out of tree build
INCLUDE(MacroOutOfSourceBuild)
MACRO_ENSURE_OUT_OF_SOURCE_BUILD("${PROJECT_NAME} requires an out source build.")

# SET(CMAKE_VERBOSE_MAKEFILE ON)
SET(CMAKE_COLOR_MAKEFILE ON)

# we use c++11
SET(CMAKE_CXX_STANDARD 11)
SET(CMAKE_CXX_STANDARD_REQUIRED ON)

# version
SET(READDY_MAJOR_VERSION 0)
SET(READDY_MINOR_VERSION 0)
SET(READDY_PATCH_VERSION 1)
SET(READDY_VERSION "${READDY_MAJOR_VERSION}.${READDY_MINOR_VERSION}.${READDY_PATCH_VERSION}")

# should run tests
SET(READDY_RUN_TESTS ON CACHE BOOL "Should run tests or not")
SET(READDY_RUN_TESTS_MEMORY_CHECK OFF CACHE BOOL "Determines if valgrind test targets should be added.")

# Compile in debug mode by default
IF (NOT CMAKE_BUILD_TYPE)
    MESSAGE(STATUS "Setting build type to 'DEBUG' as none was specified.")
    SET(CMAKE_BUILD_TYPE DEBUG)
ELSE (NOT CMAKE_BUILD_TYPE)
    MESSAGE(STATUS "Using specified '${CMAKE_BUILD_TYPE}' build type.")
ENDIF (NOT CMAKE_BUILD_TYPE)

# decide if api wrappers are being generated
SET(READDY_BUILD_PYTHON_WRAPPER ON CACHE BOOL "Build Python wrappers")

# includes
SET(READDY_GLOBAL_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/include")
LIST(APPEND READDY_HEADERS "${PROJECT_SOURCE_DIR}/include")

#
# readdy2 core
#

# include directories from least dependent to most dependent
# i.e.: A depends on B -> add subdirectory of B first
ADD_SUBDIRECTORY(readdy2/main/io)
ADD_SUBDIRECTORY(readdy2/main/plugin)
ADD_SUBDIRECTORY(readdy2/main)

#
# readdy2 kernels
#
ADD_SUBDIRECTORY(kernels)

#
# readdy2 wrappers
#
IF (READDY_BUILD_PYTHON_WRAPPER)
    ADD_SUBDIRECTORY(wrappers)
ENDIF ()

#
# tests
#
IF (EXISTS "${CMAKE_SOURCE_DIR}/libraries/googletest/CMakeLists.txt" AND READDY_RUN_TESTS)
    MESSAGE(STATUS "The submodule googletest was initialized, thus enable testing")

    INCLUDE(ExternalProject)

    SET(GOOGLETEST_ROOT ${CMAKE_SOURCE_DIR}/libraries/googletest)
    SET(GOOGLETEST_INCLUDE ${GOOGLETEST_ROOT}/googletest/include)
    SET(GOOGLETEST_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/libraries/googletest)
    SET(GOOGLETEST_LIBRARY ${GOOGLETEST_BINARY_DIR}/googletest)

    ExternalProject_Add(
            googletest
            DOWNLOAD_COMMAND ""
            SOURCE_DIR ${GOOGLETEST_ROOT}
            BINARY_DIR ${GOOGLETEST_BINARY_DIR}
            CMAKE_CACHE_ARGS "-DBUILD_GTEST:bool=ON" "-DBUILD_GMOCK:bool=OFF"
            "-Dgtest_force_shared_crt:bool=ON"
            INSTALL_COMMAND ""
    )
    ADD_DEPENDENCIES(ReaDDy2 googletest)
    ENABLE_TESTING()
    INCLUDE_DIRECTORIES(${GOOGLETEST_INCLUDE})
    ADD_SUBDIRECTORY(readdy2/test)

ELSE ()
    MESSAGE(STATUS "The submodule googletest was not initialized or READDY_RUN_TESTS was OFF, thus no testing.")
ENDIF ()