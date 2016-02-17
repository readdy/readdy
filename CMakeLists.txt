#-----------------------------------
# ReaDDy2
#-----------------------------------
PROJECT(ReaDDy2)
CMAKE_MINIMUM_REQUIRED(3.2)

# for custom cmake modules
SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules")

# we use c++11
SET(CMAKE_CXX_STANDARD 11)
SET(CMAKE_CXX_STANDARD_REQUIRED ON)

# name and version
SET(READDY_LIBRARY_NAME "readdy2")
SET(READDY_MAJOR_VERSION 0)
SET(READDY_MINOR_VERSION 0)
SET(READDY_PATCH_VERSION 1)
SET(READDY_VERSION "${READDY_MAJOR_VERSION}.${READDY_MINOR_VERSION}.${READDY_PATCH_VERSION}")

# should run tests
SET(READDY_RUN_TESTS ON CACHE BOOL "Should run tests or not")

# Compile in debug mode by default
IF(NOT CMAKE_BUILD_TYPE)
    MESSAGE(STATUS "Setting build type to 'DEBUG' as none was specified.")
    SET(CMAKE_BUILD_TYPE DEBUG) 
ELSE(NOT CMAKE_BUILD_TYPE)
    MESSAGE(STATUS "Using specified '${CMAKE_BUILD_TYPE}' build type.")
ENDIF(NOT CMAKE_BUILD_TYPE)


# target names
SET(SHARED_TARGET ${READDY_LIBRARY_NAME})
SET(STATIC_TARGET ${READDY_LIBRARY_NAME}_static)

# decide if api wrappers are being generated
SET(READDY_BUILD_PYTHON_WRAPPER ON CACHE BOOL "Build Python wrappers")
IF(READDY_BUILD_PYTHON_WRAPPER)
    ADD_SUBDIRECTORY(wrappers)
    # append sources
    SET(SOURCE_FILES ${SOURCE_FILES} wrappers/*.cpp)
ENDIF()

IF(EXISTS "${CMAKE_SOURCE_DIR}/libraries/googletest/CMakeLists.txt" AND READDY_RUN_TESTS)
    MESSAGE(STATUS "The submodule googletest was initialized, thus enable testing")
    SET(GOOGLETEST_ROOT ${CMAKE_SOURCE_DIR}/libraries/googletest)
    SET(GOOGLETEST_INCLUDE ${GOOGLETEST_ROOT}/libraries/googletest/include)
    SET(GOOGLETEST_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/libraries/googletest)
    SET(GOOGLETEST_LIBRARY ${GOOGLETEST_BINARY_DIR}/libraries/googletest)

    ExternalProject_Add(
        googletest
        DOWNLOAD_COMMAND ""
        SOURCE_DIR ${GOOGLETEST_ROOT}
        BINARY_DIR ${GOOGLETEST_BINARY_DIR}
        CMAKE_CACHE_ARGS "-DBUILD_GTEST:bool=ON" "-DBUILD_GMOCK:bool=OFF"
                         "-Dgtest_force_shared_crt:bool=ON"
        INSTALL_COMMAND ""
    )
    # TODO: ADD_DEPENDENCY(readdy googletest)
    ENABLE_TESTING()
    # TODO: ADD_SUBDIRECTORY (SUBDIRECTORIES) with tests
ELSE()
    MESSAGE(STATUS "The submodule googletest was not initialized or READDY_RUN_TESTS was OFF, thus no testing.")
ENDIF()
