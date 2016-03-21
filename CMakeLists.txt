#-----------------------------------
# ReaDDy2
#-----------------------------------
CMAKE_MINIMUM_REQUIRED(VERSION 3.2.0)

PROJECT(ReaDDy2_Main C CXX)
SET (CMAKE_CXX_STANDARD 11)

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

# build as shared lib
SET(BUILD_SHARED_LIBS ON)
# build single libraries or as a single shared library holding the whole thing
SET(READDY_BUILD_SHARED_COMBINED OFF CACHE BOOL "Determines if one single shared library or shared libraries with respect to the components should be built")

# Compile in debug mode by default
IF (NOT CMAKE_BUILD_TYPE)
    MESSAGE(STATUS "Setting build type to 'DEBUG' as none was specified.")
    SET(CMAKE_BUILD_TYPE DEBUG)
ELSE (NOT CMAKE_BUILD_TYPE)
    MESSAGE(STATUS "Using specified '${CMAKE_BUILD_TYPE}' build type.")
ENDIF (NOT CMAKE_BUILD_TYPE)

# OSX deployment target
SET(CMAKE_OSX_DEPLOYMENT_TARGET 10.9)

# decide if api wrappers are being generated
SET(READDY_BUILD_PYTHON_WRAPPER ON CACHE BOOL "Build Python wrappers")

# output directories
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# where to install
IF(WIN32)
    IF(NOT READDY_INSTALL_PREFIX)
        SET(READDY_INSTALL_PREFIX "$ENV{ProgramFiles}/readdy2")
    ENDIF(NOT READDY_INSTALL_PREFIX)
ELSE(WIN32)
    IF(NOT READDY_INSTALL_PREFIX)
        SET(READDY_INSTALL_PREFIX "/usr/local/readdy2")
    ENDIF(NOT READDY_INSTALL_PREFIX)
ENDIF(WIN32)
IF(${CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT})
    SET(CMAKE_INSTALL_PREFIX ${OPENMM_INSTALL_PREFIX} CACHE PATH "Where to install OpenMM" FORCE)
ENDIF(${CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT})

# includes
SET(READDY_GLOBAL_DIR "${PROJECT_SOURCE_DIR}")
SET(READDY_GLOBAL_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/include")
LIST(APPEND READDY_HEADERS "${PROJECT_SOURCE_DIR}/include")

# boost.dll
INCLUDE(SetupBoostDLL)

# extra compile and link flags
SET(EXTRA_COMPILE_FLAGS "-std=c++11")
IF (APPLE)
    SET(EXTRA_COMPILE_FLAGS "${EXTRA_COMPILE_FLAGS} -stdlib=libstdc++")
ENDIF(APPLE)
SET(EXTRA_LINK_FLAGS ${EXTRA_COMPILE_FLAGS})
IF (CMAKE_SYSTEM_NAME MATCHES "Linux")
    IF (NOT ANDROID)
        SET(EXTRA_LINK_FLAGS "${EXTRA_LINK_FLAGS} -Wl,--no-as-needed -lrt")
    ENDIF (NOT ANDROID)
ENDIF(CMAKE_SYSTEM_NAME MATCHES "Linux")

# link dynamically against boost log
ADD_DEFINITIONS(-DBOOST_LOG_DYN_LINK)

# main project name
SET(READDY_PROJECT_NAME "ReaDDy2")

#
# readdy2 core
#
IF(READDY_BUILD_SHARED_COMBINED)
    INCLUDE(cmake/sources_common.cmake)
    INCLUDE(cmake/sources_io.cmake)
    INCLUDE(cmake/sources_plugin.cmake)
    INCLUDE(cmake/sources_readdy.cmake)
    # create library
    ADD_LIBRARY(${READDY_PROJECT_NAME} ${READDY_ALL_SOURCES} ${READDY_GLOBAL_INCLUDE_DIR})
    TARGET_INCLUDE_DIRECTORIES(${READDY_PROJECT_NAME} PUBLIC ${READDY_INCLUDE_DIRS})
    SET_TARGET_PROPERTIES(${READDY_PROJECT_NAME} PROPERTIES
            LINK_FLAGS "${EXTRA_LINK_FLAGS}"
            COMPILE_FLAGS "${EXTRA_COMPILEFLAGS}"
    )
    TARGET_LINK_LIBRARIES(${READDY_PROJECT_NAME} ${Boost_LIBRARIES} ${CMAKE_DL_LIBS})
ELSE(READDY_BUILD_SHARED_COMBINED)
    # include directories from least dependent to most dependent
    # i.e.: A depends on B -> add subdirectory of B first
    ADD_SUBDIRECTORY(readdy2/main/common)
    ADD_SUBDIRECTORY(readdy2/main/io)
    ADD_SUBDIRECTORY(readdy2/main/plugin)
    ADD_SUBDIRECTORY(readdy2/main)
ENDIF(READDY_BUILD_SHARED_COMBINED)


#
# readdy2 kernels
#
ADD_SUBDIRECTORY(kernels)

#
# installation (for subprojects handled in their respective build files)
#
IF(READDY_BUILD_SHARED_COMBINED)
    INSTALL(TARGETS ${READDY_PROJECT_NAME}
            RUNTIME DESTINATION bin
            LIBRARY DESTINATION lib
            ARCHIVE DESTINATION lib)
    INSTALL_TARGETS(/lib RUNTIME_DIRECTORY /lib ${READDY_PROJECT_NAME})
ENDIF(READDY_BUILD_SHARED_COMBINED)
INSTALL(DIRECTORY include/ DESTINATION include)


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
    ADD_TEST(NAME readdy-tests COMMAND runUnitTests)
ELSE ()
    MESSAGE(STATUS "The submodule googletest was not initialized or READDY_RUN_TESTS was OFF, thus no testing.")
ENDIF ()