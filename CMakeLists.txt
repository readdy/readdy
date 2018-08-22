####################################################################
# Copyright © 2018 Computational Molecular Biology Group,          #
#                  Freie Universität Berlin (GER)                  #
#                                                                  #
# Redistribution and use in source and binary forms, with or       #
# without modification, are permitted provided that the            #
# following conditions are met:                                    #
#  1. Redistributions of source code must retain the above         #
#     copyright notice, this list of conditions and the            #
#     following disclaimer.                                        #
#  2. Redistributions in binary form must reproduce the above      #
#     copyright notice, this list of conditions and the following  #
#     disclaimer in the documentation and/or other materials       #
#     provided with the distribution.                              #
#  3. Neither the name of the copyright holder nor the names of    #
#     its contributors may be used to endorse or promote products  #
#     derived from this software without specific                  #
#     prior written permission.                                    #
#                                                                  #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND           #
# CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,      #
# INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF         #
# MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE         #
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR            #
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,     #
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,         #
# BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; #
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER #
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,      #
# STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)    #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF      #
# ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                       #
####################################################################

PROJECT(readdy_main C CXX)

#####################################
#                                   #
# User configurable cache variables #
#                                   #
#####################################

# test configuration
SET(READDY_CREATE_TEST_TARGET ON CACHE BOOL "Determines if the test target should be added.")
SET(READDY_CREATE_MEMORY_CHECK_TEST_TARGET OFF CACHE BOOL "Determines if valgrind test targets should be added. Requires the cache variable READDY_CREATE_TEST_TARGET to be ON.")
SET(READDY_INSTALL_UNIT_TEST_EXECUTABLE OFF CACHE BOOL "Determines if the unit test executable should be installed. Requires the cache variable READDY_CREATE_TEST_TARGET to be ON.")
SET(READDY_KERNELS_TO_TEST "SingleCPU,CPU" CACHE STRING "Variable holding which kernels should be tested in the core library.")

# build single libraries or as a single shared library holding the whole thing
SET(READDY_BUILD_SHARED_COMBINED OFF CACHE BOOL "Determines if one single shared library or shared libraries with respect to the components should be built")

# decide if api wrappers are being generated
SET(READDY_BUILD_PYTHON_WRAPPER ON CACHE BOOL "Build Python wrappers")
SET(READDY_DEBUG_PYTHON_MODULES OFF CACHE BOOL "If this flag is set to ON, the compiled python modules get installed directly into the source tree.")
SET(READDY_DEBUG_CONDA_ROOT_DIR "" CACHE PATH "This only does something, if READDY_DEBUG_PYTHON_MODULES is set to ON. Then, it compiles the readdy libraries into the @conda@/lib directory.")

# documentation target (make doc)
SET(READDY_GENERATE_DOCUMENTATION_TARGET OFF CACHE BOOL "Should generate target for documentation or not.")
SET(READDY_GENERATE_DOCUMENTATION_TARGET_ONLY OFF CACHE BOOL "If set to true, this only generates the documentation target.")

# log configuration variables of interest
SET(READDY_LOG_CMAKE_CONFIGURATION OFF CACHE BOOL "If turned on, the status of relevant cmake variables is logged at configuration time.")

#####################################
#                                   #
# Basic setup of the project        #
#                                   #
#####################################

# minimum cmake version
CMAKE_MINIMUM_REQUIRED(VERSION 3.2.0)

# extend cmake module path to include readdy custom modules
SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules")

# ensure that the build directory is not the source directory
INCLUDE(MacroOutOfSourceBuild)
MACRO_ENSURE_OUT_OF_SOURCE_BUILD("${PROJECT_NAME} requires an out source build.")

# makefile settings
SET(CMAKE_VERBOSE_MAKEFILE OFF)
SET(CMAKE_COLOR_MAKEFILE ON)

# version
set(READDY_VERSION "v0.0" CACHE STRING "The version")
mark_as_advanced(READDY_VERSION)
set(READDY_BUILD_STRING "dev" CACHE STRING "Build string")
mark_as_advanced(READDY_BUILD_STRING)

message(STATUS "====================================================")
message(STATUS "CMake configure for readdy version ${READDY_VERSION}-${READDY_BUILD_STRING}")
message(STATUS "====================================================")

# we use c++17
SET(CMAKE_CXX_STANDARD 17)
SET(CMAKE_CXX_STANDARD_REQUIRED ON)

# always build readdy as shared lib
SET(BUILD_SHARED_LIBS ON CACHE BOOL "Build shared libs" FORCE)

# dir vars
SET(READDY_GLOBAL_DIR "${PROJECT_SOURCE_DIR}")
SET(READDY_GLOBAL_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/include")

INCLUDE(CheckCXXCompilerFlag)
function(enable_cxx_compiler_flag_if_supported flag)
    string(FIND "${CMAKE_CXX_FLAGS}" "${flag}" flag_already_set)
    if(flag_already_set EQUAL -1)
        check_cxx_compiler_flag("${flag}" flag_supported)
        if(flag_supported)
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${flag}" PARENT_SCOPE)
        endif()
    endif()
endfunction()
#enable_cxx_compiler_flag_if_supported("-pedantic")
#enable_cxx_compiler_flag_if_supported("-fms-extensions")


IF(NOT READDY_GENERATE_DOCUMENTATION_TARGET_ONLY)
    # Compile in debug mode by default
    IF (NOT CMAKE_BUILD_TYPE)
        MESSAGE(STATUS "Setting build type to 'DEBUG' as none was specified.")
        SET(CMAKE_BUILD_TYPE DEBUG)
    ELSE (NOT CMAKE_BUILD_TYPE)
        MESSAGE(STATUS "Using specified '${CMAKE_BUILD_TYPE}' build type.")
    ENDIF (NOT CMAKE_BUILD_TYPE)

    string(TOLOWER "${CMAKE_BUILD_TYPE}" BUILD_TYPE_LOWER)
    if (BUILD_TYPE_LOWER STREQUAL "release")
        message(STATUS "--- release build type ---")
    else()
        message(STATUS "--- debug build type, specifying debug flag ---")
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DREADDY_DEBUG")
    endif()

    # output directories
    IF(READDY_DEBUG_PYTHON_MODULES)
        SET(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${READDY_DEBUG_CONDA_ROOT_DIR}/lib)
        SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${READDY_DEBUG_CONDA_ROOT_DIR}/lib)
        SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${READDY_DEBUG_CONDA_ROOT_DIR}/bin)
        SET(READDY_PLUGIN_OUTPUT_DIRECTORY ${READDY_DEBUG_CONDA_ROOT_DIR}/readdy/readdy_plugins)
    ELSE()
        SET(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
        SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
        SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
        SET(READDY_PLUGIN_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/readdy/readdy_plugins)
    ENDIF()

    # where to install
    IF(WIN32)
        IF(NOT READDY_INSTALL_PREFIX)
            SET(READDY_INSTALL_PREFIX "$ENV{ProgramFiles}/readdy")
        ENDIF(NOT READDY_INSTALL_PREFIX)
    ELSE(WIN32)
        IF(NOT READDY_INSTALL_PREFIX)
            SET(READDY_INSTALL_PREFIX "/usr/local/readdy")
        ENDIF(NOT READDY_INSTALL_PREFIX)
    ENDIF(WIN32)
    IF(${CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT})
        SET(CMAKE_INSTALL_PREFIX ${READDY_INSTALL_PREFIX} CACHE PATH "Where to install readdy" FORCE)
    ENDIF(${CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT})

    SET(EXTRA_COMPILE_FLAGS "")
    SET(EXTRA_LINK_FLAGS ${EXTRA_COMPILE_FLAGS})
    IF(APPLE)
        SET(CMAKE_MACOSX_RPATH 1)
        IF (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
            SET(CMAKE_XCODE_ATTRIBUTE_CLANG_CXX_LIBRARY "libc++")
            # SET(EXTRA_COMPILE_FLAGS "${EXTRA_COMPILE_FLAGS} -stdlib=libc++")
            # SET(EXTRA_LINK_FLAGS "${EXTRA_LINK_FLAGS} -v -stdlib=libc++ -lc++abi")
            if (CMAKE_BUILD_TYPE STREQUAL "Release")
                SET(EXTRA_COMPILE_FLAGS "${EXTRA_COMPILE_FLAGS}")
            endif()
        ENDIF (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    ENDIF(APPLE)
    IF (CMAKE_SYSTEM_NAME MATCHES "Linux")
        IF (NOT ANDROID)
            if (CMAKE_BUILD_TYPE STREQUAL "Release")
                SET(EXTRA_COMPILE_FLAGS "${EXTRA_COMPILE_FLAGS}") # -ffast-math
            endif()
            SET(EXTRA_LINK_FLAGS "${EXTRA_LINK_FLAGS} -Wl,--no-as-needed -lrt")
        ENDIF (NOT ANDROID)
    ENDIF(CMAKE_SYSTEM_NAME MATCHES "Linux")

    # main project name
    SET(READDY_PROJECT_NAME "readdy")

    IF (READDY_LOG_CMAKE_CONFIGURATION)
        INCLUDE(LogConfiguration)
    ENDIF (READDY_LOG_CMAKE_CONFIGURATION)

    #####################################
    #                                   #
    # setup submodules                  #
    #                                   #
    #####################################

    MESSAGE(STATUS "Initializing and updating submodules")
    EXECUTE_PROCESS(COMMAND git submodule update --init WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

    #####################################
    #                                   #
    # setup testing                     #
    #                                   #
    #####################################

    IF (READDY_CREATE_TEST_TARGET)

       if (CMAKE_VERSION VERSION_LESS 3.2)
           set(UPDATE_DISCONNECTED_IF_AVAILABLE "")
       else()
           set(UPDATE_DISCONNECTED_IF_AVAILABLE "UPDATE_DISCONNECTED 1")
       endif()

       include(cmake/Modules/DownloadProject.cmake)
       download_project(PROJ googletest
               GIT_REPOSITORY https://github.com/google/googletest.git
               GIT_TAG release-1.8.0
               ${UPDATE_DISCONNECTED_IF_AVAILABLE}
               )

       # Prevent GoogleTest from overriding our compiler/linker options
       # when building with Visual Studio
       set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

       add_subdirectory(${googletest_SOURCE_DIR} ${googletest_BINARY_DIR})

       ENABLE_TESTING()

       ADD_SUBDIRECTORY(readdy_testing)

    ELSE()
        MESSAGE(STATUS "The submodule googletest was not initialized or READDY_CREATE_TEST_TARGET was OFF, thus no testing.")
        # disable testing if there was no submodule
        SET(READDY_CREATE_TEST_TARGET OFF)
        # if tests are disabled, so are the corresponding memory checks
        SET(READDY_CREATE_MEMORY_CHECK_TEST_TARGET OFF)
        # no installation of test executable
        SET(READDY_INSTALL_UNIT_TEST_EXECUTABLE OFF)
    ENDIF()

    #####################################
    #                                   #
    # readdy core                       #
    #                                   #
    #####################################

    # include blosc
    SET(DEACTIVATE_SNAPPY ON CACHE BOOL "")
    SET(DEACTIVATE_ZSTD OFF CACHE BOOL "")
    SET(DEACTIVATE_ZLIB ON CACHE BOOL "")
    SET(BLOSC_INSTALL ON CACHE BOOL "")
    SET(BUILD_SHARED ON CACHE BOOL "")
    SET(BUILD_STATIC OFF CACHE BOOL "")
    SET(BUILD_TESTS OFF CACHE BOOL "")
    SET(BUILD_BENCHMARKS OFF CACHE BOOL "")
    ADD_SUBDIRECTORY(libraries/c-blosc)

    IF(READDY_BUILD_SHARED_COMBINED)
        INCLUDE(cmake/sources/sources_common.cmake)
        INCLUDE(cmake/sources/sources_io.cmake)
        INCLUDE(cmake/sources/sources_model.cmake)
        INCLUDE(cmake/sources/sources_plugin.cmake)
        INCLUDE(cmake/sources/kernels/singlecpu.cmake)
        INCLUDE(cmake/sources/sources_readdy.cmake)
        # create library
        MESSAGE(STATUS "Libraries ReaDDy depends on: ${READDY_DEPENDENT_LIBRARIES}")
        MESSAGE(STATUS "Include directories: ${READDY_INCLUDE_DIRS}")
        ADD_LIBRARY(${READDY_PROJECT_NAME} ${READDY_ALL_SOURCES} ${READDY_GLOBAL_INCLUDE_DIR})
        TARGET_INCLUDE_DIRECTORIES(${READDY_PROJECT_NAME} PUBLIC ${READDY_INCLUDE_DIRS})
        TARGET_LINK_LIBRARIES(${READDY_PROJECT_NAME} ${READDY_DEPENDENT_LIBRARIES})
        SET_TARGET_PROPERTIES(${READDY_PROJECT_NAME} PROPERTIES
                LINK_FLAGS "${EXTRA_LINK_FLAGS}"
                COMPILE_FLAGS "${EXTRA_COMPILE_FLAGS}")
    ELSE(READDY_BUILD_SHARED_COMBINED)
        # include directories from least dependent to most dependent
        # i.e.: A depends on B -> add subdirectory of B first
        ADD_SUBDIRECTORY(readdy/main/common)
        ADD_SUBDIRECTORY(readdy/main/io)
        ADD_SUBDIRECTORY(readdy/main/model)
        ADD_SUBDIRECTORY(readdy/main/plugin)
        ADD_SUBDIRECTORY(kernels/singlecpu)
        ADD_SUBDIRECTORY(readdy/main)
    ENDIF(READDY_BUILD_SHARED_COMBINED)

    #####################################
    #                                   #
    # readdy kernels                    #
    #                                   #
    #####################################
    ADD_SUBDIRECTORY(kernels)

    #####################################
    #                                   #
    # installation (for subprojects     #
    # handled in their respective       #
    # build files)                      #
    #                                   #
    #####################################
    IF(READDY_BUILD_SHARED_COMBINED)
        INSTALL (
                DIRECTORY ${CMAKE_SOURCE_DIR}/include/
                DESTINATION include
                FILES_MATCHING PATTERN "*.h*")
        INSTALL (
                DIRECTORY ${CMAKE_SOURCE_DIR}/kernels/cpu/include/
                DESTINATION include
                FILES_MATCHING PATTERN "*.h*"
        )
        INSTALL (
                DIRECTORY ${CMAKE_SOURCE_DIR}/libraries/h5rd/include/
                DESTINATION include
                FILES_MATCHING PATTERN "*.h*")
        INSTALL (
                DIRECTORY ${CMAKE_SOURCE_DIR}/libraries/json/include/
                DESTINATION include
                FILES_MATCHING PATTERN "*.h*"
        )
        INSTALL (
                DIRECTORY ${CMAKE_SOURCE_DIR}/libraries/spdlog/include/
                DESTINATION include
                FILES_MATCHING PATTERN "*.h*"
        )
        INSTALL (
                DIRECTORY ${CMAKE_SOURCE_DIR}/libraries/spdlog/include/
                DESTINATION include
                FILES_MATCHING PATTERN "*.cc*"
        )
        INSTALL(TARGETS ${READDY_PROJECT_NAME}
                RUNTIME DESTINATION bin
                LIBRARY DESTINATION lib
                ARCHIVE DESTINATION lib)
        INSTALL_TARGETS(/lib RUNTIME_DIRECTORY /lib ${READDY_PROJECT_NAME})
    ENDIF(READDY_BUILD_SHARED_COMBINED)

    #####################################
    #                                   #
    # readdy wrappers                   #
    #                                   #
    #####################################
    IF (READDY_BUILD_PYTHON_WRAPPER)
        ADD_SUBDIRECTORY(libraries/pybind11)
        ADD_SUBDIRECTORY(wrappers)
    ENDIF (READDY_BUILD_PYTHON_WRAPPER)

    #####################################
    #                                   #
    # tests                             #
    #                                   #
    #####################################
    IF (READDY_CREATE_TEST_TARGET)
        ADD_SUBDIRECTORY(readdy/test)
        ADD_TEST(NAME readdy-tests COMMAND runUnitTests)
        ADD_SUBDIRECTORY(kernels/singlecpu/test)
    ENDIF ()

    #####################################
    #                                   #
    # memory checks for the tests       #
    #                                   #
    #####################################
    IF(READDY_CREATE_TEST_TARGET)
        IF(READDY_CREATE_MEMORY_CHECK_TEST_TARGET)
            FIND_PACKAGE(Valgrind QUIET)
            IF (NOT VALGRIND_FOUND)
                MESSAGE(STATUS "Valgrind not found, therefore no memory checks.")
            ELSE (NOT VALGRIND_FOUND)
                MESSAGE(STATUS "Valgrind found, memory checks enabled.")
                SET(VALGRIND_OPTIONS --tool=memcheck --leak-check=yes --show-reachable=yes --num-callers=30 --leak-check=summary --leak-resolution=high)

                FUNCTION(ADD_MEMCHECK_TEST name binary)
                    SET(memcheck_command "valgrind ${VALGRIND_OPTIONS}")
                    SEPARATE_ARGUMENTS(memcheck_command)
                    ADD_TEST(NAME memcheck_${name} COMMAND ${memcheck_command} "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${binary}" ${ARGN})
                ENDFUNCTION(ADD_MEMCHECK_TEST)

                # can be called by "ctest -V -R memcheck_readdy-tests"
                ADD_MEMCHECK_TEST(readdy-tests runUnitTests)
            ENDIF (NOT VALGRIND_FOUND)
        ENDIF(READDY_CREATE_MEMORY_CHECK_TEST_TARGET)
    ENDIF(READDY_CREATE_TEST_TARGET)
ELSE(NOT READDY_GENERATE_DOCUMENTATION_TARGET_ONLY)
    SET(READDY_GENERATE_DOCUMENTATION_TARGET ON)
ENDIF(NOT READDY_GENERATE_DOCUMENTATION_TARGET_ONLY)


#####################################
#                                   #
# documentation (doxygen)           #
#                                   #
#####################################
IF(READDY_GENERATE_DOCUMENTATION_TARGET)
    FIND_PACKAGE (Doxygen REQUIRED)
    FIND_PROGRAM (DOXYFILE_MAKE make)
    MARK_AS_ADVANCED (DOXYFILE_MAKE)
    CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/docs/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
    ADD_CUSTOM_TARGET(doc
            ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
            COMMENT "Generating API documentation with Doxygen" VERBATIM)
    IF(MAKEINDEX_COMPILER AND DOXYFILE_MAKE)
        ADD_CUSTOM_COMMAND(
                TARGET doc
                POST_BUILD
                COMMAND "${DOXYFILE_MAKE}"
                COMMENT	"Running Doxygen in ${CMAKE_CURRENT_BINARY_DIR}/docs..."
                WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/docs")
    ENDIF()
ENDIF(READDY_GENERATE_DOCUMENTATION_TARGET)

if(${false})
    INCLUDE(cmake/sources/sources_common.cmake)
    INCLUDE(cmake/sources/sources_io.cmake)
    INCLUDE(cmake/sources/sources_model.cmake)
    INCLUDE(cmake/sources/sources_plugin.cmake)
    INCLUDE(cmake/sources/kernels/singlecpu.cmake)
    INCLUDE(cmake/sources/sources_readdy.cmake)
    include(cmake/Modules/clang-tidy.cmake)
    add_clang_tidy_target(clang-tidy "${READDY_ALL_SOURCES}" "${READDY_INCLUDE_DIRS}")
endif()