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

if (POLICY CMP0048)
    cmake_policy(SET CMP0048 NEW)
endif (POLICY CMP0048)
PROJECT(readdy_main VERSION 0.0 LANGUAGES C CXX)

#####################################
#                                   #
# User configurable cache variables #
#                                   #
#####################################
set(INSTALL_DIR "${CMAKE_INSTALL_PREFIX}" CACHE PATH "Installation root for anything but python bindings and files")
set(PYTHON_INSTALL_DIR "${SP_DIR}/readdy" CACHE PATH "Installation directory relative to SP_DIR (python/site-packages)")

# test configuration
set(READDY_CREATE_TEST_TARGET ON CACHE BOOL "Determines if the test target should be added.")
set(READDY_INSTALL_UNIT_TEST_EXECUTABLE OFF CACHE BOOL "Determines if the unit test executable should be installed. Requires the cache variable READDY_CREATE_TEST_TARGET to be ON.")
set(READDY_KERNELS_TO_TEST "SingleCPU,CPU,MPI" CACHE STRING "Variable holding which kernels should be tested in the core library.")
set(READDY_INSTALL_HEADERS OFF CACHE BOOL "Whether to install header files")

# decide if api wrappers are being generated
set(READDY_BUILD_PYTHON_WRAPPER ON CACHE BOOL "Build Python wrappers")
set(READDY_DEBUG_PYTHON_MODULES OFF CACHE BOOL "If this flag is set to ON, the compiled python modules get installed directly into the source tree.")
set(READDY_DEBUG_CONDA_ROOT_DIR "" CACHE PATH "This only does something, if READDY_DEBUG_PYTHON_MODULES is set to ON. Then, it compiles the readdy libraries into the @conda@/lib directory.")

# documentation target (make doc)
set(READDY_GENERATE_DOCUMENTATION_TARGET OFF CACHE BOOL "Should generate target for documentation or not.")

# log configuration variables of interest
set(READDY_LOG_CMAKE_CONFIGURATION OFF CACHE BOOL "If turned on, the status of relevant cmake variables is logged at configuration time.")

# build the mpi kernel
set(READDY_BUILD_MPI_KERNEL OFF CACHE BOOL "Whether to build the MPI kernel or not")

# build the scenarios executables for performance benchmarking
set(READDY_BUILD_SCENARIOS OFF CACHE BOOL "Whether to build the Scenarios executable or not")

option(
        IGNORE_VENDORED_DEPENDENCIES
        "Ignore libraries provided as submodules of the repository"
        OFF
)
mark_as_advanced(IGNORE_VENDORED_DEPENDENCIES)


#####################################
#                                   #
# Basic setup of the project        #
#                                   #
#####################################

# minimum cmake version
cmake_minimum_required(VERSION 3.16)

# extend cmake module path to include readdy custom modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules")

# ensure that the build directory is not the source directory
include(MacroOutOfSourceBuild)
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

# we use c++20, but c++17 on macOS to avoid SDK compatibility issues
SET(CMAKE_CXX_STANDARD 20)
SET(CMAKE_CXX_STANDARD_REQUIRED ON)
if(APPLE)
    set(CMAKE_MACOSX_RPATH TRUE)  # any shared lib being built has rpath enabled
    # Force C++17 on macOS to avoid C++20/SDK compatibility issues with Catch2 and other dependencies
    SET(CMAKE_CXX_STANDARD 17)
    SET(CMAKE_CXX_STANDARD_REQUIRED OFF)
endif()

# always build readdy as shared lib
SET(BUILD_SHARED_LIBS ON CACHE BOOL "Build shared libs" FORCE)

# dir vars
SET(READDY_GLOBAL_DIR "${PROJECT_SOURCE_DIR}")
SET(READDY_GLOBAL_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/include")

# Compile in debug mode by default
IF (NOT CMAKE_BUILD_TYPE)
    MESSAGE(STATUS "Setting build type to 'DEBUG' as none was specified.")
    SET(CMAKE_BUILD_TYPE DEBUG)
ELSE (NOT CMAKE_BUILD_TYPE)
    MESSAGE(STATUS "Using specified '${CMAKE_BUILD_TYPE}' build type.")
ENDIF (NOT CMAKE_BUILD_TYPE)

SET(EXTRA_COMPILE_FLAGS "")
SET(EXTRA_LINK_FLAGS ${EXTRA_COMPILE_FLAGS})
IF (CMAKE_SYSTEM_NAME MATCHES "Linux")
    IF (NOT ANDROID)
        if (CMAKE_BUILD_TYPE STREQUAL "Release")
            SET(EXTRA_COMPILE_FLAGS "${EXTRA_COMPILE_FLAGS}")
        endif()
        SET(EXTRA_LINK_FLAGS "${EXTRA_LINK_FLAGS} -Wl,--no-as-needed -lrt")
    ENDIF (NOT ANDROID)
ENDIF(CMAKE_SYSTEM_NAME MATCHES "Linux")

# main project name
set(READDY_PROJECT_NAME "readdy")

if (READDY_LOG_CMAKE_CONFIGURATION)
    include(LogConfiguration)
endif (READDY_LOG_CMAKE_CONFIGURATION)

#####################################
#                                   #
# setup testing                     #
#                                   #
#####################################

if (READDY_CREATE_TEST_TARGET)
    add_definitions(-DCATCH_CONFIG_NO_CPP17_UNCAUGHT_EXCEPTIONS)  # workaround for macos which chokes on some cpp17 stuff
    add_subdirectory(contrib/Catch2)
    include(${CMAKE_CURRENT_LIST_DIR}/contrib/Catch2/extras/Catch.cmake)
    add_subdirectory(readdy_testing)
else()
    # disable testing if there was no submodule
    SET(READDY_CREATE_TEST_TARGET OFF)
    # no installation of test executable
    SET(READDY_INSTALL_UNIT_TEST_EXECUTABLE OFF)
endif()

#####################################
#                                   #
# readdy core                       #
#                                   #
#####################################
if (EXISTS ${CMAKE_CURRENT_LIST_DIR}/contrib/pybind11/CMakeLists.txt)
    add_subdirectory(contrib/pybind11)
else()
    find_package(pybind11 REQUIRED)
endif()

if (NOT IGNORE_VENDORED_DEPENDENCIES)
    set(FMT_SYSTEM_HEADERS ON)
    add_subdirectory(
            ${CMAKE_SOURCE_DIR}/contrib/fmt
            EXCLUDE_FROM_ALL
    )
    set(JSON_Install OFF CACHE INTERNAL "")
    add_subdirectory(
            ${CMAKE_SOURCE_DIR}/contrib/json
            EXCLUDE_FROM_ALL
    )
    add_subdirectory(${CMAKE_SOURCE_DIR}/contrib/spdlog)
else()
    find_package(fmt REQUIRED)
    find_package(spdlog REQUIRED)
    find_package(nlohmann_json REQUIRED)
endif()
find_package(blosc REQUIRED)

set(WINDOWS_EXPORT_ALL_SYMBOLS ON)

include(cmake/sources/sources_common.cmake)
include(cmake/sources/sources_io.cmake)
include(cmake/sources/sources_model.cmake)
include(cmake/sources/sources_plugin.cmake)
include(cmake/sources/kernels/singlecpu.cmake)
include(cmake/sources/kernels/cpu.cmake)
include(cmake/sources/sources_readdy.cmake)
# create library
message(STATUS "Libraries ReaDDy depends on: ${READDY_DEPENDENT_LIBRARIES}")
message(STATUS "Include directories: ${READDY_INCLUDE_DIRS}")
add_library(${READDY_PROJECT_NAME} ${READDY_ALL_SOURCES} ${READDY_GLOBAL_INCLUDE_DIR})
target_include_directories(${READDY_PROJECT_NAME} PUBLIC ${READDY_INCLUDE_DIRS})
target_link_libraries(${READDY_PROJECT_NAME} PUBLIC ${READDY_DEPENDENT_LIBRARIES})
set_target_properties(${READDY_PROJECT_NAME} PROPERTIES
        LINK_FLAGS "${EXTRA_LINK_FLAGS}"
        COMPILE_FLAGS "${EXTRA_COMPILE_FLAGS}")
target_compile_definitions(${READDY_PROJECT_NAME} PUBLIC LIBRARY_EXPORTS=1)

#####################################
#                                   #
# readdy kernels                    #
#                                   #
#####################################
add_subdirectory(kernels)

#####################################
#                                   #
# installation (for subprojects     #
# handled in their respective       #
# build files)                      #
#                                   #
#####################################
if(READDY_INSTALL_HEADERS)
    install (
            DIRECTORY ${CMAKE_SOURCE_DIR}/include/
            DESTINATION ${INSTALL_DIR}/include
            FILES_MATCHING PATTERN "*.h*")
    install (
            DIRECTORY ${CMAKE_SOURCE_DIR}/kernels/cpu/include/
            DESTINATION ${INSTALL_DIR}/include
            FILES_MATCHING PATTERN "*.h*"
    )
    install (
            DIRECTORY ${CMAKE_SOURCE_DIR}/libraries/h5rd/include/
            DESTINATION ${INSTALL_DIR}/include
            FILES_MATCHING PATTERN "*.h*")
    install (
            DIRECTORY ${CMAKE_SOURCE_DIR}/libraries/json/include/
            DESTINATION ${INSTALL_DIR}/include
            FILES_MATCHING PATTERN "*.h*"
    )
    install (
            DIRECTORY ${CMAKE_SOURCE_DIR}/libraries/spdlog/include/
            DESTINATION ${INSTALL_DIR}/include
            FILES_MATCHING PATTERN "*.h*"
    )
    install (
            DIRECTORY ${CMAKE_SOURCE_DIR}/libraries/spdlog/include/
            DESTINATION ${INSTALL_DIR}/include
            FILES_MATCHING PATTERN "*.cc*"
    )
    install (
            DIRECTORY ${CMAKE_SOURCE_DIR}/libraries/graph/
            DESTINATION ${INSTALL_DIR}/include
            FILES_MATCHING PATTERN "*.h*"
    )
endif()
install(TARGETS ${READDY_PROJECT_NAME} LIBRARY DESTINATION "${INSTALL_DIR}/lib")

#####################################
#                                   #
# readdy wrappers                   #
#                                   #
#####################################
IF (READDY_BUILD_PYTHON_WRAPPER)
    add_subdirectory(wrappers)
ENDIF (READDY_BUILD_PYTHON_WRAPPER)

#####################################
#                                   #
# tests                             #
#                                   #
#####################################
if (READDY_CREATE_TEST_TARGET)
    add_subdirectory(readdy/test)
    add_test(NAME readdy-tests COMMAND runUnitTests)
    add_subdirectory(kernels/singlecpu/test)
    add_subdirectory(kernels/cpu/test)
endif ()

if(READDY_BUILD_SCENARIOS)
    add_subdirectory(examples/scenarios)
    if(READDY_BUILD_MPI_KERNEL)
        add_subdirectory(examples/mpi_scenarios)
    endif()
endif()


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
