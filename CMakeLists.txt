#####################################
#                                   #
# ReaDDy                           #
#                                   #
#####################################
PROJECT(ReaDDy_Main C CXX)

#####################################
#                                   #
# User configurable cache variables #
#                                   #
#####################################

# test configuration
SET(READDY_CREATE_TEST_TARGET ON CACHE BOOL "Determines if the test target should be added.")
SET(READDY_CREATE_MEMORY_CHECK_TEST_TARGET OFF CACHE BOOL "Determines if valgrind test targets should be added. Requires the cache variable READDY_CREATE_TEST_TARGET to be ON.")
SET(READDY_INSTALL_UNIT_TEST_EXECUTABLE OFF CACHE BOOL "Determines if the unit test executable should be installed. Requires the cache variable READDY_CREATE_TEST_TARGET to be ON.")

# build single libraries or as a single shared library holding the whole thing
SET(READDY_BUILD_SHARED_COMBINED OFF CACHE BOOL "Determines if one single shared library or shared libraries with respect to the components should be built")

# decide if api wrappers are being generated
SET(READDY_BUILD_PYTHON_WRAPPER ON CACHE BOOL "Build Python wrappers")
SET(READDY_DEBUG_PYTHON_MODULES OFF CACHE BOOL "If this flag is set to ON, the compiled python modules get installed directly into the source tree.")
SET(READDY_DEBUG_CONDA_ROOT_DIR "" CACHE PATH "This only does something, if READDY_DEBUG_PYTHON_MODULES is set to ON. Then, it compiles the readdy libraries into the @conda@/lib directory.")

# documentation target (make doc)
SET(READDY_GENERATE_DOCUMENTATION_TARGET OFF CACHE BOOL "Should generate target for documentation or not")

# log configuration variables of interest
SET(READDY_LOG_CMAKE_CONFIGURATION OFF CACHE BOOL "If turned on, the status of relevant cmake variables is logged at configuration time.")

# boost download cache dir
SET(BOOST_DOWNLOAD_CACHE_DIR "${CMAKE_BINARY_DIR}/download/external/boost" CACHE PATH "Path to the directory into which the boost libraries are downloaded and installed.")

#####################################
#                                   #
# Basic setup of the project        #
#                                   #
#####################################

# minimum cmake version
CMAKE_MINIMUM_REQUIRED(VERSION 3.2.0)

# extend cmake module path to include readdy custom modules
SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules")

# ensure that the build directory is not the source directory
INCLUDE(MacroOutOfSourceBuild)
MACRO_ENSURE_OUT_OF_SOURCE_BUILD("${PROJECT_NAME} requires an out source build.")

# makefile settings
SET(CMAKE_VERBOSE_MAKEFILE OFF)
SET(CMAKE_COLOR_MAKEFILE ON)

# we use c++11
SET(CMAKE_CXX_STANDARD 11)
SET(CMAKE_CXX_STANDARD_REQUIRED ON)

# version
SET(READDY_MAJOR_VERSION 0)
SET(READDY_MINOR_VERSION 0)
SET(READDY_PATCH_VERSION 1)
SET(READDY_VERSION "${READDY_MSETAJOR_VERSION}.${READDY_MINOR_VERSION}.${READDY_PATCH_VERSION}")

# always build readdy as shared lib
SET(BUILD_SHARED_LIBS ON)

# Compile in debug mode by default
IF (NOT CMAKE_BUILD_TYPE)
    MESSAGE(STATUS "Setting build type to 'DEBUG' as none was specified.")
    SET(CMAKE_BUILD_TYPE DEBUG)
ELSE (NOT CMAKE_BUILD_TYPE)
    MESSAGE(STATUS "Using specified '${CMAKE_BUILD_TYPE}' build type.")
ENDIF (NOT CMAKE_BUILD_TYPE)

# output directories
IF(READDY_DEBUG_PYTHON_MODULES)
    SET(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${READDY_DEBUG_CONDA_ROOT_DIR}/lib)
    SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${READDY_DEBUG_CONDA_ROOT_DIR}/lib)
    SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${READDY_DEBUG_CONDA_ROOT_DIR}/bin)
    SET(READDY_PLUGIN_OUTPUT_DIRECTORY ${READDY_DEBUG_CONDA_ROOT_DIR}/lib/readdy_plugins)
ELSE()
    SET(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
    SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
    SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
    SET(READDY_PLUGIN_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib/readdy_plugins)
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
    SET(CMAKE_INSTALL_PREFIX ${READDY_INSTALL_PREFIX} CACHE PATH "Where to install ReaDDy" FORCE)
ENDIF(${CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT})

# includes
SET(READDY_GLOBAL_DIR "${PROJECT_SOURCE_DIR}")
SET(READDY_GLOBAL_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/include")
LIST(APPEND READDY_HEADERS "${PROJECT_SOURCE_DIR}/include")

# download and install boost into cache dir
INCLUDE(SetupBoost)

# setup boost.dll
INCLUDE(SetupBoostDLL)

# extra compile and link flags
SET(EXTRA_COMPILE_FLAGS "-std=c++11")
SET(EXTRA_LINK_FLAGS ${EXTRA_COMPILE_FLAGS})
IF(APPLE)
    IF (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        SET(EXTRA_COMPILE_FLAGS "${EXTRA_COMPILE_FLAGS} -stdlib=libc++")
        SET(EXTRA_LINK_FLAGS "${EXTRA_LINK_FLAGS} -v -stdlib=libc++ -lc++abi")
    ENDIF (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
ENDIF(APPLE)
IF (CMAKE_SYSTEM_NAME MATCHES "Linux")
    IF (NOT ANDROID)
        SET(EXTRA_LINK_FLAGS "${EXTRA_LINK_FLAGS} -Wl,--no-as-needed -lrt")
    ENDIF (NOT ANDROID)
ENDIF(CMAKE_SYSTEM_NAME MATCHES "Linux")

# link dynamically against boost log
ADD_DEFINITIONS(-DBOOST_LOG_DYN_LINK)

# main project name
SET(READDY_PROJECT_NAME "ReaDDy")

IF (READDY_LOG_CMAKE_CONFIGURATION)
    INCLUDE(LogConfiguration)
ENDIF (READDY_LOG_CMAKE_CONFIGURATION)

#####################################
#                                   #
# setup testing                     #
#                                   #
#####################################

IF (EXISTS "${CMAKE_SOURCE_DIR}/libraries/googletest/CMakeLists.txt" AND READDY_CREATE_TEST_TARGET)

    SET(GOOGLETEST_ROOT ${CMAKE_SOURCE_DIR}/libraries/googletest)
    SET(GOOGLETEST_INCLUDE ${GOOGLETEST_ROOT}/googletest/include)
    SET(GOOGLETEST_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/libraries/googletest)
    SET(GOOGLETEST_LIBRARY ${GOOGLETEST_BINARY_DIR}/googletest)

    ADD_SUBDIRECTORY(${GOOGLETEST_ROOT})
    SET(BUILD_GTEST ON)
    SET(BUILD_GMOCK OFF)
    SET(gtest_force_shared_crt ON)
    SET(BUILD_SHARED_LIBS ON)
    SET_TARGET_PROPERTIES(gmock gmock_main gtest gtest_main PROPERTIES
            LINK_FLAGS "${EXTRA_LINK_FLAGS}"
            COMPILE_FLAGS "${EXTRA_COMPILE_FLAGS}"
            )

    ENABLE_TESTING()

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

IF(READDY_BUILD_SHARED_COMBINED)
    INCLUDE(cmake/sources/sources_common.cmake)
    INCLUDE(cmake/sources/sources_model.cmake)
    INCLUDE(cmake/sources/sources_io.cmake)
    INCLUDE(cmake/sources/sources_plugin.cmake)
    INCLUDE(cmake/sources/sources_readdy.cmake)
    # create library
    ADD_LIBRARY(${READDY_PROJECT_NAME} ${READDY_ALL_SOURCES} ${READDY_GLOBAL_INCLUDE_DIR})
    TARGET_INCLUDE_DIRECTORIES(${READDY_PROJECT_NAME} PUBLIC ${READDY_INCLUDE_DIRS})
    SET_TARGET_PROPERTIES(${READDY_PROJECT_NAME} PROPERTIES
            LINK_FLAGS "${EXTRA_LINK_FLAGS}"
            COMPILE_FLAGS "${EXTRA_COMPILE_FLAGS}"
    )
    TARGET_LINK_LIBRARIES(${READDY_PROJECT_NAME} ${READDY_DEPENDENT_LIBRARIES})
ELSE(READDY_BUILD_SHARED_COMBINED)
    # include directories from least dependent to most dependent
    # i.e.: A depends on B -> add subdirectory of B first
    ADD_SUBDIRECTORY(readdy/main/common)
    ADD_SUBDIRECTORY(readdy/main/model)
    ADD_SUBDIRECTORY(readdy/main/io)
    ADD_SUBDIRECTORY(readdy/main/plugin)
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
    INSTALL(TARGETS ${READDY_PROJECT_NAME}
            RUNTIME DESTINATION bin
            LIBRARY DESTINATION lib
            ARCHIVE DESTINATION lib)
    INSTALL_TARGETS(/lib RUNTIME_DIRECTORY /lib ${READDY_PROJECT_NAME})
ENDIF(READDY_BUILD_SHARED_COMBINED)
INSTALL(DIRECTORY include/ DESTINATION include)

#####################################
#                                   #
# readdy wrappers                   #
#                                   #
#####################################
IF (READDY_BUILD_PYTHON_WRAPPER)
    ADD_SUBDIRECTORY(wrappers)
ENDIF (READDY_BUILD_PYTHON_WRAPPER)

#####################################
#                                   #
# tests                             #
#                                   #
#####################################
IF (READDY_CREATE_TEST_TARGET)

    INCLUDE_DIRECTORIES(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})
    ADD_SUBDIRECTORY(readdy/test)
    ADD_TEST(NAME readdy-tests COMMAND runUnitTests)

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
