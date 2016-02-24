#-----------------------------------
# ReaDDy2
#-----------------------------------
CMAKE_MINIMUM_REQUIRED(VERSION 3.2.0)

PROJECT(ReaDDy2_Main C CXX)

# for custom cmake modules
SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules")

IF("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")
    MESSAGE(SEND_ERROR "In-source builds are not allowed.")
ENDIF()

SET(CMAKE_VERBOSE_MAKEFILE ON)
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

# Compile in debug mode by default
IF(NOT CMAKE_BUILD_TYPE)
    MESSAGE(STATUS "Setting build type to 'DEBUG' as none was specified.")
    SET(CMAKE_BUILD_TYPE DEBUG)
ELSE(NOT CMAKE_BUILD_TYPE)
    MESSAGE(STATUS "Using specified '${CMAKE_BUILD_TYPE}' build type.")
ENDIF(NOT CMAKE_BUILD_TYPE)

# list of libraries for easier linking
SET(READDY_LIBRARIES "")
SET(READDY_HEADERS "")

FIND_PACKAGE(Boost COMPONENTS system thread REQUIRED)

# decide if api wrappers are being generated
SET(READDY_BUILD_PYTHON_WRAPPER ON CACHE BOOL "Build Python wrappers")
IF (READDY_BUILD_PYTHON_WRAPPER)
    # add READDY_WITH_PYTHON define
    ADD_DEFINITIONS(-DREADDY_WITH_PYTHON)

    # import python
    FIND_PACKAGE(PythonInterp REQUIRED)
    FIND_PACKAGE(PythonLibs REQUIRED)

    IF (${PYTHON_VERSION_STRING} GREATER 3.0)
        MESSAGE(STATUS "--------------------- Building with python 3 ---------------------")
    ELSE ()
        MESSAGE(STATUS "--------------------- Building with python 2 ---------------------")
    ENDIF ()

    # macro that takes a list of potential library names
    MACRO(boost_python_libraries libnames)
        FOREACH (_boost_py_lib ${libnames})
            MESSAGE(STATUS "Looking for boost lib: ${_boost_py_lib}")
            FIND_PACKAGE(Boost QUIET COMPONENTS ${_boost_py_lib})
            STRING(TOUPPER ${_boost_py_lib} boost_py_lib_name)
            IF (Boost_${boost_py_lib_name}_FOUND)
                SET(Boost_PYTHON_LIBRARIES ${Boost_${boost_py_lib_name}_LIBRARIES})
                BREAK()
            ENDIF ()
        ENDFOREACH ()
    ENDMACRO(boost_python_libraries)

    # look for boost python
    IF(${PYTHON_VERSION_STRING} GREATER 3.0)
        # find boost-python3 lib
        boost_python_libraries("python3;python-py35;python-py34;python")
    ELSE ()
        # find boost-python2 lib
        boost_python_libraries("python2;python-py27;python")
    ENDIF ()
    # if no boost-python, exit
    IF (NOT Boost_PYTHON_LIBRARIES)
        MESSAGE(FATAL_ERROR "Did not find boost python libraries, exit!")
    ENDIF ()

    MESSAGE(STATUS "Python libraries: ${PYTHON_LIBRARIES}")
    MESSAGE(STATUS "Boost libraries: ${Boost_LIBRARIES}")
    MESSAGE(STATUS "Python version: ${PYTHON_VERSION_STRING}")
    MESSAGE(STATUS "Boost python library: ${Boost_PYTHON_LIBRARIES}")

    LIST(APPEND READDY_LIBRARIES ${PYTHON_LIBRARIES})
ENDIF ()
# includes
SET(READDY_GLOBAL_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/include")
LIST(APPEND READDY_HEADERS "${PROJECT_SOURCE_DIR}/include")

# append libraries
LIST(APPEND READDY_LIBRARIES ${Boost_LIBRARIES} ${Boost_PYTHON_LIBRARIES})

# include directories from least dependent to most dependent
# i.e.: A depends on B -> add subdirectory of B first
ADD_SUBDIRECTORY(readdy2/main/io)
ADD_SUBDIRECTORY(readdy2/main)

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

#
# wrappers
#
IF(READDY_BUILD_PYTHON_WRAPPER)
    ADD_SUBDIRECTORY(wrappers)
ENDIF()