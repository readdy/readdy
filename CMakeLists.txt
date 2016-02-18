#-----------------------------------
# ReaDDy2
#-----------------------------------
PROJECT(ReaDDy2)
CMAKE_MINIMUM_REQUIRED(VERSION 3.2.0)

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

# list of libraries for easier linking
SET(READDY_LIBRARIES "")

FIND_PACKAGE(Boost COMPONENTS system thread REQUIRED)

# decide if api wrappers are being generated
SET(READDY_BUILD_PYTHON_WRAPPER ON CACHE BOOL "Build Python wrappers")
IF(READDY_BUILD_PYTHON_WRAPPER)
	# add READDY_WITH_PYTHON define
    ADD_DEFINITIONS(-DREADDY_WITH_PYTHON)
	
	SET(BUILD_PY3 ON CACHE BOOL "True if should be built with python3, otherwise python2")
	IF(BUILD_PY3)
	    SET(Python_ADDITIONAL_VERSIONS 3.5)
	    MESSAGE(STATUS "--------------------- Building with python 3 ---------------------")
	ELSE()
	    SET(Python_ADDITIONAL_VERSIONS 2.7)
	    MESSAGE(STATUS "--------------------- Building with python 2 ---------------------")
	ENDIF()

	# macro that takes a list of potential library names
	MACRO(boost_python_libraries libnames)
	    FOREACH(_boost_py_lib ${libnames})
		MESSAGE(STATUS "Looking for boost lib: ${_boost_py_lib}")
	        FIND_PACKAGE(Boost QUIET COMPONENTS ${_boost_py_lib})
	        STRING(TOUPPER ${_boost_py_lib} boost_py_lib_name)
	        IF(Boost_${boost_py_lib_name}_FOUND)
	            SET(Boost_PYTHON_LIBRARIES ${Boost_${boost_py_lib_name}_LIBRARIES})
	            BREAK()
			ENDIF()
		ENDFOREACH()
	ENDMACRO(boost_python_libraries)
	# import python
	FIND_PACKAGE(PythonInterp)
	FIND_PACKAGE(PythonLibs)
	# look for boost python
	IF(BUILD_PY3)
	    # find boost-python3 lib
	    boost_python_libraries("python3;python-py35;python-py34;python")
	ELSE()
	    # find boost-python2 lib 
	    boost_python_libraries("python2;python-py27;python")
	ENDIF()
	# if no boost-python, exit
	IF(NOT Boost_PYTHON_LIBRARIES)
	    MESSAGE(FATAL_ERROR "Did not find boost python libraries, exit!" )
	ENDIF()

	MESSAGE(STATUS "Python libraries: ${PYTHON_LIBRARIES}")
	MESSAGE(STATUS "Boost libraries: ${Boost_LIBRARIES}")
	MESSAGE(STATUS "Python version: ${PYTHON_VERSION_STRING}")
	MESSAGE(STATUS "Boost python library: ${Boost_PYTHON_LIBRARIES}")
	
	INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_DIRS})
	LIST(APPEND READDY_LIBRARIES ${PYTHON_LIBRARIES})
ENDIF()
# library includes
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
# includes
INCLUDE_DIRECTORIES("${PROJECT_SOURCE_DIR}/include")

LIST(APPEND READDY_LIBRARIES ${Boost_LIBRARIES})

# this is where the sources are
ADD_SUBDIRECTORY(readdy2/main)

#
# tests
#
IF(EXISTS "${CMAKE_SOURCE_DIR}/libraries/googletest/CMakeLists.txt" AND READDY_RUN_TESTS)
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
    ADD_DEPENDENCIES(${READDY_LIBRARY_NAME} googletest)
    ENABLE_TESTING()
	INCLUDE_DIRECTORIES(${GOOGLETEST_INCLUDE})
    ADD_SUBDIRECTORY(readdy2/test)
ELSE()
    MESSAGE(STATUS "The submodule googletest was not initialized or READDY_RUN_TESTS was OFF, thus no testing.")
ENDIF()
