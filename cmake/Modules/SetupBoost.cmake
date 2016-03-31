SET(Boost_VERSION "1_59_0")
#SET(Boost_SHA1 "da39a3ee5e6b4b0d3255bfef95601890afd80709")
IF(DEFINED ENV{BOOST_DOWNLOAD_CACHE_DIR})
    SET(BOOST_DOWNLOAD_CACHE_DIR "$ENV{BOOST_DOWNLOAD_CACHE_DIR}")
ELSE()
    SET(BOOST_DOWNLOAD_CACHE_DIR "${CMAKE_BINARY_DIR}/download/external/boost")
ENDIF()

#[[
BZIP2_FOUND - system has BZip2
BZIP2_INCLUDE_DIR - the BZip2 include directory
BZIP2_LIBRARIES - Link these to use BZip2
BZIP2_NEED_PREFIX - this is set if the functions are prefixed with BZ2_
BZIP2_VERSION_STRING - the version of BZip2 found (since CMake 2.8.8)
]]
FIND_PACKAGE(BZip2 REQUIRED)
LIST(GET BZIP2_LIBRARIES 0 BZIPLIB)
GET_FILENAME_COMPONENT(BZIP2_LIBRARY_DIR ${BZIPLIB} PATH)
MESSAGE(STATUS "BZip2 Library dir: ${BZIP2_LIBRARY_DIR}")
#[[
ZLIB_INCLUDE_DIRS   - where to find zlib.h, etc.
ZLIB_LIBRARIES      - List of libraries when using zlib.
ZLIB_FOUND          - True if zlib found.
]]
FIND_PACKAGE(ZLIB REQUIRED)
LIST(GET ZLIB_LIBRARIES 0 ZLIBLIB)
GET_FILENAME_COMPONENT(ZLIB_LIBRARY_DIR ${ZLIBLIB} PATH)
MESSAGE(STATUS "ZLIB library dir: ${ZLIB_LIBRARY_DIR}")

# download into BOOST_DOWNLOAD_CACHE_DIR
IF(NOT EXISTS "${BOOST_DOWNLOAD_CACHE_DIR}/")
    FILE(MAKE_DIRECTORY ${BOOST_DOWNLOAD_CACHE_DIR})
ENDIF()
SET(BOOST_DOWNLOAD_OUT "${BOOST_DOWNLOAD_CACHE_DIR}/boost_${Boost_VERSION}.tar.gz")
IF (NOT EXISTS "${BOOST_DOWNLOAD_OUT}")
    FIND_PROGRAM(CURL curl)
    IF(NOT CURL)
        MESSAGE(FATAL_ERROR "Could not locate curl!")
    ENDIF(NOT CURL)
    MESSAGE(STATUS "Saving boost to ${BOOST_DOWNLOAD_OUT}")
    EXECUTE_PROCESS(
            COMMAND ${CURL} "-kL" "http://sourceforge.net/projects/boost/files/boost/1.59.0/boost_${Boost_VERSION}.tar.gz/download" "-o" "${BOOST_DOWNLOAD_OUT}"
            RESULT_VARIABLE Result
            OUTPUT_VARIABLE Output
    )
    IF(NOT Result EQUAL "0")
        MESSAGE(FATAL_ERROR "Download failed: ${Output}")
    ENDIF(NOT Result EQUAL "0")
    #[[FILE(
            DOWNLOAD
            ${BOOST_DOWNLOAD_OUT}
            STATUS BOOST_DOWNLOAD_STATUS
            LOG BOOST_DL_LOG
            SHOW_PROGRESS
            EXPECTED_HASH SHA1=${Boost_SHA1}
    )
    MESSAGE(STATUS "Boost DL Log:\n${BOOST_DL_LOG}")]]
ENDIF (NOT EXISTS "${BOOST_DOWNLOAD_OUT}")
# unzip
SET(BOOST_UNZIP_OUT "${BOOST_DOWNLOAD_CACHE_DIR}/boost_${Boost_VERSION}")
IF (NOT EXISTS "${BOOST_UNZIP_OUT}/")
    EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E tar xzf "${BOOST_DOWNLOAD_OUT}"
            WORKING_DIRECTORY "${BOOST_DOWNLOAD_CACHE_DIR}")
ENDIF (NOT EXISTS "${BOOST_UNZIP_OUT}/")

# boostrap
UNSET(b2Path CACHE)
FIND_PROGRAM(b2Path NAMES bjam b2 PATHS ${BOOST_UNZIP_OUT} NO_DEFAULT_PATH)
IF (NOT b2Path)
    IF (MSVC)
        FIND_PROGRAM(b2Bootstrap NAMES bootstrap.bat PATHS ${BOOST_UNZIP_OUT} NO_DEFAULT_PATH)
    ELSE ()
        FIND_PROGRAM(b2Bootstrap NAMES bootstrap.sh PATHS ${BOOST_UNZIP_OUT} NO_DEFAULT_PATH)
    ENDIF ()
    MESSAGE(STATUS "Building b2 (bjam) with ${b2Bootstrap}")
    LIST(APPEND BOOTSTRAP_ARGS "--prefix=${CMAKE_BINARY_DIR}/out/boost")
    IF (NOT MSVC)
        IF (CMAKE_CXX_COMPILER_ID MATCHES "^(Apple)?Clang$")
            LIST(APPEND BOOTSTRAP_ARGS "--with-toolset=clang")
        ELSEIF (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
            LIST(APPEND BOOTSTRAP_ARGS "--with-toolset=gcc")
        ENDIF ()
    ENDIF (NOT MSVC)
    # execute bootstrap
    EXECUTE_PROCESS(COMMAND "${b2Bootstrap}" "${BOOTSTRAP_ARGS}"
            WORKING_DIRECTORY ${BOOST_UNZIP_OUT}
            RESULT_VARIABLE Result
            OUTPUT_VARIABLE Output
            ERROR_VARIABLE Error)
    IF (NOT Result EQUAL "0")
        MESSAGE(FATAL_ERROR "Failed running ${b2Bootstrap}:\n${Output}\n${Error}\n")
    ENDIF (NOT Result EQUAL "0")
    FIND_PROGRAM(b2Path NAMES bjam b2 PATHS ${BOOST_UNZIP_OUT} NO_DEFAULT_PATH)
ENDIF (NOT b2Path)

IF (READDY_BUILD_PYTHON_WRAPPER)
    FIND_PACKAGE(PythonInterp REQUIRED)
    EXECUTE_PROCESS(
            COMMAND "${PYTHON_EXECUTABLE}" "${READDY_GLOBAL_DIR}/libraries/boost/python_include_dir.py"
            RESULT_VARIABLE Result
            OUTPUT_VARIABLE Output
            OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    IF (NOT Result EQUAL "0")
        MESSAGE(FATAL_ERROR "Failed running python_include_dir script:\n${Output}")
    ENDIF (NOT Result EQUAL "0")
    # set python include dir as environment variable
    SET(ENV{PYTHON_INCLUDE_DIR} "${Output}")
    MESSAGE(STATUS "Found python include dir \"$ENV{PYTHON_INCLUDE_DIR}\"")
ENDIF (READDY_BUILD_PYTHON_WRAPPER)

LIST(APPEND B2ARGS
        "link=shared"
        "threading=multi"
        "runtime-link=shared"
        "dll-path=${BOOST_UNZIP_OUT}/stage/lib"
        "--build-dir=Build"
        "-sBZIP2_LIBPATH=${BZIP2_LIBRARY_DIR}"
        "-sBZIP2_INCLUDE=${BZIP2_INCLUDE_DIR}"
        "-sZLIB_LIBPATH=${ZLIB_LIBRARY_DIR}"
        "-sZLIB_INCLUDE=${ZLIB_INCLUDE_DIRS}"
        "--user-config=${READDY_GLOBAL_DIR}/libraries/boost/user-config.jam"
        "stage")

IF (APPLE)
    SET(BOOST_CXX_FLAGS "-mmacosx-version-min=10.7 -std=c++11 -stdlib=libc++")
    SET(BOOST_LINK_FLAGS "-mmacosx-version-min=10.7 -stdlib=libc++")
    SET(BOOST_CXX_FLAGS "cxxflags=${BOOST_CXX_FLAGS}")
    SET(BOOST_LINK_FLAGS "linkflags=${BOOST_LINK_FLAGS}")

    LIST(APPEND B2ARGS "toolset=clang")
ELSE ()
    SET(BOOST_CXX_FLAGS "")
    SET(BOOST_LINK_FLAGS "")

    LIST(APPEND B2ARGS "toolset=gcc")
ENDIF (APPLE)

IF (MSVC)
    IF (MSVC11)
        LIST(APPEND B2ARGS toolset=msvc-11.0)
    ELSEIF (MSVC12)
        LIST(APPEND B2ARGS toolset=msvc-12.0)
    ELSEIF (MSVC14)
        LIST(APPEND B2ARGS toolset=msvc-14.0)
    ENDIF ()
    LIST(APPEND B2ARGS define=_BIND_TO_CURRENT_MFC_VERSION=1 define=_BIND_TO_CURRENT_CRT_VERSION=1 --layout=versioned)
    IF (TargetArchitecture STREQUAL "x86_64")
        LIST(APPEND B2ARGS address-model=64)
    ENDIF (TargetArchitecture STREQUAL "x86_64")
ENDIF (MSVC)

INCLUDE(ProcessorCount)
ProcessorCount(CPUS)
LIST(APPEND B2ARGS "-j${CPUS}")
IF (BOOST_CXX_FLAGS)
    LIST(APPEND B2ARGS "${BOOST_CXX_FLAGS}")
ENDIF (BOOST_CXX_FLAGS)
IF (BOOST_LINK_FLAGS)
    LIST(APPEND B2ARGS "${BOOST_LINK_FLAGS}")
ENDIF (BOOST_LINK_FLAGS)

LIST(APPEND BOOST_MODULES "system" "filesystem" "thread" "log")
IF (READDY_BUILD_PYTHON_WRAPPER)
    LIST(APPEND BOOST_MODULES "python")
ENDIF (READDY_BUILD_PYTHON_WRAPPER)
FOREACH (component ${BOOST_MODULES})
    LIST(APPEND B2ARGS "--with-${component}")
ENDFOREACH (component ${BOOST_MODULES})

MESSAGE(STATUS "UNZIP OUT=${BOOST_UNZIP_OUT}")
EXECUTE_PROCESS(
        COMMAND ${b2Path} ${B2ARGS}
        WORKING_DIRECTORY ${BOOST_UNZIP_OUT}
        RESULT_VARIABLE Result
)
IF (NOT "${Result}" EQUAL "0")
    MESSAGE(FATAL_ERROR "Could not build boost!")
ENDIF (NOT "${Result}" EQUAL "0")
SET(BOOST_INCLUDEDIR "${BOOST_UNZIP_OUT}")
SET(BOOST_LIBRARYDIR "${BOOST_UNZIP_OUT}/stage/lib")
SET(Boost_NO_SYSTEM_PATHS ON)

MESSAGE("Boost includedir: ${BOOST_INCLUDEDIR}")
MESSAGE("Boost librarydir: ${BOOST_LIBRARYDIR}")
MESSAGE("Boost no system paths: ${Boost_NO_SYSTEM_PATHS}")

SET(Boost_INCLUDE_DIRS "${BOOST_UNZIP_OUT}/boost")

INSTALL(DIRECTORY "${BOOST_UNZIP_OUT}/stage/lib/" COMPONENT Development DESTINATION "lib/")
INSTALL(DIRECTORY ${BOOST_UNZIP_OUT}/boost COMPONENT Development DESTINATION include)
