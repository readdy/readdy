SET(Boost_VERSION "1_59_0")
SET(Boost_SHA1 "5123209db194d66d69a9cfa5af8ff473d5941d97")

#[[
BZIP2_FOUND - system has BZip2
BZIP2_INCLUDE_DIR - the BZip2 include directory
BZIP2_LIBRARIES - Link these to use BZip2
BZIP2_NEED_PREFIX - this is set if the functions are prefixed with BZ2_
BZIP2_VERSION_STRING - the version of BZip2 found (since CMake 2.8.8)
]]
FIND_PACKAGE(BZip2 REQUIRED)
#[[
ZLIB_INCLUDE_DIRS   - where to find zlib.h, etc.
ZLIB_LIBRARIES      - List of libraries when using zlib.
ZLIB_FOUND          - True if zlib found.
]]
FIND_PACKAGE(ZLIB REQUIRED)

# download into libraries/external/boost (so it gets ignored by git)
SET(BOOST_DOWNLOAD_OUT "${READDY_GLOBAL_DIR}/libraries/external/boost/boost_${Boost_VERSION}.tar.gz")
IF(NOT EXISTS "${BOOST_DOWNLOAD_OUT}")
    FILE(DOWNLOAD http://sourceforge.net/projects/boost/files/boost/1.59.0/boost_${Boost_VERSION}.tar.gz/download
            ${BOOST_DOWNLOAD_OUT}
            STATUS BOOST_DOWNLOAD_STATUS
            SHOW_PROGRESS
            EXPECTED_HASH SHA1=${Boost_SHA1})
ENDIF(NOT EXISTS "${BOOST_DOWNLOAD_OUT}")
# unzip
SET(BOOST_UNZIP_OUT "${READDY_GLOBAL_DIR}/libraries/external/boost/boost_${Boost_VERSION}")
IF(NOT EXISTS "${BOOST_UNZIP_OUT}/")
    EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E tar xzf "${BOOST_DOWNLOAD_OUT}"
            WORKING_DIRECTORY "${READDY_GLOBAL_DIR}/libraries/external/boost")
ENDIF(NOT EXISTS "${BOOST_UNZIP_OUT}/")

# boostrap
UNSET(b2Path CACHE)
FIND_PROGRAM(b2Path NAMES b2 PATHS ${BOOST_UNZIP_OUT} NO_DEFAULT_PATH)
IF(NOT b2Path)
    MESSAGE(STATUS "Building b2 (bjam)")
    IF(MSVC)
        SET(b2Bootstrap "bootstrap.bat")
    ELSE(MSVC)
        SET(b2Bootstrap "./bootstrap.sh")
        IF(CMAKE_CXX_COMPILER_ID MATCHES "^(Apple)?Clang$")
            list(APPEND b2Bootstrap --with-toolset=clang)
        ELSEIF(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
            list(APPEND b2Bootstrap --with-toolset=gcc)
        ENDIF()
    ENDIF(MSVC)
    EXECUTE_PROCESS(COMMAND ${b2Bootstrap}
            WORKING_DIRECTORY ${BOOST_UNZIP_OUT}
            RESULT_VARIABLE Result
            OUTPUT_VARIABLE Output
            ERROR_VARIABLE Error)
    IF(NOT Result EQUAL "0")
        MESSAGE(FATAL_ERROR "Failed running ${b2Bootstrap}:\n${Output}\n${Error}\n")
    ENDIF(NOT Result EQUAL "0")
ENDIF(NOT b2Path)

