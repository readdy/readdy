SET(Boost_VERSION "1_59_0")
SET(Boost_SHA1 "5123209db194d66d69a9cfa5af8ff473d5941d97")

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
ENDIF()