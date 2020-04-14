#!/bin/bash

unset MACOSX_DEPLOYMENT_TARGET

echo "conda recipe build.sh ..."
echo "GIT_DESCRIBE_NUMBER ${GIT_DESCRIBE_NUMBER}"
echo "PY_VER ${PY_VER}"
echo "GIT_BUILD_STR ${GIT_BUILD_STR}"
echo "GIT_STUB ${GIT_STUB}"

#########################################################
#                                                       #
# cmake configuration flags                             #
#                                                       #
#########################################################

# install prefix
CMAKE_FLAGS="-DCMAKE_INSTALL_PREFIX=${PREFIX}"
# prefix path
CMAKE_FLAGS+=" -DCMAKE_PREFIX_PATH=${PREFIX}"
# python prefix
CMAKE_FLAGS+=" -DPYTHON_PREFIX=${PREFIX}"
# python executable
CMAKE_FLAGS+=" -DPYTHON_EXECUTABLE=${PYTHON}"
# do not generate documentation target
CMAKE_FLAGS+=" -DREADDY_GENERATE_DOCUMENTATION_TARGET:BOOL=OFF"
# build monolithic lib
CMAKE_FLAGS+=" -DREADDY_BUILD_SHARED_COMBINED:BOOL=ON"
# release compile flags
CMAKE_FLAGS+=" -DCMAKE_BUILD_TYPE=Release"
# debug cmake config
CMAKE_FLAGS+=" -DREADDY_LOG_CMAKE_CONFIGURATION:BOOL=ON"
# enable testing and install test target
CMAKE_FLAGS+=" -DREADDY_CREATE_TEST_TARGET:BOOL=ON"
CMAKE_FLAGS+=" -DREADDY_INSTALL_UNIT_TEST_EXECUTABLE:BOOL=OFF"
# hdf5 flags
CMAKE_FLAGS+=" -DHDF5_INCLUDE_DIRS=${PREFIX}/include"
# version flags
CMAKE_FLAGS+=" -DREADDY_VERSION=${PKG_VERSION}"
CMAKE_FLAGS+=" -DREADDY_BUILD_STRING=${PKG_BUILDNUM}"

#########################################################
#                                                       #
# environment variables                                 #
#   - HDF5_ROOT: is set to prefix (the find module      #
#     seems not to respect the prefix path of cmake)    #
#   - PYTHON_INCLUDE_DIR: is set w.r.t. whatever python #
#     interpreter was found first in PATH               #
#   - MACOSX_DEPLOYMENT_TARGET: fix the osx sdk         #
#                                                       #
#########################################################

export HDF5_ROOT=${PREFIX}

mkdir build || true
cd build
echo "calling cmake with flags: "
for flag in ${CMAKE_FLAGS// /$'\n'}
do
    echo "   $flag"
done

cmake .. ${CMAKE_FLAGS}
make ${MAKEFLAGS} -v
make install

export READDY_N_CORES=2

err_code=0
ret_code=0

echo "calling c++ core unit tests"
CONDA_ENV_PATH=${PREFIX} bin/runUnitTests --durations yes
err_code=$?
if [ ${err_code} -ne 0 ]; then
    ret_code=${err_code}
    echo "core unit tests failed with ${ret_code}"
fi

echo "calling c++ integration tests"
CONDA_ENV_PATH=${PREFIX} bin/runUnitTests --durations yes [integration]
err_code=$?
if [ ${err_code} -ne 0 ]; then
    ret_code=${err_code}
    echo "core unit tests failed with ${ret_code}"
fi

echo "calling c++ singlecpu unit tests"
CONDA_ENV_PATH=${PREFIX} bin/runUnitTests_singlecpu --durations yes
err_code=$?
if [ ${err_code} -ne 0 ]; then
    ret_code=${err_code}
    echo "singlecpu unit tests failed with ${ret_code}"
fi

echo "calling c++ cpu unit tests"
CONDA_ENV_PATH=${PREFIX} bin/runUnitTests_cpu --durations yes
err_code=$?
if [ ${err_code} -ne 0 ]; then
    ret_code=${err_code}
    echo "cpu unit tests failed with ${ret_code}"
fi

exit ${ret_code}
