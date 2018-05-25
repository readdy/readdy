#!/bin/bash

unset MACOSX_DEPLOYMENT_TARGET

#########################################################
#                                                       #
# cmake configuration flags                             #
#                                                       #
#########################################################

# prefix path
CMAKE_FLAGS="-DCMAKE_INSTALL_PREFIX=${PREFIX}"
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
make ${MAKEFLAGS}
make install

#if [ $(uname) = "Darwin" ]; then
#    install_name_tool -add_rpath @loader_path/../readdy/readdy_plugins/ $BUILD_PREFIX/bin/runUnitTests_singlecpu
#    install_name_tool -add_rpath @loader_path/../readdy/readdy_plugins/ $BUILD_PREFIX/bin/runUnitTests_cpu
#fi

export READDY_N_CORES=2

err_code=0
ret_code=0

echo "calling c++ core unit tests"
CONDA_ENV_PATH=${PREFIX} bin/runUnitTests
err_code=$?
if [ ${err_code} -ne 0 ]; then
    ret_code=${err_code}
    echo "core unit tests failed with ${ret_code}"
fi

echo "calling c++ singlecpu unit tests"
CONDA_ENV_PATH=${PREFIX} bin/runUnitTests_singlecpu
err_code=$?
if [ ${err_code} -ne 0 ]; then
    ret_code=${err_code}
    echo "singlecpu unit tests failed with ${ret_code}"
fi

echo "calling c++ cpu unit tests"
CONDA_ENV_PATH=${PREFIX} bin/runUnitTests_cpu
err_code=$?
if [ ${err_code} -ne 0 ]; then
    ret_code=${err_code}
    echo "cpu unit tests failed with ${ret_code}"
fi

exit ${ret_code}
