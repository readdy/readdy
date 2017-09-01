#!/bin/bash

unset MACOSX_DEPLOYMENT_TARGET

#########################################################
#                                                       #
# cmake configuration flags                             #
#                                                       #
#########################################################

# prefix path
CMAKE_FLAGS="-DCMAKE_INSTALL_PREFIX=$PREFIX"
# do not generate documentation target
CMAKE_FLAGS+=" -DREADDY_GENERATE_DOCUMENTATION_TARGET:BOOL=OFF"
# build monolithic lib
CMAKE_FLAGS+=" -DREADDY_BUILD_SHARED_COMBINED:BOOL=ON"
# release compile flags
CMAKE_FLAGS+=" -DCMAKE_BUILD_TYPE=${BUILD_TYPE}" # BUILD_TYPE="Release" if not a PR, otherwise RelWithDebInfo
# debug cmake config
CMAKE_FLAGS+=" -DREADDY_LOG_CMAKE_CONFIGURATION:BOOL=ON"
# enable testing and install test target
CMAKE_FLAGS+=" -DREADDY_CREATE_TEST_TARGET:BOOL=ON"
CMAKE_FLAGS+=" -DREADDY_INSTALL_UNIT_TEST_EXECUTABLE:BOOL=ON"
# hdf5 flags
CMAKE_FLAGS+=" -DHDF5_INCLUDE_DIR=$PREFIX/include"
# select compiler
CMAKE_FLAGS+=" -DCMAKE_C_COMPILER=${CC}"
CMAKE_FLAGS+=" -DCMAKE_CXX_COMPILER=${CXX}"

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

# cant reliably determine cpu count in a docker container,
# therefore fix this value.
if [ "$TRAVIS" == "true" ]; then CPU_COUNT=2; fi

mkdir build || true
cd build
echo "calling cmake with flags: "
for flag in ${CMAKE_FLAGS// /$'\n'}
do
    echo "   $flag"
done

cmake .. ${CMAKE_FLAGS}
make -j${CPU_COUNT}
make install &> /dev/null

if [ $(uname) = "Darwin" ]; then
    install_name_tool -add_rpath @loader_path/../readdy/readdy_plugins/ $PREFIX/bin/runUnitTests_singlecpu
    install_name_tool -add_rpath @loader_path/../readdy/readdy_plugins/ $PREFIX/bin/runUnitTests_cpu
fi