#!/bin/bash

#########################################################
#                                                       #
# fix dylib rpaths, gracefully taken from               #
# https://github.com/minrk/patch-conda-rpaths           #
#                                                       #
#########################################################

#########################################################
#                                                       #
# cmake configuration flags                             #
#                                                       #
#########################################################

# prefix path
CMAKE_FLAGS="-DCMAKE_INSTALL_PREFIX=$PREFIX"
# boost
CMAKE_FLAGS+=" -DBOOST_DOWNLOAD_CACHE_DIR:PATH=${BOOST_DOWNLOAD_CACHE_DIR}"
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
CMAKE_FLAGS+=" -DREADDY_INSTALL_UNIT_TEST_EXECUTABLE:BOOL=ON"
# hdf5 flags
CMAKE_FLAGS+=" -DHDF5_INCLUDE_DIR=$PREFIX/include"

# attempt to feed the right python library to
# FindPythonLibs (weird behavior of cmake's find module)
lib_path=""
if [ `uname` == Darwin ]; then
    if [ ${PY3K} -eq 1 ]; then
        cd ${PREFIX}/lib
        ln -s libpython${PY_VER}m.dylib libpython${PY_VER}.dylib
        cd -
    fi
    lib_path="$PREFIX/lib/libpython${PY_VER}.dylib"
fi
if [ `uname` == Linux ]; then
    if [ ${PY3K} -eq 1 ]; then
        lib_path="$PREFIX/lib/libpython${PY_VER}m.so"
    else
        lib_path="$PREFIX/lib/libpython${PY_VER}.so"
    fi
fi
CMAKE_FLAGS+=" -DPYTHON_LIBRARY:FILEPATH=${lib_path}"

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

export HDF5_DIR=${PREFIX}
# export HDF5_ROOT=${PREFIX}
export PYTHON_INCLUDE_DIR=`python -c "from __future__ import print_function; import distutils.sysconfig; print(distutils.sysconfig.get_python_inc(True))"`
if [ `uname` == Darwin ]; then
    export MACOSX_DEPLOYMENT_TARGET=10.9
fi

# cant reliably determine cpu count in a docker container,
# therefore fix this value.
if [ "$TRAVIS" == "true" ]; then CPU_COUNT=2; fi

mkdir build
cd build
echo "calling cmake with flags: "
for flag in ${CMAKE_FLAGS// /$'\n'}
do
    echo "   $flag"
done
cmake .. ${CMAKE_FLAGS}
make -j${CPU_COUNT}
make install &> /dev/null