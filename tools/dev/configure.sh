#!/bin/bash

# config
PREFIX="~/miniconda2"
PY3K=0
PY_VER="2.7"

#########################################################
#                                                       #
# cmake configuration flags                             #
#                                                       #
#########################################################

# prefix path
CMAKE_FLAGS="-DCMAKE_INSTALL_PREFIX=$PREFIX"
# do not generate documentation target
CMAKE_FLAGS+=" -DREADDY_GENERATE_DOCUMENTATION_TARGET:BOOL=ON"
# build separate libs
CMAKE_FLAGS+=" -DREADDY_BUILD_SHARED_COMBINED:BOOL=OFF"
# release compile flags
CMAKE_FLAGS+=" -DCMAKE_BUILD_TYPE=Debug"
# debug cmake config
CMAKE_FLAGS+=" -DREADDY_LOG_CMAKE_CONFIGURATION:BOOL=ON"
# enable testing and install test target
CMAKE_FLAGS+=" -DREADDY_CREATE_TEST_TARGET:BOOL=ON"
CMAKE_FLAGS+=" -DREADDY_INSTALL_UNIT_TEST_EXECUTABLE:BOOL=OFF"

if [ $1 = "clang" ]
then
    CMAKE_FLAGS+=" -DCMAKE_C_COMPILER=/usr/bin/clang"
    CMAKE_FLAGS+=" -DCMAKE_CXX_COMPILER=/usr/bin/clang++"
fi

export HDF5_ROOT=${PREFIX}
export PYTHON_INCLUDE_DIR=`python -c "from __future__ import print_function; import distutils.sysconfig; print(distutils.sysconfig.get_python_inc(True))"`

# attempt to feed the right python library to
# FindPythonLibs (weird behavior of cmake's find module)
lib_path=""
if [ `uname` == Darwin ]; then
    if ! [ $PY3K eq 1 ]; then
        cd ${PREFIX}/lib
        ln -s libpython${PY_VER}m.dylib libpython${PY_VER}.dylib
        cd -
    fi
    lib_path="$PREFIX/lib/libpython${PY_VER}.dylib"
fi
if [ `uname` == Linux ]; then
    if [ $PY3K -eq 1 ]; then
        lib_path="$PREFIX/lib/libpython${PY_VER}m.so"
    else
        lib_path="$PREFIX/lib/libpython${PY_VER}.so"
    fi
fi
CMAKE_FLAGS+=" -DPYTHON_LIBRARY:FILEPATH=${lib_path}"

# cant reliably determine cpu count in a docker container,
# therefore fix this value.
if [ "$TRAVIS" == "true" ]; then CPU_COUNT=2; fi

mkdir ~/Development/readdy/build || true
cd ~/Development/readdy/build
echo "calling cmake with flags: "
for flag in ${CMAKE_FLAGS// /$'\n'}
do
    echo "   $flag"
done
cmake ~/Development/readdy/readdy ${CMAKE_FLAGS}
cd -