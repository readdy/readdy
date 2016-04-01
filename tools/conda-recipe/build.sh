#!/bin/bash

CMAKE_FLAGS="-DCMAKE_INSTALL_PREFIX=$PREFIX -DBUILD_TESTING=OFF"
CMAKE_FLAGS+=" -DREADDY_BUILD_SHARED_COMBINED:BOOL=ON"
CMAKE_FLAGS+=" -DCMAKE_BUILD_TYPE=Release"

export PYTHON_INCLUDE_DIR=`python -c "from __future__ import print_function; import distutils.sysconfig; print(distutils.sysconfig.get_python_inc(True))"`

# cant reliably determine cpu count in a docker container,
# therefore fix this value.
if [ "$TRAVIS" == "true" ]; then CPU_COUNT=2; fi

mkdir build
cd build
cmake .. ${CMAKE_FLAGS}
make -j${CPU_COUNT}
make -j${CPU_COUNT} install &> /dev/null