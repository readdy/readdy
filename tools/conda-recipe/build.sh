#!/bin/bash

CMAKE_FLAGS="-DCMAKE_INSTALL_PREFIX=$PREFIX -DBUILD_TESTING=OFF"
CMAKE_FLAGS+=" -DREADDY_BUILD_SHARED_COMBINED:BOOL=ON"
CMAKE_FLAGS+=" -DCMAKE_BUILD_TYPE=Release"

# cant reliably determine cpu count in a docker container,
# therefore fix this value.
if [ "$TRAVIS" == "true" ]; then CPU_COUNT=2; fi

mkdir build
cd build
cmake .. ${CMAKE_FLAGS}
make -j${CPU_COUNT}
make -j${CPU_COUNT} install