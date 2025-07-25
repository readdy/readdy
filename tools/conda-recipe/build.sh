#!/bin/bash

mkdir -p build
cd build || true
rm -rf ./*

echo "conda recipe build.sh ..."
echo "GIT_DESCRIBE_NUMBER ${GIT_DESCRIBE_NUMBER}"
echo "PY_VER ${PY_VER}"
echo "GIT_BUILD_STR ${GIT_BUILD_STR}"
echo "Install prefix ${SP_DIR}"

export HDF5_ROOT=${PREFIX}

echo "Running cmake"

# On macOS, force C++17 to avoid C++20/SDK compatibility issues
if [[ "$OSTYPE" == "darwin"* ]]; then
  CXX_STANDARD="-DCMAKE_CXX_STANDARD=17"
else
  CXX_STANDARD=""
fi

cmake .. \
  -DCMAKE_INSTALL_PREFIX="${PREFIX}" \
  -DCMAKE_PREFIX_PATH="${PREFIX}" \
  -DCMAKE_OSX_SYSROOT="${CONDA_BUILD_SYSROOT}" \
  -DCMAKE_BUILD_TYPE=Release \
  -DPYTHON_EXECUTABLE="${PYTHON}" \
  -DPYTHON_PREFIX="${PREFIX}" \
  -DHDF5_INCLUDE_DIRS="${PREFIX}/include" \
  -DREADDY_LOG_CMAKE_CONFIGURATION:BOOL=ON \
  -DREADDY_CREATE_TEST_TARGET:BOOL=ON \
  -DREADDY_INSTALL_UNIT_TEST_EXECUTABLE:BOOL=OFF \
  -DREADDY_VERSION=${PKG_VERSION} \
  -DREADDY_BUILD_STRING=${PKG_BUILDNUM} \
  -DSP_DIR="${SP_DIR}" \
  -DIGNORE_VENDORED_DEPENDENCIES:BOOL=ON \
  ${CXX_STANDARD} \
  -GNinja

echo "Running ninja with makeflags ${MAKEFLAGS}"
ninja ${MAKEFLAGS}
echo "Running ninja install"
ninja install

export READDY_N_CORES=2
export READDY_PLUGIN_DIR="${PREFIX}/lib/readdy_plugins"

err_code=0
ret_code=0

echo "calling c++ core unit tests"
./bin/runUnitTests --durations yes
err_code=$?
if [ ${err_code} -ne 0 ]; then
   ret_code=${err_code}
   echo "core unit tests failed with ${ret_code}"
fi

# echo "calling c++ integration tests"
# ./readdy/test/runUnitTests --durations yes [integration]
# err_code=$?
# if [ ${err_code} -ne 0 ]; then
#    ret_code=${err_code}
#    echo "core unit tests failed with ${ret_code}"
# fi

echo "calling c++ singlecpu unit tests"
./bin/runUnitTests_singlecpu --durations yes
err_code=$?
if [ ${err_code} -ne 0 ]; then
   ret_code=${err_code}
   echo "singlecpu unit tests failed with ${ret_code}"
fi

echo "calling c++ cpu unit tests"
./bin/runUnitTests_cpu --durations yes
err_code=$?
if [ ${err_code} -ne 0 ]; then
  ret_code=${err_code}
  echo "cpu unit tests failed with ${ret_code}"
fi

exit ${ret_code}
