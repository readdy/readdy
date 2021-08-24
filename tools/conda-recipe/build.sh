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

cmake .. \
  -DCMAKE_INSTALL_PREFIX=${SP_DIR} \
  -DCMAKE_PREFIX_PATH=${PREFIX} \
  -DCMAKE_BUILD_TYPE=Release \
  -DPYTHON_EXECUTABLE=${PYTHON} \
  -DPYTHON_PREFIX=${PREFIX} \
  -DHDF5_INCLUDE_DIRS="${PREFIX}/include" \
  -DREADDY_BUILD_SHARED_COMBINED:BOOL=ON \
  -DREADDY_LOG_CMAKE_CONFIGURATION:BOOL=ON \
  -DREADDY_CREATE_TEST_TARGET:BOOL=ON \
  -DREADDY_INSTALL_UNIT_TEST_EXECUTABLE:BOOL=OFF \
  -DREADDY_VERSION=${PKG_VERSION} \
  -DREADDY_BUILD_STRING=${PKG_BUILDNUM} \
  -GNinja

ninja ${MAKEFLAGS}
ninja install

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
