#!/usr/bin/env bash

set -evx

source activate _test
if ! [[ "${TRAVIS_OS_NAME}" == "linux" ]]; then
    PREFIX=$HOME/miniconda/envs/_test runUnitTests
    PREFIX=$HOME/miniconda/envs/_test DYLD_LIBRARY_PATH=$HOME/miniconda/envs/_test/lib/readdy_plugins runUnitTests_singlecpu
fi

if [[ "${TRAVIS_OS_NAME}" == "linux" ]]; then
    PREFIX=$HOME/miniconda/envs/_test runUnitTests
    PREFIX=$HOME/miniconda/envs/_test LD_LIBRARY_PATH=$HOME/miniconda/envs/_test/lib/readdy_plugins runUnitTests_singlecpu
fi

source deactivate