#!/usr/bin/env bash

if false
then
    set -evx

    source activate _test
    if ! [[ "${TRAVIS_OS_NAME}" == "linux" ]]; then
        runUnitTests
        export DYLD_LIBRARY_PATH=${CONDA_ENV_PATH}/lib/readdy_plugins
        runUnitTests_singlecpu
    fi

    if [[ "${TRAVIS_OS_NAME}" == "linux" ]]; then
        runUnitTests
        export LD_LIBRARY_PATH=${CONDA_ENV_PATH}/lib/readdy_plugins
        runUnitTests_singlecpu
    fi

    source deactivate
fi