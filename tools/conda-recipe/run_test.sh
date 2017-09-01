#!/usr/bin/env bash

set +e

if [ $(uname) = "Linux" ]; then
    echo "installing logbt"
    mkdir ./mason/
    curl -sSfL https://github.com/mapbox/mason/archive/v0.2.0.tar.gz | tar --gunzip --extract --strip-components=1 --directory=./mason/
    ./mason/mason install gdb 7.12
    export PATH=$(./mason/mason prefix gdb 7.12)/bin:${PATH}
    curl -sSfL https://github.com/mapbox/logbt/archive/v1.6.0.tar.gz | tar --gunzip --extract --strip-components=2 --exclude="*md" --exclude="test*" --directory=.
    sudo bash -c "echo '/tmp/logbt-coredumps/core.%p.%E' > /proc/sys/kernel/core_pattern"
fi

# overall return code
ret_code=0
# current return code
err_code=0

export READDY_N_CORES=2

echo "calling c++ core unit tests"
runUnitTests
err_code=$?
if [ ${err_code} -ne 0 ]; then
    ret_code=${err_code}
    echo "core unit tests failed with ${ret_code}"
    if [ $(uname) = "Linux" ]; then
        echo "re-running with logbt"
        ./logbt runUnitTests
    fi
fi

echo "calling c++ singlecpu unit tests"
runUnitTests_singlecpu
err_code=$?
if [ ${err_code} -ne 0 ]; then
    ret_code=${err_code}
    echo "singlecpu unit tests failed with ${ret_code}"
    if [ $(uname) = "Linux" ]; then
        echo "re-running with logbt"
        ./logbt runUnitTests_singlecpu
    fi
fi

echo "calling c++ cpu unit tests"
runUnitTests_cpu
err_code=$?
if [ ${err_code} -ne 0 ]; then
    ret_code=${err_code}
    echo "cpu unit tests failed with ${ret_code}"
    if [ $(uname) = "Linux" ]; then
        echo "re-running with logbt"
        ./logbt runUnitTests_cpu
    fi
fi

echo "calling python unit tests"
nosetests readdy -s -vv --with-doctest --doctest-options=+NORMALIZE_WHITESPACE,+ELLIPSIS
err_code=$?
if [ ${err_code} -ne 0 ]; then
    ret_code=${err_code}
    echo "python unit tests failed with ${ret_code}"
fi

exit ${ret_code}
