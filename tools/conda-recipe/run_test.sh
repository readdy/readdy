#!/usr/bin/env bash

set +e

# overall return code
ret_code=0
# current return code
err_code=0

echo "calling c++ core unit tests"
runUnitTests
err_code=$?
if [ ${err_code} -ne 0 ]; then
    ret_code=${err_code}
    echo "core unit tests failed with ${ret_code}"
fi


echo "calling c++ singlecpu unit tests"
runUnitTests_singlecpu
err_code=$?
if [ ${err_code} -ne 0 ]; then
    ret_code=${err_code}
    echo "singlecpu unit tests failed with ${ret_code}"
fi

echo "calling c++ cpu unit tests"
runUnitTests_cpu
err_code=$?
if [ ${err_code} -ne 0 ]; then
    ret_code=${err_code}
    echo "cpu unit tests failed with ${ret_code}"
fi

echo "calling c++ cpu_dense unit tests"
runUnitTests_cpu_dense
err_code=$?
if [ ${err_code} -ne 0 ]; then
    ret_code=${err_code}
    echo "cpu_dense unit tests failed with ${ret_code}"
fi



echo "calling python unit tests"
nosetests readdy -vv --with-doctest --doctest-options=+NORMALIZE_WHITESPACE,+ELLIPSIS
err_code=$?
if [ ${err_code} -ne 0 ]; then
    ret_code=${err_code}
    echo "python unit tests failed with ${ret_code}"
fi

exit ${ret_code}
