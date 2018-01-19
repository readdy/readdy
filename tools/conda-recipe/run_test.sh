#!/usr/bin/env bash

set +e

# overall return code
ret_code=0
# current return code
err_code=0

# cd $BUILD_PREFIX

export READDY_N_CORES=2

echo "calling python unit tests"
nosetests readdy -s -vv --with-doctest --doctest-options=+NORMALIZE_WHITESPACE,+ELLIPSIS
err_code=$?
if [ ${err_code} -ne 0 ]; then
    ret_code=${err_code}
    echo "python unit tests failed with ${ret_code}"
fi

exit ${ret_code}
