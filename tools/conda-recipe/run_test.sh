#!/usr/bin/env bash

set -evx

echo "calling c++ core unit tests"
runUnitTests

echo "calling c++ singlecpu unit tests"
runUnitTests_singlecpu

echo "calling python unit tests"
nosetests readdy -vv --with-doctest --doctest-options=+NORMALIZE_WHITESPACE,+ELLIPSIS
