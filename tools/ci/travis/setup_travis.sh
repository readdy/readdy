#!/usr/bin/env bash

if [[ "${TRAVIS_OS_NAME}" == "linux" ]]; then
  CMAKE_URL="http://www.cmake.org/files/v3.3/cmake-3.3.2-Linux-x86_64.tar.gz"
  mkdir _cmake && wget --quiet ${CMAKE_URL} -O - | tar --strip-components=1 -xz -C _cmake
  export PATH=${DEPS_DIR}/_cmake/bin:${PATH}
else
  brew update
  brew install cmake
fi