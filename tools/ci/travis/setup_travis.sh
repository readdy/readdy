#!/usr/bin/env bash

if [[ "${TRAVIS_OS_NAME}" == "linux" ]]; then
  CMAKE_URL="http://www.cmake.org/files/v3.3/cmake-3.3.2-Linux-x86_64.tar.gz"
  mkdir -p _cmake
  wget ${CMAKE_URL} -O _cmake/__cmake.tar.gz
  for file in _cmake/*.tar.gz; do tar xzvf "${file}" -C _cmake && rm "${file}"; done
  export PATH=${DEPS_DIR}/_cmake/bin:${PATH}
else
  brew update >/dev/null
  brew unlink cmake
  brew install cmake
fi