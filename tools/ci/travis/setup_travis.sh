#!/usr/bin/env bash
set -evx
env | sort
cd $HOME
if [[ "${TRAVIS_OS_NAME}" == "linux" ]]; then
  #CMAKE_URL="http://www.cmake.org/files/v3.3/cmake-3.3.2-Linux-x86_64.tar.gz"
  #mkdir -p _cmake
  #wget ${CMAKE_URL} --no-check-certificate -q -O _cmake/__cmake.tar.gz
  #for file in _cmake/*.tar.gz; do tar -xzf "${file}" --strip=1 -C _cmake && rm "${file}"; done
  echo "skip this since cmake gets provided by conda"
else
  echo "skip this"
  #brew update >/dev/null
  # brew unlink cmake
  # brew install cmake
  if [ "$CXX" = "g++" ]; then
    brew unlink gcc || true
    brew install gcc || true
  #    brew install boost --c++11 || true
  #    brew install boost-python --c++11 || true
  fi
fi