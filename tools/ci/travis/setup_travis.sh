#!/usr/bin/env bash
set -evx
cd $HOME
if ! [[ "${TRAVIS_OS_NAME}" == "linux" ]]; then

  # if we are to use g++ on darwin, make sure it is installed in the latest (homebrew) version
  if [ "$CXX" = "g++" ]; then
    brew unlink gcc || true
    brew install gcc || true
  fi

fi
cd -
if [[ "${TRAVIS_OS_NAME}" == "linux" ]]; then
    export LLVM_VERSION=3.7.0
    export LLVM_ARCHIVE_PATH=$HOME/clang+llvm.tar.xz
    wget http://llvm.org/releases/$LLVM_VERSION/clang+llvm-$LLVM_VERSION-x86_64-linux-gnu-ubuntu-14.04.tar.xz -O $LLVM_ARCHIVE_PATH
    mkdir $HOME/clang+llvm
    tar xf $LLVM_ARCHIVE_PATH -C $HOME/clang+llvm --strip-components 1
    export PATH=$HOME/clang+llvm/bin:$PATH
fi