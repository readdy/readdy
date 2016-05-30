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