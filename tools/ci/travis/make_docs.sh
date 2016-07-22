#!/usr/bin/env bash

# exit upon failed command and use of undeclared variables
set -e -u

function set_this_up {
    if [ "$TRAVIS_PULL_REQUEST" != "false" ]
    then
        echo "This was a pull request, thus dont build docs. Exit."
        exit 0
    fi
    if [ "$TRAVIS_BRANCH" != "master" ]
    then
        echo "This commit was made against the $TRAVIS_BRANCH branch and not the master branch. Exit."
        exit 0
    fi
    if [ "$CONDA_PY" != "27" ]
    then
        echo "Only build documentation for python version 2.7. Exit."
        exit 0
    fi
    if [ "${TRAVIS_OS_NAME}" != "linux" ]
    then
        echo "Only build documentation for linux. Exit."
        exit 0
    fi
    if [ -z ${GH_TOKEN+x} ]
    then
      echo "GH_TOKEN was not set, so this is probably a fork. Exit."
      exit 0
    fi

    # install doxypypy
    yes | pip install doxypypy
}
function make_doc {
    mkdir $HOME/_readdy_docs || true
    cd $HOME/_readdy_docs
    # cant reliably determine cpu count in a docker container,
    # therefore fix this value.
    if [ "$TRAVIS" == "true" ]; then CPU_COUNT=2; fi
    cmake ${TRAVIS_BUILD_DIR} -DREADDY_GENERATE_DOCUMENTATION_TARGET_ONLY:BOOL=ON
    make doc
    cd -
}

function setup_docs_repo {
    cd $HOME/_readdy_docs/docs/html
    git init
    git config user.name "Moritz Hoffmann"
    git config user.email "clonker@users.noreply.github.com"
    git remote add upstream "https://$GH_TOKEN@github.com/readdy/readdy.git"
    git fetch upstream
    git checkout --orphan workbranch
    git reset --hard
    cd -
}

function deploy {
    # revision tag
    cd ${TRAVIS_BUILD_DIR}
    rev=$(git rev-parse --short HEAD)
    cd $HOME/_readdy_docs/docs/html
    touch .
    git add -A .
    git commit -m "github pages"
    git push -q -f upstream workbranch:gh-pages
}

set_this_up

make_doc
setup_docs_repo

deploy
