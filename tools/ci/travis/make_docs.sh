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

    # install ruby packages
    gem install jekyll bundler
    # install doxypypy
    yes | pip install doxypypy
}

function make_reference_doc {
    mkdir -p $HOME/_readdy_docs/reference || true
    cd $HOME/_readdy_docs/reference
    # cant reliably determine cpu count in a docker container,
    # therefore fix this value.
    if [ "$TRAVIS" == "true" ]; then CPU_COUNT=2; fi
    cmake ${TRAVIS_BUILD_DIR} -DREADDY_GENERATE_DOCUMENTATION_TARGET_ONLY:BOOL=ON
    make doc
    cd -
}

function make_website_doc {
    cd $HOME
    # get the jekyll source code and setup jekyll via bundler
    git clone "https://github.com/readdy/readdy_documentation.git"
    cd readdy_documentation/readdy_documentation
    bundle install
    # insert the reference documentation
    cp -r $HOME/_readdy_docs/reference/docs/html/* reference_manual/
    # insert the content/pages and move the index.md to the jekyll root
    mkdir _pages || true
    cp -r ${TRAVIS_BUILD_DIR}/docs/documentation/* _pages/
    mv _pages/index.md ./index.md
    # build
    bundle exec jekyll build
    cd _site
    rm Gemfile Gemfile.lock
    cd -
}

function setup_docs_repo {
    cd $HOME/readdy_documentation/readdy_documentation/_site
    git init
    git config user.name "Christoph Froehner"
    git config user.email "chrisfroe@users.noreply.github.com"
    git remote add upstream "https://$GH_TOKEN@github.com/readdy/readdy_documentation.git"
    git fetch upstream
    git checkout --orphan workbranch
    git reset --hard
    cd -
}

function deploy {
    # revision tag
    cd $HOME/readdy_documentation/readdy_documentation/_site/
    touch .
    git add -A .
    git commit -m "github pages"
    git push -q -f upstream workbranch:gh-pages
}

set_this_up

make_reference_doc
make_website_doc
setup_docs_repo

deploy
