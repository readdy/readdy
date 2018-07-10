#!/usr/bin/env bash

set -e -u

function set_this_up {
    ttagval=${TRAVIS_TAG:-notag}
    if [ "$TRAVIS_PULL_REQUEST" != "false" ]
    then
        echo "This was a pull request, thus dont build docs. Exit."
        exit 0
    fi
    if [ "$TRAVIS_BRANCH" != "master" ] && [ "$ttagval" == "notag" ]
    then
        echo "This commit was made against the $TRAVIS_BRANCH branch and not the master branch. Exit."
        exit 0
    fi
    if [ -z ${BINSTAR_TOKEN+x} ]
    then
      echo "BINSTAR_TOKEN was not set, so this is probably a fork. Exit."
      exit 0
    fi
}

set_this_up

CONDA_PACKAGE_FILE=$(conda build tools/conda-recipe --output | grep '.tar.bz2' | tail -1)
echo "found conda package file $CONDA_PACKAGE_FILE"

conda install anaconda-client -qy

tagval=${TRAVIS_TAG:-notag}

if [ "$tagval" == "notag" ]
then
    echo "uploading devel package"
    anaconda -t $BINSTAR_TOKEN upload -c readdy -u readdy -l dev --force $CONDA_PACKAGE_FILE
else
    echo "uploading tagged package with tag $tagval"
    anaconda -t $BINSTAR_TOKEN upload -u readdy $CONDA_PACKAGE_FILE
fi
