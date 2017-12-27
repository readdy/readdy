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

# workaround: to prevent other output than the actual file, call this once before,
# so that packages are already downloaded and extracted
conda build tools/conda-recipe --output > /dev/null 2>&1

CONDA_PACKAGE_FILE=$(conda build tools/conda-recipe --output)
echo "found conda package file $CONDA_PACKAGE_FILE"

conda install anaconda-client -qy

tagval=${TRAVIS_TAG:-notag}

if [ "$tagval" == "notag" ]
then
    echo "uploading dev package"
    anaconda -t $BINSTAR_TOKEN upload -c readdy -u readdy -l dev --force $CONDA_PACKAGE_FILE
else
    echo "uploading tagged package with tag $tagval"
    anaconda -t $BINSTAR_TOKEN upload -u readdy $CONDA_PACKAGE_FILE
fi
