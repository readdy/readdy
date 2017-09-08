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

    gem install travis
}

function trigger_doc_job {
    # obtain a travis token in exchange for our github token
    yes | travis login --org --github-token $GH_TOKEN > /dev/null 2>&1
    TRAVIS_TOKEN=$(travis token --org)

    # send the trigger request via curl using the Travis API
    body='{
    "request": {
    "branch":"master"
    }}'

    curl -s -X POST \
       -H "Content-Type: application/json" \
       -H "Accept: application/json" \
       -H "Travis-API-Version: 3" \
       -H "Authorization: token $TRAVIS_TOKEN" \
       -d "$body" \
       https://api.travis-ci.org/repo/readdy%2Freaddy_documentation/requests > /dev/null 2>&1
}

set_this_up
trigger_doc_job
