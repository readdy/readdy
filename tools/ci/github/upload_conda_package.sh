#!/usr/bin/env bash

# Uploading the package requires:
# - An already built conda package
# - A binstar token stored in envvar BINSTAR_TOKEN
# - A git ref stored in GITHUB_REF (provided by GitHub Actions)
#   pointing either to the head of a branch or to a tag

set -e -u

function resolve_github_ref() {
  tag="notag"
  branch="none"
  if [[ "${GITHUB_REF}" == "refs/heads/"* ]]; then
    branch=${GITHUB_REF:11}
    echo "GITHUB_REF ${GITHUB_REF} points to the head of branch ${branch}. There is no tag."
  elif [[ "${GITHUB_REF}" == "refs/tags/"* ]]; then
    tag=${GITHUB_REF:10}
    echo "GITHUB_REF ${GITHUB_REF} points to the tag ${tag}. There is no branchname."
  else
    echo "Could not resolve GITHUB_REF=${GITHUB_REF}. Something is wrong."
    exit 1
  fi
}

# catch the case when there is no tag and the branch is not master
function validate_this_should_run() {
    if [ -z ${PYTHON_VERSION+x} ]; then
    echo "PYTHON_VERSION was not set. Something is wrong."
    exit 1
  fi

  if [ -z ${BINSTAR_TOKEN+x} ]; then
    echo "BINSTAR_TOKEN was not set, so this is probably a fork. Exit."
    exit 0
  fi

  if [ -z ${GITHUB_REF+x} ]; then
    echo "GITHUB_REF was not set. Something is wrong."
    exit 1
  fi

  resolve_github_ref

  if [ "${tag}" == "notag" ] && [ "${branch}" != "master" ]; then
    echo "When there is no tag, the branch must be master (is ${branch}). Exit."
    exit 0
  else
    echo "O.K. There is no tag, but the branch is master."
  fi
}

# fixme removethis: this is up here for debugging the build pipeline
# resolve_github_ref
# echo "tag is ${tag}"
# echo "branch is ${branch}"

# validate_this_should_run

conda_package_file=$(conda build tools/conda-recipe --python="${PYTHON_VERSION}" --output | grep '.tar.bz2' | tail -1)
echo "Found conda package file ${conda_package_file}"

conda install anaconda-client -qy

if [ "$(Build.Reason)" == "IndividualCI" ]; then
  output=$(git describe --tags --exact-match "$(Build.SourceVersion)")
  echo "Git describe output: ${output}"
  status=$?
  if [ ${status} == 128 ]; then
    echo "Git describe status 128, not a tag."
    echo "Uploading package ${conda_package_file} to dev channel"
    anaconda -t "${BINSTAR_TOKEN}" upload -c readdy -u readdy -l dev --force "${conda_package_file}" > /dev/null 2>&1
  else
    echo "Git describe successful with tag ${output}."
    echo "Uploading tagged package ${conda_package_file} with tag ${tag} to regular channel"
    anaconda -t "${BINSTAR_TOKEN}" upload -u readdy "${conda_package_file}" > /dev/null 2>&1
  fi
else
  echo "Skipping upload because the build reason was $(Build.Reason) and not 'IndividualCI'"
fi

echo "Uploaded. Done."
