#!/usr/bin/env bash

# Uploading the package requires:
# - An already built conda package
# - A binstar token stored in envvar BINSTAR_TOKEN
# - A git ref stored in GITHUB_REF (provided by GitHub Actions)
#   pointing either to the head of a branch or to a tag

conda_package_file=$(conda build tools/conda-recipe --python="${PYTHON_VERSION}" --output | grep '.tar.bz2' | tail -1)
echo "Found conda package file ${conda_package_file}"

conda install anaconda-client -qy

if [ "${BUILD_REASON}" == "IndividualCI" ]; then
  output=$(git describe --tags --exact-match "${BUILD_SOURCE_VERSION}")
  status=$?
  echo "Git describe output: ${output}"
  if [ ${status} == 128 ] || [ -z "${output}" ]; then
    echo "Git describe status 128 (was ${status}) or output empty, not a tag."
    echo "Uploading package ${conda_package_file} to dev channel"
    anaconda -t "${BINSTAR_TOKEN}" upload -c readdy -u readdy -l dev --force "${conda_package_file}" > /dev/null 2>&1
  else
    echo "Git describe successful with tag ${output}."
    echo "Uploading tagged package ${conda_package_file} with tag ${output} to regular channel"
    anaconda -t "${BINSTAR_TOKEN}" upload -u readdy "${conda_package_file}" > /dev/null 2>&1
  fi
else
  echo "Skipping upload because the build reason was ${BUILD_REASON} and not 'IndividualCI'"
fi

echo "Uploaded. Done."
