#!/usr/bin/env bash

conda_package_file=$(conda build tools/conda-recipe --python="${CONDA_PY}" --output | grep '.tar.bz2' | tail -1)
echo "Found conda package file ${conda_package_file}"

conda install anaconda-client -qy

if [ "${BUILD_REASON}" == "IndividualCI" ]; then
  echo "Build source branch: ${BUILD_SOURCE_BRANCH}"
  if [[ "${BUILD_SOURCE_BRANCH}" == "refs/tags/"* ]]; then
    echo "Uploading tagged package ${conda_package_file} with tag ${BUILD_SOURCE_BRANCH} to regular channel"
    anaconda -t "${BINSTAR_TOKEN}" upload -u readdy "${conda_package_file}" > /dev/null 2>&1
  else
    echo "Uploading package ${conda_package_file} to dev channel"
    anaconda -t "${BINSTAR_TOKEN}" upload -c readdy -u readdy -l dev --force "${conda_package_file}" > /dev/null 2>&1
  fi
else
  echo "Skipping upload because the build reason was ${BUILD_REASON} and not 'IndividualCI'"
fi

echo "Uploaded. Done."
