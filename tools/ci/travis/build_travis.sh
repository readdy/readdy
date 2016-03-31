#!/usr/bin/env bash

conda build -q tools/conda-recipe | grep -v "\-\- Installing:" | grep -v "\-\- Up\-to\-date:"