#!/bin/bash
# This allows us to use conda run without having conda active in our root
# environment.
CONDA_ENVIRONMENT=$1
COMMAND=$2

source /sw/gaea-c5/python/3.9/anaconda-base/etc/profile.d/conda.sh
conda run --no-capture-output --live-stream -n ${CONDA_ENVIRONMENT} ${COMMAND} "${@:3}"
