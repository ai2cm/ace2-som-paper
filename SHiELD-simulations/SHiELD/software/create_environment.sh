#!/bin/bash

ENVIRONMENT_FILE=$1
ENVIRONMENT_NAME=$2
PIP_REQUIREMENTS=$3
SHiELD_RUN_DIR=$4

source /sw/gaea-c5/python/3.9/anaconda-base/etc/profile.d/conda.sh
conda env create --file ${ENVIRONMENT_FILE} --name ${ENVIRONMENT_NAME}
conda activate ${ENVIRONMENT_NAME}
pip install -r ${PIP_REQUIREMENTS}
pip install -e ${SHiELD_RUN_DIR}
