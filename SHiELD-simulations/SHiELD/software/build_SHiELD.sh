#!/bin/bash

ROOT=$(dirname "$0")
SHiELD_BUILD=${ROOT}/submodules/SHiELD_build

cd ${SHiELD_BUILD}/Build
./COMPILE clean shield 64bit repro intel
