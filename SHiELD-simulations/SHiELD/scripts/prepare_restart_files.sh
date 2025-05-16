#!/bin/bash

source=${1}
destination=${2}
timestamps=("${@:3}")

segment_end=$(basename ${source})

for timestamp in "${timestamps[@]}"
do
    timestamp_destination=${destination}/${timestamp}
    mkdir -p ${timestamp_destination}
    if [ "${timestamp}" = "${segment_end}" ]; then
        cp ${source}/[a-z]* ${timestamp_destination}/
    else
        restart_timestamp=${timestamp::-2}.000000
        cp ${source}/${restart_timestamp}.* ${timestamp_destination}/
        rename -v "${restart_timestamp}." "" ${timestamp_destination}/*
    fi
done
