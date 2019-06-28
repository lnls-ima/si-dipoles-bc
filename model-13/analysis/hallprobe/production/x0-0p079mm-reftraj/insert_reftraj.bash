#!/usr/bin/env bash

for ima in $(ls | grep BC); do

    cd $ima/M1/z-positive/ && \
    sed -i '/traj_load_filename/s/None/"trajectory-bc-pos.in"/g' trajectory.in && \
    cd ../../../;

    cd $ima/M1/z-negative/ && \
    sed -i '/traj_load_filename/s/None/"trajectory-bc-neg.in"/g' trajectory.in && \
    cd ../../../;

done
