#!/usr/bin/env bash

for ima in $(ls | grep BC); do

    cd $ima/M1/z-positive/ && \
    rm trajectory-bc-pos.in && \
    cp ../../../trajectory-bc-pos.txt . && \
    mv trajectory-bc-pos.txt trajectory-bc-pos.in &&
    cd ../../../;

    cd $ima/M1/z-negative/ && \
    rm trajectory-bc-neg.in && \
    cp ../../../trajectory-bc-neg.txt . && \
    mv trajectory-bc-neg.txt trajectory-bc-neg.in &&
    cd ../../../;

done
