#!/bin/bash
nproc=$1
file_list=$2

cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
# cmake -DCMAKE_BUILD_TYPE=Release ..
make -j$nproc

source ../../.env
../copy_libs.sh

echo "Using ANALYSIS_CONFIG_FILE: $ANALYSIS_CONFIG_FILE"
export ANALYSIS_CONFIG_FILE=${ANALYSIS_CONFIG_FILE}
./bin/KLSPM00 $file_list
