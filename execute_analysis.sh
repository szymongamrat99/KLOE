#!/bin/bash
nproc=$1
file_list=$2

cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
# cmake -DCMAKE_BUILD_TYPE=Release ..
make -j$nproc

../copy_libs.sh

./bin/KLSPM00 $file_list
