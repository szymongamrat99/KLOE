#!/bin/bash
nproc=$1
file_list=$2

cd build
cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo ..
make -j$nproc
./bin/KLSPM00 $file_list
