#!/bin/bash
nproc=$1

cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j$(nproc)
./bin/KLSPM00
