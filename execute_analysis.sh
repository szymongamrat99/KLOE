#!/bin/bash
cd build
cmake -DENABLE_PROFILING=ON ..
make -j 5
./bin/KLSPM00
