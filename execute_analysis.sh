#!/bin/bash
cd build
cmake ..
make -j 5
./bin/KLSPM00
