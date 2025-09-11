#!/bin/bash
cd build
make -j 5
taskset -c 4-7 bin/KLSPM00