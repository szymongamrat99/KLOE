#!/bin/bash
nproc=$1
file_list=$2

HOST=$(hostname)

if [ "$HOST" == "kitt" ]; 
then
  set -a
  source /data/4/users/gamrat/.env
  set +a
fi

cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j$nproc

echo "Using ANALYSIS_CONFIG_FILE: $ANALYSIS_CONFIG_FILE"
export ANALYSIS_CONFIG_FILE=${ANALYSIS_CONFIG_FILE}
./bin/KLSPM00 $file_list

if [ "$HOST" == "ui-tier1.cr.cnaf.infn.it" ]; then
  ../copy_libs.sh
fi
