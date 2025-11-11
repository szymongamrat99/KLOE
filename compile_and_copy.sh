#!/bin/bash
nproc=$1

cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j$nproc
cp bin/KLSPM00 /opt/exp_software/kloe/users/gamrat/KLSPM00/bin/KLSPM00.exe
cp Include/*/*.so /opt/exp_software/kloe/users/gamrat/KLSPM00/lib/.
cp Subanalysis/*/*.so /opt/exp_software/kloe/users/gamrat/KLSPM00/lib/.
cd ..
cp -fr Subanalysis/Properties /opt/exp_software/kloe/users/gamrat/KLSPM00/.
