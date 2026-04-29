git pull
./execute_analysis.sh 
vi ~/.bashrc 
./execute_analysis.sh 
source ~/.bashrc
./execute_analysis.sh 
vi ~/.bashrc 
source ~/.bashrc
./execute_analysis.sh 
git add .
git commit -m "Changed env variables"
git push
vi ~/.bashrc 
source ~/.bashrc
git pull
vi ~/.bashrc 
git add .
git commit -m "Fixed"
git push
git pull
./execute_analysis.sh 
cd ROOT_Macros/plots/signal/signal_vs_bcg/filtering/
python calculateCombPi0MassCut.py 
cd ../../../../
cd ../
git pull
git add .
git commit -m "Additional cuts"
git push
./execute_analysis.sh 
git add .
git commit -m "Fixed"
git push
./execute_analysis.sh 5 ../job_v26_all_phys3_1_inv_pb_1.txt 
git add .
git commit -m "Fixed uv"
git pus
git push
./execute_analysis.sh 
git add .
git commit -m "Fixed"
git push
exit
./execute_analysis.sh 
2
git add .
git commit -m "Fixed"
git push
ls
cd ..
ls
cd CNAF_Produced_Files/
ls
cd root_files/
ls
scp tier1-cnaf:/storage/gpfs_small_files/kloe/gamrat/Subanalysis/InitialAnalysis/root_files/ALL_PHYS2* .
ls
ssh-copy-id bastion-cnaf 
ssh bastion-cnaf 
ssh-keygen -f "/home/gamrat/.ssh/known_hosts" -R "bastion.cnaf.infn.it"
ssh bastion-cnaf 
scp tier1-cnaf:/storage/gpfs_small_files/kloe/gamrat/Subanalysis/InitialAnalysis/root_files/ALL_PHYS2* .
scp -r tier1-cnaf:/storage/gpfs_small_files/kloe/gamrat/Subanalysis/InitialAnalysis/root_files/ALL_PHYS2* .
root 
cd ALL_PHYS2_SIGNAL_NoSmearing/
root
cd 
exit
cd ROOT_Macros/plots/signal/signal_vs_bcg/filtering/
python rootFilesFilter_new.py 
cd ROOT_Macros/plots/signal/signal_vs_bcg/filtering/
python rootFilesFilter_new.py 
python calculateCombPi0MassCut.py 
./execute_analysis.sh 
cd ROOT_Macros/plots/signal/signal_vs_bcg/
./plots_MC_only.sh 
./plots_MC_only.sh new new
vi ~/.rootlogon.C 
cat ~/.rootlogon.C 
./plots_MC_only.sh new new
cd ROOT_Macros/plots/signal/signal_vs_bcg/
./plots_MC_only.sh 
./plots_MC_only.sh new new
cd ../CNAF_Produced_Files/root_files/
ls
scp tier1-cnaf:/storage/gpfs_small_files/kloe/gamrat/Subanalysis/InitialAnalysis/root_files/ALL_PHYS_* .
scp -r tier1-cnaf:/storage/gpfs_small_files/kloe/gamrat/Subanalysis/InitialAnalysis/root_files/ALL_PHYS_* .
scp -r tier1-cnaf:/storage/gpfs_small_files/kloe/gamrat/Subanalysis/InitialAnalysis/root_files/ALL_PHYS3_* .
ll
ls
cd ROOT_Macros/plots/signal/signal_vs_bcg/
./plots_MC_only.sh 
./plots_MC_only.sh new new
root --version
python - <<'PY'
import ROOT, glob
ROOT.EnableImplicitMT()
files = sorted(glob.glob('/data/ssd/gamrat/CNAF_Produced_Files/root_files/ALL_PHYS2_SIGNAL_NoSmearing/mk0*.root'))
rdf = ROOT.RDataFrame('h1', files)
rdf2 = rdf.Define('mctruth_int','mctruth == 0 ? 1 : mctruth')
print('count=', int(rdf2.Count().GetValue()))
print('TEST2_DONE')
PY

cd ../root_files/
root
lss
ls
cd 2025-11-25/
ls
cd ../../
cd CNAF_Produced_Files/
ls
cd root_files/
ls
root
cd ../CNAF_Produced_Files/root_files/
ls
rm -fr ALL_PHYS3*
scp -r tier1-cnaf:/storage/gpfs_small_files/kloe/gamrat/Subanalysis/InitialAnalysis/root_files/ALL_PHYS3_* .
ssh-copy-id bastion-cnaf 
ssh tier1-cnaf 
cd
du -h
df -h
mkdir .cern-root
cd .c
cd .cern-root/
mkdir root_src
mkdir root_build
mkdir root_install
cd root_src/
git clone --branch latest-stable --depth=1 --tags v6-36-10 https://github.com/root-project/root.git .
git clone --branch latest-stable  --tags v6-36-10 --depth=1 https://github.com/root-project/root.git .
git clone --branch latest-stable --depth=1 https://github.com/root-project/root.git .
ls
cd ../root_build/
pwd
cmake -DCMAKE_INSTALL_PREFIX=/home/gamrat/.cern-root/root_install -Dbuiltin_nlohmannjson=ON -Dbuiltin_veccore=ON -Ddev=ON -Dfortran=ON -Dsoversion=ON -Dxrootd=OFF -Dasan=ON -Droottest=ON -Dtesting=ON -Dfail-on-missing=ON -Dclingtest=ON  /home/gamrat/.cern-root/root_src
cd ../CNAF_Produced_Files/root_files/
ls
scp -r tier1-cnaf:/storage/gpfs_small_files/kloe/gamrat/Subanalysis/InitialAnalysis/root_files/DATA_SIGNAL_NoSmearing* .
cd DATA_SIGNAL_NoSmearing/
ls
ll
ls -hl
exit
git checkout -- .
git add .
./execute_analysis.sh 
nohup ./execute_analysis.sh 5 ../job_v26_all_phys3_1_inv_pb_1.txt < parameters.txt > nohup.log &
git add .
git commit -m "Changed floats to doubles"
git push
./execute_analysis.sh 
nohup ./execute_analysis.sh 5 ../job_v26_all_phys3_1_inv_pb_1.txt < parameters.txt > nohup.log &
git add .
git commit -m "Added SVD to kin fit"
git push
git add .
git commit -m "New histos"
git push
cd Subanalysis/InterfFunction/
./compile.sh 
vi ~/.bashrc 
./compile.sh 
vi ~/.bashrc 
source ~/.bashrc 
./compile.sh 
./hist_fits_integral_range 
vi ~/.bashrc 
source ~/.bashrc 
./hist_fits_integral_range 
vi ~/.bashrc 
source ~/.bashrc 
./hist_fits_integral_range 
source ~/.bashrc 
./compile.sh 
./hist_fits_integral_range 
./compile.sh 
./hist_fits_integral_range 
./compile.sh 
./hist_fits_integral_range 
./compile.sh 
./hist_fits_integral_range 
./compile.sh 
./hist_fits_integral_range 
./compile.sh 
git checkout -- .
./execute_analysis.sh 
./execute_analysis.sh 5 ../job_v26_all_phys3_1_inv_pb_1.txt 
./execute_analysis.sh 
cd Subanalysis/InitialAnalysis/root_files/ALL_PHYS3_SIGNAL_NoSmearing/
ls
ll
ls -hl
root mk0_initial_analysis_all_phys3_SIGNAL_NoSmearing_Signal_1.root
cd ../../../../
./execute_analysis.sh 
git commit -m "Added new things"
git push
cd Subanalysis/InterfFunction/
ls
./compile.sh 
./hist_fits_integral_range 
nohup ./hist_fits_integral_range < parameters.txt > nohup_1.log &
nohup ./hist_fits_integral_range < parameters.txt > nohup_2.log &
nohup ./hist_fits_integral_range < parameters.txt > nohup_3.log &
nohup ./hist_fits_integral_range < parameters.txt > nohup_4.log &
nohup ./hist_fits_integral_range < parameters.txt > nohup_5.log &
nohup ./hist_fits_integral_range < parameters.txt > nohup_6.log &
nohup ./hist_fits_integral_range < parameters.txt > nohup_7.log &
nohup ./hist_fits_integral_range < parameters.txt > nohup_8.log &
nohup ./hist_fits_integral_range < parameters.txt > nohup_9.log &
nohup ./hist_fits_integral_range < parameters.txt > nohup_10.log &
killall hist_fits_integral_range 
pgrep hist_fits_integral_range
nohup ./hist_fits_integral_range < parameters.txt > nohup_1.log &
nohup ./hist_fits_integral_range < parameters.txt > nohup_2.log &
nohup ./hist_fits_integral_range < parameters.txt > nohup_3.log &
nohup ./hist_fits_integral_range < parameters.txt > nohup_4.log &
nohup ./hist_fits_integral_range < parameters.txt > nohup_5.log &
nohup ./hist_fits_integral_range < parameters.txt > nohup_6.log &
nohup ./hist_fits_integral_range < parameters.txt > nohup_7.log &
nohup ./hist_fits_integral_range < parameters.txt > nohup_8.log &
nohup ./hist_fits_integral_range < parameters.txt > nohup_9.log &
nohup ./hist_fits_integral_range < parameters.txt > nohup_10.log &
nohup ./hist_fits_integral_range < parameters.txt > nohup_11.log &
nohup ./hist_fits_integral_range < parameters.txt > nohup_12.log &
nohup ./hist_fits_integral_range < parameters.txt > nohup_13.log &
nohup ./hist_fits_integral_range < parameters.txt > nohup_14.log &
nohup ./hist_fits_integral_range < parameters.txt > nohup_15.log &
nohup ./hist_fits_integral_range < parameters.txt > nohup_16.log &
nohup ./hist_fits_integral_range < parameters.txt > nohup_17.log &
nohup ./hist_fits_integral_range < parameters.txt > nohup_18.log &
nohup ./hist_fits_integral_range < parameters.txt > nohup_19.log &
nohup ./hist_fits_integral_range < parameters.txt > nohup_20.log &
nohup ./hist_fits_integral_range < parameters.txt > nohup_21.log &
nohup ./hist_fits_integral_range < parameters.txt > nohup_22.log &
nohup ./hist_fits_integral_range < parameters.txt > nohup_23.log &
nohup ./hist_fits_integral_range < parameters.txt > nohup_24.log &
nohup ./hist_fits_integral_range < parameters.txt > nohup_25.log &
nohup ./hist_fits_integral_range < parameters.txt > nohup_26.log &
nohup ./hist_fits_integral_range < parameters.txt > nohup_27.log &
nohup ./hist_fits_integral_range < parameters.txt > nohup_28.log &
nohup ./hist_fits_integral_range < parameters.txt > nohup_29.log &
lscpu
./execute_analysis.sh 
git add .
git commit -m "Simplified derivative"
git push
./execute_analysis.sh 
git add .
git commit -m "Simplified derivative"
git push
cd Subanalysis/InterfFunction/
./compile.sh 
cd img/integral_range/
rm -fr *
cd ../../
ls
chmod +x run_fit_int_range.sh 
./run_fit_int_range.sh 
./run_fit_int_range.sh 10
echo 1 | ./run_fit_int_range.sh 10
nohup echo 1 | ./run_fit_int_range.sh 10 > nohup_1.log &
killall hist_fits_integral_range 
nohup echo 1 | ./run_fit_int_range.sh 10 > nohup_1.log 2>&1 &
nohup echo 1 | ./run_fit_int_range.sh 20 > nohup_1.log 2>&1 &
nohup echo 1 | ./run_fit_int_range.sh 40 > nohup_1.log 2>&1 &
nohup echo 1 | ./run_fit_int_range.sh 60 > nohup_1.log 2>&1 &
nohup echo 1 | ./run_fit_int_range.sh 80 > nohup_1.log 2>&1 &
nohup echo 1 | ./run_fit_int_range.sh 100 > nohup_1.log 2>&1 &
nohup echo 1 | ./run_fit_int_range.sh 150 > nohup_1.log 2>&1 &
nohup echo 1 | ./run_fit_int_range.sh 200 > nohup_1.log 2>&1 &
nohup echo 1 | ./run_fit_int_range.sh 250 > nohup_1.log 2>&1 &
nohup echo 1 | ./run_fit_int_range.sh 300 > nohup_1.log 2>&1 &
nohup echo 2 | ./run_fit_int_range.sh 300 > nohup_1.log 2>&1 &
nohup echo 2 | ./run_fit_int_range.sh 250 > nohup_1.log 2>&1 &
nohup echo 2 | ./run_fit_int_range.sh 200 > nohup_1.log 2>&1 &
nohup echo 2 | ./run_fit_int_range.sh 150 > nohup_1.log 2>&1 &
nohup echo 2 | ./run_fit_int_range.sh 100 > nohup_1.log 2>&1 &
nohup echo 2 | ./run_fit_int_range.sh 80 > nohup_1.log 2>&1 &
nohup echo 2 | ./run_fit_int_range.sh 60 > nohup_1.log 2>&1 &
nohup echo 2 | ./run_fit_int_range.sh 40 > nohup_1.log 2>&1 &
nohup echo 2 | ./run_fit_int_range.sh 20 > nohup_1.log 2>&1 &
nohup echo 2 | ./run_fit_int_range.sh 10 > nohup_1.log 2>&1 &
killall hist_fits_integral_range 
./compile.sh 
vi run_fit_int_range.sh 
./run_fit_int_range.sh 
exit
cd Subanalysis/InterfFunction
ls
cd img/theoretical_plots/
ls
cd ..
cd root_files/
ls
cd exp_corrected/
ll
ls -hl
cd ..
ls
cd ..
ls
cd img/
ls -hl
cd theoretical_plots/
ls -hl
quota
ls
cd 
rm -fr .cern-root/
cd $WORKDIR
ls
cd ..
exit
ls
cd ..
mkdir .cern-root
cd .cern-root/
mkdir root_src root_install root_build
git clone --branch latest-stable --depth=1 https://github.com/root-project/root.git root_src
ls
cd root_build/
cmake -DCMAKE_INSTALL_PREFIX=../root_install -Dbuiltin_nlohmannjson=ON -Dbuiltin_veccore=ON -Ddev=ON -Dfortran=ON -Dsoversion=ON -Dxrootd=OFF -Dasan=ON -Droottest=ON -Dtesting=ON -Dfail-on-missing=ON -Dclingtest=ON  ../root_src
apt list -a libxxhash-dev
cmake --version
cd ../.cern-root/
ls
cd root_build/
cmake -DCMAKE_INSTALL_PREFIX=../root_install -Dbuiltin_nlohmannjson=ON -Dbuiltin_veccore=ON -Ddev=ON -Dfortran=ON -Dsoversion=ON -Dxrootd=OFF -Dasan=ON -Droottest=ON -Dtesting=ON -Dfail-on-missing=ON -Dclingtest=ON  ../root_src
cd ../root_src/
cat README
ls
cat requirements.txt 
ls
cd README/
ls
cat INSTALL 
cd ../../
cd root_build/
cmake -DCMAKE_INSTALL_PREFIX=../root_install -Dbuiltin_nlohmannjson=ON -Dbuiltin_veccore=ON -Ddev=ON -Dfortran=ON -Dsoversion=ON -Dxrootd=OFF -Dasan=ON -Droottest=ON -Dtesting=ON -Dfail-on-missing=ON -Dclingtest=ON  ../root_src
cd ../root_src/
ls
vi CMakeLists.txt 
llq
cd 
ls
cd
cd /data/ssd/gamrat
ls
cd Node1/
ls
cd ..
ls -a
cd .cern-root/
ls
cd root_build/
ls
cd ..
cd root_install/
ls
cd ..
exit
root snap
root
vi ~/.bashrc
exit
./execute_analysis.sh 
root
cd build
cmake -DCMAKE_PREFIX_PATH=/snap/root/current/usr/local ..
make
cd ..
./execute_analysis.sh 
cd build/
rm -fr *
cd ..
./execute_analysis.sh 
vi ~/.bashrc 
./execute_analysis.sh 
./execute_analysis.sh 
vi ~/.bashrc
echo $CMAKE_PREFIX_PATH 
vi ~/.bashrc
echo $CMAKE_PREFIX_PATH 
./execute_analysis.sh 
cd build/
rm -fr *
cd ..
./execute_analysis.sh 
cmake --version
./execute_analysis.sh 
cd build/
rm -fr *
cd ..
./execute_analysis.sh 
vi ~/.bashrc 
./execute_analysis.sh 
cd build/
rm -fr *
cd ..
./execute_analysis.sh 
vi ~/.bashrc 
source ~/.bashrc
./execute_analysis.sh 
cd /usr/
ls
cd lib
ls
cd ..
ls
cd include/
ls
cd boost/
ls
cd ..
ls
cd ..
ls
cd lib64
ls
cd ..
cd lib32
ls
cd ..
cd lib
cd ..
ls
cd local/
ls
cd lib/
ls
cd ..
cd include/
ls
cd ..
cd src/
ls
cd ..
cd bin/
ls
cd ../../
cd bin/
ls
ls boost
cd ..
ls
find --name boost
find -name boost
cd include/boost/
ls
find -name BoostConfig
find -name BoostConfig.cmake
find -name FindBoost.cmake
cd 
cd /data/ssd/gamrat/KLOE/
./execute_analysis.sh 
source ~/.bashrc
./execute_analysis.sh 
cd build/
rm -fr *
cd ..
./execute_analysis.sh 
cd build/
make VERBOSE=1
vi ~/.bashrc 
root
cd /snap/root-framework/current/lib
ls
cd ../
cd lib64
ls
cd ../../
ls
cd current
ls
cd include/
ls
cd ../lib
ls
cd python3.12/
ls
cd site-packages/
ls
cd ../../../
cd ..
ls
cd current
cd ..
cd 954/
ls
cd lib
ls
cd ..
cd share/
ls
cd ..
cd usr/
ls
cd include/
ls
cd ../
cd lib
ls
cd x86_64-linux-gnu/
ls
cd ..
ls
cd ..
ls
cd include/
ls
cd ..
ls
cd ../
ls
cd ../current
ls
cd usr/
ls
cd local/
ls
cd i
cd include/
ls
cd ../
cd include/
ls
vi ~/.bashrc
cd build/
rm -fr *
cd ..
./execute_analysis.sh 
g++ --version
vi ~/.bashrc 
source ~/.bashrc
./execute_analysis.sh 
cd build/
rm -fr *
cd ..
./execute_analysis.sh 
cd build/
rm -fr *
cd ../
./execute_analysis.sh 
ls /usr/include/boost/version.hpp
./execute_analysis.sh 
cd build/
rm -fr *
cd ..
./execute_analysis.sh 
root
vi ~/.bashrc 
source ~/.bashrc 
./execute_analysis.sh 
cd build/
rm -fr *
cd ..
./execute_analysis.sh 
cd build/
make -DCMAKE_PREFIX_PATH=/usr
make -DCMAKE_PREFIX_PATH=/usr .
make -DCMAKE_PREFIX_PATH=/usr .ls
ls
cd ..
ls
vi ~/.bashrc 
./execute_analysis.sh 
cd build/
rm -fr *
./execute_analysis.sh 
cd ..
./execute_analysis.sh 
cmake --version
cd build/bin/
cd ..
rm -fr *
cd ..
./execute_analysis.sh 
cd Subanalysis/InterfFunction/
./compile.sh 
./bttf_sample 
ulimit -n
./execute_analysis.sh 
cd Subanalysis/InterfFunction/
ls
rm nohup_*
ls
rm log_*
ls
./compile.sh 
./bttf_sample 
./bttf_sample 1
./bttf_sample 10
./compile.sh 
./bttf_sample 10
./bttf_sample 1
./compile.sh 
./bttf_sample 10
./compile.sh 
./bttf_sample 10
./compile.sh 
./bttf_sample 11
./bttf_sample 12
./compile.sh 
./bttf_sample 12
./bttf_sample 13
./compile.sh 
./bttf_sample 13
./bttf_sample
./compile.sh 
./bttf_sample
./bttf_sample 13
./compile.sh 
./bttf_sample 13
./bttf_sample 10
./compile.sh 
./bttf_sample 10
./compile.sh 
./bttf_sample 10
./compile.sh 
./bttf_sample 10
./bttf_sample 11
./bttf_sample 12
git add .
git commit -m "Fixed Interf for jobs"
git push
cd Subanalysis/InterfFunction/root_files/sampling/DoubleExponential
ls
rm -fr *
ls
cd t1Min_0_t1Max_20_t2Min_0_t2Max_20_nSamples_1000000_customParamsFlag_false_ReParam_0.00166_ImParam_-0.00198/
ls
q
/q
exit
llq
ll
g++ version
g++ --version
c++ --version
ls 
cd Subanalysis/InterfFunction/
ls
parallel -j 5 ./bttf_sample {#} ::: {1..5}
cd root_files/sampling/DoubleExponential/t1Min_0_t1Max_20_t2Min_0_t2Max_20_nSamples_1000000_customParamsFlag_false_ReParam_0.00166_ImParam_-0.00198/
ls
root
root sampling_results_job_1_1.root 
rm -fr *
ls
cd ../../../../
cd img/theoretical_plots/integralLimit_300_samples_10000000000_Re_0.00166_Im_-0.00198_tmp/
ls
ls -hl

cd ../../
cd ../
nohup parallel --results -j 5 ./bttf_sample {#} ::: {1..6} > parallel_output.log 2>&1 & 
ls
cd -j/
cd -j
rm -j
rm -fr "-j"
rm -fr ./-j
ls
nohup parallel --results log_parallel -j 5 ./bttf_sample {#} ::: {1..6} > parallel_output.log 2>&1 & 
ls
cd log_parallel/
ls
cd 1
cd ../
ls
cd ..
nohup parallel -j 5 ./bttf_sample {#} ::: {1..6} > parallel_output.log 2>&1 & 
ls
cat parallel_output.log 
nohup parallel -j 6 ./bttf_sample {#} ::: {1..6} > parallel_output.log 2>&1 & 
cat parallel_output.log 
ps -e | grep bttf
cat parallel_output.log 
killall bttf_sample
nohup parallel -j 6 ./bttf_sample {#} ::: {1..6} > parallel_output.log 2>&1 & 
cd root_files/sampling/DoubleExponential/t1Min_0_t1Max_20_t2Min_0_t2Max_20_nSamples_10000000000_customParamsFlag_false_ReParam_0.00166_ImParam_-0.00198/
ls -hl
tmux
ls -hl
top
htop
ls -hl
tmux ls
ls -hl
top
ls -hl
cd ..
ls
cd ..
ls
cd ..
ls
rm nohup*
ls
rm -r log_parallel/
ls
cat parallel_output.log 
cd 
cd /data/ssd/gamrat/PROD2ROOT/
ls
cd MC
ls
cd MK0/
ls
cd all_phys
ls
sha256sum *.root
sha256sum *.root > checksums.root
ls
vi checksums.root 
scp checksums.root tier1-cnaf:~/.
ssh tier1-cnaf 
ls
rm checksums.root 
cd /data/tape_dump/gamrat/PROD2ROOT/MC/MK0/all_phys
scp tier1-cnaf:/home/g/gamrat/sum.txt .
sha256sum -c sum.txt 
ls
ssh tier1-cnaf 
./execute_analysis.sh 
cd Subanalysis/InterfFunction/
./compile.sh 
ls
./hist_fits_integral_range
./compile.sh
./hist_fits_integral_range
./compile.sh
./hist_fits_integral_range
cd Subanalysis/InterfFunction/
./compile.sh 
cd root_files/sampling/DoubleExponential/t1Min_0_t1Max_20_t2Min_0_t2Max_20_nSamples_10000000000_customParamsFlag_false_ReParam_0.00166_ImParam_-0.00198/
ls -hl
cd Subanalysis/InterfFunction/root_files/sampling/DoubleExponential/t1Min_0_t1Max_20_t2Min_0_t2Max_20_nSamples_10000000000_customParamsFlag_false_ReParam_0.00166_ImParam_-0.00198/
ls -hl
./execute_analysis.sh 
./execute_analysis.sh 5 ../job_v26_all_phys3_1_inv_pb_1.txt 
ls /data/ssd/gamrat/PROD2ROOT/MC/MK0/all_phys3/prod2root_mk0_all_phys3_30342_v2.root
ls /data/tape_dump/gamrat/PROD2ROOT/MC/MK0/all_phys3/prod2root_mk0_all_phys3_30342_v2.root
./execute_analysis.sh 5 ../job_v26_all_phys3_1_inv_pb_1.txt 
git add .
git commit -m "Done"
git push
./execute_analysis.sh 
git add .
git commit -m "Done"
git push
cd Subanalysis/InterfFunction/root_files/sampling/DoubleExponential/t1Min_0_t1Max_20_t2Min_0_t2Max_20_nSamples_10000000000_customParamsFlag_false_ReParam_0.00166_ImParam_-0.00198/
ls -hl
du -h
df -h
pgrep bttf
ps -e
exit
./execute_analysis.sh 
git add .
git commit -m "Debug"
git push
./execute_analysis.sh 5 ../job_v26_all_phys3_1_inv_pb_1.txt 
git add .
git commit -m "Debug"
git push
git add .
git commit -m "Debug"
git push
./execute_analysis.sh 5 ../job_v26_all_phys3_1_inv_pb_1.txt 
git add .
git commit -m "Debug"
git push
./execute_analysis.sh 5 ../job_v26_all_phys3_1_inv_pb_1.txt 
git add .
git commit -m "Debug"
git push
./execute_analysis.sh 
git add .
git commit -m "Debug"
git push
git add .
git commit -m "Debug"
git push
./execute_analysis.sh 
git add .
git commit -m "Debug"
git push
./execute_analysis.sh 
git add .
git commit -m "Debug"
git push
git add .
git commit -m "Debug"
git push
./execute_analysis.sh 
git add .
git commit -m "Debug"
git push
git add .
git commit -m "Debug"
git push
./execute_analysis.sh 
git add .
git commit -m "Debug"
git push
git add .
git commit -m "Debug"
git add .
git commit -m "Debug"
git push
./execute_analysis.sh 
git add .
git commit -m "Debug"
git push
./execute_analysis.sh 
git add .
git commit -m "Debug"
git push
./execute_analysis.sh 
git add .
git commit -m "Debug"
git push
./execute_analysis.sh 
git add .
git commit -m "Debug"
git push
git add .
git commit -m "Debug"
git add .
git commit -m "Debug"
./execute_analysis.sh 
git add .
git commit -m "Debug"
git push
./execute_analysis.sh 
git add .
git commit -m "Debug"
git push
git add .
git commit -m "Debug"
git push
git add .
git commit -m "Debug"
git push
git add .
git commit -m "Debug"
git push
./execute_analysis.sh 
./execute_analysis.sh 5 ../job_v26_all_phys3_1_inv_pb_1.txt 
./execute_analysis.sh 
./execute_analysis.sh 5 ../job_v26_all_phys3_1_inv_pb_1.txt 
./execute_analysis.sh 
./execute_analysis.sh 5 ../job_v26_all_phys3_1_inv_pb_1.txt 
./execute_analysis.sh 5 ../job_v26_all_phys3_1_inv_pb_1.txt
ls
exit
cd Subanalysis/InterfFunction/
./compile.sh 
./hist_sample_draw 
nohup ./hist_sample_draw < params.txt > nohup.log &
cd Subanalysis/InterfFunction/
ls
ps -e | grep hist_sample
ps -e | grep hist_sample_draw
ps -e
cd Subanalysis/InterfFunction/
./hist_fits_integral_range 
./compile.sh 
./hist_fits_integral_range 
./compile.sh 
./hist_fits_integral_range 
./compile.sh 
./hist_fits_integral_range 
./compile.sh 
./hist_fits_integral_range 
./compile.sh 
git add .
git commit -m "Scaling"
git push
cd Subanalysis/InitialAnalysis/
cd co
cd ..
cd InterfFunction/
cd config/
cat config.json 
cd ../img/theoretical_plots/integralLimit_300_samples_10000000000_Re_0.00166_Im_-0.00198_tmp/
ls
root
ls
cd Subanalysis/InterfFunction/root_files/sampling/DoubleExponential/t1Min_0_t1Max_20_t2Min_0_t2Max_20_nSamples_10000000000_customParamsFlag_false_ReParam_0.00166_ImParam_-0.00198/
ls
rm -fr *
cd ..
rm -fr *
ls
cd ../../
cd ..
nohup parallel -j 6 ./bttf_sample {#} ::: {1..6} > parallel_output.log 2>&1 & 
cd root_files/sampling/DoubleExponential/t1Min_0_t1Max_300_t2Min_0_t2Max_300_nSamples_10000000000_customParamsFlag_false_ReParam_0.00166_ImParam_-0.00198/
ls
ls -hl
cd Subanalysis/InterfFunction/root_files/sampling/DoubleExponential/t1Min_0_t1Max_300_t2Min_0_t2Max_300_nSamples_10000000000_customParamsFlag_false_ReParam_0.00166_ImParam_-0.00198/
ll
ls -hl
cd Subanalysis/InterfFunction/
./compile.sh 
./hist_fits_integral_range 
./compile.sh 
./hist_fits_integral_range 
./compile.sh 
./bttf_histo
./hist_sample_draw 
./compile.sh 
./hist_sample_draw /data/ssd/gamrat/KLOE/Subanalysis/InterfFunction/img/theoretical_plots/integralLimit_300_samples_10000000000_Re_0.00166_Im_-0.00198_tmp/sampling_results_300.root
./compile.sh 
./hist_sample_draw /data/ssd/gamrat/KLOE/Subanalysis/InterfFunction/img/theoretical_plots/integralLimit_300_samples_10000000000_Re_0.00166_Im_-0.00198_tmp/sampling_results_300.root
./compile.sh 
./hist_sample_draw /data/ssd/gamrat/KLOE/Subanalysis/InterfFunction/img/theoretical_plots/integralLimit_300_samples_10000000000_Re_0.00166_Im_-0.00198_tmp/sampling_results_300.root
./compile.sh 
./hist_sample_draw /data/ssd/gamrat/KLOE/Subanalysis/InterfFunction/img/theoretical_plots/integralLimit_300_samples_10000000000_Re_0.00166_Im_-0.00198_tmp/sampling_results_300.root
./compile.sh 
./hist_sample_draw /data/ssd/gamrat/KLOE/Subanalysis/InterfFunction/img/theoretical_plots/integralLimit_300_samples_10000000000_Re_0.00166_Im_-0.00198_tmp/sampling_results_300.root
./compile.sh 
./hist_sample_draw /data/ssd/gamrat/KLOE/Subanalysis/InterfFunction/img/theoretical_plots/integralLimit_300_samples_10000000000_Re_0.00166_Im_-0.00198_tmp/sampling_results_300.root
cd Subanalysis/InterfFunction/root_files/filled_sampled_histograms/exp_corrected/
ls
r, histograms2D_6452077916_0_exp_corrected_sigmat1.00_t10.00_300.00_t20.00_300.00_tch0.00_300.00_tne0.00_300.00_dcbe301d9fb7d7d8.root 
rm histograms2D_6452077916_0_exp_corrected_sigmat1.00_t10.00_300.00_t20.00_300.00_tch0.00_300.00_tne0.00_300.00_dcbe301d9fb7d7d8.root 
cd ..
ls
cd exp_corrected/
ls
rm histograms2D_6452077916_0_exp_corrected_sigmat1.00_t10.00_300.00_t20.00_300.00_tch0.00_300.00_tne0.00_300.00_dcbe301d9fb7d7d8.root 
cd ../
rm -fr *
ls
cd exp_corrected/
ls
cd re0.00166_im-0.00198_sigmat1.00_t10.00_300.00_t20.00_300.00_tch0.00_300.00_tne0.00_300.00/
ls
cd ..
rm -fr *
cd ../
ls
find /data/ssd/gamrat/KLOE/Subanalysis/InterfFunction/img/theoretical_plots/integralLimit_300_samples_10000000000_Re_0.00166_Im_-0.00198_tmp -name "*.root" | sort > inputs.txt
cat inputs.txt 
cd ..
nohup parallel -a root_files/filled_sampled_histograms/inputs.txt -j 6 ./hist_sample_draw {} > nohup_sigmat_1.0.log
nohup parallel -a root_files/filled_sampled_histograms/inputs.txt -j 6 ./hist_sample_draw {} > nohup_sigmat_1.0.log &
cat nohup_sigmat_1.0.log 
pgrep nohup
pgrep hist_sample_dra 
nohup parallel -a root_files/filled_sampled_histograms/inputs.txt -j 6 ./hist_sample_draw {} 1 0.5> nohup_sigmat_0.5.log &
nohup parallel -a root_files/filled_sampled_histograms/inputs.txt -j 6 ./hist_sample_draw {} 1 0.2> nohup_sigmat_0.2.log &
nohup parallel -a root_files/filled_sampled_histograms/inputs.txt -j 6 ./hist_sample_draw {} 1 0.8> nohup_sigmat_0.8.log &
nohup parallel -a root_files/filled_sampled_histograms/inputs.txt -j 6 ./hist_sample_draw {} 1 0.1> nohup_sigmat_0.1.log &
killall hist_sample_draw 
pgrep hist_sample_dra 
./compile.sh 
nohup parallel -a root_files/filled_sampled_histograms/inputs.txt -j 6 ./hist_sample_draw {} > nohup_sigmat_1.0.log &
nohup parallel -a root_files/filled_sampled_histograms/inputs.txt -j 6 ./hist_sample_draw {} 1 0.5> nohup_sigmat_0.5.log &
nohup parallel -a root_files/filled_sampled_histograms/inputs.txt -j 6 ./hist_sample_draw {} 1 0.2> nohup_sigmat_0.2.log &
nohup parallel -a root_files/filled_sampled_histograms/inputs.txt -j 6 ./hist_sample_draw {} 1 0.8> nohup_sigmat_0.8.log &
nohup parallel -a root_files/filled_sampled_histograms/inputs.txt -j 6 ./hist_sample_draw {} 1 0.1> nohup_sigmat_0.1.log &
killall hist_sample_draw 
nohup parallel --line-buffer -a root_files/filled_sampled_histograms/inputs.txt -j 6 ./hist_sample_draw {} > nohup_sigmat_1.0.log 2>&1 &
nohup parallel --line-buffer -a root_files/filled_sampled_histograms/inputs.txt -j 6 ./hist_sample_draw {} 1 0.1 > nohup_sigmat_0.1.log 2>&1 &
nohup parallel --line-buffer -a root_files/filled_sampled_histograms/inputs.txt -j 6 ./hist_sample_draw {} 1 0.2 > nohup_sigmat_0.2.log 2>&1 &
nohup parallel --line-buffer -a root_files/filled_sampled_histograms/inputs.txt -j 6 ./hist_sample_draw {} 1 0.5 > nohup_sigmat_0.5.log 2>&1 &
nohup parallel --line-buffer -a root_files/filled_sampled_histograms/inputs.txt -j 6 ./hist_sample_draw {} 1 0.8 > nohup_sigmat_0.8.log 2>&1 &
pgrep hist_sample_dra 
cd Subanalysis/InterfFunction/root_files/sampling/DoubleExponential/t1Min_0_t1Max_300_t2Min_0_t2Max_300_nSamples_10000000000_customParamsFlag_false_ReParam_0.00166_ImParam_-0.00198/
ll
ls -hl
cd Subanalysis/InterfFunction/root_files/sampling/DoubleExponential/t1Min_0_t1Max_300_t2Min_0_t2Max_300_nSamples_10000000000_customParamsFlag_false_ReParam_0.00166_ImParam_-0.00198/
ls -hl
pgrep bttf_sample
ls -hl
root sampling_results_job_6_tmp_1_5.root
cd ../../
cd ..
ls
cd filled_sampled_histograms/
find /data/ssd/gamrat/KLOE/Subanalysis/InterfFunction/root_files/sampling/DoubleExponential/t1Min_0_t1Max_300_t2Min_0_t2Max_300_nSamples_10000000000_customParamsFlag_false_ReParam_0.00166_ImParam_-0.00198 -name "*job_1*.root" | sort > inputs1.txt
find /data/ssd/gamrat/KLOE/Subanalysis/InterfFunction/root_files/sampling/DoubleExponential/t1Min_0_t1Max_300_t2Min_0_t2Max_300_nSamples_10000000000_customParamsFlag_false_ReParam_0.00166_ImParam_-0.00198 -name "*job_2*.root" | sort > inputs2.txt
find /data/ssd/gamrat/KLOE/Subanalysis/InterfFunction/root_files/sampling/DoubleExponential/t1Min_0_t1Max_300_t2Min_0_t2Max_300_nSamples_10000000000_customParamsFlag_false_ReParam_0.00166_ImParam_-0.00198 -name "*job_3*.root" | sort > inputs3.txt
find /data/ssd/gamrat/KLOE/Subanalysis/InterfFunction/root_files/sampling/DoubleExponential/t1Min_0_t1Max_300_t2Min_0_t2Max_300_nSamples_10000000000_customParamsFlag_false_ReParam_0.00166_ImParam_-0.00198 -name "*job_4*.root" | sort > inputs4.txt
find /data/ssd/gamrat/KLOE/Subanalysis/InterfFunction/root_files/sampling/DoubleExponential/t1Min_0_t1Max_300_t2Min_0_t2Max_300_nSamples_10000000000_customParamsFlag_false_ReParam_0.00166_ImParam_-0.00198 -name "*job_5*.root" | sort > inputs5.txt
find /data/ssd/gamrat/KLOE/Subanalysis/InterfFunction/root_files/sampling/DoubleExponential/t1Min_0_t1Max_300_t2Min_0_t2Max_300_nSamples_10000000000_customParamsFlag_false_ReParam_0.00166_ImParam_-0.00198 -name "*job_6*.root" | sort > inputs6.txt
ls
cat inputs6.txt 
cd ../../
nohup parallel -a root_files/filled_sampled_histograms/inputs1.txt -j 6 ./hist_sample_draw {} > nohup_sigmat_1.0_1.log &
nohup parallel -a root_files/filled_sampled_histograms/inputs2.txt -j 6 ./hist_sample_draw {} > nohup_sigmat_1.0_2.log &
nohup parallel -a root_files/filled_sampled_histograms/inputs3.txt -j 6 ./hist_sample_draw {} > nohup_sigmat_1.0_3.log &
nohup parallel -a root_files/filled_sampled_histograms/inputs4.txt -j 6 ./hist_sample_draw {} > nohup_sigmat_1.0_4.log &
nohup parallel -a root_files/filled_sampled_histograms/[Btograms/inputs5.txt -j 6 ./hist_sample_draw {} > nohup_sigmat_1.0_6.log &
nohup parallel -a root_files/filled_sampled_histograms/inputs6.txt -j 6 ./hist_sample_draw {} > nohup_sigmat_1.0_6.log &
cd root_files/filled_sampled_histograms/exp_corrected/re0.00166_im-0.00198_sigmat1.00_t10.00_300.00_t20.00_300.00_tch0.00_300.00_tne0.00_300.00/
ls
ls -hl
pgrep hist_sample_draw
pgrep hist_sample
pgrep hist_sample_
pgrep hist_sample_d
pgrep hist_sample_dr
pgrep hist_sample_dra
pgrep hist_sample_draw
ls
ls -hl
cd Subanalysis/InterfFunction/root_files/filled_sampled_histograms/exp_corrected/re0.00166_im-0.00198_sigmat
cd Subanalysis/InterfFunction/root_files/filled_sampled_histograms/exp_corrected/re0.00166_im-0.00198_sigmat1.00_t10.00_300.00_t20.00_300.00_tch0.00_300.00_tne0.00_300.00/
ls
ls -hl
pgrep hist_sample_d
ls -hl
pgrep hist_sample_d
cd ../../../
cd ..
nohup parallel -a root_files/filled_sampled_histograms/inputs1.txt -j 6 ./hist_sample_draw {} 1 0.1 > nohup_sigmat_0.1_1.log &
nohup parallel -a root_files/filled_sampled_histograms/inputs2.txt -j 6 ./hist_sample_draw {} 1 0.1 > nohup_sigmat_0.1_2.log &
nohup parallel -a root_files/filled_sampled_histograms/inputs3.txt -j 6 ./hist_sample_draw {} 1 0.1 > nohup_sigmat_0.1_3.log &
nohup parallel -a root_files/filled_sampled_histograms/inputs4.txt -j 6 ./hist_sample_draw {} 1 0.1 > nohup_sigmat_0.1_4.log &
nohup parallel -a root_files/filled_sampled_histograms/inputs5.txt -j 6 ./hist_sample_draw {} 1 0.1 > nohup_sigmat_0.1_5.log &
nohup parallel -a root_files/filled_sampled_histograms/inputs6.txt -j 6 ./hist_sample_draw {} 1 0.1 > nohup_sigmat_0.1_6.log &
ls
cd root_files/filled_sampled_histograms/exp_corrected/re0.00166_im-0.00198_sigmat0.10_t10.00_300.00_t20.00_300.00_tch0.00_300.00_tne0.00_300.00/
ls
ls -hl
ls
cd Subanalysis/InterfFunction/root_files/filled_sampled_histograms/exp_corrected/re0.00166_im-0.00198_sigmat1.00_t10.00_300.00_t20.00_300.00_tch0.00_300.00_tne0.00_300.00/
ls
hadd -f merged_2DHist.root *.root
root merged_2DHist.root 
cd ../../../../
./compile.sh 
./hist_fits_integral_range 
cd root_files/filled_sampled_histograms/exp_corrected/re0.00166_im-0.00198_sigmat0.10_t10.00_300.00_t20.00_300.00_tch0.00_300.00_tne0.00_300.00/
ls
ls -hl
hadd -f merged_2DHist.root *.root
cd ../../../
cd ..
./compile.sh 
./hist_fits_integral_range 
./compile.sh 
nohup parallel -a root_files/filled_sampled_histograms/inputs1.txt -j 6 ./hist_sample_draw {} 1 0.5 > nohup_sigmat_0.5_1.log &
nohup parallel -a root_files/filled_sampled_histograms/inputs2.txt -j 6 ./hist_sample_draw {} 1 0.5 > nohup_sigmat_0.5_2.log &
nohup parallel -a root_files/filled_sampled_histograms/inputs3.txt -j 6 ./hist_sample_draw {} 1 0.5 > nohup_sigmat_0.5_3.log &
nohup parallel -a root_files/filled_sampled_histograms/inputs4.txt -j 6 ./hist_sample_draw {} 1 0.5 > nohup_sigmat_0.5_4.log &
nohup parallel -a root_files/filled_sampled_histograms/inputs5.txt -j 6 ./hist_sample_draw {} 1 0.5 > nohup_sigmat_0.5_5.log &
nohup parallel -a root_files/filled_sampled_histograms/inputs6.txt -j 6 ./hist_sample_draw {} 1 0.5 > nohup_sigmat_0.5_6.log &
ls -hl root_files/filled_sampled_histograms/exp_corrected/re0.00166_im-0.00198_sigmat0.50_t10.00_300.00_t20.00_300.00_tch0.00_300.00_tne0.00_300.00/
ls -hl root_files/filled_sampled_histograms/exp_corrected/
nohup parallel -a root_files/filled_sampled_histograms/inputs1.txt -j 6 ./hist_sample_draw {} 1 0.8 > nohup_sigmat_0.8_1.log &
nohup parallel -a root_files/filled_sampled_histograms/inputs2.txt -j 6 ./hist_sample_draw {} 1 0.8 > nohup_sigmat_0.8_2.log &
nohup parallel -a root_files/filled_sampled_histograms/inputs3.txt -j 6 ./hist_sample_draw {} 1 0.8 > nohup_sigmat_0.8_3.log &
nohup parallel -a root_files/filled_sampled_histograms/inputs4.txt -j 6 ./hist_sample_draw {} 1 0.8 > nohup_sigmat_0.8_4.log &
nohup parallel -a root_files/filled_sampled_histograms/inputs5.txt -j 6 ./hist_sample_draw {} 1 0.8 > nohup_sigmat_0.8_5.log &
nohup parallel -a root_files/filled_sampled_histograms/inputs6.txt -j 6 ./hist_sample_draw {} 1 0.8 > nohup_sigmat_0.8_6.log &
git add .
git commit -m "Test"
git push
./execute_analysis.sh 
git add .
git commit -m "Test"
git push
./execute_analysis.sh 
git add .
git commit -m "Test"
git push
./execute_analysis.sh 
git add .
git push
git commit -m "Test"
git push
./execute_analysis.sh 
git add .
git commit -m "Test"
git push
git add .
git commit -m "Test"
git push
./execute_analysis.sh 
ls
cat job_v26_all_phys3_1_inv_pb_1.txt 
git add .
git commit -m "Test"
git push
cd build/bin/
./KLSPM00 ../../job_v26_all_phys3_1_inv_pb_1.txt 
cd ../../
./execute_analysis.sh 
cd build/bin/
./KLSPM00 ../../job_v26_all_phys3_1_inv_pb_1.txt 
pgrep sample
pgrep bttf
pgrep hist
git add .
git commit -m "Test"
cd ../../
cd ../
cd KLOE/
git add .
git commit -m "Test"
git push
nohup parallel -a root_files/filled_sampled_histograms/inputs1.txt -j 6 ./hist_sample_draw {} 1 1.5 > nohup_sigmat_1.5_1.log &
cd Subanalysis/InterfFunction/
./compile.sh 
nohup parallel -a root_files/filled_sampled_histograms/inputs1.txt -j 6 ./hist_sample_draw {} 1 1.5 > nohup_sigmat_1.5_1.log &
nohup parallel -a root_files/filled_sampled_histograms/inputs2.txt -j 6 ./hist_sample_draw {} 1 1.5 > nohup_sigmat_1.5_2.log &
nohup parallel -a root_files/filled_sampled_histograms/inputs3.txt -j 6 ./hist_sample_draw {} 1 1.5 > nohup_sigmat_1.5_3.log &
nohup parallel -a root_files/filled_sampled_histograms/inputs4.txt -j 6 ./hist_sample_draw {} 1 1.5 > nohup_sigmat_1.5_4.log &
nohup parallel -a root_files/filled_sampled_histograms/inputs5.txt -j 6 ./hist_sample_draw {} 1 1.5 > nohup_sigmat_1.5_5.log &
nohup parallel -a root_files/filled_sampled_histograms/inputs6.txt -j 6 ./hist_sample_draw {} 1 1.5 > nohup_sigmat_1.5_6.log &
htop
./compile.sh 
ls
./hist_fits_integral_range.exe 
./hist_fits_integral_range.exe exp_corrected re0.00166_im-0.00198_sigmat1.00_t10.00_300.00_t20.00_300.00_tch0.00_300.00_tne0.00_300.00
./compile.sh 
./hist_fits_integral_range.exe exp_corrected re0.00166_im-0.00198_sigmat1.00_t10.00_300.00_t20.00_300.00_tch0.00_300.00_tne0.00_300.00
./compile.sh 
./hist_fits_integral_range.exe exp_corrected re0.00166_im-0.00198_sigmat1.00_t10.00_300.00_t20.00_300.00_tch0.00_300.00_tne0.00_300.00
./hist_fits_integral_range.exe exp_corrected re0.00166_im-0.00198_sigmat1.00_t10.00_300.00_t20.00_300.00_tch0.00_300.00_tne0.00_300.00 test
./compile.sh 
./hist_fits_integral_range.exe exp_corrected re0.00166_im-0.00198_sigmat1.00_t10.00_300.00_t20.00_300.00_tch0.00_300.00_tne0.00_300.00 test
nohup parallel -a root_files/filled_sampled_histograms/inputs5.txt -j 6 ./hist_sample_draw {} 1 1.5 > nohup_sigmat_1.5_5.log &
nohup parallel -a root_files/filled_sampled_histograms/inputs.txt -j 6 ./hist_sample_draw.exe {} 1 1.5 > nohup_sigmat_1.5.log &
ls
chmod +x exec_int_range.sh 
./exec_int_range.sh 
./exec_int_range.sh exp_corrected re0.00166_im-0.00198_sigmat1.00_t10.00_300.00_t20.00_300.00_tch0.00_300.00_tne0.00_300.00 150_points_both_free_sigma_1.0 0 300 0.5 300 150
cat nohup_150_points_both_free_sigma_1.0.log 
./exec_int_range.sh exp_corrected re0.00166_im-0.00198_sigmat1.00_t10.00_300.00_t20.00_300.00_tch0.00_300.00_tne0.00_300.00 150_points_re_fixed_sigma_1.0 0 300 0.5 300 150
./exec_int_range.sh exp_corrected re0.00166_im-0.00198_sigmat1.00_t10.00_300.00_t20.00_300.00_tch0.00_300.00_tne0.00_300.00 150_points_im_fixed_sigma_1.0 0 300 0.5 300 150
./exec_int_range.sh exp_corrected re0.00166_im-0.00198_sigmat0.10_t10.00_300.00_t20.00_300.00_tch0.00_300.00_tne0.00_300.00 150_points_im_fixed_sigma_0.1 0 300 0.5 300 150
./exec_int_range.sh exp_corrected re0.00166_im-0.00198_sigmat0.10_t10.00_300.00_t20.00_300.00_tch0.00_300.00_tne0.00_300.00 150_points_re_fixed_sigma_0.1 0 300 0.5 300 150
./exec_int_range.sh exp_corrected re0.00166_im-0.00198_sigmat0.10_t10.00_300.00_t20.00_300.00_tch0.00_300.00_tne0.00_300.00 150_points_both_free_sigma_0.1 0 300 0.5 300 150
./exec_int_range.sh exp_corrected re0.00166_im-0.00198_sigmat0.50_t10.00_300.00_t20.00_300.00_tch0.00_300.00_tne0.00_300.00 150_points_both_free_sigma_0.5 0 300 0.5 300 150
./exec_int_range.sh exp_corrected re0.00166_im-0.00198_sigmat0.80_t10.00_300.00_t20.00_300.00_tch0.00_300.00_tne0.00_300.00 150_points_both_free_sigma_0.8 0 300 0.5 300 150
./exec_int_range.sh exp_corrected re0.00166_im-0.00198_sigmat0.80_t10.00_300.00_t20.00_300.00_tch0.00_300.00_tne0.00_300.00 150_points_re_fixed_sigma_0.8 0 300 0.5 300 150
./exec_int_range.sh exp_corrected re0.00166_im-0.00198_sigmat0.50_t10.00_300.00_t20.00_300.00_tch0.00_300.00_tne0.00_300.00 150_points_re_fixed_sigma_0.5 0 300 0.5 300 150
./exec_int_range.sh exp_corrected re0.00166_im-0.00198_sigmat0.50_t10.00_300.00_t20.00_300.00_tch0.00_300.00_tne0.00_300.00 150_points_im_fixed_sigma_0.5 0 300 0.5 300 150
./exec_int_range.sh exp_corrected re0.00166_im-0.00198_sigmat0.80_t10.00_300.00_t20.00_300.00_tch0.00_300.00_tne0.00_300.00 150_points_im_fixed_sigma_0.8 0 300 0.5 300 150
cd img/integral_range/150_points_both_free_sigma_0.8/
ls
cd ../../
cd ..
cat nohup_150_points_both_free_sigma_0.5.log 
cd root_files/filled_sampled_histograms/
ls
cd exp_corrected/
ls
cd re0.00166_im-0.00198_sigmat0.50_t10.00_300.00_t20.00_300.00_tch0.00_300.00_tne0.00_300.00/
hadd merged_2DHist.root *.root
cd ../re0.00166_im-0.00198_sigmat0.80_t10.00_300.00_t20.00_300.00_tch0.00_300.00_tne0.00_300.00/
hadd merged_2DHist.root *.root
cd ../../
cd ../../../
cd InterfFunction/
./exec_int_range.sh exp_corrected re0.00166_im-0.00198_sigmat0.50_t10.00_300.00_t20.00_300.00_tch0.00_300.00_tne0.00_300.00 150_points_im_fixed_sigma_0.5 0 300 0.5 300 150
./exec_int_range.sh exp_corrected re0.00166_im-0.00198_sigmat0.80_t10.00_300.00_t20.00_300.00_tch0.00_300.00_tne0.00_300.00 150_points_im_fixed_sigma_0.8 0 300 0.5 300 150
./exec_int_range.sh exp_corrected re0.00166_im-0.00198_sigmat0.50_t10.00_300.00_t20.00_300.00_tch0.00_300.00_tne0.00_300.00 150_points_re_fixed_sigma_0.5 0 300 0.5 300 150
./exec_int_range.sh exp_corrected re0.00166_im-0.00198_sigmat0.80_t10.00_300.00_t20.00_300.00_tch0.00_300.00_tne0.00_300.00 150_points_re_fixed_sigma_0.8 0 300 0.5 300 150
./exec_int_range.sh exp_corrected re0.00166_im-0.00198_sigmat0.50_t10.00_300.00_t20.00_300.00_tch0.00_300.00_tne0.00_300.00 150_points_both_free_sigma_0.5 0 300 0.5 300 150
./exec_int_range.sh exp_corrected re0.00166_im-0.00198_sigmat0.80_t10.00_300.00_t20.00_300.00_tch0.00_300.00_tne0.00_300.00 150_points_both_free_sigma_0.8 0 300 0.5 300 150
cd Subanalysis/InterfFunction/root_files/filled_sampled_histograms/exp_corrected/re0.00166_im-0.00198_sigmat1.50_t10.00_300.00_t20.00_300.00_tch0.00_300.00_tne0.00_300.00/
ls
ls -hl
hadd merged_2DHist.root *.root
ls -hl
cd ../
cd ../../
./exec_int_range.sh exp_corrected re0.00166_im-0.00198_sigmat1.50_t10.00_300.00_t20.00_300.00_tch0.00_300.00_tne0.00_300.00 150_points_both_free_sigma_1.5 0 300 0.5 300 150
./exec_int_range.sh exp_corrected re0.00166_im-0.00198_sigmat1.50_t10.00_300.00_t20.00_300.00_tch0.00_300.00_tne0.00_300.00 150_points_re_fixed_sigma_1.5 0 300 0.5 300 150
./exec_int_range.sh exp_corrected re0.00166_im-0.00198_sigmat1.50_t10.00_300.00_t20.00_300.00_tch0.00_300.00_tne0.00_300.00 150_points_im_fixed_sigma_1.5 0 300 0.5 300 150
ssh tier1-cnaf 
cd build/bin/
gdb KLSPM00 
valgrind KLSPM00_EXECUTE 
valgrind KLSPM00
cd ../../
./execute_analysis.sh 
llq
ll
find /usr/lib64 /usr/lib -name "libasan.so*"
cd build/
rm -fr *
cd ..
./execute_analysis.sh 10
./execute_analysis.sh
./execute_analysis.sh 
cd build/bin/
gdb ./KLSPM00 
cd ../../
./execute_analysis.sh 
export ASAN_OPTIONS=print_stacktrace=1
./execute_analysis.sh 
git add .
git commit -m "Fix"
git push
./execute_analysis.sh 
git checkout -- .
./execute_analysis.sh 
export ASAN_OPTIONS=symbolize=1
./execute_analysis.sh 
cd build/bin/
valgrind --tool=memcheck          --leak-check=full          --show-leak-kinds=all          --track-origins=yes          --verbose          ./KLSPM00
valgrind --tool=memcheck          --suppressions=$ROOTSYS/etc/valgrind-root.supp          --leak-check=full          ./KLSPM00
cd ../../
./execute_analysis.sh 
git add .
git commit -m "Fix"
git pus
git push
./execute_analysis.sh 
cd Subanalysis/Properties/
code properties.json 
code analysis_config.json 
cd ../../
./execute_analysis.sh 
cat ~/.bashrc 
cd ..
cd root_files/
ls
cd kitt/
ls
exit
cd ../root_files/kitt/
ls
cd ../hal/
ls
cd ../kitt/
ls
cd ALL_PHYS3_SIGNAL_NoSmearing/
ls
ls -hl
cd ..
rm -fr *
cd ..
cd kitt/
ls
rm -fr *
cd ..
cd kitt/
ls
ALL_PHYS3_SIGNAL_NoSmearing/
ls
cd ALL_PHYS3_SIGNAL_NoSmearing/
ls
ls -hl
git add .
git commit -m "Fix"
git push
git pull
vi ~/.bashrc 
source ~/.bashrc 
./execute_analysis.sh 
git add .
git commit -m "Custom root folder"
git push
vi ~/.bashrc 
ls
cd ..
ls
cd root_files/
ls
cd ..
cd DBV-26/
ls
cd all_phys
ls
cd 20260315/
ls
cd ../
ls
cd ..
ls
cd ..
ls
./lumi_batch.sh 
vi lumi_batch.sh 
./lumi_batch.sh PROD2ROOT/MC/MK0 
./lumi_batch.sh PROD2ROOT/MC/MK0 26 all_phys3 4 lumi_file/lumi_per_ev_nb.log 
ls
cd root_files/
ls
cd hal/
ls
cd 2026-04-14/
ls
cd ..
cd kitt/
ls
cd 2026-04-14/
ls
cd ..
rm -fr *
cd ../hal/
rm -fr *
ll
ls -hl
cd ..
ls
cd kitt/
ls
cd ..
cd ../
cd DBV-26/all_phys3/
ls
cd 20260414
ls
cd ../
cd ../../
cd DBV-26/all_phys3/20260414/
ls
cat job_v26_all_phys3_4_inv_pb_1.txt 
cd ../../
cd ../
./lumi_batch.sh /data/ssd/gamrat/PROD2ROOT/MC/MK0 all_phys3 
./lumi_batch.sh /data/ssd/gamrat/PROD2ROOT/MC/MK0 26 all_phys3 4 lumi_file/lumi_per_ev_nb.log 
cd DBV-26/
cd all_phys3/
rm -r 20260414/
cd ../../
./lumi_batch.sh /data/ssd/gamrat/PROD2ROOT/MC/MK0 26 all_phys3 4 lumi_file/lumi_per_ev_nb.log 
cd DBV-26/all_phys3/20260414/
ls
cd ../../
ls
cd ../
./lumi_batch.sh /data/ssd/gamrat/PROD2ROOT/MC/MK0 26 all_phys2 4 lumi_file/lumi_per_ev_nb.log 
cd KLOE/
./run_parallel.sh 
./run_parallel.sh all_phys2 20260415 1 5
ls
cd ..
llq
htop
cd root_files/hal/
ls
cd ..
cd KLOE/
rm -r para
rm -r parallel_logs/
./run_parallel.sh all_phys2 20260415 1 5
cd ../
htop
cd root_files/hal/
ls
cd ..
cd kitt/
ls
cd ../
cd ../KLOE/parallel_logs/
ls
cd 1/1/
ls
cat seq
cat stderr
source ~/.bashrc
cd ../../
cd ..
rm -r parallel_logs/
./run_parallel.sh all_phys2 20260415 1 5
cd parallel_logs/1/1
cat stderr 
cd ../../
pgrep KLSPM00
cd ..
rm -r parallel_logs/
vi run_parallel.sh 
./run_parallel.sh all_phys2 20260415 1 5
cd parallel_logs/1/1/
cat stderr 
cd ../../
cd ../
cd ..
ls
cd root_files/hal/ALL_PHYS2_SIGNAL_NoSmearing/
ls
ls -hl
cd ../../
cd kitt/ALL_PHYS3_SIGNAL_NoSmearing/
ls -hl
cd ../../
cd ../
cd KLOE/Subanalysis/Properties/
cat analysis_config.json 
cd ../../
cd ../
./lumi_batch.sh /data/ssd/gamrat/PROD2ROOT/MC/MK0 26 all_phys 4 lumi_file/lumi_per_ev_nb.log 
cd root_files/kitt/
ls
cd ALL_PHYS
cd ALL_PHYS_SIGNAL_NoSmearing/
ls -hl
htop
ls /data/tape_dump/gamrat/PROD2ROOT/MC/MK0/
./lumi_batch.sh /data/tape_dump/gamrat/PROD2ROOT/MC/MK0 26 all_phys 4 lumi_file/lumi_per_ev_nb.log 
cd ../../
cd ..
./lumi_batch.sh /data/tape_dump/gamrat/PROD2ROOT/MC/MK0 26 all_phys 4 lumi_file/lumi_per_ev_nb.log 
nproc
lscpu
cd KLOE/
nproc --all
vi run_parallel.sh 
./run_parallel.sh all_phys 20260415 2 80 
cd ../root_files/hal/ALL_PHYS
cd ../root_files/hal/ALL_PHYS_SIGNAL_NoSmearing/
ls
llq
ll
cd ../
htop
cd ..
cd ../KLOE/
ls
cd parallel_logs_all_phys_2_80/
ls
cd 1/1
ls
cd 1
ls
cd 10/
ls
cat stdout
exit
ls -hl
cd parallel_logs_all_phys_2_80/
ls
cd 1/
ls
cd 10
ls
cat stdout 
exit
cd ../root_files/hal/
ls
cd ALL_PHYS2_SIGNAL_NoSmearing/
ls -hl
cd ../ALL_PHYS_SIGNAL_NoSmearing/
ls -hl
llq
cd ../../../KLOE/
cd parallel_logs_all_phys_2_80/
cd 1/13/
ls
cat stdout 
cat stderr 
cd ../../../
cd ../
cd root_files/
ls
cd hal/
ls -hl
cd ALL_PHYS_SIGNAL_NoSmearing/
ls
ls -hl
exit
vi run_parallel.sh 
./execute_analysis.sh 
./execute_analysis.sh 8 ../DBV-26/all_phys/20260415/job_v26_all_phys_4_inv_pb_13.txt 
./execute_analysis.sh 8 ../../../DBV-26/all_phys/20260415/job_v26_all_phys_4_inv_pb_13.txt 
./execute_analysis.sh 8 ../../DBV-26/all_phys/20260415/job_v26_all_phys_4_inv_pb_13.txt 
./execute_analysis.sh 8 ../../DBV-26/all_phys/20260415/job_v26_all_phys_4_inv_pb_54.txt 
./execute_analysis.sh 8 ../../DBV-26/all_phys/20260415/job_v26_all_phys_4_inv_pb_55.txt 
./execute_analysis.sh 8 ../../DBV-26/all_phys/20260415/job_v26_all_phys_4_inv_pb_34.txt 
./execute_analysis.sh 
./execute_analysis.sh 8 ../../DBV-26/all_phys/20260415/job_v26_all_phys_4_inv_pb_55.txt 
killall KLSPM00
cd ../root_files/
ls
cd hal/
rm -fr *
cd ../kitt/
rm -fr *
cd ..
cd ../KLOE/
./execute_analysis.sh 
git add .
git commit -m "Fixed fast(free)"
git restore --staged .
git reset 
git rm -r --cached parallel_logs*
git add .
git commit -m "Fixed fast(free)"
git push
rm -r parallel_logs*
./run_parallel.sh all_phys2 20260415 1 5
pgrep KLSPM00
cd ../root_files/hal/
ls
cd ALL_PHYS2_SIGNAL_NoSmearing/
ls -hl
cd ../
cd ../../
cd KLOE/
./run_parallel.sh all_phys 20260415 1 100
htop
killall KLSPM00
cd ../
cd ki
cd root_files/kitt/
ls -hl
cd ../
cd hal/
rm -fr *
cd ../../KLOE/
./execute_analysis.sh 
./run_parallel.sh all_phys2 20260415 1 5
./run_parallel.sh all_phys2 20260415 1 5 3
./run_parallel.sh all_phys 20260415 1 100 50
htop
cd ../root_files/hal/
killall KLSPM00
rm -fr *
cd ../../
cd KLOE/
./run_parallel.sh all_phys2 20260415 1 5 3
htop
pgrep KLSPM00
killall parallel
killall nohup
ps -e
killall execute_analysis
killall KLSPM00
htop
killall KLSPM00
cd ../
cd hal
cd root_files/hal/
ls
rm -fr *
ls
cd ..
cd kitt/
ls
cd ../../
cd KLOE/
ls
./run_parallel.sh all_phys2 20260415 1 5 5
killall KLSPM00
killall execute_analysis
rm parallel_logs_all_phys*
rm -r parallel_logs_all_phys*
./run_parallel.sh all_phys2 20260415 1 5 5
./run_parallel.sh all_phys 20260415 1 100 45
htop
git add .
git commit -m "Possibility to choose nproc"
git push
cd ../root_files/kitt/
ls
cd ALL_PHYS3_SIGNAL_NoSmearing/
ls
ls -hl
cd ../../hal/
ls
cd ALL_PHYS
cd ALL_PHYS_SIGNAL_NoSmearing/
ls -hl
vi ~/.bashrc
ls -hl
cd ../root_files/hal/ALL_PHYS_SIGNAL_NoSmearing/
ls -hl
htop
cd../../../
cd ../
cd ../../
cd KLOE/
cd parallel_logs_all_phys_1_100/
cd 1/9/
cat stderr
cat stdout 
exit
ls -hl
vi ~/.bashrc 
ls Subanalysis/Properties/
htop
source ~/.bashrc 
vi ~/.bashrc 
./execute_analysis.sh 
./execute_analysis.sh 5
git add .
git commit -m "Customizable analysis config"
git push
htop
vi ~/.bashrc
cat ~/.bashrc 
ls
cd ..
ls
cd CNAF_Produced_Files/
ls
cd root_files/
ls
rm -fr *
ls
cd ..
rm -fr root_files/
scp -r tier1-cnaf:~/root_files .
llq
ll
ls -hl
cd root_files/
ls -hl
cd 2026-04-17/
ls
cd Signal/
ls
cd ALL_PHYS_SIGNAL_NoSmearing/
ls
exit
cd ../CNAF_Produced_Files/root_files/2026-04-17/Signal/
ls
cd DATA_SIGNAL_NoSmearing/
ls
cd ../
ls
cd ALL_PHYS_SIGNAL_NoSmearing/
ls -l
cd ..
cd ALL_PHYS2_SIGNAL_NoSmearing/
ls -l
cd ../ALL_PHYS3_SIGNAL_NoSmearing/
ls
cd ..
ls
cd ../CNAF_Produced_Files/root_files/2026-04-17/
ls
cd Signal/
ls
cd ALL_PHYS_SIGNAL_NoSmearing/
ls
cd ..
cd log/
ls
cd all_phys
ls
cd 1
ls
cat cut.analysis.log 
cd ../
ls
cd ..
ls
vi calculate_purity_eff.sh
chmod +x calculate_purity_eff.sh 
./calculate_purity_eff.sh 
vi calculate_purity_eff.sh
rm calculate_purity_eff.sh 
vi calculate.sh
chmod +x calculate.sh 
./calculate.sh 
cd all_phys/1
cat analysis.config.log 
cat cut.analysis.log 
cd ../../
rm calculate.sh 
vi calculate.sh
chmod +x calculate.sh 
./calculate.sh 
rm calculate.sh 
vi calculate.sh
chmod +x calculate.sh 
./calculate.sh 
vi calculate.sh
rm calculate.sh 
vi calculate.sh
chmod +x calculate.sh 
./calculate.sh 
cd data/1/
cat cut.analysis.log 
cd ../../
ls
cd data/
ls
cd 1
ls
cat cut.analysis.log 
cd ../../
vi calculate_total.sh
chmod +x calculate_total.sh 
./calculate_total.sh 
./calculate.sh 
ls
cp *.sh ../../4pi/log/.
cd ../../4pi/log/
./calculate
./calculate.sh 
vi calculate.sh 
vi calculate_total.sh 
./calculate.sh 
./execute_analysis.sh 
./execute_analysis.sh 
root
./execute_analysis.sh 
