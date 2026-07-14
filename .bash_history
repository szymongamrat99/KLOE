git add .
git commit -m "New approach"
git push
git checkout -b feat/mixed-channel-constants
git push --set-upstream origin feat/mixed-channel-constants 
./execute_analysis.sh 
git add .
git commit -m "Dynamic variables for fitting in the cp_fit_final"
git add .
git commit -m "change"
./execute_analysis.sh 
git add .
git commit -m "Changes for interference"
./execute_analysis.sh 
git checkout -- .
./execute_analysis.sh 
git add .
git commit -m "Flexible normalization"
git push
./execute_analysis.sh 
cd ../root_files/hal/2026-05-05/
ls
cd ALL_PHYS3_SEMILEPTONIC_NoSmearing/
ls
ls -hl
root
cd ../ALL_PHYS3_THREE_PI0_NoSmearing/
ls
root
ls
cd
ls
cd /data/ssd/gamrat/root_files/hal/2026-05-05/ALL_PHYS3_THREE_PI0_NoSmearing/
ls
ls -hl
rm mk0_initial_analysis_all_phys3_THREE_PI0_NoSmearing_6.root
rm mk0_initial_analysis_all_phys3_THREE_PI0_NoSmearing_7.root
rm mk0_initial_analysis_all_phys3_THREE_PI0_NoSmearing_8.root
rm mk0_initial_analysis_all_phys3_THREE_PI0_NoSmearing_9.root
rm mk0_initial_analysis_all_phys3_THREE_PI0_NoSmearing_10.root
rm mk0_initial_analysis_all_phys3_THREE_PI0_NoSmearing_11.root
rm mk0_initial_analysis_all_phys3_THREE_PI0_NoSmearing_12.root
ls -hl
root
git add .
git checkout feat/mixed-channel-constants 
git add .
git commit -m "Added new hypotheses"
./execute_analysis.sh 
nohup ./execute_analysis.sh < 3 1 2 > nohup.log &
nohup ./execute_analysis.sh < '3 1 2' > nohup.log &
nohup ./execute_analysis.sh < "3 1 2" > nohup.log &
ls
cat parameters.txt 
code parameters.txt 
nohup ./execute_analysis.sh < parameters.txt > nohup.log &
cat nohup
cat nohup.log 
ps -e
pgrep KLSPM00
kill 30879
nohup ./execute_analysis.sh < parameters.txt > nohup.log &
cat nohup.log 
nohup ./execute_analysis.sh < parameters.txt > nohup1.log &
cat nohup1.log 
cat nohup.log 
git add .
git commit -m "Prepared for control sample choice"
git push
git checkout mast
git checkout maste
git checkout main 
git merge feat/mixed-channel-constants 
git push
cd build/
rm -fr *
cd ..
./execute_analysis.sh 
cd log/
ls
cd 2026-05-05
ls
cat error.log 
./execute_analysis.sh 
git add .
git commit -m "Final version"
git push
git checkout main
git checkout -- .
git checkout main 
git pull
./execute_analysis.sh 
ssh tier1-cnaf 
exit
./execute_analysis.sh 
git add .
git checkout -b feat/scaling-factors-regeneration
git push
git push --set-upstream origin feat/scaling-factors-regeneration
git add .
git commit -m "Preparation for scaling fctors"
git push
./execute_analysis.sh 
git add .
git commit -m "Scaling factors added to fitting"
git push
./execute_analysis.sh 
git add .
git commit -m "Scaling factors propagated along with their errors"
git push
./execute_analysis.sh 
.q
./execute_analysis.sh 
./e
./execute_analysis.sh 
cd ../root_files/
ls
cd kitt/
ls
rm -fr ALL_PHYS*
cd ..
ls
cd ..
ls
cd gamrat/
ls
cd DBV-26/all_phys
ls
cd 20260415/
ls
cd ../../
ls
cd ../
cd root_files/
ls
cd kitt/
ls
cd 2026-04-15/
ls
cd ..
cd ../../
ls
cd DBV-26/all_phys2/
ls
cd 20260415/
ls -hl
cd ..
ks
ls
cd ../
ls
cd all_phys3/
ls
cd 20260414/
ls
cd ../../DK0/
ls
cd 20260315/
ls -hl
cd ../
ls
cd ..
ls
cd ..
ls
cd DBV-26/
ls
cd all_phys
ls
cd 20260315/
ls
cd ..
ls
cd ../root_files/kitt/
ls
cd 2026-04-15/
ls
cd ALL_PHYS3_SIGNAL_NoSmearing/
ls -hl
cd ..
ls
cd ..
cd h
cd ../
cd ..
ls
./lumi_
./lumi_batch.sh 
./lumi_batch.sh /data/tape_dump/gamrat/PROD2ROOT/MC/MK0 all_phys 4 lumi_file/lumi_per_ev_nb.log 
./lumi_batch.sh /data/tape_dump/gamrat/PROD2ROOT/MC/MK0 26 all_phys 4 lumi_file/lumi_per_ev_nb.log 
ls /data/tape_dump/gamrat/PROD2ROOT/MC/MK0/all_phys2/
ls /data/tape_dump/gamrat/PROD2ROOT/MC/MK0/
./lumi_batch.sh /data/tape_dump/gamrat/PROD2ROOT/MC/MK0 26 all_phys2 4 lumi_file/lumi_per_ev_nb.log 
./lumi_batch.sh /data/tape_dump/gamrat/PROD2ROOT/MC/MK0 26 all_phys3 4 lumi_file/lumi_per_ev_nb.log 
./lumi_batch.sh /data/ssd/gamrat/PROD2ROOT 26 DK0 4 lumi_file/lumi_per_ev_nb.log 
./lumi_batch.sh /data/ssd/gamrat/PROD2ROOT/DK0 26 DK0 4 lumi_file/lumi_per_ev_nb.log 
vi lumi_batch.sh 
./lumi_batch.sh /data/ssd/gamrat/PROD2ROOT 26 DK0 4 lumi_file/lumi_per_ev_nb.log 
vi lumi_batch.sh 
cd DBV-26/all_phys/20260516/
ls
cd ../
cd ../
ls PROD2ROOT/MC/MK0/all_phys
ls PROD2ROOT/MC/MK0/all_phys2
ls PROD2ROOT/MC/MK0/all_phys3
ls PROD2ROOT/DK0/
cd DBV-26/DK0/
ls
cd 20260516/
ls
cd ..
cd ../
cd all_phys
ls
cd 20260516/
ls
cd ../
cd ../DK0/
ls
cd 20260516/
ls
cat luminosity_report.log 
cd ../../../PROD2ROOT/DK0/
ls
cd ../DBV-26/all_phys2/
ls -hl
cd 20260516/
ls
cd ../../
cd all_phys3/
ls
cd 20260516/
ls -hl
ls
./execute_analysis.sh 
htop
git checkout 
git add .
git commit -m "Scalingfactors almost done"
git push
git checkout 
git checkout main 
git pull
./execute_analysis.sh 
htop
./run_parallel.sh 
./run_parallel.sh all_phys 20260516 1 206 15
htop
./run_parallel.sh all_phys2 20260516 1 71 15
htop
./run_parallel.sh all_phys3 20260516 1 76 15
htop
ls
pgrep parallel
pgrep run
pgrep nohup
ps -e
killall parallel
killall nohup
killall KLSPM00
killall execute_analysis.sh
ps -e
killall bash
ps -e
killall KLSPM00
killall execute_analysis.sh
killall KLSPM00
killall execute_analysis.sh
killall KLSPM00
killall execute_analysis.sh
killall KLSPM00
killall execute_analysis.sh
killall KLSPM00
killall execute_analysis.sh
killall KLSPM00
rm -fr parallel_logs_all_phys*
rm -fr nohup_all_phys*
htop
./run_parallel.sh all_phys 20260516 1 206 15
./run_parallel.sh all_phys2 20260516 1 71 15
./run_parallel.sh all_phys3 20260516 1 76 15
cd ../root_files/hal/
ls
cd ALL_PHYS_SEMILEPTONIC_NoSmearing/
ls
ls -hl
cd ..
ls
htop
cd ../root_files/hal/ALL_PHYS_SEMILEPTONIC_NoSmearing/
ls
ls -hl
htop
du -h
df -h
exit
df -h
quota
quota --human-readable 
cd ../root_files/hal
ls
cd ALL_PHYS_SEMILEPTONIC_NoSmearing/
ls
ls -hl
htop
cd ../
cd ALL_PHYS3_SEMILEPTONIC_NoSmearing/
ls -hl
ll
ls -hl
htop
cd ../../../
cd KLOE/
ls -hl
cat nohup_all_phys_1_206.log 
ls
cd ../root_files/hal/
ls
cd ALL_PHYS_SEMILEPTONIC_NoSmearing/
ls
ls -hl
ls
cat file_lumi_ALL_PHYS_SEMILEPTONIC_NoSmearing.log 
root 
root mk0_initial_analysis_all_phys_SEMILEPTONIC_NoSmearing_12.root
htop
cd ../root_files/hal/
ls
cd ALL_PHYS2_SEMILEPTONIC_NoSmearing/
ls
ls -hl
;s
ls
cd ..
ls
cd ALL_PHYS3_SEMILEPTONIC_NoSmearing/
ls
ls -hl
cd ../
cd ALL_PHYS_SEMILEPTONIC_NoSmearing/
ls
ls -hl
ls
ls -l
ls
ls | wc -l
ls -hl
cd ../
cd ..
cd ../KLOE/
ls
cd parallel_logs_all_phys_1_206/
ls
cd 1
ls
cd 1
ls
cat stdout 
ls
cd ..
ls
cd ..
cd parallel_logs_all_phys2_1_71/
ls
cd 1
ls
cd ../
ls
cd ../log/
ls
cd 2026-05-16
ls
cd ..
cd all_phys
ls
cd 1
ls
cat cut.analysis.log 
cd ..
ls
cd 110
ls
cat error.log 
cat general.log 
ls
cat analysis.config.log 
exit
cd ../root_files/hal/ALL_PHYS_SEMILEPTONIC_NoSmearing/
ls
ls -hl
cd ../ALL_PHYS2_SEMILEPTONIC_NoSmearing/
root
git checkout -b feat/semileptonic-reconstruction-method
cd ../PROD2ROOT/MC/
ls
cd MK0/
ls
cd all_phys
ls
root prod2root_mk0_all_phys_32508_v2.root 
cd ../root_files/hal/
ls
cd ALL_PHYS2_SEMILEPTONIC_NoSmearing/
ls -h
ls -hl
cp ../root_files/initial_analysis/lib/HistFactory.h Include/Codes/inc/.
cp ../root_files/initial_analysis/lib/HistFactory.cpp Include/Codes/src/.
./execute_analysis.sh 
htop
./execute_analysis.sh 
pgrep KLSPM00
killall KLSPM00
pgrep KLSPM00
killall KLSPM00
rm -fr parallel_logs_all_phys*
cd ../
cd root_files/
ls
cd hal/
ls
rm -fr *
cd ../../
ls
cd KLOE/
ls
rm nohup_all_phys*
./run_parallel.sh 
./run_parallel.sh all_phys 1 206 15
./run_parallel.sh all_phys 20260516 1 206 15
./run_parallel.sh all_phys2 20260516 1 71 15
./run_parallel.sh all_phys3 20260516 1 76 15
cd initial_analysis/src/
ls
root
ls
cd ../root_files/hal/ALL_PHYS
cd ../root_files/hal/ALL_PHYS_SEMILEPTONIC_NoSmearing/
ls
ls -hl
htop
cd ../root_files/
ls
cd kitt/
ls
cd 2026-05-18/
ls
cd ALL_PHYS3_THREE_PI0_NoSmearing/
ls
cd ..
rm -fr *
cd ..
htop
cd kitt/
ls
cd ..
ls
cd hal/
ls
cd ALL_PHYS
ls
cd ALL_PHYS_SEMILEPTONIC_NoSmearing/
ls
ls -hl
htop
ls -hl
cat input_luminosity_ALL_PHYS_SEMILEPTONIC_NoSmearing.log 
cat file_lumi_ALL_PHYS_SEMILEPTONIC_NoSmearing.log 
git add .
git commit -m "Reconstruction method adjusted"
git push
git push --set-upstream origin feat/semileptonic-reconstruction-method
ls ../root_files/hal/ALL_PHYS_SEMILEPTONIC_NoSmearing/
ls -hl ../root_files/hal/ALL_PHYS_SEMILEPTONIC_NoSmearing/
cd ../root_files/hal/
ls
cd 2026-05-18/
ls
cd ALL_PHYS3_SEMILEPTONIC_NoSmearing/
ls
ls -hl
root
./execute_analysis.sh 
./execute_analysis.sh 
cd ../root_files/hal/2026-05-19/
ls
cd ALL_PHYS3_SEMILEPTONIC_NoSmearing/
root
ls
ls -hl
root 
htop
cd ../../../
cd hal/ALL_PHYS2_SEMILEPTONIC_NoSmearing/
ls -hl
git add .
git commit -m "New semileptonic method done"
git push
ls ../root_files/kitt/
ls
ls ../root_files/hal/ALL_PHYS_SEMILEPTONIC_NoSmearing/
exit
cd Subanalysis/CPFit/config/
ls
ssh tier1-cnaf 
scp tier1-cnaf:/home/g/gamrat/init_analysis/scripts/regeneration_correction_functions.root .
scp tier1-cnaf:/home/g/gamrat/init_analysis/scripts/regeneration_fit_results.txt .
scp tier1-cnaf:/home/g/gamrat/init_analysis/scripts/regeneration_correction_functions.root .
scp tier1-cnaf:/home/g/gamrat/init_analysis/scripts/regeneration_fit_results.txt .
scp tier1-cnaf:/home/g/gamrat/init_analysis/scripts/regeneration_correction_functions.root .
scp tier1-cnaf:/home/g/gamrat/init_analysis/scripts/regeneration_fit_results.txt .
scp tier1-cnaf:/home/g/gamrat/init_analysis/scripts/regeneration_correction_functions.root .
scp tier1-cnaf:/home/g/gamrat/init_analysis/scripts/regeneration_fit_results.txt .
scp tier1-cnaf:/home/g/gamrat/init_analysis/scripts/regeneration_correction_functions.root .
scp tier1-cnaf:/home/g/gamrat/init_analysis/scripts/regeneration_fit_results.txt .
scp tier1-cnaf:/home/g/gamrat/init_analysis/scripts/regeneration_correction_functions.root .
scp tier1-cnaf:/home/g/gamrat/init_analysis/scripts/regeneration_fit_results.txt .
scp tier1-cnaf:/home/g/gamrat/init_analysis/scripts/regeneration_correction_functions.root .
scp tier1-cnaf:/home/g/gamrat/init_analysis/scripts/regeneration_fit_results.txt .
scp tier1-cnaf:/home/g/gamrat/init_analysis/scripts/regeneration_correction_functions.root .
scp tier1-cnaf:/home/g/gamrat/init_analysis/scripts/regeneration_fit_results.txt .
scp tier1-cnaf:/home/g/gamrat/init_analysis/scripts/regeneration_correction_functions.root .
scp tier1-cnaf:/home/g/gamrat/init_analysis/scripts/regeneration_fit_results.txt .
scp tier1-cnaf:/home/g/gamrat/init_analysis/scripts/regeneration_correction_functions.root .
scp tier1-cnaf:/home/g/gamrat/init_analysis/scripts/regeneration_fit_results.txt .
scp tier1-cnaf:/home/g/gamrat/init_analysis/scripts/regeneration_correction_functions.root .
scp tier1-cnaf:/home/g/gamrat/init_analysis/scripts/regeneration_fit_results.txt .
exit
htop
cd ../root_files/hal/
ls
cd ALL_PHYS2_SEMILEPTONIC_NoSmearing/
ls
ls -hl
cd ../../../
cd KLOE/
git add .
git commit -m "New function"
git push
./execute_analysis.sh 
cd Subanalysis/CPFit/img/
ls
ls -l
cd ../../../
./execute_analysis.sh 
.q
./execute_analysis.sh 
.q
./execute_analysis.sh 
cd Subanalysis/CPFit/config/
scp tier1-cnaf:/home/g/gamrat/init_analysis/scripts/regeneration_correction_functions.root .
scp tier1-cnaf:/home/g/gamrat/init_analysis/scripts/regeneration_fit_results.txt .
git pull
cd ..
git clone https://github.com/szymongamrat99/python-kloe-analysis .
git clone https://github.com/szymongamrat99/python-kloe-analysis
cd python-kloe-analysis/
ls
cd scripts/
ls
cd results/
ls
cd ..
root
ls
root
root-config --incdir
root
ls ../../CNAF_Produced_Files/root_files/2026-04-17/Signal/
root
cd ..
git add .
git commit -m "Created weighting matrix for the regeneration"
git push
cd scripts/
root
git add .
cd ..
git commit -m "Added new"
git push
cd scripts/
root
root-config
root-config --incdir
root
cd Subanalysis/CPFit/config/
cp /data/ssd/gamrat/python-kloe-analysis/scripts/results/* .
cd../../
cd ..
cd ../../\
cd ../../
./execute_analysis.sh 
cd ../../
cd gamrat/KLOE/
./execute_analysis.sh 
.q
./execute_analysis.sh 
cd scripts/
root
cd scripts/
ls
root
cp results/* ../../KLOE/Subanalysis/CPFit/config/.
root
.q
root
cp results/* ../../KLOE/Subanalysis/CPFit/config/.
root
cp results/* ../../KLOE/Subanalysis/CPFit/config/.
root
cp results/* ../../KLOE/Subanalysis/CPFit/config/.
root
cp results/* ../../KLOE/Subanalysis/CPFit/config/.
root
cp results/* ../../KLOE/Subanalysis/CPFit/config/.
root
cp results/* ../../KLOE/Subanalysis/CPFit/config/.
root
./execute_analysis.sh 
.q
./execute_analysis.sh 
git add .
git commit -m "Single regeneration scaling"
git push
./execute_analysis.sh 
git add .
git commit -m "Single regeneration scaling"
git push
./execute_analysis.sh 
rot
root
./execute_analysis.sh 
root
./execute_analysis.sh 
cd scripts/
root
cp /data/ssd/gamrat/python-kloe-analysis/scripts/results/* .
root
./execute_analysis.sh 
1
./execute_analysis.sh 
cd build/
rm -fr *
cd ..
cd build/
cmake -DCMAKE_BUILD_TYPE=Debug ..
make -j8
./execute_analysis.sh
make -j8
./execute_analysis.sh
git add .
git commit -m "Last version"
git push
git checkout -b feat/1d-corrections-regeneration
./execute_analysis.sh 
cd ../root_files/
ls
cd hal/
ls
cd ..
cd kitt/
ls
cd ALL_PHYS2_THREE_PI0_NoSmearing/
ls -hl
root 
cd ../../
ls
cd ..
ls
cd DBV-26/
ls
cd DK0/
ls
cd 20260516/
ls
cd ..
ls
cd all_phys
ls
cd ..
ls
cd all_phys
ls
cd 20260516/
ls
cd ../../
cd ../
ls
./lumi_

./lumi_batch.sh PROD2ROOT/DK0 26 data 4 lumi_file/lumi_per_ev_nb.log 
./lumi_batch.sh PROD2ROOT/DK0 26 DK0 4 lumi_file/lumi_per_ev_nb.log 
./lumi_batch.sh PROD2ROOT 26 DK0 4 lumi_file/lumi_per_ev_nb.log 
cat lumi_batch.sh 
vi lumi_batch.sh 
cp lumi_batch.sh lumi_batch_backup.sh 
vi lumi_batch
vi lumi_batch.sh 
./lumi_batch.sh PROD2ROOT 26 DK0 4 lumi_file/lumi_per_ev_nb.log 
cd DBV-26/
ls
mv DK0/ data/
ls
./lumi_batch.sh /data/ssd/gamrat/PROD2ROOT 26 DK0 4 lumi_file/lumi_per_ev_nb.log 
cd ..
./lumi_batch.sh /data/ssd/gamrat/PROD2ROOT 26 DK0 4 lumi_file/lumi_per_ev_nb.log 
cd DBV-26/
ls
rm data/
rm -r data/
mv DK0/ data/
cd ../root_files/hal/
ls
cd DATA_SEMILEPTONIC_NoSmearing/
ls
rm -fr *
ls
ls -hl
cd ../../kitt/
ls
cd DATA_THREE_PI0_NoSmearing/
ls
ls -hl
cd ..
cd ../hal/
ls
cd ALL_PHYS_SEMILEPTONIC_NoSmearing/
ls
root mk0_initial_analysis_all_phys_SEMILEPTONIC_NoSmearing_Semileptonic_99.root
git add .
git commit -m "Middle point of the three pi0 hypothesis"
git push
git checkout master
git checkout main 
ls
./execute_analysis.sh 
ls
cat run_parallel.sh 
vi parameters.txt 
cat para
cat parameters.txt 
./run_parallel.sh data 20260627 1 421 100
cd parallel_logs_data_1_421/
ls
cd 1
ls
cd 1
ls
cat stderr 
cd ..
ls
cd ..
cd 1/1
ls
cat seq 
cat stdout 
cd ../../../
./execute_analysis.sh /data/ssd/gamrat/DBV-26/data/20260627/job_v26_data_4_inv_pb_1.txt
./execute_analysis.sh 10 /data/ssd/gamrat/DBV-26/data/20260627/job_v26_data_4_inv_pb_1.txt
htop
./run_parallel.sh data 20260627 1 421 60
htop
ls
cd ../root_files/hal/
ls
cd DATA_SEMILEPTONIC_NoSmearing/
ls
ls -hl
htop
root dk0_initial_analysis_SEMILEPTONIC_NoSmearing_13.root
cd ../../kitt/DATA_THREE_PI0_NoSmearing/
ls
ls -hl
htop
cd ../root_files/
ls
cd hal/
ls
cd DATA_SEMILEPTONIC_NoSmearing/
ls
ls -hl
cd ../../kitt/
ls
cd DATA_THREE_PI0_NoSmearing/
ls
ls -hl
exit
cd ../root_files/kitt/
ls
cd DATA_THREE_PI0_NoSmearing/
ls
cd ..
cd ALL_PHYS_THREE_PI0_NoSmearing/
ls
cd ../
cd ..
cd kitt/DATA_THREE_PI0_NoSmearing/
ls
cat file_lumi_DATA_THREE_PI0_NoSmearing.log 
cd ..
ls
cd ../
ls
cd ../CNAF_Produced_Files/
ls
cd root_files/
ls
cd 2026-04-17/
ls
cd Signal/
ls
cd log/
ls
cat calculate
cat calculate.sh 
cat calculate_total.sh 
./calculate_total.sh 
cp calculate* ../../../../../KLOE/log/.
cp calculate* /data/4/users/gamrat/KLOE/log/.
cd ../../../../../root_files/kitt/
ls
cd ALL_PHYS_THREE_PI0_NoSmearing/
ls
cat file_lumi_ALL_PHYS_THREE_PI0_NoSmearing.log 
vi calculate_lumi.sh
chmod +x calculate_lumi.sh file_lumi_ALL_PHYS_THREE_PI0_NoSmearing.log 
./calculate_lumi.sh file_lumi_ALL_PHYS_THREE_PI0_NoSmearing.log 
cp calculate_lumi.sh ../.
cd ..
./calculate_lumi.sh ALL_PHYS_THREE_PI0_NoSmearing/file_lumi_ALL_PHYS_THREE_PI0_NoSmearing.log 
./calculate_lumi.sh ALL_PHYS2_THREE_PI0_NoSmearing/file_lumi_ALL_PHYS2_THREE_PI0_NoSmearing.log 
./calculate_lumi.sh ALL_PHYS3_THREE_PI0_NoSmearing/file_lumi_ALL_PHYS3_THREE_PI0_NoSmearing.log 
./calculate_lumi.sh DATA_THREE_PI0_NoSmearing/file_lumi_DATA_THREE_PI0_NoSmearing.log 
cd scripts/
root
root\
root
htop
root
cd scripts/
root
./execute_analysis.sh 
cd build/
rm -fr *
cd ..
./execute_analysis.sh 
q
./execute_analysis.sh 
.q
./execute_analysis.sh 
nohup ./execute_analysis.sh < parameters.txt > nohup.log &
cd ../root_files/hal/
ls
cd 2026-07-01/
ls
cd ALL_PHYS3_THREE_PI0_NoSmearing/
ls -hl
exit
cd scripts/
root
git add .
git checkout -b feature/three_pi0_additional_functions
git add .
git commit -m "Addition of good clusters check for six gammas"
git push
git push --set-upstream origin feature/three_pi0_additional_functions
./execute_analysis.sh 
cd ../root_files/hal/2026-07-01/
ls
cd ALL_PHYS3_THREE_PI0_NoSmearing/
ls -hl
root
cd ../root_files/hal/2026-07-02/ALL_PHYS3_THREE_PI0_NoSmearing/
ls
root
./execute_analysis.sh 
cd scripts/
root
./execute_analysis.sh 
git checkout main 
git checkout -- .
git checkout main
git pull
git checkout -b feature/semileptonic-analysis
git add .
git commit -m "Single vtx for semileptonic events"
git push
git add .
git commit -m "Added special configs"
git push
./run_parallel.sh 
./run_parallel.sh all_phys 20260516 1 206 15 analysis_config_MC_Semileptonic.json
cd ../root_files/hal/ALL_PHYS
cd ../root_files/hal/ALL_PHYS_SEMILEPTONIC_NoSmearing/
ls
cd ..
ls
rm -fr *
htop
killall KLSPM00
killall
killall -u gamrat
htop
./execute_analysis.sh 
killall
ps -e gamrat
kill
pkill -u gamrat
./execute_analysis.sh 
cd build/
rm -fr *
cd ..
./execute_analysis.sh 
./run_parallel.sh all_phys 20260516 1 206 15 analysis_config_MC_Semileptonic.json
./run_parallel.sh all_phys2 20260516 1 71 15 analysis_config_MC_Semileptonic.json
./run_parallel.sh all_phys3 20260516 1 76 15 analysis_config_MC_Semileptonic.json
./run_parallel.sh data 20260627 1 421 15 analysis_config_Data_Semileptonic.json
killall KLSPM00
./run_parallel.sh data 20260627 1 421 15 analysis_config_Data_Semileptonic.json
killall KLSPM00
killall parallel
killall KLSPM00
./run_parallel.sh all_phys 20260516 1 206 15 analysis_config_MC_Semileptonic.json
./run_parallel.sh all_phys2 20260516 1 71 15 analysis_config_MC_Semileptonic.json
./run_parallel.sh all_phys3 20260516 1 76 15 analysis_config_MC_Semileptonic.json
./run_parallel.sh data 20260627 1 421 15 analysis_config_Data_Semileptonic.json
killall KLSPM00
git add .
git commit -m "Elastic analysis config file"
git push
./run_parallel.sh all_phys 20260516 1 206 15 analysis_config_MC_Semileptonic.json
./run_parallel.sh all_phys2 20260516 1 71 15 analysis_config_MC_Semileptonic.json
killall KLSPM00
./run_parallel.sh all_phys3 20260516 1 76 15 analysis_config_MC_Semileptonic.json
./run_parallel.sh data 20260627 1 421 15 analysis_config_Data_Semileptonic.json
htop
./execute_analysis.sh 
kill $(pgrep -u gamrat | grep -v -e $$ -e $PPID)
pkill -u gamrat -t $(basename $(tty))
htop
pkill -u gamrat -t $(basename $(tty))
htop
cd /data/ssd/gamrat/root_files/hal/
;s
ls
cd ../../DBV-26/all_phys2/20260516/
ls
cd ../../data/
ls
cd 20260627/
ls -hl
cd ../../../
cd root_files/hal/
ls
pgrep -u gamrat parallel
pkill -9 -P $(pgrep -u gamrat parallel)
pkill -u gamrat parallel
ls
cd ALL_PHYS_THREE_PI0_NoSmearing/
ls -hl
cd ..
rm -fr *
ls
cd ALL_PHYS_SEMILEPTONIC_NoSmearing/
ls
ls -hl
cd ..
ls
cd ALL_PHYS2_SEMILEPTONIC_NoSmearing/
ls
ls -hl
htop
cd ../root_files/kitt/ALL_PHYS_THREE_PI0_NoSmearing/
ls -hl
cd ../root_files/kitt/ALL_PHYS2_THREE_PI0_NoSmearing/
cd ../ALL_PHYS2_THREE_PI0_NoSmearing/
ls
ls -hl
cd ../ALL_PHYS3_THREE_PI0_NoSmearing/
ls -hl
cd ../DATA_THREE_PI0_NoSmearing/
ls -hl
cd ../../hal/ALL_PHYS_SEMILEPTONIC_NoSmearing/
ls -hl
cd ../
ls
cd ../
ls
cd hal/
ls
exit
ls
cd ../root_files/hal/ALL_PHYS_SEMILEPTONIC_NoSmearing/
ls 0hl
ls -hl
cd ..
cd ../kitt/DATA_THREE_PI0_NoSmearing/
ls -hl
cd ../root_files/hal/ALL_PHYS_SEMILEPTONIC_NoSmearing/
ls
ls -hl
wc -l
l -hl wc -l
l -hl | wc -l
l -hl | wc
l -hl | wl
ls -hl | wc -l
ls -hl
root
ls
cat file_lumi_ALL_PHYS_SEMILEPTONIC_NoSmearing.log 
root
exit
cd ../python-kloe-analysis/
root
cd scripts/
root
htop
cd ../root_files/hal/DATA_SEMILEPTONIC_NoSmearing/
ls -hl
ls
cd ../ALL_PHYS
cd ../ALL_PHYS_SEMILEPTONIC_NoSmearing/
ls -hl
cd ../ALL_PHYS2_SEMILEPTONIC_NoSmearing/
ls -hl
cd ../
ls
mkdir ROOT-in-docker
docker ps
cd KLe
cd KLOE/
cd ../python-kloe-analysis/
cd scripts/
root
htop
cd ../root_files/
ls
cd hal/
ls
cd../ki
cd ../kitt
ls
cd ALL_PHYS_THREE_PI0_NoSmearing/
ls
root
ls
root
./execute_analysis.sh 
cd build/
rm -fr *
cd ..
./execute_analysis.sh 
.q
./execute_analysis.sh 
cd so
cd scripts/
root
.q
root
.q
root
git checkout -b feature/other-flags
git add .
git commit -m "Flag check for other channel"
git push
git push --set-upstream origin feature/other-flags
root
cd img/other_flags/
ls
cd ..
root
cd results/detailed_other_flags/
ls
root flags.root 
cd ../../
root
cd results/detailed_other_flags/
ls
root flags.root 
cd ../../
htop
cd ../root_files/
ls
cd hal
ls
cd 2026-07-10/
ls -hl
cd ALL_PHYS3_THREE_PI0_NoSmearing/
ls
root 
cd build/
rm -fr *
cd ..
./execute_analysis.sh 
cd ../
cd root_files/
ls -hl
cd hal/
ls
rm -fr *
cd ../kitt/
rm -fr *
ls
cd ..
htop
cd ../h
cd hal/
ls
killall KLSPM00
htop
killall KLSPM00
ls
rm -fr *
cd ../
ls
cd kitt/
ls -hl
cd ../../
ls
cd KLOE/
git add .
git commit -m "Other flag + isr"
git push
git checkout main 
git pull
./execute_analysis.sh 
cd log/
ls
rm -fr all_phys*
rm -fr data*
ls
cd ../
rm -fr nohup*
rm -fr parallel_logs_*
./run_parallel.sh all_phys 20260516 1 206 15 "analysis_config_MC_Semileptonic.json"
./run_parallel.sh all_phys2 20260516 1 71 15 "analysis_config_MC_Semileptonic.json"
./run_parallel.sh all_phys3 20260516 1 76 25 "analysis_config_MC
./run_parallel.sh all_phys3 20260516 1 76 15 "analysis_config_MC_Semileptonic.json"
./run_parallel.sh data 20260627 1 421 15 "analysis_config_Data_Semileptonic.json"
ls
cd ../root_files/hal/
ls
cd ALL_PHYS
cd ALL_PHYS_SEMILEPTONIC_NoSmearing/
ls -hl
cd scripts/
root
cd scripts/
root
cd scripts/
root
cd ../../root_files/kitt/DATA_THREE_PI0_NoSmearing/
;s
ls
root
cd 
cd /data/ssd/gamrat/python-kloe-analysis/s
cd /data/ssd/gamrat/python-kloe-analysis/scripts/
root
htop
cd scripts/
root
cd ../../root_files/kitt/ALL_PHYS_THREE_PI0_NoSmearing/
root
cd scripts/
root
