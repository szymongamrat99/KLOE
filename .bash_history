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
cd scripts/
root
cd scripts/
root
./execute_analysis.sh 
cd build/
git checkout -b feature/three_pi0_additional_functions 
git checkout feature/three_pi0_additional_functions 
git add .
git commit -m "Bug fixed"
cd ..
git add .
git commit -m "Bug fixed"
git checkout feature/three_pi0_additional_functions 
git push
git merge master
git merge main
git push
./execute_analysis.sh 
cd build/
rm -fr *
cd ..
./execute_analysis.sh 
nohup ./execute_analysis.sh < par.txt > nohup.log &
cd ../root_files/hal/
ls
cd 2026-07-14/
ls
cd ALL_PHYS3_THREE_PI0_NoSmearing/
ls
root 
root
cd ../root_files/hal/2026-07-14/
ls
cd ALL_PHYS3_THREE_PI0_NoSmearing/
root
exit
./run_parallel.sh 
./execute_analysis.sh 
git checkout feature/three_pi0_additional_functions 
git add .
git commit -m "Modified charged"
git push
git checkout feature/three_pi0_additional_functions 
ls
git checkout feature/control-samples-semileptonic-3pi0 
git checkout feature/three_pi0_additional_functions 
./execute_analysis.sh 
./e
./execute_analysis.sh 
cd scripts/results/
ls
cd control_sample_corr_factors/
pwd
cd scripts/
./analyses_corr_fact.sh 
./analyses_corr_fact.sh kitt THREE_PI0 20 0.0
pkill -9 root.exe 
./analyses_corr_fact.sh kitt THREE_PI0 20 0.0
./analyses_corr_fact.sh kitt THREE_PI0 450 0.0
./analyses_corr_fact.sh hal SEMILEPTONIC 450
ps -e root
pgrep root
ps -e "root.exe"
ps -e | grep root.exe
kill 15244
./analyses_corr_fact.sh hal SEMILEPTONIC 450
./execute_analysis.sh 
git add .
git commit -m "Corrections"
git push
htop
cd ../root_files/
ls
cd hal/
ls
ls -hl
cd DATA_SEMILEPTONIC_NoSmearing/
ls -hl
cd ../
cd ALL_PHYS
cd ALL_PHYS_SEMILEPTONIC_NoSmearing/
root
cd ../
ls
cd ../
ls
cd CNAF_Produced_Files/
ls
cd root_files/
ls
rsync -avz tier1-cnaf:/home/g/gamrat/root_files/2026-07-26 .
ssh tier1-cnaf 
rsync -avz tier1-cnaf:/home/g/gamrat/root_files/2026-07-19 .
cd scripts/
ls
./analyses_corr_fact.sh hal SEMILEPTONIC 20 1 0.0
./analyses_corr_fact.sh hal SEMILEPTONIC 20 1
pkill -9 root.exe
./analyses_corr_fact.sh hal SEMILEPTONIC 20 1
./analyses_corr_fact.sh hal SEMILEPTONIC 20 5
./analyses_corr_fact.sh hal SEMILEPTONIC 20 2
./analyses_corr_fact.sh hal SEMILEPTONIC 20 3
./analyses_corr_fact.sh hal SEMILEPTONIC 20 4
./analyses_corr_fact.sh hal SEMILEPTONIC 20 5
pkill -9 root.exe
./analyses_corr_fact.sh hal SEMILEPTONIC 20 1
./analyses_corr_fact.sh hal SEMILEPTONIC 20 2
./analyses_corr_fact.sh hal SEMILEPTONIC 20 3
./analyses_corr_fact.sh hal SEMILEPTONIC 20 4
./analyses_corr_fact.sh hal SEMILEPTONIC 20 5
./analyses_corr_fact.sh hal SEMILEPTONIC 20 1
./analyses_corr_fact.sh hal SEMILEPTONIC 20 2
kill -9 root.exe
pkill -9 root.exe
./analyses_corr_fact.sh hal SEMILEPTONIC 20 1
./analyses_corr_fact.sh hal SEMILEPTONIC 20 2
./analyses_corr_fact.sh hal SEMILEPTONIC 20 3
./analyses_corr_fact.sh hal SEMILEPTONIC 20 4
./analyses_corr_fact.sh hal SEMILEPTONIC 20 4 true
./analyses_corr_fact.sh hal SEMILEPTONIC 450 1
./analyses_corr_fact.sh hal SEMILEPTONIC 450 2
./analyses_corr_fact.sh hal SEMILEPTONIC 450 3
./analyses_corr_fact.sh hal SEMILEPTONIC 450 4 true
./analyses_corr_fact.sh kitt THREE_PI0 450 1 0.0
./analyses_corr_fact.sh kitt THREE_PI0 450 2 0.0
./analyses_corr_fact.sh kitt THREE_PI0 450 3 0.0 true
./analyses_corr_fact.sh hal SEMILEPTONIC 450 4 true
cd ../root_files/kitt/ALL_PHYS_SIGNAL_NoSmearing/
ls -hl
cd ../../../KLOE/
./execute_analysis.sh 
cd scripts/
ls
./analyses_corr_fact.sh hal SEMILEPTONIC 450 1
./analyses_corr_fact.sh hal SEMILEPTONIC 450 2
./analyses_corr_fact.sh hal SEMILEPTONIC 450 3
./analyses_corr_fact.sh kitt THREE_PI0 450 2 0.0
./analyses_corr_fact.sh kitt THREE_PI0 450 1 0.0
./analyses_corr_fact.sh kitt THREE_PI0 450 3 0.0 true
./analyses_corr_fact.sh hal SEMILEPTONIC 450 4 true
cd ../CNAF_Produced_Files/
ls
cd root_files/
ls
cd 2026-07-19/
ls
cd signal/
ls
cd ALL_PHYS
cd ALL_PHYS_SIGNAL_NoSmearing/
ls
scp tier1-cnaf:/home/g/gamrat/root_files/2026-07-19/signal/ALL_PHYS_SIGNAL_NoSmearing/mk0_initial_analysis_all_phys_SIGNAL_NoSmearing_173.root .
cd
ext
exit
cd ../
cd CNAF_Produced_Files/
ls
cd root_files/
ls
cd 2026-0
cd 2026-04-17/
ls
cd Signal/
ls
mv ALL_PHYS_SIGNAL_NoSmearing/ ALL_PHYS_SIGNAL_NoSmearing_backup
mv ALL_PHYS2_SIGNAL_NoSmearing/ ALL_PHYS2_SIGNAL_NoSmearing_backup
mv ALL_PHYS3_SIGNAL_NoSmearing/ ALL_PHYS3_SIGNAL_NoSmearing_backup
ls
cp -r ../../2026-07-19/signal/ALL_PHYS_SIGNAL_NoSmearing .
cp -r ../../2026-07-19/signal/ALL_PHYS2_SIGNAL_NoSmearing .
cp -r ../../2026-07-19/signal/ALL_PHYS3_SIGNAL_NoSmearing .
ls -hl
cd ..
cd ../../
cd KLOE/
ls
./execute_analysis.sh 
cd log/
ls
./calculate_total.sh 
./calculate_total.sh all_phys2
./calculate_total.sh all_phys3
./calculate_total.sh all_phys
./calculate_total.sh
./calculate.sh all_phys
./calculate.sh all_phys2
./calculate.sh all_phys3
./calculate.sh all_phys3 .
./calculate.sh all_phys .
./calculate.sh all_phys2 .
./calculate.sh all_phys3 .
./execute_analysis.sh 
cd ..
cd CNAF_Produced_Files/root_files/2026-07-19/
ls
ls ../../../root_files/kitt/
cd ../../../
cd root_files/
mkdir cnaf
cd cnaf/
ls
cd ..
ls
cd hal/
ls
cd ../
cd cnaf/
ls
cp -fr ../../CNAF_Produced_Files/root_files/2026-07-19/4pi/* .
cd ..
ls
cd cnaf
ls
cd DATA_FOUR_PI_NoSmearing/
ls -hl
exit
cd scripts/
root
./analyses_corr_fact.sh 
./analyses_corr_fact.sh hal SEMILEPTONIC 1 4
./analyses_corr_fact.sh hal SEMILEPTONIC 10 4
./analyses_corr_fact.sh hal SEMILEPTONIC 420 4 true
./analyses_corr_fact.sh hal SEMILEPTONIC 10 4 true
./analyses_corr_fact.sh kitt THREE_PI0 10 4 0.0 true
./analyses_corr_fact.sh kitt THREE_PI0 450 3 0.0 true
./analyses_corr_fact.sh hal SEMILEPTONIC 450 4 true
./analyses_corr_fact.sh cnaf DOUBLE_PI 450 3 0.0 true
./analyses_corr_fact.sh cnaf FOUR_PI 450 3 0.0 true
./analyses_corr_fact.sh cnaf FOUR_PI 450 1 0.0 true
pkill -9 root.exe
./analyses_corr_fact.sh cnaf FOUR_PI 450 1 0.0 true
./analyses_corr_fact.sh cnaf FOUR_PI 450 2 0.0 true
./analyses_corr_fact.sh cnaf FOUR_PI 450 3 0.0 true
./analyses_corr_fact.sh hal SEMILEPTONIC 450 4 true
./analyses_corr_fact.sh kitt THREE_PI0 450 3 0.0 true
git add .
git commit -m "Inclusion of additional things for Regeneration weight calculation"
./execute_analysis.sh 
git add .
git commit -m "Inclusion of additional things for Regeneration weight calculation"
git push
git checkout -b feature/regeneration-weight-with-fit
./execute_analysis.sh 
git checkout feature/three_pi0_additional_functions 
git checkout -- .
./execute_analysis.sh 
cd ../CNAF_Produced_Files/
ls
cd root_files/2026-04-17/
ls
cd Signal/
ls
cd ALL_PHYS2_SIGNAL_NoSmearing
ls
cd ../
cd ALL_PHYS3_SIGNAL_NoSmearing
ls
cd ..
ls
cd DATA_SIGNAL_NoSmearing/
ls
cd ..
ls
cd ..
cd ../../../
cd KLOE/
./execute_analysis.sh 
cd scripts/
root
git add .
git commit -m "Preparation of signal weights"
git push
./execute_analysis.sh 
git add .
git commit -m "Determination of corr factors for signal - approach 1"
git push
cd scripts/
./analyses_corr_fact.sh cnaf FOUR_PI 450 3 0.0 true
./analyses_corr_fact.sh kitt THREE_PI0 450 3 0.0 true
./analyses_corr_fact.sh hal SEMILEPTONIC 450 4 true
root
cd ../CNAF_Produced_Files/
ls
cd root_files/
ls
cd 2026-04-17/
ls
cd SI
cd Signal/
ls
cd log/
ls
./calculate.sh 
cd lo
cd ..
ls
cd DATA_SIGNAL_NoSmearing/
ls
cd ../
ls
cd DATA_SIGNAL_NoSmearing/
cat file_lumi_DATA_SIGNAL_NoSmearing.log 
cd ..
vi lumi.sh
chmod +x lumi.sh 
./lumi.sh DATA_SIGNAL_NoSmearing/file_lumi_DATA_SIGNAL_NoSmearing.log 
./lumi.sh ALL_PHYS_SIGNAL_NoSmearing/file_lumi_ALL_PHYS_SIGNAL_NoSmearing.log 
./lumi.sh ALL_PHYS2_SIGNAL_NoSmearing/file_lumi_ALL_PHYS2_SIGNAL_NoSmearing.log 
./lumi.sh ALL_PHYS3_SIGNAL_NoSmearing/file_lumi_ALL_PHYS3_SIGNAL_NoSmearing.log 
exit
cd ..
cd CNAF_Produced_Files/root_files/2026-07-19/
ls
rsync -avz tier1-cnaf:/home/g/gamrat/root_files/2026-07-19/log .
ls
cd log/
ls
cd all_phys
ls
cd 1
ls
cat cut.analysis.log 
cd ..
cp ../../../KLOE/log/calculate* .
ls
./calculate.sh log/all_phys
mv calculate* log/.
cd log/
ls
./calculate.sh all_phys
vi calculate
vi calculate.sh 
./calculate.sh all_phys .
cat calculate.sh 
cd all_phys/1
ls
cat cut.analysis.log 
cd ../../
ls
vi calculate.sh 
touch calculate_new.sh
vi calculate_new.sh 
chmod +x calculate_new.sh 
./calculate_new.sh 
./calculate_new.sh all_phys .
./calculate_new.sh all_phys2 .
./calculate_new.sh all_phys3 .
ls
cd ..
ls
cd signal/
ls
cd ..
cd 4pi/
ls
cd ..
cd ../2026-04-17/
ls
cd si
cd Signal/
ls
./lumi.sh ALL_PHYS_SIGNAL_NoSmearing/file_lumi_ALL_PHYS_SIGNAL_NoSmearing.log 
./lumi.sh ALL_PHYS2_SIGNAL_NoSmearing/file_lumi_ALL_PHYS2_SIGNAL_NoSmearing.log 
./lumi.sh ALL_PHYS3_SIGNAL_NoSmearing/file_lumi_ALL_PHYS3_SIGNAL_NoSmearing.log 
./lumi.sh DATA_SIGNAL_NoSmearing/file_lumi_DATA_SIGNAL_NoSmearing.log 
htop
cd ../../2026-07-19/
ls
cd log/
ls
./calculate_new.sh all_phys .
./calculate_new.sh all_phys2 .
./calculate_new.sh all_phys3 .
ls
cd all_phys
ls
cd 1
ls
cat error.log 
cd scripts/
root
./analyses_corr_fact.sh cnaf FOUR_PI 450 3 0.0 true
./analyses_corr_fact.sh hal SEMILEPTONIC 450 4 true
./analyses_corr_fact.sh kitt THREE_PI0 450 3 0.0 true
root
./analyses_corr_fact.sh hal SEMILEPTONIC 450 4 true
./analyses_corr_fact.sh kitt THREE_PI0 450 3 0.0 true
root
./execute_analysis.sh 
./e
./execute_analysis.sh 
git add .
git commit -m "New version of the corr factor check"
git push
root
ls
cd scripts/
ls
./analyses_corr_fact.sh 
./analyses_corr_fact.sh hal SEMILEPTONIC 450 4 true
./analyses_corr_fact.sh kitt THREE_PI0 450 3 0.0 true
./execute_analysis.sh 
cd ../root_files/
ls
cd hal/
ls
cd ALL_PHYS_SEMILEPTONIC_NoSmearing/
ls
root mk0_initial_analysis_all_phys_SEMILEPTONIC_NoSmearing_99.root
cd ..
ls
mv ALL_PHYS_SEMILEPTONIC_NoSmearing ALL_PHYS_SEMILEPTONIC_NoSmearing_backup
ls
cd ..
ls
cd hal/
ls
cd ALL_PHYS_SEMILEPTONIC_NoSmearing
ls -hl
ls
ls -hl
root
cd log/
ls
mv all_phys all_phys_bkp
mv all_phys2 all_phys2_bkp
mv all_phys3 all_phys3_bkp
mv data data_bkp
cd ..
./run_parallel.sh all_phys 20260516 1 5 5 "analysis_config_MC_Semileptonic.json"
htop
ls
htop
./run_parallel.sh all_phys 20260516 6 20 15 "analysis_config_MC_Semileptonic.json"
cd ../root_files/
ls
cd hal/
ls
cd ALL_PHYS_SEMILEPTONIC_NoSmearing
ls
ls -hl
root
ssh kloerec
exit
ls
./execute_analysis.sh 5 job_v26_all_phys3_1_inv_pb_1.txt 
./execute_analysis.sh 5 ../../job_v26_all_phys3_1_inv_pb_1.txt 
./execute_analysis.sh 5 ../job_v26_all_phys3_1_inv_pb_1.txt 
root
ls
scp tier1-cnaf:/home/g/gamrat/geanfi_test/sample.root .
ls
root sample.root 
ls
./execute_analysis.sh 
cd ../root_files/hal/
ls
cd ALL_PHYS3_SIGNAL_NoSmearing/
ls
ls -hl
root mk0_initial_analysis_all_phys3_SIGNAL_NoSmearing_1.root 
exit
cd ../root_files/
ls
cd ../DBV-26/all_phys
ls
cd 20260516/
ls
cd ../../
cd ../
ls
cd PROD2ROOT/
ls
cd MC/MK0/all_phys3/
ls
root prod2root_mk0_all_phys3_41228_v2.root
cd ../../../../
ls
cd root_files/
ls
cd hal/
ls
cd ALL_PHYS_SEMILEPTONIC_NoSmearing
ls
root mk0_initial_analysis_all_phys_SEMILEPTONIC_NoSmearing_9.root 
cd ../../../PROD2ROOT/
ls
cd MC/MK0/all_phys3/
root prod2root_mk0_all_phys3_41129_v2.root
scp tier1-cnaf:~/geanfi_test/sample.root .
ls
./execute_analysis.sh 
./execute_analysis.sh 5 ../job_v26_all_phys3_1_inv_pb_1.txt
ls
root
cd ../
ls
cd KLOE/
cd ../root_files/hal/ALL_PHYS3_SIGNAL_NoSmearing/
ls
cat file_lumi_ALL_PHYS3_SIGNAL_NoSmearing.log 
ls
cat input_luminosity_ALL_PHYS3_SIGNAL_NoSmearing.log 
root mk0_initial_analysis_all_phys3_SIGNAL_NoSmearing_1.root 
exit
scp tier1-cnaf:~/geanfi_test/sample.root .
./execute_analysis.sh 5 ../job_v26_all_phys3_1_inv_pb_1.txt
root sample.root 
exit
cd scripts/semileptonic_cuts_optimal/
root
cd scripts/semileptonic_cuts_optimal/
ls
root
ls
cd ../root_files/hal/
ls
htop
./execute_analysis.sh 
cd scripts/semileptonic_cuts_optimal/
root
cd ../
cd results/
ls -hl
cd control_sample_corr_factors/
ls -hl
docker
docker ps
exit
docker ps
exit
ls
cd ..
ls
cd ksoft-container/
ls
docker compose up -d --build
cat docker-compose.yaml 
exit
ssh-add -l
cd ../ksoft-container/
ls
docker compose up -d --build
vi docker-compose.yaml 
vi Dockerfile 
docker compose up -d --build
docker ps
docker exec -it 83ca398ed6a6 bash
cd ../
cd ksoft-container/
ls
vi docker-compose.yaml 
vi Dockerfile 
vi docker-compose.yaml 
cd ..
mkdir geanfi-workspace
cd ksoft-container/
vi docker-compose.yaml 
docker exec -it 83ca398ed6a6 bash
docker compose up -d --build
ls
cd ..
cd geanfi-workspace/
ls
docker exec -it 83ca398ed6a6 bash
docker ps
docker exec -it 1a02f3c2501f bash
ls
rm geanfi-workspace/*
ls
rm *
ls
rm -fr geanfi-workspace/
ls
docker exec -it 1a02f3c2501f bash
exit
cd ..
scp fibm15:/kbackup/db2.tar .
scp fibm15:/kbackup/DB2/db2.tar .
ls
exit
htop
wget https://public.dhe.ibm.com/ibmdl/export/pub/software/data/db2/drivers/odbc_cli/linuxx64_odbc_cli.tar.gz 
ls
rm linuxx64_odbc_cli.tar.gz 
cd ..
rm linuxx64_odbc_cli.tar.gz 
wget https://public.dhe.ibm.com/ibmdl/export/pub/software/data/db2/drivers/odbc_cli/linuxx64_odbc_cli.tar.gz 
ls
tar -zxvf linuxx64_odbc_cli.tar.gz
ls
docker ps
docker exec -it 29511853c50a bash
docker exec -ti ksoft_container bash
docker exec -i db2server su - db2inst1 -c "db2 connect to kloemaps user dbread using HASLO > /dev/null && db2 \"select layer, prob_k1, prob_k2 from cndrng.mc_kch_noise where (Start_of_validity<=26759) and (26759<End_of_validity) order by layer\""
docker exec -i db2server su - dbread -c "db2 connect to kloemaps user dbread using HASLO > /dev/null && db2 \"select layer, prob_k1, prob_k2 from cndrng.mc_kch_noise where (Start_of_validity<=26759) and (26759<End_of_validity) order by layer\""
docker exec -i db2server su - db2inst2 -c "db2 connect to kloemaps user dbread using HASLO > /dev/null && db2 \"select layer, prob_k1, prob_k2 from cndrng.mc_kch_noise where (Start_of_validity<=26759) and (26759<End_of_validity) order by layer\""
ls
cp clidriver/ ksoft-container/.
cp -r clidriver ksoft-container/.
cd ksoft-container/
ls
rm last.kumac paw.metafile 
ls
cd ../
ls
cd db2-container/
ls
cd database/
ls
cd DB2/
ls
cat LOAD.out 
llq
ls
cd..
cd ..
ls
docker login
cd db2-container/
ls
cat docker-compose.yaml 
docker login
ls
rm -fr docker-compose.yaml 
cd db2_data/
ls
cd ..
docker login 
docker pull icr.io/db2_community/db2
ls
dbdir=`pwd`
#fill .env_list according to IBM docs (DB2INSTANCE=db2inst2,leave  DBNAME= empty, set your DB2INST2_PASSWORD=[password])
echo -e "LICENSE=accept\nDB2INSTANCE=db2inst2\nDB2INST1_PASSWORD=megapass\nDBNAME=\n\
BLU=false\nENABLE_ORACLE_COMPATIBILITY=false\nUPDATEAVAIL=NO\nTO_CREATE_SAMPLEDB=false\n\
REPODB=false\nIS_OSXFS=false\nPERSISTENT_HOME=true\nHADR_ENABLED=false\n\
ls
echo -e "LICENSE=accept\nDB2INSTANCE=db2inst2\nDB2INST1_PASSWORD=megapass\nDBNAME=\n\
BLU=false\nENABLE_ORACLE_COMPATIBILITY=false\nUPDATEAVAIL=NO\nTO_CREATE_SAMPLEDB=false\n\
REPODB=false\nIS_OSXFS=false\nPERSISTENT_HOME=true\nHADR_ENABLED=false\n\
vi fill_env.sh
bash fill_env.sh 
ls -a
cat .env_list 
ls
rm fill_env.sh 
rm -fr db2_data/
ls
mkdir database
docker run -h db2server --name db2server --restart=always --detach --privileged=true -p 50000:50000 --env-file .env_list -v ${dbdir}/database:/database icr.io/db2_community/db2
docker ps
cd database/
ls
cd ../../
mkdir DB2
mv db2.tar DB2/.
cd DB2
ls
tar -xvf nazwa_pliku.tar
tar -xvf db2.tar
ls
mv db2.tar ../.
cd ..
cp -a DB2 ${dbdir}/database/
cd ${dbdir}/database/DB2
ls
rm *.msg
ls
kloeddl=kloemaps_ddl_fix.sql
cp kloemaps_ddl.sql $kloeddl
sed -i 's|/DB2data/db2inst2/db2inst2/NODE0000/SQL00001/SQLT000DMS.0|SQLT000DMS.0|g' $kloeddl
sed -i 's|/DB2data/db2inst2/db2inst2/NODE0000/SQL00001/SYSTOOLSPACE|SYSTOOLSPACE|g' $kloeddl
sed -i 's|/DB2data/db2inst2/db2inst2/NODE0000/SQL00001/SYSTOOLSTMPSPACE|SYSTOOLSTMPSPACE|g' $kloeddl
sed -i 's|/DB2data/db2inst2/db2inst2|/database/data/db2inst2|g' $kloeddl
sed -i -e "s|SQLSTATE 'KLTRGM1'|SQLSTATE 'KLTM1'|g"  -e "s|SQLSTATE 'KLTRGM3'|SQLSTATE 'KLTM3'|g" $kloeddl
sed -i "s|^create trigger DESCRIPT.FSP_FARM_ins|--#SET TERMINATOR @\ncreate trigger DESCRIPT.FSP_FARM_ins|g" $kloeddl
sed -i "s|^create trigger DESCRIPT.FSP_MC_ins|--#SET TERMINATOR @\ncreate trigger DESCRIPT.FSP_MC_ins|g" $kloeddl
sed -i "s|^create trigger DESCRIPT.FSP_OFFLINE_ins|--#SET TERMINATOR @\ncreate trigger DESCRIPT.FSP_OFFLINE_ins|g" $kloeddl
sed -i "s|^END;|END@\n--#SET TERMINATOR ;|g" $kloeddl
sed -i "s|LOW2KEY='charged rad DST'|LOW2KEY='Bhabha Monte Carlo DST'|g" $kloeddl
sed -i "s|LOW2KEY=X'636F736D69635F6D75'|LOW2KEY='Kl_all'|g" $kloeddl
sed -i "s|LOW2KEY=X'652B652D'|LOW2KEY='K+'|g" $kloeddl
sed -i "s|LOW2KEY=X'636F736D69635F6D756F6E73'|LOW2KEY='K+'|g" $kloeddl
sed -i 's|RUNSTATS ON TABLE "SYSTOOLS"."HMON_ATM_INFO"|--RUNSTATS ON TABLE "SYSTOOLS"."HMON_ATM_INFO"|g' $kloeddl
sed -i 's|RUNSTATS ON TABLE "SYSTOOLS"."HMON_COLLECTION"|--RUNSTATS ON TABLE "SYSTOOLS"."HMON_COLLECTION"|g' $kloeddl
sed -i '/ALTER TABLE "DB2INST2"."PROVABINNL"/,+3 s/^/--/' $kloeddl
sed -i '/ALTER TABLE "HEPDB   "."PROVABINNL"/,+3 s/^/--/' $kloeddl
sed -i '/CREATE UNIQUE INDEX "DB2INST2"."INDID"/,+3 s/;/;\nALTER TABLE "DB2INST2"."PROVABINNL" ADD PRIMARY KEY ("ID");/' $kloeddl
sed -i '/CREATE UNIQUE INDEX "DB2INST2"."INDBIN"/,+3 s/;/;\nALTER TABLE "HEPDB   "."PROVABINNL" ADD PRIMARY KEY ("ID");/' $kloeddl
sed -i '/PUBBLIC/s/^/--/' $kloeddl
(cat $kloeddl | sed -n '/This CLP file was created/,/Statements for foreign keys/p'; echo -e "COMMIT WORK;\nCONNECT RESET;\nTERMINATE;") > ${kloeddl}_tables
(echo "CONNECT TO KLOEMAPS;"; cat $kloeddl | sed -n '/Statements for foreign keys/,/TERMINATE/p' ) > ${kloeddl}_foreigns
cp -a db2move.lst db2move.lst_bckp
sed -i '/"SYSTOOLS"."HMON/d' db2move.lst
ls -hl
declare -A tbl=(["HEPDB.BANKS"]="REF" ["HEPDB.BANKS1"]="REF" ["HEPDB.BANKSA"]="REF" ["LOGGER.FILE_CNAF"]="FID" ["LOGGER.DATAREC_LOG"]="ID"  ["LOGGER.MC_LOG"]="ID"  ["LOGGER.RAW_LOG"]="ID")
(echo "connect to kloemaps;"; for t in ${!tbl[@]}; do ( echo "ALTER TABLE ${t} ALTER COLUMN ${tbl[$t]} SET GENERATED BY DEFAULT;"; ) done; echo -e "CONNECT RESET;\nTERMINATE;") > generatedcols.sql
(echo "connect to kloemaps;"; for t in ${!tbl[@]}; do ( echo "ALTER TABLE ${t} ALTER COLUMN ${tbl[$t]} SET GENERATED ALWAYS;"; ) done; echo -e "CONNECT RESET;\nTERMINATE;") > generatedcols_revert.sql
echo -e 'connect to kloemaps;\nRUNSTATS ON TABLE "DESCRIPT"."MCSTREAM_DESCRIPT";\nRUNSTATS ON TABLE "DESCRIPT"."MCSTREAM_GROUP";\nRUNSTATS ON TABLE "DESCRIPT"."STREAM_OFFLINE";\nCONNECT RESET;\nTERMINATE;' > updatestat.sql
docker exec -ti db2server bash -c "su - db2inst2"
ls -l
cd ..
ls -l
cd DB2/
ls
docker exec -ti db2server bash
ls
cd ..
chmod -R 666 DB2/
cd DB2
ls
ls -l
docker exec -ti db2server bash
docker exec -ti db2server bash -c "su - db2inst2"
guid
id -u
id -g
cd DB2
ls -l
chown 777 DB2
chmod 777 DB2
ls
cd DB2/
ls
ls -l
docker exec -ti db2server bash -c "su - db2inst2"
docker exec -ti db2server bash -c "useradd -r -s /sbin/nologin --no-create-home dbread; passwd dbread" 
docker exec -ti db2server bash -c "su - db2inst2"
cd ksoft-container/
docker compose up -d --build
docker ps
docker exec -it a5dc57b8f0ad bash
tar -zxf v12.1.5_linuxx64_client.tar.gz 
cd client/
./db2setup -f sysreq
cd ..
rm -fr client/
docker compose up -d --build
docker logs
docker ps
docker logs 29511853c50a
docker compose up -d --build
tar -zxf v12.1.5_linuxx64_client.tar.gz 
cd client/
./db2_install --help
cd ..
rm -fr client/
docker compose up -d --build
tar -zxf v12.1.5_linuxx64_client.tar.gz 
./db2_install -b /root/sqllib -p CLIENT -f sysreq
cd client/
./db2_install -b ~/sqllib -p CLIENT -f sysreq
./db2setup -f sysreq
./db2setup
cd ..
rm -fr client/
docker compose up -d --build
du -h
df -h
docker compose up -d --build
tar -zxf v12.1.5_linuxx64_client.tar.gz 
docker compose up -d --build
docker system prune -a --volumes
docker compose up -d --build
docker system prune -a --volumes
docker compose up -d --build
docker build --progress=plain -t ksoft-cernlib:latest .
docker compose up --build --progress=plain
BUILDKIT_PROGRESS=plain docker compose up --build
docker compose up -d --build
du -h
df -h
docker compose up -d --build
docker images
docker builder prune -f
docker compose up -d --build

docker images
cd scripts/
ls
root
./analyses_corr_fact.sh 
nohup ./analyses_corr_fact.sh hal semileptonic 500 4 true > logs.semi.4 &
nohup ./analyses_corr_fact.sh hal SEMILEPTONIC 500 4 true > logs.semi.4 &
nohup ./analyses_corr_fact.sh kitt THREE_PI0 500 4 true > logs.three.3 &
cd ../../
cd root_files/hal/
ls
mv ALL_PHYS_SEMILEPTONIC_NoSmearing ALL_PHYS_SEMILEPTONIC_NoSmearing_backup2
mv ALL_PHYS_SEMILEPTONIC_NoSmearing_backup ALL_PHYS_SEMILEPTONIC_NoSmearing
ls
cd ../../python-kloe-analysis/scripts/
nohup ./analyses_corr_fact.sh hal SEMILEPTONIC 500 4 true > logs.semi.4 &
nohup ./analyses_corr_fact.sh kitt THREE_PI0 500 3 true > logs.three.3 &
nohup ./analyses_corr_fact.sh kitt THREE_PI0 1 3 true > logs.three.3 &
cd scripts/
nohup ./analyses_corr_fact.sh kitt THREE_PI0 1 3 true > logs.three.3 &
nohup ./analyses_corr_fact.sh kitt THREE_PI0 1 3 0.0 true > logs.three.3 &
nohup ./analyses_corr_fact.sh kitt THREE_PI0 500 3 0.0 true > logs.three.3 &
nohup ./analyses_corr_fact.sh hal SEMILEPTONIC 500 5 true > logs.semi.4 &
pkill -9 root.exe
nohup ./analyses_corr_fact.sh hal SEMILEPTONIC 500 5 > logs.semi.4 &
nohup ./analyses_corr_fact.sh hal SEMILEPTONIC 500 6 true > logs.semi.4 &
nohup ./analyses_corr_fact.sh kitt THREE_PI0 500 4 0.0 > logs.three.3 &
nohup ./analyses_corr_fact.sh kitt THREE_PI0 500 5 0.0 > logs.three.3 &
nohup ./analyses_corr_fact.sh kitt THREE_PI0 500 6 0.0 true > logs.three.3 &
nohup ./analyses_corr_fact.sh hal SEMILEPTONIC 500 5 true > logs.semi.4 &
pkill -9 root.exe
nohup ./analyses_corr_fact.sh hal SEMILEPTONIC 500 5 > logs.semi.4 &
nohup ./analyses_corr_fact.sh hal SEMILEPTONIC 500 6 true > logs.semi.4 &
nohup ./analyses_corr_fact.sh kitt THREE_PI0 500 4 0.0 > logs.three.3 &
nohup ./analyses_corr_fact.sh kitt THREE_PI0 500 5 0.0 > logs.three.3 &
nohup ./analyses_corr_fact.sh kitt THREE_PI0 500 6 0.0 true > logs.three.3 &
./execute_analysis.sh 
.q
./execute_analysis.sh 
cd scripts/results/
ls
cd control_sample_corr_factors/
ls
root control_sample_corr_factors_backup.root 
pip install uproot
python venv .venv
python -m pyenv .venv
python -m venv .venv
cd ..
python3 plotter.py 
python plotter.py 
./execute_analysis.sh
ls
cd ../root_files/hal/
ls
cd 2026-07-20/
ls
cd ..
ls
exit
