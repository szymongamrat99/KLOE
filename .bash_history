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
