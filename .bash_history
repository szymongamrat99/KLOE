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
