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
