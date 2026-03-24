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
