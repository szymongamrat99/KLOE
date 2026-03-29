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
