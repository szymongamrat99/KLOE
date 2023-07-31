#include <iostream>
#include <TFile.h>
#include <TTree.h>

gROOT->ProcessLine(".x load.C");

void time_diff_cm()
{
    TFile file_friend(data_files + "time_diff_cm.root", "RECREATE");
    TTree *tree = new TTree("h1", "Friend tree for time difference by proper times of kaons.");

    

}