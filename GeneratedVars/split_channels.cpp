#include "../../Include/const.h"
#include <TMath.h>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>

Int_t split_channels()
{
    TFile file("mctruth.root", "recreate");
    TTree *tree = new TTree("h1", "Mctruth for all channels");

    Int_t mctruth_new = 0;

    tree->Branch("mctruth", &mctruth_new, "mctruth/I");

    Int_t ntmc, nvtxmc;
    UChar_t pidmc[200], vtxmc[200], mother[200], mctruth, mcflag;

    chain.SetBranchAddress("ntmc", &ntmc);
    chain.SetBranchAddress("nvtxmc", &nvtxmc);

    chain.SetBranchAddress("pidmc", pidmc);
    chain.SetBranchAddress("vtxmc", vtxmc);
    chain.SetBranchAddress("mother", mother);

    chain.SetBranchAddress("mcflag", &mcflag);
    chain.SetBranchAddress("mctruth", &mctruth);

    Int_t nentries = (Int_t)chain.GetEntries();

    UInt_t Ks = 0, Kl = 0, Ksregen = 0, piplusks = 0, pipluskl = 0, piminusks = 0, piminuskl = 0, 
           muonplusks = 0, muonpluskl = 0, muonminusks = 0, muonminuskl = 0, electronks = 0, electronkl = 0, 
           positronks = 0, positronkl = 0, pi0ks = 0, pi0kl = 0, pi0phi = 0, piplusphi = 0, piminusphi = 0, otherphi = 0,
           otherkl = 0, otherks = 0, gammaphi = 0;

    Bool_t signal_cond, regen_cond, omega_cond, three_cond, semi_cond, other_cond;

    for(Int_t i = 0; i < nentries; i++)
    {
        chain.GetEntry(i);

        if(mcflag == 1 && mctruth != 0 && mctruth != 2)
        {
            for(Int_t j = 0; j < ntmc; j++)
            {
                if(mother[vtxmc[j] - 1] == 50)
                {
                    if(pidmc[j] == 10) Kl++;
                    else if(pidmc[j] == 16) Ks++;
                    else if(pidmc[j] == 7) pi0phi++;
                    else if(pidmc[j] == 8) piplusphi++;
                    else if(pidmc[j] == 9) piminusphi++;
                    else if(pidmc[j] == 1) gammaphi++;
                    else otherphi++;
                }
                else if(mother[vtxmc[j] - 1] == 10)
                {
                    if(pidmc[j] == 16) Ksregen++;
                    else if(pidmc[j] == 7) pi0kl++;
                    else if(pidmc[j] == 8) pipluskl++;
                    else if(pidmc[j] == 9) piminuskl++;
                    else if(pidmc[j] == 5) muonpluskl++;
                    else if(pidmc[j] == 6) muonminuskl++;
                    else if(pidmc[j] == 2) positronkl++;
                    else if(pidmc[j] == 3) electronkl++;
                    else otherkl++;
                }
                else if(mother[vtxmc[j] - 1] == 16)
                {
                    if(pidmc[j] == 7) pi0ks++;
                    else if(pidmc[j] == 8) piplusks++;
                    else if(pidmc[j] == 9) piminusks++;
                    else if(pidmc[j] == 5) muonplusks++;
                    else if(pidmc[j] == 6) muonminusks++;
                    else if(pidmc[j] == 2) positronks++;
                    else if(pidmc[j] == 3) electronks++;
                    else otherks++;
                }  
            }

            signal_cond = (pi0phi == 0 && piplusphi == 0 && piminusphi == 0 && otherphi == 0 && otherks == 0 && otherkl == 0 &&
                           positronkl + positronks == 0 && electronkl + electronks == 0 && muonminuskl + muonminusks == 0 &&
                           muonpluskl + muonplusks == 0 && Ksregen == 0 && 
                           Ks == 1 && Kl == 1 && 
                           ((pi0ks == 2 && pipluskl == 1 && piminuskl == 1 && pi0kl == 0 && piplusks == 0 && piminusks == 0 ) || 
                            (pi0kl == 2 && piplusks == 1 && piminusks == 1 && pi0ks == 0 && pipluskl == 0 && piminuskl == 0) ));

            regen_cond = (Ksregen == 1 && Ks == 1 && Kl == 1);

            omega_cond = (pi0phi == 2 && piplusphi == 1 && piminusphi == 1 && otherphi == 0 && otherks == 0 && otherkl == 0 &&
                           positronkl + positronks == 0 && electronkl + electronks == 0 && muonminuskl + muonminusks == 0 &&
                           muonpluskl + muonplusks == 0 && Ksregen == 0 && 
                           Ks == 0 && Kl == 0 && pi0ks == 0 && pi0kl == 0 && pipluskl + piplusks == 0 && piminuskl + piminusks == 0);

            three_cond = (pi0phi == 0 && piplusphi == 0 && piminusphi == 0 && otherphi == 0 && otherks == 0 && otherkl == 0 &&
                           positronkl + positronks == 0 && electronkl + electronks == 0 && muonminuskl + muonminusks == 0 &&
                           muonpluskl + muonplusks == 0 && Ksregen == 0 && 
                           Ks == 1 && Kl == 1 && (pi0kl == 3 && piplusks == 1 && piminusks == 1 && pi0ks == 0 && pipluskl == 0 && piminuskl == 0));

            semi_cond = (pi0phi == 0 && piplusphi == 0 && piminusphi == 0 && otherphi == 0 && otherks == 0 && otherkl == 0 &&
                           Ksregen == 0 && Ks == 1 && Kl == 1 && 
                           ((pi0ks == 2 && positronkl == 1 && piminuskl == 1 && pi0kl == 0 ) ||
                            (pi0ks == 2 && pipluskl == 1 && electronkl == 1 && pi0kl == 0 ) ||
                            (pi0ks == 2 && pipluskl == 1 && muonminuskl == 1 && pi0kl == 0 ) ||
                            (pi0ks == 2 && piminuskl == 1 && muonpluskl == 1 && pi0kl == 0 ) ||
                            (pi0kl == 2 && positronks == 1 && piminusks == 1 && pi0ks == 0 ) || 
                            (pi0kl == 2 && piplusks == 1 && electronks == 1 && pi0ks == 0 ) ||
                            (pi0kl == 2 && piplusks == 1 && muonminusks == 1 && pi0ks == 0 ) ||
                            (pi0kl == 2 && piminusks == 1 && muonplusks == 1 && pi0ks == 0 ) ));

            if( signal_cond ) mctruth_new = 1;
            else if( regen_cond ) mctruth_new = 3;
            else if( omega_cond ) mctruth_new = 4;
            else if( three_cond ) mctruth_new = 5;
            else if( semi_cond ) mctruth_new = 6;
            else mctruth_new = 7;

        }
        else mctruth_new = 0;

        tree->Fill();

        Ks = 0;
        Kl = 0;
        Ksregen = 0;
        pi0kl = 0;
        pi0ks = 0;
        pi0phi = 0;
        pipluskl = 0;
        piplusks = 0;
        piplusphi = 0;
        piminuskl = 0;
        piminusks = 0;
        piminusphi = 0;
        electronkl = 0;
        electronks = 0;
        positronkl = 0;
        positronks = 0;
        muonminuskl = 0;
        muonpluskl = 0;
        muonminusks = 0;
        muonplusks = 0;
        otherkl = 0;
        otherks = 0;
        otherphi = 0;

    }

    file.Write();
    file.Close();

    return 0;
}