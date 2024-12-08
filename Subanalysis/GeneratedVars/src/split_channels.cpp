#include <TROOT.h>

#include "../inc/genvars.hpp"

Int_t split_channels(Int_t firstFile, Int_t lastFile)
{
    TFile
        *file;
    TTree
        *tree;
    TChain
        *chain;

    BaseKinematics baseKin;

    chain = new TChain("INTERF/h1");
    chain_init(chain, firstFile, lastFile);

    chain->SetBranchAddress("ntmc", &baseKin.ntmc);
    chain->SetBranchAddress("nvtxmc", &baseKin.nvtxmc);

    chain->SetBranchAddress("pidmc", baseKin.pidmc);
    chain->SetBranchAddress("vtxmc", baseKin.vtxmc);
    chain->SetBranchAddress("mother", baseKin.mother);

    chain->SetBranchAddress("mcflag", &baseKin.mcflag);
    chain->SetBranchAddress("mctruth", &baseKin.mctruth);

    Int_t nentries = (Int_t)chain->GetEntries();

    // =================================================================================

    TString name = mctruth_filename + std::to_string(firstFile) + "_" + std::to_string(lastFile) + ext_root;

    file = new TFile(gen_vars_dir + root_files_dir + name, "recreate");
    tree = new TTree(gen_vars_tree, "Mctruth for all channels");

    tree->Branch("mctruth", &baseKin.mctruth_int, "mctruth/I");

    // =================================================================================

    UInt_t Ks = 0, Kl = 0, Ksregen = 0, piplusks = 0, pipluskl = 0, piminusks = 0, piminuskl = 0,
           muonplusks = 0, muonpluskl = 0, muonminusks = 0, muonminuskl = 0, electronks = 0, electronkl = 0,
           positronks = 0, positronkl = 0, pi0ks = 0, pi0kl = 0, pi0phi = 0, piplusphi = 0, piminusphi = 0, otherphi = 0,
           otherkl = 0, otherks = 0, gammaphi = 0;

    Bool_t signal_cond, regen_cond, omega_cond, three_cond, semi_cond, other_cond;

    for (Int_t i = 0; i < nentries; i++)
    {
        chain->GetEntry(i);

        if (baseKin.mcflag == 1 && baseKin.mctruth != 0 && baseKin.mctruth != 2)
        {
            for (Int_t j = 0; j < baseKin.ntmc; j++)
            {
                if (baseKin.mother[baseKin.vtxmc[j] - 1] == 50)
                {
                    if (baseKin.pidmc[j] == 10)
                        Kl++;
                    else if (baseKin.pidmc[j] == 16)
                        Ks++;
                    else if (baseKin.pidmc[j] == 7)
                        pi0phi++;
                    else if (baseKin.pidmc[j] == 8)
                        piplusphi++;
                    else if (baseKin.pidmc[j] == 9)
                        piminusphi++;
                    else if (baseKin.pidmc[j] == 1)
                        gammaphi++;
                    else
                        otherphi++;
                }
                else if (baseKin.mother[baseKin.vtxmc[j] - 1] == 10)
                {
                    if (baseKin.pidmc[j] == 16)
                        Ksregen++;
                    else if (baseKin.pidmc[j] == 7)
                        pi0kl++;
                    else if (baseKin.pidmc[j] == 8)
                        pipluskl++;
                    else if (baseKin.pidmc[j] == 9)
                        piminuskl++;
                    else if (baseKin.pidmc[j] == 5)
                        muonpluskl++;
                    else if (baseKin.pidmc[j] == 6)
                        muonminuskl++;
                    else if (baseKin.pidmc[j] == 2)
                        positronkl++;
                    else if (baseKin.pidmc[j] == 3)
                        electronkl++;
                    else
                        otherkl++;
                }
                else if (baseKin.mother[baseKin.vtxmc[j] - 1] == 16)
                {
                    if (baseKin.pidmc[j] == 7)
                        pi0ks++;
                    else if (baseKin.pidmc[j] == 8)
                        piplusks++;
                    else if (baseKin.pidmc[j] == 9)
                        piminusks++;
                    else if (baseKin.pidmc[j] == 5)
                        muonplusks++;
                    else if (baseKin.pidmc[j] == 6)
                        muonminusks++;
                    else if (baseKin.pidmc[j] == 2)
                        positronks++;
                    else if (baseKin.pidmc[j] == 3)
                        electronks++;
                    else
                        otherks++;
                }
            }

            signal_cond = (pi0phi == 0 && piplusphi == 0 && piminusphi == 0 && otherphi == 0 && otherks == 0 && otherkl == 0 &&
                           positronkl + positronks == 0 && electronkl + electronks == 0 && muonminuskl + muonminusks == 0 &&
                           muonpluskl + muonplusks == 0 && Ksregen == 0 &&
                           Ks == 1 && Kl == 1 &&
                           ((pi0ks == 2 && pipluskl == 1 && piminuskl == 1 && pi0kl == 0 && piplusks == 0 && piminusks == 0) ||
                            (pi0kl == 2 && piplusks == 1 && piminusks == 1 && pi0ks == 0 && pipluskl == 0 && piminuskl == 0)));

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
                         ((pi0ks == 2 && positronkl == 1 && piminuskl == 1 && pi0kl == 0) ||
                          (pi0ks == 2 && pipluskl == 1 && electronkl == 1 && pi0kl == 0) ||
                          (pi0ks == 2 && pipluskl == 1 && muonminuskl == 1 && pi0kl == 0) ||
                          (pi0ks == 2 && piminuskl == 1 && muonpluskl == 1 && pi0kl == 0) ||
                          (pi0kl == 2 && positronks == 1 && piminusks == 1 && pi0ks == 0) ||
                          (pi0kl == 2 && piplusks == 1 && electronks == 1 && pi0ks == 0) ||
                          (pi0kl == 2 && piplusks == 1 && muonminusks == 1 && pi0ks == 0) ||
                          (pi0kl == 2 && piminusks == 1 && muonplusks == 1 && pi0ks == 0)));

            if (signal_cond)
                baseKin.mctruth_int = 1;
            else if (regen_cond)
                baseKin.mctruth_int = 3;
            else if (omega_cond)
                baseKin.mctruth_int = 4;
            else if (three_cond)
                baseKin.mctruth_int = 5;
            else if (semi_cond)
                baseKin.mctruth_int = 6;
            else
                baseKin.mctruth_int = 7;
        }
        else if (baseKin.mcflag == 1 && baseKin.mctruth == 2)
            baseKin.mctruth_int = 2;
        else if (baseKin.mcflag == 0 || (baseKin.mcflag == 1 && baseKin.mctruth == 0))
            baseKin.mctruth_int = 0;

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

    file->Write();
    file->Close();

    return 0;
}