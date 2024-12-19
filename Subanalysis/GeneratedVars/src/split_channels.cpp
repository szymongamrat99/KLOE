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

    chain = new TChain("h1");
    chain_init(chain, firstFile, lastFile);

    chain->SetBranchAddress("nTMC", &interfcommon_.ntmc);
    chain->SetBranchAddress("nVtxMC", &interfcommon_.nvtxmc);

    chain->SetBranchAddress("PidMC", interfcommon_.pidmc);
    chain->SetBranchAddress("VtxMC", interfcommon_.vtxmc);
    chain->SetBranchAddress("Mother", interfcommon_.mother);
    // chain->SetBranchAddress("mctruth", &baseKin.mctruth);

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

    Bool_t signal_cond, regen_cond, omega_cond, three_cond, semi_cond, other_cond, pipi_cond;

    for (Int_t i = 0; i < nentries; i++)
    {
        chain->GetEntry(i);

        interfcommon_.mcflag = 1;

        if (interfcommon_.mcflag == 1 /*&& baseKin.mctruth != 0 && baseKin.mctruth != 2*/)
        {
            for (Int_t j = 0; j < interfcommon_.ntmc; j++)
            {
                if (interfcommon_.mother[interfcommon_.vtxmc[j] - 1] == 50)
                {
                    if (interfcommon_.pidmc[j] == 10)
                        Kl++;
                    else if (interfcommon_.pidmc[j] == 16)
                        Ks++;
                    else if (interfcommon_.pidmc[j] == 7)
                        pi0phi++;
                    else if (interfcommon_.pidmc[j] == 8)
                        piplusphi++;
                    else if (interfcommon_.pidmc[j] == 9)
                        piminusphi++;
                    else if (interfcommon_.pidmc[j] == 1)
                        gammaphi++;
                    else
                        otherphi++;
                }
                else if (interfcommon_.mother[interfcommon_.vtxmc[j] - 1] == 10)
                {
                    if (interfcommon_.pidmc[j] == 16)
                        Ksregen++;
                    else if (interfcommon_.pidmc[j] == 7)
                        pi0kl++;
                    else if (interfcommon_.pidmc[j] == 8)
                        pipluskl++;
                    else if (interfcommon_.pidmc[j] == 9)
                        piminuskl++;
                    else if (interfcommon_.pidmc[j] == 5)
                        muonpluskl++;
                    else if (interfcommon_.pidmc[j] == 6)
                        muonminuskl++;
                    else if (interfcommon_.pidmc[j] == 2)
                        positronkl++;
                    else if (interfcommon_.pidmc[j] == 3)
                        electronkl++;
                    else
                        otherkl++;
                }
                else if (interfcommon_.mother[interfcommon_.vtxmc[j] - 1] == 16)
                {
                    if (interfcommon_.pidmc[j] == 7)
                        pi0ks++;
                    else if (interfcommon_.pidmc[j] == 8)
                        piplusks++;
                    else if (interfcommon_.pidmc[j] == 9)
                        piminusks++;
                    else if (interfcommon_.pidmc[j] == 5)
                        muonplusks++;
                    else if (interfcommon_.pidmc[j] == 6)
                        muonminusks++;
                    else if (interfcommon_.pidmc[j] == 2)
                        positronks++;
                    else if (interfcommon_.pidmc[j] == 3)
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

            pipi_cond = (pi0phi == 0 && piplusphi == 0 && piminusphi == 0 && otherphi == 0 && otherks == 0 && otherkl == 0 &&
                         positronkl + positronks == 0 && electronkl + electronks == 0 && muonminuskl + muonminusks == 0 &&
                         muonpluskl + muonplusks == 0 && Ksregen == 0 &&
                         Ks == 1 && Kl == 1 &&
                         ((piminusks == 1 && piplusks == 1 && pipluskl == 1 && piminuskl == 1 && pi0kl == 0 && pi0ks == 0)));

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
            else if (pipi_cond)
                baseKin.mctruth_int = 8;
            else
                baseKin.mctruth_int = 7;
        }
        else if (interfcommon_.mcflag == 0)
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