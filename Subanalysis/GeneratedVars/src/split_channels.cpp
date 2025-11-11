#include "../inc/genvars.hpp"

Int_t split_channels(TChain &chain, Controls::DataType &data_type, ErrorHandling::ErrorLogs &logger, KLOE::pm00 &Obj)
{
    TFile
        *file;
    TTree
        *tree;

    KLOE::BaseKinematics baseKin;

    Int_t ntmc, nvtxmc;
    UChar_t pidmcOld[50], vtxmcOld[50], motherOld[50], mctruth, mcflag;

    chain.SetBranchAddress("ntmc", &ntmc);
    chain.SetBranchAddress("nvtxmc", &nvtxmc);

    chain.SetBranchAddress("pidmcOld", pidmcOld);
    chain.SetBranchAddress("vtxmcOld", vtxmcOld);
    chain.SetBranchAddress("motherOld", motherOld);
    chain.SetBranchAddress("mctruth", &mctruth);
    chain.SetBranchAddress("mcflag", &mcflag);

    Int_t nentries = (Int_t)chain.GetEntries();

    // =================================================================================

    // Creation of filename for the analysis step
    std::string
        datestamp = Obj.getCurrentDate(),
        name = "";

    name = Paths::base_path + Paths::root_files_dir + Filenames::mctruth_filename + datestamp + "_" + int(data_type) + Paths::ext_root;
    // -----------------------------------------------------------------------------------------

    file = new TFile(name.c_str(), "recreate");
    tree = new TTree(Filenames::gen_vars_tree, "Mctruth for all channels");

    tree->Branch("mctruth", &baseKin.mctruth_int, "mctruth/I");

    // =================================================================================

    Bool_t signal_cond, regen_cond, omega_cond, three_cond, semi_cond, other_cond, pipi_cond;

    for (Int_t i = 0; i < nentries; i++)
    {
        chain.GetEntry(i);

        UInt_t Ks = 0, Kl = 0, Ksregen = 0, piplusks = 0, pipluskl = 0, piminusks = 0, piminuskl = 0,
               muonplusks = 0, muonpluskl = 0, muonminusks = 0, muonminuskl = 0, electronks = 0, electronkl = 0,
               positronks = 0, positronkl = 0, pi0ks = 0, pi0kl = 0, pi0phi = 0, piplusphi = 0, piminusphi = 0, otherphi = 0,
               otherkl = 0, otherks = 0, gammaphi = 0;

        if (mcflag == 1 && mctruth != 0)
        {
            for (Int_t j = 0; j < ntmc; j++)
            {
                if (motherOld[vtxmcOld[j] - 1] == 50)
                {
                    if (pidmcOld[j] == 10)
                        Kl++;
                    else if (pidmcOld[j] == 16)
                        Ks++;
                    else if (pidmcOld[j] == 7)
                        pi0phi++;
                    else if (pidmcOld[j] == 8)
                        piplusphi++;
                    else if (pidmcOld[j] == 9)
                        piminusphi++;
                    else if (pidmcOld[j] == 1)
                        gammaphi++;
                    else
                        otherphi++;
                }
                else if (motherOld[vtxmcOld[j] - 1] == 10)
                {
                    if (pidmcOld[j] == 16)
                        Ksregen++;
                    else if (pidmcOld[j] == 7)
                        pi0kl++;
                    else if (pidmcOld[j] == 8)
                        pipluskl++;
                    else if (pidmcOld[j] == 9)
                        piminuskl++;
                    else if (pidmcOld[j] == 5)
                        muonpluskl++;
                    else if (pidmcOld[j] == 6)
                        muonminuskl++;
                    else if (pidmcOld[j] == 2)
                        positronkl++;
                    else if (pidmcOld[j] == 3)
                        electronkl++;
                    else
                        otherkl++;
                }
                else if (motherOld[vtxmcOld[j] - 1] == 16)
                {
                    if (pidmcOld[j] == 7)
                        pi0ks++;
                    else if (pidmcOld[j] == 8)
                        piplusks++;
                    else if (pidmcOld[j] == 9)
                        piminusks++;
                    else if (pidmcOld[j] == 5)
                        muonplusks++;
                    else if (pidmcOld[j] == 6)
                        muonminusks++;
                    else if (pidmcOld[j] == 2)
                        positronks++;
                    else if (pidmcOld[j] == 3)
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

            if (signal_cond && mctruth == 1)
                baseKin.mctruth_int = 1;
            else if (signal_cond && mctruth == 2)
                baseKin.mctruth_int = 2;
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
        else if (mcflag == 0)
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

    Utils::properties["variables"]["tree"]["filename"]["mctruth"] = (std::string)name;
    Utils::properties["variables"]["tree"]["treename"]["mctruth"] = (std::string)Filenames::gen_vars_tree;

    Utils::properties["lastScript"] = "Mctruth division of files.";
    Utils::properties["lastUpdate"] = Obj.getCurrentTimestamp();

    std::ofstream outfile(Paths::propName);
    outfile << Utils::properties.dump(4);
    outfile.close();

    return 0;
}
