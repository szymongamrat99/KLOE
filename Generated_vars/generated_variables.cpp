#include <TTree.h>
#include <TFile.h>


#include "../../Include/const.h"

void generated_variables(UInt_t filenumber = 1, TString directory = "230531_data", TString rootname = "data_stream42_")
{
    TString treedir = "INTERF", treename = "h1", fulltree = "",
            extension = ".root", fullname = "";

    fullname = "../../ROOT_files/" + directory + "/backup/" + rootname + filenumber + extension;
    fulltree = treedir + "/" + treename;

    TFile *file = new TFile(fullname, "UPDATE");
    TTree *tree = (TTree*)file->Get(fulltree);

    //Branches' addresses
    //Bhabha vars
    Int_t ntmc, nvtxmc;
    UInt_t pidmc[200], vtxmc[200], mother[200], mctruth = 0, mcflag = 0;
    Float_t pos_mc[3][200], mom_mc[3][200];

    tree->SetBranchAddress("ntmc", &ntmc);
    tree->SetBranchAddress("nvtxmc", &nvtxmc);

    tree->SetBranchAddress("pidmc", pidmc);
    tree->SetBranchAddress("vtxmc", vtxmc);
    tree->SetBranchAddress("mother", mother);

    tree->SetBranchAddress("xvmc", pos_mc[0]);
    tree->SetBranchAddress("yvmc", pos_mc[1]);
    tree->SetBranchAddress("zvmc", pos_mc[2]);

    tree->SetBranchAddress("pxmc", mom_mc[0]);
    tree->SetBranchAddress("pymc", mom_mc[1]);
    tree->SetBranchAddress("pzmc", mom_mc[2]);

    tree->SetBranchAddress("mctruth", &mctruth);
    tree->SetBranchAddress("mcflag", &mcflag);

    Int_t nentries = (Int_t)tree->GetEntries();

    Float_t Ks_three[9], Kl_three[9], Ks_semi[9], Kl_semi[9], Ks_regen[9], Kl_regen[9];

    Float_t ipmc_three[3], ipmc_semi[3], ipmc_regen[3], ipmc_omega[3], ipmc_other[3],
            Knemc_three[9], Knemc_semi[9], Knemc_regen[9], Kchmc_three[9], Kchmc_semi[9], Kchmc_regen[9];

    TBranch *b_ipmcthree = tree->Branch("ipmc_three", ipmc_three, "ipmc_three[3]/F");
    TBranch *b_ipmcsemi = tree->Branch("ipmc_semi", ipmc_semi, "ipmc_semi[3]/F");
    TBranch *b_ipmcregen = tree->Branch("ipmc_regen", ipmc_regen, "ipmc_regen[3]/F");
    TBranch *b_ipmcomega = tree->Branch("ipmc_omega", ipmc_omega, "ipmc_omega[3]/F");
    TBranch *b_ipmcother = tree->Branch("ipmc_other", ipmc_other, "ipmc_other[3]/F");

    TBranch *b_Knemcthree = tree->Branch("Knemc_three", Knemc_three, "Knemc_three[9]/F");
    TBranch *b_Knemcsemi = tree->Branch("Knemc_semi", Knemc_semi, "Knemc_semi[9]/F");
    TBranch *b_Knemcregen = tree->Branch("Knemc_regen", Knemc_regen, "Knemc_regen[9]/F");

    TBranch *b_Kchmcthree = tree->Branch("Kchmc_three", Kchmc_three, "Kchmc_three[9]/F");
    TBranch *b_Kchmcsemi = tree->Branch("Kchmc_semi", Kchmc_semi, "Kchmc_semi[9]/F");
    TBranch *b_Kchmcregen = tree->Branch("Kchmc_regen", Kchmc_regen, "Kchmc_regen[9]/F");
    for(Int_t i = 0; i < nentries; i++)
    {
        tree->GetEntry(i);


        for(Int_t j = 0; j < 3; j++) ipmc_three[j] = -999.;
        for(Int_t j = 0; j < 3; j++) ipmc_semi[j] = -999.;
        for(Int_t j = 0; j < 3; j++) ipmc_regen[j] = -999.;
        for(Int_t j = 0; j < 3; j++) ipmc_omega[j] = -999.;
        for(Int_t j = 0; j < 3; j++) ipmc_other[j] = -999.;

        for(Int_t j = 0; j < 9; j++) Knemc_three[j] = -999.;
        for(Int_t j = 0; j < 9; j++) Knemc_semi[j] = -999.;
        for(Int_t j = 0; j < 9; j++) Knemc_regen[j] = -999.;

        for(Int_t j = 0; j < 9; j++) Kchmc_three[j] = -999.;
        for(Int_t j = 0; j < 9; j++) Kchmc_semi[j] = -999.;
        for(Int_t j = 0; j < 9; j++) Kchmc_regen[j] = -999.;

        if(mcflag == 1)
        {
            if(mctruth == 4)
            {
                for(Int_t j = 0; j < ntmc; j++)
                {
                    if(mother[vtxmc[j] - 1] == 50)
                    {
                        ipmc_omega[0] = pos_mc[0][vtxmc[j] - 1];
                        ipmc_omega[1] = pos_mc[1][vtxmc[j] - 1];
                        ipmc_omega[2] = pos_mc[2][vtxmc[j] - 1];
                    }
                }
            }

            if(mctruth == 7)
            {
                for(Int_t j = 0; j < ntmc; j++)
                {
                    if(mother[vtxmc[j] - 1] == 50)
                    {
                        ipmc_other[0] = pos_mc[0][vtxmc[j] - 1];
                        ipmc_other[1] = pos_mc[1][vtxmc[j] - 1];
                        ipmc_other[2] = pos_mc[2][vtxmc[j] - 1];
                    }
                }
            }

            if(mctruth == 3)
            {
                for(Int_t j = 0; j < ntmc; j++)
                {
                    if(mother[vtxmc[j] - 1] == 50)
                    {
                        ipmc_regen[0] = pos_mc[0][vtxmc[j] - 1];
                        ipmc_regen[1] = pos_mc[1][vtxmc[j] - 1];
                        ipmc_regen[2] = pos_mc[2][vtxmc[j] - 1];
                    }

                    if(mother[vtxmc[j] - 1] == 50 && pidmc[j] == 10)
                    {
                        Kl_regen[0] = mom_mc[0][j];
                        Kl_regen[1] = mom_mc[1][j];
                        Kl_regen[2] = mom_mc[2][j];
                    }

                    if(mother[vtxmc[j] - 1] == 50 && pidmc[j] == 16)
                    {
                        Ks_regen[0] = mom_mc[0][j];
                        Ks_regen[1] = mom_mc[1][j];
                        Ks_regen[2] = mom_mc[2][j];
                    }

                    if(mother[vtxmc[j] - 1] == 10 && pidmc[j] == 7)
                    {
                        Knemc_regen[0] = Kl_regen[0];
                        Knemc_regen[1] = Kl_regen[1];
                        Knemc_regen[2] = Kl_regen[2];
                        Knemc_regen[3] = sqrt(pow(Kl_regen[0],2) + pow(Kl_regen[1],2) + pow(Kl_regen[2],2) + pow(m_k0,2));
                        Knemc_regen[4] = sqrt(pow(Kl_regen[0],2) + pow(Kl_regen[1],2) + pow(Kl_regen[2],2));
                        Knemc_regen[5] = m_k0;
                        Knemc_regen[6] = pos_mc[0][vtxmc[j] - 1];
                        Knemc_regen[7] = pos_mc[1][vtxmc[j] - 1];
                        Knemc_regen[8] = pos_mc[2][vtxmc[j] - 1];
                    }

                    if(mother[vtxmc[j] - 1] == 16 && pidmc[j] == 7)
                    {
                        Knemc_regen[0] = Ks_regen[0];
                        Knemc_regen[1] = Ks_regen[1];
                        Knemc_regen[2] = Ks_regen[2];
                        Knemc_regen[3] = sqrt(pow(Ks_regen[0],2) + pow(Ks_regen[1],2) + pow(Ks_regen[2],2) + pow(m_k0,2));
                        Knemc_regen[4] = sqrt(pow(Ks_regen[0],2) + pow(Ks_regen[1],2) + pow(Ks_regen[2],2));
                        Knemc_regen[5] = m_k0;
                        Knemc_regen[6] = pos_mc[0][vtxmc[j] - 1];
                        Knemc_regen[7] = pos_mc[1][vtxmc[j] - 1];
                        Knemc_regen[8] = pos_mc[2][vtxmc[j] - 1];
                    }

                    if(mother[vtxmc[j] - 1] == 10 && (pidmc[j] == 8 || pidmc[j] == 9))
                    {
                        Kchmc_regen[0] = Kl_regen[0];
                        Kchmc_regen[1] = Kl_regen[1];
                        Kchmc_regen[2] = Kl_regen[2];
                        Kchmc_regen[3] = sqrt(pow(Kl_regen[0],2) + pow(Kl_regen[1],2) + pow(Kl_regen[2],2) + pow(m_k0,2));
                        Kchmc_regen[4] = sqrt(pow(Kl_regen[0],2) + pow(Kl_regen[1],2) + pow(Kl_regen[2],2));
                        Kchmc_regen[5] = m_k0;
                        Kchmc_regen[6] = pos_mc[0][vtxmc[j] - 1];
                        Kchmc_regen[7] = pos_mc[1][vtxmc[j] - 1];
                        Kchmc_regen[8] = pos_mc[2][vtxmc[j] - 1];
                    }

                    if(mother[vtxmc[j] - 1] == 16 && (pidmc[j] == 8 || pidmc[j] == 9))
                    {
                        Kchmc_regen[0] = Ks_regen[0];
                        Kchmc_regen[1] = Ks_regen[1];
                        Kchmc_regen[2] = Ks_regen[2];
                        Kchmc_regen[3] = sqrt(pow(Ks_regen[0],2) + pow(Ks_regen[1],2) + pow(Ks_regen[2],2) + pow(m_k0,2));
                        Kchmc_regen[4] = sqrt(pow(Ks_regen[0],2) + pow(Ks_regen[1],2) + pow(Ks_regen[2],2));
                        Kchmc_regen[5] = m_k0;
                        Kchmc_regen[6] = pos_mc[0][vtxmc[j] - 1];
                        Kchmc_regen[7] = pos_mc[1][vtxmc[j] - 1];
                        Kchmc_regen[8] = pos_mc[2][vtxmc[j] - 1];
                    }
                }
            }

            if(mctruth == 5)
            {
                for(Int_t j = 0; j < ntmc; j++)
                {
                    if(mother[vtxmc[j] - 1] == 50)
                    {
                        ipmc_three[0] = pos_mc[0][j];
                        ipmc_three[1] = pos_mc[1][j];
                        ipmc_three[2] = pos_mc[2][j];
                    }

                    if(mother[vtxmc[j] - 1] == 50 && pidmc[j] == 10)
                    {
                        Kl_three[0] = mom_mc[0][j];
                        Kl_three[1] = mom_mc[1][j];
                        Kl_three[2] = mom_mc[2][j];
                    }

                    if(mother[vtxmc[j] - 1] == 50 && pidmc[j] == 16)
                    {
                        Ks_three[0] = mom_mc[0][j];
                        Ks_three[1] = mom_mc[1][j];
                        Ks_three[2] = mom_mc[2][j];
                    }

                    if(mother[vtxmc[j] - 1] == 10 && pidmc[j] == 7)
                    {
                        Knemc_three[0] = Kl_three[0];
                        Knemc_three[1] = Kl_three[1];
                        Knemc_three[2] = Kl_three[2];
                        Knemc_three[3] = sqrt(pow(Kl_three[0],2) + pow(Kl_three[1],2) + pow(Kl_three[2],2) + pow(m_k0,2));
                        Knemc_three[4] = sqrt(pow(Kl_three[0],2) + pow(Kl_three[1],2) + pow(Kl_three[2],2));
                        Knemc_three[5] = m_k0;
                        Knemc_three[6] = pos_mc[0][vtxmc[j] - 1];
                        Knemc_three[7] = pos_mc[1][vtxmc[j] - 1];
                        Knemc_three[8] = pos_mc[2][vtxmc[j] - 1];

                    }

                    if(mother[vtxmc[j] - 1] == 16 && pidmc[j] == 7)
                    {
                        Knemc_three[0] = Ks_three[0];
                        Knemc_three[1] = Ks_three[1];
                        Knemc_three[2] = Ks_three[2];
                        Knemc_three[3] = sqrt(pow(Ks_three[0],2) + pow(Ks_three[1],2) + pow(Ks_three[2],2) + pow(m_k0,2));
                        Knemc_three[4] = sqrt(pow(Ks_three[0],2) + pow(Ks_three[1],2) + pow(Ks_three[2],2));
                        Knemc_three[5] = m_k0;
                        Knemc_three[6] = pos_mc[0][vtxmc[j] - 1];
                        Knemc_three[7] = pos_mc[1][vtxmc[j] - 1];
                        Knemc_three[8] = pos_mc[2][vtxmc[j] - 1];

                    }

                    if(mother[vtxmc[j] - 1] == 10 && (pidmc[j] == 8 || pidmc[j] == 9))
                    {
                        Kchmc_three[0] = Kl_three[0];
                        Kchmc_three[1] = Kl_three[1];
                        Kchmc_three[2] = Kl_three[2];
                        Kchmc_three[3] = sqrt(pow(Kl_three[0],2) + pow(Kl_three[1],2) + pow(Kl_three[2],2) + pow(m_k0,2));
                        Kchmc_three[4] = sqrt(pow(Kl_three[0],2) + pow(Kl_three[1],2) + pow(Kl_three[2],2));
                        Kchmc_three[5] = m_k0;
                        Kchmc_three[6] = pos_mc[0][vtxmc[j] - 1];
                        Kchmc_three[7] = pos_mc[1][vtxmc[j] - 1];
                        Kchmc_three[8] = pos_mc[2][vtxmc[j] - 1];
                    }

                    if(mother[vtxmc[j] - 1] == 16 && (pidmc[j] == 8 || pidmc[j] == 9))
                    {
                        Kchmc_three[0] = Ks_three[0];
                        Kchmc_three[1] = Ks_three[1];
                        Kchmc_three[2] = Ks_three[2];
                        Kchmc_three[3] = sqrt(pow(Ks_three[0],2) + pow(Ks_three[1],2) + pow(Ks_three[2],2) + pow(m_k0,2));
                        Kchmc_three[4] = sqrt(pow(Ks_three[0],2) + pow(Ks_three[1],2) + pow(Ks_three[2],2));
                        Kchmc_three[5] = m_k0;
                        Kchmc_three[6] = pos_mc[0][vtxmc[j] - 1];
                        Kchmc_three[7] = pos_mc[1][vtxmc[j] - 1];
                        Kchmc_three[8] = pos_mc[2][vtxmc[j] - 1];
                    }
                }
            }

            if(mctruth == 6)
            {
                for(Int_t j = 0; j < ntmc; j++)
                {
                    if(mother[vtxmc[j] - 1] == 50)
                    {
                        ipmc_semi[0] = pos_mc[0][vtxmc[j] - 1];
                        ipmc_semi[1] = pos_mc[1][vtxmc[j] - 1];
                        ipmc_semi[2] = pos_mc[2][vtxmc[j] - 1];
                    }

                    if(mother[vtxmc[j] - 1] == 50 && pidmc[j] == 10)
                    {
                        Kl_semi[0] = mom_mc[0][j];
                        Kl_semi[1] = mom_mc[1][j];
                        Kl_semi[2] = mom_mc[2][j];
                    }

                    if(mother[vtxmc[j] - 1] == 50 && pidmc[j] == 16)
                    {
                        Ks_semi[0] = mom_mc[0][j];
                        Ks_semi[1] = mom_mc[1][j];
                        Ks_semi[2] = mom_mc[2][j];
                    }

                    if(mother[vtxmc[j] - 1] == 10 && pidmc[j] == 7)
                    {
                        Knemc_semi[0] = Kl_semi[0];
                        Knemc_semi[1] = Kl_semi[1];
                        Knemc_semi[2] = Kl_semi[2];
                        Knemc_semi[3] = sqrt(pow(Kl_semi[0],2) + pow(Kl_semi[1],2) + pow(Kl_semi[2],2) + pow(m_k0,2));
                        Knemc_semi[4] = sqrt(pow(Kl_semi[0],2) + pow(Kl_semi[1],2) + pow(Kl_semi[2],2));
                        Knemc_semi[5] = m_k0;
                        Knemc_semi[6] = pos_mc[0][vtxmc[j] - 1];
                        Knemc_semi[7] = pos_mc[1][vtxmc[j] - 1];
                        Knemc_semi[8] = pos_mc[2][vtxmc[j] - 1];
                    }

                    if(mother[vtxmc[j] - 1] == 16 && pidmc[j] == 7)
                    {
                        Knemc_semi[0] = Ks_semi[0];
                        Knemc_semi[1] = Ks_semi[1];
                        Knemc_semi[2] = Ks_semi[2];
                        Knemc_semi[3] = sqrt(pow(Ks_semi[0],2) + pow(Ks_semi[1],2) + pow(Ks_semi[2],2) + pow(m_k0,2));
                        Knemc_semi[4] = sqrt(pow(Ks_semi[0],2) + pow(Ks_semi[1],2) + pow(Ks_semi[2],2));
                        Knemc_semi[5] = m_k0;
                        Knemc_semi[6] = pos_mc[0][vtxmc[j] - 1];
                        Knemc_semi[7] = pos_mc[1][vtxmc[j] - 1];
                        Knemc_semi[8] = pos_mc[2][vtxmc[j] - 1];
                    }

                    if(mother[vtxmc[j] - 1] == 10 && (pidmc[j] == 8 || pidmc[j] == 9))
                    {
                        Kchmc_semi[0] = Kl_semi[0];
                        Kchmc_semi[1] = Kl_semi[1];
                        Kchmc_semi[2] = Kl_semi[2];
                        Kchmc_semi[3] = sqrt(pow(Kl_semi[0],2) + pow(Kl_semi[1],2) + pow(Kl_semi[2],2) + pow(m_k0,2));
                        Kchmc_semi[4] = sqrt(pow(Kl_semi[0],2) + pow(Kl_semi[1],2) + pow(Kl_semi[2],2));
                        Kchmc_semi[5] = m_k0;
                        Kchmc_semi[6] = pos_mc[0][vtxmc[j] - 1];
                        Kchmc_semi[7] = pos_mc[1][vtxmc[j] - 1];
                        Kchmc_semi[8] = pos_mc[2][vtxmc[j] - 1];
                    }

                    if(mother[vtxmc[j] - 1] == 16 && (pidmc[j] == 8 || pidmc[j] == 9))
                    {
                        Kchmc_semi[0] = Ks_semi[0];
                        Kchmc_semi[1] = Ks_semi[1];
                        Kchmc_semi[2] = Ks_semi[2];
                        Kchmc_semi[3] = sqrt(pow(Ks_semi[0],2) + pow(Ks_semi[1],2) + pow(Ks_semi[2],2) + pow(m_k0,2));
                        Kchmc_semi[4] = sqrt(pow(Ks_semi[0],2) + pow(Ks_semi[1],2) + pow(Ks_semi[2],2));
                        Kchmc_semi[5] = m_k0;
                        Kchmc_semi[6] = pos_mc[0][vtxmc[j] - 1];
                        Kchmc_semi[7] = pos_mc[1][vtxmc[j] - 1];
                        Kchmc_semi[8] = pos_mc[2][vtxmc[j] - 1];
                    }
                }
            }

        }

        b_ipmcthree->Fill();
        b_ipmcsemi->Fill();
        b_ipmcregen->Fill();
        b_ipmcomega->Fill();
        b_ipmcother->Fill();

        b_Knemcthree->Fill();
        b_Knemcsemi->Fill();
        b_Knemcregen->Fill();

        b_Kchmcthree->Fill();
        b_Kchmcsemi->Fill();
        b_Kchmcregen->Fill();
    }

    file->cd(treedir);
    tree->Write();
    file->Close();
    delete file;

}
