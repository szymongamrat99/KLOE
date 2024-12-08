#include "../../Include/const.h"

void pipi_find(UInt_t filenumber = 1, TString directory = "230531_data", TString rootname = "data_stream42_")
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

    tree->SetBranchAddress("mcflag", &mcflag);
    tree->SetBranchAddress("pidmc", pidmc);
    tree->SetBranchAddress("vtxmc", vtxmc);
    tree->SetBranchAddress("mother", mother);

    tree->SetBranchAddress("xvmc", pos_mc[0]);
    tree->SetBranchAddress("yvmc", pos_mc[1]);
    tree->SetBranchAddress("zvmc", pos_mc[2]);

    tree->SetBranchAddress("pxmc", mom_mc[0]);
    tree->SetBranchAddress("pymc", mom_mc[1]);
    tree->SetBranchAddress("pzmc", mom_mc[2]);

    Int_t nentries = (Int_t)tree->GetEntries();

    Int_t kspip, kspim, ksothers, klpip, klpim, klothers, phiks, phikl, phiothers;

    Int_t mctruth_pipi = 0;

    TBranch *b_mctruth_pipi = tree->Branch("mctruth_pipi", &mctruth_pipi, "mctruth_pipi/I");

    for(Int_t i = 0; i < nentries; i++)
    {
        tree->GetEntry(i);

        kspip = 0;
        kspim = 0;
        ksothers = 0;
        klpip = 0;
        klpim = 0;
        klothers = 0;
        phikl = 0;
        phiks = 0;
        phiothers = 0;

        mctruth_pipi = 0;
	
	if(mcflag == 1){
        for(Int_t j = 0; j < ntmc; j++)
        {
                if(mother[vtxmc[j] - 1] == 16)
                {
                    if(pidmc[j] == 8) kspip++;
                        else if(pidmc[j] == 9) kspim++;
                            else ksothers++;
                }
                else if(mother[vtxmc[j] - 1] == 10)
                {
                    if(pidmc[j] == 8) klpip++;
                        else if(pidmc[j] == 9) klpim++;
                            else klothers++;
                }
                else if(mother[vtxmc[j] - 1] == 50)
                {
                    if(pidmc[j] == 16) phiks++;
                        else if(pidmc[j] == 10) phikl++;
                            else phiothers++;
                }
            
        }

        if( phiks == 1 && phikl == 1 && kspim == 1 && kspip == 1 && klpim == 1 && klpip == 1 
            && klothers == 0 && ksothers == 0 && phiothers == 0  ) mctruth_pipi = 1;
        else mctruth_pipi = 0;
	}

        b_mctruth_pipi->Fill();

    }

    file->cd(treedir);
    tree->Write();
    file->Close();
    delete file;

}
