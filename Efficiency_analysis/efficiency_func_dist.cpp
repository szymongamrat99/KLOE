#include "../../Include/const.h"

void efficiency_func_dist()
{
    TString fullname = "",
            dirnamemc = "230528_mc", dirnamedata = "230528_data",
            filenamemc = "mc_stream62_mccard2_", filenamedata = "data_stream42_",
            extension = ".root";

    TChain chain("INTERF/h1");

    for(Int_t i = 1; i <= 11; i++)
    {
        
        fullname = "../../ROOT_files/" + dirnamemc + "/" + filenamemc + i + extension;
        chain.Add(fullname);

        //fullname = "../../ROOT_files/" + dirnamedata + "/" + filenamedata + i + extension;
        //chain.Add(fullname);
    }

    TEfficiency *efficiency;

    Float_t dt_max = 100.0;
    Float_t dt_min = -abs(dt_max);

    UInt_t nbins = (dt_max - dt_min) + 1;

    efficiency = new TEfficiency("efficiency", "", nbins, dt_min, dt_max);

    Float_t Dtboostlor, minv4gam;
    UChar_t mctruth;

    chain.SetBranchAddress("Dtboostlor", &Dtboostlor);
    chain.SetBranchAddress("mctruth", &mctruth);
    chain.SetBranchAddress("minv4gam", &minv4gam);

    Int_t nentries = chain.GetEntries();

    Bool_t passed_cond;

    for(Int_t i = 0; i < nentries; i++)
    {
        chain.GetEntry(i);

        if(mctruth == 6)
        {
            passed_cond = abs(minv4gam - m_k0) < 76;
            efficiency->Fill(passed_cond, Dtboostlor);
        }
    }

    efficiency->SetStatisticOption(TEfficiency::kBUniform);
    efficiency->Draw();

}