#include <iostream>
#include <fstream>

#include <TROOT.h>
#include <TTree.h>
#include <TFile.h>
#include <TStyle.h>
#include <TChain.h>
#include <Math/Minimizer.h>
#include <Math/Factory.h>
#include <Math/Functor.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TRatioPlot.h>
#include <TF1.h>

#include "../../Include/Codes/interference.h"
#include "../../Include/Codes/kloe_class.h"
#include "../../Include/const.h"
#include "chain_init.C"

using namespace KLOE;

void resolutions()
{
    TFile file("../../../old_root_files/mctruth.root");
    TTree *tree = (TTree*)file.Get("h1");

    Int_t mctruth = 0;

    tree->SetBranchAddress("mctruth", &mctruth);

    TChain *chain = new TChain("INTERF/h1");

    chain_init(chain);

    UChar_t mcflag;
    Float_t Dtboostlor, Dtmc, ip[3], phi_mom[4], Kchboost[9], Knerec[9];

    chain->SetBranchAddress("mcflag", &mcflag);

    chain->SetBranchAddress("Dtboostlor", &Dtboostlor);
    chain->SetBranchAddress("Dtmc", &Dtmc);

    chain->SetBranchAddress("Broots", &phi_mom[0]);
    chain->SetBranchAddress("Bpx", &phi_mom[1]);
    chain->SetBranchAddress("Bpy", &phi_mom[2]);
    chain->SetBranchAddress("Bpz", &phi_mom[3]);

    chain->SetBranchAddress("ip", ip);

    chain->SetBranchAddress("Kchboost", Kchboost);
    chain->SetBranchAddress("Knerec", Knerec);

    chain->AddFriend("h = h1", "../../../old_root_files/mctruth.root");

    UInt_t nentries = chain->GetEntries();

    UInt_t nbins = 400;
    Double_t x_min = -50.0, x_max = 50.0;
    
    TH1 *res_dt = new TH1D("", ";#Deltat [#tau_{S}];Counts", nbins, x_min, x_max);
    
    for(UInt_t i = 0; i < nentries; i++)
    {
        chain->GetEntry(i);
        tree->GetEntry(i);

        if(mcflag == 1)
        {
            if(mctruth == 1)
            {
                res_dt->Fill(Dtmc - Dtboostlor);
            }
        }
    }

    gStyle->SetOptFit(1011);

    res_dt->Fit("gaus");
    TF1 *fit = (TF1*)res_dt->GetFunction("gaus");

    Double_t deltat_chi2 = fit->GetChisquare();
    Double_t deltat_res = fit->GetParameter(2);
    Double_t deltat_res_err = fit->GetParError(2);


    std::ofstream myfile_dt;
    myfile_dt.open ("resolutions.csv");
    myfile_dt << "Variable (MC - REC),Resolution,Std deviation\n";
    myfile_dt << "\u0394t," << deltat_res << "," << deltat_res_err << "," << deltat_chi2 << ",\n";
    myfile_dt.close();

    TCanvas *c1 = new TCanvas("c1", "", 790, 790);
    c1->cd();

    res_dt->SetLineWidth(3);
    res_dt->SetLineColor(kBlack);
    res_dt->Draw();

    c1->Print("deltat_res.png");

}