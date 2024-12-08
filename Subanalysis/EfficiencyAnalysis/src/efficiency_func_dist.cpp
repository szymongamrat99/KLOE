#include "TEfficiency.h"
#include <TFile.h>
#include <TGraphAsymmErrors.h>
#include <TH1F.h>
#include <TCanvas.h>

#include <const.h>
#include <lorentz_transf.h>

void efficiency_func_dist(UInt_t first_file, UInt_t last_file)
{
	TChain *chain = new TChain("INTERF/h1");
	chain_init(chain, first_file, last_file);

    TFile file_tri("../Neutrec/neuvtx_tri_kin_fit_1_56_10_9_1.root");
    TTree *tree_tri = (TTree*)file_tri.Get("h_tri_kin_fit");

    TFile file_corr("correction_factor.root","recreate");

    Int_t done4 = 0, g4takentri_kinfit[4];
    Float_t fourKnetri[10] = {0.}, iptri_kinfit[3];

    tree_tri->SetBranchAddress("done4_kinfit", &done4);
    tree_tri->SetBranchAddress("fourKnetri_kinfit", fourKnetri);
    tree_tri->SetBranchAddress("g4takentri_kinfit", g4takentri_kinfit);
    tree_tri->SetBranchAddress("iptri_kinfit", iptri_kinfit);

    TEfficiency *efficiency;

    Float_t dt_max = 90.0;
    Float_t dt_min = -dt_max;

    UInt_t nbins = floor((dt_max - dt_min)/2. + 1);

    efficiency = new TEfficiency("efficiency", "", nbins, dt_min, dt_max);

    Float_t Dtboostlor;
    UChar_t mctruth;
    Float_t Chi2, minv4gam, Kchrec[9], Qmiss, ip[3], Kchboost[9], phi_mom[4], Xcl[50], Ycl[50], Zcl[50], Tcl[50];

    chain->SetBranchAddress("mctruth", &mctruth);
    chain->SetBranchAddress("minv4gam", &minv4gam);
    chain->SetBranchAddress("Chi2", &Chi2);
    chain->SetBranchAddress("Kchboost", Kchboost);
    chain->SetBranchAddress("Qmiss", &Qmiss);

    chain->SetBranchAddress("Xcl", &Xcl);
    chain->SetBranchAddress("Ycl", &Ycl);
    chain->SetBranchAddress("Zcl", &Zcl);
    chain->SetBranchAddress("Tcl", &Tcl);

    chain->SetBranchAddress("ip", ip);

    chain->SetBranchAddress("Bpx", &phi_mom[0]);
    chain->SetBranchAddress("Bpy", &phi_mom[1]);
    chain->SetBranchAddress("Bpz", &phi_mom[2]);
    chain->SetBranchAddress("Broots", &phi_mom[3]);

    chain->AddFriend(tree_tri);

    Int_t nentries = chain->GetEntries();

    Bool_t passed_cond;

    Double_t velocity_kch, velocity_kne, tch_LAB, tch_CM, tne_LAB, tne_CM, k_path00_tri, trcv_sum, TRCV[4];
    Float_t Kch_LAB[4], Kne_LAB[4], Kch_CM[4], Kne_CM[4], Kch_CMCM[4], Kne_CMCM[4], Kne_boost[3], Kch_boost[3], Phi_boost[3], Kchmom_LAB[4], Knemom_LAB[4], Kchmom_CM[4], Knemom_CM[4];

    TH1 *trcv_hist = new TH1F("trcv", "", 100, -100.0, 50.0);

    for(Int_t i = 0; i < nentries; i++)
    {
        chain->GetEntry(i);

        if((mctruth == 1 || mctruth == 2) && done4 == 1)
        {
            velocity_kch = cVel*sqrt(pow(Kchboost[0],2) + pow(Kchboost[1],2) + pow(Kchboost[2],2))/Kchboost[3];

            velocity_kne = fourKnetri[4]/fourKnetri[3];

            k_path00_tri = sqrt(pow(fourKnetri[6] - iptri_kinfit[0],2) + pow(fourKnetri[7] - iptri_kinfit[1],2) + pow(fourKnetri[8] - iptri_kinfit[2],2));
            
            tch_LAB = sqrt(pow(Kchboost[6] - ip[0],2) + pow(Kchboost[7] - ip[1],2) + pow(Kchboost[8] - ip[2],2))/velocity_kch;
            tne_LAB = fourKnetri[9];

            Kch_LAB[0] = Kchboost[6] - ip[0];
            Kch_LAB[1] = Kchboost[7] - ip[1];
            Kch_LAB[2] = Kchboost[8] - ip[2];
            Kch_LAB[3] = tch_LAB * cVel;

            Kchmom_LAB[0] = Kchboost[0];
            Kchmom_LAB[1] = Kchboost[1];
            Kchmom_LAB[2] = Kchboost[2];
            Kchmom_LAB[3] = Kchboost[3];

            Kne_LAB[0] = fourKnetri[6] - iptri_kinfit[0];
            Kne_LAB[1] = fourKnetri[7] - iptri_kinfit[1];
            Kne_LAB[2] = fourKnetri[8] - iptri_kinfit[2];
            Kne_LAB[3] = tne_LAB * cVel;

            Knemom_LAB[0] = fourKnetri[0];
            Knemom_LAB[1] = fourKnetri[1];
            Knemom_LAB[2] = fourKnetri[2];
            Knemom_LAB[3] = fourKnetri[3];

            Phi_boost[0] = -phi_mom[0] / phi_mom[3];
            Phi_boost[1] = -phi_mom[1] / phi_mom[3];
            Phi_boost[2] = -phi_mom[2] / phi_mom[3]; 

            lorentz_transf(Phi_boost, Kch_LAB, Kch_CM);
            lorentz_transf(Phi_boost, Kne_LAB, Kne_CM);
            lorentz_transf(Phi_boost, Kchmom_LAB, Kchmom_CM);
            lorentz_transf(Phi_boost, Knemom_LAB, Knemom_CM);

            Kch_boost[0] = -Kchmom_CM[0] / Kchmom_CM[3];
            Kch_boost[1] = -Kchmom_CM[1] / Kchmom_CM[3];
            Kch_boost[2] = -Kchmom_CM[2] / Kchmom_CM[3];

            Kne_boost[0] = Kchmom_CM[0] / Kchmom_CM[3];
            Kne_boost[1] = Kchmom_CM[1] / Kchmom_CM[3];
            Kne_boost[2] = Kchmom_CM[2] / Kchmom_CM[3];

            // std::cout << Kch_boost[0] << " " << Kne_boost[0] << std::endl; 
            // std::cout << Kch_boost[1] << " " << Kne_boost[1] << std::endl; 
            // std::cout << Kch_boost[2] << " " << Kne_boost[2] << std::endl << std::endl; 

            lorentz_transf(Kch_boost, Kch_CM, Kch_CMCM);
            lorentz_transf(Kch_boost, Kne_CM, Kne_CMCM);

            Dtboostlor = (Kch_CMCM[3] - Kne_CMCM[3])/(cVel * tau_S_nonCPT);

            for(Int_t i = 0; i < 4; i++) TRCV[i] = Tcl[g4takentri_kinfit[i]-1] - (sqrt(pow(Xcl[g4takentri_kinfit[i]-1] - fourKnetri[6],2) + pow(Ycl[g4takentri_kinfit[i]-1] - fourKnetri[7],2) + pow(Zcl[g4takentri_kinfit[i]-1] - fourKnetri[8],2))/cVel) - fourKnetri[9];

            trcv_sum = (TRCV[0] + TRCV[1] + TRCV[2] + TRCV[3]);

            passed_cond = (minv4gam - mK0) < 76 && abs(Kchboost[5] - mK0) < 1.2 && Qmiss < 3.75;
            efficiency->Fill(passed_cond, Dtboostlor);

            trcv_hist->Fill(trcv_sum);
        }
    }
    TCanvas *c1 = new TCanvas("c1", "", 750, 750);

    trcv_hist->Draw();
    c1->Print("trcv_hist.png");

    efficiency->SetStatisticOption(TEfficiency::kBUniform);
    efficiency->Draw();

    gPad->Update();

    efficiency->GetPaintedGraph()->Write("correction_factor");

    file_corr.Write();
    file_corr.Close();

}