#include <iostream>

#include <TMath.h>
#include <TLorentzVector.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLatex.h>

#include "chain_init.C"
#include "../../Include/const.h"
#include "../../Include/Codes/kloe_class.h"
#include "../../Include/Codes/interf_function.h"
#include "../../Include/Codes/chi2_dist.h"

int comp_of_methods()
{
  TChain *chain = new TChain("INTERF/h1");
  chain_init(chain, 1, 56);

  TFile *file = new TFile("neuvtx_tri_kin_fit_1_56.root");
  TTree *tree = (TTree *)file->Get("h_tri_kin_fit");

  Float_t Kchboost[9], Knereclor[9], Knerec[9],
      Kchmc[9], Knemc[9], ip[3], ipmc[3], phi_mom[4], Dtmc;
  UChar_t mctruth;

  chain->SetBranchAddress("Kchboost", Kchboost);
  chain->SetBranchAddress("Knereclor", Knereclor);
  chain->SetBranchAddress("Knerec", Knerec);
  chain->SetBranchAddress("Kchmc", Kchmc);
  chain->SetBranchAddress("Knemc", Knemc);

  chain->SetBranchAddress("ip", ip);
  chain->SetBranchAddress("ipmc", ipmc);

  chain->SetBranchAddress("Bpx", &phi_mom[0]);
  chain->SetBranchAddress("Bpy", &phi_mom[1]);
  chain->SetBranchAddress("Bpz", &phi_mom[2]);
  chain->SetBranchAddress("Broots", &phi_mom[3]);

  chain->SetBranchAddress("mctruth", &mctruth);

  chain->SetBranchAddress("Dtmc", &Dtmc);

  Float_t gamma_kinfit[4][8], ip_kinfit[3], Knetri_kinfit[10], chi2min;
  Int_t done_kinfit;

  tree->SetBranchAddress("fourgamma1tri_kinfit", gamma_kinfit[0]);
  tree->SetBranchAddress("fourgamma2tri_kinfit", gamma_kinfit[1]);
  tree->SetBranchAddress("fourgamma3tri_kinfit", gamma_kinfit[2]);
  tree->SetBranchAddress("fourgamma4tri_kinfit", gamma_kinfit[3]);

  tree->SetBranchAddress("fourKnetri_kinfit", Knetri_kinfit);

  tree->SetBranchAddress("iptri_kinfit", ip_kinfit);
  tree->SetBranchAddress("done4_kinfit", &done_kinfit);

  tree->SetBranchAddress("chi2min", &chi2min);

  chain->AddFriend(tree);

  UInt_t nentries = (UInt_t)chain->GetEntries();

  Float_t lengthneu_mc, lengthneu_tri, lengthneu_rec, lengthch_mc, lengthch_rec,
          v_Kneutri, v_Kchrec, v_Kneurec, v_Kneumc, v_Kchmc,
          t_chmc, t_neumc, t_chrec, t_neutri, t_neurec;

  TString id_hist;
  TH1 *dttri_hist;
  TH2 *neu_vtx_corr[4];

  for(Int_t i = 0; i < 3; i++)
  {
    id_hist = "Coor" + std::to_string(i);
    neu_vtx_corr[i] = new TH2F(id_hist, "", 100, -50, 50, 100, -50, 50);
  }

  neu_vtx_corr[3] = new TH2F("Lengths", "", 201, 0, 100, 201, 0, 100);

  dttri_hist = new TH1F("Dttri", "", 100, 0, 50);


  for (Int_t i = 0; i < nentries; i++)
  {
    chain->GetEntry(i);

    if (done_kinfit == 1)
    {
      neu_vtx_corr[0]->Fill(Knemc[6], Knetri_kinfit[6]);
      neu_vtx_corr[1]->Fill(Knemc[7], Knetri_kinfit[7]);
      neu_vtx_corr[2]->Fill(Knemc[8], Knetri_kinfit[8]);

      lengthneu_mc = sqrt(pow(Knemc[6] - ipmc[0],2) + pow(Knemc[7] - ipmc[1],2) + pow(Knemc[8] - ipmc[2],2));
      lengthneu_tri = sqrt(pow(Knetri_kinfit[6] - ip_kinfit[0],2) + pow(Knetri_kinfit[7] - ip_kinfit[1],2) + pow(Knetri_kinfit[8] - ip_kinfit[2],2));
      lengthneu_rec = sqrt(pow(Knereclor[6] - ip[0],2) + pow(Knereclor[7] - ip[1],2) + pow(Knereclor[8] - ip[2],2));

      lengthch_mc = sqrt(pow(Kchmc[6] - ipmc[0],2) + pow(Kchmc[7] - ipmc[1],2) + pow(Kchmc[8] - ipmc[2],2));
      lengthch_rec = sqrt(pow(Kchboost[6] - ip[0],2) + pow(Kchboost[7] - ip[1],2) + pow(Kchboost[8] - ip[2],2));

      v_Kneumc = c_vel * Knemc[4]/Knemc[3];
      v_Kneutri = c_vel * Knetri_kinfit[4]/Knetri_kinfit[3];
      v_Kneurec = c_vel * Knereclor[5]/Knereclor[3];

      v_Kchmc = c_vel * Kchmc[4]/Kchmc[3];
      v_Kchrec = c_vel * Kchboost[4]/Kchboost[3];

      t_chrec = lengthch_rec / v_Kchrec;
      t_chmc = lengthch_mc / v_Kchmc;

      t_neurec = lengthneu_rec / v_Kneurec;
      t_neutri = lengthneu_tri / v_Kneutri;
      t_neumc = lengthneu_mc / v_Kneumc;

      neu_vtx_corr[3]->Fill(t_chrec/tau_S_nonCPT, t_neutri/tau_S_nonCPT, interf_function((t_chmc - t_neumc)/tau_S_nonCPT));

      dttri_hist->Fill(chi2min);

    }
  }

  gStyle->SetOptStat(0);

  TCanvas *c[5];
  TString id_canva;
  Int_t width = 750, height = 750;

  for(Int_t i = 0; i < 4; i++)
  {
    id_canva = "Canva" + std::to_string(i);
    c[i] = new TCanvas(id_canva, "", width, height);
    c[i]->SetRightMargin(0.15);
    c[i]->cd();

    neu_vtx_corr[i]->GetXaxis()->SetTitle("V^{neu}_{gen} [cm]");
    neu_vtx_corr[i]->GetYaxis()->SetTitle("V^{neu}_{tri} [cm]");
    neu_vtx_corr[i]->Draw("COLZ");

    c[i]->Print(id_canva + ".png");
  }

  TF1 *func = new TF1("func",chi2dist,0,50,2);
  func->SetParameters(dttri_hist->Integral(0,100)/10.,6);
  func->SetParNames("Normalization","Degrees of freedom");

  id_canva = "Canva" + std::to_string(4);
  c[4] = new TCanvas(id_canva, "", width, height);
  c[4]->SetRightMargin(0.15);
  c[4]->cd();

  //dttri_hist->Fit(func);
  dttri_hist->Draw();
  func->Draw("SAMES");
  c[4]->Print(id_canva + ".png");




  return 0;
}