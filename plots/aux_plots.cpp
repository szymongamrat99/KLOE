#include <iostream>

#include <TH1.h>
#include <THStack.h>
#include <TCanvas.h>

#include "../../Include/const.h"
#include "../../Include/Codes/chain_init.cpp"
#include "../../Include/Codes/kloe_class.h"

using namespace std;

int main()
{
  TChain *chain = new TChain("INTERF/h1");
  chain_init(chain, 1, 56);

  KLOE::pm00 ev;

  Float_t Kchrec[9] = {0.}, Kchboost[9] = {0.}, ipBhabha[3] = {0.}, pBhabha[4] = {0.}, ip[3] = {0.}, Knerec[9] = {0.}, Knereclor[9] = {0.};
  UChar_t mctruth = 0;

  chain->SetBranchAddress("Kchrec", Kchrec);
  chain->SetBranchAddress("Kchboost", Kchboost);

  chain->SetBranchAddress("Knerec", Knerec);
  chain->SetBranchAddress("Knereclor", Knereclor);

  chain->SetBranchAddress("Bx", &ipBhabha[0]);
  chain->SetBranchAddress("By", &ipBhabha[1]);
  chain->SetBranchAddress("Bz", &ipBhabha[2]);

  chain->SetBranchAddress("Bpx", &pBhabha[0]);
  chain->SetBranchAddress("Bpy", &pBhabha[1]);
  chain->SetBranchAddress("Bpz", &pBhabha[2]);
  chain->SetBranchAddress("Broots", &pBhabha[3]);

  chain->SetBranchAddress("ip", ip);

  chain->SetBranchAddress("mctruth", &mctruth);

  Bool_t sig_only = false;

  Int_t nentries = chain->GetEntries();
  const size_t ch_chann = 3;

  TString hist_name = "";
  array<TH1 *, ch_chann> kch_hist, kne_hist;

  for (Int_t i = 0; i < kch_hist.max_size(); i++)
  {
    hist_name = "Kch_inv_mass_" + to_string(i);
    kch_hist[i] = new TH1D(hist_name, hist_name, 50, mK0 - 10., mK0 + 10.);
  }

  for (Int_t i = 0; i < kne_hist.max_size(); i++)
  {
    hist_name = "Kne_inv_mass_" + to_string(i);
    kne_hist[i] = new TH1D(hist_name, hist_name, 50, mK0 - 100., mK0 + 100.);
  }

  for (Int_t i = 0; i < nentries; i++)
  {
    chain->GetEntry(i);

    sig_only = (mctruth != 0);

    if (sig_only)
    {
      for (Int_t i = 0; i < kch_hist.max_size(); i++)
        kch_hist[i]->Fill(Kchrec[5]);

      for (Int_t i = 0; i < kne_hist.max_size(); i++)
        kne_hist[i]->Fill(Knerec[5]);
    }
  }

  TCanvas *c1 = new TCanvas("c1","c1",790,790);
  TCanvas *c2 = new TCanvas("c2","c2",790,790);

  THStack *kch_all_hist = new THStack("Kch_all", "Kch_all;m_{K_{ch},inv} [MeV];Counts");
  THStack *kne_all_hist = new THStack("Kne_all", "Kne_all;m_{K_{ne},inv} [MeV];Counts");
  Color_t color[ch_chann] = {kRed, kBlue, kGreen};

  c1->cd();
  for (Int_t i = 0; i < kch_hist.max_size(); i++)
  {
    kch_hist[i]->SetLineColor(color[i]);
    kch_hist[i]->SetLineWidth(3);
    kch_all_hist->Add(kch_hist[i]);
  }
  
  kch_all_hist->Draw();
  c1->Print("test1.svg");

  c2->cd();
  for (Int_t i = 0; i < kne_hist.max_size(); i++)
  {
    kne_hist[i]->SetLineColor(color[i]);
    kne_hist[i]->SetLineWidth(3);
    kne_all_hist->Add(kne_hist[i]);
  }
  
  kne_all_hist->Draw();
  c2->Print("test2.svg");


  return 0;
}
