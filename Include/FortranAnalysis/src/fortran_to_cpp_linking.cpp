// Author: Szymon Gamrat
// Date of last update: .05.2024

#include <iostream>

#include <TH1.h>
#include <TCanvas.h>

#include <fort_common.h> // Linking of FORTRAN common block interfcommon
#include <fort_func.h>   // Linking of FORTRAN klspm00_lib library

#include <const.h> // Linking of const.h with standard definitions for klspm00 analysis

using namespace std;

int main()
{
  clearstruct_(); // Clear of the structure

  TChain *chain = new TChain("h1");
  chain->Add(prod2root_path_v26 + "/*" + ext_root);

  chain->SetBranchAddress("Bx", &interfcommon_.Bx);
  chain->SetBranchAddress("By", &interfcommon_.By);
  chain->SetBranchAddress("Bz", &interfcommon_.Bz);
  chain->SetBranchAddress("QuaLv", interfcommon_.qualv);
  chain->SetBranchAddress("nV", &interfcommon_.nv);
  chain->SetBranchAddress("nTv", &interfcommon_.ntv);
  chain->SetBranchAddress("iV", interfcommon_.ivOld);
  chain->SetBranchAddress("CurV", interfcommon_.CurV);
  chain->SetBranchAddress("PhiV", interfcommon_.PhiV);
  chain->SetBranchAddress("CoTv", interfcommon_.CotV);
  chain->SetBranchAddress("xV", interfcommon_.xvOld);
  chain->SetBranchAddress("yV", interfcommon_.yvOld);
  chain->SetBranchAddress("zV", interfcommon_.zvOld);

  Int_t nentries = chain->GetEntries();

  std::vector<TH1*> massKaon;
  TString massKaon_name = "";
  
  for (Int_t i = 0; i < 4; i++)
  {
    massKaon_name = "massKaon_" + i;
    massKaon.push_back(new TH1D(massKaon_name, "", 100, 0.0 ,200.0));
  };

  for (Int_t i = 0; i < nentries/100.; i++)
  {
    chain->GetEntry(i);

    // Finding Ks from KSKL->pi+pi-pi+pi-
    int
      findKl = 0,
      findKs = 1,
      findClose = 0;

    int
      last_vtx = 0;

    find_kchrec_(&findKs, &findKl, &last_vtx, &findClose, &interfcommon_.Bx, &interfcommon_.By, &interfcommon_.Bz, interfcommon_.qualv, &interfcommon_.nv, &interfcommon_.ntv, interfcommon_.ivOld, interfcommon_.CurV, interfcommon_.PhiV, interfcommon_.CotV, interfcommon_.xvOld, interfcommon_.yvOld, interfcommon_.zvOld, interfcommon_.vtakenks, interfcommon_.KchRecKS, interfcommon_.trk1KS, interfcommon_.trk2KS, &interfcommon_.cosTrkKS);

    // Finding Kl from KSKL->pi+pi-pi+pi-
    findKs = 0;
    findKl = 1;

    last_vtx = interfcommon_.vtakenks[0];

     cos_pipi_kaoncm_(&interfcommon_.Bpx, &interfcommon_.Bpy, &interfcommon_.Bpz, &interfcommon_.Broots, interfcommon_.ip_plane, interfcommon_.KchRec, interfcommon_.trk1, interfcommon_.trk2, &interfcommon_.cosTrkCM);

    find_kchrec_(&findKs, &findKl, &last_vtx, &findClose, &interfcommon_.Bx, &interfcommon_.By, &interfcommon_.Bz, interfcommon_.qualv, &interfcommon_.nv, &interfcommon_.ntv, interfcommon_.ivOld, interfcommon_.CurV, interfcommon_.PhiV, interfcommon_.CotV, interfcommon_.xvOld, interfcommon_.yvOld, interfcommon_.zvOld, interfcommon_.vtakenkl, interfcommon_.KchRecKL, interfcommon_.trk1KL, interfcommon_.trk2KL, &interfcommon_.cosTrkKL);

    // Finding Kch from KSKL->pi+pi-pi0pi0
    findKs = 0;
    findKl = 0;

    find_kchrec_(&findKs, &findKl, &last_vtx, &findClose, &interfcommon_.Bx, &interfcommon_.By, &interfcommon_.Bz, interfcommon_.qualv, &interfcommon_.nv, &interfcommon_.ntv, interfcommon_.ivOld, interfcommon_.CurV, interfcommon_.PhiV, interfcommon_.CotV, interfcommon_.xvOld, interfcommon_.yvOld, interfcommon_.zvOld, interfcommon_.vtaken, interfcommon_.KchRec, interfcommon_.trk1, interfcommon_.trk2, &interfcommon_.cosTrk);

    // Finding vtx closest to IP without any other hypothesis
    findClose = 1;

    find_kchrec_(&findKs, &findKl, &last_vtx, &findClose, &interfcommon_.Bx, &interfcommon_.By, &interfcommon_.Bz, interfcommon_.qualv, &interfcommon_.nv, &interfcommon_.ntv, interfcommon_.ivOld, interfcommon_.CurV, interfcommon_.PhiV, interfcommon_.CotV, interfcommon_.xvOld, interfcommon_.yvOld, interfcommon_.zvOld, interfcommon_.vtakenclose, interfcommon_.KchRecClose, interfcommon_.trk1Close, interfcommon_.trk2Close, &interfcommon_.cosTrkClose);

    Float_t dist[4];

    dist[0] = sqrt(pow(interfcommon_.KchRec[6] - interfcommon_.Bx,2) + 
                   pow(interfcommon_.KchRec[7] - interfcommon_.By,2));

    dist[1] = sqrt(pow(interfcommon_.KchRecKS[6] - interfcommon_.Bx,2) + 
                   pow(interfcommon_.KchRecKS[7] - interfcommon_.By,2));
    
    dist[2] = sqrt(pow(interfcommon_.KchRecKL[6] - interfcommon_.Bx,2) + 
                   pow(interfcommon_.KchRecKL[7] - interfcommon_.By,2));

    dist[3] = sqrt(pow(interfcommon_.KchRecClose[6] - interfcommon_.Bx,2) + 
                   pow(interfcommon_.KchRecClose[7] - interfcommon_.By,2));

    massKaon[0]->Fill(dist[0]);
    massKaon[1]->Fill(dist[1]);
    massKaon[2]->Fill(dist[2]);
    massKaon[3]->Fill(dist[3]);
  }

  TCanvas *c1 = new TCanvas("canva", "", 790, 790);

  c1->SetLogy(1);

  massKaon[0]->SetLineColor(kBlack);
  massKaon[1]->SetLineColor(kRed);
  massKaon[2]->SetLineColor(kBlue);
  massKaon[3]->SetLineColor(kGreen);

  massKaon[1]->GetYaxis()->SetRangeUser(10,1E6);

  massKaon[1]->Draw();
  massKaon[0]->Draw("SAMES");
  massKaon[2]->Draw("SAMES");
  massKaon[3]->Draw("SAMES");

  c1->Print("massKaon.png");

  return 0;
}