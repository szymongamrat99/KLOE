// Author: Szymon Gamrat
// Date of last update: 07.12.2024

// Check my .bashrc to find out, which paths are defined in the standard library path

#include <iostream>
#include <stdlib.h>
#include <fstream>

#include <json.hpp> // Nlohmann library for json: https://json.nlohmann.me/

#include <TChain.h>
#include <TTree.h>
#include <TBranch.h>
#include <TFile.h>
#include <TMath.h>
#include <TH1.h>
#include <TCanvas.h>

#include <const.h>

#include <fort_func.h>   // Library of functions from fortran linked to C++
#include <fort_common.h> // Common block from fortran linked to C++

using json = nlohmann::json;

int main(int argc, const char *argv[])
{
  // The paths are defined in .bashrc file, can be changed to literal
  std::string
      kloedataPath = getenv("KLOE_DBV26_DK0"),
      workdirPath = getenv("WORKDIR"),
      chainFiles = kloedataPath + "/*.root",
      pdgConstFilePath = workdirPath + "/scripts/Scripts/Properties/pdg_api/pdg_const.json";

  // Chain opened for all files in the directory
  TChain *chain = new TChain("h1");
  chain_init(chain, 40000, 40003);

  chain->SetBranchAddress("Bx", &interfcommon_.Bx);
  chain->SetBranchAddress("By", &interfcommon_.By);
  chain->SetBranchAddress("Bz", &interfcommon_.Bz);
  chain->SetBranchAddress("BPx", &interfcommon_.Bpx);
  chain->SetBranchAddress("BPy", &interfcommon_.Bpy);
  chain->SetBranchAddress("BPz", &interfcommon_.Bpz);
  chain->SetBranchAddress("Broots", &interfcommon_.Broots);
  chain->SetBranchAddress("nV", &interfcommon_.nv);
  chain->SetBranchAddress("nClu", &interfcommon_.nclu);
  chain->SetBranchAddress("nTv", &interfcommon_.ntv);
  chain->SetBranchAddress("iV", interfcommon_.iv);
  chain->SetBranchAddress("CurV", interfcommon_.CurV);
  chain->SetBranchAddress("PhiV", interfcommon_.PhiV);
  chain->SetBranchAddress("CoTv", interfcommon_.CotV);
  chain->SetBranchAddress("xV", interfcommon_.xv);
  chain->SetBranchAddress("yV", interfcommon_.yv);
  chain->SetBranchAddress("zV", interfcommon_.zv);
  chain->SetBranchAddress("EneCl", interfcommon_.EneCl);

  TFile
      *file_mctruth,
      *file_omega,
      *file_triangle;

  TTree
      *tree_mctruth,
      *tree_omega,
      *tree_triangle;

  TString
      mctruth_name = gen_vars_dir + root_files_dir + mctruth_filename + "40000_40003" + ext_root;

  file_mctruth = new TFile(mctruth_name);
  tree_mctruth = (TTree *)file_mctruth->Get(gen_vars_tree);

  tree_mctruth->SetBranchAddress("mctruth", &interfcommon_.mctruth);

  chain->AddFriend(tree_mctruth);

  // Creation of a new
  TFile *file = new TFile("invmass_test_pdg.root", "RECREATE");
  TTree *tree = new TTree("h1", "Test of pdg const file");

  Float_t
      invMass = 0.;

  TBranch *b = tree->Branch("invMass", &invMass, "invMass/F");

  // Parsing the JSON file with PDG values updated with API
  // std::ifstream PDGConst(pdgConstFilePath);
  // json constants = json::parse(PDGConst);

  const Double_t massKaonConst = 497.611; // constants["values"]["/S011M"];

  std::cout << "Mass of neutral Kaon: " << massKaonConst << std::endl;

  Int_t nentries = chain->GetEntries();

  std::vector<TH1 *> massKaon;
  TString massKaon_name = "";

  for (Int_t i = 0; i < 4; i++)
  {
    massKaon_name = "massKaon_" + std::to_string(i);
    massKaon.push_back(new TH1D(massKaon_name, "", 50, -1.0, 400.0));
  };

  for (Int_t i = 0; i < nentries; i++)
  {
    chain->GetEntry(i);

    // Finding Ks from KSKL->pi+pi-pi+pi-
    int
        findKl = 0,
        findKs = 1,
        findClose = 0;

    int
        last_vtx = 0;

    if (interfcommon_.nclu >= 4 && interfcommon_.nv > 1)
    {
      find_kchrec_(&findKs, &findKl, &last_vtx, &findClose, &interfcommon_.Bx, &interfcommon_.By, &interfcommon_.Bz, interfcommon_.qualv, &interfcommon_.nv, &interfcommon_.ntv, interfcommon_.iv, interfcommon_.CurV, interfcommon_.PhiV, interfcommon_.CotV, interfcommon_.xv, interfcommon_.yv, interfcommon_.zv, interfcommon_.vtakenks, interfcommon_.KchRecKS, interfcommon_.trk1KS, interfcommon_.trk2KS, &interfcommon_.cosTrkKS);

      // Finding Kl from KSKL->pi+pi-pi+pi-
      findKs = 0;
      findKl = 1;

      last_vtx = interfcommon_.vtakenks[0];

      find_kchrec_(&findKs, &findKl, &last_vtx, &findClose, &interfcommon_.Bx, &interfcommon_.By, &interfcommon_.Bz, interfcommon_.qualv, &interfcommon_.nv, &interfcommon_.ntv, interfcommon_.iv, interfcommon_.CurV, interfcommon_.PhiV, interfcommon_.CotV, interfcommon_.xv, interfcommon_.yv, interfcommon_.zv, interfcommon_.vtakenkl, interfcommon_.KchRecKL, interfcommon_.trk1KL, interfcommon_.trk2KL, &interfcommon_.cosTrkKL);

      // Finding Kch from KSKL->pi+pi-pi0pi0
      findKs = 0;
      findKl = 0;

      find_kchrec_(&findKs, &findKl, &last_vtx, &findClose, &interfcommon_.Bx, &interfcommon_.By, &interfcommon_.Bz, interfcommon_.qualv, &interfcommon_.nv, &interfcommon_.ntv, interfcommon_.iv, interfcommon_.CurV, interfcommon_.PhiV, interfcommon_.CotV, interfcommon_.xv, interfcommon_.yv, interfcommon_.zv, interfcommon_.vtaken, interfcommon_.KchRec, interfcommon_.trk1, interfcommon_.trk2, &interfcommon_.cosTrk);

      cor_ip_boost_(interfcommon_.KchRec, &interfcommon_.Bpx, &interfcommon_.Bpy, &interfcommon_.Bpz, &interfcommon_.Bx, &interfcommon_.By, &interfcommon_.Bz, &interfcommon_.Broots, interfcommon_.trk1, interfcommon_.trk2, interfcommon_.KchBoost, interfcommon_.ip_closest, interfcommon_.ip_plane, &interfcommon_.chdist, &interfcommon_.Qmiss);

      cos_pipi_kaoncm_(&interfcommon_.Bpx, &interfcommon_.Bpy, &interfcommon_.Bpz, &interfcommon_.Broots, interfcommon_.ip_plane, interfcommon_.KchRec, interfcommon_.trk1, interfcommon_.trk2, &interfcommon_.cosTrkCM);

      // Finding vtx closest to IP without any other hypothesis
      findClose = 1;

      find_kchrec_(&findKs, &findKl, &last_vtx, &findClose, &interfcommon_.Bx, &interfcommon_.By, &interfcommon_.Bz, interfcommon_.qualv, &interfcommon_.nv, &interfcommon_.ntv, interfcommon_.iv, interfcommon_.CurV, interfcommon_.PhiV, interfcommon_.CotV, interfcommon_.xv, interfcommon_.yv, interfcommon_.zv, interfcommon_.vtakenclose, interfcommon_.KchRecClose, interfcommon_.trk1Close, interfcommon_.trk2Close, &interfcommon_.cosTrkClose);

      // Kch from KSKL->pi+pi-pi0pi0 is saved to tree
      if (interfcommon_.KchRec[5] < 485.)
      {
        invMass = interfcommon_.KchRec[5];
        // Comparison of reconstructions to be shown on the histogram
        if(interfcommon_.mctruth == 6)
          massKaon[0]->Fill(interfcommon_.Qmiss);

        if(interfcommon_.mctruth == 5) 
          massKaon[1]->Fill(interfcommon_.Qmiss);

        if(interfcommon_.mctruth == 3)
          massKaon[2]->Fill(interfcommon_.Qmiss);

        if(interfcommon_.mctruth == 7)
          massKaon[3]->Fill(interfcommon_.Qmiss);
      }
    }

    tree->Fill();
  }

  TCanvas *c1 = new TCanvas("canva", "", 790, 790);

  massKaon[0]->SetLineColor(kBlack);
  massKaon[1]->SetLineColor(kRed);
  massKaon[2]->SetLineColor(kBlue);
  massKaon[3]->SetLineColor(kGreen);

   massKaon[3]->Draw();
  massKaon[0]->Draw("SAMES");
  massKaon[1]->Draw("SAMES");
  massKaon[2]->Draw("SAMES");

  massKaon[3]->GetYaxis()->SetRangeUser(0.0, 1.5*massKaon[3]->GetMaximum());

  c1->Print("massKaon.svg");

  tree->Print();

  file->Write();
  file->Close();

  return 0;
}