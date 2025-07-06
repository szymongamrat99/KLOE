#include <iostream>
#include <fstream>

#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TH2.h>
#include <TFitResult.h>
#include <TGraphAsymmErrors.h>

#include <fort_common.h>
#include <kloe_class.h>
#include <const.h>

#include "../inc/regenrejec.hpp"

using namespace std;

int plots(Int_t firstFile, Int_t lastFile)
{
  TChain *chain = new TChain("INTERF/h1");
  chain_init(chain, firstFile, lastFile);

  // Tree names and filenames to be analyzed
  std::map<TString, TString>
      treeNames,
      fileNames;

  std::vector<std::string>
      propertyName;

  propertyName.push_back("trianglefinal");
  propertyName.push_back("mctruth");

  treeNames[propertyName[0]] = (std::string)properties["variables"]["tree"]["treename"][propertyName[0]];
  fileNames[propertyName[0]] = (std::string)properties["variables"]["tree"]["filename"][propertyName[0]];
  treeNames[propertyName[1]] = (std::string)properties["variables"]["tree"]["treename"][propertyName[1]];
  fileNames[propertyName[1]] = (std::string)properties["variables"]["tree"]["filename"][propertyName[1]];

  //

  TFile file_tri(fileNames[propertyName[0]]);
  TTree *tree_tri = (TTree *)file_tri.Get(treeNames[propertyName[0]]);

  TFile file_mctruth(fileNames[propertyName[1]]);
  TTree *tree_mctruth = (TTree *)file_mctruth.Get(treeNames[propertyName[1]]);

  chain->SetBranchAddress("mcflag", &interfcommon_.mcflag);
  chain->SetBranchAddress("Kchboost", interfcommon_.KchBoost);
  chain->SetBranchAddress("Knereclor", interfcommon_.KneRecLor);
  chain->SetBranchAddress("ip", interfcommon_.ip);

  tree_tri->SetBranchAddress("fourKnetriangle", interfcommon_.KneRecTriangle);
  tree_tri->SetBranchAddress("done_triangle", &interfcommon_.done4);

  tree_mctruth->SetBranchAddress("mctruth", &interfcommon_.mctruth);

  chain->AddFriend(tree_tri);
  chain->AddFriend(tree_mctruth);

  Int_t nentries = chain->GetEntries();

  map<Int_t, TString> region_dict, method_dict;

  region_dict[0] = "R_blue";
  region_dict[1] = "#rho_blue";
  region_dict[2] = "R_green";
  region_dict[3] = "#rho_green";

  method_dict[0] = "z-coor_section";
  method_dict[1] = "#rho_section";

  Double_t
      width_rho_triangle = properties["variables"]["Resolutions"]["rhoNeutral"]["triTriangle"],
      width_R_triangle = properties["variables"]["Resolutions"]["pathNeutral"]["triTriangle"],
      width_rho_pure_triangle = properties["variables"]["Resolutions"]["rhoNeutral"]["triTriangle"],
      width_R_pure_triangle = properties["variables"]["Resolutions"]["pathNeutral"]["triTriangle"],
      width_rho_charged = properties["variables"]["Resolutions"]["rhoCharged"],
      width_R_charged = properties["variables"]["Resolutions"]["pathCharged"];

  Double_t
      max_path = 50.0,
      min_path = 0.0;

  Int_t
      bins_rho_pure_triangle = Int_t((max_path - min_path) / width_rho_pure_triangle),
      bins_path_pure_triangle = Int_t((max_path - min_path) / width_R_pure_triangle),
      bins_rho_triangle = Int_t((max_path - min_path) / width_rho_triangle),
      bins_path_triangle = Int_t((max_path - min_path) / width_R_triangle),
      bins_rho_charged = Int_t((max_path - min_path) / width_rho_charged),
      bins_path_charged = Int_t((max_path - min_path) / width_R_charged);

  std::vector<TCanvas *> deltaTcanva;
  TString deltaTcanva_name = "";

  for (Int_t i = 0; i < 7; i++)
  {
    deltaTcanva_name = "deltaTcanva_" + std::to_string(i);
    deltaTcanva.push_back(new TCanvas(deltaTcanva_name, deltaTcanva_name));
  };

  std::vector<TH1 *> deltaT;
  TString deltaT_name = "";

  Double_t
      dt_min = properties["variables"]["CPFit"]["histoResults"]["rangeX"][0],
      dt_max = properties["variables"]["CPFit"]["histoResults"]["rangeX"][1],
      dtWidth = properties["variables"]["Resolutions"]["deltaT"];

  Int_t bin_num = floor((dt_max - dt_min) / dtWidth) + 1;

  for (Int_t i = 0; i < 7; i++)
  {
    deltaT_name = "delta_" + std::to_string(i);
    deltaT.push_back(new TH1D(deltaT_name, "", bin_num, dt_min, dt_max));
  };

  std::vector<std::string>
      method;

  method.push_back("methodB");
  method.push_back("methodA");

  Double_t
      methodA_lower = properties["variables"]["RegenRejection"]["boundaries"][method[1]][0],
      methodA_higher = properties["variables"]["RegenRejection"]["boundaries"][method[1]][1],
      methodB_bound = properties["variables"]["RegenRejection"]["boundaries"][method[0]][0],
      sigma = properties["variables"]["RegenRejection"]["sigma"];

  Long64_t count_tot = 0, count_neg_reg = 0, count_pos_reg = 0, count_sig = 0, count_sig_neg = 0, count_sig_pos = 0;

  KLOE::pm00 auxilliary;

  Double_t
        meanRadiusChHigher = properties["variables"]["RegenRejection"]["results"]["methodA"]["charged"]["spherical"]["mean"][1],
        errorRadiusChHigher = properties["variables"]["RegenRejection"]["results"]["methodA"]["charged"]["spherical"]["width"][1],
        meanRadiusChLower = properties["variables"]["RegenRejection"]["results"]["methodA"]["charged"]["cylindrical"]["mean"][0],
        errorRadiusChLower = properties["variables"]["RegenRejection"]["results"]["methodA"]["charged"]["cylindrical"]["width"][0],
        
        meanRadiusTriangleHigher = properties["variables"]["RegenRejection"]["results"]["methodA"]["triangle"]["spherical"]["mean"][1],
        errorRadiusTriangleHigher = properties["variables"]["RegenRejection"]["results"]["methodA"]["triangle"]["spherical"]["width"][1],
        meanRadiusTriangleLower = properties["variables"]["RegenRejection"]["results"]["methodA"]["triangle"]["cylindrical"]["mean"][0],
        errorRadiusTriangleLower = properties["variables"]["RegenRejection"]["results"]["methodA"]["triangle"]["cylindrical"]["width"][0],
        
        meanRadiusPureTriangleHigher = properties["variables"]["RegenRejection"]["results"]["methodA"]["pureTriangle"]["spherical"]["mean"][1],
        errorRadiusPureTriangleHigher = properties["variables"]["RegenRejection"]["results"]["methodA"]["pureTriangle"]["spherical"]["width"][1],
        meanRadiusPureTriangleLower = properties["variables"]["RegenRejection"]["results"]["methodA"]["pureTriangle"]["cylindrical"]["mean"][0],
        errorRadiusPureTriangleLower = properties["variables"]["RegenRejection"]["results"]["methodA"]["pureTriangle"]["cylindrical"]["width"][0];

  for (Int_t i = 0; i < nentries; i++)
  {
    chain->GetEntry(i);

    Bool_t cuts[6] = {false, false, false, false, false, false};
    Double_t radius[2] = {0., 0.}, radius_ch[2] = {0., 0.}, radius_tri[2] = {0., 0.};

    for (Int_t i = 0; i < 3; i++)
    {
      // Spherical
      radius[0] += pow(interfcommon_.KneRecLor[6 + i] - interfcommon_.ip[i], 2);
      radius_ch[0] += pow(interfcommon_.KchBoost[6 + i] - interfcommon_.ip[i], 2);
      radius_tri[0] += pow(interfcommon_.KneRecTriangle[6 + i] - interfcommon_.ip[i], 2);

      if (i < 2)
      {
        // Cylindrical
        radius[1] += pow(interfcommon_.KneRecLor[6 + i] - interfcommon_.ip[i], 2);
        radius_ch[1] += pow(interfcommon_.KchBoost[6 + i] - interfcommon_.ip[i], 2);
        radius_tri[1] += pow(interfcommon_.KneRecTriangle[6 + i] - interfcommon_.ip[i], 2);
      }
    }

    radius[0] = sqrt(radius[0]);
    radius[1] = sqrt(radius[1]);

    radius_ch[0] = sqrt(radius_ch[0]);
    radius_ch[1] = sqrt(radius_ch[1]);

    radius_tri[0] = sqrt(radius_tri[0]);
    radius_tri[1] = sqrt(radius_tri[1]);

    count_tot++;

    if (interfcommon_.mctruth == 3)
    {
      // Neutral vtx triangle method

      TLorentzVector *momKch = new TLorentzVector(interfcommon_.KchBoost[0], interfcommon_.KchBoost[1], interfcommon_.KchBoost[2], interfcommon_.KchBoost[3]);

      TVector3 boostCh = momKch->BoostVector() * cVel;
      Double_t tCh = radius_ch[0] / (boostCh.Mag());

      TLorentzVector *posKch = new TLorentzVector(interfcommon_.KchBoost[6] - interfcommon_.ip[0], interfcommon_.KchBoost[7] - interfcommon_.ip[1], interfcommon_.KchBoost[8] - interfcommon_.ip[2], tCh * cVel);

      TLorentzVector *momKne = new TLorentzVector(interfcommon_.KneRecTriangle[0], interfcommon_.KneRecTriangle[1], interfcommon_.KneRecTriangle[2], interfcommon_.KneRecTriangle[3]);

      TVector3 boostNeu = momKne->BoostVector() * cVel;
      Double_t tNeu = radius_tri[0] / (boostNeu.Mag());

      TLorentzVector *posKne = new TLorentzVector(interfcommon_.KneRecTriangle[6] - interfcommon_.ip[0], interfcommon_.KneRecTriangle[7] - interfcommon_.ip[1], interfcommon_.KneRecTriangle[8] - interfcommon_.ip[2], tNeu * cVel);

      interfcommon_.DtBoostRec = auxilliary.DeltaT(momKch, posKch, momKne, posKne);

      cuts[0] = abs(radius[0] - meanRadiusPureTriangleHigher) > errorRadiusPureTriangleHigher * sigma;      // && radius[1] > 8;
      cuts[1] = abs(radius[1] - meanRadiusPureTriangleLower) > errorRadiusPureTriangleLower * sigma;     // && radius[1] <= 8;
      cuts[2] = abs(radius_ch[0] - meanRadiusChHigher) > errorRadiusChHigher * sigma; // && radius_ch[1] > 8;
      cuts[3] = abs(radius_ch[1] - meanRadiusChLower) > errorRadiusChLower * sigma;  // && radius_ch[1] <= 8;
      cuts[4] = abs(radius_tri[0] - meanRadiusTriangleHigher) > errorRadiusTriangleHigher * sigma; // && radius_tri[1] > 8;
      cuts[5] = abs(radius_tri[1] - meanRadiusTriangleLower) > errorRadiusTriangleLower * sigma;         // && radius_tri[1] <= 8;

      deltaT[0]->Fill(interfcommon_.DtBoostRec);
      if (cuts[0])
        deltaT[1]->Fill(interfcommon_.DtBoostRec);

      if (cuts[1])
        deltaT[2]->Fill(interfcommon_.DtBoostRec);

      if (cuts[2])
        deltaT[3]->Fill(interfcommon_.DtBoostRec);

      if (cuts[3])
        deltaT[4]->Fill(interfcommon_.DtBoostRec);

      if (cuts[4])
        deltaT[5]->Fill(interfcommon_.DtBoostRec);

      if (cuts[5])
        deltaT[6]->Fill(interfcommon_.DtBoostRec);
    }
  }

  TString img_name;

  for (Int_t i = 1; i < 7; i++)
  {
    img_name = SystemPath::img_dir + "RegenerationAnalysis/deltaT" + i + ext_img;
    deltaTcanva[i - 1]->cd();

    deltaT[0]->SetLineColor(kGreen);
    deltaT[0]->GetYaxis()->SetMaxDigits(3);

    deltaT[0]->Draw();
    deltaT[i]->Draw("SAME");

    deltaTcanva[i - 1]->Print(img_name);
  }

  return 0;
}