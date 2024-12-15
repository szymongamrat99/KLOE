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
#include <const.h>

#include "../inc/regenrejec.hpp"

using namespace std;

int regenrejec(Int_t firstFile, Int_t lastFile)
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

  TCanvas *canvas[100];
  TH1 *h_radius[2][4], *h_radius_tri[2][4], *h_radius_ch[2][4];

  TString hist_name, canvas_name;

  for (Int_t i = 0; i < 100; i++)
  {
    canvas_name = "Canva_" + to_string(i);
    canvas[i] = new TCanvas(canvas_name, canvas_name);
    canvas[i]->SetLeftMargin(0.2);
    canvas[i]->SetBottomMargin(0.2);
  }

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

  for (Int_t i = 0; i < 2; i++)
    for (Int_t j = 0; j < 4; j++)
    {
      hist_name = "Hist_radius_triangle_" + method_dict[i] + "_" + region_dict[j];
      if (j == 0 || j == 2)
        h_radius[i][j] = new TH1D(hist_name, hist_name, bins_path_pure_triangle, min_path, max_path);
      else
        h_radius[i][j] = new TH1D(hist_name, hist_name, bins_rho_pure_triangle, min_path, max_path);

      hist_name = "Hist_radius_tri_" + method_dict[i] + "_" + region_dict[j];
      if (j == 0 || j == 2)
        h_radius_tri[i][j] = new TH1D(hist_name, hist_name, bins_path_triangle, min_path, max_path);
      else
        h_radius_tri[i][j] = new TH1D(hist_name, hist_name, bins_rho_triangle, min_path, max_path);

      hist_name = "Hist_radius_ch_" + method_dict[i] + "_" + region_dict[j];
      if (j == 0 || j == 2)
        h_radius_ch[i][j] = new TH1D(hist_name, hist_name, bins_path_charged, min_path, max_path);
      else
        h_radius_ch[i][j] = new TH1D(hist_name, hist_name, bins_rho_charged, min_path, max_path);

      h_radius[i][j]->GetYaxis()->SetMaxDigits(3);
      h_radius_tri[i][j]->GetYaxis()->SetMaxDigits(3);
      h_radius_ch[i][j]->GetYaxis()->SetMaxDigits(3);
    }

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

      if (abs(interfcommon_.KneRecLor[8]) > methodB_bound)
      {
        h_radius[0][0]->Fill(radius[0]);
        h_radius[0][1]->Fill(radius[1]);
      }
      else
      {
        h_radius[0][2]->Fill(radius[0]);
        h_radius[0][3]->Fill(radius[1]);
      }

      if (radius[1] > methodA_lower && radius[1] < methodA_higher)
      {
        h_radius[1][0]->Fill(radius[0]);
        h_radius[1][1]->Fill(radius[1]);
      }
      else if (radius[1] < methodA_lower)
      {
        h_radius[1][2]->Fill(radius[0]);
        h_radius[1][3]->Fill(radius[1]);
      }

      // Charged vtx

      if (abs(interfcommon_.KchBoost[8]) > methodB_bound)
      {
        h_radius_ch[0][0]->Fill(radius_ch[0]);
        h_radius_ch[0][1]->Fill(radius_ch[1]);
      }
      else
      {
        h_radius_ch[0][2]->Fill(radius_ch[0]);
        h_radius_ch[0][3]->Fill(radius_ch[1]);
      }

      if (radius_ch[1] > methodA_lower && radius_ch[1] < methodA_higher)
      {
        h_radius_ch[1][0]->Fill(radius_ch[0]);
        h_radius_ch[1][1]->Fill(radius_ch[1]);
      }
      else if (radius_ch[1] < methodA_lower)
      {
        h_radius_ch[1][2]->Fill(radius_ch[0]);
        h_radius_ch[1][3]->Fill(radius_ch[1]);
      }

      // Neutral vtx trilateration method

      if (interfcommon_.done4 == 1)
      {
        if (abs(interfcommon_.KneRecTriangle[8]) > methodB_bound)
        {
          h_radius_tri[0][0]->Fill(radius_tri[0]);
          h_radius_tri[0][1]->Fill(radius_tri[1]);
        }
        else
        {
          h_radius_tri[0][2]->Fill(radius_tri[0]);
          h_radius_tri[0][3]->Fill(radius_tri[1]);
        }

        if (radius_tri[1] > methodA_lower && radius_tri[1] < methodA_higher)
        {
          h_radius_tri[1][0]->Fill(radius_tri[0]);
          h_radius_tri[1][1]->Fill(radius_tri[1]);
        }
        else if (radius_tri[1] < methodA_lower)
        {
          h_radius_tri[1][2]->Fill(radius_tri[0]);
          h_radius_tri[1][3]->Fill(radius_tri[1]);
        }
      }
    }
  }

  TFitResultPtr r;
  Int_t fitStatus = 0;

  TString img_name;

  Int_t canva_num = 0;

  for (Int_t i = 0; i < 2; i++)
  {
    for (Int_t j = 0; j < 4; j++)
    {
      img_name = img_dir + "RegenerationAnalysis/radius_" + method_dict[i] + "_" + region_dict[j] + ext_img;
      canvas[canva_num]->cd();
      canvas[canva_num]->SetLogy();

      if (j == 0 || j == 2)
      {
        r = h_radius[i][j]->Fit("gaus", "SM", "", 8., 13.);
        fitStatus = r;
        if (fitStatus >= 0)
        {
          if (j == 2)
          {
            properties["variables"]["RegenRejection"]["results"][method[i]]["pureTriangle"]["spherical"]["mean"][0] = r->Parameter(1);
            properties["variables"]["RegenRejection"]["results"][method[i]]["pureTriangle"]["spherical"]["width"][0] = r->Parameter(2);
            properties["variables"]["RegenRejection"]["errors"][method[i]]["pureTriangle"]["spherical"]["mean"][0] = r->Error(1);
            properties["variables"]["RegenRejection"]["errors"][method[i]]["pureTriangle"]["spherical"]["width"][0] = r->Error(2);
          }
          else
          {
            properties["variables"]["RegenRejection"]["results"][method[i]]["pureTriangle"]["spherical"]["mean"][1] = r->Parameter(1);
            properties["variables"]["RegenRejection"]["results"][method[i]]["pureTriangle"]["spherical"]["width"][1] = r->Parameter(2);
            properties["variables"]["RegenRejection"]["errors"][method[i]]["pureTriangle"]["spherical"]["mean"][1] = r->Error(1);
            properties["variables"]["RegenRejection"]["errors"][method[i]]["pureTriangle"]["spherical"]["width"][1] = r->Error(2);
          }

          cout << img_name << ": (" << r->Parameter(1) << "+-" << r->Error(1) << ")" << endl;
          cout << img_name << ": (" << r->Parameter(2) << "+-" << r->Error(2) << ")" << endl;
        }
        else
        {
          if (j == 2)
          {
            properties["variables"]["RegenRejection"]["results"][method[i]]["pureTriangle"]["spherical"]["mean"][0] = nullptr;
            properties["variables"]["RegenRejection"]["results"][method[i]]["pureTriangle"]["spherical"]["width"][0] = nullptr;
            properties["variables"]["RegenRejection"]["errors"][method[i]]["pureTriangle"]["spherical"]["mean"][0] = nullptr;
            properties["variables"]["RegenRejection"]["errors"][method[i]]["pureTriangle"]["spherical"]["width"][0] = nullptr;
          }
          else
          {
            properties["variables"]["RegenRejection"]["results"][method[i]]["pureTriangle"]["spherical"]["mean"][1] = nullptr;
            properties["variables"]["RegenRejection"]["results"][method[i]]["pureTriangle"]["spherical"]["width"][1] = nullptr;
            properties["variables"]["RegenRejection"]["errors"][method[i]]["pureTriangle"]["spherical"]["mean"][1] = nullptr;
            properties["variables"]["RegenRejection"]["errors"][method[i]]["pureTriangle"]["spherical"]["width"][1] = nullptr;
          }
        }

        h_radius[i][j]->GetXaxis()->SetTitle("R [cm]");
        h_radius[i][j]->GetYaxis()->SetTitle(Form("Counts/(%.2f cm)", width_R_pure_triangle));
      }
      else if (j == 1 || j == 3)
      {
        r = h_radius[i][j]->Fit("gaus", "SM", "", 4., 6.);
        fitStatus = r;
        if (fitStatus >= 0)
        {
          if (j == 3)
          {
            properties["variables"]["RegenRejection"]["results"][method[i]]["pureTriangle"]["cylindrical"]["mean"][0] = r->Parameter(1);
            properties["variables"]["RegenRejection"]["results"][method[i]]["pureTriangle"]["cylindrical"]["width"][0] = r->Parameter(2);
            properties["variables"]["RegenRejection"]["errors"][method[i]]["pureTriangle"]["cylindrical"]["mean"][0] = r->Error(1);
            properties["variables"]["RegenRejection"]["errors"][method[i]]["pureTriangle"]["cylindrical"]["width"][0] = r->Error(2);
          }
          else
          {
            properties["variables"]["RegenRejection"]["results"][method[i]]["pureTriangle"]["cylindrical"]["mean"][1] = r->Parameter(1);
            properties["variables"]["RegenRejection"]["results"][method[i]]["pureTriangle"]["cylindrical"]["width"][1] = r->Parameter(2);
            properties["variables"]["RegenRejection"]["errors"][method[i]]["pureTriangle"]["cylindrical"]["mean"][1] = r->Error(1);
            properties["variables"]["RegenRejection"]["errors"][method[i]]["pureTriangle"]["cylindrical"]["width"][1] = r->Error(2);
          }

          cout << img_name << ": (" << r->Parameter(1) << "+-" << r->Error(1) << ")" << endl;
          cout << img_name << ": (" << r->Parameter(2) << "+-" << r->Error(2) << ")" << endl;
        }
        else
        {
          if (j == 3)
          {
            properties["variables"]["RegenRejection"]["results"][method[i]]["pureTriangle"]["cylindrical"]["mean"][0] = nullptr;
            properties["variables"]["RegenRejection"]["results"][method[i]]["pureTriangle"]["cylindrical"]["width"][0] = nullptr;
            properties["variables"]["RegenRejection"]["errors"][method[i]]["pureTriangle"]["cylindrical"]["mean"][0] = nullptr;
            properties["variables"]["RegenRejection"]["errors"][method[i]]["pureTriangle"]["cylindrical"]["width"][0] = nullptr;
          }
          else
          {
            properties["variables"]["RegenRejection"]["results"][method[i]]["pureTriangle"]["cylindrical"]["mean"][1] = nullptr;
            properties["variables"]["RegenRejection"]["results"][method[i]]["pureTriangle"]["cylindrical"]["width"][1] = nullptr;
            properties["variables"]["RegenRejection"]["errors"][method[i]]["pureTriangle"]["cylindrical"]["mean"][1] = nullptr;
            properties["variables"]["RegenRejection"]["errors"][method[i]]["pureTriangle"]["cylindrical"]["width"][1] = nullptr;
          }
        }
      }

      h_radius[i][j]->GetXaxis()->SetTitle("#rho [cm]");
      h_radius[i][j]->GetYaxis()->SetTitle(Form("Counts/(%.2f cm)", width_rho_pure_triangle));

      h_radius[i][j]->Draw();
      canvas[canva_num]->Print(img_name);

      canva_num++;
    }
  }

  for (Int_t i = 0; i < 2; i++)
  {
    for (Int_t j = 0; j < 4; j++)
    {
      img_name = img_dir + "RegenerationAnalysis/radius_ch_" + method_dict[i] + "_" + region_dict[j] + ext_img;
      canvas[canva_num]->cd();
      canvas[canva_num]->SetLogy();

      if (j == 0 || j == 2)
      {
        r = h_radius_ch[i][j]->Fit("gaus", "SM", "", 8., 12.);
        fitStatus = r;
        if (fitStatus >= 0)
        {
          if (j == 2)
          {
            properties["variables"]["RegenRejection"]["results"][method[i]]["charged"]["spherical"]["mean"][0] = r->Parameter(1);
            properties["variables"]["RegenRejection"]["results"][method[i]]["charged"]["spherical"]["width"][0] = r->Parameter(2);
            properties["variables"]["RegenRejection"]["errors"][method[i]]["charged"]["spherical"]["mean"][0] = r->Error(1);
            properties["variables"]["RegenRejection"]["errors"][method[i]]["charged"]["spherical"]["width"][0] = r->Error(2);
          }
          else
          {
            properties["variables"]["RegenRejection"]["results"][method[i]]["charged"]["spherical"]["mean"][1] = r->Parameter(1);
            properties["variables"]["RegenRejection"]["results"][method[i]]["charged"]["spherical"]["width"][1] = r->Parameter(2);
            properties["variables"]["RegenRejection"]["errors"][method[i]]["charged"]["spherical"]["mean"][1] = r->Error(1);
            properties["variables"]["RegenRejection"]["errors"][method[i]]["charged"]["spherical"]["width"][1] = r->Error(2);
          }

          cout << img_name << ": (" << r->Parameter(1) << "+-" << r->Error(1) << ")" << endl;
          cout << img_name << ": (" << r->Parameter(2) << "+-" << r->Error(2) << ")" << endl;
        }
        else
        {
          if (j == 2)
          {
            properties["variables"]["RegenRejection"]["results"][method[i]]["charged"]["spherical"]["mean"][0] = nullptr;
            properties["variables"]["RegenRejection"]["results"][method[i]]["charged"]["spherical"]["width"][0] = nullptr;
            properties["variables"]["RegenRejection"]["errors"][method[i]]["charged"]["spherical"]["mean"][0] = nullptr;
            properties["variables"]["RegenRejection"]["errors"][method[i]]["charged"]["spherical"]["width"][0] = nullptr;
          }
          else
          {
            properties["variables"]["RegenRejection"]["results"][method[i]]["charged"]["spherical"]["mean"][1] = nullptr;
            properties["variables"]["RegenRejection"]["results"][method[i]]["charged"]["spherical"]["width"][1] = nullptr;
            properties["variables"]["RegenRejection"]["errors"][method[i]]["charged"]["spherical"]["mean"][1] = nullptr;
            properties["variables"]["RegenRejection"]["errors"][method[i]]["charged"]["spherical"]["width"][1] = nullptr;
          }
        }
        h_radius_ch[i][j]->GetXaxis()->SetTitle("R [cm]");
        h_radius_ch[i][j]->GetYaxis()->SetTitle(Form("Counts/(%.2f cm)", width_R_charged));
      }
      else if (j == 1 || j == 3)
      {
        r = h_radius_ch[i][j]->Fit("gaus", "SM", "", 3.5, 6.5);
        fitStatus = r;
        if (fitStatus >= 0)
        {
          if (j == 3)
          {
            properties["variables"]["RegenRejection"]["results"][method[i]]["charged"]["cylindrical"]["mean"][0] = r->Parameter(1);
            properties["variables"]["RegenRejection"]["results"][method[i]]["charged"]["cylindrical"]["width"][0] = r->Parameter(2);
            properties["variables"]["RegenRejection"]["errors"][method[i]]["charged"]["cylindrical"]["mean"][0] = r->Error(1);
            properties["variables"]["RegenRejection"]["errors"][method[i]]["charged"]["cylindrical"]["width"][0] = r->Error(2);
          }
          else
          {
            properties["variables"]["RegenRejection"]["results"][method[i]]["charged"]["cylindrical"]["mean"][1] = r->Parameter(1);
            properties["variables"]["RegenRejection"]["results"][method[i]]["charged"]["cylindrical"]["width"][1] = r->Parameter(2);
            properties["variables"]["RegenRejection"]["errors"][method[i]]["charged"]["cylindrical"]["mean"][1] = r->Error(1);
            properties["variables"]["RegenRejection"]["errors"][method[i]]["charged"]["cylindrical"]["width"][1] = r->Error(2);
          }

          cout << img_name << ": (" << r->Parameter(1) << "+-" << r->Error(1) << ")" << endl;
          cout << img_name << ": (" << r->Parameter(2) << "+-" << r->Error(2) << ")" << endl;
        }
        else
        {
          if (j == 3)
          {
            properties["variables"]["RegenRejection"]["results"][method[i]]["charged"]["cylindrical"]["mean"][0] = nullptr;
            properties["variables"]["RegenRejection"]["results"][method[i]]["charged"]["cylindrical"]["width"][0] = nullptr;
            properties["variables"]["RegenRejection"]["errors"][method[i]]["charged"]["cylindrical"]["mean"][0] = nullptr;
            properties["variables"]["RegenRejection"]["errors"][method[i]]["charged"]["cylindrical"]["width"][0] = nullptr;
          }
          else
          {
            properties["variables"]["RegenRejection"]["results"][method[i]]["charged"]["cylindrical"]["mean"][1] = nullptr;
            properties["variables"]["RegenRejection"]["results"][method[i]]["charged"]["cylindrical"]["width"][1] = nullptr;
            properties["variables"]["RegenRejection"]["errors"][method[i]]["charged"]["cylindrical"]["mean"][1] = nullptr;
            properties["variables"]["RegenRejection"]["errors"][method[i]]["charged"]["cylindrical"]["width"][1] = nullptr;
          }
        }

        h_radius_ch[i][j]->GetXaxis()->SetTitle("#rho [cm]");
        h_radius_ch[i][j]->GetYaxis()->SetTitle(Form("Counts/(%.2f cm)", width_rho_charged));
      }

      h_radius_ch[i][j]->Draw();
      canvas[canva_num]->Print(img_name);

      canva_num++;
    }
  }

  for (Int_t i = 0; i < 2; i++)
    for (Int_t j = 0; j < 4; j++)
    {
      img_name = img_dir + "RegenerationAnalysis/radius_tri_" + method_dict[i] + "_" + region_dict[j] + ext_img;
      canvas[canva_num]->cd();
      canvas[canva_num]->SetLogy();

      if (j == 0 || j == 2)
      {

        if (j == 0)
          r = h_radius_tri[i][j]->Fit("gaus", "SM", "", 8., 13.);
        else
          r = h_radius_tri[i][j]->Fit("gaus", "SM", "", 4., 6.);

        fitStatus = r;

        if (fitStatus >= 0)
        {
          if (j == 2)
          {
            properties["variables"]["RegenRejection"]["results"][method[i]]["triangle"]["spherical"]["mean"][0] = r->Parameter(1);
            properties["variables"]["RegenRejection"]["results"][method[i]]["triangle"]["spherical"]["width"][0] = r->Parameter(2);
            properties["variables"]["RegenRejection"]["errors"][method[i]]["triangle"]["spherical"]["mean"][0] = r->Error(1);
            properties["variables"]["RegenRejection"]["errors"][method[i]]["triangle"]["spherical"]["width"][0] = r->Error(2);
          }
          else
          {
            properties["variables"]["RegenRejection"]["results"][method[i]]["triangle"]["spherical"]["mean"][1] = r->Parameter(1);
            properties["variables"]["RegenRejection"]["results"][method[i]]["triangle"]["spherical"]["width"][1] = r->Parameter(2);
            properties["variables"]["RegenRejection"]["errors"][method[i]]["triangle"]["spherical"]["mean"][1] = r->Error(1);
            properties["variables"]["RegenRejection"]["errors"][method[i]]["triangle"]["spherical"]["width"][1] = r->Error(2);
          }

          cout << img_name << ": (" << r->Parameter(1) << "+-" << r->Error(1) << ")" << endl;
          cout << img_name << ": (" << r->Parameter(2) << "+-" << r->Error(2) << ")" << endl;
        }
        else
        {
          if (j == 2)
          {
            properties["variables"]["RegenRejection"]["results"][method[i]]["triangle"]["spherical"]["mean"][0] = nullptr;
            properties["variables"]["RegenRejection"]["results"][method[i]]["triangle"]["spherical"]["width"][0] = nullptr;
            properties["variables"]["RegenRejection"]["errors"][method[i]]["triangle"]["spherical"]["mean"][0] = nullptr;
            properties["variables"]["RegenRejection"]["errors"][method[i]]["triangle"]["spherical"]["width"][0] = nullptr;
          }
          else
          {
            properties["variables"]["RegenRejection"]["results"][method[i]]["triangle"]["spherical"]["mean"][1] = nullptr;
            properties["variables"]["RegenRejection"]["results"][method[i]]["triangle"]["spherical"]["width"][1] = nullptr;
            properties["variables"]["RegenRejection"]["errors"][method[i]]["triangle"]["spherical"]["mean"][1] = nullptr;
            properties["variables"]["RegenRejection"]["errors"][method[i]]["triangle"]["spherical"]["width"][1] = nullptr;
          }
        }
        h_radius_tri[i][j]->GetXaxis()->SetTitle("R [cm]");
        h_radius_tri[i][j]->GetYaxis()->SetTitle(Form("Counts/(%.2f cm)", width_R_triangle));
      }
      else if (j == 1 || j == 3)
      {
        if (j == 1)
          r = h_radius_tri[i][j]->Fit("gaus", "SM", "", 8., 13.);
        else
          r = h_radius_tri[i][j]->Fit("gaus", "SM", "", 4., 6.);

        fitStatus = r;

        if (fitStatus >= 0)
        {
          if (j == 3)
          {
            properties["variables"]["RegenRejection"]["results"][method[i]]["triangle"]["cylindrical"]["mean"][0] = r->Parameter(1);
            properties["variables"]["RegenRejection"]["results"][method[i]]["triangle"]["cylindrical"]["width"][0] = r->Parameter(2);
            properties["variables"]["RegenRejection"]["errors"][method[i]]["triangle"]["cylindrical"]["mean"][0] = r->Error(1);
            properties["variables"]["RegenRejection"]["errors"][method[i]]["triangle"]["cylindrical"]["width"][0] = r->Error(2);
          }
          else
          {
            properties["variables"]["RegenRejection"]["results"][method[i]]["triangle"]["cylindrical"]["mean"][1] = r->Parameter(1);
            properties["variables"]["RegenRejection"]["results"][method[i]]["triangle"]["cylindrical"]["width"][1] = r->Parameter(2);
            properties["variables"]["RegenRejection"]["errors"][method[i]]["triangle"]["cylindrical"]["mean"][1] = r->Error(1);
            properties["variables"]["RegenRejection"]["errors"][method[i]]["triangle"]["cylindrical"]["width"][1] = r->Error(2);
          }

          cout << img_name << ": (" << r->Parameter(1) << "+-" << r->Error(1) << ")" << endl;
          cout << img_name << ": (" << r->Parameter(2) << "+-" << r->Error(2) << ")" << endl;
        }
        else
        {
          if (j == 3)
          {
            properties["variables"]["RegenRejection"]["results"][method[i]]["triangle"]["cylindrical"]["mean"][0] = nullptr;
            properties["variables"]["RegenRejection"]["results"][method[i]]["triangle"]["cylindrical"]["width"][0] = nullptr;
            properties["variables"]["RegenRejection"]["errors"][method[i]]["triangle"]["cylindrical"]["mean"][0] = nullptr;
            properties["variables"]["RegenRejection"]["errors"][method[i]]["triangle"]["cylindrical"]["width"][0] = nullptr;
          }
          else
          {
            properties["variables"]["RegenRejection"]["results"][method[i]]["triangle"]["cylindrical"]["mean"][1] = nullptr;
            properties["variables"]["RegenRejection"]["results"][method[i]]["triangle"]["cylindrical"]["width"][1] = nullptr;
            properties["variables"]["RegenRejection"]["errors"][method[i]]["triangle"]["cylindrical"]["mean"][1] = nullptr;
            properties["variables"]["RegenRejection"]["errors"][method[i]]["triangle"]["cylindrical"]["width"][1] = nullptr;
          }
        }

        h_radius_tri[i][j]->GetXaxis()->SetTitle("#rho [cm]");
        h_radius_tri[i][j]->GetYaxis()->SetTitle(Form("Counts/(%.2f cm)", width_rho_triangle));
      }

      h_radius_tri[i][j]->Draw();
      canvas[canva_num]->Print(img_name);

      canva_num++;
    }

  std::ofstream result(propName);
  result << properties.dump(4);
  result.close();

  // count_tot = 982142857;

  // std::cout << "Regen branching ratio for negative: " << count_neg_reg / (Float_t)count_tot << std::endl;
  // std::cout << "Regen branching ratio for positive: " << count_pos_reg / (Float_t)count_tot << std::endl;
  // std::cout << "Regen total branching: " << (count_neg_reg / (Float_t)count_tot) + (count_pos_reg / (Float_t)count_tot) << std::endl;
  // std::cout << "Regen events neg: " << count_neg_reg << std::endl;
  // std::cout << "Regen events pos: " << count_pos_reg << std::endl;
  // std::cout << "Events total: " << count_tot << std::endl;

  // std::cout << "Comparison of regen branching ratio for positive: " << count_pos_reg / (Float_t)(count_pos_reg + count_neg_reg) << std::endl;
  // std::cout << "Comparison of regen branching ratio for negative: " << count_neg_reg / (Float_t)(count_pos_reg + count_neg_reg) << std::endl;

  // std::cout << "Count sig total branching ratio: " << count_sig / (Float_t)(count_tot) << std::endl;
  // std::cout << "Count sig negative branching ratio: " << count_sig_neg / (Float_t)(count_tot) << " pm " << sqrt(pow(sqrt(count_sig_neg) / (Float_t)(count_tot), 2) + pow(sqrt(count_tot) * count_sig_neg / pow((Float_t)(count_tot), 2), 2)) << std::endl;
  // std::cout << "Count sig positive branching ratio: " << count_sig_pos / (Float_t)(count_tot) << " pm " << sqrt(pow(sqrt(count_sig_pos) / (Float_t)(count_tot), 2) + pow(sqrt(count_tot) * count_sig_pos / pow((Float_t)(count_tot), 2), 2)) << std::endl;

  return 0;
}