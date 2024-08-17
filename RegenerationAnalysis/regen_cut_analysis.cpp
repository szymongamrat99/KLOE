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

#include "chain_init.C"

using namespace std;

void regen_plots()
{
  TChain *chain = new TChain("INTERF/h1");
  chain_init(chain);

  Float_t Dtboostlor = 0., Dtmc = 0., Kchboost[9] = {0.}, Knereclor[9] = {0.};
  UChar_t mcflag = 0;

  chain->SetBranchAddress("Dtboostlor", &Dtboostlor);
  chain->SetBranchAddress("Dtmc", &Dtmc);
  chain->SetBranchAddress("Kchboost", Kchboost);
  chain->SetBranchAddress("Knereclor", Knereclor);
  chain->SetBranchAddress("mcflag", &mcflag);


  TFile file_tri("../Neutrec/root_files/neuvtx_tri_kin_fit_1_56_10_9_1_2_bunch_corr_automatic.root");
  TTree *tree_tri = (TTree *)file_tri.Get("h_tri_kin_fit");

  Int_t done4 = 0, g4takentri_kinfit[4];
  Float_t fourKnetri[10] = {0.};

  tree_tri->SetBranchAddress("fourKnetri_kinfit", fourKnetri);
  tree_tri->SetBranchAddress("done4_kinfit", &done4);

  TFile file("../../../old_root_files/mctruth.root");
  TTree *tree = (TTree *)file.Get("h1");

  Int_t mctruth = 0;

  tree->SetBranchAddress("mctruth", &mctruth);

  Int_t nentries = chain->GetEntries();

  chain->AddFriend(tree_tri);
  chain->AddFriend(tree);

  Double_t height = 750, width = 750;
  TString dir = "img/", ext = ".svg";

  TCanvas *canvas[100];
  TH1 *h_radius[2][4], *h_radius_tri[2][4], *h_radius_ch[2][4], *h_dt[7];

  TString hist_name, canvas_name;

  for (Int_t i = 0; i < 100; i++)
  {
    canvas_name = "Canva_" + to_string(i);
    canvas[i] = new TCanvas(canvas_name, canvas_name, width, height);
    canvas[i]->SetLeftMargin(0.2);
    canvas[i]->SetBottomMargin(0.2);
  }

  Int_t nbins = 100;
  Double_t minx = 0.0, maxx = 50.0, one_bin_val = (maxx - minx) / (Double_t)nbins;

  map<Int_t, TString> region_dict, method_dict;

  region_dict[0] = "R_blue";
  region_dict[1] = "#rho_blue";
  region_dict[2] = "R_green";
  region_dict[3] = "#rho_green";

  method_dict[0] = "z-coor_section";
  method_dict[1] = "#rho_section";

  Int_t bins_rho_triangle = Int_t(50. / 1.29), bins_path_triangle = Int_t(50. / 1.14), bins_rho_tri = Int_t(50. / 2.58), bins_path_tri = Int_t(50. / 2.11);  

  for (Int_t i = 0; i < 2; i++)
    for (Int_t j = 0; j < 4; j++)
    {
      hist_name = "Hist_radius_triangle_" + method_dict[i] + "_" + region_dict[j];
      if(j == 0 || j == 2)
        h_radius[i][j] = new TH1D(hist_name, hist_name, bins_path_triangle, 0, 50.0);
      else
        h_radius[i][j] = new TH1D(hist_name, hist_name, bins_rho_triangle, 0, 50.0);

      hist_name = "Hist_radius_tri_" + method_dict[i] + "_" + region_dict[j];
      if(j == 0 || j == 2)
        h_radius_tri[i][j] = new TH1D(hist_name, hist_name, bins_path_tri, 0, 50.0);
      else
        h_radius_tri[i][j] = new TH1D(hist_name, hist_name, bins_rho_tri, 0, 50.0);

      hist_name = "Hist_radius_ch_" + method_dict[i] + "_" + region_dict[j];
      h_radius_ch[i][j] = new TH1D(hist_name, hist_name, 100, 0, 50.0);

      h_radius[i][j]->GetYaxis()->SetMaxDigits(3);
      h_radius_tri[i][j]->GetYaxis()->SetMaxDigits(3);
      h_radius_ch[i][j]->GetYaxis()->SetMaxDigits(3);
    }

  h_dt[0] = new TH1D("Delta T bef", "Delta T before / after cut;#Deltat [#tau_{S}];Counts / 2.1#tau_{S}", 91, -90.0, 90.0);

  h_dt[1] = new TH1D("Delta T aft1", "Delta T after cut 1;#Deltat [#tau_{S}];Counts / 2.1#tau_{S}", 91, -90.0, 90.0);

  h_dt[2] = new TH1D("Delta T aft2", "Delta T after cut 2;#Deltat [#tau_{S}];Counts / 2.1#tau_{S}", 91, -90.0, 90.0);

  h_dt[3] = new TH1D("Delta T aft3", "Delta T after cut 3;#Deltat [#tau_{S}];Counts / 2.1#tau_{S}", 91, -90.0, 90.0);

  h_dt[4] = new TH1D("Delta T aft4", "Delta T after cut 4;#Deltat [#tau_{S}];Counts / 2.1#tau_{S}", 91, -90.0, 90.0);

  h_dt[5] = new TH1D("Delta T aft5", "Delta T after cut 5;#Deltat [#tau_{S}];Counts / 2.1#tau_{S}", 91, -90.0, 90.0);

  h_dt[6] = new TH1D("Delta T aft6", "Delta T after cut 6;#Deltat [#tau_{S}];Counts / 2.1#tau_{S}", 91, -90.0, 90.0);

  Double_t sphere_bound = 10, bp_bound = 4.4;

  Double_t sigma = 1.5;

  Long64_t count_tot = 0, count_neg_reg = 0, count_pos_reg = 0, count_sig = 0, count_sig_neg = 0, count_sig_pos = 0;

  for (Int_t i = 0; i < nentries; i++)
  {
    chain->GetEntry(i);

    Bool_t cuts[6] = {false, false, false, false, false, false};
    Double_t radius[2] = {0., 0.}, radius_ch[2] = {0., 0.}, radius_tri[2] = {0., 0.};

    for (Int_t i = 0; i < 3; i++)
    {
      radius[0] += pow(Knereclor[6 + i], 2);
      radius_ch[0] += pow(Kchboost[6 + i], 2);
      radius_tri[0] += pow(fourKnetri[6 + i], 2);

      if (i < 2)
      {
        radius[1] += pow(Knereclor[6 + i], 2);
        radius_ch[1] += pow(Kchboost[6 + i], 2);
        radius_tri[1] += pow(fourKnetri[6 + i], 2);
      }
    }

    radius[0] = sqrt(radius[0]);
    radius[1] = sqrt(radius[1]);

    radius_ch[0] = sqrt(radius_ch[0]);
    radius_ch[1] = sqrt(radius_ch[1]);

    radius_tri[0] = sqrt(radius_tri[0]);
    radius_tri[1] = sqrt(radius_tri[1]);

    if(mcflag == 1 && (mctruth == 0 || mctruth == 2 || mctruth == 1))
    {
      count_sig++;

      if(Dtmc < 0) count_sig_neg++;
      else count_sig_pos++;
    } 
    
    count_tot++;


    if (mctruth == 3)
    {
      // Neutral vtx triangle method

      if (abs(Knereclor[8]) > 10)
      {
        h_radius[0][0]->Fill(radius[0]);
        h_radius[0][1]->Fill(radius[1]);
      }
      else
      {
        h_radius[0][2]->Fill(radius[0]);
        h_radius[0][3]->Fill(radius[1]);
      }

      if (radius[1] > 8 && radius[1] < 20)
      {
        h_radius[1][0]->Fill(radius[0]);
        h_radius[1][1]->Fill(radius[1]);
      }
      else if (radius[1] < 8)
      {
        h_radius[1][2]->Fill(radius[0]);
        h_radius[1][3]->Fill(radius[1]);
      }

      // Charged vtx

      if (abs(Kchboost[8]) > 10)
      {
        h_radius_ch[0][0]->Fill(radius_ch[0]);
        h_radius_ch[0][1]->Fill(radius_ch[1]);
      }
      else
      {
        h_radius_ch[0][2]->Fill(radius_ch[0]);
        h_radius_ch[0][3]->Fill(radius_ch[1]);
      }

      if (radius_ch[1] > 8 && radius_ch[1] < 20)
      {
        h_radius_ch[1][0]->Fill(radius_ch[0]);
        h_radius_ch[1][1]->Fill(radius_ch[1]);
      }
      else if (radius_ch[1] < 8)
      {
        h_radius_ch[1][2]->Fill(radius_ch[0]);
        h_radius_ch[1][3]->Fill(radius_ch[1]);
      }

      // Neutral vtx trilateration method

      if (done4 == 1)
      {
        if (abs(fourKnetri[8]) > 10)
        {
          h_radius_tri[0][0]->Fill(radius_tri[0]);
          h_radius_tri[0][1]->Fill(radius_tri[1]);
        }
        else
        {
          h_radius_tri[0][2]->Fill(radius_tri[0]);
          h_radius_tri[0][3]->Fill(radius_tri[1]);
        }

        if (radius_tri[1] > 8 && radius_tri[1] < 20)
        {
          h_radius_tri[1][0]->Fill(radius_tri[0]);
          h_radius_tri[1][1]->Fill(radius_tri[1]);
        }
        else if (radius_tri[1] < 8)
        {
          h_radius_tri[1][2]->Fill(radius_tri[0]);
          h_radius_tri[1][3]->Fill(radius_tri[1]);
        }
      }

      cuts[0] = abs(radius[0] - 11.2386) > 1.0456 * sigma;// && radius[1] > 8;
      cuts[1] = abs(radius[1] - 5.05747) > 1.54588 * sigma;// && radius[1] <= 8;
      cuts[2] = abs(radius_ch[0] - 10.5840) > 0.681656 * sigma;// && radius_ch[1] > 8;
      cuts[3] = abs(radius_ch[1] - 4.45238) > 1.29837 * sigma;// && radius_ch[1] <= 8;
      cuts[4] = abs(radius_tri[0] - 11.0739) > 1.36243 * sigma;// && radius_tri[1] > 8;
      cuts[5] = abs(radius_tri[1] - 4.5) > 1.0 * sigma;// && radius_tri[1] <= 8;

      h_dt[0]->Fill(Dtboostlor);
      if (cuts[0])
        h_dt[1]->Fill(Dtboostlor);

      if (cuts[1])
        h_dt[2]->Fill(Dtboostlor);

      if (cuts[2])
        h_dt[3]->Fill(Dtboostlor);

      if (cuts[3])
        h_dt[4]->Fill(Dtboostlor);

      if (cuts[4])
        h_dt[5]->Fill(Dtboostlor);
        
      if (cuts[5])
        h_dt[6]->Fill(Dtboostlor);
    }
  }


  TFitResultPtr r;
  Int_t fitStatus = 0;

  TString img_name;

  Int_t canva_num = 0;

  ofstream myfile;
  myfile.open("fit_results.txt");

  for (Int_t i = 0; i < 2; i++)
    for (Int_t j = 0; j < 4; j++)
    {
      img_name = dir + "radius_" + method_dict[i] + "_" + region_dict[j] + ext;
      canvas[canva_num]->cd();
      canvas[canva_num]->SetLogy();

      if (j == 0 || j == 2)
      {
        r = h_radius[i][j]->Fit("gaus", "SM", "", 8., 13.);
        fitStatus = r;
        if (fitStatus >= 0)
        {
          myfile << img_name << ": (" << r->Parameter(1) << "+-" << r->Error(1) << ")" << endl;
          myfile << img_name <<": (" << r->Parameter(2) << "+-" << r->Error(2) << ")" << endl;

          cout << img_name << ": (" << r->Parameter(1) << "+-" << r->Error(1) << ")" << endl;
          cout << img_name <<": (" << r->Parameter(2) << "+-" << r->Error(2) << ")" << endl;
        }

        h_radius[i][j]->GetXaxis()->SetTitle("R [cm]");
        h_radius[i][j]->GetYaxis()->SetTitle(Form("Counts/(%.2f cm)", one_bin_val));
      }
      else if (j == 1 || j == 3)
      {
        r = h_radius[i][j]->Fit("gaus", "SM", "", 4., 7.);
        fitStatus = r;
        if (fitStatus >= 0)
        {
          myfile << img_name << ": (" << r->Parameter(1) << "+-" << r->Error(1) << ")" << endl;
          myfile << img_name <<": (" << r->Parameter(2) << "+-" << r->Error(2) << ")" << endl;
        }

        h_radius[i][j]->GetXaxis()->SetTitle("#rho [cm]");
        h_radius[i][j]->GetYaxis()->SetTitle(Form("Counts/(%.2f cm)", one_bin_val));
      }

      h_radius[i][j]->Draw();
      canvas[canva_num]->Print(img_name);

      canva_num++;
    }

  for (Int_t i = 0; i < 2; i++)
    for (Int_t j = 0; j < 4; j++)
    {
      img_name = dir + "radius_ch_" + method_dict[i] + "_" + region_dict[j] + ext;
      canvas[canva_num]->cd();
      canvas[canva_num]->SetLogy();

      if (j == 0 || j == 2)
      {
        r = h_radius_ch[i][j]->Fit("gaus", "SM", "", 9., 12.);
        fitStatus = r;
        if (fitStatus >= 0)
        {
          myfile << img_name << ": (" << r->Parameter(1) << "+-" << r->Error(1) << ")" << endl;
          myfile << img_name <<": (" << r->Parameter(2) << "+-" << r->Error(2) << ")" << endl;
        }

        h_radius_ch[i][j]->GetXaxis()->SetTitle("R [cm]");
        h_radius_ch[i][j]->GetYaxis()->SetTitle(Form("Counts/(%.2f cm)", one_bin_val));
      }
      else if (j == 1 || j == 3)
      {
        r = h_radius_ch[i][j]->Fit("gaus", "SM", "", 3, 6.);
        fitStatus = r;
        if (fitStatus >= 0)
        {
          myfile << img_name << ": (" << r->Parameter(1) << "+-" << r->Error(1) << ")" << endl;
          myfile << img_name <<": (" << r->Parameter(2) << "+-" << r->Error(2) << ")" << endl;
        }

        h_radius_ch[i][j]->GetXaxis()->SetTitle("#rho [cm]");
        h_radius_ch[i][j]->GetYaxis()->SetTitle(Form("Counts/(%.2f cm)", one_bin_val));
      }

      h_radius_ch[i][j]->Draw();
      canvas[canva_num]->Print(img_name);

      canva_num++;
    }

  for (Int_t i = 0; i < 2; i++)
    for (Int_t j = 0; j < 4; j++)
    {
      img_name = dir + "radius_tri_" + method_dict[i] + "_" + region_dict[j] + ext;
      canvas[canva_num]->cd();
      canvas[canva_num]->SetLogy();

      if (j == 0 || j == 2)
      {

        if (j == 0)
          r = h_radius_tri[i][j]->Fit("gaus", "SM", "", 8., 13.);
        else
          r = h_radius_tri[i][j]->Fit("gaus", "SM", "", 4., 5.);

        fitStatus = r;

        if (fitStatus >= 0)
        {
          myfile << img_name << ": (" << r->Parameter(1) << "+-" << r->Error(1) << ")" << endl;
          myfile << img_name <<": (" << r->Parameter(2) << "+-" << r->Error(2) << ")" << endl;
        }
        h_radius_tri[i][j]->GetXaxis()->SetTitle("R [cm]");
        h_radius_tri[i][j]->GetYaxis()->SetTitle(Form("Counts/(%.2f cm)", one_bin_val));
      }
      else if (j == 1 || j == 3)
      {
        if (j == 1)
          r = h_radius_tri[i][j]->Fit("gaus", "SM", "", 8., 13.);
        else
          r = h_radius_tri[i][j]->Fit("gaus", "SM", "", 4., 5.);

        fitStatus = r;

        if (fitStatus >= 0)
        {
          myfile << img_name << ": (" << r->Parameter(1) << "+-" << r->Error(1) << ")" << endl;
          myfile << img_name <<": (" << r->Parameter(2) << "+-" << r->Error(2) << ")" << endl;
        }

        h_radius_tri[i][j]->GetXaxis()->SetTitle("#rho [cm]");
        h_radius_tri[i][j]->GetYaxis()->SetTitle(Form("Counts/(%.2f cm)", one_bin_val));
      }

      h_radius_tri[i][j]->Draw();
      canvas[canva_num]->Print(img_name);

      canva_num++;
    }

  myfile.close();

  for(Int_t i = 1; i < 7; i++)
  {
    img_name = dir + "deltaT" + i + ext;
    canvas[canva_num]->cd();

    h_dt[0]->SetLineColor(kGreen);
    h_dt[0]->GetYaxis()->SetMaxDigits(3);

    h_dt[0]->Draw();
    h_dt[i]->Draw("SAME");

    canvas[canva_num]->Print(img_name);

    canva_num++;
  }

  count_tot = 982142857;

  std::cout << "Regen branching ratio for negative: "<< count_neg_reg / (Float_t)count_tot << std::endl;
  std::cout << "Regen branching ratio for positive: "<< count_pos_reg / (Float_t)count_tot << std::endl; 
  std::cout << "Regen total branching: "<< (count_neg_reg / (Float_t)count_tot) + (count_pos_reg / (Float_t)count_tot) << std::endl;
  std::cout << "Regen events neg: "<< count_neg_reg << std::endl; 
  std::cout << "Regen events pos: "<< count_pos_reg << std::endl; 
  std::cout << "Events total: "<< count_tot << std::endl;

  std::cout << "Comparison of regen branching ratio for positive: "<< count_pos_reg / (Float_t)(count_pos_reg + count_neg_reg) << std::endl;
  std::cout << "Comparison of regen branching ratio for negative: "<< count_neg_reg / (Float_t)(count_pos_reg + count_neg_reg) << std::endl;

  std::cout << "Count sig total branching ratio: "<< count_sig / (Float_t)(count_tot) << std::endl;
  std::cout << "Count sig negative branching ratio: "<< count_sig_neg / (Float_t)(count_tot) << " pm " << sqrt(pow(sqrt(count_sig_neg)/(Float_t)(count_tot),2) + pow(sqrt(count_tot)*count_sig_neg/pow((Float_t)(count_tot),2),2)) << std::endl; 
  std::cout << "Count sig positive branching ratio: "<< count_sig_pos / (Float_t)(count_tot) << " pm " << sqrt(pow(sqrt(count_sig_pos)/(Float_t)(count_tot),2) + pow(sqrt(count_tot)*count_sig_pos/pow((Float_t)(count_tot),2),2)) << std::endl;    
}