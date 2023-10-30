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
#include <TFitResult.h>

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
      Kchmc[9], Knemc[9], ip[3], ipmc[3], phi_mom[4], Dtmc, Tcl[50];
  UChar_t mctruth;

  chain->SetBranchAddress("Kchboost", Kchboost);
  chain->SetBranchAddress("Knereclor", Knereclor);
  chain->SetBranchAddress("Knerec", Knerec);
  chain->SetBranchAddress("Kchmc", Kchmc);
  chain->SetBranchAddress("Knemc", Knemc);

  chain->SetBranchAddress("ip", ip);
  chain->SetBranchAddress("ipmc", ipmc);

  chain->SetBranchAddress("Tcl", Tcl);

  chain->SetBranchAddress("Bpx", &phi_mom[0]);
  chain->SetBranchAddress("Bpy", &phi_mom[1]);
  chain->SetBranchAddress("Bpz", &phi_mom[2]);
  chain->SetBranchAddress("Broots", &phi_mom[3]);

  chain->SetBranchAddress("mctruth", &mctruth);

  chain->SetBranchAddress("Dtmc", &Dtmc);

  Float_t gamma_kinfit[4][8], ip_kinfit[3], Knetri_kinfit[10], chi2min;
  Int_t done_kinfit, g4taken_kinfit[4];

  tree->SetBranchAddress("fourgamma1tri_kinfit", gamma_kinfit[0]);
  tree->SetBranchAddress("fourgamma2tri_kinfit", gamma_kinfit[1]);
  tree->SetBranchAddress("fourgamma3tri_kinfit", gamma_kinfit[2]);
  tree->SetBranchAddress("fourgamma4tri_kinfit", gamma_kinfit[3]);

  tree->SetBranchAddress("fourKnetri_kinfit", Knetri_kinfit);

  tree->SetBranchAddress("iptri_kinfit", ip_kinfit);
  tree->SetBranchAddress("done4_kinfit", &done_kinfit);

  tree->SetBranchAddress("g4taken_kinfit", g4taken_kinfit);

  tree->SetBranchAddress("chi2min", &chi2min);

  chain->AddFriend(tree);

  UInt_t nentries = (UInt_t)chain->GetEntries();

  Float_t lengthneu_mc, lengthneu_tri, lengthneu_rec, lengthch_mc, lengthch_rec,
      v_Kneutri, v_Kchrec, v_Kneurec, v_Kneumc, v_Kchmc,
      t_chmc, t_neumc, t_chrec, t_neutri, t_neurec, Rtneu_rec, Rtneu_mc, Rtneu_tri,
      angle_path_mom, cos_path_mom;

  TString id_hist;
  TH1 *dttri_hist, *prob_hist;
  TH2 *neu_vtx_corr[4], *neu_mom[4], *ip_coor[3];

  TH2 *sigmas_std[5], *sigmas_tri[5], *sigmas_tri_kin_fit[5];
  TH1 *neu_std_hist[5], *neu_tri_hist[5], *neu_tri_kin_fit_hist[5];
  TH1 *res_std_hist[5], *res_tri_hist[5], *res_tri_kin_fit_hist[5];

  for (Int_t i = 0; i < 3; i++)
  {
    id_hist = "Coor" + std::to_string(i);
    neu_vtx_corr[i] = new TH2F(id_hist, "", 100, -50, 50, 100, -50, 50);

    id_hist = "Mom" + std::to_string(i);
    neu_mom[i] = new TH2F(id_hist, "", 200, -300, 300, 200, -300, 300);

    id_hist = "IP coor" + std::to_string(i);
    ip_coor[i] = new TH2F(id_hist, "", 100, -100, 100, 100, -100, 100);
  }

  neu_vtx_corr[3] = new TH2F("Lengths", "", 100, 0, 20, 100, 0, 20);
  neu_mom[3] = new TH2F("Energies", "", 100, 504, 520, 100, 480, 600);

  dttri_hist = new TH1F("Dttri", "", 100, 0, 100);
  prob_hist = new TH1F("Prob", "", 100, 0, 1);

  TCanvas *canvas[30];

  const UInt_t number_of_points = 21;

  TString name[4] = {"x", "y", "z", "length"};

  for (Int_t i = 0; i < 4; i++)
  {
    if (i < 2)
    {
      neu_std_hist[i] = new TH1F(name[i] + " std coor", "", 100, -200, 200);
      neu_tri_hist[i] = new TH1F(name[i] + " tri coor", "", 100, -200, 200);
      neu_tri_kin_fit_hist[i] = new TH1F(name[i] + " kinfit coor", "", 100, -165, 165);
      res_std_hist[i] = new TH1F(name[i] + " std res", "", number_of_points, 0, 50);
      res_tri_hist[i] = new TH1F(name[i] + " tri res", "", number_of_points, 0, 50);
      res_tri_kin_fit_hist[i] = new TH1F(name[i] + " kinfit res", "", number_of_points, 0, 50);
      sigmas_std[i] = new TH2F(name[i] + " sigmas std", "", number_of_points, 0, 50, 100, -100, 100);
      sigmas_tri[i] = new TH2F(name[i] + " sigmas tri", "", number_of_points, 0, 50, 100, -100, 100);
      sigmas_tri_kin_fit[i] = new TH2F(name[i] + " sigmas kinfit", "", number_of_points, 0, 50, 100, -100, 100);
    }
    else if (i == 2)
    {
      neu_std_hist[i] = new TH1F(name[i] + " std", "", 100, 0, 5);
      neu_tri_hist[i] = new TH1F(name[i] + " tri", "", 100, 0, 5);
      neu_tri_kin_fit_hist[i] = new TH1F(name[i] + " kinfit", "", 100, 0, 5);
      res_std_hist[i] = new TH1F(name[i] + " std res", "", number_of_points, 0, 50);
      res_tri_hist[i] = new TH1F(name[i] + " tri res", "", number_of_points, 0, 50);
      res_tri_kin_fit_hist[i] = new TH1F(name[i] + " kinfit res", "", number_of_points, 0, 50);
      sigmas_std[i] = new TH2F(name[i] + " sigmas std", "", number_of_points, 0, 50, 100, -100, 100);
      sigmas_tri[i] = new TH2F(name[i] + " sigmas tri", "", number_of_points, 0, 50, 100, -100, 100);
      sigmas_tri_kin_fit[i] = new TH2F(name[i] + " sigmas kinfit", "", number_of_points, 0, 50, 100, -100, 100);
    }
    else if (i == 3)
    {
      neu_std_hist[i] = new TH1F(name[i] + " std", "", 100, 0, 5);
      neu_tri_hist[i] = new TH1F(name[i] + " tri", "", 100, 0, 5);
      neu_tri_kin_fit_hist[i] = new TH1F(name[i] + " kinfit", "", 100, 0, 5);
      res_std_hist[i] = new TH1F(name[i] + " std res", "", number_of_points, 0, 90);
      res_tri_hist[i] = new TH1F(name[i] + " tri res", "", number_of_points, 0, 90);
      res_tri_kin_fit_hist[i] = new TH1F(name[i] + " kinfit res", "", number_of_points, 0, 90);
      sigmas_std[i] = new TH2F(name[i] + " sigmas std", "", number_of_points, 0, 90, 200, -100, 100);
      sigmas_tri[i] = new TH2F(name[i] + " sigmas tri", "", number_of_points, 0, 90, 200, -100, 100);
      sigmas_tri_kin_fit[i] = new TH2F(name[i] + " sigmas kinfit", "", number_of_points, 0, 90, 200, -100, 100);
    }

  }

  Float_t gamma_path[4], cluster_time[4];

  for (Int_t i = 0; i < 30; i++)
    canvas[i] = new TCanvas(("Canvas" + std::to_string(i)).c_str(), "", 750, 750);

  for (Int_t i = 0; i < nentries; i++)
  {
    chain->GetEntry(i);

    if (done_kinfit == 1 && chi2min < 100)
    {
      ip_coor[0]->Fill(ipmc[0], ip_kinfit[0]);
      ip_coor[1]->Fill(ipmc[1], ip_kinfit[1]);
      ip_coor[2]->Fill(ipmc[2], ip_kinfit[2]);

      neu_vtx_corr[0]->Fill(Knemc[6], Knetri_kinfit[6]);
      neu_vtx_corr[1]->Fill(Knemc[7], Knetri_kinfit[7]);
      neu_vtx_corr[2]->Fill(Knemc[8], Knetri_kinfit[8]);

      neu_mom[0]->Fill(Knemc[0], Knetri_kinfit[0]);
      neu_mom[1]->Fill(Knemc[1], Knetri_kinfit[1]);
      neu_mom[2]->Fill(Knemc[2], Knetri_kinfit[2]);
      neu_mom[3]->Fill(Knemc[3], Knetri_kinfit[3]);

      lengthneu_mc = sqrt(pow(Knemc[6] - ipmc[0], 2) + pow(Knemc[7] - ipmc[1], 2) + pow(Knemc[8] - ipmc[2], 2));
      lengthneu_tri = sqrt(pow(Knetri_kinfit[6] - ip_kinfit[0], 2) + pow(Knetri_kinfit[7] - ip_kinfit[1], 2) + pow(Knetri_kinfit[8] - ip_kinfit[2], 2));
      lengthneu_rec = sqrt(pow(Knereclor[6] - ip[0], 2) + pow(Knereclor[7] - ip[1], 2) + pow(Knereclor[8] - ip[2], 2));

      lengthch_mc = sqrt(pow(Kchmc[6] - ipmc[0], 2) + pow(Kchmc[7] - ipmc[1], 2) + pow(Kchmc[8] - ipmc[2], 2));
      lengthch_rec = sqrt(pow(Kchboost[6] - ip[0], 2) + pow(Kchboost[7] - ip[1], 2) + pow(Kchboost[8] - ip[2], 2));

      Rtneu_mc = sqrt(pow(Knemc[6] - ipmc[0], 2) + pow(Knemc[7] - ipmc[1], 2));
      Rtneu_tri = sqrt(pow(Knetri_kinfit[6] - ip_kinfit[0], 2) + pow(Knetri_kinfit[7] - ip_kinfit[1], 2));
      Rtneu_rec = sqrt(pow(Knereclor[6] - ip[0], 2) + pow(Knereclor[7] - ip[1], 2));

      v_Kneumc = c_vel * Knemc[4] / Knemc[3];
      v_Kneutri = c_vel * Knetri_kinfit[4] / Knetri_kinfit[3];
      v_Kneurec = c_vel * Knereclor[5] / Knereclor[3];

      v_Kchmc = c_vel * Kchmc[4] / Kchmc[3];
      v_Kchrec = c_vel * Kchboost[4] / Kchboost[3];

      t_chrec = lengthch_rec / v_Kchrec;
      t_chmc = lengthch_mc / v_Kchmc;

      t_neurec = lengthneu_rec / v_Kneurec;
      t_neutri = lengthneu_tri / v_Kneutri;
      t_neumc = lengthneu_mc / v_Kneumc;

      cos_path_mom = ((Knetri_kinfit[6] - ip_kinfit[0])*Knetri_kinfit[0] + (Knetri_kinfit[7] - ip_kinfit[1])*Knetri_kinfit[1] + (Knetri_kinfit[8] - ip_kinfit[2])*Knetri_kinfit[2])/ (lengthneu_tri*sqrt(pow(Knetri_kinfit[0],2) + pow(Knetri_kinfit[1],2) + pow(Knetri_kinfit[2],2)));
      angle_path_mom = M_PI*acos(cos_path_mom)/180.;

      neu_vtx_corr[3]->Fill(t_neumc, Knetri_kinfit[9]);

      cluster_time[0] = Tcl[g4taken_kinfit[0]];
      cluster_time[1] = Tcl[g4taken_kinfit[1]];
      cluster_time[2] = Tcl[g4taken_kinfit[2]];
      cluster_time[3] = Tcl[g4taken_kinfit[3]];

      cluster_time[0] = Tcl[g4taken_kinfit[0]];
      cluster_time[1] = Tcl[g4taken_kinfit[1]];
      cluster_time[2] = Tcl[g4taken_kinfit[2]];
      cluster_time[3] = Tcl[g4taken_kinfit[3]];

      dttri_hist->Fill(Knetri_kinfit[9]/tau_S_nonCPT);
      prob_hist->Fill(TMath::Prob(chi2min, 6));

      sigmas_std[0]->Fill(abs(Knemc[6] - ipmc[0]), Knetri_kinfit[6] - Knemc[6]);
      sigmas_std[1]->Fill(abs(Knemc[7] - ipmc[1]), Knetri_kinfit[7] - Knemc[7]);
      sigmas_std[2]->Fill(abs(Knemc[8] - ipmc[2]), Knetri_kinfit[8] - Knemc[8]);
      sigmas_std[3]->Fill(t_neumc/tau_S_nonCPT, (Knetri_kinfit[9] - t_neumc)/tau_S_nonCPT);
    }
  }

  gStyle->SetOptStat("iMR");
  gStyle->SetStatX(0.85);
  gStyle->SetStatY(0.9);

  TCanvas *c[20];
  TString id_canva, x_title, y_title;
  Int_t width = 750, height = 750;

  for (Int_t i = 0; i < 3; i++)
  {
    id_canva = "coor" + std::to_string(i + 1);
    c[i] = new TCanvas(id_canva, "", width, height);
    c[i]->SetRightMargin(0.15);
    c[i]->cd();

    x_title = Form("#vec{X}_{neu,%d}^{gen} [cm]", i + 1);
    y_title = Form("#vec{X}_{neu,%d}^{tri} [cm]", i + 1);

    neu_vtx_corr[i]->GetXaxis()->SetTitle(x_title);
    neu_vtx_corr[i]->GetYaxis()->SetTitle(y_title);
    neu_vtx_corr[i]->Draw("COLZ");

    c[i]->Print(id_canva + ".png");

    id_canva = "ip_coor" + std::to_string(i + 1);
    c[i + 6] = new TCanvas(id_canva, "", width, height);
    c[i + 6]->SetRightMargin(0.15);
    c[i + 6]->cd();

    x_title = Form("#vec{X}_{IP,%d}^{gen} [cm]", i + 1);
    y_title = Form("#vec{X}_{IP,%d}^{tri} [cm]", i + 1);

    ip_coor[i]->GetXaxis()->SetTitle(x_title);
    ip_coor[i]->GetYaxis()->SetTitle(y_title);
    ip_coor[i]->Draw("COLZ");

    c[i + 6]->Print(id_canva + ".png");

    id_canva = "mom" + std::to_string(i + 1);
    c[i + 9] = new TCanvas(id_canva, "", width, height);
    c[i + 9]->SetRightMargin(0.15);
    c[i + 9]->cd();

    x_title = Form("#vec{p}_{K#rightarrow#pi^{0}#pi^{0},%d}^{gen} [MeV/c]", i + 1);
    y_title = Form("#vec{p}_{K#rightarrow#pi^{0}#pi^{0},%d}^{tri} [MeV/c]", i + 1);

    neu_mom[i]->GetXaxis()->SetTitle(x_title);
    neu_mom[i]->GetYaxis()->SetTitle(y_title);
    neu_mom[i]->Draw("COLZ");

    c[i + 9]->Print(id_canva + ".png");
  }

  id_canva = "time_neutral";
  c[3] = new TCanvas(id_canva, "", width, height);
  c[3]->SetRightMargin(0.15);
  c[3]->cd();

  x_title = Form("t_{neu}^{gen} [ns]");
  y_title = Form("t_{neu}^{tri} [ns]");

  neu_vtx_corr[3]->GetXaxis()->SetTitle(x_title);
  neu_vtx_corr[3]->GetYaxis()->SetTitle(y_title);
  neu_vtx_corr[3]->Draw("COLZ");

  c[3]->Print(id_canva + ".png");

  id_canva = "energy";
  c[12] = new TCanvas(id_canva, "", width, height);
  c[12]->SetRightMargin(0.15);
  c[12]->cd();

  x_title = Form("E_{K#rightarrow#pi^{0}#pi^{0}}^{gen} [MeV]");
  y_title = Form("E_{K#rightarrow#pi^{0}#pi^{0}}^{tri} [MeV]");

  neu_mom[3]->GetXaxis()->SetTitle(x_title);
  neu_mom[3]->GetYaxis()->SetTitle(y_title);
  neu_mom[3]->Draw("COLZ");

  c[12]->Print(id_canva + ".png");

  TF1 *func1 = new TF1("func", chi2dist, 0, 50, 2);
  func1->SetParameters(dttri_hist->Integral(0, 100) / 10., 6);
  func1->SetParNames("Normalization", "Degrees of freedom");

  id_canva = "Canva" + std::to_string(4);
  c[4] = new TCanvas(id_canva, "", width, height);
  c[4]->SetRightMargin(0.15);
  c[4]->cd();

  dttri_hist->GetYaxis()->SetMaxDigits(3);

  dttri_hist->GetXaxis()->CenterTitle();
  dttri_hist->GetYaxis()->CenterTitle();
  dttri_hist->GetXaxis()->SetTitle("#chi^{2}_{tri}");
  dttri_hist->GetYaxis()->SetTitle("Counts");
  dttri_hist->Draw();
  c[4]->Print(id_canva + ".png");

  id_canva = "Canva" + std::to_string(5);
  c[5] = new TCanvas(id_canva, "", width, height);
  c[5]->SetRightMargin(0.15);
  c[5]->cd();

  prob_hist->GetYaxis()->SetMaxDigits(3);

  prob_hist->GetXaxis()->CenterTitle();
  prob_hist->GetYaxis()->CenterTitle();
  prob_hist->GetXaxis()->SetTitle("Prob(#chi^{2},6)");
  prob_hist->GetYaxis()->SetTitle("Counts");
  prob_hist->Draw();
  c[5]->Print(id_canva + ".png");

  TFitResultPtr r;
  Float_t parameter[3];

  TF1 *func = new TF1("f1", "[0]*exp(-0.5*pow((x-[1])/[2],2))", -200, 200);
  func->SetParNames("Constant", "Mean", "Sigma");

  Color_t res_color[4] = {kRed, kBlue, kGreen, kBlack};

  for (Int_t j = 0; j < 4; j++)
  {
    for (Int_t i = 1; i <= number_of_points; i++)
    {
      parameter[0] = sigmas_std[j]->ProjectionY("_py", i, i)->Integral();
      parameter[1] = sigmas_std[j]->ProjectionY("_py", i, i)->GetBinCenter(sigmas_std[j]->ProjectionY("_py", i, i)->GetMaximumBin());
      parameter[2] = sigmas_std[j]->ProjectionY("_py", i, i)->GetStdDev();

      func->SetParameters(parameter[0], parameter[1], parameter[2]);

      func->SetParLimits(2,0.0001,100.);

      r = sigmas_std[j]->ProjectionY("_py", i, i)->Fit("f1", "S");
      res_std_hist[j]->SetBinContent(i, r->Parameter(2));
      res_std_hist[j]->SetBinError(i, r->ParError(2));
    }

    res_std_hist[j]->SetLineColor(res_color[j]);
  }

  TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9, "Axis");

  legend->AddEntry(res_std_hist[0],"x", "le");
  legend->AddEntry(res_std_hist[1],"y", "le");
  legend->AddEntry(res_std_hist[2],"z", "le");

  TString title_x[2] = {"|#vec{X}^{gen}_{K#rightarrow#pi^{0}#pi^{0},i} - #vec{X}^{gen}_{#Phi,i}| [cm]",
                        "t_{K#rightarrow#pi^{0}#pi^{0}}^{gen} [#tau_{S}]"};

  TString title_y[2] = {"#sigma(#vec{X}^{rec}_{K#rightarrow#pi^{0}#pi^{0},i} - #vec{X}^{gen}_{K#rightarrow#pi^{0}#pi^{0},i}) [cm]",
                        "#sigma(t_{K#rightarrow#pi^{0}#pi^{0}}^{rec} - t_{K#rightarrow#pi^{0}#pi^{0}}^{gen}) [#tau_{S}]"};

  gStyle->SetOptStat(0);

  canvas[0]->cd();
  res_std_hist[0]->SetXTitle(title_x[0]);
  res_std_hist[0]->SetYTitle(title_y[0]);
  res_std_hist[0]->GetYaxis()->SetRangeUser(0.0, 15.0);
  res_std_hist[0]->Draw("PE1");
  res_std_hist[1]->Draw("PE1SAME");
  res_std_hist[2]->Draw("PE1SAME");
  legend->Draw();
  canvas[0]->Print(("sigmas_std" + std::to_string(0) + ".png").c_str());

  canvas[1]->cd();
  res_std_hist[3]->SetXTitle(title_x[1]);
  res_std_hist[3]->SetYTitle(title_y[1]);
  res_std_hist[3]->GetYaxis()->SetRangeUser(0.0, 15.0);
  res_std_hist[3]->Draw("PE1");
  canvas[1]->Print(("sigmas_std" + std::to_string(2) + ".png").c_str());

  return 0;
}