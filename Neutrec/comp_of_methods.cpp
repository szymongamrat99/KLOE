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
#include <TMatrixD.h>
#include <TVectorD.h>
#include <TPaveText.h>

#include "chain_init.C"
#include "../../Include/const.h"
#include "../../Include/Codes/kloe_class.h"
#include "../../Include/Codes/interf_function.h"
#include "../../Include/Codes/chi2_dist.h"
#include "../../Include/Codes/uncertainties.h"
#include "../../Include/Codes/lorentz_transf.h"
#include "../../Include/Codes/triple_gaus.h"

const Int_t M = 5;

int comp_of_methods()
{
  TChain *chain = new TChain("INTERF/h1");
  chain_init(chain, 1, 56);

  TFile *file = new TFile("neuvtx_tri_kin_fit_1_56_100_5.root");
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

  Float_t gamma_kinfit[4][8], ip_kinfit[3], Knetri_kinfit[10], chi2min, bhabha_mom_err[4];
  Int_t done_kinfit, g4taken_kinfit[4];
  TMatrixD *cov_matrix = new TMatrixD(27,27);
  TVectorD *min_const = new TVectorD(M);
  TVectorD *lag_mult = new TVectorD(M), *X_min = new TVectorD(27), *X_init = new TVectorD(27);

  tree->SetBranchAddress("fourgamma1tri_kinfit", gamma_kinfit[0]);
  tree->SetBranchAddress("fourgamma2tri_kinfit", gamma_kinfit[1]);
  tree->SetBranchAddress("fourgamma3tri_kinfit", gamma_kinfit[2]);
  tree->SetBranchAddress("fourgamma4tri_kinfit", gamma_kinfit[3]);

  chain->SetBranchAddress("Bwidpx", &bhabha_mom_err[0]);
	chain->SetBranchAddress("Bwidpy", &bhabha_mom_err[1]);
	chain->SetBranchAddress("Bwidpz", &bhabha_mom_err[2]);
	chain->SetBranchAddress("Brootserr", &bhabha_mom_err[3]);

  tree->SetBranchAddress("fourKnetri_kinfit", Knetri_kinfit);

  tree->SetBranchAddress("iptri_kinfit", ip_kinfit);
  tree->SetBranchAddress("done4_kinfit", &done_kinfit);

  //tree->SetBranchAddress("g4taken_kinfit", g4taken_kinfit);

  tree->SetBranchAddress("chi2min", &chi2min);
  tree->SetBranchAddress("min_cov", &cov_matrix);
  tree->SetBranchAddress("const_min", &min_const);
  tree->SetBranchAddress("lag_mult", &lag_mult);

  tree->SetBranchAddress("min_vars", &X_min);
  tree->SetBranchAddress("init_vars", &X_init);

  chain->AddFriend(tree);

  UInt_t nentries = (UInt_t)chain->GetEntries();

  Float_t lengthneu_mc, lengthneu_tri, lengthneu_rec, lengthch_mc, lengthch_rec,
      v_Kneutri, v_Kchrec, v_Kneurec, v_Kneumc, v_Kchmc,
      t_chmc, t_neumc, t_chrec, t_neutri, t_neurec, Rtneu_rec, Rtneu_mc, Rtneu_tri,
      angle_path_mom, cos_path_mom, sigmas_init[24];

  TString id_hist;
  TH1 *dttri_hist, *prob_hist;
  TH2 *neu_vtx_corr[4], *neu_mom[4], *ip_coor[3];

  TH2 *sigmas_std[5], *sigmas_tri[5], *sigmas_tri_kin_fit[5];
  TH1 *neu_std_hist[5], *neu_tri_hist[5], *neu_tri_kin_fit_hist[5];
  TH1 *res_std_hist[5], *res_tri_hist[5], *res_tri_kin_fit_hist[5];
  TH1 *pulls[8];

  for (Int_t i = 0; i < 3; i++)
  {
    id_hist = "Coor" + std::to_string(i);
    neu_vtx_corr[i] = new TH2F(id_hist, "", 100, -50, 50, 100, -50, 50);

    id_hist = "Mom" + std::to_string(i);
    neu_mom[i] = new TH2F(id_hist, "", 200, -300, 300, 200, -300, 300);

    id_hist = "IP coor" + std::to_string(i);
    ip_coor[i] = new TH2F(id_hist, "", 100, -100, 100, 100, -100, 100);
  }

  for(Int_t i = 0; i < 24; i++)
  {
    id_hist = "Var" + std::to_string(i);
    if(i == 0)
      pulls[i] = new TH1F(id_hist, ";#vec{p}^{tri}_{neu,1} - #vec{p}^{gen}_{neu,1} [MeV/c];Counts", 50, -300, 300);
    else if(i == 1)
      pulls[i] = new TH1F(id_hist, ";#vec{p}^{tri}_{neu,2} - #vec{p}^{gen}_{neu,2} [MeV/c];Counts", 50, -300, 300);
    else if(i == 2)
      pulls[i] = new TH1F(id_hist, ";#vec{p}^{tri}_{neu,3} - #vec{p}^{gen}_{neu,3} [MeV/c];Counts", 50, -300, 300);
    else if(i == 3)
      pulls[i] = new TH1F(id_hist, ";E^{tri}_{neu} - E^{gen}_{neu,3} [MeV];Counts", 50, -30, 30);
    else if(i == 4)
      pulls[i] = new TH1F(id_hist, ";#vec{X}^{tri}_{neu,1} - #vec{X}^{gen}_{neu,1} [cm];Counts", 50, -30, 30);
    else if(i == 5)
      pulls[i] = new TH1F(id_hist, ";#vec{X}^{tri}_{neu,2} - #vec{X}^{gen}_{neu,2} [cm];Counts", 50, -30, 30);
    else if(i == 6)
      pulls[i] = new TH1F(id_hist, ";#vec{X}^{tri}_{neu,3} - #vec{X}^{gen}_{neu,3} [cm];Counts", 50, -30, 30);
    else if(i == 7)
      pulls[i] = new TH1F(id_hist, ";t^{tri}_{neu} - t^{gen}_{neu} [#tau_{S}];Counts", 50, -20, 20);
  }

  neu_vtx_corr[3] = new TH2F("Lengths", "", 100, 0.0, 1.0, 100, 0, 20.0);
  neu_mom[3] = new TH2F("Energies", "", 100, 504, 520, 100, 507, 512);

  dttri_hist = new TH1F("Dttri", "", 200, 0.0, 50.0);
  prob_hist = new TH1F("Prob", "", 200, 0.0, 1.0);

  TH2 *corr_matrix_plot = new TH2F("corr_matrix_plot", "", 24, 0, 24, 24, 0, 24);

  TCanvas *canvas[30];

  const UInt_t number_of_points = 2;

  TString name[4] = {"x", "y", "z", "length"};

  for (Int_t i = 0; i < 4; i++)
  {
    if (i < 2)
    {
      neu_std_hist[i] = new TH1F(name[i] + " std coor", "", 100, -200, 200);
      neu_tri_hist[i] = new TH1F(name[i] + " tri coor", "", 100, -200, 200);
      neu_tri_kin_fit_hist[i] = new TH1F(name[i] + " kinfit coor", "", 100, -165, 165);
      res_std_hist[i] = new TH1F(name[i] + " std res", "", number_of_points, 0, 10);
      res_tri_hist[i] = new TH1F(name[i] + " tri res", "", number_of_points, 0, 10);
      res_tri_kin_fit_hist[i] = new TH1F(name[i] + " kinfit res", "", number_of_points, 0, 10);
      sigmas_std[i] = new TH2F(name[i] + " sigmas std", "", number_of_points, 0, 10, 25, -10, 10);
      sigmas_tri[i] = new TH2F(name[i] + " sigmas tri", "", number_of_points, 0, 10, 25, -10, 10);
      sigmas_tri_kin_fit[i] = new TH2F(name[i] + " sigmas kinfit", "", number_of_points, 0, 10, 25, -10, 10);
    }
    else if (i == 2)
    {
      neu_std_hist[i] = new TH1F(name[i] + " std", "", 100, 0, 5);
      neu_tri_hist[i] = new TH1F(name[i] + " tri", "", 100, 0, 5);
      neu_tri_kin_fit_hist[i] = new TH1F(name[i] + " kinfit", "", 100, 0, 5);
      res_std_hist[i] = new TH1F(name[i] + " std res", "", number_of_points, 0, 10);
      res_tri_hist[i] = new TH1F(name[i] + " tri res", "", number_of_points, 0, 10);
      res_tri_kin_fit_hist[i] = new TH1F(name[i] + " kinfit res", "", number_of_points, 0, 10);
      sigmas_std[i] = new TH2F(name[i] + " sigmas std", "", number_of_points, 0, 10, 25, -10, 10);
      sigmas_tri[i] = new TH2F(name[i] + " sigmas tri", "", number_of_points, 0, 10, 25, -10, 10);
      sigmas_tri_kin_fit[i] = new TH2F(name[i] + " sigmas kinfit", "", number_of_points, 0, 10, 25, -10, 10);
    }
    else if (i == 3)
    {
      neu_std_hist[i] = new TH1F(name[i] + " std", "", 100, 0, 5);
      neu_tri_hist[i] = new TH1F(name[i] + " tri", "", 100, 0, 5);
      neu_tri_kin_fit_hist[i] = new TH1F(name[i] + " kinfit", "", 100, 0, 5);
      res_std_hist[i] = new TH1F(name[i] + " std res", "", number_of_points, 0, 10);
      res_tri_hist[i] = new TH1F(name[i] + " tri res", "", number_of_points, 0, 10);
      res_tri_kin_fit_hist[i] = new TH1F(name[i] + " kinfit res", "", number_of_points, 0, 10);
      sigmas_std[i] = new TH2F(name[i] + " sigmas std", "", number_of_points, 0, 10, 50, -25, 25);
      sigmas_tri[i] = new TH2F(name[i] + " sigmas tri", "", number_of_points, 0, 10, 50, -25, 25);
      sigmas_tri_kin_fit[i] = new TH2F(name[i] + " sigmas kinfit", "", number_of_points, 0, 10, 50, -25, 25);
    }

  }

  Float_t gamma_path[4], cluster_time[4], mean_corr, corr;
  Int_t counts = 0, counts_all = 0;

  TCanvas *canvas_pulls[24];

  for (Int_t i = 0; i < 30; i++)
  {
    canvas[i] = new TCanvas(("Canvas" + std::to_string(i)).c_str(), "", 750, 750);

    if(i < 24)
      canvas_pulls[i] = new TCanvas(("Canvas_pulls" + std::to_string(i)).c_str(), "", 750, 750);
  }

  Float_t vec_before[4], vec_after[4], boost[3], ip_before[4], ip_after[4];

  for (Int_t i = 0; i < nentries; i++)
  {
    chain->GetEntry(i);
    if((mctruth == 1 || mctruth == 2)) counts_all++;

    if ((mctruth == 1 || mctruth == 2) && done_kinfit == 1 && chi2min < 40)
    {
      counts++;

      for(Int_t j = 0; j < 4; j++)
      {
        sigmas_init[j * 5] = pow(1.2,2);
        sigmas_init[j * 5 + 1] = pow(1.2,2);
        sigmas_init[j * 5 + 2] = pow(1.2 / sqrt((*X_init)(j * 5 + 4) / 1000.), 2);
        sigmas_init[j * 5 + 3] = pow(clu_time_error((*X_init)(j * 5 + 4)), 2);			// ns
        sigmas_init[j * 5 + 4] = pow(clu_ene_error((*X_init)(j * 5 + 4)), 2);			// ns
      }

      sigmas_init[20] = bhabha_mom_err[0];
      sigmas_init[21] = bhabha_mom_err[1];
      sigmas_init[22] = bhabha_mom_err[2];
      sigmas_init[23] = bhabha_mom_err[3];

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

      t_neutri = (Knetri_kinfit[8] - ip_kinfit[2]) / v_Kneutri;
      t_neumc = lengthneu_mc / v_Kneumc;

      cos_path_mom = ((Knetri_kinfit[6] - ip_kinfit[0])*Knetri_kinfit[0] + (Knetri_kinfit[7] - ip_kinfit[1])*Knetri_kinfit[1] + (Knetri_kinfit[8] - ip_kinfit[2])*Knetri_kinfit[2])/ (lengthneu_tri*sqrt(pow(Knetri_kinfit[0],2) + pow(Knetri_kinfit[1],2) + pow(Knetri_kinfit[2],2)));

      angle_path_mom = M_PI*acos(cos_path_mom)/180.;

      boost[0] = -(*X_min)(20)/(*X_min)(23);
      boost[1] = -(*X_min)(21)/(*X_min)(23);
      boost[2] = -(*X_min)(22)/(*X_min)(23);

      vec_before[0] = Knetri_kinfit[6];
      vec_before[1] = Knetri_kinfit[7];
      vec_before[2] = Knetri_kinfit[8];
      vec_before[3] = Knetri_kinfit[9];

      ip_before[0] = ip_kinfit[0];
      ip_before[1] = ip_kinfit[1];
      ip_before[2] = ip_kinfit[2];
      ip_before[3] = 0.;
      
      lorentz_transf(boost, vec_before, vec_after);
      lorentz_transf(boost, ip_before, ip_after);

      lengthneu_tri = sqrt(pow(vec_after[0] - ip_after[0], 2) + pow(vec_after[1] - ip_after[1], 2) + pow(vec_after[2] - ip_after[2], 2));

      //if((lengthneu_tri/(c_vel * vec_after[3])) > 0.25)
      neu_vtx_corr[3]->Fill(TMath::Prob(chi2min,M), Knetri_kinfit[9]);
      //else
      //neu_vtx_corr[3]->Fill(t_neumc, Knetri_kinfit[9]);
      // cluster_time[0] = Tcl[g4taken_kinfit[0]];
      // cluster_time[1] = Tcl[g4taken_kinfit[1]];
      // cluster_time[2] = Tcl[g4taken_kinfit[2]];
      // cluster_time[3] = Tcl[g4taken_kinfit[3]];

      // cluster_time[0] = Tcl[g4taken_kinfit[0]];
      // cluster_time[1] = Tcl[g4taken_kinfit[1]];
      // cluster_time[2] = Tcl[g4taken_kinfit[2]];
      // cluster_time[3] = Tcl[g4taken_kinfit[3]];

      dttri_hist->Fill(chi2min/(Float_t)M);
      prob_hist->Fill(lengthneu_tri/(c_vel * vec_after[3]));

      sigmas_std[0]->Fill(abs(Knemc[6] - ipmc[0]), Knetri_kinfit[6] - Knemc[6]);
      sigmas_std[1]->Fill(abs(Knemc[7] - ipmc[1]), Knetri_kinfit[7] - Knemc[7]);
      sigmas_std[2]->Fill(abs(Knemc[8] - ipmc[2]), Knetri_kinfit[8] - Knemc[8]);

      //if(lengthneu_tri/(c_vel * Knetri_kinfit[9]) > 0.9 && lengthneu_tri/(c_vel * Knetri_kinfit[9]) < 1.0)
      //sigmas_std[3]->Fill(TMath::Prob(chi2min,M), Knetri_kinfit[9]);
      //else
      sigmas_std[3]->Fill(t_neumc/tau_S_nonCPT, (Knetri_kinfit[9] - t_neumc)/tau_S_nonCPT);

      //sigmas_std[3]->Fill(t_neumc/tau_S_nonCPT, (Knetri_kinfit[9] - t_neumc)/tau_S_nonCPT);

      pulls[0]->Fill(Knetri_kinfit[0] - Knemc[0]);
      pulls[1]->Fill(Knetri_kinfit[1] - Knemc[1]);
      pulls[2]->Fill(Knetri_kinfit[2] - Knemc[2]);
      pulls[3]->Fill(Knetri_kinfit[3] - Knemc[3]);

      pulls[4]->Fill(Knetri_kinfit[6] - Knemc[6]);
      pulls[5]->Fill(Knetri_kinfit[7] - Knemc[7]);
      pulls[6]->Fill(Knetri_kinfit[8] - Knemc[8]);
      pulls[7]->Fill((Knetri_kinfit[9] - t_neumc)/tau_S_nonCPT);

      /*for(Int_t k = 0; k < 24; k++)
      {
        for(Int_t l = 0; l < 24; l++)
        {
          corr = (*cov_matrix)(k,l)/(sqrt((*cov_matrix)(k,k)) * sqrt((*cov_matrix)(l,l)));

          corr_matrix_plot->SetBinContent(corr_matrix_plot->GetBin(k,l), corr);
        }

        pulls[k]->Fill( ((*X_init)(k) - (*X_min)(k)) / sqrt(sigmas_init[k] - (*cov_matrix)(k,k)) );
      }*/
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

  x_title = Form("Prob(#chi^{2},%d)", M);
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
  dttri_hist->GetYaxis()->SetRangeUser(0.0, 1.2 * dttri_hist->GetMaximum());
  dttri_hist->Draw();
  c[4]->Print(id_canva + ".png");

  TF1 *func2 = new TF1("chi2dist", chi2dist, 0, 0.1, 2);

	func2->SetParameters(70000.0, (Double_t)M);
	func2->SetParNames("Norm", "Degrees of Freedom");

	func2->SetParLimits(0, 0.001, 100000.0);
	func2->SetParLimits(1, 1.0, 10.0);
	//prob_hist->Fit(func2);

  id_canva = "Canva" + std::to_string(5);
  c[5] = new TCanvas(id_canva, "", width, height);
  c[5]->SetRightMargin(0.15);
  c[5]->cd();
  //c[5]->SetLogy();

  prob_hist->GetYaxis()->SetMaxDigits(3);

  TString prob_title = Form("#beta^{CM}_{neu}");

  prob_hist->GetXaxis()->CenterTitle();
  prob_hist->GetYaxis()->CenterTitle();
  prob_hist->GetXaxis()->SetTitle(prob_title);
  prob_hist->GetYaxis()->SetTitle("Counts");
  prob_hist->GetYaxis()->SetRangeUser(0.0, 1.2 * prob_hist->GetMaximum());
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
                        "#beta^{CM}_{neu}"};//"t_{K#rightarrow#pi^{0}#pi^{0}}^{gen} [#tau_{S}]"};

  TString title_y[2] = {"#sigma(#vec{X}^{rec}_{K#rightarrow#pi^{0}#pi^{0},i} - #vec{X}^{gen}_{K#rightarrow#pi^{0}#pi^{0},i}) [cm]",
                        "t_{neu}^{tri} [ns]"};//"#sigma(t_{K#rightarrow#pi^{0}#pi^{0}}^{rec} - t_{K#rightarrow#pi^{0}#pi^{0}}^{gen}) [#tau_{S}]"};

  gStyle->SetOptStat(0);

  canvas[0]->cd();
  res_std_hist[0]->SetXTitle(title_x[0]);
  res_std_hist[0]->SetYTitle(title_y[0]);
  res_std_hist[0]->GetYaxis()->SetRangeUser(0.0, 8.0);
  res_std_hist[0]->Draw("PE1");
  res_std_hist[1]->Draw("PE1SAME");
  res_std_hist[2]->Draw("PE1SAME");
  legend->Draw();
  canvas[0]->Print(("sigmas_std" + std::to_string(0) + ".png").c_str());

  canvas[1]->cd();
  res_std_hist[3]->SetXTitle(title_x[1]);
  res_std_hist[3]->SetYTitle(title_y[1]);
  res_std_hist[3]->GetYaxis()->SetRangeUser(0.0, 3.0);
  res_std_hist[3]->Draw("PE1");
  canvas[1]->Print(("sigmas_std" + std::to_string(2) + ".png").c_str());

  TString fit_stats[3];
  TPaveText *fit_text = new TPaveText(0.7, 0.7, 0.9, 0.9, "NDC");

  TF1 *triple_fit = new TF1("triple_gaus", triple_gaus, -1500.0, 1500.0, 9, 1);
  triple_fit->SetParNames("Norm1", "Avg1", "Std1", "Norm2", "Avg2", "Std2", "Norm3", "Avg3", "Std3");

  TFitResultPtr result;

  for(Int_t i = 0; i < 8; i++)
  {
    canvas_pulls[i]->cd();

    //triple_fit->SetParameters(2E5, 0.0, 10.0, 2E5, 0.0, 10.0, 2E5, 0.0, 10.0);

    triple_fit->SetParLimits(0, 0.0, 5E6);
    triple_fit->SetParLimits(3, 0.0, 5E6);
    triple_fit->SetParLimits(6, 0.0, 5E6);
    triple_fit->SetParLimits(2, 0.001, 30.0);
    triple_fit->SetParLimits(5, 0.001, 30.0);
    triple_fit->SetParLimits(8, 0.001, 30.0);

    result = pulls[i]->Fit(triple_fit, "SF");

    std::cout << comb_std_dev(result->GetParams()) << std::endl;

    fit_stats[0] = Form("#chi^{2}/%.2d = %.2f", result->Ndf(), result->Chi2()/Double_t(result->Ndf()));
    fit_stats[1] = Form("Mean = %.2f#pm%.2f", result->Parameter(1), result->Error(1));
    fit_stats[2] = Form("Std dev = %.2f#pm%.2f", result->Parameter(2), result->Error(2));

    fit_text->AddText(fit_stats[0]);
    fit_text->AddText(fit_stats[1]);
    fit_text->AddText(fit_stats[2]);

    pulls[i]->GetYaxis()->CenterTitle(1);
    pulls[i]->SetLineWidth(5);

    pulls[i]->GetXaxis()->CenterTitle(1);

    pulls[i]->GetYaxis()->SetMaxDigits(3);
    pulls[i]->GetYaxis()->CenterTitle(1);

    pulls[i]->Draw();
    fit_text->Draw();

    canvas_pulls[i]->Print(("pulls" + std::to_string(i + 1) +".png").c_str());

    fit_text->Clear();
  }

  std::cout << "Events all: " << counts_all << " " << "Events: " << counts << std::endl;
  std::cout << "Efficiency: " << counts/(Float_t)counts_all << std::endl;

  return 0;
}
