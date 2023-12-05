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
const TString ext = ".png";

int comp_of_methods(Int_t file_num = 0, Double_t cut_prob = 0.0)
{
  TChain *chain = new TChain("INTERF/h1");
  chain_init(chain, 1, 56);

  TString file_name[3];

  file_name[0] = "neuvtx_tri_kin_fit_1_56_100_5_no_time.root";
  file_name[1] = "neuvtx_tri_kin_fit_1_56_30_5_time.root";
  file_name[2] = "neuvtx_tri_rec_1_56.root";

  TFile *file = new TFile(file_name[file_num]);
  TTree *tree;

  if(file_num == 0 || file_num == 1)
    tree = (TTree *)file->Get("h_tri_kin_fit");
  else
    tree = (TTree *)file->Get("h_tri");

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
  Int_t done_kinfit, g4taken_kinfit[4], bunchnum;
  TMatrixD *cov_matrix = new TMatrixD(27, 27);
  TVectorD *min_const = new TVectorD(M);
  TVectorD *lag_mult = new TVectorD(M), *X_min = new TVectorD(27), *X_init = new TVectorD(27);

  if(file_num == 0 || file_num == 1)
  {
    tree->SetBranchAddress("fourgamma1tri_kinfit", gamma_kinfit[0]);
    tree->SetBranchAddress("fourgamma2tri_kinfit", gamma_kinfit[1]);
    tree->SetBranchAddress("fourgamma3tri_kinfit", gamma_kinfit[2]);
    tree->SetBranchAddress("fourgamma4tri_kinfit", gamma_kinfit[3]);

    tree->SetBranchAddress("fourKnetri_kinfit", Knetri_kinfit);

    tree->SetBranchAddress("iptri_kinfit", ip_kinfit);
    tree->SetBranchAddress("done4_kinfit", &done_kinfit);

    if(file_num == 1)
      tree->SetBranchAddress("bunchnum", &bunchnum);

    tree->SetBranchAddress("chi2min", &chi2min);
    tree->SetBranchAddress("min_cov", &cov_matrix);
    tree->SetBranchAddress("const_min", &min_const);
    tree->SetBranchAddress("lag_mult", &lag_mult);

    tree->SetBranchAddress("min_vars", &X_min);
    tree->SetBranchAddress("init_vars", &X_init);
  }
  else
  {
    tree->SetBranchAddress("fourgamma1tri", gamma_kinfit[0]);
    tree->SetBranchAddress("fourgamma2tri", gamma_kinfit[1]);
    tree->SetBranchAddress("fourgamma3tri", gamma_kinfit[2]);
    tree->SetBranchAddress("fourgamma4tri", gamma_kinfit[3]);

    tree->SetBranchAddress("fourKnetri", Knetri_kinfit);

    tree->SetBranchAddress("iptri", ip_kinfit);
    tree->SetBranchAddress("done4", &done_kinfit);
  }

  chain->AddFriend(tree);

  UInt_t nentries = (UInt_t)chain->GetEntries();

  Float_t lengthneu_mc, lengthneu_tri, lengthneu_tri_CM, lengthneu_rec, lengthch_mc, lengthch_rec,
      v_Kneutri, v_Kneutri_CM, v_Kchrec, v_Kneurec, v_Kneumc, v_Kchmc,
      t_chmc, t_neumc, t_chrec, t_neutri, t_neurec, Rtneu_rec, Rtneu_mc, Rtneu_tri,
      angle_path_mom, cos_path_mom, sigmas_init[24];

  TString id_hist, id_canva;

  TCanvas *canvas[2][30], *canvas_pulls[2][24];
  Int_t width = 750, height = 750;

  TH1 *chi2_hist[2], *prob_hist[2];
  TH1 *beta_hist[2];
  TH2 *neu_vtx_corr[2][4], *neu_mom[2][4], *ip_coor[2][3];
  TH2 *beta_time[2][2];

  TH2 *sigmas_std[2][5], *sigmas_tri[2][5], *sigmas_tri_kin_fit[2][5];
  TH1 *neu_std_hist[2][5], *neu_tri_hist[2][5], *neu_tri_kin_fit_hist[2][5];
  TH1 *res_std_hist[2][5], *res_tri_hist[2][5], *res_tri_kin_fit_hist[2][5];
  TH1 *pulls[2][8];
  TH1 *bunch[2];

  TString name[5] = {"x", "y", "z", "time", "length"};

  const UInt_t number_of_points = 10;

  for (Int_t j = 0; j < 2; j++)
  {
    for (Int_t i = 0; i < 30; i++)
    {
      id_canva = "Canva" + std::to_string(j) + std::to_string(i);
      canvas[j][i] = new TCanvas(id_canva, "", width, height);

      if (i < 24)
      {
        id_canva = "Canva_pull" + std::to_string(j) + std::to_string(i);
        canvas_pulls[j][i] = new TCanvas(id_canva, "", width, height);
      }
    }

    for (Int_t i = 0; i < 3; i++)
    {
      id_hist = "Coor" + std::to_string(j) + std::to_string(i);
      neu_vtx_corr[j][i] = new TH2F(id_hist, "", 100, -50, 50, 100, -50, 50);

      id_hist = "Mom" + std::to_string(j) + std::to_string(i);
      neu_mom[j][i] = new TH2F(id_hist, "", 200, -300, 300, 200, -300, 300);

      id_hist = "IP coor" + std::to_string(j) + std::to_string(i);
      ip_coor[j][i] = new TH2F(id_hist, "", 100, -10, 10, 100, -10, 10);
    }

    id_hist = "Coor" + std::to_string(j) + std::to_string(3);
    neu_vtx_corr[j][3] = new TH2F(id_hist, "", 100, 0.0, 20.0, 100, 0.0, 20.0);

    id_hist = "Mom" + std::to_string(j) + std::to_string(3);
    neu_mom[j][3] = new TH2F(id_hist, "", 100, 504, 520, 100, 507, 512);

    for (Int_t i = 0; i < 8; i++)
    {
      id_hist = "Var" + std::to_string(j) + std::to_string(i);
      if (i == 0)
        pulls[j][i] = new TH1F(id_hist, ";#vec{p}^{tri}_{neu,1} - #vec{p}^{gen}_{neu,1} [MeV/c];Counts", 101, -400, 400);
      else if (i == 1)
        pulls[j][i] = new TH1F(id_hist, ";#vec{p}^{tri}_{neu,2} - #vec{p}^{gen}_{neu,2} [MeV/c];Counts", 101, -400, 400);
      else if (i == 2)
        pulls[j][i] = new TH1F(id_hist, ";#vec{p}^{tri}_{neu,3} - #vec{p}^{gen}_{neu,3} [MeV/c];Counts", 101, -400, 400);
      else if (i == 3)
        pulls[j][i] = new TH1F(id_hist, ";E^{tri}_{neu} - E^{gen}_{neu,3} [MeV];Counts", 201, -100, 100);
      else if (i == 4)
        pulls[j][i] = new TH1F(id_hist, ";#vec{X}^{tri}_{neu,1} - #vec{X}^{gen}_{neu,1} [cm];Counts", 201, -100, 100);
      else if (i == 5)
        pulls[j][i] = new TH1F(id_hist, ";#vec{X}^{tri}_{neu,2} - #vec{X}^{gen}_{neu,2} [cm];Counts", 201, -100, 100);
      else if (i == 6)
        pulls[j][i] = new TH1F(id_hist, ";#vec{X}^{tri}_{neu,3} - #vec{X}^{gen}_{neu,3} [cm];Counts", 201, -100, 100);
      else if (i == 7)
        pulls[j][i] = new TH1F(id_hist, ";t^{tri}_{neu} - t^{gen}_{neu} [#tau_{S}];Counts", 101, -10, 10);
    }

    id_hist = "Chi2" + std::to_string(j);
    chi2_hist[j] = new TH1F(id_hist, "", 200, 0.0, 100.0);

    id_hist = "Prob" + std::to_string(j);
    prob_hist[j] = new TH1F(id_hist, "", 100, 0.0, 1.0);

    id_hist = "Beta" + std::to_string(j);
    beta_hist[j] = new TH1F(id_hist, ";#beta^{CM}_{K#rightarrow#pi^{0}#pi^{0}};Counts", 100, 0.1, 0.4);

    id_hist = "Bunch_num" + std::to_string(j);
    bunch[j] = new TH1I(id_hist, ";Number of bunch correction;Counts", 101, -10, 10);

    for (Int_t i = 0; i < 5; i++)
    {
      id_hist = Form(name[i] + "%d%d", j + 1, i + 1);
      if (i < 2)
      {
        neu_std_hist[j][i] = new TH1F(id_hist + " std coor", "", 100, -200, 200);
        neu_tri_hist[j][i] = new TH1F(id_hist + " tri coor", "", 100, -200, 200);
        neu_tri_kin_fit_hist[j][i] = new TH1F(id_hist + " kinfit coor", "", 100, -165, 165);
        res_std_hist[j][i] = new TH1F(id_hist + " std res", "", number_of_points, 0, 50);
        res_tri_hist[j][i] = new TH1F(id_hist + " tri res", "", number_of_points, 0, 50);
        res_tri_kin_fit_hist[j][i] = new TH1F(id_hist + " kinfit res", "", number_of_points, 0, 50);
        sigmas_std[j][i] = new TH2F(id_hist + " sigmas std", "", number_of_points, 0, 50, 25, -10, 10);
        sigmas_tri[j][i] = new TH2F(id_hist + " sigmas tri", "", number_of_points, 0, 50, 25, -10, 10);
        sigmas_tri_kin_fit[j][i] = new TH2F(id_hist + " sigmas kinfit", "", number_of_points, 0, 50, 25, -10, 10);
      }
      else if (i == 2)
      {
        neu_std_hist[j][i] = new TH1F(id_hist + " std", "", 100, 0, 5);
        neu_tri_hist[j][i] = new TH1F(id_hist + " tri", "", 100, 0, 5);
        neu_tri_kin_fit_hist[j][i] = new TH1F(id_hist + " kinfit", "", 100, 0, 5);
        res_std_hist[j][i] = new TH1F(id_hist + " std res", "", number_of_points, 0, 50);
        res_tri_hist[j][i] = new TH1F(id_hist + " tri res", "", number_of_points, 0, 50);
        res_tri_kin_fit_hist[j][i] = new TH1F(id_hist + " kinfit res", "", number_of_points, 0, 50);
        sigmas_std[j][i] = new TH2F(id_hist + " sigmas std", "", number_of_points, 0, 50, 25, -10, 10);
        sigmas_tri[j][i] = new TH2F(id_hist + " sigmas tri", "", number_of_points, 0, 50, 25, -10, 10);
        sigmas_tri_kin_fit[j][i] = new TH2F(id_hist + " sigmas kinfit", "", number_of_points, 0, 50, 25, -10, 10);
      }
      else if (i == 3)
      {
        neu_std_hist[j][i] = new TH1F(id_hist + " std", "", 100, 0, 5);
        neu_tri_hist[j][i] = new TH1F(id_hist + " tri", "", 100, 0, 5);
        neu_tri_kin_fit_hist[j][i] = new TH1F(id_hist + " kinfit", "", 100, 0, 5);
        res_std_hist[j][i] = new TH1F(id_hist + " std res", "", number_of_points, 0, 50);
        res_tri_hist[j][i] = new TH1F(id_hist + " tri res", "", number_of_points, 0, 50);
        res_tri_kin_fit_hist[j][i] = new TH1F(id_hist + " kinfit res", "", number_of_points, 0, 50);
        sigmas_std[j][i] = new TH2F(id_hist + " sigmas std", "", number_of_points, 0, 50, 50, -25, 25);
        sigmas_tri[j][i] = new TH2F(id_hist + " sigmas tri", "", number_of_points, 0, 50, 50, -25, 25);
        sigmas_tri_kin_fit[j][i] = new TH2F(id_hist + " sigmas kinfit", "", number_of_points, 0, 50, 50, -25, 25);
      }
      else if (i == 4)
      {
        neu_std_hist[j][i] = new TH1F(id_hist + " std", "", 100, 0, 5);
        neu_tri_hist[j][i] = new TH1F(id_hist + " tri", "", 100, 0, 5);
        neu_tri_kin_fit_hist[j][i] = new TH1F(id_hist + " kinfit", "", 100, 0, 5);
        res_std_hist[j][i] = new TH1F(id_hist + " std res", "", number_of_points, 0, 50);
        res_tri_hist[j][i] = new TH1F(id_hist + " tri res", "", number_of_points, 0, 50);
        res_tri_kin_fit_hist[j][i] = new TH1F(id_hist + " kinfit res", "", number_of_points, 0, 50);
        sigmas_std[j][i] = new TH2F(id_hist + " sigmas std", "", number_of_points, 0, 50, 50, -25, 25);
        sigmas_tri[j][i] = new TH2F(id_hist + " sigmas tri", "", number_of_points, 0, 50, 50, -25, 25);
        sigmas_tri_kin_fit[j][i] = new TH2F(id_hist + " sigmas kinfit", "", number_of_points, 0, 50, 50, -25, 25);
      }
    }
  }

  Int_t counts = 0, counts_all = 0;
  Float_t vec_before[4], vec_after[4], boost[3], ip_before[4], ip_after[4];

  Bool_t cut = true;

  for (Int_t i = 0; i < nentries; i++)
  {
    chain->GetEntry(i);

    if ((mctruth == 1 || mctruth == 2))
    {
      counts_all++;

      // Kaon path length
      lengthch_mc = sqrt(pow(Kchmc[6] - ipmc[0], 2) + pow(Kchmc[7] - ipmc[1], 2) + pow(Kchmc[8] - ipmc[2], 2));
      lengthch_rec = sqrt(pow(Kchboost[6] - ip[0], 2) + pow(Kchboost[7] - ip[1], 2) + pow(Kchboost[8] - ip[2], 2));

      lengthneu_mc = sqrt(pow(Knemc[6] - ipmc[0], 2) + pow(Knemc[7] - ipmc[1], 2) + pow(Knemc[8] - ipmc[2], 2));
      lengthneu_tri = sqrt(pow(Knetri_kinfit[6] - ip_kinfit[0], 2) + pow(Knetri_kinfit[7] - ip_kinfit[1], 2) + pow(Knetri_kinfit[8] - ip_kinfit[2], 2));
      lengthneu_rec = sqrt(pow(Knereclor[6] - ip[0], 2) + pow(Knereclor[7] - ip[1], 2) + pow(Knereclor[8] - ip[2], 2));
      //
      // Kaon transverse radius
      Rtneu_mc = sqrt(pow(Knemc[6] - ipmc[0], 2) + pow(Knemc[7] - ipmc[1], 2));
      Rtneu_tri = sqrt(pow(Knetri_kinfit[6] - ip_kinfit[0], 2) + pow(Knetri_kinfit[7] - ip_kinfit[1], 2));
      Rtneu_rec = sqrt(pow(Knereclor[6] - ip[0], 2) + pow(Knereclor[7] - ip[1], 2));
      //
      // Kaon velocity
      v_Kchmc = c_vel * Kchmc[4] / Kchmc[3];
      v_Kchrec = c_vel * Kchboost[4] / Kchboost[3];

      v_Kneumc = c_vel * Knemc[4] / Knemc[3];
      v_Kneutri = c_vel * Knetri_kinfit[4] / Knetri_kinfit[3];
      v_Kneurec = c_vel * Knereclor[4] / Knereclor[3];
      //
      // Kaon flight times
      t_chrec = lengthch_rec / v_Kchrec;
      t_chmc = lengthch_mc / v_Kchmc;

      t_neurec = lengthneu_rec / v_Kneurec;
      t_neutri = lengthneu_tri / v_Kneutri;
      t_neumc = lengthneu_mc / v_Kneumc;
      //
      // Angle between vector between IP and neuvtx and Kaon's momentum vector
      cos_path_mom = ((Knetri_kinfit[6] - ip_kinfit[0]) * Knetri_kinfit[0] + (Knetri_kinfit[7] - ip_kinfit[1]) * Knetri_kinfit[1] + (Knetri_kinfit[8] - ip_kinfit[2]) * Knetri_kinfit[2]) / (lengthneu_tri * sqrt(pow(Knetri_kinfit[0], 2) + pow(Knetri_kinfit[1], 2) + pow(Knetri_kinfit[2], 2)));

      angle_path_mom = M_PI * acos(cos_path_mom) / 180.;
      //
      //! Lorentz transformation to get the beta in Phi's CM
      if(file_num == 0 || file_num == 1)
      {
        boost[0] = -(*X_min)(20) / (*X_min)(23);
        boost[1] = -(*X_min)(21) / (*X_min)(23);
        boost[2] = -(*X_min)(22) / (*X_min)(23);
      }
      else
      {
        boost[0] = -phi_mom[0] / phi_mom[3];
        boost[1] = -phi_mom[1] / phi_mom[3];
        boost[2] = -phi_mom[2] / phi_mom[3];
      }

      vec_before[0] = Knetri_kinfit[6] - ip_kinfit[0];
      vec_before[1] = Knetri_kinfit[7] - ip_kinfit[1];
      vec_before[2] = Knetri_kinfit[8] - ip_kinfit[2];
      vec_before[3] = Knetri_kinfit[9];

      lorentz_transf(boost, vec_before, vec_after);

      lengthneu_tri_CM = sqrt(pow(vec_after[0], 2) + pow(vec_after[1], 2) + pow(vec_after[2], 2));
      v_Kneutri_CM = lengthneu_tri_CM / vec_after[3];
      //!

      for (Int_t j = 0; j < 2; j++)
      {
        if (j == 0)
          cut = TMath::Prob(chi2min, M) > 0.0;
        else
          cut = TMath::Prob(chi2min, M) > cut_prob;

        if (cut && done_kinfit == 1)
        {
          if (j != 0)
            counts++;

          chi2_hist[j]->Fill(chi2min);
          prob_hist[j]->Fill(TMath::Prob(chi2min, M));

          beta_hist[j]->Fill(v_Kneutri_CM / c_vel);

          bunch[j]->Fill(bunchnum);

          for (Int_t k = 0; k < 3; k++)
            ip_coor[j][k]->Fill(ipmc[k], ip_kinfit[k]);

          for (Int_t k = 0; k < 5; k++)
          {
            if (k < 3)
            {
              if(abs(Knemc[6 + k] - Knetri_kinfit[6 + k]) < 10000)
                neu_vtx_corr[j][k]->Fill(Knemc[6 + k], Knetri_kinfit[6 + k]);

              sigmas_std[j][k]->Fill(abs(Knemc[6 + k] - ipmc[k]), Knetri_kinfit[6 + k] - Knemc[6 + k]);
              pulls[j][4 + k]->Fill(Knetri_kinfit[6 + k] - Knemc[6 + k]);
            }
            else if (k == 3)
            {
              if(abs(Knemc[6] - Knetri_kinfit[6]) < 10000 && abs(Knemc[7] - Knetri_kinfit[7]) < 10000 && abs(Knemc[8] - Knetri_kinfit[8]) < 10000)
                neu_vtx_corr[j][3]->Fill(t_neumc, Knetri_kinfit[6 + k]);

              sigmas_std[j][3]->Fill(lengthneu_mc, (Knetri_kinfit[9] - t_neumc) / tau_S_nonCPT);
              pulls[j][4 + k]->Fill((Knetri_kinfit[6 + k] - t_neumc) / tau_S_nonCPT);
            }
            else
            {
                sigmas_std[j][4]->Fill(lengthneu_mc, (lengthneu_tri - lengthneu_mc));
            }

            neu_mom[j][k]->Fill(Knemc[k], Knetri_kinfit[k]);

            pulls[j][k]->Fill(Knetri_kinfit[k] - Knemc[k]);
            pulls[j][4 + k]->Fill(Knetri_kinfit[k] - Knemc[k]);
          }
        }
      }
    }
  }

  TString fit_stats[3];
  TPaveText *fit_text = new TPaveText(0.7, 0.7, 0.9, 0.9, "NDC");

  TF1 *triple_fit = new TF1("triple_gaus", triple_gaus, -400.0, 400.0, 9, 1);
  triple_fit->SetParNames("Norm1", "Avg1", "Std1", "Norm2", "Avg2", "Std2", "Norm3", "Avg3", "Std3");
  triple_fit->SetLineWidth(4);

  TFitResultPtr result;

  TString x_title, y_title;

  Float_t parameter[3];

  TF1 *func = new TF1("f1", "[0]*exp(-0.5*pow((x-[1])/[2],2))", -200, 200);
  func->SetParNames("Constant", "Mean", "Sigma");

  Color_t res_color[5] = {kRed, kBlue, kGreen, kBlack, kBlack};
  TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9, "Axis");

  for (Int_t j = 0; j < 2; j++)
  {
    gStyle->SetOptStat("iMROU");
    gStyle->SetStatX(0.85);
    gStyle->SetStatY(0.9);

    for (Int_t i = 0; i < 3; i++)
    {
      canvas[j][i]->SetRightMargin(0.15);
      canvas[j][i]->SetLeftMargin(0.17);
      canvas[j][i]->cd();

      id_canva = "coor" + std::to_string(j + 1) + std::to_string(i + 1);
      x_title = Form("#vec{X}_{neu,%d}^{gen} [cm]", i + 1);
      y_title = Form("#vec{X}_{neu,%d}^{tri} [cm]", i + 1);

      neu_vtx_corr[j][i]->GetXaxis()->SetTitle(x_title);
      neu_vtx_corr[j][i]->GetYaxis()->SetTitle(y_title);
      neu_vtx_corr[j][i]->Draw("COLZ");

      canvas[j][i]->Print(id_canva + ext);
      //!
      canvas[j][i + 4]->SetRightMargin(0.15);
      canvas[j][i + 4]->SetLeftMargin(0.17);
      canvas[j][i + 4]->cd();

      id_canva = "ip_coor" + std::to_string(j + 1) + std::to_string(i + 1);
      x_title = Form("#vec{X}_{IP,%d}^{gen} [cm]", i + 1);
      y_title = Form("#vec{X}_{IP,%d}^{tri} [cm]", i + 1);

      ip_coor[j][i]->GetXaxis()->SetTitle(x_title);
      ip_coor[j][i]->GetYaxis()->SetTitle(y_title);
      ip_coor[j][i]->Draw("COLZ");

      canvas[j][i + 4]->Print(id_canva + ext);
      //!
      canvas[j][i + 7]->SetRightMargin(0.15);
      canvas[j][i + 7]->SetLeftMargin(0.17);
      canvas[j][i + 7]->cd();

      id_canva = "mom" + std::to_string(j + 1) + std::to_string(i + 1);
      x_title = Form("#vec{p}_{K#rightarrow#pi^{0}#pi^{0},%d}^{gen} [MeV/c]", i + 1);
      y_title = Form("#vec{p}_{K#rightarrow#pi^{0}#pi^{0},%d}^{tri} [MeV/c]", i + 1);

      neu_mom[j][i]->GetXaxis()->SetTitle(x_title);
      neu_mom[j][i]->GetYaxis()->SetTitle(y_title);
      neu_mom[j][i]->Draw("COLZ");

      canvas[j][i + 7]->Print(id_canva + ext);
    }
    //!
    canvas[j][3]->SetRightMargin(0.15);
    canvas[j][3]->SetLeftMargin(0.17);
    canvas[j][3]->cd();

    id_canva = "time_neutral" + std::to_string(j + 1);
    x_title = Form("t_{neu}^{gen} [ns]");
    y_title = Form("t_{neu}^{tri} [ns]");

    neu_vtx_corr[j][3]->GetXaxis()->SetTitle(x_title);
    neu_vtx_corr[j][3]->GetYaxis()->SetTitle(y_title);
    neu_vtx_corr[j][3]->Draw("COLZ");

    canvas[j][3]->Print(id_canva + ext);
    //!
    canvas[j][10]->SetRightMargin(0.15);
    canvas[j][10]->SetLeftMargin(0.17);
    canvas[j][10]->cd();

    id_canva = "energy" + std::to_string(j + 1);
    x_title = Form("E_{K#rightarrow#pi^{0}#pi^{0}}^{gen} [MeV]");
    y_title = Form("E_{K#rightarrow#pi^{0}#pi^{0}}^{tri} [MeV]");

    neu_mom[j][3]->GetXaxis()->SetTitle(x_title);
    neu_mom[j][3]->GetYaxis()->SetTitle(y_title);
    neu_mom[j][3]->Draw("COLZ");

    canvas[j][10]->Print(id_canva + ext);
    //!

    gStyle->SetOptStat(0);
    for (Int_t i = 0; i < 8; i++)
    {
      canvas[j][i + 13]->cd();

      triple_fit->SetParameters(pulls[j][i]->Integral(), pulls[j][i]->GetMean(), pulls[j][i]->GetStdDev(), pulls[j][i]->Integral(), pulls[j][i]->GetMean(), pulls[j][i]->GetStdDev(), pulls[j][i]->Integral(), pulls[j][i]->GetMean(), pulls[j][i]->GetStdDev());

      if (i < 3)
      {
        triple_fit->SetParLimits(0, 0.01 * pulls[j][i]->Integral(), 1000.0 * pulls[j][i]->Integral());
        triple_fit->SetParLimits(3, 0.01 * pulls[j][i]->Integral(), 1000.0 * pulls[j][i]->Integral());
        triple_fit->SetParLimits(6, 0.01 * pulls[j][i]->Integral(), 1000.0 * pulls[j][i]->Integral());

        triple_fit->SetParLimits(2, 0.01 * pulls[j][i]->GetStdDev(), 1000.0 * pulls[j][i]->GetStdDev());
        triple_fit->SetParLimits(5, 0.01 * pulls[j][i]->GetStdDev(), 1000.0 * pulls[j][i]->GetStdDev());
        triple_fit->SetParLimits(8, 0.01 * pulls[j][i]->GetStdDev(), 1000.0 * pulls[j][i]->GetStdDev());
      }
      else if (i >= 3 && i < 7)
      {
        triple_fit->SetParLimits(0, 0.01 * pulls[j][i]->Integral(), 1000.0 * pulls[j][i]->Integral());
        triple_fit->SetParLimits(3, 0.01 * pulls[j][i]->Integral(), 1000.0 * pulls[j][i]->Integral());
        triple_fit->SetParLimits(6, 0.01 * pulls[j][i]->Integral(), 1000.0 * pulls[j][i]->Integral());

        triple_fit->SetParLimits(2, 0.01 * pulls[j][i]->GetStdDev(), 100.0 * pulls[j][i]->GetStdDev());
        triple_fit->SetParLimits(5, 0.01 * pulls[j][i]->GetStdDev(), 100.0 * pulls[j][i]->GetStdDev());
        triple_fit->SetParLimits(8, 0.01 * pulls[j][i]->GetStdDev(), 100.0 * pulls[j][i]->GetStdDev());
      }
      else
      {
        triple_fit->SetParLimits(0, 0.01 * pulls[j][i]->Integral(), 100.0 * pulls[j][i]->Integral());
        triple_fit->SetParLimits(3, 0.01 * pulls[j][i]->Integral(), 100.0 * pulls[j][i]->Integral());
        triple_fit->SetParLimits(6, 0.01 * pulls[j][i]->Integral(), 100.0 * pulls[j][i]->Integral());

        triple_fit->SetParLimits(2, 0.01 * pulls[j][i]->GetStdDev(), 5.0 * pulls[j][i]->GetStdDev());
        triple_fit->SetParLimits(5, 0.01 * pulls[j][i]->GetStdDev(), 5.0 * pulls[j][i]->GetStdDev());
        triple_fit->SetParLimits(8, 0.01 * pulls[j][i]->GetStdDev(), 5.0 * pulls[j][i]->GetStdDev());
      }

      result = pulls[j][i]->Fit(triple_fit, "SF");

      fit_stats[2] = Form("Width = %.2f", comb_std_dev(result->GetParams()));

      fit_text->AddText(fit_stats[2]);

      pulls[j][i]->SetLineWidth(5);

      pulls[j][i]->GetXaxis()->CenterTitle(1);

      pulls[j][i]->GetYaxis()->SetMaxDigits(3);
      pulls[j][i]->GetYaxis()->CenterTitle(1);
      pulls[j][i]->GetYaxis()->SetRangeUser(0.0, 1.3 * pulls[j][i]->GetMaximum());

      pulls[j][i]->Draw();
      fit_text->Draw();

      id_canva = "pull" + std::to_string(j + 1) + std::to_string(i + 1);
      canvas[j][i + 13]->Print(id_canva + ext);

      fit_text->Clear();
    }
    //!
    canvas[j][21]->cd();
    bunch[j]->GetYaxis()->SetRangeUser(0.0, 1.3 * bunch[j]->GetMaximum());
    bunch[j]->SetLineColor(kBlack);
    bunch[j]->Draw();

    id_canva = "bunchnum" + std::to_string(j + 1);
    canvas[j][21]->Print(id_canva + ext);
    //!
    for (Int_t k = 0; k < 5; k++)
    {
      for (Int_t i = 1; i <= number_of_points; i++)
      {
        parameter[0] = sigmas_std[j][k]->ProjectionY("_py", i, i)->Integral();
        parameter[1] = sigmas_std[j][k]->ProjectionY("_py", i, i)->GetBinCenter(sigmas_std[j][k]->ProjectionY("_py", i, i)->GetMaximumBin());
        parameter[2] = sigmas_std[j][k]->ProjectionY("_py", i, i)->GetStdDev();

        func->SetParameters(parameter[0], parameter[1], parameter[2]);

        func->SetParLimits(2, 0.0001, 100.);

        result = sigmas_std[j][k]->ProjectionY("_py", i, i)->Fit("f1", "S");
        res_std_hist[j][k]->SetBinContent(i, result->Parameter(2));
        res_std_hist[j][k]->SetBinError(i, result->ParError(2));
      }

      res_std_hist[j][k]->SetLineColor(res_color[k]);
    }

    legend->AddEntry(res_std_hist[j][0], "x", "le");
    legend->AddEntry(res_std_hist[j][1], "y", "le");
    legend->AddEntry(res_std_hist[j][2], "z", "le");
    //!
    id_canva = "sigmas_std_coordinates" + std::to_string(j + 1);
    x_title = "|#vec{X}^{gen}_{K#rightarrow#pi^{0}#pi^{0},i} - #vec{X}^{gen}_{#Phi,i}| [cm]";
    y_title = "#sigma(#vec{X}^{rec}_{K#rightarrow#pi^{0}#pi^{0},i} - #vec{X}^{gen}_{K#rightarrow#pi^{0}#pi^{0},i}) [cm]";

    canvas[j][22]->cd();
    res_std_hist[j][0]->SetXTitle(x_title);
    res_std_hist[j][0]->SetYTitle(y_title);
    res_std_hist[j][0]->GetYaxis()->SetRangeUser(0.0, 15.0);
    res_std_hist[j][0]->Draw("PE1");
    res_std_hist[j][1]->Draw("PE1SAME");
    res_std_hist[j][2]->Draw("PE1SAME");
    legend->Draw();
    canvas[j][22]->Print(id_canva + ext);
    //!
    id_canva = "sigmas_std_times" + std::to_string(j + 1);
    x_title = "K#rightarrow#pi^{0}#pi^{0} path generated [cm]";
    y_title = "#sigma(t_{K#rightarrow#pi^{0}#pi^{0}}^{rec} - t_{K#rightarrow#pi^{0}#pi^{0}}^{gen}) [#tau_{S}]";
    canvas[j][23]->cd();
    res_std_hist[j][3]->SetXTitle(x_title);
    res_std_hist[j][3]->SetYTitle(y_title);
    res_std_hist[j][3]->GetYaxis()->SetRangeUser(0.0, 6);
    res_std_hist[j][3]->Draw("PE1");
    canvas[j][23]->Print(id_canva + ext);
    //!
    id_canva = "sigmas_std_length" + std::to_string(j + 1);
    x_title = "K#rightarrow#pi^{0}#pi^{0} path generated [cm]";
    y_title = "#sigma(l_{K#rightarrow#pi^{0}#pi^{0}}^{rec} - l_{K#rightarrow#pi^{0}#pi^{0}}^{gen}) [#tau_{S}]";
    canvas[j][24]->cd();
    res_std_hist[j][4]->SetXTitle(x_title);
    res_std_hist[j][4]->SetYTitle(y_title);
    res_std_hist[j][4]->GetYaxis()->SetRangeUser(0.0, 40.0);
    res_std_hist[j][4]->Draw("PE1");
    canvas[j][24]->Print(id_canva + ext);
    //!

    canvas[j][25]->SetRightMargin(0.15);
    canvas[j][25]->SetLeftMargin(0.17);
    canvas[j][25]->cd();

    beta_hist[j]->GetYaxis()->SetMaxDigits(3);

    id_canva = "beta_CM" + std::to_string(j + 1);
    x_title = "#beta^{CM}_{K#rightarrow#pi^{0}#pi^{0}}";
    y_title = "Counts";

    beta_hist[j]->GetXaxis()->CenterTitle();
    beta_hist[j]->GetYaxis()->CenterTitle();
    beta_hist[j]->GetXaxis()->SetTitle(x_title);
    beta_hist[j]->GetYaxis()->SetTitle(y_title);
    beta_hist[j]->GetYaxis()->SetRangeUser(0.0, 1.3 * beta_hist[j]->GetMaximum());
    beta_hist[j]->SetLineColor(kBlack);
    beta_hist[j]->Draw();
    canvas[j][25]->Print(id_canva + ext);

    legend->Clear();
  }

  //! Chi2 and Prob histos

  canvas[0][11]->SetRightMargin(0.15);
  canvas[0][11]->SetLeftMargin(0.17);
  canvas[0][11]->cd();

  chi2_hist[0]->GetYaxis()->SetMaxDigits(3);

  id_canva = "chi2_dist";
  x_title = "#chi^{2}_{tri}";
  y_title = "Counts";

  chi2_hist[0]->GetXaxis()->CenterTitle();
  chi2_hist[0]->GetYaxis()->CenterTitle();
  chi2_hist[0]->GetXaxis()->SetTitle(x_title);
  chi2_hist[0]->GetYaxis()->SetTitle(y_title);
  chi2_hist[0]->GetYaxis()->SetRangeUser(0.0, 1.3 * chi2_hist[0]->GetMaximum());
  chi2_hist[0]->SetLineColor(kBlack);
  chi2_hist[1]->SetLineColor(kRed);
  chi2_hist[0]->Draw();
  chi2_hist[1]->Draw("SAME");
  canvas[0][11]->Print(id_canva + ext);
  //!
  canvas[0][12]->SetRightMargin(0.15);
  canvas[0][12]->SetLeftMargin(0.17);
  canvas[0][12]->cd();
  canvas[0][12]->SetLogy();

  prob_hist[0]->GetYaxis()->SetMaxDigits(3);

  id_canva = "prob";
  x_title = Form("Prob(#chi^{2}_{tri},%d)", M);
  y_title = "Counts";

  prob_hist[0]->GetXaxis()->CenterTitle();
  prob_hist[0]->GetYaxis()->CenterTitle();
  prob_hist[0]->GetXaxis()->SetTitle(x_title);
  prob_hist[0]->GetYaxis()->SetTitle("Counts");
  if (canvas[0][12]->GetLogy())
    prob_hist[0]->GetYaxis()->SetRangeUser(1.0, 10. * prob_hist[0]->GetMaximum());
  else
    prob_hist[0]->GetYaxis()->SetRangeUser(0.0, 1.3 * prob_hist[0]->GetMaximum());
  prob_hist[0]->SetLineColor(kBlack);
  prob_hist[1]->SetLineColor(kRed);
  prob_hist[0]->Draw();
  prob_hist[1]->Draw("SAME");
  canvas[0][12]->Print(id_canva + ext);
  //!

  std::cout << "Events all: " << counts_all << " "
            << "Events: " << counts << std::endl;
  std::cout << "Efficiency: " << counts / (Float_t)counts_all << std::endl;

  return 0;
}
