#include <iostream>
#include <fstream>

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

#include "src/chain_init.C"
#include "../../Include/const.h"
#include "../../Include/Codes/kloe_class.h"
#include "../../Include/Codes/interf_function.h"
#include "../../Include/Codes/chi2_dist.h"
#include "../../Include/Codes/uncertainties.h"
#include "../../Include/Codes/lorentz_transf.h"
#include "../../Include/Codes/triple_gaus.h"
#include "../../Include/Codes/arrays_equality.h"

const TString ext = ".png";

int comp_of_methods(Int_t first, Int_t last, Int_t loopcount, const Int_t M, Int_t file_num = 0, Double_t cut_prob = 0.0, Double_t time_cut = 100.0)
{
  TChain *chain = new TChain("INTERF/h1");
  chain_init(chain, first, last);

  TString file_name[3];

  file_name[0] = "neuvtx_tri_kin_fit_1_56_100_5_no_time.root";
  file_name[1] = "neuvtx_tri_kin_fit_" + std::to_string(first) + "_" + std::to_string(last) + "_" + std::to_string(loopcount) + "_" + std::to_string(M) + ".root"; //"neuvtx_tri_kin_fit_1_56_30_5_time.root";
  file_name[2] = "neuvtx_tri_rec_1_56.root";

  TFile *file = new TFile(file_name[file_num]);
  TTree *tree;

  if (file_num == 0 || file_num == 1)
    tree = (TTree *)file->Get("h_tri_kin_fit");
  else
    tree = (TTree *)file->Get("h_tri");

  TFile *file_gen = new TFile("../Generated_vars/gen_vars_1_56.root");
  TTree *tree_gen = (TTree *)file_gen->Get("h_gen_vars");

  Float_t Kchboost[9], Knereclor[9], Knerec[9],
      Kchmc[9], Knemc[9], ip[3], ipmc[3], phi_mom[4], Dtmc, Tcl[50],
      cluster[5][200], bhabha_vtx[3], T0step1;
  UChar_t mctruth, mcisr, g4taken[4];
  UChar_t pidmc[200], vtxmc[200], mother[200];
  Int_t ntmc, nvtxmc, nclu;

  chain->SetBranchAddress("Kchboost", Kchboost);
  chain->SetBranchAddress("Knereclor", Knereclor);
  chain->SetBranchAddress("Knerec", Knerec);
  chain->SetBranchAddress("Kchmc", Kchmc);
  chain->SetBranchAddress("Knemc", Knemc);

  chain->SetBranchAddress("nclu", &nclu);
  chain->SetBranchAddress("Xcl", cluster[0]);
  chain->SetBranchAddress("Ycl", cluster[1]);
  chain->SetBranchAddress("Zcl", cluster[2]);
  chain->SetBranchAddress("Tcl", cluster[3]);
  chain->SetBranchAddress("Enecl", cluster[4]);

  chain->SetBranchAddress("Bx", &bhabha_vtx[0]);
  chain->SetBranchAddress("By", &bhabha_vtx[1]);
  chain->SetBranchAddress("Bz", &bhabha_vtx[2]);

  chain->SetBranchAddress("ip", ip);
  chain->SetBranchAddress("ipmc", ipmc);

  chain->SetBranchAddress("ntmc", &ntmc);
  chain->SetBranchAddress("nvtxmc", &nvtxmc);
  chain->SetBranchAddress("pidmc", pidmc);
  chain->SetBranchAddress("vtxmc", vtxmc);
  chain->SetBranchAddress("mother", mother);

  chain->SetBranchAddress("mcisr", &mcisr);
  chain->SetBranchAddress("T0step1", &T0step1);

  chain->SetBranchAddress("Bpx", &phi_mom[0]);
  chain->SetBranchAddress("Bpy", &phi_mom[1]);
  chain->SetBranchAddress("Bpz", &phi_mom[2]);
  chain->SetBranchAddress("Broots", &phi_mom[3]);

  chain->SetBranchAddress("mctruth", &mctruth);

  chain->SetBranchAddress("Dtmc", &Dtmc);

  chain->SetBranchAddress("g4taken", g4taken);

  Float_t gamma_kinfit[4][8], ip_kinfit[3], Knetri_kinfit[10], chi2min, bhabha_mom_err[4], sol1[4], sol2[4], sol1err, sol2err;
  Int_t done_kinfit, g4taken_kinfit[4], bunchnum, chosen, region[4], clusindgood[4];
  TMatrixD *cov_matrix = new TMatrixD(27, 27);
  TVectorD *min_const = new TVectorD(M);
  TVectorD *lag_mult = new TVectorD(M), *X_min = new TVectorD(27), *X_init = new TVectorD(27);

  if (file_num == 0 || file_num == 1)
  {
    tree->SetBranchAddress("fourgamma1tri_kinfit", gamma_kinfit[0]);
    tree->SetBranchAddress("fourgamma2tri_kinfit", gamma_kinfit[1]);
    tree->SetBranchAddress("fourgamma3tri_kinfit", gamma_kinfit[2]);
    tree->SetBranchAddress("fourgamma4tri_kinfit", gamma_kinfit[3]);

    tree->SetBranchAddress("fourKnetri_kinfit", Knetri_kinfit);

    tree->SetBranchAddress("iptri_kinfit", ip_kinfit);
    tree->SetBranchAddress("done4_kinfit", &done_kinfit);

    tree->SetBranchAddress("g4takentri_kinfit", g4taken_kinfit);

    if (file_num == 1)
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

    tree->SetBranchAddress("sol1", sol1);
    tree->SetBranchAddress("sol2", sol2);
    tree->SetBranchAddress("sol1err", &sol1err);
    tree->SetBranchAddress("sol2err", &sol2err);

    tree->SetBranchAddress("chosen", &chosen);

    tree->SetBranchAddress("fourg4taken", g4taken_kinfit);

    tree->SetBranchAddress("iptri", ip_kinfit);
    tree->SetBranchAddress("done4", &done_kinfit);
  }

  tree_gen->SetBranchAddress("region", region);
  tree_gen->SetBranchAddress("clusindgood", clusindgood);

  chain->AddFriend(tree);
  chain->AddFriend(tree_gen);

  UInt_t nentries = (UInt_t)chain->GetEntries();

  Float_t lengthneu_mc, lengthneu_tri, lengthneu_tri_CM, lengthneu_rec, lengthch_mc, lengthch_rec,
      v_Kneutri, v_Kneutri_CM, v_Kchrec, v_Kneurec, v_Kneumc, v_Kchmc,
      t_chmc, t_neumc, t_chrec, t_neutri, t_neurec, Rtneu_rec, Rtneu_mc, Rtneu_tri,
      angle_path_mom, cos_path_mom, sigmas_init[24];

  TString id_hist, id_canva;

  TCanvas *canvas[2][100], *canvas_pulls[2][24];
  Int_t width = 750, height = 750;

  TH1 *chi2_hist[2], *prob_hist[2], *clusenergy_hist[2];
  TH1 *beta_hist[2];
  TH2 *neu_vtx_corr[2][4], *neu_mom[2][4], *ip_coor[2][3];
  TH2 *beta_time[2][2], *cluscorr_hist[2][2];

  TH2 *sigmas_std[2][5], *sigmas_tri[2][5], *sigmas_tri_kin_fit[2][5], *angle_vs_time[2][2], *dist_vs_time[2][2], *sol_err_hist[2];
  TH1 *neu_std_hist[2][5], *neu_tri_hist[2][5], *neu_tri_kin_fit_hist[2][5], *chosen_hist[2];
  TH1 *res_std_hist[2][5], *res_tri_hist[2][5], *res_tri_kin_fit_hist[2][5];
  TH1 *pulls[2][8], *first_clus_hist[2];
  TH1 *bunch[2], *pidmc_hist[2], *mcisr_hist[2], *tcl_hist[2];

  TString name[5] = {"x", "y", "z", "time", "length"};

  const UInt_t number_of_points = 10;

  for (Int_t j = 0; j < 2; j++)
  {
    for (Int_t i = 0; i < 100; i++)
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
    neu_vtx_corr[j][3] = new TH2F(id_hist, "", 200, 0.0, 20.0, 200, -10.0, 30.0);

    id_hist = "Mom" + std::to_string(j) + std::to_string(3);
    neu_mom[j][3] = new TH2F(id_hist, "", 100, 490, 550, 100, 490, 550);

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
        pulls[j][i] = new TH1F(id_hist, ";#vec{X}^{tri}_{neu,1} - #vec{X}^{gen}_{neu,1} [cm];Counts", 201, -300, 300);
      else if (i == 5)
        pulls[j][i] = new TH1F(id_hist, ";#vec{X}^{tri}_{neu,2} - #vec{X}^{gen}_{neu,2} [cm];Counts", 201, -300, 300);
      else if (i == 6)
        pulls[j][i] = new TH1F(id_hist, ";#vec{X}^{tri}_{neu,3} - #vec{X}^{gen}_{neu,3} [cm];Counts", 201, -300, 300);
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

    id_hist = "Chosen" + std::to_string(j);
    chosen_hist[j] = new TH1I(id_hist, ";Chosen solution;Counts", 101, -1, 3);

    id_hist = "Pidmc" + std::to_string(j);
    pidmc_hist[j] = new TH1I(id_hist, ";PID;Counts", 101, 0, 20);

    id_hist = "MCISR" + std::to_string(j);
    mcisr_hist[j] = new TH1I(id_hist, ";ISR MC Flag;Counts", 101, -1, 2);

    id_hist = "First_clus" + std::to_string(j);
    first_clus_hist[j] = new TH1F(id_hist, ";T_{cl} - #frac{d_{cl}}{c} [ns];Counts", 101, -6.0, 6.0);

    id_hist = "Tcl_diff" + std::to_string(j);
    tcl_hist[j] = new TH1F(id_hist, ";T_{cl} - #frac{d_{cl}}{c} - t_{neu} [ns];Counts", 201, -100.0, 50.0);

    id_hist = "Energyofcl" + std::to_string(j);
    clusenergy_hist[j] = new TH1F(id_hist, ";Inv mass [MeV/c^{2}];Counts", 201, 0.0, 1000.0);

    id_hist = "Angle 0 vs time" + std::to_string(j);
    angle_vs_time[j][0] = new TH2F(id_hist, ";min(#angle(#vec{p}_{#gamma,i},#vec{p}_{#gamma,j})) [#circ];t^{tri}_{neu} [ns]", 100, 0, 180, 100, -10, 25);

    id_hist = "Angle 180 vs time" + std::to_string(j);
    angle_vs_time[j][1] = new TH2F(id_hist, ";min(180#circ - #angle(#vec{p}_{#gamma,i},#vec{p}_{#gamma,j})) [#circ];t^{tri}_{neu} [ns]", 100, 0, 180, 100, -10, 25);

    id_hist = "Min dist vs angle 0" + std::to_string(j);
    dist_vs_time[j][0] = new TH2F(id_hist, ";min(dist(#vec{X}_{cl,i}),dist(#vec{X}_{cl,j})) [cm];min(#angle(#vec{p}_{#gamma,i},#vec{p}_{#gamma,j})) [#circ]", 100, 0, 500, 100, 0, 180);

    id_hist = "Min dist vs angle 180" + std::to_string(j);
    dist_vs_time[j][1] = new TH2F(id_hist, ";min(dist(#vec{X}_{cl,i}),dist(#vec{X}_{cl,j})) [cm];min(180#circ - #angle(#vec{p}_{#gamma,i},#vec{p}_{#gamma,j})) [#circ]", 100, 0, 500, 100, 0, 180);

    id_hist = "Sol error" + std::to_string(j);
    sol_err_hist[j] = new TH2F(id_hist, ";c#beta_{K^{0},1} t_{neu,1} - d_{neu,1} [cm];c#beta_{K^{0},2} t_{neu,2} - d_{neu,2} [cm]", 100, -200, 100, 100, -50, 500);

    id_hist = "clus_corr_x_y" + std::to_string(j);
    cluscorr_hist[j][0] = new TH2F(id_hist, ";x-axis [cm];y-axis [cm]", 100, -250, 250, 100, -250, 250);

    id_hist = "clus_corr_x_z" + std::to_string(j);
    cluscorr_hist[j][1] = new TH2F(id_hist, ";x-axis [cm];z-axis [cm]", 100, -250, 250, 100, -200, 200);

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

  Int_t count_angle = 0, min_ind[6], min_ind_180[6], min_ind_dist[6], isEqual = 0, bad_clus = 0;
  Float_t vec_before[4], vec_after[4], boost[3], ip_before[4], ip_after[4], d_cl, angle_min_0, angle_min_180, angle[6], angle_180[6], cos_tmp, num, den, dist[6], dist_min, tcl;

  Bool_t cut = true;

  Int_t counts_cut = 0, counts_all = 0, counts_bad = 0, counts_good = 0, g4taken_corr[4], counts_bad_smaller = 0, counts_bad_bigger = 0, counts_bad1 = 0, counts_bad2 = 0, counts_bad3 = 0, counts_bad4 = 0;

  for (Int_t i = 0; i < nentries; i++)
  {
    chain->GetEntry(i);

    count_angle = 0;

    if ((mctruth == 1 || mctruth == 2))
    {

      // Kaon path length
      lengthch_mc = sqrt(pow(Kchmc[6] - ipmc[0], 2) + pow(Kchmc[7] - ipmc[1], 2) + pow(Kchmc[8] - ipmc[2], 2));
      lengthch_rec = sqrt(pow(Kchboost[6] - ip[0], 2) + pow(Kchboost[7] - ip[1], 2) + pow(Kchboost[8] - ip[2], 2));

      lengthneu_mc = sqrt(pow(Knemc[6] - ipmc[0], 2) + pow(Knemc[7] - ipmc[1], 2) + pow(Knemc[8] - ipmc[2], 2));
      lengthneu_tri = sqrt(pow(Knetri_kinfit[6] - ip_kinfit[0], 2) + pow(Knetri_kinfit[7] - ip_kinfit[1], 2) + pow(Knetri_kinfit[8] - ip_kinfit[2], 2));
      lengthneu_rec = sqrt(pow(Knerec[6] - ip[0], 2) + pow(Knerec[7] - ip[1], 2) + pow(Knerec[8] - ip[2], 2));
      //
      // Kaon transverse radius
      Rtneu_mc = sqrt(pow(Knemc[6] - ipmc[0], 2) + pow(Knemc[7] - ipmc[1], 2));
      Rtneu_tri = sqrt(pow(Knetri_kinfit[6] - ip_kinfit[0], 2) + pow(Knetri_kinfit[7] - ip_kinfit[1], 2));
      Rtneu_rec = sqrt(pow(Knerec[6] - ip[0], 2) + pow(Knerec[7] - ip[1], 2));
      //
      // Kaon velocity
      v_Kchmc = c_vel * Kchmc[4] / Kchmc[3];
      v_Kchrec = c_vel * Kchboost[4] / Kchboost[3];

      v_Kneumc = c_vel * Knemc[4] / Knemc[3];
      v_Kneutri = c_vel * Knetri_kinfit[4] / Knetri_kinfit[3];
      v_Kneurec = c_vel * Knerec[4] / Knerec[3];
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
      if (file_num == 0 || file_num == 1)
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
          cut = TMath::Prob(chi2min, M) > 0;
        else
          cut = TMath::Prob(chi2min, M) > cut_prob;

        if (cut && done_kinfit == 1)
        {

          for (Int_t k1 = 0; k1 < 3; k1++)
            for (Int_t k2 = k1 + 1; k2 < 4; k2++)
            {
              num = (gamma_kinfit[k1][0] * gamma_kinfit[k2][0] +
                     gamma_kinfit[k1][1] * gamma_kinfit[k2][1] +
                     gamma_kinfit[k1][2] * gamma_kinfit[k2][2]);

              den = sqrt(pow(gamma_kinfit[k1][0], 2) +
                         pow(gamma_kinfit[k1][1], 2) +
                         pow(gamma_kinfit[k1][2], 2)) *
                    sqrt(pow(gamma_kinfit[k2][0], 2) +
                         pow(gamma_kinfit[k2][1], 2) +
                         pow(gamma_kinfit[k2][2], 2));

              cos_tmp = num / den;

              angle[count_angle] = 180. * acos(cos_tmp) / M_PI;
              angle_180[count_angle] = 180 - angle[count_angle];

              dist[count_angle] = sqrt(pow(gamma_kinfit[k1][4] - gamma_kinfit[k2][4], 2) +
                                       pow(gamma_kinfit[k1][5] - gamma_kinfit[k2][5], 2) +
                                       pow(gamma_kinfit[k1][6] - gamma_kinfit[k2][6], 2));
              count_angle++;
            }

          TMath::Sort(count_angle, angle, min_ind, kFALSE);
          TMath::Sort(count_angle, angle_180, min_ind_180, kFALSE);
          TMath::Sort(count_angle, dist, min_ind_dist, kFALSE);

          count_angle = 0;

          angle_min_0 = angle[min_ind[0]];
          angle_min_180 = angle_180[min_ind_180[0]];
          dist_min = dist[min_ind_dist[0]];

          if (j == 0)
          {
            counts_all++;

            if (Knetri_kinfit[9] > time_cut)
              counts_cut++;

            //! Check how many cluster sets are found good and which bad

            for (Int_t k = 0; k < 4; k++)
              g4taken_corr[k] = g4taken[k] - 1;

            isEqual = arr_eq(clusindgood, g4taken_kinfit);

            bad_clus = isEqual;

            if (isEqual == 0)
              counts_good++;
            else
            {
              counts_bad++;

              if (gamma_kinfit[0][7] > Knetri_kinfit[9] && gamma_kinfit[1][7] > Knetri_kinfit[9] &&
                  gamma_kinfit[2][7] > Knetri_kinfit[9] && gamma_kinfit[3][7] > Knetri_kinfit[9])
                counts_bad_bigger++;
              else
                counts_bad_smaller++;

              if (bad_clus == 1)
                counts_bad1++;
              if (bad_clus == 2)
                counts_bad2++;
              if (bad_clus == 3)
                counts_bad3++;
              if (bad_clus == 4)
                counts_bad4++;
            }
          }

          if (isEqual >= 0)
          {
            chi2_hist[j]->Fill(chi2min);
            prob_hist[j]->Fill(TMath::Prob(chi2min, M));

            beta_hist[j]->Fill(v_Kneutri / c_vel);

            bunch[j]->Fill(bunchnum);
          }

          for (Int_t k = 0; k < 3; k++)
            ip_coor[j][k]->Fill(ipmc[k], ip_kinfit[k]);

          for (Int_t k = 0; k < 5; k++)
          {
            if (k < 3)
            {
              if (isEqual >= 0)
              {
                neu_vtx_corr[j][k]->Fill(Knemc[6 + k], Knetri_kinfit[6 + k]);
                sigmas_std[j][k]->Fill(abs(Knemc[6 + k] - ipmc[k]), Knetri_kinfit[6 + k] - Knemc[6 + k]);
                pulls[j][4 + k]->Fill(Knetri_kinfit[6 + k] - Knemc[6 + k]);
                neu_mom[j][k]->Fill(Knemc[k], Knetri_kinfit[k]);
              }
            }
            else if (k == 3)
            {
              if (isEqual >= 0)
              {
                neu_vtx_corr[j][3]->Fill(t_neumc, Knetri_kinfit[9]); // Knetri_kinfit[6 + k]);
                sigmas_std[j][3]->Fill(lengthneu_mc, (Knetri_kinfit[9] - t_neurec) / tau_S_nonCPT);
                pulls[j][4 + k]->Fill((Knetri_kinfit[9] - t_neumc) / tau_S_nonCPT);
                neu_mom[j][k]->Fill(Knemc[k], Knetri_kinfit[k]);
              }
            }
            else
            {
              if (isEqual >= 0)
                sigmas_std[j][4]->Fill(lengthneu_mc, (lengthneu_tri - lengthneu_mc));
            }

            if (isEqual >= 0)
            {
              pulls[j][k]->Fill(Knetri_kinfit[k] - Knemc[k]);
              pulls[j][4 + k]->Fill(Knetri_kinfit[4 + k] - Knemc[4 + k]);
            }
          }

          d_cl = sqrt(pow(cluster[0][0] - bhabha_vtx[0], 2) + pow(cluster[1][0] - bhabha_vtx[1], 2) + pow(cluster[2][0] - bhabha_vtx[2], 2));

          if (isEqual >= 0)
          {

            for (Int_t k = 0; k < ntmc; k++)
              pidmc_hist[j]->Fill(pidmc[k]);

            angle_vs_time[j][0]->Fill(angle_min_0, Knetri_kinfit[9]);

            angle_vs_time[j][1]->Fill(angle_min_180, Knetri_kinfit[9]);
            dist_vs_time[j][0]->Fill(dist_min, angle_min_0);
            dist_vs_time[j][1]->Fill(dist_min, angle_min_180);

            sol_err_hist[j]->Fill(sol1err, sol2err);
            chosen_hist[j]->Fill(chosen);

            for (Int_t k = 0; k < 4; k++)
            {
              tcl = gamma_kinfit[k][7] - (sqrt(pow(gamma_kinfit[k][4] - Knetri_kinfit[6], 2) + pow(gamma_kinfit[k][5] - Knetri_kinfit[7], 2) + pow(gamma_kinfit[k][6] - Knetri_kinfit[8], 2)) / c_vel) - Knetri_kinfit[9];
              tcl_hist[j]->Fill(tcl);

              cluscorr_hist[j][0]->Fill(gamma_kinfit[0][4], gamma_kinfit[1][4]);
              cluscorr_hist[j][1]->Fill(gamma_kinfit[0][5], gamma_kinfit[1][5]);
            }

            clusenergy_hist[j]->Fill(Knetri_kinfit[5]);

            first_clus_hist[j]->Fill(cluster[3][0] - (d_cl / c_vel));

            mcisr_hist[j]->Fill(mcisr);
          }
        }
      }
    }
  }

  std::ofstream log_file;
  log_file.open("tri.log");

  log_file << "Events all: " << counts_all << std::endl;
  log_file << "Events after tneu < 9 ns: " << counts_cut << std::endl;
  log_file << "Efficiency after tneu < 9 ns: " << counts_cut / (Float_t)counts_all << std::endl;
  log_file << "Good cluster sets eff: " << counts_good / (Float_t)counts_all << std::endl;
  log_file << "Bad cluster sets: " << counts_bad << std::endl;
  log_file << "Bad clusters 1: " << counts_bad1 / (Float_t)counts_bad << std::endl;
  log_file << "Bad clusters 2: " << counts_bad2 / (Float_t)counts_bad << std::endl;
  log_file << "Bad clusters 3: " << counts_bad3 / (Float_t)counts_bad << std::endl;
  log_file << "Bad clusters 4: " << counts_bad4 / (Float_t)counts_bad << std::endl;
  log_file << "Bad cluster smaller: " << counts_bad_smaller << std::endl;
  log_file << "Bad cluster bigger: " << counts_bad_bigger << std::endl;
  log_file << "Bad cluster sets eff: " << counts_bad / (Float_t)counts_all << std::endl;
  log_file << "Bad cluster smaller eff: " << counts_bad_smaller / (Float_t)counts_bad << std::endl;
  log_file << "Bad cluster bigger eff: " << counts_bad_bigger / (Float_t)counts_bad << std::endl;

  log_file.close();

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

  //! sol1 and sol2 errors
  canvas[0][74]->SetRightMargin(0.15);
  canvas[0][74]->SetLeftMargin(0.17);
  canvas[0][74]->cd();

  chi2_hist[0]->GetYaxis()->SetMaxDigits(3);

  id_canva = "chosen_solutions";

  chosen_hist[0]->GetXaxis()->CenterTitle();
  chosen_hist[0]->GetYaxis()->CenterTitle();
  chosen_hist[0]->GetYaxis()->SetRangeUser(0.0, 1.3 * chosen_hist[0]->GetMaximum());
  chosen_hist[0]->SetLineColor(kBlack);
  chosen_hist[1]->SetLineColor(kRed);
  chosen_hist[0]->Draw();
  chosen_hist[1]->Draw("SAME");
  canvas[0][74]->Print(id_canva + ext);

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

  //!
  canvas[0][26]->SetRightMargin(0.15);
  canvas[0][26]->SetLeftMargin(0.17);
  canvas[0][26]->cd();
  canvas[0][26]->SetLogy();

  pidmc_hist[0]->GetYaxis()->SetMaxDigits(3);

  id_canva = "pidmc";

  pidmc_hist[0]->GetXaxis()->CenterTitle();
  pidmc_hist[0]->GetYaxis()->CenterTitle();
  // pidmc_hist[0]->GetXaxis()->SetTitle(x_title);
  pidmc_hist[0]->GetYaxis()->SetTitle("Counts");
  if (canvas[0][26]->GetLogy())
    pidmc_hist[0]->GetYaxis()->SetRangeUser(1.0, 10. * pidmc_hist[0]->GetMaximum());
  else
    pidmc_hist[0]->GetYaxis()->SetRangeUser(0.0, 1.3 * pidmc_hist[0]->GetMaximum());
  pidmc_hist[0]->SetLineColor(kBlack);
  pidmc_hist[1]->SetLineColor(kRed);
  pidmc_hist[0]->Draw();
  pidmc_hist[1]->Draw("SAME");
  canvas[0][26]->Print(id_canva + ext);
  //!
  canvas[0][27]->SetRightMargin(0.15);
  canvas[0][27]->SetLeftMargin(0.17);
  canvas[0][27]->cd();
  canvas[0][27]->SetLogy();

  mcisr_hist[0]->GetYaxis()->SetMaxDigits(3);

  id_canva = "mcisr";

  mcisr_hist[0]->GetXaxis()->CenterTitle();
  mcisr_hist[0]->GetYaxis()->CenterTitle();
  // mcisr_hist[0]->GetXaxis()->SetTitle(x_title);
  mcisr_hist[0]->GetYaxis()->SetTitle("Counts");
  if (canvas[0][27]->GetLogy())
    mcisr_hist[0]->GetYaxis()->SetRangeUser(1.0, 10. * mcisr_hist[0]->GetMaximum());
  else
    mcisr_hist[0]->GetYaxis()->SetRangeUser(0.0, 1.3 * mcisr_hist[0]->GetMaximum());
  mcisr_hist[0]->SetLineColor(kBlack);
  mcisr_hist[1]->SetLineColor(kRed);
  mcisr_hist[0]->Draw();
  mcisr_hist[1]->Draw("SAME");
  canvas[0][27]->Print(id_canva + ext);
  //!
  canvas[0][28]->SetRightMargin(0.15);
  canvas[0][28]->SetLeftMargin(0.17);
  canvas[0][28]->cd();
  canvas[0][28]->SetLogy();

  first_clus_hist[0]->GetYaxis()->SetMaxDigits(3);

  id_canva = "first_clus";

  first_clus_hist[0]->GetXaxis()->CenterTitle();
  first_clus_hist[0]->GetYaxis()->CenterTitle();
  // first_clus_hist[0]->GetXaxis()->SetTitle("T0step1 [ns]");
  first_clus_hist[0]->GetYaxis()->SetTitle("Counts");
  if (canvas[0][28]->GetLogy())
    first_clus_hist[0]->GetYaxis()->SetRangeUser(1.0, 10. * first_clus_hist[0]->GetMaximum());
  else
    first_clus_hist[0]->GetYaxis()->SetRangeUser(0.0, 1.3 * first_clus_hist[0]->GetMaximum());
  first_clus_hist[0]->SetLineColor(kBlack);
  first_clus_hist[1]->SetLineColor(kRed);
  first_clus_hist[0]->Draw();
  first_clus_hist[1]->Draw("SAME");
  canvas[0][28]->Print(id_canva + ext);
  //!
  //!
  canvas[0][29]->SetRightMargin(0.15);
  canvas[0][29]->SetLeftMargin(0.17);
  canvas[0][29]->cd();
  canvas[0][29]->SetLogy();

  tcl_hist[0]->GetYaxis()->SetMaxDigits(3);

  id_canva = "tcl";

  tcl_hist[0]->GetXaxis()->CenterTitle();
  tcl_hist[0]->GetYaxis()->CenterTitle();
  // tcl_hist[0]->GetXaxis()->SetTitle("T0step1 [ns]");
  tcl_hist[0]->GetYaxis()->SetTitle("Counts");
  if (canvas[0][29]->GetLogy())
    tcl_hist[0]->GetYaxis()->SetRangeUser(1.0, 10. * tcl_hist[0]->GetMaximum());
  else
    tcl_hist[0]->GetYaxis()->SetRangeUser(0.0, 1.3 * tcl_hist[0]->GetMaximum());
  tcl_hist[0]->SetLineColor(kBlack);
  tcl_hist[1]->SetLineColor(kRed);
  tcl_hist[0]->Draw();
  tcl_hist[1]->Draw("SAME");
  canvas[0][29]->Print(id_canva + ext);
  //!

  for (Int_t j = 0; j < 2; j++)
  {
    gStyle->SetOptStat("iMROU");
    gStyle->SetStatX(0.85);
    gStyle->SetStatY(0.9);

    id_canva = "Angle0vstime" + std::to_string(j + 1);

    canvas[j][70]->SetRightMargin(0.15);
    canvas[j][70]->SetLeftMargin(0.17);
    canvas[j][70]->cd();

    angle_vs_time[j][0]->Draw("COLZ");

    canvas[j][70]->Print(id_canva + ext);

    id_canva = "Angle180vstime" + std::to_string(j + 1);

    canvas[j][71]->SetRightMargin(0.15);
    canvas[j][71]->SetLeftMargin(0.17);
    canvas[j][71]->cd();

    angle_vs_time[j][1]->Draw("COLZ");

    canvas[j][71]->Print(id_canva + ext);

    id_canva = "Distvsangle0" + std::to_string(j + 1);

    canvas[j][72]->SetRightMargin(0.15);
    canvas[j][72]->SetLeftMargin(0.17);
    canvas[j][72]->cd();

    dist_vs_time[j][0]->Draw("COLZ");

    canvas[j][72]->Print(id_canva + ext);

    id_canva = "Distvsangle180" + std::to_string(j + 1);

    canvas[j][74]->SetRightMargin(0.15);
    canvas[j][74]->SetLeftMargin(0.17);
    canvas[j][74]->cd();

    dist_vs_time[j][1]->Draw("COLZ");

    canvas[j][74]->Print(id_canva + ext);

    id_canva = "Solerr" + std::to_string(j + 1);

    canvas[j][73]->SetRightMargin(0.15);
    canvas[j][73]->SetLeftMargin(0.17);
    canvas[j][73]->cd();

    sol_err_hist[j]->Draw("COLZ");

    canvas[j][73]->Print(id_canva + ext);

    id_canva = "Clusenergy" + std::to_string(j + 1);

    canvas[j][75]->SetRightMargin(0.15);
    canvas[j][75]->SetLeftMargin(0.17);
    canvas[j][75]->cd();

    clusenergy_hist[j]->Draw();

    canvas[j][75]->Print(id_canva + ext);

    id_canva = "Cluscorr_x_y_" + std::to_string(j + 1);

    canvas[j][76]->SetRightMargin(0.15);
    canvas[j][76]->SetLeftMargin(0.17);
    canvas[j][76]->cd();

    cluscorr_hist[j][0]->Draw("COLZ");

    canvas[j][76]->Print(id_canva + ext);

    id_canva = "Cluscorr_x_z_" + std::to_string(j + 1);

    canvas[j][77]->SetRightMargin(0.15);
    canvas[j][77]->SetLeftMargin(0.17);
    canvas[j][77]->cd();

    cluscorr_hist[j][1]->Draw("COLZ");

    canvas[j][77]->Print(id_canva + ext);

    for (Int_t i = 0; i < 3; i++)
    {
      canvas[j][i]->SetRightMargin(0.15);
      canvas[j][i]->SetLeftMargin(0.17);
      canvas[j][i]->cd();

      id_canva = "coor" + std::to_string(j + 1) + std::to_string(i + 1);
      x_title = Form("#vec{X}_{neu,%d}^{gen}[cm]", i + 1);
      y_title = Form("#vec{X}_{neu,%d}^{tri}[cm]", i + 1);

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
    x_title = Form("t_{neu}^{gen} - t_{neu}^{gen} [ns]");
    y_title = Form("t_{neu}^{tri} - t_{neu}^{gen} [ns]");

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
    y_title = "#sigma(#vec{X}^{tri}_{K#rightarrow#pi^{0}#pi^{0},i} - #vec{X}^{gen}_{K#rightarrow#pi^{0}#pi^{0},i}) [cm]";

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
    y_title = "#sigma(t_{K#rightarrow#pi^{0}#pi^{0}}^{tri} - t_{K#rightarrow#pi^{0}#pi^{0}}^{gen}) [#tau_{S}]";
    canvas[j][23]->cd();
    res_std_hist[j][3]->SetXTitle(x_title);
    res_std_hist[j][3]->SetYTitle(y_title);
    res_std_hist[j][3]->GetYaxis()->SetRangeUser(0.0, 6);
    res_std_hist[j][3]->Draw("PE1");
    canvas[j][23]->Print(id_canva + ext);
    //!
    id_canva = "sigmas_std_length" + std::to_string(j + 1);
    x_title = "K#rightarrow#pi^{0}#pi^{0} path generated [cm]";
    y_title = "#sigma(l_{K#rightarrow#pi^{0}#pi^{0}}^{tri} - l_{K#rightarrow#pi^{0}#pi^{0}}^{gen}) [#tau_{S}]";
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

  return 0;
}
