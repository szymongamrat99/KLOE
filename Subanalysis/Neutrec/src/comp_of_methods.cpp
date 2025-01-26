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
#include <TLine.h>
#include <TBranch.h>
#include <TChain.h>
#include <TLegend.h>

#include "kloe_class.h"
#include "interf_function.h"
#include "chi2_dist.h"
#include "uncertainties.h"
#include "lorentz_transf.h"
#include "triple_gaus.h"

#include "../inc/trilateration.hpp"

int comp_of_methods(int first_file, int last_file, Controls::DataType data_type)
{
  ErrorHandling::ErrorLogs errLogger;

  Double_t cut_prob = 0.0, time_cut = 100.0;

  int
      loopcount = (Int_t)properties["variables"]["KinFit"]["Trilateration"]["loopCount"],
      M = (Int_t)properties["variables"]["KinFit"]["Trilateration"]["numOfConstraints"],
      jmin = (Int_t)properties["variables"]["KinFit"]["Trilateration"]["bunchMin"],
      jmax = (Int_t)properties["variables"]["KinFit"]["Trilateration"]["bunchMax"],
      range = Int_t(jmax - jmin) + 1,
      bunch_int = 9;

  TChain *chain = new TChain("INTERF/h1");
  chain_init(chain, first_file, last_file);

  TString
      filename_gen = gen_vars_dir + root_files_dir + gen_vars_filename + first_file + "_" + last_file + ext_root,
      filename_trilateration = neutrec_dir + root_files_dir + gen_vars_filename + first_file + "_" + last_file + ext_root,
      filename_trilateration_kin_fit = neutrec_dir + root_files_dir + neu_trilateration_kin_fit_filename + first_file + "_" + last_file + "_" + loopcount + "_" + M + "_" + range + "_" + int(data_type) + ext_root,
      filename_triangle = neutrec_dir + root_files_dir + neu_triangle_filename + first_file + "_" + last_file + "_" + loopcount + "_" + M + "_" + range + "_" + int(data_type) + ext_root;

  std::vector<TString> file_name, tree_name;

  file_name.push_back(filename_trilateration);
  file_name.push_back(filename_trilateration_kin_fit);
  file_name.push_back(filename_triangle);

  tree_name.push_back(neutrec_tri_tree);
  tree_name.push_back(neutrec_kin_fit_tree);
  tree_name.push_back(neutrec_triangle_tree);

  Controls::Menu menu(3, file_name);

  menu.ShowOpt();
  menu.EndMenu();

  Int_t
      file_num;
  TFile
      *file,
      *file_gen;
  TTree
      *tree,
      *tree_gen;

  try
  {
    std::cin >> file_num;

    if (!std::cin)
      throw(ErrorHandling::ErrorCodes::DATA_TYPE);
    else if (file_num < 1 || file_num > file_name.size())
      throw(ErrorHandling::ErrorCodes::MENU_RANGE);

    file = new TFile(file_name[file_num]);
    file_gen = new TFile(filename_gen);

    if (file->IsZombie() || file_gen->IsZombie())
      throw(ErrorHandling::ErrorCodes::FILE_NOT_EXIST);

    tree = (TTree *)file->Get(tree_name[file_num]);
    tree_gen = (TTree *)file_gen->Get(gen_vars_tree);

    if (tree->IsZombie() || tree_gen->IsZombie())
      throw(ErrorHandling::ErrorCodes::TREE_NOT_EXIST);
  }
  catch (ErrorHandling::ErrorCodes err)
  {
    errLogger.setErrCount(err);
    errLogger.getErrLog(err);
    return int(err);
  }

  std::string interface;

  if (file_num == 1)
    interface = "trilateration";
  if (file_num == 2)
    interface = "pureTriangle";

  BaseKinematics baseKin;

  UChar_t errId, cutId;

  chain->SetBranchAddress("Kchboost", baseKin.Kchboost);
  chain->SetBranchAddress("Knereclor", baseKin.Knereclor);
  chain->SetBranchAddress("Knerec", baseKin.Knerec);
  chain->SetBranchAddress("Kchrec", baseKin.Kchrec);
  chain->SetBranchAddress("Kchmc", baseKin.Kchmc);
  chain->SetBranchAddress("Knemc", baseKin.Knemc);

  chain->SetBranchAddress("trk1", baseKin.trk[0]);
  chain->SetBranchAddress("trk2", baseKin.trk[1]);

  chain->SetBranchAddress("nclu", &baseKin.nclu);
  chain->SetBranchAddress("Xcl", baseKin.cluster[0]);
  chain->SetBranchAddress("Ycl", baseKin.cluster[1]);
  chain->SetBranchAddress("Zcl", baseKin.cluster[2]);
  chain->SetBranchAddress("Tcl", baseKin.cluster[3]);
  chain->SetBranchAddress("Enecl", baseKin.cluster[4]);

  chain->SetBranchAddress("Bx", &baseKin.bhabha_vtx[0]);
  chain->SetBranchAddress("By", &baseKin.bhabha_vtx[1]);
  chain->SetBranchAddress("Bz", &baseKin.bhabha_vtx[2]);

  chain->SetBranchAddress("ip", baseKin.ip);
  chain->SetBranchAddress("ipmc", baseKin.ipmc);

  chain->SetBranchAddress("ntmc", &baseKin.ntmc);
  chain->SetBranchAddress("nvtxmc", &baseKin.nvtxmc);
  chain->SetBranchAddress("nv", &baseKin.nv);
  chain->SetBranchAddress("pidmc", baseKin.pidmc);
  chain->SetBranchAddress("vtxmc", baseKin.vtxmc);
  chain->SetBranchAddress("mother", baseKin.mother);

  chain->SetBranchAddress("mcisr", &baseKin.mcisr);
  chain->SetBranchAddress("T0step1", &baseKin.T0step1);

  chain->SetBranchAddress("Bpx", &baseKin.phi_mom[0]);
  chain->SetBranchAddress("Bpy", &baseKin.phi_mom[1]);
  chain->SetBranchAddress("Bpz", &baseKin.phi_mom[2]);
  chain->SetBranchAddress("Broots", &baseKin.phi_mom[3]);

  chain->SetBranchAddress("mctruth", &baseKin.mctruth);
  chain->SetBranchAddress("mcflag", &baseKin.mcflag);

  chain->SetBranchAddress("Dtmc", &baseKin.Dtmc);

  chain->SetBranchAddress("g4taken", baseKin.g4taken);
  chain->SetBranchAddress("ncll", baseKin.ncll);

  chain->SetBranchAddress("Errid", &errId);
  chain->SetBranchAddress("Cutid", &cutId);

  Float_t gamma_kinfit[4][8], ip_kinfit[3], Knetri_kinfit[10], chi2min, bhabha_mom_err[4], sol1[4], sol2[4], sol1err, sol2err, trc[4], trcsum, minv4gam;
  Int_t done_kinfit = 10, g4taken_kinfit[4], bunchnum, chosen, region[4], clusindgood[4];
  TMatrixD *cov_matrix = new TMatrixD(27, 27);
  TVectorD *min_const = new TVectorD(M);
  TVectorD *lag_mult = new TVectorD(M), *X_min = new TVectorD(27), *X_init = new TVectorD(27);

  Float_t gamma_kinfit_no_bunch[4][8], ip_kinfit_no_bunch[3], Knetri_kinfit_no_bunch[10], chi2min_no_bunch, bhabha_mom_err_no_bunch[4], sol1_no_bunch[4], sol2_no_bunch[4], sol1err_no_bunch, sol2err_no_bunch;
  Int_t done_kinfit_no_bunch, g4taken_kinfit_no_bunch[4], bunchnum_no_bunch, chosen_no_bunch, region_no_bunch[4], clusindgood_no_bunch[4];

  if (file_num == 0)
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
  else if (file_num == 1)
  {
    tree->SetBranchAddress("fourgamma1tri_kinfit", gamma_kinfit[0]);
    tree->SetBranchAddress("fourgamma2tri_kinfit", gamma_kinfit[1]);
    tree->SetBranchAddress("fourgamma3tri_kinfit", gamma_kinfit[2]);
    tree->SetBranchAddress("fourgamma4tri_kinfit", gamma_kinfit[3]);

    tree->SetBranchAddress("fourKnetri_kinfit", Knetri_kinfit);

    tree->SetBranchAddress("iptri_kinfit", ip_kinfit);
    tree->SetBranchAddress("done4_kinfit", &done_kinfit);

    tree->SetBranchAddress("g4takentri_kinfit", g4taken_kinfit);

    tree->SetBranchAddress("bunchnum", &bunchnum);

    tree->SetBranchAddress("chi2min", &chi2min);
    tree->SetBranchAddress("min_cov", &cov_matrix);
    tree->SetBranchAddress("const_min", &min_const);
    tree->SetBranchAddress("lag_mult", &lag_mult);

    tree->SetBranchAddress("min_vars", &X_min);
    tree->SetBranchAddress("init_vars", &X_init);
  }
  else if (file_num == 2)
  {
    tree->SetBranchAddress("fourgamma1triangle", gamma_kinfit[0]);
    tree->SetBranchAddress("fourgamma2triangle", gamma_kinfit[1]);
    tree->SetBranchAddress("fourgamma3triangle", gamma_kinfit[2]);
    tree->SetBranchAddress("fourgamma4triangle", gamma_kinfit[3]);

    tree->SetBranchAddress("fourKnetriangle", Knetri_kinfit);

    tree->SetBranchAddress("g4taken_triangle", g4taken_kinfit);

    tree->SetBranchAddress("chi2min", &chi2min);

    tree->SetBranchAddress("trc", trc);
    tree->SetBranchAddress("trcsum", &trcsum);

    tree->SetBranchAddress("minv4gam", &minv4gam);

    tree->SetBranchAddress("iptriangle", ip_kinfit);
    tree->SetBranchAddress("done_triangle", &done_kinfit);
  }

  Float_t Knemc_bcg[9], Kchmc_bcg[9], ipmc_bcg[3];

  tree_gen->SetBranchAddress("region", region);
  tree_gen->SetBranchAddress("clusindgood", clusindgood);

  chain->AddFriend(tree);
  chain->AddFriend(tree_gen);

  UInt_t nentries = (UInt_t)chain->GetEntries();

  Float_t lengthneu_mc, lengthneu_tri, lengthneu_tri_CM, lengthneu_rec, lengthch_mc, lengthch_rec,
      v_Kneutri, v_Kneutri_CM, v_Kchrec, v_Kneurec, v_Kneumc, v_Kchmc,
      t_chmc, t_neumc, t_chrec, t_neutri, t_neurec, Rtneu_rec, Rtneu_mc, Rtneu_tri, Rtch_rec, Rtch_mc,
      angle_path_mom, cos_path_mom, sigmas_init[24];

  TString id_hist, id_canva;

  TCanvas *canvas[2][100], *canvas_pulls[2][24], *canvas_std[2][6];
  Int_t width = 750, height = 750;

  TH1 *chi2_hist[2], *chi2_no_bunch_hist[2], *chi2add[2], *chi2subst[2], *chi2div[2], *prob_hist[2], *clusenergy_hist[2];
  TH1 *betadt_hist[2], *betadtover90_hist[2], *betapm_hist[2], *betapmover90_hist[2];
  TH2 *neu_vtx_corr[2][4], *neu_mom[2][4], *ip_coor[2][3];
  TH2 *beta_time[2][2], *cluscorr_hist[2][2];

  TH2 *sigmas_std[2][5], *sigmas_tri[2][5], *sigmas_tri_kin_fit[2][5], *angle_vs_time[2][2], *dist_vs_time[2][2], *sol_err_hist[2], *beta1_hist[2], *beta2_hist[2];
  TH1 *neu_std_hist[2][5], *neu_tri_hist[2][5], *neu_tri_kin_fit_hist[2][5], *chosen_hist[2];
  TH1 *res_std_hist[2][5], *res_tri_hist[2][5], *res_tri_kin_fit_hist[2][5];
  TH1 *pulls[2][10], *first_clus_hist[2];
  TH1 *bunch[2], *pidmc_hist[2], *mcisr_hist[2], *tcl_hist[2];

  TH1 *chi2_bunch_hist[2][4];

  TString name[5] = {"x", "y", "z", "time", "length"};

  const UInt_t number_of_points = 10;

  TCanvas *canvas_proj[2][6][5];

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

      if (i < 6)
      {
        id_canva = "Canva_std" + std::to_string(j) + std::to_string(i);
        canvas_std[j][i] = new TCanvas(id_canva, "", width, height);

        for (Int_t k = 0; k < 5; k++)
        {
          id_canva = "Canva_proj" + std::to_string(j) + std::to_string(i) + std::to_string(k);
          canvas_proj[j][i][k] = new TCanvas(id_canva, "", width, height);
        }
      }
    }

    for (Int_t i = 0; i < 3; i++)
    {
      id_hist = "Coor" + std::to_string(j) + std::to_string(i);
      neu_vtx_corr[j][i] = new TH2F(id_hist, "", 100, -200, 200, 100, -200, 200);

      id_hist = "Mom" + std::to_string(j) + std::to_string(i);
      neu_mom[j][i] = new TH2F(id_hist, "", 200, -300, 300, 200, -300, 300);

      id_hist = "IP coor" + std::to_string(j) + std::to_string(i);
      ip_coor[j][i] = new TH2F(id_hist, "", 100, -10, 10, 100, -10, 10);
    }

    id_hist = "Beta1" + std::to_string(j);
    beta1_hist[j] = new TH2F(id_hist, "", 200, 0.15, 0.3, 200, 0.15, 0.3);

    id_hist = "Beta2" + std::to_string(j);
    beta2_hist[j] = new TH2F(id_hist, "", 200, 0.15, 0.3, 200, 0.15, 0.3);

    id_hist = "Coor" + std::to_string(j) + std::to_string(3);
    neu_vtx_corr[j][3] = new TH2F(id_hist, "", 200, 0.0, 50.0, 200, -5.0, 50.0);

    id_hist = "Mom" + std::to_string(j) + std::to_string(3);
    neu_mom[j][3] = new TH2F(id_hist, "", 100, 490, 550, 100, 490, 550);

    for (Int_t i = 0; i < 10; i++)
    {
      id_hist = "Var" + std::to_string(j) + std::to_string(i);
      if (i == 0)
        pulls[j][i] = new TH1F(id_hist, ";#vec{p}^{tri}_{neu,1} - #vec{p}^{gen}_{neu,1} [MeV/c];Counts", 101, -40, 40);
      else if (i == 1)
        pulls[j][i] = new TH1F(id_hist, ";#vec{p}^{tri}_{neu,2} - #vec{p}^{gen}_{neu,2} [MeV/c];Counts", 101, -40, 40);
      else if (i == 2)
        pulls[j][i] = new TH1F(id_hist, ";#vec{p}^{tri}_{neu,3} - #vec{p}^{gen}_{neu,3} [MeV/c];Counts", 101, -40, 40);
      else if (i == 3)
        pulls[j][i] = new TH1F(id_hist, ";E^{tri}_{neu} - E^{gen}_{neu,3} [MeV];Counts", 201, -6, 6);
      else if (i == 4)
        pulls[j][i] = new TH1F(id_hist, ";#vec{X}^{tri}_{neu,1} - #vec{X}^{gen}_{neu,1} [cm];Counts", 201, -8, 8);
      else if (i == 5)
        pulls[j][i] = new TH1F(id_hist, ";#vec{X}^{tri}_{neu,2} - #vec{X}^{gen}_{neu,2} [cm];Counts", 201, -8, 8);
      else if (i == 6)
        pulls[j][i] = new TH1F(id_hist, ";#vec{X}^{tri}_{neu,3} - #vec{X}^{gen}_{neu,3} [cm];Counts", 201, -20, 20);
      else if (i == 7)
        pulls[j][i] = new TH1F(id_hist, ";t^{tri}_{neu} - t^{gen}_{neu} [#tau_{S}];Counts", 101, -5, 10);
      else if (i == 8)
        pulls[j][i] = new TH1F(id_hist, ";#rho^{tri}_{neu} - #rho^{gen}_{neu} [cm];Counts", 101, -5, 15);
      else if (i == 9)
        pulls[j][i] = new TH1F(id_hist, ";d^{tri}_{neu} - d^{gen}_{neu} [cm];Counts", 101, -5, 15);
    }

    id_hist = "Chi2" + std::to_string(j);
    chi2_hist[j] = new TH1F(id_hist, "", 200, -1.0, 100.0);

    for (Int_t i = 0; i < 4; i++)
    {
      id_hist = "Chi2" + std::to_string(j) + std::to_string(i);
      chi2_bunch_hist[j][i] = new TH1F(id_hist, "", 200, -1.0, 100.0);
    }

    id_hist = "Chi2_no_bunch" + std::to_string(j);
    chi2_no_bunch_hist[j] = new TH1F(id_hist, "", 200, -1.0, 100.0);

    id_hist = "Chi2_add" + std::to_string(j);
    chi2add[j] = new TH1F(id_hist, "", 200, -1.0, 100.0);

    id_hist = "Chi2_subst" + std::to_string(j);
    chi2subst[j] = new TH1F(id_hist, ";#chi^{2}_{test} - #chi^{2}_{notest};Counts", 200, -1.0, 100.0);

    id_hist = "Chi2_div" + std::to_string(j);
    chi2div[j] = new TH1F(id_hist, ";#frac{#chi^{2}_{test} - #chi^{2}_{notest}}{#chi^{2}_{test} + #chi^{2}_{notest}};Counts", 200, -1.0, 100.0);

    id_hist = "Prob" + std::to_string(j);
    prob_hist[j] = new TH1F(id_hist, "", 100, 0.0, 1.0);

    id_hist = "Betadt" + std::to_string(j);
    betadt_hist[j] = new TH1F(id_hist, ";#beta_{K#rightarrow#pi^{0}#pi^{0}};Counts", 100, 0.15, 0.6);

    id_hist = "Betadtover90" + std::to_string(j);
    betadtover90_hist[j] = new TH1F(id_hist, ";#beta_{K#rightarrow#pi^{0}#pi^{0}};Counts", 100, 0.15, 0.6);

    id_hist = "Betapm" + std::to_string(j);
    betapm_hist[j] = new TH1F(id_hist, ";#beta_{K#rightarrow#pi^{0}#pi^{0}};Counts", 100, 0.15, 0.6);

    id_hist = "Betapmover90" + std::to_string(j);
    betapmover90_hist[j] = new TH1F(id_hist, ";#beta_{K#rightarrow#pi^{0}#pi^{0}};Counts", 100, 0.15, 0.6);

    id_hist = "Bunch_num" + std::to_string(j);
    bunch[j] = new TH1I(id_hist, ";Number of bunch correction;Counts", 11, -5, 5);

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
    cluscorr_hist[j][0] = new TH2F(id_hist, ";x-axis [cm];y-axis [cm]", 100, -50, 50, 100, -50, 50);

    id_hist = "clus_corr_x_z" + std::to_string(j);
    cluscorr_hist[j][1] = new TH2F(id_hist, ";x-axis [cm];z-axis [cm]", 100, -50, 50, 100, -50, 50);

    for (Int_t i = 0; i < 5; i++)
    {
      id_hist = Form(name[i] + "%d%d", j + 1, i + 1);
      if (i < 2)
      {
        neu_std_hist[j][i] = new TH1F(id_hist + " std coor", "", 100, -200, 200);
        neu_tri_hist[j][i] = new TH1F(id_hist + " tri coor", "", 100, -200, 200);
        neu_tri_kin_fit_hist[j][i] = new TH1F(id_hist + " kinfit coor", "", 100, -200, 200);
        res_std_hist[j][i] = new TH1F(id_hist + " std res", "", number_of_points, 0, 200);
        res_tri_hist[j][i] = new TH1F(id_hist + " tri res", "", number_of_points, 0, 200);
        res_tri_kin_fit_hist[j][i] = new TH1F(id_hist + " kinfit res", "", number_of_points, 0, 200);
        sigmas_std[j][i] = new TH2F(id_hist + " sigmas std", "", number_of_points, 0, 200, 200, -8, 8);
        sigmas_tri[j][i] = new TH2F(id_hist + " sigmas tri", "", number_of_points, 0, 200, 200, -8, 8);
        sigmas_tri_kin_fit[j][i] = new TH2F(id_hist + " sigmas kinfit", "", number_of_points, 0, 200, 200, -8, 8);
      }
      else if (i == 2)
      {
        neu_std_hist[j][i] = new TH1F(id_hist + " std", "", 100, 0, 5);
        neu_tri_hist[j][i] = new TH1F(id_hist + " tri", "", 100, 0, 5);
        neu_tri_kin_fit_hist[j][i] = new TH1F(id_hist + " kinfit", "", 100, 0, 5);
        res_std_hist[j][i] = new TH1F(id_hist + " std res", "", number_of_points, 0, 200);
        res_tri_hist[j][i] = new TH1F(id_hist + " tri res", "", number_of_points, 0, 200);
        res_tri_kin_fit_hist[j][i] = new TH1F(id_hist + " kinfit res", "", number_of_points, 0, 200);
        sigmas_std[j][i] = new TH2F(id_hist + " sigmas std", "", number_of_points, 0, 200, 200, -8, 8);
        sigmas_tri[j][i] = new TH2F(id_hist + " sigmas tri", "", number_of_points, 0, 200, 200, -8, 8);
        sigmas_tri_kin_fit[j][i] = new TH2F(id_hist + " sigmas kinfit", "", number_of_points, 0, 200, 50, -8, 8);
      }
      else if (i == 3)
      {
        neu_std_hist[j][i] = new TH1F(id_hist + " std", "", 100, 0, 5);
        neu_tri_hist[j][i] = new TH1F(id_hist + " tri", "", 100, 0, 5);
        neu_tri_kin_fit_hist[j][i] = new TH1F(id_hist + " kinfit", "", 100, 0, 5);
        res_std_hist[j][i] = new TH1F(id_hist + " std res", "", number_of_points, 0, 200);
        res_tri_hist[j][i] = new TH1F(id_hist + " tri res", "", number_of_points, 0, 200);
        res_tri_kin_fit_hist[j][i] = new TH1F(id_hist + " kinfit res", "", number_of_points, 0, 200);
        sigmas_std[j][i] = new TH2F(id_hist + " sigmas std", "", number_of_points, 0, 200, 200, -5, 10);
        sigmas_tri[j][i] = new TH2F(id_hist + " sigmas tri", "", number_of_points, 0, 200, 200, -5, 10);
        sigmas_tri_kin_fit[j][i] = new TH2F(id_hist + " sigmas kinfit", "", number_of_points, 0, 200, 50, -5, 10);
      }
      else if (i == 4)
      {
        neu_std_hist[j][i] = new TH1F(id_hist + " std", "", 100, 0, 5);
        neu_tri_hist[j][i] = new TH1F(id_hist + " tri", "", 100, 0, 5);
        neu_tri_kin_fit_hist[j][i] = new TH1F(id_hist + " kinfit", "", 100, 0, 5);
        res_std_hist[j][i] = new TH1F(id_hist + " std res", "", number_of_points, 0, 200);
        res_tri_hist[j][i] = new TH1F(id_hist + " tri res", "", number_of_points, 0, 200);
        res_tri_kin_fit_hist[j][i] = new TH1F(id_hist + " kinfit res", "", number_of_points, 0, 200);
        sigmas_std[j][i] = new TH2F(id_hist + " sigmas std", "", number_of_points, 0, 200, 200, -5, 5);
        sigmas_tri[j][i] = new TH2F(id_hist + " sigmas tri", "", number_of_points, 0, 200, 200, -5, 5);
        sigmas_tri_kin_fit[j][i] = new TH2F(id_hist + " sigmas kinfit", "", number_of_points, 0, 200, 50, -5, 5);
      }
    }
  }

  Int_t count_angle = 0, min_ind[6], min_ind_180[6], min_ind_dist[6], isEqual = 0, bad_clus = 0;
  Float_t vec_before[4], vec_after[4], boost[3], ip_before[4], ip_after[4], d_cl, angle_min_0, angle_min_180, angle[6], angle_180[6], cos_tmp, num, den, dist[6], dist_min, tcl, cos_boost_kaon_mom, angle_boost_kaon_mom;

  Bool_t cut = true, bunch_0 = true, bunch_plus1 = true, bunch_minus1 = true, blob = true, chosen_bunch = true;

  Int_t counts_cut = 0, counts_all = 0, counts_bad = 0, counts_good = 0, counts_goodzero = 0, counts_goodpone = 0, counts_goodmone = 0, g4taken_corr[4], counts_bad_smaller = 0, counts_bad_bigger = 0, counts_bad1 = 0, counts_bad2 = 0, counts_bad3 = 0, counts_bad4 = 0, counts_blob = 0, counts_rej = 0, counts_badblob = 0, counts_badzero = 0, counts_badpone = 0, counts_badmone = 0;

  Float_t bunch_corr = 0;

  // General statistics
  Int_t
      counts_sig_all = 0,
      counts_sig_err[17] = {0},
      counts_sig_cut[6] = {0},
      counts_sig_sel = 0;

  KLOE::pm00 event;

  for (Int_t i = 0; i < nentries; i++)
  {
    chain->GetEntry(i);

    if (baseKin.mcflag == 1 && (baseKin.mctruth == 1 || baseKin.mctruth == 2 || baseKin.mctruth == 0))
    {
      counts_sig_all++;

      if (errId != 2)
      {
        counts_sig_err[1]++;
        if (errId != 4)
        {
          counts_sig_err[3]++;
          if (errId != 5)
          {
            counts_sig_err[4]++;
            if (errId != 6)
            {
              counts_sig_err[5]++;
              if (errId != 7)
              {
                counts_sig_err[6]++;
                if (errId != 11)
                {
                  counts_sig_err[10]++;
                  if (cutId != 1)
                  {
                    counts_sig_cut[0]++;
                    if (cutId != 2)
                    {
                      counts_sig_cut[1]++;
                      if (cutId != 3)
                      {
                        counts_sig_cut[2]++;
                        if (cutId != 4)
                        {
                          counts_sig_cut[3]++;
                          if (cutId != 5)
                          {
                            counts_sig_cut[4]++;
                            if (cutId != 6)
                            {
                              counts_sig_cut[5]++;
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }

    count_angle = 0;

    if (baseKin.mcflag == 1 && (baseKin.mctruth == 1 || baseKin.mctruth == 2 || baseKin.mctruth == 0))
    {
      // Kaon path length


      if (baseKin.mctruth == 1 || baseKin.mctruth == 2 || baseKin.mctruth == 0)
      {

        lengthch_mc = sqrt(pow(baseKin.Kchmc[6] - baseKin.ipmc[0], 2) + pow(baseKin.Kchmc[7] - baseKin.ipmc[1], 2) + pow(baseKin.Kchmc[8] - baseKin.ipmc[2], 2));
        Rtch_mc = sqrt(pow(baseKin.Kchmc[6] - baseKin.ipmc[0], 2) + pow(baseKin.Kchmc[7] - baseKin.ipmc[1], 2));
        lengthneu_mc = sqrt(pow(baseKin.Knemc[6] - baseKin.ipmc[0], 2) + pow(baseKin.Knemc[7] - baseKin.ipmc[1], 2) + pow(baseKin.Knemc[8] - baseKin.ipmc[2], 2));
        Rtneu_mc = sqrt(pow(baseKin.Knemc[6] - baseKin.ipmc[0], 2) + pow(baseKin.Knemc[7] - baseKin.ipmc[1], 2));
      }
      else if (baseKin.mctruth == 6)
      {
        lengthch_mc = sqrt(pow(Kchmc_bcg[6] - ipmc_bcg[0], 2) + pow(Kchmc_bcg[7] - ipmc_bcg[1], 2) + pow(Kchmc_bcg[8] - ipmc_bcg[2], 2));
        lengthneu_mc = sqrt(pow(Knemc_bcg[6] - ipmc_bcg[0], 2) + pow(Knemc_bcg[7] - ipmc_bcg[1], 2) + pow(Knemc_bcg[8] - ipmc_bcg[2], 2));
        Rtneu_mc = sqrt(pow(Knemc_bcg[6] - ipmc_bcg[0], 2) + pow(Knemc_bcg[7] - ipmc_bcg[1], 2));
      }

      lengthch_rec = sqrt(pow(baseKin.Kchrec[6] - baseKin.ip[0], 2) + pow(baseKin.Kchrec[7] - baseKin.ip[1], 2) + pow(baseKin.Kchrec[8] - baseKin.ip[2], 2));

      lengthneu_tri = sqrt(pow(Knetri_kinfit[6] - ip_kinfit[0], 2) + pow(Knetri_kinfit[7] - ip_kinfit[1], 2) + pow(Knetri_kinfit[8] - ip_kinfit[2], 2));
      lengthneu_rec = sqrt(pow(baseKin.Knereclor[6] - baseKin.ip[0], 2) + pow(baseKin.Knereclor[7] - baseKin.ip[1], 2) + pow(baseKin.Knereclor[8] - baseKin.ip[2], 2));

      //
      // Kaon transverse radius
      // Rtneu_mc = sqrt(pow(baseKin.Knemc[6] - baseKin.ipmc[0], 2) + pow(baseKin.Knemc[7] - baseKin.ipmc[1], 2));
      Rtneu_tri = sqrt(pow(Knetri_kinfit[6] - ip_kinfit[0], 2) + pow(Knetri_kinfit[7] - ip_kinfit[1], 2));
      Rtneu_rec = sqrt(pow(baseKin.Knereclor[6] - baseKin.ip[0], 2) + pow(baseKin.Knereclor[7] - baseKin.ip[1], 2));

      Rtch_rec = sqrt(pow(baseKin.Kchrec[6] - baseKin.ip[0], 2) + pow(baseKin.Kchrec[7] - baseKin.ip[1], 2));
      //
      // Kaon velocity
      if (baseKin.mctruth == 1 || baseKin.mctruth == 2 || baseKin.mctruth == 0)
        v_Kchmc = cVel * baseKin.Kchmc[4] / baseKin.Kchmc[3];
      else
        v_Kchmc = cVel * Kchmc_bcg[4] / Kchmc_bcg[3];

      v_Kchrec = cVel * baseKin.Kchboost[4] / baseKin.Kchboost[3];

      if (baseKin.mctruth == 1 || baseKin.mctruth == 2 || baseKin.mctruth == 0)
        v_Kneumc = cVel * baseKin.Knemc[4] / baseKin.Knemc[3];
      else
        v_Kneumc = cVel * Knemc_bcg[4] / Knemc_bcg[3];

      v_Kneutri = cVel * sqrt(pow(Knetri_kinfit[0], 2) + pow(Knetri_kinfit[1], 2) + pow(Knetri_kinfit[2], 2)) / Knetri_kinfit[3];
      v_Kneurec = cVel * baseKin.Knereclor[4] / baseKin.Knereclor[3];
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
        boost[0] = -baseKin.phi_mom[0] / baseKin.phi_mom[3];
        boost[1] = -baseKin.phi_mom[1] / baseKin.phi_mom[3];
        boost[2] = -baseKin.phi_mom[2] / baseKin.phi_mom[3];
      }

      vec_before[0] = Knetri_kinfit[6] - ip_kinfit[0];
      vec_before[1] = Knetri_kinfit[7] - ip_kinfit[1];
      vec_before[2] = Knetri_kinfit[8] - ip_kinfit[2];
      vec_before[3] = Knetri_kinfit[9];

      lorentz_transf(boost, vec_before, vec_after);

      lengthneu_tri_CM = sqrt(pow(vec_after[0], 2) + pow(vec_after[1], 2) + pow(vec_after[2], 2));
      v_Kneutri_CM = lengthneu_tri_CM / vec_after[3];
      //!

      // Boost - kaon momentum angle

      cos_boost_kaon_mom = (baseKin.phi_mom[0] * Knetri_kinfit[0] + baseKin.phi_mom[1] * Knetri_kinfit[1] + baseKin.phi_mom[2] * Knetri_kinfit[2]) / (sqrt(pow(baseKin.phi_mom[0], 2) + pow(baseKin.phi_mom[1], 2) + pow(baseKin.phi_mom[2], 2)) * sqrt(pow(Knetri_kinfit[0], 2) + pow(Knetri_kinfit[1], 2) + pow(Knetri_kinfit[2], 2)));

      angle_boost_kaon_mom = 180. * acos(cos_boost_kaon_mom) / M_PI;

      for (Int_t j = 0; j < 2; j++)
      {
        if (j == 0)
          cut = chi2min < 10000000000.; // && abs(Knetri_kinfit[9] - t_neumc) < 2.375/2.;
        else
          cut = chi2min < 40.0; // && abs(Knetri_kinfit[9] - t_neumc) < 2.375/2.;

        if (cut && (errId > 3 || errId == 0) && done_kinfit == 1)
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

          Float_t
              trk_before[2][4],
              trk_after[2][4];

          for (Int_t ll = 0; ll < 2; ll++)
          {
            trk_before[ll][0] = baseKin.trk[ll][0];
            trk_before[ll][1] = baseKin.trk[ll][1];
            trk_before[ll][2] = baseKin.trk[ll][2];
            trk_before[ll][3] = baseKin.trk[ll][3];
          }

          boost[0] = -baseKin.Kchboost[0] / baseKin.Kchboost[3];
          boost[1] = -baseKin.Kchboost[1] / baseKin.Kchboost[3];
          boost[2] = -baseKin.Kchboost[2] / baseKin.Kchboost[3];

          lorentz_transf(boost, trk_before[0], trk_after[0]);
          lorentz_transf(boost, trk_before[1], trk_after[1]);

          Float_t 
                costrk_num = trk_after[0][0]*trk_after[1][0] + trk_after[0][1]*trk_after[1][1] + trk_after[0][2]*trk_after[1][2],
                costrk_den = sqrt(pow(trk_after[0][0], 2) + pow(trk_after[0][1], 2) + pow(trk_after[0][2], 2)) * sqrt(pow(trk_after[1][0], 2) + pow(trk_after[1][1], 2) + pow(trk_after[1][2], 2)),
                costrk = costrk_num / costrk_den;

          Bool_t
              cut_trcsum   = 1,//= trcsum > -1., 
              cut_Kneuminv = 1,//= (minv4gam - mK0) < 76,
              cut_Kchminv = 1,//= (baseKin.Kchrec[5] - mK0) < 1.2,
              cut_Qmiss = 1,//= sqrt(pow(baseKin.Kchboost[3] - baseKin.Kchrec[3], 2) +
                              //  pow(baseKin.Kchboost[0] - baseKin.Kchrec[0], 2) +
                              //  pow(baseKin.Kchboost[1] - baseKin.Kchrec[1], 2) +
                              //  pow(baseKin.Kchboost[2] - baseKin.Kchrec[2], 2)) < 3.75,
              cut_costrk = 1;//= costrk < -0.8;

          if (j == 0 /*&& (lengthneu_tri < 50 && lengthch_rec < 50)*/ && (errId > 3 || errId == 0) && done_kinfit == 1 && cut_trcsum && cut_Kneuminv && cut_Kchminv && cut_Qmiss && cut_costrk)
          { 
            counts_all++;

            // std::cout << counts_all << std::endl;

            if (Knetri_kinfit[9] > time_cut)
              counts_cut++;

            //! Check how many baseKin.cluster sets are found good and which bad

            // std::cout << int(baseKin.mctruth) << std::endl;

            for (Int_t k = 0; k < 4; k++)
            {
              g4taken_corr[k] = baseKin.ncll[g4taken_kinfit[k]] - 1;
              // std::cout << g4taken_corr[k] << " ";
            }

            // std::cout << std::endl;

            isEqual = event.ArrayEquality(clusindgood, g4taken_corr, 4);

            // for(Int_t k = 0; k < 4; k++)
            // {
            //   std::cout << clusindgood[k] << " " << g4taken_corr[k] << std::endl;
            // }

            bad_clus = isEqual;

            bunch_0 = abs(Knetri_kinfit[9] - t_neumc) < 2.715 / 2.;
            bunch_plus1 = Knetri_kinfit[9] - t_neumc < 3. * 2.715 / 2. && Knetri_kinfit[9] - t_neumc > 2.715 / 2.;
            bunch_minus1 = Knetri_kinfit[9] - t_neumc < -2.715 / 2. && Knetri_kinfit[9] - t_neumc > -3. * 2.715 / 2.;
            blob = Knetri_kinfit[9] - t_neumc > 3. * 2.715 / 2. && Knetri_kinfit[9] - t_neumc < 200.;

            if (isEqual == 0)
            {
              counts_good++;

              if (bunch_0)
                counts_goodzero++;
              else if (bunch_plus1)
                counts_goodpone++;
              else if (bunch_minus1)
                counts_goodmone++;
              else if (blob)
                counts_blob++;
            }
            else
            {
              counts_bad++;

              if (bunch_0)
                counts_badzero++;
              else if (bunch_plus1)
                counts_badpone++;
              else if (bunch_minus1)
                counts_badmone++;
              else if (blob)
                counts_badblob++;

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

          if (bunch_int == 0)
            chosen_bunch = bunch_0;
          else if (bunch_int == 1)
            chosen_bunch = bunch_plus1;
          else if (bunch_int == -1)
            chosen_bunch = bunch_minus1;
          else if (bunch_int == 2)
            chosen_bunch = blob;
          else
            chosen_bunch = true;

          if (isEqual >= 0 && chosen_bunch /*&& (lengthneu_tri < 50 && lengthch_rec < 50)*/ && (errId > 3 || errId == 0) && done_kinfit == 1 && cut_trcsum && cut_Kneuminv && cut_Kchminv && cut_Qmiss && cut_costrk)
          {
            chi2_hist[j]->Fill(chi2min);
            prob_hist[j]->Fill(TMath::Prob(chi2min, M));

            if (bunch_0)
              chi2_bunch_hist[j][0]->Fill(chi2min);
            else if (bunch_minus1)
              chi2_bunch_hist[j][1]->Fill(chi2min);
            else if (bunch_plus1)
              chi2_bunch_hist[j][2]->Fill(chi2min);
            else
              chi2_bunch_hist[j][3]->Fill(chi2min);

            chi2_no_bunch_hist[j]->Fill(chi2min_no_bunch);

            if (angle_boost_kaon_mom < 90.)
            {
              betadt_hist[j]->Fill(lengthneu_tri / (cVel * t_neurec));
              betapm_hist[j]->Fill(v_Kneutri / cVel);
            }
            else if (angle_boost_kaon_mom > 90.)
            {
              betadtover90_hist[j]->Fill(lengthneu_tri / (cVel * t_neurec));
              betapmover90_hist[j]->Fill(v_Kneutri / cVel);
            }

            beta1_hist[j]->Fill(v_Kneumc / cVel, v_Kneutri / cVel);

            beta2_hist[j]->Fill(v_Kneumc / cVel, lengthneu_tri / (cVel * t_neurec));

            bunch[j]->Fill(bunchnum);
          }

          for (Int_t k = 0; k < 3; k++)
            ip_coor[j][k]->Fill(baseKin.ipmc[k], baseKin.ip[k]);

          for (Int_t k = 0; k < 5; k++)
          {
            if (k < 3)
            {
              if (isEqual >= 0 && chosen_bunch /*&& (lengthneu_tri < 50 && lengthch_rec < 50)*/ && (errId > 3 || errId == 0) && done_kinfit == 1 && cut_trcsum && cut_Kneuminv && cut_Kchminv && cut_Qmiss && cut_costrk)
              {
                if (baseKin.mctruth == 1 || baseKin.mctruth == 2 || baseKin.mctruth == 0)
                {
                  neu_vtx_corr[j][k]->Fill(baseKin.Knemc[6 + k], Knetri_kinfit[6 + k]);
                  sigmas_std[j][k]->Fill(lengthneu_mc, Knetri_kinfit[6 + k] - baseKin.Knemc[6 + k]);
                  pulls[j][4 + k]->Fill(Knetri_kinfit[6 + k] - baseKin.Knemc[6 + k]);
                  pulls[j][k]->Fill(Knetri_kinfit[k] - baseKin.Knemc[k]);
                  neu_mom[j][k]->Fill(baseKin.Knemc[k], Knetri_kinfit[k]);
                }
                else
                {
                  neu_vtx_corr[j][k]->Fill(Knemc_bcg[6 + k], Knetri_kinfit[6 + k]);
                  sigmas_std[j][k]->Fill(lengthneu_mc, Knetri_kinfit[6 + k] - Knemc_bcg[6 + k]);
                  pulls[j][4 + k]->Fill(Knetri_kinfit[6 + k] - Knemc_bcg[6 + k]);
                  pulls[j][k]->Fill(Knetri_kinfit[k] - baseKin.Knemc[k]);
                  neu_mom[j][k]->Fill(Knemc_bcg[k], Knetri_kinfit[k]);
                }
              }
            }
            else if (k == 3)
            {
              if (isEqual >= 0 && chosen_bunch /*&& (lengthneu_tri < 50 && lengthch_rec < 50)*/ && (errId > 3 || errId == 0) && done_kinfit == 1 && cut_trcsum && cut_Kneuminv && cut_Kchminv && cut_Qmiss && cut_costrk)
              {
                if (baseKin.mctruth == 1 || baseKin.mctruth == 2 || baseKin.mctruth == 0)
                {
                  neu_vtx_corr[j][3]->Fill(t_neumc, t_neutri); // Knetri_kinfit[6 + k]);
                  sigmas_std[j][3]->Fill(lengthneu_mc, (t_neutri - t_neumc) / tau_S_nonCPT);
                  pulls[j][4 + k]->Fill((t_neutri - t_neumc) / tau_S_nonCPT);
                  pulls[j][k]->Fill(Knetri_kinfit[k] - baseKin.Knemc[k]);
                  neu_mom[j][k]->Fill(baseKin.Knemc[k], Knetri_kinfit[k]);
                }
                else
                {
                  neu_vtx_corr[j][3]->Fill(t_neumc, t_neutri); // Knetri_kinfit[6 + k]);
                  sigmas_std[j][3]->Fill(lengthneu_mc, (t_neutri - t_neumc) / tau_S_nonCPT);
                  pulls[j][4 + k]->Fill((t_neutri - t_neumc) / tau_S_nonCPT);
                  pulls[j][k]->Fill(Knetri_kinfit[k] - baseKin.Knemc[k]);
                  neu_mom[j][k]->Fill(Knemc_bcg[k], Knetri_kinfit[k]);
                }
              }
            }
            else
            {
              if (isEqual >= 0 && chosen_bunch /*&& (lengthneu_tri < 50 && lengthch_rec < 50)*/ && (errId > 3 || errId == 0) && done_kinfit == 1 && cut_trcsum && cut_Kneuminv && cut_Kchminv && cut_Qmiss && cut_costrk)
              {
                sigmas_std[j][4]->Fill(lengthneu_mc, (lengthneu_tri - lengthneu_mc));
              }
            }
          }

          if (isEqual >= 0 && chosen_bunch /*&& (lengthneu_tri < 50 && lengthch_rec < 50)*/ && (errId > 3 || errId == 0) && done_kinfit == 1 && cut_trcsum && cut_Kneuminv && cut_Kchminv && cut_Qmiss && cut_costrk)
          {
            pulls[j][8]->Fill(Rtneu_tri - Rtneu_mc);
            pulls[j][9]->Fill(lengthneu_tri - lengthneu_mc);
          }

          d_cl = sqrt(pow(baseKin.cluster[0][0] - baseKin.bhabha_vtx[0], 2) + pow(baseKin.cluster[1][0] - baseKin.bhabha_vtx[1], 2) + pow(baseKin.cluster[2][0] - baseKin.bhabha_vtx[2], 2));

          if (isEqual >= 0 && chosen_bunch /*&& (lengthneu_tri < 50 && lengthch_rec < 50)*/ && (errId > 3 || errId == 0) && done_kinfit == 1 && cut_trcsum && cut_Kneuminv && cut_Kchminv && cut_Qmiss && cut_costrk)
          {

            for (Int_t k = 0; k < baseKin.ntmc; k++)
              pidmc_hist[j]->Fill(baseKin.pidmc[k]);

            angle_vs_time[j][0]->Fill(angle_min_0, t_neutri);

            angle_vs_time[j][1]->Fill(angle_min_180, t_neutri);
            dist_vs_time[j][0]->Fill(dist_min, angle_min_0);
            dist_vs_time[j][1]->Fill(dist_min, angle_min_180);

            sol_err_hist[j]->Fill(sol1err, sol2err);
            chosen_hist[j]->Fill(chosen);

            for (Int_t k = 0; k < 4; k++)
            {
              tcl = gamma_kinfit[k][7] - (sqrt(pow(gamma_kinfit[k][4] - Knetri_kinfit[6], 2) + pow(gamma_kinfit[k][5] - Knetri_kinfit[7], 2) + pow(gamma_kinfit[k][6] - Knetri_kinfit[8], 2)) / cVel) - t_neutri;
              tcl_hist[j]->Fill(tcl);

              cluscorr_hist[j][0]->Fill(gamma_kinfit[0][4], gamma_kinfit[1][4]);
              cluscorr_hist[j][1]->Fill(gamma_kinfit[0][5], gamma_kinfit[1][5]);
            }

            clusenergy_hist[j]->Fill(Knetri_kinfit[5]);

            first_clus_hist[j]->Fill(baseKin.cluster[3][0] - (d_cl / cVel));

            mcisr_hist[j]->Fill(baseKin.mcisr);
          }
        }
        else
        {
          counts_rej++;
        }
      }
    }
    else
    {
      bunch_corr = 999.;
    }
  }

  std::ofstream log_file;
  log_file.open("tri.log");

  log_file << "All events: " << counts_sig_all << std::endl;

  for (Int_t i = 0; i < 17; i++)
  {
    log_file << "Events without err " << i + 1 << ": " << counts_sig_err[i] / (Float_t)counts_sig_all << std::endl;
  }

  for (Int_t i = 0; i < 6; i++)
  {
    log_file << "Events without cut " << i + 1 << ": " << counts_sig_cut[i] / (Float_t)counts_sig_all << std::endl;
  }

  log_file << "Events at the starting point of analysis: " << counts_sig_sel / (Float_t)counts_sig_all << std::endl
           << std::endl;

  log_file << "Rejected events: " << counts_rej << std::endl
           << std::endl;

  log_file << "Total signal events: " << counts_all << std::endl
           << std::endl;

  log_file << "Good baseKin.cluster sets: " << counts_good << std::endl;
  log_file << "Good baseKin.cluster sets for 0 bunch: " << counts_goodzero << std::endl;
  log_file << "Good baseKin.cluster sets for +1 bunch: " << counts_goodpone << std::endl;
  log_file << "Good baseKin.cluster sets for -1 bunch: " << counts_goodmone << std::endl;
  log_file << "Good baseKin.cluster sets for blob: " << counts_blob << std::endl;
  log_file << "Good baseKin.cluster sets eff: " << counts_good / (Float_t)counts_all << std::endl;

  log_file << "Bad baseKin.cluster sets: " << counts_bad << std::endl;
  log_file << "Bad baseKin.cluster sets for 0 bunch: " << counts_badzero << std::endl;
  log_file << "Bad baseKin.cluster sets for +1 bunch: " << counts_badpone << std::endl;
  log_file << "Bad baseKin.cluster sets for -1 bunch: " << counts_badmone << std::endl;
  log_file << "Bad baseKin.cluster sets for blob: " << counts_badblob << std::endl;
  log_file << "Bad clusters 1: " << counts_bad1 / (Float_t)counts_all << std::endl;
  log_file << "Bad clusters 2: " << counts_bad2 / (Float_t)counts_all << std::endl;
  log_file << "Bad clusters 3: " << counts_bad3 / (Float_t)counts_all << std::endl;
  log_file << "Bad clusters 4: " << counts_bad4 / (Float_t)counts_all << std::endl;
  log_file << "Bad baseKin.cluster smaller: " << counts_bad_smaller << std::endl;
  log_file << "Bad baseKin.cluster bigger: " << counts_bad_bigger << std::endl;
  log_file << "Bad baseKin.cluster sets eff: " << counts_bad / (Float_t)counts_all << std::endl;
  log_file << "Bad baseKin.cluster smaller eff: " << counts_bad_smaller / (Float_t)counts_all << std::endl;
  log_file << "Bad baseKin.cluster bigger eff: " << counts_bad_bigger / (Float_t)counts_all << std::endl;

  log_file.close();

  TString fit_stats[3];
  TPaveText *fit_text = new TPaveText(0.7, 0.7, 0.9, 0.9, "NDC");

  TF1 *triple_fit = new TF1("triple_gaus", triple_gaus, -400.0, 400.0, 9, 1);
  triple_fit->SetParNames("Norm1", "Avg1", "Std1", "Norm2", "Avg2", "Std2", "Norm3", "Avg3", "Std3");
  triple_fit->SetLineWidth(4);

  TFitResultPtr result;

  TString x_title, y_title, x_title_list[5] = {"X^{tri}_{neu,1} - X^{gen}_{neu,1} [cm]", "X^{tri}_{neu,2} - X^{gen}_{neu,2} [cm]", "X^{tri}_{neu,3} - X^{gen}_{neu,3} [cm]", "t_{K#rightarrow#pi^{0}#pi^{0}}^{tri} - t_{K#rightarrow#pi^{0}#pi^{0}}^{gen} [#tau_{S}]", "d^{tri}_{neu} - d^{gen}_{neu} [cm]"}, hist_title[5] = {"d^{gen}_{neu}#in[0,10] cm", "d^{gen}_{neu}#in(10,20] cm", "d^{gen}_{neu}#in(20,30] cm", "d^{gen}_{neu}#in(30,40] cm", "d^{gen}_{neu}#in(40,50] cm"};

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
  canvas[0][74]->Print(img_dir + id_canva + ext_img);

  //! Chi2 and Prob histos

  canvas[0][11]->SetRightMargin(0.15);
  canvas[0][11]->SetLeftMargin(0.17);
  canvas[0][11]->cd();

  chi2_hist[0]->GetYaxis()->SetMaxDigits(3);

  id_canva = "chi2_dist";
  x_title = "#chi^{2}_{test}";
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
  canvas[0][11]->Print(img_dir + id_canva + ext_img);
  //!
  //! Chi2 and Prob histos

  canvas[0][89]->SetRightMargin(0.15);
  canvas[0][89]->SetLeftMargin(0.17);
  canvas[0][89]->cd();

  chi2_bunch_hist[0][0]->GetYaxis()->SetMaxDigits(3);

  id_canva = "chi2_dist_bunches";
  x_title = "#chi^{2} [-]";
  y_title = "Counts";

  legend->AddEntry(chi2_bunch_hist[0][0], "Bunch 0", "l");
  legend->AddEntry(chi2_bunch_hist[0][1], "Bunch -1", "l");
  legend->AddEntry(chi2_bunch_hist[0][2], "Bunch 1", "l");
  legend->AddEntry(chi2_bunch_hist[0][3], "Others", "l");

  chi2_bunch_hist[0][0]->GetXaxis()->CenterTitle();
  chi2_bunch_hist[0][0]->GetYaxis()->CenterTitle();
  chi2_bunch_hist[0][0]->GetXaxis()->SetTitle(x_title);
  chi2_bunch_hist[0][0]->GetYaxis()->SetTitle(y_title);
  chi2_bunch_hist[0][0]->GetYaxis()->SetRangeUser(0.0, 1.3 * chi2_bunch_hist[0][0]->GetMaximum());
  chi2_bunch_hist[0][0]->SetLineColor(kGreen);
  chi2_bunch_hist[0][1]->SetLineColor(kBlue);
  chi2_bunch_hist[0][2]->SetLineColor(kRed);
  chi2_bunch_hist[0][3]->SetLineColor(kBlack);
  chi2_bunch_hist[0][0]->Draw();
  chi2_bunch_hist[0][1]->Draw("SAMES");
  chi2_bunch_hist[0][2]->Draw("SAMES");
  chi2_bunch_hist[0][3]->Draw("SAMES");

  legend->Draw();

  canvas[0][89]->Print(img_dir + id_canva + ext_img);

  legend->Clear();
  //!

  canvas[0][12]->SetRightMargin(0.15);
  canvas[0][12]->SetLeftMargin(0.17);
  canvas[0][12]->cd();
  canvas[0][12]->SetLogy();

  prob_hist[0]->GetYaxis()->SetMaxDigits(3);

  id_canva = "prob";
  x_title = Form("Prob(#chi^{2}_{test},%d)", M);
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
  canvas[0][12]->Print(img_dir + id_canva + ext_img);
  //!

  //!
  canvas[0][26]->SetRightMargin(0.15);
  canvas[0][26]->SetLeftMargin(0.17);
  canvas[0][26]->cd();
  canvas[0][26]->SetLogy();

  pidmc_hist[0]->GetYaxis()->SetMaxDigits(3);

  id_canva = "baseKin.pidmc";

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
  canvas[0][26]->Print(img_dir + id_canva + ext_img);
  //!
  canvas[0][27]->SetRightMargin(0.15);
  canvas[0][27]->SetLeftMargin(0.17);
  canvas[0][27]->cd();
  canvas[0][27]->SetLogy();

  mcisr_hist[0]->GetYaxis()->SetMaxDigits(3);

  id_canva = "baseKin.mcisr";

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
  canvas[0][27]->Print(img_dir + id_canva + ext_img);
  //!
  canvas[0][28]->SetRightMargin(0.15);
  canvas[0][28]->SetLeftMargin(0.17);
  canvas[0][28]->cd();
  canvas[0][28]->SetLogy();

  first_clus_hist[0]->GetYaxis()->SetMaxDigits(3);

  id_canva = "first_clus";

  first_clus_hist[0]->GetXaxis()->CenterTitle();
  first_clus_hist[0]->GetYaxis()->CenterTitle();
  // first_clus_hist[0]->GetXaxis()->SetTitle("baseKin.T0step1 [ns]");
  first_clus_hist[0]->GetYaxis()->SetTitle("Counts");
  if (canvas[0][28]->GetLogy())
    first_clus_hist[0]->GetYaxis()->SetRangeUser(1.0, 10. * first_clus_hist[0]->GetMaximum());
  else
    first_clus_hist[0]->GetYaxis()->SetRangeUser(0.0, 1.3 * first_clus_hist[0]->GetMaximum());
  first_clus_hist[0]->SetLineColor(kBlack);
  first_clus_hist[1]->SetLineColor(kRed);
  first_clus_hist[0]->Draw();
  first_clus_hist[1]->Draw("SAME");
  canvas[0][28]->Print(img_dir + id_canva + ext_img);
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
  // tcl_hist[0]->GetXaxis()->SetTitle("baseKin.T0step1 [ns]");
  tcl_hist[0]->GetYaxis()->SetTitle("Counts");
  if (canvas[0][29]->GetLogy())
    tcl_hist[0]->GetYaxis()->SetRangeUser(1.0, 10. * tcl_hist[0]->GetMaximum());
  else
    tcl_hist[0]->GetYaxis()->SetRangeUser(0.0, 1.3 * tcl_hist[0]->GetMaximum());
  tcl_hist[0]->SetLineColor(kBlack);
  tcl_hist[1]->SetLineColor(kRed);
  tcl_hist[0]->Draw();
  tcl_hist[1]->Draw("SAME");
  canvas[0][29]->Print(img_dir + id_canva + ext_img);
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

    canvas[j][70]->Print(img_dir + id_canva + ext_img);

    id_canva = "Angle180vstime" + std::to_string(j + 1);

    canvas[j][71]->SetRightMargin(0.15);
    canvas[j][71]->SetLeftMargin(0.17);
    canvas[j][71]->cd();

    angle_vs_time[j][1]->Draw("COLZ");

    canvas[j][71]->Print(img_dir + id_canva + ext_img);

    id_canva = "Distvsangle0" + std::to_string(j + 1);

    canvas[j][72]->SetRightMargin(0.15);
    canvas[j][72]->SetLeftMargin(0.17);
    canvas[j][72]->cd();

    dist_vs_time[j][0]->Draw("COLZ");

    canvas[j][72]->Print(img_dir + id_canva + ext_img);

    id_canva = "Distvsangle180" + std::to_string(j + 1);

    canvas[j][74]->SetRightMargin(0.15);
    canvas[j][74]->SetLeftMargin(0.17);
    canvas[j][74]->cd();

    dist_vs_time[j][1]->Draw("COLZ");

    canvas[j][74]->Print(img_dir + id_canva + ext_img);

    id_canva = "Solerr" + std::to_string(j + 1);

    canvas[j][73]->SetRightMargin(0.15);
    canvas[j][73]->SetLeftMargin(0.17);
    canvas[j][73]->cd();

    sol_err_hist[j]->Draw("COLZ");

    canvas[j][73]->Print(img_dir + id_canva + ext_img);

    id_canva = "Clusenergy" + std::to_string(j + 1);

    canvas[j][75]->SetRightMargin(0.15);
    canvas[j][75]->SetLeftMargin(0.17);
    canvas[j][75]->cd();

    clusenergy_hist[j]->Draw();

    canvas[j][75]->Print(img_dir + id_canva + ext_img);

    id_canva = "Cluscorr_x_y_" + std::to_string(j + 1);

    canvas[j][76]->SetRightMargin(0.15);
    canvas[j][76]->SetLeftMargin(0.17);
    canvas[j][76]->cd();

    cluscorr_hist[j][0]->Draw("COLZ");

    canvas[j][76]->Print(img_dir + id_canva + ext_img);

    id_canva = "Cluscorr_x_z_" + std::to_string(j + 1);

    canvas[j][77]->SetRightMargin(0.15);
    canvas[j][77]->SetLeftMargin(0.17);
    canvas[j][77]->cd();

    cluscorr_hist[j][1]->Draw("COLZ");

    canvas[j][77]->Print(img_dir + id_canva + ext_img);

    id_canva = "Time_projection" + std::to_string(j + 1);

    canvas[j][78]->SetRightMargin(0.15);
    canvas[j][78]->SetLeftMargin(0.17);
    canvas[j][78]->cd();

    neu_vtx_corr[j][3]->ProjectionY()->Draw("HIST");

    canvas[j][78]->Print(img_dir + id_canva + ext_img);

    id_canva = "Chi2_substraction" + std::to_string(j + 1);

    canvas[j][79]->SetRightMargin(0.15);
    canvas[j][79]->SetLeftMargin(0.17);
    canvas[j][79]->cd();

    chi2subst[j]->Add(chi2_no_bunch_hist[j], chi2_hist[j], 1.0, -1.0);

    chi2subst[j]->Draw("HIST");

    canvas[j][79]->Print(img_dir + id_canva + ext_img);

    id_canva = "Chi2_division" + std::to_string(j + 1);

    canvas[j][80]->SetRightMargin(0.15);
    canvas[j][80]->SetLeftMargin(0.17);
    canvas[j][80]->cd();

    chi2add[j]->Add(chi2_no_bunch_hist[j], chi2_hist[j], 1.0, 1.0);

    chi2div[j]->Divide(chi2subst[j], chi2add[j]);

    chi2div[j]->Draw("HIST");

    canvas[j][80]->Print(img_dir + id_canva + ext_img);

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

      canvas[j][i]->Print(img_dir + id_canva + ext_img);
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

      canvas[j][i + 4]->Print(img_dir + id_canva + ext_img);
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

      canvas[j][i + 7]->Print(img_dir + id_canva + ext_img);
    }
    //!
    gStyle->SetOptStat(0);
    canvas[j][3]->SetRightMargin(0.15);
    canvas[j][3]->SetLeftMargin(0.17);
    canvas[j][3]->cd();

    TLine line_1(0., -2.715 - (2.715 / 2.), 15., 15. - 2.715 - (2.715 / 2.));
    TLine line_2(0., -(2.715 / 2.), 15., 15. - (2.715 / 2.));
    TLine line_3(0., (2.715 / 2.), 15., 15. + (2.715 / 2.));
    TLine line_4(0., 2.715 + (2.715 / 2.), 15., 15. + 2.715 + (2.715 / 2.));

    line_1.SetLineWidth(2);
    line_2.SetLineWidth(2);
    line_3.SetLineWidth(2);
    line_4.SetLineWidth(2);

    id_canva = "time_neutral" + std::to_string(j + 1);
    x_title = Form("t_{neu}^{gen} [ns]");
    y_title = Form("t_{neu}^{tri} [ns]");

    neu_vtx_corr[j][3]->GetXaxis()->SetTitle(x_title);
    neu_vtx_corr[j][3]->GetYaxis()->SetTitle(y_title);
    neu_vtx_corr[j][3]->Draw("COLZ");
    /*line_1.Draw("SAME");
    line_2.Draw("SAME");
    line_3.Draw("SAME");
    line_4.Draw("SAME");*/

    canvas[j][3]->Print(img_dir + id_canva + ext_img);

    gStyle->SetOptStat(1);
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

    canvas[j][10]->Print(img_dir + id_canva + ext_img);
    //!

    gStyle->SetOptStat(0);
    for (Int_t i = 0; i < 10; i++)
    {
      canvas[j][i + 13]->cd();

      parameter[0] = pulls[j][i]->GetEntries();
      parameter[1] = pulls[j][i]->GetBinCenter(pulls[j][i]->GetMaximumBin());
      parameter[2] = pulls[j][i]->GetStdDev();

      triple_fit->SetParameters(0.15 * parameter[0], parameter[1], parameter[2], 0.15 * parameter[0], parameter[1] - 1.0, parameter[2], parameter[0], parameter[1] + 1.0, parameter[2]);

      if (i < 3)
      {
        triple_fit->SetParLimits(0, 0.0, 100.0 * parameter[0]);
        triple_fit->SetParLimits(3, 0.0, 100.0 * parameter[0]);
        triple_fit->SetParLimits(6, 0.0, 100.0 * parameter[0]);

        triple_fit->SetParLimits(1, parameter[1] - 1.0, parameter[1] + 1.0);
        triple_fit->SetParLimits(4, parameter[1] - 1.0, parameter[1] + 1.0);
        triple_fit->SetParLimits(7, parameter[1] - 1.0, parameter[1] + 1.0);

        triple_fit->SetParLimits(2, 0.01 * parameter[2], 5.0 * parameter[2]);
        triple_fit->SetParLimits(5, 0.01 * parameter[2], 5.0 * parameter[2]);
        triple_fit->SetParLimits(8, 0.01 * parameter[2], 5.0 * parameter[2]);
      }
      else if (i >= 3 && i < 7)
      {
        triple_fit->SetParLimits(0, 0.0, 100.0 * parameter[0]);
        triple_fit->SetParLimits(3, 0.0, 100.0 * parameter[0]);
        triple_fit->SetParLimits(6, 0.0, 100.0 * parameter[0]);

        triple_fit->SetParLimits(1, parameter[1] - 1.0, parameter[1] + 1.0);
        triple_fit->SetParLimits(4, parameter[1] - 1.0, parameter[1] + 1.0);
        triple_fit->SetParLimits(7, parameter[1] - 1.0, parameter[1] + 1.0);

        triple_fit->SetParLimits(2, 0.01 * parameter[2], 5.0 * parameter[2]);
        triple_fit->SetParLimits(5, 0.01 * parameter[2], 5.0 * parameter[2]);
        triple_fit->SetParLimits(8, 0.01 * parameter[2], 5.0 * parameter[2]);
      }
      else
      {
        triple_fit->SetParLimits(0, 0.0, 100.0 * parameter[0]);
        triple_fit->SetParLimits(3, 0.0, 100.0 * parameter[0]);
        triple_fit->SetParLimits(6, 0.0, 100.0 * parameter[0]);

        triple_fit->SetParLimits(1, parameter[1] - 1.0, parameter[1] + 1.0);
        triple_fit->SetParLimits(4, parameter[1] - 1.0, parameter[1] + 1.0);
        triple_fit->SetParLimits(7, parameter[1] - 1.0, parameter[1] + 1.0);

        triple_fit->SetParLimits(2, 0.01 * parameter[2], 5.0 * parameter[2]);
        triple_fit->SetParLimits(5, 0.01 * parameter[2], 5.0 * parameter[2]);
        triple_fit->SetParLimits(8, 0.01 * parameter[2], 5.0 * parameter[2]);
      }

      result = pulls[j][i]->Fit(triple_fit, "SF");

      if (result == 0 && j == 0)
      {
        fit_stats[1] = Form("Mean = %.2f#pm%.2f", comb_mean(result->GetParams(), result->GetErrors()), comb_mean_err(result->GetParams(), result->GetErrors()));
        fit_stats[2] = Form("Width = %.2f#pm%.2f", comb_std_dev(result->GetParams(), result->GetErrors()), comb_std_dev_err(result->GetParams(), result->GetErrors()));
        fit_text->AddText(fit_stats[1]);
        fit_text->AddText(fit_stats[2]);

        // if(i >= 0 && i < 4)
        //   properties["variables"]["Resolutions"]["momCharged"][i] = comb_std_dev(result->GetParams(), result->GetErrors());
        // else if(i >= 4 && i < 8)
        //   properties["variables"]["Resolutions"]["vtxCharged"][i - 4] = comb_std_dev(result->GetParams(), result->GetErrors());
        // else if( i == 8 )
        //   properties["variables"]["Resolutions"]["rhoCharged"] = comb_std_dev(result->GetParams(), result->GetErrors());
        // else if( i == 9 )
        //   properties["variables"]["Resolutions"]["pathCharged"] = comb_std_dev(result->GetParams(), result->GetErrors());
        if (i >= 0 && i < 4)
          properties["variables"]["Resolutions"]["momNeutral"][interface][i] = comb_std_dev(result->GetParams(), result->GetErrors());
        else if (i >= 4 && i < 8)
          properties["variables"]["Resolutions"]["vtxNeutral"][interface][i - 4] = comb_std_dev(result->GetParams(), result->GetErrors());
        else if (i == 8)
          properties["variables"]["Resolutions"]["rhoNeutral"][interface] = comb_std_dev(result->GetParams(), result->GetErrors());
        else if (i == 9)
          properties["variables"]["Resolutions"]["pathNeutral"][interface] = comb_std_dev(result->GetParams(), result->GetErrors());
      }
      else if (result == 1 && j == 0)
      {
        // if(i >= 0 && i < 4)
        //   properties["variables"]["Resolutions"]["momCharged"][i] = nullptr;
        // else if(i >= 4 && i < 8)
        //   properties["variables"]["Resolutions"]["vtxCharged"][i - 4] = nullptr;
        // else if( i == 8 )
        //   properties["variables"]["Resolutions"]["rhoCharged"] = nullptr;
        // else if( i == 9 )
        //   properties["variables"]["Resolutions"]["pathCharged"] = nullptr;
        if (i >= 0 && i < 4)
          properties["variables"]["Resolutions"]["momNeutral"][interface][i] = nullptr;
        else if (i >= 4 && i < 8)
          properties["variables"]["Resolutions"]["vtxNeutral"][interface][i - 4] = nullptr;
        else if (i == 8)
          properties["variables"]["Resolutions"]["rhoNeutral"][interface] = nullptr;
        else if (i == 9)
          properties["variables"]["Resolutions"]["pathNeutral"][interface] = nullptr;
      }

      pulls[j][i]->SetLineWidth(5);

      pulls[j][i]->GetXaxis()->CenterTitle(1);

      pulls[j][i]->GetYaxis()->SetMaxDigits(3);
      pulls[j][i]->GetYaxis()->CenterTitle(1);
      pulls[j][i]->GetYaxis()->SetRangeUser(0.0, 1.3 * pulls[j][i]->GetMaximum());

      pulls[j][i]->Draw();
      fit_text->Draw();

      id_canva = "pull" + std::to_string(j + 1) + std::to_string(i + 1);
      canvas[j][i + 13]->Print(img_dir + id_canva + ext_img);

      fit_text->Clear();
    }
    //!
    canvas[j][23]->cd();
    canvas[j][23]->SetLogy(1);
    bunch[j]->GetYaxis()->SetRangeUser(10.0, 10.3 * bunch[j]->GetMaximum());
    bunch[j]->SetLineColor(kBlack);
    bunch[j]->Draw();

    id_canva = "bunchnum" + std::to_string(j + 1);
    canvas[j][23]->Print(img_dir + id_canva + ext_root);

    canvas[j][25]->SetRightMargin(0.15);
    canvas[j][25]->SetLeftMargin(0.17);
    canvas[j][25]->cd();

    betadt_hist[j]->GetYaxis()->SetMaxDigits(3);

    id_canva = "betadt" + std::to_string(j + 1);
    x_title = "#beta^{p/E}_{K#rightarrow#pi^{0}#pi^{0}} [-]";
    y_title = "Counts";

    legend->AddEntry(betadt_hist[j], "#theta < 90^{#circ}", "l");
    legend->AddEntry(betadtover90_hist[j], "#theta > 90^{#circ}", "l");

    betadt_hist[j]->GetXaxis()->CenterTitle();
    betadt_hist[j]->GetYaxis()->CenterTitle();
    betadt_hist[j]->GetXaxis()->SetTitle(x_title);
    betadt_hist[j]->GetYaxis()->SetTitle(y_title);
    betadt_hist[j]->GetYaxis()->SetRangeUser(0.0, 1.3 * betadt_hist[j]->GetMaximum());
    betadt_hist[j]->SetLineColor(kBlack);
    betadtover90_hist[j]->SetLineColor(kRed);
    betadt_hist[j]->Draw();
    betadtover90_hist[j]->Draw("SAME");
    legend->Draw();
    canvas[j][25]->Print(img_dir + id_canva + ext_img);

    legend->Clear();

    canvas[j][81]->SetRightMargin(0.15);
    canvas[j][81]->SetLeftMargin(0.17);
    canvas[j][81]->cd();

    betapm_hist[j]->GetYaxis()->SetMaxDigits(3);

    id_canva = "betapm" + std::to_string(j + 1);
    x_title = "#beta^{p/E}_{K#rightarrow#pi^{0}#pi^{0}} [-]";
    y_title = "Counts";

    legend->AddEntry(betapm_hist[j], "#theta < 90^{#circ}", "l");
    legend->AddEntry(betapmover90_hist[j], "#theta > 90^{#circ}", "l");

    betapm_hist[j]->GetXaxis()->CenterTitle();
    betapm_hist[j]->GetYaxis()->CenterTitle();
    betapm_hist[j]->GetXaxis()->SetTitle(x_title);
    betapm_hist[j]->GetYaxis()->SetTitle(y_title);
    betapm_hist[j]->GetYaxis()->SetRangeUser(0.0, 1.3 * betapm_hist[j]->GetMaximum());
    betapm_hist[j]->SetLineColor(kBlack);
    betapmover90_hist[j]->SetLineColor(kRed);
    betapm_hist[j]->Draw();
    betapmover90_hist[j]->Draw("SAME");
    legend->Draw();
    canvas[j][81]->Print(img_dir + id_canva + ext_img);

    legend->Clear();

    canvas[j][26]->SetRightMargin(0.15);
    canvas[j][26]->SetLeftMargin(0.17);
    canvas[j][26]->cd();

    beta1_hist[j]->GetYaxis()->SetMaxDigits(3);

    id_canva = "beta1" + std::to_string(j + 1);
    x_title = "#beta^{gen}_{K#rightarrow#pi^{0}#pi^{0}} [-]";
    y_title = "#beta^{rec, p/E}_{K#rightarrow#pi^{0}#pi^{0}} [-]";

    beta1_hist[j]->GetXaxis()->CenterTitle();
    beta1_hist[j]->GetYaxis()->CenterTitle();
    beta1_hist[j]->GetXaxis()->SetTitle(x_title);
    beta1_hist[j]->GetYaxis()->SetTitle(y_title);
    beta1_hist[j]->Draw("COLZ");
    canvas[j][26]->Print(img_dir + id_canva + ext_img);

    canvas[j][26]->SetRightMargin(0.15);
    canvas[j][26]->SetLeftMargin(0.17);
    canvas[j][26]->cd();

    id_canva = "beta2" + std::to_string(j + 1);
    x_title = "#beta^{gen}_{K#rightarrow#pi^{0}#pi^{0}} [-]";
    y_title = "#beta^{rec, p/E}_{K#rightarrow#pi^{0}#pi^{0}} [-]";

    beta2_hist[j]->GetXaxis()->CenterTitle();
    beta2_hist[j]->GetYaxis()->CenterTitle();
    beta2_hist[j]->GetXaxis()->SetTitle(x_title);
    beta2_hist[j]->GetYaxis()->SetTitle(y_title);
    beta2_hist[j]->Draw("COLZ");
    canvas[j][26]->Print(img_dir + id_canva + ext_img);

    //!
    for (Int_t k = 0; k < 5; k++)
    {
      for (Int_t i = 1; i <= number_of_points; i++)
      {
        parameter[0] = sigmas_std[j][k]->ProjectionY("_py", i, i)->GetEntries();
        parameter[1] = sigmas_std[j][k]->ProjectionY("_py", i, i)->GetBinCenter(sigmas_std[j][k]->ProjectionY("_py", i, i)->GetMaximumBin());
        parameter[2] = sigmas_std[j][k]->ProjectionY("_py", i, i)->GetStdDev();

        triple_fit->SetParameters(0.15 * parameter[0], parameter[1], parameter[2], 0.15 * parameter[0], parameter[1] - 1.0, parameter[2], parameter[0], parameter[1] + 1.0, parameter[2]);

        triple_fit->SetParLimits(0, 0.0, 100.0 * parameter[0]);
        triple_fit->SetParLimits(3, 0.0, 100.0 * parameter[0]);
        triple_fit->SetParLimits(6, 0.0, 100.0 * parameter[0]);

        triple_fit->SetParLimits(1, parameter[1] - 2.0, parameter[1] + 2.0);
        triple_fit->SetParLimits(4, parameter[1] - 2.0, parameter[1] + 2.0);
        triple_fit->SetParLimits(7, parameter[1] - 2.0, parameter[1] + 2.0);

        triple_fit->SetParLimits(2, 0.1 * parameter[2], 6.0 * parameter[2]);
        triple_fit->SetParLimits(5, 0.1 * parameter[2], 6.0 * parameter[2]);
        triple_fit->SetParLimits(8, 0.1 * parameter[2], 6.0 * parameter[2]);

        TH1 *hist_temp = sigmas_std[j][k]->ProjectionY("_py", i, i);
        result = hist_temp->Fit(triple_fit, "LMES");

        res_std_hist[j][k]->SetBinContent(i, comb_std_dev(result->GetParams(), result->GetErrors()));
        res_std_hist[j][k]->SetBinError(i, comb_std_dev_err(result->GetParams(), result->GetErrors()));

        canvas_proj[j][k][i - 1]->SetRightMargin(0.15);
        canvas_proj[j][k][i - 1]->SetLeftMargin(0.17);
        canvas_proj[j][k][i - 1]->cd();

        id_canva = "std_devs_proj" + std::to_string(j + 1) + std::to_string(k + 1) + std::to_string(i);
        x_title = x_title_list[k];
        y_title = "Counts";

        hist_temp->GetYaxis()->SetMaxDigits(3);
        hist_temp->GetXaxis()->CenterTitle();
        hist_temp->GetYaxis()->CenterTitle();
        hist_temp->SetTitle(hist_title[i - 1]);
        hist_temp->GetXaxis()->SetTitle(x_title);
        hist_temp->GetYaxis()->SetTitle(y_title);

        hist_temp->Draw();
        canvas_proj[j][k][i - 1]->Print(img_dir + id_canva + ext_img);
      }

      canvas_std[j][k]->SetRightMargin(0.15);
      canvas_std[j][k]->SetLeftMargin(0.17);
      canvas_std[j][k]->cd();

      sigmas_std[j][k]->GetYaxis()->SetMaxDigits(3);

      id_canva = "std_devs" + std::to_string(j + 1) + std::to_string(k + 1);
      x_title = "#beta_{K#rightarrow#pi^{0}#pi^{0}} [-]";
      y_title = "Counts";

      sigmas_std[j][k]->GetXaxis()->CenterTitle();
      sigmas_std[j][k]->GetYaxis()->CenterTitle();
      sigmas_std[j][k]->GetXaxis()->SetTitle(x_title);
      sigmas_std[j][k]->GetYaxis()->SetTitle(y_title);
      sigmas_std[j][k]->GetYaxis()->SetRangeUser(-10.0, 10.0);
      sigmas_std[j][k]->Draw("COLZ");
      canvas_std[j][k]->Print(img_dir + id_canva + ext_img);

      res_std_hist[j][k]->SetLineColor(res_color[k]);
    }

    legend->AddEntry(res_std_hist[j][0], "x", "le");
    legend->AddEntry(res_std_hist[j][1], "y", "le");
    legend->AddEntry(res_std_hist[j][2], "z", "le");
    //!
    id_canva = "sigmas_std_coordinates" + std::to_string(j + 1);
    x_title = "K#rightarrow#pi^{0}#pi^{0} path generated [cm]";
    y_title = "#sigma(#vec{X}^{tri}_{K#rightarrow#pi^{0}#pi^{0},i} - #vec{X}^{gen}_{K#rightarrow#pi^{0}#pi^{0},i}) [cm]";

    canvas[j][22]->cd();
    canvas_std[j][22]->SetLeftMargin(0.2);
    res_std_hist[j][0]->SetXTitle(x_title);
    res_std_hist[j][0]->SetYTitle(y_title);
    res_std_hist[j][0]->GetYaxis()->SetRangeUser(0.0, 10.0);
    res_std_hist[j][0]->Draw("PE1");
    res_std_hist[j][1]->Draw("PE1SAME");
    res_std_hist[j][2]->Draw("PE1SAME");
    legend->Draw();
    canvas[j][22]->Print(img_dir + id_canva + ext_img);
    //!
    id_canva = "sigmas_std_times" + std::to_string(j + 1);
    x_title = "K#rightarrow#pi^{0}#pi^{0} path generated [cm]";
    y_title = "#sigma(t_{K#rightarrow#pi^{0}#pi^{0}}^{tri} - t_{K#rightarrow#pi^{0}#pi^{0}}^{gen}) [#tau_{S}]";
    canvas[j][99]->cd();
    canvas_std[j][99]->SetLeftMargin(0.2);
    res_std_hist[j][3]->SetXTitle(x_title);
    res_std_hist[j][3]->SetYTitle(y_title);
    res_std_hist[j][3]->GetYaxis()->SetMaxDigits(3);
    res_std_hist[j][3]->GetYaxis()->SetRangeUser(0.0, 5.0);
    res_std_hist[j][3]->Draw("PE1");
    canvas[j][99]->Print(img_dir + id_canva + ext_img);
    //!
    id_canva = "sigmas_std_length" + std::to_string(j + 1);
    x_title = "K#rightarrow#pi^{0}#pi^{0} path generated [cm]";
    y_title = "#sigma(d_{neu}^{tri} - d_{neu}^{gen}) [cm]";
    canvas[j][24]->cd();
    canvas_std[j][24]->SetLeftMargin(0.2);
    res_std_hist[j][4]->SetXTitle(x_title);
    res_std_hist[j][4]->SetYTitle(y_title);
    res_std_hist[j][4]->GetYaxis()->SetRangeUser(0.0, 5.0);
    res_std_hist[j][4]->Draw("PE1");
    canvas[j][24]->Print(img_dir + id_canva + ext_img);
    //!

    legend->Clear();
  }

  std::ofstream outfile(propName);
  outfile << properties.dump(4);
  outfile.close();

  return 0;
}
