#include <const.h>

#include <stdlib.h>
#include <ctime>
#include <fstream>

#include <iostream>
#include <sstream>

#include <TH1.h>
#include <TH2.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TString.h>
#include <TPrincipal.h>

#include <PdgManager.h>

#include <TGaxis.h>

using json = nlohmann::json;

namespace Paths
{
  std::string kloedataPath = getenv("KLOE_DBV26_DK0");
  std::string kloeMCPath = getenv("KLOE_DBV26_MK0");
  std::string workdirPath = getenv("WORKDIR");
  std::string chainDataFiles = kloedataPath + "/*.root";
  std::string chainMCFiles = kloeMCPath + "/*.root";
  std::string pdgConstFilePath = (std::string)getenv("PDGAPI") + "/pdg_const.json";
  std::string propertiesPath = getenv("PROPERTIESKLOE");
  std::string histogramConfigDir = propertiesPath + "/histogram_conf";
  std::string histogramConfig1DPath = histogramConfigDir + "/histogram1D.csv";
  std::string histogramConfig2DPath = histogramConfigDir + "/histogram2D.csv";
  std::string propName = propertiesPath + "/properties.json";
  std::string analysisConfigPath = propertiesPath + "/analysis_config.json";
  std::string rootfilesName = propertiesPath + "/root-files.json";
  std::string cutlimitsName = propertiesPath + "/cut-limits.json";
  std::string reportConfigPath = propertiesPath + "/report-config.json";

  TString base_path = workdirPath + "/KLOE/";
  TString path_tmp = "";
  TString path_cs = "";
  TString prod2root_path_v26 = "/data/k2/DBV-26/DK0";
  TString ext_root = ".root";
  TString ext_img = ".pdf";
  TString ext_csv = ".csv";

  TString gen_vars_dir = base_path + "Subanalysis/GeneratedVars/";
  TString neutrec_dir = base_path + "Subanalysis/Neutrec/";
  TString cpfit_dir = base_path + "Subanalysis/CPFit/";
  TString covmatrix_dir = base_path + "Subanalysis/CovarianceMatrix/";
  TString initialanalysis_dir = base_path + "Subanalysis/InitialAnalysis/";
  TString omegarec_dir = base_path + "Subanalysis/OmegaRec/";
  TString efficiency_dir = base_path + "Subanalysis/EfficiencyAnalysis/";
  TString charged_dir = base_path + "Subanalysis/KchRec/";
  TString plots_dir = base_path + "Subanalysis/Plots/";
  TString root_files_dir = "root_files/";
  TString input_dir = "input/";
  TString logs_dir = "log/";
  TString result_dir = "results/";
  TString img_dir = "img/";
}

namespace Filenames
{
  const TString gen_vars_filename = "gen_vars_";
  const TString mctruth_filename = "mctruth_";
  const TString neu_tri_filename = "neuvtx_tri_rec_";
  const TString neu_triangle_filename = "neuvtx_triangle_rec_";
  const TString neu_trilateration_kin_fit_filename = "neuvtx_tri_kin_fit_";
  const TString omega_rec_filename = "omega_rec_";
  const TString omega_rec_kin_fit_filename = "omega_rec_kin_fit_";
  const TString cut_vars_filename = "cut_vars_";

  TString gen_vars_tree;
  TString neutrec_triangle_tree;
  TString neutrec_tri_tree = "h_tri";
  TString neutrec_kin_fit_tree = "h_tri_kin_fit";
  TString omegarec_tree;
  TString omegarec_kin_fit_tree;
  TString cut_vars_tree = "h_cut_vars";

  TString omegaRecPath;
  TString mctruthPath;
  TString genvarsPath;
  TString trianglePath;
}

namespace PhysicsConstants
{
  // Constants used in the analysis
  // Fundamental quantities
  Double_t cVel = 29.9792458;       // cm/ns
  Double_t hBar = 6.582119569E-34;  // MeV*s
  Double_t eleCh = 1.602176634E-19; // C

  // Particles' masses
  Double_t mPhi = 1019.461;     // MeV/c^2
  Double_t mK0 = 4;             // MeV/c^2
  Double_t mPi0 = 134.9768;     // MeV/c^2
  Double_t mPiCh = 139.57039;   // MeV/c^2
  Double_t mMuon = 105.6583755; // MeV/c^2
  Double_t mElec = 0.510998950; // MeV/c^2
  Double_t mOmega = 782.66;     // MeV/c^2

  // Branching ratios
  // Phi
  Double_t br_phi_kskl = 0.339;
  Double_t br_phi_omegapi0 = 4.7E-5;

  // K-short
  Double_t br_ks_pi0pi0 = 0.3069;
  Double_t br_ks_pippim = 0.6920;
  Double_t br_ks_pippimgamma = 1.79E-3;
  Double_t br_ks_piele = 7.04E-4;
  Double_t br_ks_pimu = 4.56E-4;

  // K-long
  Double_t br_kl_pi0pi0 = 8.64E-4;
  Double_t br_kl_pippim = 1.967E-3;
  Double_t br_kl_pippimpi0 = 0.1254;
  Double_t br_kl_3pi0 = 0.1952;
  Double_t br_kl_piele = 0.4055;
  Double_t br_kl_pimu = 0.2704;

  // Kaons' properties and CPV
  Double_t tau_S_nonCPT = 0.89564E-1;     // ns
  Double_t tau_S_CPT = 0.8954E-1;         // ns
  Double_t tau_L = 51.16;                 // ns
  Double_t delta_mass_nonCPT = 0.5289E10; // hbar s^-1
  Double_t delta_mass_CPT = 0.5293E10;    // hbar s^-1
  Double_t mod_epsilon = 2.228E-3;
  Double_t Re = 1.66E-3;
  Double_t Im_nonCPT = -0.11 * (M_PI / 180.); // rad
  Double_t Im_CPT = -0.002;                   // deg
  Double_t phi_pm_nonCPT = 43.4;              // deg
  Double_t phi_pm_CPT = 43.51;                // deg
  Double_t phi_00_nonCPT = 43.7;              // deg
  Double_t phi_00_CPT = 43.52;                // deg
}

namespace KLOE
{
  // General
  const Double_t T0 = 2.715; // ns
  const UInt_t MaxNumtrkv = 200;
  const UInt_t MIN_CLU_ENE = 20;

  Int_t firstFileMax = 1;
  Int_t lastFileMax = 1;
  Int_t numOfThreads = 1;

  const std::map<Int_t, TString> channName = {
      {0, "Data"},
      {1, "Signal"},
      {2, "Regeneration"},
      {3, "Omega"},
      {4, "3pi0"},
      {5, "Semileptonic"},
      {6, "Other"},
      {7, "pi+pi-pi+pi-"},
      {8, "MC sum"}};

  const std::map<TString, TString> channTitle = {
      {"Data", "Data"},
      {"Signal", "K_{S}K_{L}#rightarrow#pi^{+}#pi^{-}#pi^{0}#pi^{0}"},
      {"Regeneration", "Regeneration"},
      {"Omega", "#omega#pi^{0}#rightarrow#pi^{+}#pi^{-}#pi^{0}#pi^{0}"},
      {"3pi0", "K_{S}K_{L}#rightarrow#pi^{+}#pi^{-}3#pi^{0}"},
      {"Semileptonic", "K_{S}K_{L}#rightarrow#pi^{#pm}l^{#mp}#nu#pi^{0}#pi^{0}"},
      {"Other", "Other bcg"},
      {"pi+pi-pi+pi-", "K_{S}K_{L}#rightarrow#pi^{+}#pi^{-}#pi^{+}#pi^{-}"},
      {"MC sum", "MC sum"}};

  const std::map<TString, Color_t> channColor = {
      {"Data", kBlack},
      {"Signal", kRed},
      {"Regeneration", kGreen},
      {"Omega", kViolet},
      {"3pi0", kCyan},
      {"Semileptonic", kBlue},
      {"Other", kGreen - 1},
      {"pi+pi-pi+pi-", kYellow},
      {"MC sum", kOrange}};

  std::map<Int_t, Int_t> channEventCount = {
      {0, 0},
      {1, 0},
      {2, 0},
      {3, 0},
      {4, 0},
      {5, 0},
      {6, 0},
      {7, 0},
      {8, 0}};

  Int_t TotalCountMC()
  {
    Int_t total = 0;
    for (const auto &pair : channEventCount)
    {
      if (pair.first != 0 && pair.first != 8) // Pomijamy dane
        total += pair.second;
    }
    return total;
  }

  namespace Histograms
  {
    const std::vector<TString> varNames = {
        "px_Pi1", "py_Pi1", "pz_Pi1", "Energy_Pi1",
        "px_Pi2", "py_Pi2", "pz_Pi2", "Energy_Pi2",
        "px_Kch", "py_Kch", "pz_Kch", "Energy_Kch",
        "px_Kne", "py_Kne", "pz_Kne", "Energy_Kne",
        "px_phi", "py_phi", "pz_phi", "Energy_phi",
        "mass_Kch", "mass_Kne", "mass_phi", "mass_pi01", "mass_pi02", "combined_mass_pi0", "mass_omega", "mass_omega_rec",
        "vKne", "deltaPhiv", "deltaPhivFit", "Qmiss",
        "chi2_signalKinFit", "chi2_trilaterationKinFit", "prob_signal",
        "curv1", "phiv1", "cotv1",
        "curv2", "phiv2", "cotv2",
        "vtxNeu_x", "vtxNeu_y", "vtxNeu_z",
        "vtxNeu_x_Fit", "vtxNeu_y_Fit", "vtxNeu_z_Fit",
        "phi_vtx_x", "phi_vtx_y", "phi_vtx_z",
        "time_neutral_MC", "delta_t",
        "pull1", "pull2", "pull3", "pull4", "pull5", "pull6", "pull7", "pull8", "pull9", "pull10", "pull11", "pull12", "pull13", "pull14", "pull15", "pull16", "pull17", "pull18", "pull19", "pull20", "pull21", "pull22", "pull23", "pull24", "pull25", "pull26", "pull27", "pull28", "pull29", "pull30", "pull31", "pull32", "pull33", "pull34", "pull35", "pull36",
        "openingAngleCharged", "openingAngleNeutral",
        "nev", "nrun", "TransvRadius", "T0Omega"};

    // Konfiguracje histogramów 1D
    const std::map<TString, HistConfig1D> histConfigs1D = {
        // Pędy cząstek
        {"px_Pi1", {100, -100., 100., "Pion 1 p_{x}", "p_{x} [MeV/c]", "Counts"}},
        {"py_Pi1", {100, -100., 100., "Pion 1 p_{y}", "p_{y} [MeV/c]", "Counts"}},
        {"pz_Pi1", {100, -100., 100., "Pion 1 p_{z}", "p_{z} [MeV/c]", "Counts"}},
        {"Energy_Pi1", {100, -50.0, 50.0, "Pion 1 Energy", "E [MeV]", "Counts"}},

        {"px_Pi2", {100, -100., 100., "Pion 2 p_{x}", "p_{x} [MeV/c]", "Counts"}},
        {"py_Pi2", {100, -100., 100., "Pion 2 p_{y}", "p_{y} [MeV/c]", "Counts"}},
        {"pz_Pi2", {100, -100., 100., "Pion 2 p_{z}", "p_{z} [MeV/c]", "Counts"}},
        {"Energy_Pi2", {100, -50.0, 50.0, "Pion 2 Energy", "E [MeV]", "Counts"}},

        {"px_Kch", {100, -20., 20., "Kch p_{x}", "p_{x} [MeV/c]", "Counts"}},
        {"py_Kch", {100, -20., 20., "Kch p_{y}", "p_{y} [MeV/c]", "Counts"}},
        {"pz_Kch", {100, -30., 30., "Kch p_{z}", "p_{z} [MeV/c]", "Counts"}},
        {"Energy_Kch", {100, -10.0, 10.0, "Kch Energy", "E [MeV]", "Counts"}},

        {"px_Kne", {100, -100., 100., "Kne p_{x}", "p_{x} [MeV/c]", "Counts"}},
        {"py_Kne", {100, -100., 100., "Kne p_{y}", "p_{y} [MeV/c]", "Counts"}},
        {"pz_Kne", {100, -100., 100., "Kne p_{z}", "p_{z} [MeV/c]", "Counts"}},
        {"Energy_Kne", {100, -50.0, 50.0, "Kne Energy", "E [MeV]", "Counts"}},

        {"px_Phi", {100, -100., 100., "Phi p_{x}", "p_{x} [MeV/c]", "Counts"}},
        {"py_Phi", {100, -100., 100., "Phi p_{y}", "p_{y} [MeV/c]", "Counts"}},
        {"pz_Phi", {100, -100., 100., "Phi p_{z}", "p_{z} [MeV/c]", "Counts"}},
        {"Energy_Phi", {100, -50.0, 50.0, "Phi Energy", "E [MeV]", "Counts"}},

        // Masy
        {"mass_Kch", {100, -5., 5., "Kaon Mass", "m^{inv}_{#pi^{+}#pi^{-}} - m_{K^{0}} [MeV/c^{2}]", "Counts"}},
        {"mass_Kne", {100, -200., 200., "Kaon Mass", "m^{inv}_{4#gamma} - m_{K^{0}} [MeV/c^{2}]", "Counts"}},
        {"mass_phi", {100, 1010., 1030., "#phi Meson Mass", "m_{#phi} [MeV/c^{2}]", "Counts"}},
        {"mass_pi01", {100, -20., 20., "#pi^{0} Mass 1", "m^{inv}_{2#gamma,1} - m_{#pi^{0}} [MeV/c^{2}]", "Counts"}},
        {"mass_pi02", {100, -20., 20., "#pi^{0} Mass 2", "m^{inv}_{2#gamma,2} - m_{#pi^{0}} [MeV/c^{2}]", "Counts"}},
        {"combined_mass_pi0", {100, 0., 100., "#pi^{0} Combined", "Comb^{#pi^{0}}_{err} [MeV/c^{2}]", "Counts"}},
        {"mass_omega", {100, -200., 100., "#omega Meson Mass", "m^{inv}_{#pi^{+}#pi^{-}#pi^{0}} [MeV/c^{2}]", "Counts"}},
        {"mass_omega_rec", {100, -200., 50., "#omega Meson Mass Rec", "m^{inv}_{#pi^{+}#pi^{-}#pi^{0}} [MeV/c^{2}]", "Counts"}},

        // Chi-square
        {"chi2_signalKinFit", {100, 0., 10., "Signal Kinematic Fit #chi^{2}", "#chi^{2}_{signal}", "Counts"}},
        {"prob_signal", {100, 0., 1., "Signal Kinematic Fit Probability", "Prob_{signal}", "Counts"}},
        {"chi2_trilaterationKinFit", {100, 0., 15., "Trilateration Fit #chi^{2}", "#chi^{2}_{#omega#pi^{0}}", "Counts"}},

        // Phi vertex
        {"phi_vtx_x", {100, -1.0, 1.0, "#phi Vertex x", "V^{#phi}_{x} [cm]", "Counts"}},
        {"phi_vtx_y", {100, -0.02, 0.02, "#phi Vertex y", "V^{#phi}_{y} [cm]", "Counts"}},
        {"phi_vtx_z", {100, -5.0, 5.0, "#phi Vertex z", "V^{#phi}_{z} [cm]", "Counts"}},

        // Vertex
        {"vtxNeu_x", {100, -4., 4., "Neutral Vertex x", "x [cm]", "Counts"}},
        {"vtxNeu_y", {100, -4., 4., "Neutral Vertex y", "y [cm]", "Counts"}},
        {"vtxNeu_z", {100, -10., 10., "Neutral Vertex z", "z [cm]", "Counts"}},

        // Vertex
        {"vtxNeu_x_Fit", {100, -4., 4., "Neutral Vertex x Fit", "V^{K_{ne}}_{x} [cm]", "Counts"}},
        {"vtxNeu_y_Fit", {100, -4., 4., "Neutral Vertex y Fit", "V^{K_{ne}}_{y} [cm]", "Counts"}},
        {"vtxNeu_z_Fit", {100, -10., 10., "Neutral Vertex z Fit", "V^{K_{ne}}_{z} [cm]", "Counts"}},

        // Path i Radius
        {"TransvRadius", {100, 0., 50., "Transversal Radius of Kaons", "R [cm]", "Counts"}},

        // T0 Omega
        {"T0Omega", {100, 0., 350., "T0 from Omega", "T0 [MeV]", "Counts"}},

        // Pull Signal Fit
        {"pull1", {100, -5., 5., "Pull 1 Signal Fit", "Pull", "Counts"}},
        {"pull2", {100, -5., 5., "Pull 2 Signal Fit", "Pull", "Counts"}},
        {"pull3", {100, -5., 5., "Pull 3 Signal Fit", "Pull", "Counts"}},
        {"pull4", {100, -5., 5., "Pull 4 Signal Fit", "Pull", "Counts"}},
        {"pull5", {100, -5., 5., "Pull 5 Signal Fit", "Pull", "Counts"}},
        {"pull6", {100, -5., 5., "Pull 6 Signal Fit", "Pull", "Counts"}},
        {"pull7", {100, -5., 5., "Pull 7 Signal Fit", "Pull", "Counts"}},
        {"pull8", {100, -5., 5., "Pull 8 Signal Fit", "Pull", "Counts"}},
        {"pull9", {100, -5., 5., "Pull 9 Signal Fit", "Pull", "Counts"}},
        {"pull10", {100, -5., 5., "Pull 10 Signal Fit", "Pull", "Counts"}},
        {"pull11", {100, -5., 5., "Pull 11 Signal Fit", "Pull", "Counts"}},
        {"pull12", {100, -5., 5., "Pull 12 Signal Fit", "Pull", "Counts"}},
        {"pull13", {100, -5., 5., "Pull 13 Signal Fit", "Pull", "Counts"}},
        {"pull14", {100, -5., 5., "Pull 14 Signal Fit", "Pull", "Counts"}},
        {"pull15", {100, -5., 5., "Pull 15 Signal Fit", "Pull", "Counts"}},
        {"pull16", {100, -5., 5., "Pull 16 Signal Fit", "Pull", "Counts"}},
        {"pull17", {100, -5., 5., "Pull 17 Signal Fit", "Pull", "Counts"}},
        {"pull18", {100, -5., 5., "Pull 18 Signal Fit", "Pull", "Counts"}},
        {"pull19", {100, -5., 5., "Pull 19 Signal Fit", "Pull", "Counts"}},
        {"pull20", {100, -5., 5., "Pull 20 Signal Fit", "Pull", "Counts"}},
        {"pull21", {100, -5., 5., "Pull 21 Signal Fit", "Pull", "Counts"}},
        {"pull22", {100, -5., 5., "Pull 22 Signal Fit", "Pull", "Counts"}},
        {"pull23", {100, -5., 5., "Pull 23 Signal Fit", "Pull", "Counts"}},
        {"pull24", {100, -5., 5., "Pull 24 Signal Fit", "Pull", "Counts"}},
        {"pull25", {100, -5., 5., "Pull 25 Signal Fit", "Pull", "Counts"}},
        {"pull26", {100, -5., 5., "Pull 26 Signal Fit", "Pull", "Counts"}},
        {"pull27", {100, -5., 5., "Pull 27 Signal Fit", "Pull", "Counts"}},
        {"pull28", {100, -5., 5., "Pull 28 Signal Fit", "Pull", "Counts"}},
        {"pull29", {100, -5., 5., "Pull 29 Signal Fit", "Pull", "Counts"}},
        {"pull30", {100, -5., 5., "Pull 30 Signal Fit", "Pull", "Counts"}},
        {"pull31", {100, -5., 5., "Pull 31 Signal Fit", "Pull", "Counts"}},
        {"pull32", {100, -5., 5., "Pull 32 Signal Fit", "Pull", "Counts"}},
        {"pull33", {100, -5., 5., "Pull 33 Signal Fit", "Pull", "Counts"}},
        {"pull34", {100, -5., 5., "Pull 34 Signal Fit", "Pull", "Counts"}},
        {"pull35", {100, -5., 5., "Pull 35 Signal Fit", "Pull", "Counts"}},
        {"pull36", {100, -5., 5., "Pull 36 Signal Fit", "Pull", "Counts"}},

        // Czasy
        {"time_neutral_MC", {100, -1., 1., "Neutral Kaon Time (MC)", "Trc_{sum} [ns]", "Counts"}},
        {"delta_t", {161, -40., 40., "Time Difference", "#Deltat [#tau_{S}]", "Counts"}},

        // Kąty i zmienne kinematyczne
        {"openingAngleCharged", {100, 0., 180., "Charged Particles Opening Angle", "#alpha [deg]", "Counts"}},
        {"openingAngleNeutral", {100, 0., 180., "Neutral Particles Opening Angle", "#alpha [deg]", "Counts"}},
        {"Qmiss", {100, 0, 15., "Missing 4-momentum", "Q_{miss} [MeV]", "Counts"}},

        // Track parameters
        {"curv1", {100, -2.0, 2.0, "Track 1 Curvature", "1/p_{T} [1/MeV]", "Counts"}},
        {"phiv1", {100, -0.5, 0.5, "Track 1 #phi", "#phi [rad]", "Counts"}},
        {"cotv1", {100, -1., 1., "Track 1 cot(#theta)", "cot(#theta)", "Counts"}},
        {"curv2", {100, -2.0, 2.0, "Track 2 Curvature", "1/p_{T} [1/MeV]", "Counts"}},
        {"phiv2", {100, -0.5, 0.5, "Track 2 #phi", "#phi [rad]", "Counts"}},
        {"cotv2", {100, -1., 1., "Track 2 cot(#theta)", "cot(#theta)", "Events"}},
        {"deltaPhiv", {100, 2., 4., "#phi_{+} - #phi_{-}", "#Delta#phi_{+-} [rad]", "Counts"}},
        {"deltaPhivFit", {100, 2., 4., "#phi_{+} - #phi_{-} Fit", "#Delta#phi^{fit}_{+-} [rad]", "Counts"}}};

    // Konfiguracje histogramów 2D
    const std::map<TString, HistConfig2D> histConfigs2D = {
        // Korelacje mas
        {"mass_pi01_vs_mass_pi02", {50, 50, 120., 150., 120., 150., "#pi^{0} Mass Correlation", "m_{#pi^{0}_{1}} [MeV/c^{2}]", "m_{#pi^{0}_{2}} [MeV/c^{2}]", ""}},

        {"mass_Kch_vs_mass_Kne", {50, 50, 480., 520., 480., 520., "Kaon Mass Correlation", "m_{K^{#pm}} [MeV/c^{2}]", "m_{K^{0}} [MeV/c^{2}]", ""}},

        // Vertex pozycje
        {"vtxNeu_x_vs_vtxNeu_y", {100, 100, -20., 20., -20., 20., "Neutral Vertex Position", "x [cm]", "y [cm]", ""}},

        {"vtxNeu_x_vs_vtxNeu_z", {100, 100, -20., 20., -50., 50., "Neutral Vertex x vs z", "x [cm]", "z [cm]", ""}},

        // Pędy vs masy
        {"Energy_Pi1_vs_mass_pi01", {50, 50, 0., 1000., 120., 150., "Pion Energy vs #pi^{0} Mass", "E_{#pi^{#pm}} [MeV]", "m_{#pi^{0}} [MeV/c^{2}]", ""}},

        // Chi-square korelacje
        {"chi2_signalKinFit_vs_chi2_trilaterationKinFit", {50, 50, 0., 100., 0., 100., "#chi^{2} Correlation", "#chi^{2}_{signal}", "#chi^{2}_{tri}", ""}},

        // Track korelacje
        {"curv1_vs_curv2", {50, 50, -0.01, 0.01, -0.01, 0.01, "Track Curvature Correlation", "1/p_{T1} [1/MeV]", "1/p_{T2} [1/MeV]", ""}},

        {"phiv1_vs_phiv2", {50, 50, -TMath::Pi(), TMath::Pi(), -TMath::Pi(), TMath::Pi(), "Track #phi Correlation", "#phi_{1} [rad]", "#phi_{2} [rad]", ""}},

        // Czas vs Czas MC
        {"delta_t_vs_delta_t_mc", {161, 81, -40., 40., -20., 20., "#Deltat_{rec} vs. #Deltat_{mc}", "#Deltat_{MC} [#tau_{S}]", "(#Deltat_{rec} - #Deltat_{MC}) [#tau_{S}]", ""}},
        {"delta_t_fit_vs_delta_t_mc", {161, 81, -40., 40., -20., 20., "#Deltat_{fit} vs. #Deltat_{mc}", "#Deltat_{MC} [#tau_{S}]", "(#Deltat_{fit} - #Deltat_{MC}) [#tau_{S}]", ""}},
        {"t00_tri_vs_t00_mc", {100, 100, -5., 300., 0., 300., "t^{00}_{tri} vs. t^{00}_{MC}", "t^{00}_{MC} [#tau_{S}]", "t^{00}_{tri} [#tau_{S}]", ""}},
        {"t00_triangle_vs_t00_mc", {100, 100, -5., 300., 0., 300., "t^{00}_{triangle} vs. t^{00}_{MC}", "t^{00}_{MC} [#tau_{S}]", "t^{00}_{triangle} [#tau_{S}]", ""}},
        {"t00_triangle_vs_t00_tri", {100, 100, -5., 300., 0., 300., "t^{00}_{triangle} vs. t^{00}_{tri}", "t^{00}_{triangle} [#tau_{S}]", "t^{00}_{tri} [#tau_{S}]", ""}}};

    // Mapowanie nazw 2D na pary zmiennych (dla automatycznego tworzenia)
    const std::map<TString, std::pair<TString, TString>> histConfigs2D_Variables = {
        {"mass_correlation", {"mass_pi01", "mass_pi02"}},
        {"kaon_masses", {"mass_Kch", "mass_Kne"}},
        {"vertex_xy", {"vtxNeu_x", "vtxNeu_y"}},
        {"vertex_xz", {"vtxNeu_x", "vtxNeu_z"}},
        {"energy_mass", {"Energy_Pi1", "mass_pi01"}},
        {"chi2_correlation", {"chi2_signalKinFit", "chi2_trilaterationKinFit"}},
        {"track_curvature", {"curv1", "curv2"}},
        {"track_phi", {"phiv1", "phiv2"}},
        {"delta_t_vs_delta_t_mc", {"delta_t", "delta_t_mc"}},
        {"delta_t_fit_vs_delta_t_mc", {"delta_t_fit", "delta_t_mc"}}};

    // Funkcja pomocnicza do usuwania cudzysłowów
    TString RemoveQuotes(const TString &str)
    {
      TString result = str;
      if (result = "")
        return result;
      // Usuń cudzysłowy z początku
      if (result.BeginsWith("\""))
        result.Remove(0, 1);
      // Usuń cudzysłowy z końca
      if (result.EndsWith("\""))
        result.Remove(result.Length() - 1, 1);
      return result;
    }

    std::map<TString, HistConfig1D> LoadHistogramConfigs1D(const TString &filePath)
    {
      std::map<TString, HistConfig1D> configs;
      std::ifstream file(filePath.Data());
      if (!file.is_open())
      {
        std::cerr << "ERROR: Cannot open histogram configuration file: " << filePath << std::endl;
        return configs;
      }

      TString line;
      while (line.ReadLine(file))
      {
        if (line.BeginsWith("#") || line.IsWhitespace())
          continue;

        TObjArray *tokens = line.Tokenize(";");

        if (tokens->GetEntries() < 7)
        {
          std::cerr << "WARNING: Invalid line format: " << line << std::endl;
          delete tokens;
          continue;
        }

        // Pobierz tokeny, usuń białe znaki i cudzysłowy
        TString name = TString(((TObjString *)tokens->At(0))->String()).Strip(TString::kBoth);
        name.ReplaceAll("\"", "");

        TString nBinsStr = TString(((TObjString *)tokens->At(1))->String()).Strip(TString::kBoth);
        TString xMinStr = TString(((TObjString *)tokens->At(2))->String()).Strip(TString::kBoth);
        TString xMaxStr = TString(((TObjString *)tokens->At(3))->String()).Strip(TString::kBoth);

        TString title = TString(((TObjString *)tokens->At(4))->String()).Strip(TString::kBoth);
        title.ReplaceAll("\"", "");

        TString xLabel = TString(((TObjString *)tokens->At(5))->String()).Strip(TString::kBoth);
        xLabel.ReplaceAll("\"", "");

        TString yLabel = TString(((TObjString *)tokens->At(6))->String()).Strip(TString::kBoth);
        yLabel.ReplaceAll("\"", "");

        Int_t nBins = nBinsStr.Atoi();
        Double_t xMin = xMinStr.Atof();
        Double_t xMax = xMaxStr.Atof();

        configs[name] = {nBins, xMin, xMax, title, xLabel, yLabel};

        delete tokens;
      }

      file.close();
      return configs;
    }

    std::map<TString, HistConfig2D> LoadHistogramConfigs2D(const TString &filePath)
    {
      std::map<TString, HistConfig2D> configs;
      std::ifstream file(filePath.Data());
      if (!file.is_open())
      {
        std::cerr << "ERROR: Cannot open histogram configuration file: " << filePath << std::endl;
        return configs;
      }

      TString line;
      while (line.ReadLine(file))
      {
        if (line.BeginsWith("#") || line.IsWhitespace())
          continue;

        TObjArray *tokens = line.Tokenize(";");

        if (tokens->GetEntries() < 10)
        {
          std::cerr << "WARNING: Invalid line format: " << line << std::endl;
          delete tokens;
          continue;
        }

        // Pobierz tokeny, usuń białe znaki i cudzysłowy
        TString name = TString(((TObjString *)tokens->At(0))->String()).Strip(TString::kBoth);
        name.ReplaceAll("\"", "");
        TString nBinsXStr = TString(((TObjString *)tokens->At(1))->String()).Strip(TString::kBoth);
        TString nBinsYStr = TString(((TObjString *)tokens->At(2))->String()).Strip(TString::kBoth);
        TString xMinStr = TString(((TObjString *)tokens->At(3))->String()).Strip(TString::kBoth);
        TString xMaxStr = TString(((TObjString *)tokens->At(4))->String()).Strip(TString::kBoth);
        TString yMinStr = TString(((TObjString *)tokens->At(5))->String()).Strip(TString::kBoth);
        TString yMaxStr = TString(((TObjString *)tokens->At(6))->String()).Strip(TString::kBoth);
        TString title = TString(((TObjString *)tokens->At(7))->String()).Strip(TString::kBoth);
        title.ReplaceAll("\"", "");
        TString xLabel = TString(((TObjString *)tokens->At(8))->String()).Strip(TString::kBoth);
        xLabel.ReplaceAll("\"", "");
        TString yLabel = TString(((TObjString *)tokens->At(9))->String()).Strip(TString::kBoth);
        yLabel.ReplaceAll("\"", "");
        TString zLabel = "";
        if (tokens->GetEntries() > 10)
        {
          zLabel = TString(((TObjString *)tokens->At(10))->String()).Strip(TString::kBoth);
          zLabel.ReplaceAll("\"", "");
        }

        Int_t nBinsX = nBinsXStr.Atoi();
        Int_t nBinsY = nBinsYStr.Atoi();
        Double_t xMin = xMinStr.Atof();
        Double_t xMax = xMaxStr.Atof();
        Double_t yMin = yMinStr.Atof();
        Double_t yMax = yMaxStr.Atof();

        configs[name] = {nBinsX, nBinsY, xMin, xMax, yMin, yMax, title, xLabel, yLabel, zLabel};

        delete tokens;
      }

      file.close();
      return configs;
    }

    // Helper functions implementacja
    TH1F *CreateHist1D(const TString &varName, const TString &histName)
    {
      auto it = histConfigs1D.find(varName);
      if (it == histConfigs1D.end())
      {
        // std::cerr << "WARNING: No 1D config found for variable: " << varName << std::endl;
        return new TH1F(histName.IsNull() ? varName : histName, varName, 100, -10., 10.);
      }

      const auto &config = it->second;
      TString name = histName.IsNull() ? ("h_" + varName) : histName;
      TH1F *hist = new TH1F(name, config.title, config.nBins, config.xMin, config.xMax);
      hist->GetXaxis()->SetTitle(config.xLabel);
      hist->GetYaxis()->SetTitle(config.yLabel);
      return hist;
    }

    TH1F *CreateHist1D(const TString &varName, const TString &channName, const TString &histName)
    {
      auto it = histConfigs1D.find(varName);
      if (it == histConfigs1D.end())
      {
        // std::cerr << "WARNING: No 1D config found for variable: " << varName << std::endl;
        return new TH1F(histName.IsNull() ? varName : histName, varName, 100, -10., 10.);
      }

      const auto &config = it->second;
      TString name = histName.IsNull() ? ("h_" + varName + "_" + channName) : histName;
      TH1F *hist = new TH1F(name, config.title, config.nBins, config.xMin, config.xMax);
      hist->GetXaxis()->SetTitle(config.xLabel);
      hist->GetYaxis()->SetTitle(config.yLabel);
      return hist;
    }

    TH1F *CreateHist1D(const TString &histName, const HistConfig1D &config)
    {
      // Tworzenie histogramu na podstawie konfiguracji
      TH1F *hist = new TH1F(histName, config.title, config.nBins, config.xMin, config.xMax);
      hist->GetXaxis()->SetTitle(config.xLabel);
      hist->GetYaxis()->SetTitle(config.yLabel);
      return hist;
    }

    TH2F *CreateHist2D(const TString &var1, const TString &var2, const TString &histName)
    {
      auto config1 = histConfigs1D.find(var1);
      auto config2 = histConfigs1D.find(var2);

      if (config1 == histConfigs1D.end() || config2 == histConfigs1D.end())
      {
        std::cerr << "WARNING: Config not found for variables: " << var1 << ", " << var2 << std::endl;
        return new TH2F(histName.IsNull() ? (var1 + "_vs_" + var2) : histName,
                        var1 + " vs " + var2, 100, -10., 10., 100, -10., 10.);
      }

      const auto &c1 = config1->second;
      const auto &c2 = config2->second;
      TString name = histName.IsNull() ? ("h_" + var1 + "_vs_" + var2) : histName;
      TString title = c1.title + " vs " + c2.title;

      TH2F *hist = new TH2F(name, title, c1.nBins, c1.xMin, c1.xMax,
                            c2.nBins, c2.xMin, c2.xMax);
      hist->GetXaxis()->SetTitle(c1.xLabel);
      hist->GetYaxis()->SetTitle(c2.xLabel);
      hist->GetZaxis()->SetTitle("Events");
      return hist;
    }

    TH2F *CreateHist2D(const TString &configName, const TString &histName)
    {
      auto it = histConfigs2D.find(configName);
      if (it == histConfigs2D.end())
      {
        std::cerr << "WARNING: No 2D config found for: " << configName << std::endl;
        return new TH2F(histName.IsNull() ? configName : histName, configName,
                        100, -10., 10., 100, -10., 10.);
      }

      const auto &config = it->second;
      TString name = histName.IsNull() ? ("h_" + configName) : histName;
      TH2F *hist = new TH2F(name, config.title, config.nBinsX, config.xMin, config.xMax,
                            config.nBinsY, config.yMin, config.yMax);
      hist->GetXaxis()->SetTitle(config.xLabel);
      hist->GetYaxis()->SetTitle(config.yLabel);
      hist->GetZaxis()->SetTitle(config.zLabel);
      return hist;
    }

    TH2F *CreateHist2D(const TString &histName, const HistConfig2D &config)
    {
      // Tworzenie histogramu na podstawie konfiguracji
      TH2F *hist = new TH2F(histName, config.title, config.nBinsX, config.xMin, config.xMax,
                            config.nBinsY, config.yMin, config.yMax);
      hist->GetXaxis()->SetTitle(config.xLabel);
      hist->GetYaxis()->SetTitle(config.yLabel);
      hist->GetZaxis()->SetTitle(config.zLabel);
      return hist;
    }
  }

  void setGlobalStyle()
  {
    // Global Style of histograms, pads, etc.

    gStyle->SetStatX(0.85); // prawy brzeg w NDC
    gStyle->SetStatY(0.90); // górny brzeg w NDC
    gStyle->SetStatW(0.23); // szerokość
    gStyle->SetStatH(0.18); // wysokość

    gStyle->SetOptStat("iouMn");

    gStyle->SetFitFormat("6.2g");
    gStyle->SetStatFormat("6.2g");

    gStyle->SetCanvasDefH(750);
    gStyle->SetCanvasDefW(750);

    gStyle->SetPadLeftMargin(0.15);
    gStyle->SetPadRightMargin(0.15);
    gStyle->SetPadBottomMargin(0.15);

    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);

    gStyle->SetOptLogz(1);
    gStyle->SetPalette(1);

    gStyle->SetHistLineWidth(3);
    gStyle->SetLineWidth(2);
    gStyle->SetFrameLineWidth(2);
    gStyle->SetTitleFont(62, "XYZ");

    gStyle->SetLabelSize(0.04, "X");
    gStyle->SetLabelSize(0.04, "Y");
    gStyle->SetLabelSize(0.04, "Z");
    gStyle->SetLabelOffset(0.02, "X");
    gStyle->SetLabelOffset(0.02, "Y");
    gStyle->SetLabelOffset(0.02, "Z");
    gStyle->SetLabelFont(62, "X");
    gStyle->SetLabelFont(62, "Y");
    gStyle->SetLabelFont(62, "Z");

    gStyle->SetTitleSize(0.05, "X");
    gStyle->SetTitleSize(0.05, "Y");
    gStyle->SetTitleSize(0.05, "Z");
    gStyle->SetTitleOffset(1.2, "X");
    gStyle->SetTitleOffset(1.5, "Y");
    gStyle->SetTitleOffset(1.2, "Z");
    gStyle->SetTitleFont(62, "X");
    gStyle->SetTitleFont(62, "Y");
    gStyle->SetTitleFont(62, "Z");

    gStyle->SetTitle("");

    TGaxis::SetMaxDigits(3);

    gStyle->cd();
  }
}

namespace Utils
{
  json properties;
  json constants;
  json paths;

  TString elapsedTimeHMS(double totalSeconds)
  {
    int elapsedMinutes, elapsedHours;
    double elapsedSeconds;
    TString elapsedHMS;

    elapsedHours = int(totalSeconds / 3600.);
    elapsedMinutes = int(((totalSeconds / 3600.) - elapsedHours) * 60.);
    elapsedSeconds = ((((totalSeconds / 3600.) - elapsedHours) * 60.) - elapsedMinutes) * 60.;

    elapsedHMS = std::to_string(elapsedHours) + "h " + std::to_string(elapsedMinutes) + "min " + std::to_string(elapsedSeconds) + "s";

    return elapsedHMS;
  }

  void InitializeVariables()
  {

    // Parsing of paths
    std::string pathsFilePath = Paths::propertiesPath + "/paths-extensions.json";
    std::ifstream fpath(pathsFilePath.c_str());

    if (fpath.is_open())
    {
      paths = json::parse(fpath);

      Paths::propName = (std::string)paths["analysisProperties"]["properties"];
      Paths::pdgConstFilePath = (std::string)paths["analysisProperties"]["PDGConst"];
      Paths::analysisConfigPath = (std::string)paths["analysisProperties"]["analysisConfig"];

      Paths::ext_img = (std::string)paths["extensions"]["img"];
      Paths::ext_root = (std::string)paths["extensions"]["root"];
      Paths::cutlimitsName = (std::string)paths["cutLimits"];

      std::cout << "Paths initialized from: " << pathsFilePath << std::endl;
    }

    // Parsing of constants from PDG JSON file
    std::ifstream fconst(Paths::pdgConstFilePath.c_str());

    auto &pdg = PdgManager::getInstance();
    std::map<TString, PDGProvider::SummaryEdition *> summaryMap;

    summaryMap["K0"] = &pdg.getParticleData("S011", 2025); // K0 summary data for 2025 PDG
    summaryMap["KS"] = &pdg.getParticleData("S012", 2025); // K0 summary data for 2025 PDG
    summaryMap["KL"] = &pdg.getParticleData("S013", 2025); // K0 summary data for 2025 PDG

    auto setConstant = [](boost::optional<std::vector<PDGProvider::Property>> &properties_vec, const TString &pdgid, Double_t &constant, CPTStatus CPTOrNotCPT = CPTStatus::UNDEFINED, Double_t multiplier = 1.)
    {
      for (const auto &prop : *properties_vec)
      {
        auto pdgid_opt = prop.get_pdgid();
        auto value_desc = prop.get_description();

        if (*pdgid_opt == pdgid)
        {
          constant = PdgManager::getBestValue(prop, 1.0) * multiplier;
          std::cout << "Set " << *value_desc << " to: " << constant << std::endl;
        }
      }
    };

    for (const auto &entry : summaryMap)
    {
      auto summaries_opt = entry.second->get_summaries();

      const auto &summaries = *summaries_opt;
      auto properties_opt = summaries.get_properties();

      const auto &properties_vec = *properties_opt;

      setConstant(properties_opt, "S011M/2025", PhysicsConstants::mK0);
      setConstant(properties_opt, "S012T/2025", PhysicsConstants::tau_S_nonCPT, CPTStatus::NON_CPT);
      setConstant(properties_opt, "S012T/2025", PhysicsConstants::tau_S_CPT, CPTStatus::CPT);
      setConstant(properties_opt, "S013T/2025", PhysicsConstants::tau_L, CPTStatus::NON_CPT, 1E9);  
      setConstant(properties_opt, "S013D/2025", PhysicsConstants::delta_mass_nonCPT, CPTStatus::NON_CPT);
      setConstant(properties_opt, "S013EP/2025", PhysicsConstants::mod_epsilon);
      setConstant(properties_opt, "S013EPS/2025", PhysicsConstants::Re);
      setConstant(properties_opt, "S013EPI/2025", PhysicsConstants::Im_CPT, CPTStatus::CPT);
      setConstant(properties_opt, "S013EPI/2025", PhysicsConstants::Im_nonCPT, CPTStatus::NON_CPT);
      setConstant(properties_opt, "S013F+-/2025", PhysicsConstants::phi_pm_nonCPT, CPTStatus::NON_CPT);
      setConstant(properties_opt, "S013FOO/2025", PhysicsConstants::phi_00_nonCPT, CPTStatus::NON_CPT);
    }

    // if (fconst.is_open())
    // {
    //   constants = json::parse(fconst);

    //   PhysicsConstants::mK0 = constants.at("values").at("/S011M").get<Double_t>();
    //   PhysicsConstants::tau_S_nonCPT = constants.at("values").at("/S012T").get<Double_t>() * 1E9;
    //   PhysicsConstants::tau_L = constants.at("values").at("/S013T").get<Double_t>() * 1E9;
    //   PhysicsConstants::delta_mass_nonCPT = constants.at("values").at("/S013D").get<Double_t>();
    //   PhysicsConstants::mod_epsilon = constants.at("values").at("/S013EP").get<Double_t>();
    //   PhysicsConstants::Re = constants.at("values").at("/S013EPS").get<Double_t>();
    //   PhysicsConstants::Im_nonCPT = constants.at("values").at("/S013EPI").get<Double_t>() * (M_PI / 180.);
    //   PhysicsConstants::phi_pm_nonCPT = constants.at("values").at("/S013F+-").get<Double_t>();
    //   PhysicsConstants::phi_00_nonCPT = constants.at("values").at("/S013FOO").get<Double_t>();
    // }

    // Parsing of the KLOE properties
    std::ifstream fprop(Paths::propName.c_str());
    if (fprop.is_open())
    {
      properties = json::parse(fprop);

      Paths::path_tmp = (std::string)properties["variables"]["rootFiles"]["path"];
      Paths::path_cs = (std::string)properties["variables"]["rootFiles"]["pathControlSample"];

      KLOE::firstFileMax = properties["variables"]["rootFiles"]["firstFileMax"];
      KLOE::lastFileMax = properties["variables"]["rootFiles"]["lastFileMax"];
      KLOE::numOfThreads = properties["variables"]["parallelization"]["numOfThreads"];

      Filenames::gen_vars_tree = (std::string)properties["variables"]["tree"]["treename"]["mctruth"];
      Filenames::neutrec_triangle_tree = (std::string)properties["variables"]["tree"]["treename"]["trianglefinal"];
      Filenames::omegarec_tree = (std::string)properties["variables"]["tree"]["treename"]["omegarec"];
      Filenames::omegarec_kin_fit_tree = (std::string)properties["variables"]["tree"]["treename"]["omegarec"];

      Filenames::omegaRecPath = (std::string)properties["variables"]["tree"]["filename"]["omegarec"];
      Filenames::mctruthPath = (std::string)properties["variables"]["tree"]["filename"]["mctruth"];
      Filenames::genvarsPath = (std::string)properties["variables"]["tree"]["filename"]["generatedvars"];
      Filenames::trianglePath = (std::string)properties["variables"]["tree"]["filename"]["trianglefinal"];
    }
  }

}