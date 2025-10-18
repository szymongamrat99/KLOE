#include <const.h>

#include <stdlib.h>
#include <ctime>
#include <fstream>

using json = nlohmann::json;

namespace Paths
{
  const std::string kloedataPath = getenv("KLOE_DBV26_DK0");
  const std::string kloeMCPath = getenv("KLOE_DBV26_MK0");
  const std::string workdirPath = getenv("WORKDIR");
  const std::string chainDataFiles = kloedataPath + "/*.root";
  const std::string chainMCFiles = kloeMCPath + "/*.root";
  const std::string pdgConstFilePath = (std::string)getenv("PDGAPI") + "/pdg_const.json";
  const std::string propertiesPath = getenv("PROPERTIESKLOE");
  const std::string propName = propertiesPath + "/properties.json";
  const std::string rootfilesName = propertiesPath + "/root-files.json";
  const std::string cutlimitsName = propertiesPath + "/cut-limits.json";

  TString base_path = workdirPath + "/KLOE/";
  TString path_tmp = "";
  TString path_cs = "";
  TString prod2root_path_v26 = "/data/k2/DBV-26/DK0";
  TString ext_root = ".root";
  TString ext_img = ".png";
  TString ext_csv = ".csv";

  const TString gen_vars_dir = base_path + "Subanalysis/GeneratedVars/";
  const TString neutrec_dir = base_path + "Subanalysis/Neutrec/";
  const TString cpfit_dir = base_path + "Subanalysis/CPFit/";
  const TString covmatrix_dir = base_path + "Subanalysis/CovarianceMatrix/";
  const TString initialanalysis_dir = base_path + "Subanalysis/InitialAnalysis/";
  const TString omegarec_dir = base_path + "Subanalysis/OmegaRec/";
  const TString efficiency_dir = base_path + "Subanalysis/EfficiencyAnalysis/";
  const TString charged_dir = base_path + "Subanalysis/KchRec/";
  const TString plots_dir = base_path + "Subanalysis/Plots/";
  const TString root_files_dir = "root_files/";
  const TString input_dir = "input/";
  const TString logs_dir = "log/";
  const TString result_dir = "results/";
  const TString img_dir = "img/";
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
  Double_t mK0 = 497.611;       // MeV/c^2
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

  void setGlobalStyle()
  {
    // Global Style of histograms, pads, etc.

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

    gStyle->cd();
  }
}

namespace Utils
{
  json properties;
  json constants;

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
  };

  void InitializeVariables()
  {
    // Parsing of constants from PDG JSON file
    std::ifstream fconst(Paths::pdgConstFilePath.c_str());
    if (fconst.is_open())
    {
      constants = json::parse(fconst);

      // PhysicsConstants::mK0 = constants.at("values").at("/S011M").get<Double_t>();
      // PhysicsConstants::tau_S_nonCPT = constants.at("values").at("/S012T").get<Double_t>() * 1E9;
      // PhysicsConstants::tau_L = constants.at("values").at("/S013T").get<Double_t>() * 1E9;
      // PhysicsConstants::delta_mass_nonCPT = constants.at("values").at("/S013D").get<Double_t>();
      // PhysicsConstants::mod_epsilon = constants.at("values").at("/S013EP").get<Double_t>();
      // PhysicsConstants::Re = constants.at("values").at("/S013EPS").get<Double_t>();
      // PhysicsConstants::Im_nonCPT = constants.at("values").at("/S013EPI").get<Double_t>() * (M_PI / 180.);
      // PhysicsConstants::phi_pm_nonCPT = constants.at("values").at("/S013F+-").get<Double_t>();
      // PhysicsConstants::phi_00_nonCPT = constants.at("values").at("/S013FOO").get<Double_t>();
    }

    // Parsing of the KLOE properties
    std::ifstream fprop(Paths::propName.c_str());
    if (fprop.is_open())
    {
      properties = json::parse(fprop);

      // Paths::path_tmp = (std::string)properties["variables"]["rootFiles"]["path"];
      // Paths::path_cs = (std::string)properties["variables"]["rootFiles"]["pathControlSample"];

      // KLOE::firstFileMax = properties["variables"]["rootFiles"]["firstFileMax"];
      // KLOE::lastFileMax = properties["variables"]["rootFiles"]["lastFileMax"];
      // KLOE::numOfThreads = properties["variables"]["parallelization"]["numOfThreads"];

      // Filenames::gen_vars_tree = (std::string)properties["variables"]["tree"]["treename"]["mctruth"];
      // Filenames::neutrec_triangle_tree = (std::string)properties["variables"]["tree"]["treename"]["trianglefinal"];
      // Filenames::omegarec_tree = (std::string)properties["variables"]["tree"]["treename"]["omegarec"];
      // Filenames::omegarec_kin_fit_tree = (std::string)properties["variables"]["tree"]["treename"]["omegarec"];

      // Filenames::omegaRecPath = (std::string)properties["variables"]["tree"]["filename"]["omegarec"];
      // Filenames::mctruthPath = (std::string)properties["variables"]["tree"]["filename"]["mctruth"];
      // Filenames::genvarsPath = (std::string)properties["variables"]["tree"]["filename"]["generatedvars"];
      // Filenames::trianglePath = (std::string)properties["variables"]["tree"]["filename"]["trianglefinal"];
    }
  }

}