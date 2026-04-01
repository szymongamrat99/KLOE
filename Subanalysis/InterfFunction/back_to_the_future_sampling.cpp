#include <TCanvas.h>
#include <TF1.h>
#include <TF2.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TStyle.h>
#include <TAxis.h>
#include <TMath.h>
#include <TLegend.h>
#include <TFile.h>
#include <TDirectory.h>
#include <TRandom3.h>
#include <TBranch.h>
#include <TTree.h>
#include <nlohmann/json.hpp>

#include <boost/filesystem.hpp>
#include <boost/progress.hpp>

#include <interf_function.h>
#include <const.h>

namespace fs = boost::filesystem;

Double_t evTick = 10000.0;

void createDirIfNotExists(const TString &path)
{
  boost::filesystem::path dirPath(path.Data());

  if (!boost::filesystem::exists(dirPath))
    boost::filesystem::create_directories(dirPath);
}

void renameFile(std::string directory = "./", Int_t jobNum = 0)
{
  std::vector<fs::path> rootFiles;
  std::string pattern = std::to_string(jobNum) + "_tmp";

  // 1. Zbieramy wskaźniki (ścieżki) do plików spełniających warunek
  for (const auto &entry : fs::directory_iterator(directory))
  {
    std::string filename = entry.path().filename().string();

    // Szukamy plików name_tmp_1.root, name_tmp_2.root itd.
    if (filename.find(pattern) != std::string::npos && entry.path().extension() == ".root")
    {
      rootFiles.push_back(entry.path());
    }
  }

  // Opcjonalne sortowanie (żeby pliki szły po kolei _1, _2...)
  std::sort(rootFiles.begin(), rootFiles.end());

  // 2. Przetwarzamy zebraną listę
  for (const auto &oldPath : rootFiles)
  {
    std::string newName = oldPath.filename().string();

    size_t pos = newName.find("_tmp");
    if (pos != std::string::npos)
    {
      newName.erase(pos, 4); // usuwa "_tmp"
    }

    fs::path newPath = oldPath.parent_path() / newName;

    try
    {
      fs::rename(oldPath, newPath);
      std::cout << "Zmieniono: " << oldPath.filename() << " -> " << newName << std::endl;
    }
    catch (const fs::filesystem_error &e)
    {
      std::cerr << "Błąd przy pliku " << oldPath << ": " << e.what() << std::endl;
    }
  }
}

void GenerateHitOrMissSample(TH2D *hist, TF2 *func, Long64_t nSamplesTarget, Double_t &t1Hit, Double_t &t2Hit, TTree *tree)
{
  boost::progress_display progress(nSamplesTarget / evTick);

  Double_t t1min, t1max, t2min, t2max;
  func->GetRange(t1min, t2min, t1max, t2max);
  Double_t maxFuncVal = func->GetMaximum();

  Long64_t nSamples = 0;
  TRandom3 randGen(0); // Random generator with seed 0

  while (nSamples < nSamplesTarget)
  {
    // Generate random (t1, t2) within histogram range
    Double_t t1 = randGen.Uniform(t1min, t1max);
    Double_t t2 = randGen.Uniform(t2min, t2max);

    // Evaluate function value at (t1, t2)
    Double_t funcValue = func->Eval(t1, t2);
    // Generate a random number for hit-or-miss
    Double_t randValue = randGen.Uniform(0, maxFuncVal);

    // Accept or reject the sample
    if (randValue <= funcValue)
    {
      hist->Fill(t1, t2);
      nSamples++;

      t1Hit = t1;
      t2Hit = t2;

      tree->Fill();

      if (nSamples % Int_t(evTick) == 0)
        ++progress;
    }
  }
}

double draw_from_mixed(double ts, double tl)
{
  if (gRandom->Uniform() < 0.5)
  {
    return gRandom->Exp(ts); // Losuj z KS
  }
  else
  {
    return gRandom->Exp(tl); // Losuj z KL
  }
}

void GenerateImportanceSample(TH2D *hist, TF2 *func, Long64_t nSamplesTarget)
{
  boost::progress_display progress(nSamplesTarget / evTick);

  Double_t t1min, t1max, t2min, t2max;
  func->GetRange(t1min, t2min, t1max, t2max);
  Double_t maxFuncVal = func->GetMaximum();

  Long64_t nSamples = 0;
  TRandom3 randGen(0); // Random generator with seed 0

  while (nSamples < nSamplesTarget)
  {
    // Generate random (t1, t2) within histogram range
    Double_t t1 = draw_from_mixed(1, PhysicsConstants::tau_L / PhysicsConstants::tau_S_nonCPT);
    Double_t t2 = draw_from_mixed(1, PhysicsConstants::tau_L / PhysicsConstants::tau_S_nonCPT);

    if (t1 < t1min || t1 > t1max || t2 < t2min || t2 > t2max)
      continue;

    auto pdf = [&](double t)
    {
      Double_t tau_L = PhysicsConstants::tau_L / PhysicsConstants::tau_S_nonCPT;

      return 0.5 * (exp(-t)) + 0.5 * (1.0 / tau_L * exp(-t / tau_L));
    };

    // Evaluate function value at (t1, t2)
    Double_t funcValue = func->Eval(t1, t2);
    Double_t weight = pdf(t1) * pdf(t2); // Jacobian for exponential sampling
    // Generate a random number for hit-or-miss
    Double_t M = 1000.0;

    Double_t randValue = randGen.Uniform(0, M);

    // Accept or reject the sample
    if (randValue <= funcValue / weight)
    {
      hist->Fill(t1, t2);
      nSamples++;

      if (nSamples % Int_t(evTick) == 0)
        ++progress;
    }
  }
}

void GenerateRandom2(TH2D *hist, TF2 *func, Long64_t nSamplesTarget)
{
  boost::progress_display progress(nSamplesTarget / evTick);

  Int_t granularity = 1000;

  func->SetNpx(granularity);
  func->SetNpy(granularity);

  Long64_t nSamples = 0;
  Double_t t1, t2;

  TRandom3 *randGen = new TRandom3(0); // Random generator with seed 0

  while (nSamples < nSamplesTarget)
  {
    // Generate random (t1, t2) within histogram range
    func->GetRandom2(t1, t2, randGen);
    hist->Fill(t1, t2); // Weight to account for equal probabilities
    nSamples++;

    if (nSamples % Int_t(evTick) == 0)
      ++progress;
  }
}

void GenerateRandom21D(TH1 *hist, TF2 *func, Long64_t nSamplesTarget)
{
  boost::progress_display progress(nSamplesTarget / evTick);

  Int_t granularity = 1000;

  func->SetNpx(granularity);
  func->SetNpy(granularity);

  Long64_t nSamples = 0;
  Double_t t1, t2;

  TRandom3 *randGen = new TRandom3(0); // Random generator with seed 0

  while (nSamples < nSamplesTarget)
  {
    // Generate random (t1, t2) within histogram range
    func->GetRandom2(t1, t2, randGen);
    hist->Fill(t1); // Weight to account for equal probabilities
    nSamples++;

    if (nSamples % Int_t(evTick) == 0)
      ++progress;
  }
}

static TRandom3 gRand(0);

void generate_t1_t2(double &t1, double &t2, bool equalProb, double t1Min, double t1Max, double t2Min, double t2Max)
{
  Double_t weightA = PhysicsConstants::br_ks_pi0pi0 * PhysicsConstants::br_kl_pippim,
           weightB = PhysicsConstants::br_ks_pippim * PhysicsConstants::br_kl_pi0pi0,
           totalWeight = weightA + weightB,
           probA = weightA / totalWeight;

  if (equalProb)
    probA = 0.5;

  Double_t tauS = 1.0; // Wartość tau_S do ustawienia zakresów
  Double_t tauL = PhysicsConstants::tau_L / PhysicsConstants::tau_S_nonCPT;

  Double_t u1MinS = 1 - exp(-t1Min / tauS);
  Double_t u1MaxS = 1 - exp(-t1Max / tauS);
  Double_t u2MinS = 1 - exp(-t2Min / tauS);
  Double_t u2MaxS = 1 - exp(-t2Max / tauS);

  Double_t u1MinL = 1 - exp(-t1Min / tauL);
  Double_t u1MaxL = 1 - exp(-t1Max / tauL);
  Double_t u2MinL = 1 - exp(-t2Min / tauL);
  Double_t u2MaxL = 1 - exp(-t2Max / tauL);

  Double_t u = gRand.Uniform();
  Double_t u1 = 0.0, u2 = 0.0;

  if (u < probA)
  {
    u1 = gRand.Uniform(u1MinS, u1MaxS);
    u2 = gRand.Uniform(u2MinL, u2MaxL);

    t1 = -TMath::Log(1 - u1);
    t2 = -PhysicsConstants::tau_L / PhysicsConstants::tau_S_nonCPT * TMath::Log(1 - u2);
  }
  else
  {
    u1 = gRand.Uniform(u1MinL, u1MaxL);
    u2 = gRand.Uniform(u2MinS, u2MaxS);

    t1 = -PhysicsConstants::tau_L / PhysicsConstants::tau_S_nonCPT * TMath::Log(1 - u1);
    t2 = -TMath::Log(1 - u2);
  }
}

int main(int argc, char *argv[])
{
  Int_t jobNum = 0;

  if (argc > 1 && std::atoi(argv[1]) >= 0)
    jobNum = std::atoi(argv[1]);

  KLOE::setGlobalStyle();
  ErrorHandling::ErrorLogs logger("log/");
  Utils::InitializeVariables(logger);

  TString configPath = "config/config.json";
  std::ifstream configFile(configPath);
  if (!configFile.is_open())
  {
    std::cerr << "Error opening config file: " << configPath << std::endl;
    return 1;
  }

  std::unordered_map<int, TString> methodNames = {
      {1, "HitOrMiss"},
      {2, "GetRandom2D"},
      {3, "DoubleExponential"}};

  // Parameters from config file with default values

  Double_t t1Min = 0.0, t1Max = 300.0, t2Min = 0.0, t2Max = 300.0;
  Int_t nSamplesTarget = 100000, samplingMethod = 3;
  Bool_t customParamsFlag = false;
  Double_t reParam = PhysicsConstants::Re, imParam = PhysicsConstants::Im_nonCPT;

  nlohmann::json config;
  config = nlohmann::json::parse(configFile);

  Utils::JsonFieldLookupDouble(config, "sampling/t1Min", t1Min, logger);
  Utils::JsonFieldLookupDouble(config, "sampling/t1Max", t1Max, logger);
  Utils::JsonFieldLookupDouble(config, "sampling/t2Min", t2Min, logger);
  Utils::JsonFieldLookupDouble(config, "sampling/t2Max", t2Max, logger);

  Utils::JsonFieldLookupInt(config, "sampling/nSamplesTarget", nSamplesTarget, logger);
  Utils::JsonFieldLookupInt(config, "sampling/samplingMethod", samplingMethod, logger);

  Utils::JsonFieldLookupBool(config, "sampling/customParameters/flag", customParamsFlag, logger);

  if (customParamsFlag)
  {
    Utils::JsonFieldLookupDouble(config, "sampling/customParameters/parameters/ReParam", reParam, logger);
    Utils::JsonFieldLookupDouble(config, "sampling/customParameters/parameters/ImParam", imParam, logger);
  }

  // -----------------------------------------------------------------------------------

  // Sampling functions to create 2D histograms
  const Float_t sigmaT = 1; // Time resolution in tau_S units
  const Int_t nBinst1 = (t1Max - t1Min) / sigmaT, nBinst2 = (t2Max - t2Min) / sigmaT;

  TF2 *func_00pm = new TF2("I(#pi^{0}#pi^{0},t_{1},#pi^{+}#pi^{-},t_{2});t_{1} [#tau_{S}]; t_{2} [#tau_{S}]", &interf_function_00pm, t1Min, t1Max, t2Min, t2Max, 2);
  func_00pm->SetParameters(reParam, imParam);

  TF2 *func_pm00 = new TF2("I(#pi^{+}#pi^{-},t_{1},#pi^{0}#pi^{0},t_{2});t_{1} [#tau_{S}]; t_{2} [#tau_{S}]", &interf_function_pm00, t1Min, t1Max, t2Min, t2Max, 2);
  func_pm00->SetParameters(reParam, imParam);

  TF2 *func_pmpm = new TF2("I(#pi^{+}#pi^{-},t_{1},#pi^{+}#pi^{-},t_{2});t_{1} [#tau_{S}]; t_{2} [#tau_{S}]", &interf_function_pmpm, t1Min, t1Max, t2Min, t2Max, 2);
  func_pmpm->SetParameters(reParam, imParam);

  Double_t maxIntegral = 300.0;

  Long64_t nSamples = 0;

  TH2D *hist_00pm = new TH2D("hist_00pm", "I(#pi^{0}#pi^{0},t_{1},#pi^{+}#pi^{-},t_{2});t_{1} [#tau_{S}];t_{2} [#tau_{S}]", nBinst1, t1Min, t1Max, nBinst2, t2Min, t2Max);
  TH2D *hist_pm00 = new TH2D("hist_pm00", "I(#pi^{+}#pi^{-},t_{1},#pi^{0}#pi^{0},t_{2});t_{1} [#tau_{S}];t_{2} [#tau_{S}]", nBinst1, t1Min, t1Max, nBinst2, t2Min, t2Max);
  TH2D *hist_pmpm = new TH2D("hist_pmpm", "I(#pi^{+}#pi^{-},t_{1},#pi^{+}#pi^{-},t_{2});t_{1} [#tau_{S}];t_{2} [#tau_{S}]", nBinst1, t1Min, t1Max, nBinst2, t2Min, t2Max);

  // Creation of the output ROOT folder if it does not exist
  TString
      timeOrientedDir = "t1Min_" + TString::Format("%.0f", t1Min) + "_t1Max_" + TString::Format("%.0f", t1Max) + "_t2Min_" + TString::Format("%.0f", t2Min) + "_t2Max_" + TString::Format("%.0f", t2Max) + "_nSamples_" + TString::Format("%lld", nSamplesTarget) + "_customParamsFlag_" + (customParamsFlag ? "true" : "false"),
      root_dir = Paths::workdirPath + "/Subanalysis/InterfFunction/root_files/sampling/" + methodNames[samplingMethod] + "/" + timeOrientedDir + "/";

  createDirIfNotExists(root_dir);

  TString rootFileName = "";

  if (argc > 1 && std::atoi(argv[1]) >= 0)
    rootFileName = root_dir + "sampling_results_job_" + TString::Format("%d", jobNum) + "_tmp" + "_1.root";
  else
    rootFileName = root_dir + "sampling_results_1.root";

  Double_t tne, tch, t1, t2;

  TFile *outputFile = new TFile(rootFileName, "RECREATE");
  outputFile->cd();

  TTree *tree_00pm = new TTree("00pm", "Tree with sampled events");
  TTree *tree_pmpm = new TTree("pmpm", "Tree with sampled events");

  TBranch *branch_t1_00pm = tree_00pm->Branch("tne", &tne, "tne/D");
  TBranch *branch_t2_00pm = tree_00pm->Branch("tch", &tch, "tch/D");

  TBranch *branch_t1_pmpm = tree_pmpm->Branch("t1", &t1, "t1/D");
  TBranch *branch_t2_pmpm = tree_pmpm->Branch("t2", &t2, "t2/D");

  switch (samplingMethod)
  {
  case 1:
  {
    GenerateHitOrMissSample(hist_00pm, func_00pm, nSamplesTarget, tne, tch, tree_00pm);
    GenerateHitOrMissSample(hist_pm00, func_pm00, nSamplesTarget, tch, tne, tree_00pm);
    GenerateHitOrMissSample(hist_pmpm, func_pmpm, nSamplesTarget, t1, t2, tree_pmpm);
    break;
  }
  case 2:
  {
    GenerateRandom2(hist_00pm, func_00pm, nSamplesTarget);
    GenerateRandom2(hist_pm00, func_pm00, nSamplesTarget);
    GenerateRandom2(hist_pmpm, func_pmpm, nSamplesTarget);

    break;
  }
  case 3:
  {

    Double_t inverseWeight = PhysicsConstants::br_ks_pi0pi0 * PhysicsConstants::br_ks_pippim / (PhysicsConstants::br_ks_pippim * PhysicsConstants::br_ks_pippim);

    Double_t gammaS = 1.0; // Wartość gamma_S do ustawienia zakresów
    Double_t gammaL = PhysicsConstants::tau_S_nonCPT / PhysicsConstants::tau_L;

    Double_t weightA = PhysicsConstants::br_ks_pi0pi0 * PhysicsConstants::br_kl_pippim,
             weightB = PhysicsConstants::br_ks_pippim * PhysicsConstants::br_kl_pi0pi0,
             weightC = PhysicsConstants::br_ks_pippim * PhysicsConstants::br_kl_pippim,
             totalWeight = weightA + weightB,
             probA = weightA / totalWeight,
             normFactS = (1 - exp(-gammaS * maxIntegral)),
             normFactL = (1 - exp(-gammaL * maxIntegral));

    auto double_exp = [&](Double_t t1, Double_t t2)
    {
      return (0.5 / (normFactS * normFactL)) * (gammaS * gammaL * exp(-gammaS * t1 - gammaL * t2) + gammaS * gammaL * exp(-gammaL * t1 - gammaS * t2));
    };

    auto double_exp_mixed_corr_00pm = [&](Double_t t1, Double_t t2)
    {
      // Każdy człon musi być podzielony przez swoją całkę w zakresie [0, Tmax]
      Double_t term1 = (gammaS * exp(-gammaS * t1) / normFactS) * (gammaL * exp(-gammaL * t2) / normFactL);
      Double_t term2 = (gammaL * exp(-gammaL * t1) / normFactL) * (gammaS * exp(-gammaS * t2) / normFactS);

      return probA * term1 + (1 - probA) * term2;
    };

    auto double_exp_mixed_corr_pm00 = [&](Double_t t1, Double_t t2)
    {
      // Każdy człon musi być podzielony przez swoją całkę w zakresie [0, Tmax]
      Double_t term1 = (gammaS * exp(-gammaS * t2) / normFactS) * (gammaL * exp(-gammaL * t1) / normFactL);
      Double_t term2 = (gammaL * exp(-gammaL * t2) / normFactL) * (gammaS * exp(-gammaS * t1) / normFactS);

      return probA * term1 + (1 - probA) * term2;
    };

    auto double_exp_mixed_00pm = [&](Double_t t1, Double_t t2)
    {
      // Każdy człon musi być podzielony przez swoją całkę w zakresie [0, Tmax]
      Double_t term1 = (gammaS * exp(-gammaS * t1)) * (gammaL * exp(-gammaL * t2));
      Double_t term2 = (gammaL * exp(-gammaL * t1)) * (gammaS * exp(-gammaS * t2));

      return probA * term1 + (1 - probA) * term2;
    };

    auto double_exp_mixed_pm00 = [&](Double_t t1, Double_t t2)
    {
      // Każdy człon musi być podzielony przez swoją całkę w zakresie [0, Tmax]
      Double_t term1 = (gammaS * exp(-gammaS * t2)) * (gammaL * exp(-gammaL * t1));
      Double_t term2 = (gammaL * exp(-gammaL * t2)) * (gammaS * exp(-gammaS * t1));

      return probA * term1 + (1 - probA) * term2;
    };

    boost::progress_display progress00pm(nSamplesTarget / evTick);

    Double_t weightSamples = (PhysicsConstants::br_ks_pippim * PhysicsConstants::br_ks_pippim) / (PhysicsConstants::br_ks_pi0pi0 * PhysicsConstants::br_ks_pippim);

    std::cout << "Weight for samples: " << weightSamples << std::endl;

    Long64_t nSamples = 0;
    while (nSamples < nSamplesTarget)
    {
      generate_t1_t2(tne, tch, false, t1Min, t1Max, t2Min, t2Max);

      if (tne <= maxIntegral && tch <= maxIntegral)
      {
        Double_t weight00pm = func_00pm->Eval(tne, tch) / (double_exp_mixed_corr_00pm(tne, tch));
        Double_t weightpm00 = func_pm00->Eval(tch, tne) / (double_exp_mixed_corr_pm00(tch, tne));

        hist_00pm->Fill(tne, tch, weight00pm);
        hist_pm00->Fill(tch, tne, weightpm00);

        tree_00pm->Fill();

        nSamples++;

        if (nSamples % Int_t(evTick) == 0)
          ++progress00pm;
      }
    }

    boost::progress_display progresspmpm(weightSamples * nSamplesTarget / evTick);
    nSamples = 0;

    while (nSamples < weightSamples * nSamplesTarget)
    {
      generate_t1_t2(t1, t2, true, t1Min, t1Max, t2Min, t2Max);

      if (t1 <= maxIntegral && t2 <= maxIntegral)
      {
        Double_t weight = func_pmpm->Eval(t1, t2) / (double_exp(t1, t2) * weightSamples); // Dodatkowy weightSamples, aby uwzględnić fakt, że generujemy więcej próbek

        hist_pmpm->Fill(t1, t2, weight);

        tree_pmpm->Fill();

        nSamples++;

        if (nSamples % Int_t(evTick) == 0)
          ++progresspmpm;
      }
    }
    break;
  }
  }

  tree_00pm->Write();
  tree_pmpm->Write();

  outputFile->Close();

  renameFile(root_dir.Data(), jobNum);
}