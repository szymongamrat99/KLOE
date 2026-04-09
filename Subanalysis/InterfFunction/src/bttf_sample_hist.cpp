#include <const.h>
#include <interf_function.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TF2.h>
#include <TCanvas.h>
#include <TChain.h>
#include <iostream>
#include <TKey.h>
#include <TROOT.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <boost/filesystem.hpp>
#include <boost/progress.hpp>
#include <TRandom3.h>

namespace fs = boost::filesystem;

std::vector<std::string> GetDirectories(const std::string &path)
{
  std::vector<std::string> directories;

  if (!fs::exists(path) || !fs::is_directory(path))
  {
    std::cerr << "Błąd: ścieżka nie istnieje lub nie jest folderem!" << std::endl;
    return directories;
  }

  // Iteruj po zawartości folderu
  fs::directory_iterator end_iter;
  for (fs::directory_iterator dir_itr(path); dir_itr != end_iter; ++dir_itr)
  {
    if (fs::is_directory(dir_itr->status()))
    {
      directories.push_back(dir_itr->path().filename().string());
    }
  }

  return directories;
}

int SelectDirectory(const std::vector<std::string> &directories)
{
  if (directories.empty())
  {
    std::cerr << "Brak folderów do wyboru!" << std::endl;
    return -1;
  }

  // Wyświetl listę
  std::cout << "\nDostępne foldery:\n"
            << std::endl;
  for (size_t i = 0; i < directories.size(); ++i)
  {
    std::cout << "  [" << i << "] " << directories[i] << std::endl;
  }

  // Pobierz wybór użytkownika
  std::cout << "\nWybierz folder (podaj numer): ";
  int choice;
  std::cin >> choice;

  if (choice < 0 || choice >= static_cast<int>(directories.size()))
  {
    std::cerr << "Błąd: nieprawidłowy wybór!" << std::endl;
    return -1;
  }

  return choice;
}

void ListHistogramNames(TFile *file)
{
  if (!file || file->IsZombie())
  {
    std::cerr << "Błąd: nieprawidłowy plik ROOT!" << std::endl;
    return;
  }

  std::cout << "\n=== Histogramy w pliku ===" << std::endl;

  TIter next(file->GetListOfKeys());
  TKey *key;

  while ((key = (TKey *)next()))
  {
    TString className = key->GetClassName();

    // Wyświetl tylko histogramy
    if (className.Contains("TH1") || className.Contains("TH2"))
    {
      std::cout << "  " << key->GetName() << std::endl;
    }
  }
}

int main()
{
  ROOT::EnableImplicitMT(8);

  ErrorHandling::ErrorLogs logger("log/");

  // // Ścieżka do folderu głównego
  // std::string basePath = "/data/ssd/gamrat/KLOE/Subanalysis/InterfFunction/img/theoretical_plots";

  // // Pobierz listę folderów
  // std::vector<std::string> directories = GetDirectories(basePath);

  // // Pozwól użytkownikowi wybrać
  // int selectedIndex = SelectDirectory(directories);

  // if (selectedIndex == -1)
  // {
  //   return 1;
  // }

  // // Zbuduj pełną ścieżkę do wybranego folderu
  // std::string selectedPath = (fs::path(basePath) / directories[selectedIndex]).string();
  // std::cout << "\nWybrany folder: " << selectedPath << std::endl;

  // Teraz możesz użyć selectedPath do otwarcia pliku ROOT
  TString rootFilePathMain = "/data/ssd/gamrat/KLOE/Subanalysis/InterfFunction/img/theoretical_plots/integralLimit_300_samples_10000000000_Re_0.00166_Im_-0.00198_tmp/*.root",
          rootFilePathAux = "/data/ssd/gamrat/KLOE/Subanalysis/InterfFunction/root_files/sampling/DoubleExponential/t1Min_0_t1Max_20_t2Min_0_t2Max_20_nSamples_10000000000_customParamsFlag_false_ReParam_0.00166_ImParam_-0.00198/*.root";

  TChain *chain_00pm = new TChain("00pm");
  TChain *chain_pmpm = new TChain("pmpm");

  chain_00pm->Add(rootFilePathMain);
  chain_pmpm->Add(rootFilePathMain);

  chain_00pm->Add(rootFilePathAux);
  chain_pmpm->Add(rootFilePathAux);

  KLOE::setGlobalStyle();
  Utils::InitializeVariables(logger);

  std::cout << "Do You want to set custom parameters? (1 - yes, 0 - no): ";
  int customParams = 0;
  std::cin >> customParams;

  Double_t reParam = PhysicsConstants::Re;
  Double_t imParam = PhysicsConstants::Im_nonCPT;

  if (customParams)
  {
    std::cout << "Set Re(epsilon'/epsilon) (default " << PhysicsConstants::Re << "): ";
    std::cin >> reParam;

    std::cout << "Set Im(epsilon'/epsilon) (default " << PhysicsConstants::Im_nonCPT << "): ";
    std::cin >> imParam;
  }

  std::cout << "Do You want to choose the number of events? (1 - yes, 0 - no): ";
  int customEvents = 0;
  std::cin >> customEvents;

  Long64_t nEvents = 0;

  if (customEvents)
  {
    std::cout << "Set number of events to process: ";
    std::cin >> nEvents;
  }
  else
  {
    nEvents = chain_00pm->GetEntries();
    std::cout << "Chosen default number of events: " << nEvents << std::endl;
  }

  TString fileMethodAddition[3] = {"exp_not_corrected", "exp_corrected", "uniform_weighted"};

  std::cout << "Choose method: " << std::endl;
  std::cout << "  [0] Use exponentially generated times (no initial weight correction)." << std::endl;
  std::cout << "  [1] Use exponentially generated times (corrected weight)." << std::endl;
  std::cout << "  [2] Use unifomly generated times weighted with interference." << std::endl;
  int methodChoice = 0;
  std::cin >> methodChoice;

  if (methodChoice < 0 || methodChoice > 2)
  {
    std::cerr << "Błąd: nieprawidłowy wybór metody!" << std::endl;
    return 1;
  }

  TTreeReader reader_00pm(chain_00pm);
  TTreeReaderValue<Double_t> tch(reader_00pm, "tch");
  TTreeReaderValue<Double_t> tne(reader_00pm, "tne");

  TTreeReader reader_pmpm(chain_pmpm);
  TTreeReaderValue<Double_t> t1(reader_pmpm, "t1");
  TTreeReaderValue<Double_t> t2(reader_pmpm, "t2");

  Double_t t1Min = 0.0, t1Max = 300.0, t2Min = 0.0, t2Max = 300.0;
  Double_t tchMin = 0.0, tchMax = 300.0, tneMin = 0.0, tneMax = 300.0;
  Double_t sigmat = 1.0; // Wartość sigma_t do ustawienia zakresów
  Int_t nBins = (t1Max - t1Min) / sigmat;

  TF2 *func_00pm = new TF2("I(#pi^{0}#pi^{0},t_{1},#pi^{+}#pi^{-},t_{2});t_{1} [#tau_{S}]; t_{2} [#tau_{S}]", &interf_function_00pm, t1Min, t1Max, t2Min, t2Max, 2);
  func_00pm->SetParameters(reParam, imParam);

  TF2 *func_pm00 = new TF2("I(#pi^{0}#pi^{0},t_{1},#pi^{+}#pi^{-},t_{2});t_{1} [#tau_{S}]; t_{2} [#tau_{S}]", &interf_function_pm00, t1Min, t1Max, t2Min, t2Max, 2);
  func_pm00->SetParameters(reParam, imParam);

  TF2 *func_pmpm = new TF2("I(#pi^{+}#pi^{-},t_{1},#pi^{+}#pi^{-},t_{2});t_{1} [#tau_{S}]; t_{2} [#tau_{S}]", &interf_function_pmpm, t1Min, t1Max, t2Min, t2Max, 2);
  func_pmpm->SetParameters(reParam, imParam);

  TH2 *h_00pm = new TH2D("h_00pm", "(00,+-);t_{00} [#tau_{S}];t_{+-} [#tau_{S}]", nBins, tchMin, tchMax, nBins, tneMin, tneMax);
  TH2 *h_pm00 = new TH2D("h_pm00", "(+-,00);t_{+-} [#tau_{S}];t_{00} [#tau_{S}]", nBins, tneMin, tneMax, nBins, tchMin, tchMax);
  TH2 *h_pmpm = new TH2D("h_pmpm", "(+-,+-);t_{1} [#tau_{S}];t_{2} [#tau_{S}]", nBins, t1Min, t1Max, nBins, t2Min, t2Max);

  TH2 *h_00pm_not_weighted = new TH2D("h_00pm_not_weighted", "(00,+-);t_{00} [#tau_{S}];t_{+-} [#tau_{S}]", nBins, tchMin, tchMax, nBins, tneMin, tneMax);
  TH2 *h_pm00_not_weighted = new TH2D("h_pm00_not_weighted", "(+-,00);t_{+-} [#tau_{S}];t_{00} [#tau_{S}]", nBins, tneMin, tneMax, nBins, tchMin, tchMax);
  TH2 *h_pmpm_not_weighted = new TH2D("h_pmpm_not_weighted", "(+-,+-);t_{1} [#tau_{S}];t_{2} [#tau_{S}]", nBins, t1Min, t1Max, nBins, t2Min, t2Max);

  TH1 *h_00_not_weighted = new TH1D("h_00_not_weighted", "t_{00}; t_{00} [#tau_{S}];Counts", nBins, tchMin, tchMax);
  TH1 *h_pm_not_weighted = new TH1D("h_pm_not_weighted", "t_{+-}; t_{+-} [#tau_{S}];Counts", nBins, tneMin, tneMax);
  TH1 *h_pmpm1_not_weighted = new TH1D("h_pmpm1_not_weighted", "t_{1}; t_{1} [#tau_{S}];Counts", nBins, t1Min, t1Max);
  TH1 *h_pmpm2_not_weighted = new TH1D("h_pmpm2_not_weighted", "t_{2}; t_{2} [#tau_{S}];;Counts", nBins, t2Min, t2Max);

  const Long64_t evTick = 1000000;
  Long64_t currentEvent = 0;

  Double_t inverseWeight = PhysicsConstants::br_ks_pi0pi0 * PhysicsConstants::br_ks_pippim / (PhysicsConstants::br_ks_pippim * PhysicsConstants::br_ks_pippim);

  Double_t gammaS = 1.0; // Wartość gamma_S do ustawienia zakresów
  Double_t gammaL = PhysicsConstants::tau_S_nonCPT / PhysicsConstants::tau_L;

  Double_t weightA = PhysicsConstants::br_ks_pi0pi0 * PhysicsConstants::br_kl_pippim,
           weightB = PhysicsConstants::br_ks_pippim * PhysicsConstants::br_kl_pi0pi0,
           weightC = PhysicsConstants::br_ks_pippim * PhysicsConstants::br_kl_pippim,
           totalWeight = weightA + weightB,
           probA = weightA / totalWeight,
           normFactS = (1 - exp(-gammaS * 300.0)),
           normFactL = (1 - exp(-gammaL * 300.0));

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

  TRandom3 randGen(0); // Inicjalizacja generatora losowego z losowym seedem

  while (reader_00pm.Next())
  {
    // if (currentEvent >= nEvents)
    //   break;

    Double_t weight00pm;
    Double_t weightpm00;

    Double_t tchValue = *tch, tneValue = *tne;

    if (methodChoice == 0)
    {
      weight00pm = func_00pm->Eval(tneValue, tchValue);
      weightpm00 = func_pm00->Eval(tchValue, tneValue);
    }
    else if (methodChoice == 1)
    {
      weight00pm = func_00pm->Eval(tneValue, tchValue) / (double_exp_mixed_corr_00pm(tneValue, tchValue));
      weightpm00 = func_pm00->Eval(tchValue, tneValue) / (double_exp_mixed_corr_pm00(tchValue, tneValue));
    }
    else // methodChoice == 2
    {
      tchValue = randGen.Uniform(0.0, 300.0); // Generowanie losowej wartości czasu t1
      tneValue = randGen.Uniform(0.0, 300.0); // Generowanie losowej wartości czasu t2

      weight00pm = func_00pm->Eval(tneValue, tchValue);
      weightpm00 = func_pm00->Eval(tchValue, tneValue);
    }

    h_00pm->Fill(tneValue, tchValue, weight00pm);
    h_pm00->Fill(tchValue, tneValue, weightpm00);

    h_00pm_not_weighted->Fill(tneValue, tchValue);
    h_pm00_not_weighted->Fill(tchValue, tneValue);

    h_00_not_weighted->Fill(tneValue);
    h_pm_not_weighted->Fill(tchValue);

    if (currentEvent % evTick == 0)
      std::cout << "Processing entry: " << currentEvent << "\r";

    currentEvent++;
  }

  currentEvent = 0;

  while (reader_pmpm.Next())
  {
    // if (currentEvent >= nEvents)
    //   break;

    Double_t weight = 1.0;

    Double_t t1Value = *t1, t2Value = *t2;

    if (methodChoice == 0)
    {
      weight = func_pmpm->Eval(t1Value, t2Value);
    }
    else if (methodChoice == 1)
    {
      weight = inverseWeight * func_pmpm->Eval(t1Value, t2Value) / (double_exp(t1Value, t2Value));
    }
    else // methodChoice == 2
    {
      t1Value = randGen.Uniform(0.0, 300.0); // Generowanie losowej wartości czasu t1
      t2Value = randGen.Uniform(0.0, 300.0); // Generowanie losowej wartości czasu t2

      weight = func_pmpm->Eval(t1Value, t2Value);
    }

    h_pmpm->Fill(t1Value, t2Value, weight);
    h_pmpm_not_weighted->Fill(t1Value, t2Value);

    h_pmpm1_not_weighted->Fill(t1Value);
    h_pmpm2_not_weighted->Fill(t2Value);

    if (currentEvent % evTick == 0)
      std::cout << "Processing entry: " << currentEvent << "\r";

    currentEvent++;
  }

  // h_00pm->Scale(1.0 / h_00pm->GetEntries());
  // h_pm00->Scale(1.0 / h_pm00->GetEntries());
  // h_pmpm->Scale(1.0 / h_pmpm->GetEntries());

  TString rootFolderPath = "root_files/" + fileMethodAddition[methodChoice];

  if (!fs::exists((std::string)rootFolderPath))
  {
    fs::create_directories((std::string)rootFolderPath);
    std::cout << "Created directory: " << rootFolderPath << std::endl;
  }

  TString rootFileName = Form("%s/histograms2D_%lld_%.5f_%.5f.root", rootFolderPath.Data(), nEvents, reParam, imParam);

  TFile *outFile = new TFile(rootFileName, "RECREATE");

  outFile->cd();
  h_00pm->Write();
  h_pm00->Write();
  h_pmpm->Write();
  h_00pm_not_weighted->Write();
  h_pm00_not_weighted->Write();
  h_pmpm_not_weighted->Write();
  h_00_not_weighted->Write();
  h_pm_not_weighted->Write();
  h_pmpm1_not_weighted->Write();
  h_pmpm2_not_weighted->Write();
  outFile->Close();

  return 0;
}