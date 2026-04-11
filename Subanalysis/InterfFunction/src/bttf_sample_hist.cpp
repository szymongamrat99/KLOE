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

#include <sstream>

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

int main(Int_t argc, char *argv[])
{
  ErrorHandling::ErrorLogs logger("log/");
  KLOE::setGlobalStyle();
  Utils::InitializeVariables(logger);

  if (argc < 2)
  {
    std::cerr << "Usage: " << argv[0]
              << " input.root [method] [t_res] [re] [im] [t1Min] [t1Max] [t2Min] [t2Max] [tchMin] [tchMax] [tneMin] [tneMax]" << std::endl;
    return 1;
  }

  TString inputFilePath = argv[1];
  Long64_t customEvents = 0;
  Double_t reParam = PhysicsConstants::Re,
           imParam = PhysicsConstants::Im_nonCPT,
           sigmat = 1.0,
           t1Min = 0.0, t1Max = 300.0,
           t2Min = 0.0, t2Max = 300.0,
           tchMin = 0.0, tchMax = 300.0,
           tneMin = 0.0, tneMax = 300.0;

  Int_t chosenMethod = 1;

  if (argc > 2)
  {
    chosenMethod = std::atoi(argv[2]);

    if (chosenMethod < 0 || chosenMethod > 2)
    {
      std::cerr << "Error: Improper method!" << std::endl;
      std::cout << "Choose method: " << std::endl;
      std::cout << "  [0] Use exponentially generated times (no initial weight correction)." << std::endl;
      std::cout << "  [1] Use exponentially generated times (corrected weight)." << std::endl;
      std::cout << "  [2] Use unifomly generated times weighted with interference." << std::endl;
      return 1;
    }
    std::cout << "Chosen method: " << chosenMethod << std::endl;
  }

  if (argc > 3)
  {
    sigmat = std::stod(argv[3]);
    if (sigmat <= 0)
    {
      std::cerr << "Error: sigma_t must be positive!" << std::endl;
      return 1;
    }

    std::cout << "Chosen sigma_t: " << sigmat << std::endl;
  }

  if (argc > 4)
  {
    reParam = std::stod(argv[4]);
    std::cout << "Chosen Re(epsilon'/epsilon): " << reParam << std::endl;
  }

  if (argc > 5)
  {
    imParam = std::stod(argv[5]);
    std::cout << "Chosen Im(epsilon'/epsilon): " << imParam << std::endl;
  }

  if (argc > 6)
  {
    sigmat = std::stod(argv[6]);
    if (sigmat <= 0)
    {
      std::cerr << "Error: sigma_t must be positive!" << std::endl;
      return 1;
    }

    std::cout << "Chosen sigma_t: " << sigmat << std::endl;
  }

  if (argc > 7)
  {
    t1Min = std::stod(argv[7]);
    if (t1Min < 0)
    {
      std::cerr << "Error: t1Min must be non-negative!" << std::endl;
      return 1;
    }
    std::cout << "Chosen t1Min: " << t1Min << std::endl;
  }

  if (argc > 8)
  {
    t1Max = std::stod(argv[8]);
    if (t1Max < 0)
    {
      std::cerr << "Error: t1Max must be non-negative!" << std::endl;
      return 1;
    }
    std::cout << "Chosen t1Max: " << t1Max << std::endl;
  }

  if (argc > 9)
  {
    t2Min = std::stod(argv[9]);
    if (t2Min < 0)
    {
      std::cerr << "Error: t2Min must be non-negative!" << std::endl;
      return 1;
    }
    std::cout << "Chosen t2Min: " << t2Min << std::endl;
  }

  if (argc > 10)
  {
    t2Max = std::stod(argv[10]);
    if (t2Max < 0)
    {
      std::cerr << "Error: t2Max must be non-negative!" << std::endl;
      return 1;
    }
    std::cout << "Chosen t2Max: " << t2Max << std::endl;
  }

  if (argc > 11)
  {
    tchMin = std::stod(argv[11]);
    if (tchMin < 0)
    {
      std::cerr << "Error: tchMin must be non-negative!" << std::endl;
      return 1;
    }

    std::cout << "Chosen tchMin: " << tchMin << std::endl;
  }

  if (argc > 12)
  {
    tchMax = std::stod(argv[12]);
    if (tchMax < 0)
    {
      std::cerr << "Error: tchMax must be non-negative!" << std::endl;
      return 1;
    }
    std::cout << "Chosen tchMax: " << tchMax << std::endl;
  }

  if (argc > 13)
  {
    tneMin = std::stod(argv[13]);
    if (tneMin < 0)
    {
      std::cerr << "Error: tneMin must be non-negative!" << std::endl;
      return 1;
    }
    std::cout << "Chosen tneMin: " << tneMin << std::endl;
  }

  if (argc > 14)
  {
    tneMax = std::stod(argv[14]);
    if (tneMax < 0)
    {
      std::cerr << "Error: tneMax must be non-negative!" << std::endl;
      return 1;
    }
    std::cout << "Chosen tneMax: " << tneMax << std::endl;
  }

  TChain *chain_00pm = new TChain("00pm");
  TChain *chain_pmpm = new TChain("pmpm");

  chain_00pm->Add(inputFilePath);
  chain_pmpm->Add(inputFilePath);

  Long64_t nEntries_00pm = chain_00pm->GetEntries(),
           nEntries_pmpm = chain_pmpm->GetEntries();

  TString fileMethodAddition[3] = {"exp_not_corrected", "exp_corrected", "uniform_weighted"};

  // Creation of file
  TString rootFolderPath = "root_files/filled_sampled_histograms/" + fileMethodAddition[chosenMethod] + Form("/re%.5f_im%.5f_sigmat%.2f_t1%.2f_%.2f_t2%.2f_%.2f_tch%.2f_%.2f_tne%.2f_%.2f",
                         reParam, imParam, sigmat, t1Min, t1Max, t2Min, t2Max, tchMin, tchMax, tneMin, tneMax);

  if (!fs::exists((std::string)rootFolderPath))
  {
    fs::create_directories((std::string)rootFolderPath);
    std::cout << "Created directory: " << rootFolderPath << std::endl;
  }

  fs::path inputPath(argv[1]);
  std::string key = inputPath.string();
  std::size_t h = std::hash<std::string>{}(key);
  std::ostringstream oss;
  oss << std::hex << h;
  std::string hashStr = oss.str();

  TString outName = Form("%s/histograms2D_%lld_%lld_%s.root",
                         rootFolderPath.Data(), nEntries_00pm, nEntries_pmpm, hashStr.c_str());

  TFile *outFile = new TFile(outName, "RECREATE");
  // ----

  TTreeReader reader_00pm(chain_00pm);
  TTreeReaderValue<Double_t> tch(reader_00pm, "tch");
  TTreeReaderValue<Double_t> tne(reader_00pm, "tne");

  TTreeReader reader_pmpm(chain_pmpm);
  TTreeReaderValue<Double_t> t1(reader_pmpm, "t1");
  TTreeReaderValue<Double_t> t2(reader_pmpm, "t2");

  Int_t nBinst1 = (t1Max - t1Min) / sigmat,
        nBinst2 = (t2Max - t2Min) / sigmat,
        nBinsch = (tchMax - tchMin) / sigmat,
        nBinsne = (tneMax - tneMin) / sigmat;

  TF2 *func_00pm = new TF2("I(#pi^{0}#pi^{0},t_{1},#pi^{+}#pi^{-},t_{2});t_{1} [#tau_{S}]; t_{2} [#tau_{S}]", &interf_function_00pm, t1Min, t1Max, t2Min, t2Max, 2);
  func_00pm->SetParameters(reParam, imParam);

  TF2 *func_pm00 = new TF2("I(#pi^{0}#pi^{0},t_{1},#pi^{+}#pi^{-},t_{2});t_{1} [#tau_{S}]; t_{2} [#tau_{S}]", &interf_function_pm00, t1Min, t1Max, t2Min, t2Max, 2);
  func_pm00->SetParameters(reParam, imParam);

  TF2 *func_pmpm = new TF2("I(#pi^{+}#pi^{-},t_{1},#pi^{+}#pi^{-},t_{2});t_{1} [#tau_{S}]; t_{2} [#tau_{S}]", &interf_function_pmpm, t1Min, t1Max, t2Min, t2Max, 2);
  func_pmpm->SetParameters(reParam, imParam);

  TH2 *h_00pm = new TH2D("h_00pm", "(00,+-);t_{00} [#tau_{S}];t_{+-} [#tau_{S}]", nBinsch, tchMin, tchMax, nBinsne, tneMin, tneMax);
  TH2 *h_pm00 = new TH2D("h_pm00", "(+-,00);t_{+-} [#tau_{S}];t_{00} [#tau_{S}]", nBinsne, tneMin, tneMax, nBinsch, tchMin, tchMax);
  TH2 *h_pmpm = new TH2D("h_pmpm", "(+-,+-);t_{1} [#tau_{S}];t_{2} [#tau_{S}]", nBinst1, t1Min, t1Max, nBinst2, t2Min, t2Max);

  TH2 *h_00pm_not_weighted = new TH2D("h_00pm_not_weighted", "(00,+-);t_{00} [#tau_{S}];t_{+-} [#tau_{S}]", nBinsch, tchMin, tchMax, nBinsne, tneMin, tneMax);
  TH2 *h_pm00_not_weighted = new TH2D("h_pm00_not_weighted", "(+-,00);t_{+-} [#tau_{S}];t_{00} [#tau_{S}]", nBinsne, tneMin, tneMax, nBinsch, tchMin, tchMax);
  TH2 *h_pmpm_not_weighted = new TH2D("h_pmpm_not_weighted", "(+-,+-);t_{1} [#tau_{S}];t_{2} [#tau_{S}]", nBinst1, t1Min, t1Max, nBinst2, t2Min, t2Max);

  TH1 *h_00_not_weighted = new TH1D("h_00_not_weighted", "t_{00}; t_{00} [#tau_{S}];Counts", nBinsch, tchMin, tchMax);
  TH1 *h_pm_not_weighted = new TH1D("h_pm_not_weighted", "t_{+-}; t_{+-} [#tau_{S}];Counts", nBinsne, tneMin, tneMax);
  TH1 *h_pmpm1_not_weighted = new TH1D("h_pmpm1_not_weighted", "t_{1}; t_{1} [#tau_{S}];Counts", nBinst1, t1Min, t1Max);
  TH1 *h_pmpm2_not_weighted = new TH1D("h_pmpm2_not_weighted", "t_{2}; t_{2} [#tau_{S}];;Counts", nBinst2, t2Min, t2Max);

  const Long64_t evTick = 10000000; // co ile zdarzeń wypisywać postęp
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
    Double_t weight00pm;
    Double_t weightpm00;

    Double_t tchValue = *tch, tneValue = *tne;

    if (chosenMethod == 0)
    {
      weight00pm = func_00pm->Eval(tneValue, tchValue);
      weightpm00 = func_pm00->Eval(tchValue, tneValue);
    }
    else if (chosenMethod == 1)
    {
      weight00pm = func_00pm->Eval(tneValue, tchValue) / (double_exp_mixed_corr_00pm(tneValue, tchValue));
      weightpm00 = func_pm00->Eval(tchValue, tneValue) / (double_exp_mixed_corr_pm00(tchValue, tneValue));
    }
    else // chosenMethod == 2
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
    Double_t weight = 1.0;

    Double_t t1Value = *t1, t2Value = *t2;

    if (chosenMethod == 0)
    {
      weight = func_pmpm->Eval(t1Value, t2Value);
    }
    else if (chosenMethod == 1)
    {
      weight = inverseWeight * func_pmpm->Eval(t1Value, t2Value) / (double_exp(t1Value, t2Value));
    }
    else // chosenMethod == 2
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