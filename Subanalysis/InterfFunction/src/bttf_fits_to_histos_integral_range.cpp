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
#include <TLegend.h>
#include <TPaveText.h>
#include <TLine.h>
#include <TPaveStats.h>

#include <TFitResult.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <cmath>

namespace fs = boost::filesystem;

// Przykład 1: Podstawowa iteracja po plikach
void listFilesBasic(const std::string &folderPath)
{
  fs::path dir(folderPath);

  // Guard clause - sprawdź czy folder istnieje
  if (!fs::exists(dir) || !fs::is_directory(dir))
  {
    std::cerr << "ERROR: Directory does not exist: " << folderPath << std::endl;
    return;
  }

  fs::directory_iterator end_iter;
  for (fs::directory_iterator it(dir); it != end_iter; ++it)
  {
    if (fs::is_regular_file(it->status()))
    {
      std::string filename = it->path().filename().string();
      std::cout << "File: " << filename << std::endl;
    }
  }
}

// Przykład 2: Zebranie plików do vectora do dalszej obróbki
std::vector<std::string> collectRootFiles(const std::string &folderPath)
{
  std::vector<std::string> rootFiles;
  fs::path dir(folderPath);

  if (!fs::exists(dir) || !fs::is_directory(dir))
    return rootFiles;

  fs::directory_iterator end_iter;
  for (fs::directory_iterator it(dir); it != end_iter; ++it)
  {
    if (fs::is_regular_file(it->status()))
    {
      fs::path filePath = it->path();

      // Filtrowanie po rozszerzeniu
      if (filePath.extension() == ".root")
      {
        rootFiles.push_back(filePath.string()); // Pełna ścieżka
        // lub: rootFiles.push_back(filePath.filename().string()); // Tylko nazwa
      }
    }
  }

  return rootFiles;
}

// Przykład 3: Parsowanie nazw plików (np. "histograms2D_1000000_0.00166_0.00000.root")
struct FileParams
{
  Long64_t nEvents;
  Double_t reParam;
  Double_t imParam;
  std::string fullPath;
};

std::map<Long64_t, TString> parseHistogramFiles(const std::string &folderPath)
{
  std::map<Long64_t, TString> fileParams;
  fs::path dir(folderPath);

  if (!fs::exists(dir) || !fs::is_directory(dir))
    return fileParams;

  fs::directory_iterator end_iter;
  for (fs::directory_iterator it(dir); it != end_iter; ++it)
  {
    if (fs::is_regular_file(it->status()))
    {
      std::string filename = it->path().filename().string();

      // Parsowanie nazwy pliku
      if (filename.find("histograms2D_") != std::string::npos)
      {
        FileParams params;
        params.fullPath = it->path().string();

        // Parsowanie parametrów z nazwy
        int parsed = sscanf(filename.c_str(), "histograms2D_%lld_%lf_%lf.root",
                            &params.nEvents, &params.reParam, &params.imParam);

        Double_t reParamTemplate = 0.00166;
        Double_t imParamTemplate = -0.00192;
        Double_t tolerance = 1e-5;

        if (parsed == 3)
        {
          if (abs(params.reParam - reParamTemplate) < tolerance &&
              abs(params.imParam - imParamTemplate) < tolerance)
          {
            fileParams[params.nEvents] = filename;
            std::cout << "Found: " << filename
                      << " (events=" << params.nEvents
                      << ", Re=" << params.reParam
                      << ", Im=" << params.imParam << ")" << std::endl;
          }
        }
      }
    }
  }

  return fileParams;
}

int main()
{
  std::cout << "What addition for the folder name do you want to use? (e.g. 'integral_range' or 'fit_range')" << std::endl;
  std::string folderAddition = "";
  std::getline(std::cin, folderAddition);

  TString imgFolderPath = "img/integral_range/" + folderAddition;

  if (folderAddition == "")
  {
    imgFolderPath = "img/integral_range";
  }

  if (!fs::exists((std::string)imgFolderPath))
  {
    fs::create_directories((std::string)imgFolderPath);
    std::cout << "Created directory: " << imgFolderPath << std::endl;
  }

  ROOT::EnableImplicitMT(8);

  auto formatExpTitle = [](Long64_t n) -> TString
  {
    if (n <= 0)
      return "0";
    int exp10 = static_cast<int>(std::round(std::log10(static_cast<double>(n))));
    return Form("10^{%d}", exp10);
  };

  auto formatExpFile = [](Long64_t n) -> TString
  {
    if (n <= 0)
      return "0";
    int exp10 = static_cast<int>(std::round(std::log10(static_cast<double>(n))));
    return Form("1e%d", exp10);
  };

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

  std::string rootFolder = "root_files/" + (std::string)fileMethodAddition[methodChoice];
  std::map<Long64_t, TString> files = parseHistogramFiles(rootFolder);

  // Guard clause - sprawdź czy są pliki
  if (files.empty())
  {
    std::cerr << "ERROR: No files found!" << std::endl;
    return 1;
  }

  // Weź plik z największą liczbą zliczeń
  auto maxEntry = *files.rbegin();
  Long64_t maxEvents = maxEntry.first;
  TString maxFileName = maxEntry.second;
  TString maxEventsTitle = formatExpTitle(maxEvents);
  TString maxEventsFile = formatExpFile(maxEvents);

  std::cout << "Processing file with max events: " << maxFileName
            << " (" << maxEvents << " events)" << std::endl;

  // Wczytaj histogramy 2D z pliku z max zliczeniami
  TString filePath = Form("%s/%s", rootFolder.c_str(), maxFileName.Data());
  TFile *file = TFile::Open(filePath);
  if (!file || file->IsZombie())
  {
    std::cerr << "Error opening file: " << filePath << std::endl;
    return 1;
  }

  TH2 *hist_00pm2D_base = nullptr;
  TH2 *hist_pm002D_base = nullptr;
  TH2 *hist_pmpm2D_base = nullptr;

  TH2 *hist_00pm2D_base_nw = nullptr;
  TH2 *hist_pm002D_base_nw = nullptr;
  TH2 *hist_pmpm2D_base_nw = nullptr;

  TH1 *hist_00 = nullptr;
  TH1 *hist_pm = nullptr;
  TH1 *hist_1 = nullptr;
  TH1 *hist_2 = nullptr;

  TH1 *deltaTpm00_not_weighted = nullptr;
  TH1 *deltaTpm00_weighted = nullptr;

  TH1 *deltaT12_not_weighted = nullptr;
  TH1 *deltaT12_weighted = nullptr;

  file->GetObject("h_00pm", hist_00pm2D_base);
  file->GetObject("h_pm00", hist_pm002D_base);
  file->GetObject("h_pmpm", hist_pmpm2D_base);

  file->GetObject("h_00pm_not_weighted", hist_00pm2D_base_nw);
  file->GetObject("h_pm00_not_weighted", hist_pm002D_base_nw);
  file->GetObject("h_pmpm_not_weighted", hist_pmpm2D_base_nw);

  file->GetObject("h_00_not_weighted", hist_00);
  file->GetObject("h_pm_not_weighted", hist_pm);
  file->GetObject("h_pmpm1_not_weighted", hist_1);
  file->GetObject("h_pmpm2_not_weighted", hist_2);

  file->GetObject("deltaT_pm00_not_weighted", deltaTpm00_not_weighted);
  file->GetObject("deltaT_pm00", deltaTpm00_weighted);
  file->GetObject("deltaT_pmpm_not_weighted", deltaT12_not_weighted);
  file->GetObject("deltaT_pmpm", deltaT12_weighted);

  if (!hist_00pm2D_base || !hist_pm002D_base || !hist_pmpm2D_base)
  {
    std::cerr << "ERROR: Could not load histograms from file!" << std::endl;
    return 1;
  }

  hist_00pm2D_base->Sumw2();
  hist_pm002D_base->Sumw2();
  hist_pmpm2D_base->Sumw2();

  Double_t t2MaxMin = 0.0;
  Double_t t2MaxMax = 300.0;
  Int_t numSteps = 10;
  std::cout << "Integral range min: ";
  std::cin >> t2MaxMin;
  std::cout << "Integral range max: ";
  std::cin >> t2MaxMax;
  std::cout << "Number of steps: ";
  std::cin >> numSteps;

  // Definiuj zakresy t2Max do przeskanowania
  std::vector<Double_t> t2MaxValues;
  Double_t step = (t2MaxMax - t2MaxMin) / (Double_t)numSteps;
  for (Int_t i = 0; i <= numSteps; ++i)
  {
    t2MaxValues.push_back(t2MaxMin + i * step);
  }

  std::map<Double_t, TH1 *> hist_00pm2D_projX;
  std::map<Double_t, TH1 *> hist_pm002D_projX;
  std::map<Double_t, TH1 *> hist_pmpmRA2D_projX;
  std::map<Double_t, TH1 *> hist_pmpmRB2D_projX;

  std::map<Double_t, TH1 *> hist_RA;
  std::map<Double_t, TH1 *> hist_RB;
  std::map<Double_t, TH1 *> hist_RC;

  // Lambda for projection
  auto make_proj_x = [](TH2 *hist2D, const TString &name, Double_t t2Min, Double_t t2Max) -> TH1 *
  {
    Int_t lowerBin = hist2D->GetYaxis()->FindBin(t2Min);
    Int_t upperBin = hist2D->GetYaxis()->FindBin(t2Max);
    return hist2D->ProjectionX(name, lowerBin, upperBin);
  };

  ErrorHandling::ErrorLogs logger("log/");
  Utils::InitializeVariables(logger);
  KLOE::setGlobalStyle();

  // Pętla po zakresach t2Max
  for (Double_t t2Max : t2MaxValues)
  {
    std::cout << "Processing t2Max = " << t2Max << std::endl;

    hist_00pm2D_projX[t2Max] = make_proj_x(hist_00pm2D_base, Form("h_00pm_projX_%.0f", t2Max), 0, t2Max);
    hist_pm002D_projX[t2Max] = make_proj_x(hist_pm002D_base, Form("h_pm00_projX_%.0f", t2Max), 0, t2Max);
    hist_pmpmRA2D_projX[t2Max] = make_proj_x(hist_pmpm2D_base, Form("h_pmpmRA_projX_%.0f", t2Max), 0, t2Max);
    hist_pmpmRB2D_projX[t2Max] = make_proj_x(hist_pmpm2D_base, Form("h_pmpmRB_projX_%.0f", t2Max), 0, t2Max);

    hist_RA[t2Max] = (TH1 *)hist_pmpmRA2D_projX[t2Max]->Clone(Form("hist_RA_%.0f", t2Max));
    hist_RA[t2Max]->Divide(hist_00pm2D_projX[t2Max]);

    hist_RB[t2Max] = (TH1 *)hist_pmpmRB2D_projX[t2Max]->Clone(Form("hist_RB_%.0f", t2Max));
    hist_RB[t2Max]->Divide(hist_pm002D_projX[t2Max]);

    hist_RC[t2Max] = (TH1 *)hist_RA[t2Max]->Clone(Form("hist_RC_%.0f", t2Max));
    hist_RC[t2Max]->Divide(hist_RB[t2Max]);
  }

  // Rysuj oryginalne histogramy 2D (tylko raz)
  TCanvas *c2D = new TCanvas("c2D", "2D Histograms", 1800, 600);
  c2D->Divide(3, 1);
  c2D->cd(1);
  hist_00pm2D_base->SetTitle(Form("h_00pm - %s events", maxEventsTitle.Data()));
  hist_00pm2D_base->GetXaxis()->SetRangeUser(0, 300.0);
  hist_00pm2D_base->GetYaxis()->SetRangeUser(0, 300.0);
  hist_00pm2D_base->Draw("COLZ");
  c2D->cd(2);
  hist_pm002D_base->SetTitle(Form("h_pm00 - %s events", maxEventsTitle.Data()));
  hist_pm002D_base->GetXaxis()->SetRangeUser(0, 300.0);
  hist_pm002D_base->GetYaxis()->SetRangeUser(0, 300.0);
  hist_pm002D_base->Draw("COLZ");
  c2D->cd(3);
  hist_pmpm2D_base->SetTitle(Form("h_pmpm - %s events", maxEventsTitle.Data()));
  hist_pmpm2D_base->GetXaxis()->SetRangeUser(0, 300.0);
  hist_pmpm2D_base->GetYaxis()->SetRangeUser(0, 300.0);
  hist_pmpm2D_base->Draw("COLZ");
  c2D->SaveAs(Form(imgFolderPath + "/2D_histograms_%s.svg", maxEventsFile.Data()));

  // Rysuj projekcje 1D dla różnych t2Max
  TCanvas *cProjX = new TCanvas("cProjX", "1D Projections", 1800, 600);
  cProjX->Divide(3, 1);
  for (Double_t t2Max : t2MaxValues)
  {
    cProjX->cd(1);
    gPad->SetLogy();
    hist_00pm2D_projX[t2Max]->SetTitle(Form("h_00pm projX - t2Max=%.0f", t2Max));
    hist_00pm2D_projX[t2Max]->SetLineColor(kBlue);
    hist_00pm2D_projX[t2Max]->GetXaxis()->SetRangeUser(0, 20.0);
    hist_00pm2D_projX[t2Max]->Draw();

    cProjX->cd(2);
    gPad->SetLogy();
    hist_pm002D_projX[t2Max]->SetTitle(Form("h_pm00 projX - t2Max=%.0f", t2Max));
    hist_pm002D_projX[t2Max]->SetLineColor(kBlue);
    hist_pm002D_projX[t2Max]->GetXaxis()->SetRangeUser(0, 20.0);
    hist_pm002D_projX[t2Max]->Draw();

    cProjX->cd(3);
    gPad->SetLogy();
    hist_pmpmRA2D_projX[t2Max]->SetTitle(Form("h_pmpmRA projX - t2Max=%.0f", t2Max));
    hist_pmpmRA2D_projX[t2Max]->SetLineColor(kBlue);
    hist_pmpmRA2D_projX[t2Max]->GetXaxis()->SetRangeUser(0, 20.0);
    hist_pmpmRA2D_projX[t2Max]->Draw();

    cProjX->SaveAs(Form(imgFolderPath + "/1DprojX_histograms_t2Max%.0f.svg", t2Max));
  }

  /////////////////////////////////////////////////////////////////
  // Fitting ratios R_{A}, R_{B}, R_{C}

  Double_t reParam = PhysicsConstants::Re;
  Double_t imParam = PhysicsConstants::Im_nonCPT;

  TF2 *func_00pm = new TF2("I(#pi^{0}#pi^{0},t_{1},#pi^{+}#pi^{-},t_{2});t_{1} [#tau_{S}]; t_{2} [#tau_{S}]", &interf_function_00pm, 0.0, 300, 0.0, 300, 2);
  func_00pm->SetParameters(reParam, imParam);

  TF2 *func_pm00 = new TF2("I(#pi^{+}#pi^{-},t_{1},#pi^{0}#pi^{0},t_{2});t_{1} [#tau_{S}]; t_{2} [#tau_{S}]", &interf_function_pm00, 0.0, 300, 0.0, 300, 2);
  func_pm00->SetParameters(reParam, imParam);

  TF2 *func_pmpm = new TF2("I(#pi^{+}#pi^{-},t_{1},#pi^{+}#pi^{-},t_{2});t_{1} [#tau_{S}]; t_{2} [#tau_{S}]", &interf_function_pmpm, 0.0, 300, 0.0, 300, 2);
  func_pmpm->SetParameters(reParam, imParam);

  auto func_00pm_normalized = [&](Double_t *x, Double_t *par)
  {
    // Pas poziomy: y jest małe (0-6), x może być dowolne (0-20)
    Bool_t inHorizontalStrip = (x[1] > 0.0 && x[1] < 6.0);

    // Pas pionowy: x jest małe (0-6), y może być dowolne (0-20)
    Bool_t inVerticalStrip = (x[0] > 0.0 && x[0] < 6.0);

    // Jeśli punkt nie leży w ŻADNYM z tych pasów, odrzuć go
    if (!(inHorizontalStrip || inVerticalStrip))
    {
      TF2::RejectPoint();
      return 0.0;
    }

    return par[2] * interf_function_00pm(x, par);
  };

  auto func_pm00_normalized = [&](Double_t *x, Double_t *par)
  {
    // Pas poziomy: y jest małe (0-6), x może być dowolne (0-20)
    Bool_t inHorizontalStrip = (x[1] > 0.0 && x[1] < 6.0);

    // Pas pionowy: x jest małe (0-6), y może być dowolne (0-20)
    Bool_t inVerticalStrip = (x[0] > 0.0 && x[0] < 6.0);

    // Jeśli punkt nie leży w ŻADNYM z tych pasów, odrzuć go
    if (!(inHorizontalStrip || inVerticalStrip))
    {
      TF2::RejectPoint();
      return 0.0;
    }

    return par[2] * interf_function_pm00(x, par);
  };

  auto func_pmpm_normalized = [&](Double_t *x, Double_t *par)
  {
    // Pas poziomy: y jest małe (0-6), x może być dowolne (0-20)
    Bool_t inHorizontalStrip = (x[1] > 0.0 && x[1] < 6.0);

    // Pas pionowy: x jest małe (0-6), y może być dowolne (0-20)
    Bool_t inVerticalStrip = (x[0] > 0.0 && x[0] < 6.0);

    // Jeśli punkt nie leży w ŻADNYM z tych pasów, odrzuć go
    if (!(inHorizontalStrip || inVerticalStrip))
    {
      TF2::RejectPoint();
      return 0.0;
    }

    return par[1] * interf_function_pmpm(x, par);
  };

  TF2 *func_00pm_norm = new TF2("I(#pi^{0}#pi^{0},t_{1},#pi^{+}#pi^{-},t_{2});t_{1} [#tau_{S}]; t_{2} [#tau_{S}]", &func_00pm_normalized, 0.0, 300, 0.0, 300, 3);
  func_00pm_norm->SetParameters(reParam, imParam, 100.0);

  func_00pm_norm->FixParameter(0, reParam);
  func_00pm_norm->FixParameter(1, imParam);

  TF2 *func_pm00_norm = new TF2("I(#pi^{+}#pi^{-},t_{1},#pi^{0}#pi^{0},t_{2});t_{1} [#tau_{S}]; t_{2} [#tau_{S}]", &func_pm00_normalized, 0.0, 300, 0.0, 300, 3);
  func_pm00_norm->SetParameters(reParam, imParam, 100.0);

  func_pm00_norm->FixParameter(0, reParam);
  func_pm00_norm->FixParameter(1, imParam);

  TF2 *func_pmpm_norm = new TF2("I(#pi^{+}#pi^{-},t_{1},#pi^{+}#pi^{-},t_{2});t_{1} [#tau_{S}]; t_{2} [#tau_{S}]", &func_pmpm_normalized, 0.0, 300, 0.0, 300, 2);
  func_pmpm_norm->SetParameters(reParam, 100.0);

  func_pmpm_norm->FixParameter(0, reParam);

  std::cout << "Mock function fitting:" << std::endl;
  TFitResultPtr fitResult_00pm = hist_00pm2D_base->Fit(func_00pm_norm, "RSL");
  TFitResultPtr fitResult_pm00 = hist_pm002D_base->Fit(func_pm00_norm, "RSL");
  TFitResultPtr fitResult_pmpm = hist_pmpm2D_base->Fit(func_pmpm_norm, "RSL");

  if (fitResult_00pm->IsValid())
    std::cout << "Fit 00pm successful: Re = " << func_00pm_norm->GetParameter(0) << ", Im = " << func_00pm_norm->GetParameter(1) << ", Chi2: " << fitResult_00pm->Chi2() << std::endl;
  else
    std::cout << "Fit 00pm failed!" << std::endl;

  if (fitResult_pm00->IsValid())
    std::cout << "Fit pm00 successful: Re = " << func_pm00_norm->GetParameter(0) << ", Im = " << func_pm00_norm->GetParameter(1) << ", Chi2: " << fitResult_pm00->Chi2() << std::endl;
  else
    std::cout << "Fit pm00 failed!" << std::endl;

  if (fitResult_pmpm->IsValid())
    std::cout << "Fit pmpm successful: Re = " << func_pmpm_norm->GetParameter(0) << ", Chi2: " << fitResult_pmpm->Chi2() << std::endl;
  else
    std::cout << "Fit pmpm failed!" << std::endl;

  TH1D *hResiduals_00pm = new TH1D("hResiduals_00pm", "Residuals (00pm) (Data-Fit)/#sigma;Residual;Counts", 120, -6, 6);
  TH1D *hResiduals_pm00 = new TH1D("hResiduals_pm00", "Residuals (pm00) (Data-Fit)/#sigma;Residual;Counts", 120, -6, 6);
  TH1D *hResiduals_pmpm = new TH1D("hResiduals_pmpm", "Residuals (pmpm) (Data-Fit)/#sigma;Residual;Counts", 120, -6, 6);

  hResiduals_00pm->Sumw2();
  hResiduals_pm00->Sumw2();
  hResiduals_pmpm->Sumw2();

  auto fillResiduals2D = [&](TH2 *hist2D, TH2 *hist2DAux, TF2 *func, TH1D *outHist)
  {
    if (!hist2D || !func || !outHist)
      return;

    const int nx = hist2D->GetNbinsX();
    const int ny = hist2D->GetNbinsY();

    const double zMax = hist2D->GetMaximum();
    const double fitMax = func->GetMaximum(0.0, 0.0); // Zakładamy, że maksimum funkcji jest w pobliżu (0,0)

    for (int ix = 1; ix <= nx; ++ix)
    {
      for (int iy = 1; iy <= ny; ++iy)
      {
        const double x = hist2D->GetXaxis()->GetBinCenter(ix);
        const double y = hist2D->GetYaxis()->GetBinCenter(iy);

        const int entries = hist2DAux->GetBinContent(ix, iy);

        if (entries < 20)
          continue;

        Bool_t inStripH = (y < 6.0) && (y > 0.0);
        Bool_t inStripV = (x < 6.0) && (x > 0.0);

        if (!(inStripH || inStripV))
          continue;

        const double data = hist2D->GetBinContent(ix, iy);
        const double err = hist2D->GetBinError(ix, iy);

        const double fit = func->Eval(x, y);
        const double resid = (err > 0) ? (data - fit) / err : 0.0;

        outHist->Fill(resid);
      }
    }
  };

  fillResiduals2D(hist_00pm2D_base, hist_00pm2D_base_nw, func_00pm_norm, hResiduals_00pm);
  fillResiduals2D(hist_pm002D_base, hist_pm002D_base_nw, func_pm00_norm, hResiduals_pm00);
  fillResiduals2D(hist_pmpm2D_base, hist_pmpm2D_base_nw, func_pmpm_norm, hResiduals_pmpm);

  gStyle->SetOptStat(1111); // Entries, Mean, RMS (Entries będą widoczne)
  gStyle->SetOptFit(1);

  TF1 *fitRes00pm = new TF1("fitRes00pm", "gaus", -6, 6);
  TF1 *fitRespm00 = new TF1("fitRespm00", "gaus", -6, 6);
  TF1 *fitRespmpm = new TF1("fitRespmpm", "gaus", -6, 6);

  hResiduals_00pm->Fit(fitRes00pm, "RS");
  hResiduals_pm00->Fit(fitRespm00, "RS");
  hResiduals_pmpm->Fit(fitRespmpm, "RS");

  TCanvas *cRes00pm = new TCanvas("cRes00pm", "Residuals 00pm", 800, 600);
  hResiduals_00pm->Draw();
  gPad->Update();
  if (auto *st = (TPaveStats *)hResiduals_00pm->FindObject("stats"))
  {
    st->SetX1NDC(0.72);
    st->SetX2NDC(0.93);
    st->SetY1NDC(0.77);
    st->SetY2NDC(0.93);
  }
  cRes00pm->SaveAs(Form(imgFolderPath + "/residuals_00pm_%s.svg", maxEventsFile.Data()));

  TCanvas *cRespm00 = new TCanvas("cRespm00", "Residuals pm00", 800, 600);
  hResiduals_pm00->Draw();
  gPad->Update();
  if (auto *st = (TPaveStats *)hResiduals_pm00->FindObject("stats"))
  {
    st->SetX1NDC(0.72);
    st->SetX2NDC(0.93);
    st->SetY1NDC(0.77);
    st->SetY2NDC(0.93);
  }
  cRespm00->SaveAs(Form(imgFolderPath + "/residuals_pm00_%s.svg", maxEventsFile.Data()));

  TCanvas *cRespmpm = new TCanvas("cRespmpm", "Residuals pmpm", 800, 600);
  hResiduals_pmpm->Draw();
  gPad->Update();
  if (auto *st = (TPaveStats *)hResiduals_pmpm->FindObject("stats"))
  {
    st->SetX1NDC(0.72);
    st->SetX2NDC(0.93);
    st->SetY1NDC(0.77);
    st->SetY2NDC(0.93);
  }
  cRespmpm->SaveAs(Form(imgFolderPath + "/residuals_pmpm_%s.svg", maxEventsFile.Data()));

  auto func_00pm_1D = [&](Double_t *x, Double_t *par)
  {
    Double_t current_x = x[0];

    func_00pm->SetParameters(par[0], par[1]);

    auto temp_f1_lambda = [&](Double_t *t, Double_t *par)
    {
      return func_00pm->Eval(current_x, t[0]);
    };

    TF1 temp_f1("temp_f1_integral", temp_f1_lambda, 0.0, par[2], 0);

    Double_t integral_result = temp_f1.Integral(0.0, par[2]);

    return integral_result;
  };

  auto func_pm00_1D = [&](Double_t *x, Double_t *par)
  {
    Double_t current_x = x[0];

    func_pm00->SetParameters(par[0], par[1]);

    auto temp_f1_lambda = [&](Double_t *t, Double_t *par)
    {
      return func_pm00->Eval(current_x, t[0]);
    };

    // Utwórz tymczasowy TF1 z funkcji temp_f1_lambda
    // Uwaga: Nowe obiekty TF1 powinny mieć unikalne nazwy!
    TF1 temp_f1("temp_f1_integral", temp_f1_lambda, 0.0, par[2], 0);

    // Użyj TF1::Integral do obliczenia całki po y (od y_min do y_max)
    Double_t integral_result = temp_f1.Integral(0.0, par[2]);

    return integral_result;
  };

  auto func_pmpm_1D = [&](Double_t *x, Double_t *par)
  {
    Double_t current_x = x[0];

    func_pmpm->SetParameters(par[0], par[1]);

    auto temp_f1_lambda = [&](Double_t *t, Double_t *par)
    {
      // F_xy->GetParameter(i) to wartość p[i]
      // par[0] to wartość x - używamy, gdyby x było w p
      return func_pmpm->Eval(current_x, t[0]);
    };

    // Utwórz tymczasowy TF1 z funkcji temp_f1_lambda
    // Uwaga: Nowe obiekty TF1 powinny mieć unikalne nazwy!
    TF1 temp_f1("temp_f1_integral", temp_f1_lambda, 0.0, par[2], 0);

    // Użyj TF1::Integral do obliczenia całki po y (od y_min do y_max)
    Double_t integral_result = temp_f1.Integral(0.0, par[2]);

    return integral_result;
  };

  auto func_pmpm_00pm_1D = [&](Double_t *x, Double_t *par)
  {
    return func_pmpm_1D(x, par) / func_00pm_1D(x, par);
  };

  auto func_pmpm_pm00_1D = [&](Double_t *x, Double_t *par)
  {
    return func_pmpm_1D(x, par) / func_pm00_1D(x, par);
  };

  auto func_RA_RB_1D = [&](Double_t *x, Double_t *par)
  {
    return func_pmpm_00pm_1D(x, par) / func_pmpm_pm00_1D(x, par);
  };

  Double_t imLimitPercentage = 1.0; // 50% od wartości bezwzględnej imParam
  Double_t reLimitPercentage = 0.1; // 10% od wartości reParam

  std::map<std::string, Double_t *> fitLimits = {
      {"Re", new Double_t[2]{reParam - reLimitPercentage * reParam, reParam + reLimitPercentage * reParam}},
      {"Im", new Double_t[2]{imParam - imLimitPercentage * abs(imParam), imParam + imLimitPercentage * abs(imParam)}},
      {"Range", new Double_t[2]{0.0, 330.0}}};

  TF1 *fitFunc_RA = new TF1("fitFunc_RA", func_pmpm_00pm_1D, 0.0, 20.0, 3);
  fitFunc_RA->SetParameters(reParam, imParam, 300.0);
  fitFunc_RA->SetParNames("Re", "Im", "Range");
  fitFunc_RA->SetParLimits(0, fitLimits["Re"][0], fitLimits["Re"][1]);
  fitFunc_RA->SetParLimits(1, fitLimits["Im"][0], fitLimits["Im"][1]);
  fitFunc_RA->SetParLimits(2, fitLimits["Range"][0], fitLimits["Range"][1]);

  TF1 *fitFunc_RB = new TF1("fitFunc_RB", func_pmpm_pm00_1D, 0.0, 20.0, 3);
  fitFunc_RB->SetParameters(reParam, imParam, 300.0);
  fitFunc_RB->SetParNames("Re", "Im", "Range");
  fitFunc_RB->SetParLimits(0, fitLimits["Re"][0], fitLimits["Re"][1]);
  fitFunc_RB->SetParLimits(1, fitLimits["Im"][0], fitLimits["Im"][1]);
  fitFunc_RB->SetParLimits(2, fitLimits["Range"][0], fitLimits["Range"][1]);

  TF1 *fitFunc_RC = new TF1("fitFunc_RC", func_RA_RB_1D, 0.0, 20.0, 3);
  fitFunc_RC->SetParameters(reParam, imParam, 300.0);
  fitFunc_RC->SetParNames("Re", "Im", "Range");
  fitFunc_RC->SetParLimits(0, fitLimits["Re"][0], fitLimits["Re"][1]);
  fitFunc_RC->SetParLimits(1, fitLimits["Im"][0], fitLimits["Im"][1]);
  fitFunc_RC->SetParLimits(2, fitLimits["Range"][0], fitLimits["Range"][1]);

  // Prepare plots for fitted parameters vs t2Max
  TGraphErrors *graphRe_RA = new TGraphErrors();
  graphRe_RA->SetTitle("Fitted Re parameter for R_{A} vs t_{2}^{max};t_{2}^{max} [#tau_{S}];Fitted Re");
  TGraphErrors *graphIm_RA = new TGraphErrors();
  graphIm_RA->SetTitle("Fitted Im parameter for R_{A} vs t_{2}^{max};t_{2}^{max} [#tau_{S}];Fitted Im");
  TGraph *graphRange_RA = new TGraph();
  graphRange_RA->SetTitle("Fitted Range parameter for R_{A} vs t_{2}^{max};t_{2}^{max} [#tau_{S}];Fitted Range");

  TGraphErrors *graphReError_RA = new TGraphErrors();
  graphReError_RA->SetTitle("Fitted Re Error for R_{A} vs t_{2}^{max};t_{2}^{max} [#tau_{S}];Error on Re");
  TGraphErrors *graphImError_RA = new TGraphErrors();
  graphImError_RA->SetTitle("Fitted Im Error for R_{A} vs t_{2}^{max};t_{2}^{max} [#tau_{S}];Error on Im");

  TGraphErrors *graphRe_RB = new TGraphErrors();
  graphRe_RB->SetTitle("Fitted Re parameter for R_{B} vs t_{2}^{max};t_{2}^{max} [#tau_{S}];Fitted Re");
  TGraphErrors *graphIm_RB = new TGraphErrors();
  graphIm_RB->SetTitle("Fitted Im parameter for R_{B} vs t_{2}^{max};t_{2}^{max} [#tau_{S}];Fitted Im");
  TGraph *graphReError_RB = new TGraph();
  graphReError_RB->SetTitle("Fitted Re Error for R_{B} vs t_{2}^{max};t_{2}^{max} [#tau_{S}];Error on Re");
  TGraph *graphImError_RB = new TGraph();
  graphImError_RB->SetTitle("Fitted Im Error for R_{B} vs t_{2}^{max};t_{2}^{max} [#tau_{S}];Error on Im");

  TGraphErrors *graphRe_RC = new TGraphErrors();
  graphRe_RC->SetTitle("Fitted Re parameter for R_{C} vs t_{2}^{max};t_{2}^{max} [#tau_{S}];Fitted Re");
  TGraphErrors *graphIm_RC = new TGraphErrors();
  graphIm_RC->SetTitle("Fitted Im parameter for R_{C} vs t_{2}^{max};t_{2}^{max} [#tau_{S}];Fitted Im");
  TGraph *graphReError_RC = new TGraph();
  graphReError_RC->SetTitle("Fitted Re Error for R_{C} vs t_{2}^{max};t_{2}^{max} [#tau_{S}];Error on Re");
  TGraph *graphImError_RC = new TGraph();
  graphImError_RC->SetTitle("Fitted Im Error for R_{C} vs t_{2}^{max};t_{2}^{max} [#tau_{S}];Error on Im");

  TCanvas *cR = new TCanvas("cR", "Ratios RA, RB, RC", 1800, 600);
  cR->Divide(3, 1);

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  TLine *lineZero = new TLine(30, 0, 320., 0);
  lineZero->SetLineColor(kBlack);
  lineZero->SetLineStyle(2);

  std::map<std::string, TLine *> linesRe = {
      {"Min", new TLine(30, fitLimits["Re"][0] - reParam, 320., fitLimits["Re"][0] - reParam)},
      {"Max", new TLine(30, fitLimits["Re"][1] - reParam, 320., fitLimits["Re"][1] - reParam)}};

  linesRe["Min"]->SetLineColor(kRed);
  linesRe["Min"]->SetLineStyle(3);
  linesRe["Max"]->SetLineColor(kRed);
  linesRe["Max"]->SetLineStyle(3);

  std::map<std::string, TLine *> linesIm = {
      {"Min", new TLine(30, fitLimits["Im"][0] - imParam, 320., fitLimits["Im"][0] - imParam)},
      {"Max", new TLine(30, fitLimits["Im"][1] - imParam, 320., fitLimits["Im"][1] - imParam)}};

  linesIm["Min"]->SetLineColor(kRed);
  linesIm["Min"]->SetLineStyle(3);
  linesIm["Max"]->SetLineColor(kRed);
  linesIm["Max"]->SetLineStyle(3);

  TFitResultPtr fitResultRA, fitResultRB, fitResultRC;
  for (Double_t t2Max : t2MaxValues)
  {
    // Ustaw parametr Range w funkcji fitowej zgodnie z aktualnym t2Max (jako hint dla fitowania)
    fitFunc_RA->SetParameter(2, t2Max);
    fitFunc_RB->SetParameter(2, t2Max);
    fitFunc_RC->SetParameter(2, t2Max);

    cR->cd(1);
    hist_RA[t2Max]->SetTitle(Form("R_{A} - t_{2}^{max}=%.0f #tau_{S}", t2Max));
    hist_RA[t2Max]->GetXaxis()->SetRangeUser(0, 20);
    fitResultRA = hist_RA[t2Max]->Fit(fitFunc_RA, "RSEM");
    hist_RA[t2Max]->SetMarkerStyle(20);
    hist_RA[t2Max]->SetLineWidth(2);
    hist_RA[t2Max]->SetMarkerColor(kBlack);
    hist_RA[t2Max]->SetLineColor(kBlack);

    hist_RA[t2Max]->Draw();
    if (fitResultRA->IsValid())
    {
      graphRe_RA->SetPoint(graphRe_RA->GetN(), t2Max, fitResultRA->Parameter(0) - reParam);
      graphRe_RA->SetPointError(graphRe_RA->GetN() - 1, 0, fitResultRA->Error(0));
      graphIm_RA->SetPoint(graphIm_RA->GetN(), t2Max, fitResultRA->Parameter(1) - imParam);
      graphIm_RA->SetPointError(graphIm_RA->GetN() - 1, 0, fitResultRA->Error(1));
      graphRange_RA->SetPoint(graphRange_RA->GetN(), t2Max, fitResultRA->Parameter(2) - t2Max);

      graphReError_RA->SetPoint(graphReError_RA->GetN(), t2Max, fitResultRA->Error(0));
      graphImError_RA->SetPoint(graphImError_RA->GetN(), t2Max, fitResultRA->Error(1));

      // Dodaj box z parametrami w prawym dolnym rogu
      TPaveText *ptRA = new TPaveText(0.55, 0.15, 0.85, 0.35, "NDC");
      ptRA->SetLineWidth(1);
      ptRA->SetFillColor(kWhite);
      ptRA->SetBorderSize(1);
      ptRA->AddText(Form("Re = %.2g #pm %.2g", fitResultRA->Parameter(0), fitResultRA->Error(0)));
      ptRA->AddText(Form("Im = %.2g #pm %.2g", fitResultRA->Parameter(1), fitResultRA->Error(1)));
      ptRA->AddText(Form("Range = %.2g #pm %.2g", fitResultRA->Parameter(2), fitResultRA->Error(2)));
      ptRA->Draw();
    }
    else
    {
      std::cout << "Fit for R_{A} (t2Max=" << t2Max << ") was not valid." << std::endl;
    }

    cR->cd(2);
    hist_RB[t2Max]->SetTitle(Form("R_{B} - t_{2}^{max}=%.0f #tau_{S}", t2Max));
    hist_RB[t2Max]->GetXaxis()->SetRangeUser(0, 20);
    fitResultRB = hist_RB[t2Max]->Fit(fitFunc_RB, "RSEM");
    hist_RB[t2Max]->SetMarkerStyle(20);
    hist_RB[t2Max]->SetLineWidth(2);
    hist_RB[t2Max]->SetMarkerColor(kBlack);
    hist_RB[t2Max]->SetLineColor(kBlack);
    hist_RB[t2Max]->Draw();

    if (fitResultRB->IsValid())
    {
      graphRe_RB->SetPoint(graphRe_RB->GetN(), t2Max, fitResultRB->Parameter(0) - reParam);
      graphRe_RB->SetPointError(graphRe_RB->GetN() - 1, 0, fitResultRB->Error(0));
      graphIm_RB->SetPoint(graphIm_RB->GetN(), t2Max, fitResultRB->Parameter(1) - imParam);
      graphIm_RB->SetPointError(graphIm_RB->GetN() - 1, 0, fitResultRB->Error(1));

      graphReError_RB->SetPoint(graphReError_RB->GetN(), t2Max, fitResultRB->Error(0));
      graphImError_RB->SetPoint(graphImError_RB->GetN(), t2Max, fitResultRB->Error(1));

      // Dodaj box z parametrami w prawym górnym rogu
      TPaveText *ptRB = new TPaveText(0.55, 0.70, 0.85, 0.90, "NDC");
      ptRB->SetLineWidth(1);
      ptRB->SetFillColor(kWhite);
      ptRB->SetBorderSize(1);
      ptRB->AddText(Form("Re = %.2g #pm %.2g", fitResultRB->Parameter(0), fitResultRB->Error(0)));
      ptRB->AddText(Form("Im = %.2g #pm %.2g", fitResultRB->Parameter(1), fitResultRB->Error(1)));
      ptRB->AddText(Form("Range = %.2g #pm %.2g", fitResultRB->Parameter(2), fitResultRB->Error(2)));
      ptRB->Draw();
    }
    else
    {
      std::cout << "Fit for R_{B} (t2Max=" << t2Max << ") was not valid." << std::endl;
    }

    cR->cd(3);
    hist_RC[t2Max]->SetTitle(Form("R_{C} - t_{2}^{max}=%.0f #tau_{S}", t2Max));
    hist_RC[t2Max]->GetXaxis()->SetRangeUser(0, 20);
    fitResultRC = hist_RC[t2Max]->Fit(fitFunc_RC, "RSEM");
    hist_RC[t2Max]->SetMarkerStyle(20);
    hist_RC[t2Max]->SetLineWidth(2);
    hist_RC[t2Max]->SetMarkerColor(kBlack);
    hist_RC[t2Max]->SetLineColor(kBlack);
    hist_RC[t2Max]->Draw();

    if (fitResultRC->IsValid())
    {
      graphRe_RC->SetPoint(graphRe_RC->GetN(), t2Max, fitResultRC->Parameter(0) - reParam);
      graphRe_RC->SetPointError(graphRe_RC->GetN() - 1, 0, fitResultRC->Error(0));
      graphIm_RC->SetPoint(graphIm_RC->GetN(), t2Max, fitResultRC->Parameter(1) - imParam);
      graphIm_RC->SetPointError(graphIm_RC->GetN() - 1, 0, fitResultRC->Error(1));

      graphReError_RC->SetPoint(graphReError_RC->GetN(), t2Max, fitResultRC->Error(0));
      graphImError_RC->SetPoint(graphImError_RC->GetN(), t2Max, fitResultRC->Error(1));

      // Dodaj box z parametrami w prawym dolnym rogu
      TPaveText *ptRC = new TPaveText(0.55, 0.15, 0.85, 0.35, "NDC");
      ptRC->SetLineWidth(1);
      ptRC->SetFillColor(kWhite);
      ptRC->SetBorderSize(1);
      ptRC->AddText(Form("Re = %.2g #pm %.2g", fitResultRC->Parameter(0), fitResultRC->Error(0)));
      ptRC->AddText(Form("Im = %.2g #pm %.2g", fitResultRC->Parameter(1), fitResultRC->Error(1)));
      ptRC->AddText(Form("Range = %.2g #pm %.2g", fitResultRC->Parameter(2), fitResultRC->Error(2)));
      ptRC->Draw();
    }
    else
    {
      std::cout << "Fit for R_{C} (t2Max=" << t2Max << ") was not valid." << std::endl;
    }

    cR->SaveAs(Form(imgFolderPath + "/R_histograms_t2Max%.0f.svg", t2Max));

    // Reset parametrów do wartości początkowych
    fitFunc_RA->SetParameter(0, reParam);
    fitFunc_RA->SetParameter(1, imParam);
    fitFunc_RB->SetParameter(0, reParam);
    fitFunc_RB->SetParameter(1, imParam);
    fitFunc_RC->SetParameter(0, reParam);
    fitFunc_RC->SetParameter(1, imParam);
  }

  auto getGraphLimit = [](Bool_t maxOrMin, TGraph *graph1, TGraph *graph2, TGraph *graph3) -> Double_t
  {
    Double_t extreme = 0.0;
    Double_t values[3] = {0.0, 0.0, 0.0};

    if (maxOrMin)
    {
      for (int i = 0; i < graph1->GetN(); ++i)
      {
        Double_t y1, y2, y3;
        y1 = graph1->GetPointY(i);
        y2 = graph2->GetPointY(i);
        y3 = graph3->GetPointY(i);

        values[0] = std::max(values[0], y1);
        values[1] = std::max(values[1], y2);
        values[2] = std::max(values[2], y3);
      }

      extreme = *std::max_element(values, values + 3);

      return extreme + 0.15 * abs(extreme); // Dodaj 15% marginesu powyżej maksimum
    }
    else
    {
      for (int i = 0; i < graph1->GetN(); ++i)
      {
        Double_t y1, y2, y3;
        y1 = graph1->GetPointY(i);
        y2 = graph2->GetPointY(i);
        y3 = graph3->GetPointY(i);

        values[0] = std::min(values[0], y1);
        values[1] = std::min(values[1], y2);
        values[2] = std::min(values[2], y3);
      }

      extreme = *std::min_element(values, values + 3);

      return extreme - 0.15 * abs(extreme); // Dodaj 15% marginesu poniżej minimum
    }
  };

  // Ustaw style dla wykresów
  graphRe_RA->SetMarkerColor(kBlue);
  graphRe_RA->SetMarkerSize(1.2);
  graphRe_RA->SetMarkerStyle(20);
  graphRe_RA->SetLineColor(kBlue);

  graphIm_RA->SetMarkerColor(kRed);
  graphIm_RA->SetMarkerSize(1.2);
  graphIm_RA->SetMarkerStyle(20);
  graphIm_RA->SetLineColor(kRed);

  graphRe_RB->SetMarkerColor(kGreen + 2);
  graphRe_RB->SetMarkerSize(1.2);
  graphRe_RB->SetMarkerStyle(21);
  graphRe_RB->SetLineColor(kGreen + 2);

  graphIm_RB->SetMarkerColor(kOrange + 1);
  graphIm_RB->SetMarkerSize(1.2);
  graphIm_RB->SetMarkerStyle(21);
  graphIm_RB->SetLineColor(kOrange + 1);

  graphRe_RC->SetMarkerColor(kMagenta + 1);
  graphRe_RC->SetMarkerSize(1.2);
  graphRe_RC->SetMarkerStyle(22);
  graphRe_RC->SetLineColor(kMagenta + 1);

  graphIm_RC->SetMarkerColor(kCyan + 1);
  graphIm_RC->SetMarkerSize(1.2);
  graphIm_RC->SetMarkerStyle(22);
  graphIm_RC->SetLineColor(kCyan + 1);

  graphReError_RA->SetMarkerColor(kBlue);
  graphReError_RA->SetMarkerSize(1.2);
  graphReError_RA->SetMarkerStyle(20);
  graphReError_RA->SetLineColor(kBlue);

  graphImError_RA->SetMarkerColor(kRed);
  graphImError_RA->SetMarkerSize(1.2);
  graphImError_RA->SetMarkerStyle(20);
  graphImError_RA->SetLineColor(kRed);

  graphReError_RB->SetMarkerColor(kGreen + 2);
  graphReError_RB->SetMarkerSize(1.2);
  graphReError_RB->SetMarkerStyle(21);
  graphReError_RB->SetLineColor(kGreen + 2);

  graphImError_RB->SetMarkerColor(kOrange + 1);
  graphImError_RB->SetMarkerSize(1.2);
  graphImError_RB->SetMarkerStyle(21);
  graphImError_RB->SetLineColor(kOrange + 1);

  graphReError_RC->SetMarkerColor(kMagenta + 1);
  graphReError_RC->SetMarkerSize(1.2);
  graphReError_RC->SetMarkerStyle(22);
  graphReError_RC->SetLineColor(kMagenta + 1);

  graphImError_RC->SetMarkerColor(kCyan + 1);
  graphImError_RC->SetMarkerSize(1.2);
  graphImError_RC->SetMarkerStyle(22);
  graphImError_RC->SetLineColor(kCyan + 1);

  // Rysuj wyniki fitów vs t2Max - pojedyncze wykresy dla każdego ratio
  TCanvas *cGraphRe_RA = new TCanvas("cGraphRe_RA", "Fitted Re parameter for R_{A} vs t2Max", 800, 600);
  graphRe_RA->SetMinimum(-1e-4);
  graphRe_RA->SetMaximum(1e-4);
  graphRe_RA->GetXaxis()->SetLimits(30, 320);
  graphRe_RA->Draw("AP");
  lineZero->Draw("SAME");
  linesRe["Min"]->Draw("SAME");
  linesRe["Max"]->Draw("SAME");
  cGraphRe_RA->SaveAs(Form(imgFolderPath + "/Fitted_Re_RA_vs_t2Max.svg", maxEventsFile.Data()));

  TCanvas *cGraphRe_RB = new TCanvas("cGraphRe_RB", "Fitted Re parameter for R_{B} vs t2Max", 800, 600);
  graphRe_RB->SetMinimum(-1e-4);
  graphRe_RB->SetMaximum(1e-4);
  graphRe_RB->GetXaxis()->SetLimits(30, 320);
  graphRe_RB->Draw("AP");
  lineZero->Draw("SAME");
  linesRe["Min"]->Draw("SAME");
  linesRe["Max"]->Draw("SAME");
  cGraphRe_RB->SaveAs(Form(imgFolderPath + "/Fitted_Re_RB_vs_t2Max.svg", maxEventsFile.Data()));

  TCanvas *cGraphRe_RC = new TCanvas("cGraphRe_RC", "Fitted Re parameter for R_{C} vs t2Max", 800, 600);
  graphRe_RC->SetMinimum(-1e-4);
  graphRe_RC->SetMaximum(1e-4);
  graphRe_RC->GetXaxis()->SetLimits(30, 320);
  graphRe_RC->Draw("AP");
  lineZero->Draw("SAME");
  linesRe["Min"]->Draw("SAME");
  linesRe["Max"]->Draw("SAME");
  cGraphRe_RC->SaveAs(Form(imgFolderPath + "/Fitted_Re_RC_vs_t2Max.svg", maxEventsFile.Data()));

  TCanvas *cGraphIm_RA = new TCanvas("cGraphIm_RA", "Fitted Im parameter for R_{A} vs t2Max", 800, 600);
  graphIm_RA->SetMinimum(-0.005);
  graphIm_RA->SetMaximum(0.005);
  graphIm_RA->GetXaxis()->SetLimits(30, 320);
  graphIm_RA->Draw("AP");
  lineZero->Draw("SAME");
  linesIm["Min"]->Draw("SAME");
  linesIm["Max"]->Draw("SAME");
  cGraphIm_RA->SaveAs(Form(imgFolderPath + "/Fitted_Im_RA_vs_t2Max.svg", maxEventsFile.Data()));

  TCanvas *cGraphIm_RB = new TCanvas("cGraphIm_RB", "Fitted Im parameter for R_{B} vs t2Max", 800, 600);
  graphIm_RB->SetMinimum(-0.005);
  graphIm_RB->SetMaximum(0.005);
  graphIm_RB->GetXaxis()->SetLimits(30, 320);
  graphIm_RB->Draw("AP");
  lineZero->Draw("SAME");
  linesIm["Min"]->Draw("SAME");
  linesIm["Max"]->Draw("SAME");
  cGraphIm_RB->SaveAs(Form(imgFolderPath + "/Fitted_Im_RB_vs_t2Max.svg", maxEventsFile.Data()));

  TCanvas *cGraphIm_RC = new TCanvas("cGraphIm_RC", "Fitted Im parameter for R_{C} vs t2Max", 800, 600);
  graphIm_RC->SetMinimum(-0.005);
  graphIm_RC->SetMaximum(0.005);
  graphIm_RC->GetXaxis()->SetLimits(30, 320);
  graphIm_RC->Draw("AP");
  lineZero->Draw("SAME");
  linesIm["Min"]->Draw("SAME");
  linesIm["Max"]->Draw("SAME");
  cGraphIm_RC->SaveAs(Form(imgFolderPath + "/Fitted_Im_RC_vs_t2Max.svg", maxEventsFile.Data()));

  //
  Double_t reErrorMax = getGraphLimit(kTRUE, graphReError_RA, graphReError_RB, graphReError_RC);
  Double_t reErrorMin = getGraphLimit(kFALSE, graphReError_RA, graphReError_RB, graphReError_RC);
  Double_t imErrorMax = getGraphLimit(kTRUE, graphImError_RA, graphImError_RB, graphImError_RC);
  Double_t imErrorMin = getGraphLimit(kFALSE, graphImError_RA, graphImError_RB, graphImError_RC);
  //

  TCanvas *cGraphReError_RA = new TCanvas("cGraphReError_RA", "Error on Re parameter for R_{A} vs t2Max", 800, 600);
  graphReError_RA->SetMinimum(reErrorMin);
  graphReError_RA->SetMaximum(reErrorMax);
  graphReError_RA->GetXaxis()->SetLimits(30, 320);
  graphReError_RA->Draw("AP");
  cGraphReError_RA->SaveAs(Form(imgFolderPath + "/Fitted_ReError_RA_vs_t2Max.svg", maxEventsFile.Data()));

  TCanvas *cGraphReError_RB = new TCanvas("cGraphReError_RB", "Error on Re parameter for R_{B} vs t2Max", 800, 600);
  graphReError_RB->SetMinimum(reErrorMin);
  graphReError_RB->SetMaximum(reErrorMax);
  graphReError_RB->GetXaxis()->SetLimits(30, 320);
  graphReError_RB->Draw("AP");
  cGraphReError_RB->SaveAs(Form(imgFolderPath + "/Fitted_ReError_RB_vs_t2Max.svg", maxEventsFile.Data()));

  TCanvas *cGraphReError_RC = new TCanvas("cGraphReError_RC", "Error on Re parameter for R_{C} vs t2Max", 800, 600);
  graphReError_RC->SetMinimum(reErrorMin);
  graphReError_RC->SetMaximum(reErrorMax);
  graphReError_RC->GetXaxis()->SetLimits(30, 320);
  graphReError_RC->Draw("AP");
  cGraphReError_RC->SaveAs(Form(imgFolderPath + "/Fitted_ReError_RC_vs_t2Max.svg", maxEventsFile.Data()));

  TCanvas *cGraphImError_RA = new TCanvas("cGraphImError_RA", "Error on Im parameter for R_{A} vs t2Max", 800, 600);
  graphImError_RA->SetMinimum(imErrorMin);
  graphImError_RA->SetMaximum(imErrorMax);
  graphImError_RA->GetXaxis()->SetLimits(30, 320);
  graphImError_RA->Draw("AP");
  cGraphImError_RA->SaveAs(Form(imgFolderPath + "/Fitted_ImError_RA_vs_t2Max.svg", maxEventsFile.Data()));

  TCanvas *cGraphImError_RB = new TCanvas("cGraphImError_RB", "Error on Im parameter for R_{B} vs t2Max", 800, 600);
  graphImError_RB->SetMinimum(imErrorMin);
  graphImError_RB->SetMaximum(imErrorMax);
  graphImError_RB->GetXaxis()->SetLimits(30, 320);
  graphImError_RB->Draw("AP");
  cGraphImError_RB->SaveAs(Form(imgFolderPath + "/Fitted_ImError_RB_vs_t2Max.svg", maxEventsFile.Data()));

  TCanvas *cGraphImError_RC = new TCanvas("cGraphImError_RC", "Error on Im parameter for R_{C} vs t2Max", 800, 600);
  graphImError_RC->SetMinimum(imErrorMin);
  graphImError_RC->SetMaximum(imErrorMax);
  graphImError_RC->GetXaxis()->SetLimits(30, 320);
  graphImError_RC->Draw("AP");
  cGraphImError_RC->SaveAs(Form(imgFolderPath + "/Fitted_ImError_RC_vs_t2Max.svg", maxEventsFile.Data()));

  // Wykresy zbiorowe - wszystkie Re i Im na jednym wykresie
  TCanvas *cGraphRe_All = new TCanvas("cGraphRe_All", "Fitted Re parameter vs t2Max", 1200, 800);
  graphRe_RA->SetTitle("Fitted Re parameter vs t_{2}^{max};t_{2}^{max} [#tau_{S}];Fitted Re");
  graphRe_RA->Draw("AP");
  graphRe_RB->Draw("P SAME");
  graphRe_RC->Draw("P SAME");

  TLegend *legRe = new TLegend(0.7, 0.7, 0.9, 0.9);
  legRe->AddEntry(graphRe_RA, "R_{A}", "lp");
  legRe->AddEntry(graphRe_RB, "R_{B}", "lp");
  legRe->AddEntry(graphRe_RC, "R_{C}", "lp");
  legRe->Draw();
  cGraphRe_All->SaveAs(Form(imgFolderPath + "/Fitted_Re_All_vs_t2Max.svg", maxEventsFile.Data()));

  TCanvas *cGraphIm_All = new TCanvas("cGraphIm_All", "Fitted Im parameter vs t2Max", 1200, 800);
  graphIm_RA->SetTitle("Fitted Im parameter vs t_{2}^{max};t_{2}^{max} [#tau_{S}];Fitted Im");
  graphIm_RA->Draw("AP");
  graphIm_RB->Draw("P SAME");
  graphIm_RC->Draw("P SAME");

  TLegend *legIm = new TLegend(0.7, 0.7, 0.9, 0.9);
  legIm->AddEntry(graphIm_RA, "R_{A}", "lp");
  legIm->AddEntry(graphIm_RB, "R_{B}", "lp");
  legIm->AddEntry(graphIm_RC, "R_{C}", "lp");
  legIm->Draw();
  cGraphIm_All->SaveAs(Form(imgFolderPath + "/Fitted_Im_All_vs_t2Max.svg", maxEventsFile.Data()));

  // Wykresy zbiorowe - wszystkie błędy Re na jednym wykresie
  TCanvas *cGraphReError_All = new TCanvas("cGraphReError_All", "Error on Re parameter vs t2Max", 1200, 800);
  graphReError_RA->SetTitle("Error on Re parameter vs t_{2}^{max};t_{2}^{max} [#tau_{S}];Error on Re");
  graphReError_RA->Draw("AP");
  graphReError_RB->Draw("P SAME");
  graphReError_RC->Draw("P SAME");

  TLegend *legReError = new TLegend(0.7, 0.7, 0.9, 0.9);
  legReError->AddEntry(graphReError_RA, "R_{A}", "lp");
  legReError->AddEntry(graphReError_RB, "R_{B}", "lp");
  legReError->AddEntry(graphReError_RC, "R_{C}", "lp");
  legReError->Draw();
  cGraphReError_All->SaveAs(Form(imgFolderPath + "/Fitted_ReError_All_vs_t2Max.svg", maxEventsFile.Data()));

  // Wykresy zbiorowe - wszystkie błędy Im na jednym wykresie
  TCanvas *cGraphImError_All = new TCanvas("cGraphImError_All", "Error on Im parameter vs t2Max", 1200, 800);
  graphImError_RA->SetTitle("Error on Im parameter vs t_{2}^{max};t_{2}^{max} [#tau_{S}];Error on Im");
  graphImError_RA->Draw("AP");
  graphImError_RB->Draw("P SAME");
  graphImError_RC->Draw("P SAME");

  TLegend *legImError = new TLegend(0.7, 0.7, 0.9, 0.9);
  legImError->AddEntry(graphImError_RA, "R_{A}", "lp");
  legImError->AddEntry(graphImError_RB, "R_{B}", "lp");
  legImError->AddEntry(graphImError_RC, "R_{C}", "lp");
  legImError->Draw();
  cGraphImError_All->SaveAs(Form(imgFolderPath + "/Fitted_ImError_All_vs_t2Max.svg", maxEventsFile.Data()));

  TCanvas *cSingleTimes00pm = new TCanvas("cSingleTimes00pm", "Single Times Histograms", 1200, 600);
  cSingleTimes00pm->Divide(2, 1);
  cSingleTimes00pm->cd(1);
  gPad->SetLogy();
  hist_00->SetTitle(Form("t_{00} - %s events", maxEventsTitle.Data()));
  hist_00->GetYaxis()->SetRangeUser(1, 100.0 * hist_00->GetMaximum());
  hist_00->GetXaxis()->SetRangeUser(0, 50.0);
  hist_00->SetLineColor(kRed);
  hist_00->Draw();
  cSingleTimes00pm->cd(2);
  gPad->SetLogy();
  hist_pm->SetTitle(Form("t_{+-} - %s events", maxEventsTitle.Data()));
  hist_pm->GetYaxis()->SetRangeUser(1, 100.0 * hist_pm->GetMaximum());
  hist_pm->GetXaxis()->SetRangeUser(0, 50.0);
  hist_pm->SetLineColor(kRed);
  hist_pm->Draw();
  cSingleTimes00pm->SaveAs(Form(imgFolderPath + "/single_times_00pm_%s.svg", maxEventsFile.Data()));

  TCanvas *cSingleTimes12 = new TCanvas("cSingleTimes12", "Single Times Histograms", 1200, 600);
  cSingleTimes12->Divide(2, 1);
  cSingleTimes12->cd(1);
  gPad->SetLogy();
  hist_1->SetTitle(Form("t_{1} - %s events", maxEventsTitle.Data()));
  hist_1->GetYaxis()->SetRangeUser(1, 100.0 * hist_1->GetMaximum());
  hist_1->GetXaxis()->SetRangeUser(0, 50.0);
  hist_1->SetLineColor(kRed);
  hist_1->Draw();
  cSingleTimes12->cd(2);
  gPad->SetLogy();
  hist_2->SetTitle(Form("t_{2} - %s events", maxEventsTitle.Data()));
  hist_2->GetYaxis()->SetRangeUser(1, 100.0 * hist_2->GetMaximum());
  hist_2->GetXaxis()->SetRangeUser(0, 50.0);
  hist_2->SetLineColor(kRed);
  hist_2->Draw();
  cSingleTimes12->SaveAs(Form(imgFolderPath + "/single_times_12_%s.svg", maxEventsFile.Data()));

  return 0;
}