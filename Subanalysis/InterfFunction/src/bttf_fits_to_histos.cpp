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

#include <TFitResult.h>
#include <TGraph.h>

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
  ROOT::EnableImplicitMT(8);

  std::string rootFolder = "root_files";
  std::map<Long64_t, TString> files = parseHistogramFiles(rootFolder);

  std::map<Long64_t, TH2 *> hist_00pm2D;
  std::map<Long64_t, TH2 *> hist_pm002D;
  std::map<Long64_t, TH2 *> hist_pmpm2D;

  std::map<Long64_t, TH1 *> hist_00pm2D_projX;
  std::map<Long64_t, TH1 *> hist_pm002D_projX;
  std::map<Long64_t, TH1 *> hist_pmpmRA2D_projX;
  std::map<Long64_t, TH1 *> hist_pmpmRB2D_projX;

  std::map<Long64_t, TH1 *> hist_RA;
  std::map<Long64_t, TH1 *> hist_RB;
  std::map<Long64_t, TH1 *> hist_RC;

  // Lambda for projection
  auto make_proj_x = [](TH2 *hist2D, const TString &name, Double_t t2Min, Double_t t2Max) -> TH1 *
  {
    Int_t lowerBin = hist2D->GetYaxis()->FindBin(t2Min);
    Int_t upperBin = hist2D->GetYaxis()->FindBin(t2Max);
    return hist2D->ProjectionX(name, lowerBin, upperBin);
  };

  KLOE::setGlobalStyle();

  for (const auto &entry : files)
  {
    TString filePath = Form("%s/%s", rootFolder.c_str(), entry.second.Data());

    TFile *file = TFile::Open(filePath);
    if (!file || file->IsZombie())
    {
      std::cerr << "Error opening file: " << filePath << std::endl;
      continue;
    }
    file->cd();

    file->GetObject("h_00pm", hist_00pm2D[entry.first]);
    file->GetObject("h_pm00", hist_pm002D[entry.first]);
    file->GetObject("h_pmpm", hist_pmpm2D[entry.first]);

    hist_00pm2D[entry.first]->Sumw2();
    hist_pm002D[entry.first]->Sumw2();
    hist_pmpm2D[entry.first]->Sumw2();

    hist_00pm2D_projX[entry.first] = make_proj_x(hist_00pm2D[entry.first], Form("h_00pm_projX_%lld", entry.first), 0, 300.0);
    hist_pm002D_projX[entry.first] = make_proj_x(hist_pm002D[entry.first], Form("h_pm00_projX_%lld", entry.first), 0, 300.0);
    hist_pmpmRA2D_projX[entry.first] = make_proj_x(hist_pmpm2D[entry.first], Form("h_pmpmRA_projX_%lld", entry.first), 0, 300.0);
    hist_pmpmRB2D_projX[entry.first] = make_proj_x(hist_pmpm2D[entry.first], Form("h_pmpmRB_projX_%lld", entry.first), 0, 300.0);

    hist_RA[entry.first] = (TH1 *)hist_pmpmRA2D_projX[entry.first]->Clone(Form("hist_RA_%lld", entry.first));
    hist_RA[entry.first]->Divide(hist_00pm2D_projX[entry.first]);

    hist_RB[entry.first] = (TH1 *)hist_pmpmRB2D_projX[entry.first]->Clone(Form("hist_RB_%lld", entry.first));
    hist_RB[entry.first]->Divide(hist_pm002D_projX[entry.first]);

    hist_RC[entry.first] = (TH1 *)hist_RA[entry.first]->Clone(Form("hist_RC_%lld", entry.first));
    hist_RC[entry.first]->Divide(hist_RB[entry.first]);
  }

  TCanvas *c2D = new TCanvas("c2D", "2D Histograms", 1800, 600);
  c2D->Divide(3, 1);
  for (const auto &entry : files)
  {
    c2D->cd(1);
    hist_00pm2D[entry.first]->SetTitle(Form("h_00pm - %lld events", entry.first));
    hist_00pm2D[entry.first]->GetXaxis()->SetRangeUser(0, 20);
    hist_00pm2D[entry.first]->GetYaxis()->SetRangeUser(0, 300);
    hist_00pm2D[entry.first]->Draw("COLZ");

    c2D->cd(2);
    hist_pm002D[entry.first]->SetTitle(Form("h_pm00 - %lld events", entry.first));
    hist_pm002D[entry.first]->GetXaxis()->SetRangeUser(0, 20);
    hist_pm002D[entry.first]->GetYaxis()->SetRangeUser(0, 300);
    hist_pm002D[entry.first]->Draw("COLZ");

    c2D->cd(3);
    hist_pmpm2D[entry.first]->SetTitle(Form("h_pmpm - %lld events", entry.first));
    hist_pmpm2D[entry.first]->GetXaxis()->SetRangeUser(0, 20);
    hist_pmpm2D[entry.first]->GetYaxis()->SetRangeUser(0, 300);
    hist_pmpm2D[entry.first]->Draw("COLZ");

    c2D->SaveAs(Form("img/2D_histograms_%lld.pdf", entry.first));
  }

  TCanvas *cProjX = new TCanvas("cProjX", "1D Projections", 1800, 600);
  cProjX->SetLogy();
  cProjX->Divide(3, 1);
  for (const auto &entry : files)
  {
    cProjX->cd(1);
    gPad->SetLogy();
    hist_00pm2D_projX[entry.first]->SetTitle(Form("h_00pm projX - %lld events", entry.first));
    hist_00pm2D_projX[entry.first]->GetXaxis()->SetRangeUser(0, 20);
    hist_00pm2D_projX[entry.first]->Draw();
    cProjX->cd(2);
    gPad->SetLogy();
    hist_pm002D_projX[entry.first]->SetTitle(Form("h_pm00 projX - %lld events", entry.first));
    hist_pm002D_projX[entry.first]->GetXaxis()->SetRangeUser(0, 20);
    hist_pm002D_projX[entry.first]->Draw();
    cProjX->cd(3);
    gPad->SetLogy();
    hist_pmpmRA2D_projX[entry.first]->SetTitle(Form("h_pmpmRA projX - %lld events", entry.first));
    hist_pmpmRA2D_projX[entry.first]->GetXaxis()->SetRangeUser(0, 20);
    hist_pmpmRA2D_projX[entry.first]->Draw();
    cProjX->SaveAs(Form("img/1DprojX_histograms_%lld.pdf", entry.first));
  }

  /////////////////////////////////////////////////////////////////
  // Fitting ratios R_A, R_B, R_C

  Double_t reParam = PhysicsConstants::Re;
  Double_t imParam = PhysicsConstants::Im_nonCPT;

  TF2 *func_00pm = new TF2("I(#pi^{0}#pi^{0},t_{1},#pi^{+}#pi^{-},t_{2});t_{1} [#tau_{S}]; t_{2} [#tau_{S}]", &interf_function_00pm, 0.0, 300, 0.0, 300, 2);

  TF2 *func_pm00 = new TF2("I(#pi^{+}#pi^{-},t_{1},#pi^{0}#pi^{0},t_{2});t_{1} [#tau_{S}]; t_{2} [#tau_{S}]", &interf_function_pm00, 0.0, 300, 0.0, 300, 2);

  TF2 *func_pmpm = new TF2("I(#pi^{+}#pi^{-},t_{1},#pi^{+}#pi^{-},t_{2});t_{1} [#tau_{S}]; t_{2} [#tau_{S}]", &interf_function_pmpm, 0.0, 300, 0.0, 300, 2);

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

  TF1 *fitFunc_RA = new TF1("fitFunc_RA", func_pmpm_00pm_1D, 0.0, 20.0, 3);
  fitFunc_RA->SetParameters(reParam, imParam, 300.0);
  fitFunc_RA->SetParNames("Re", "Im", "Range");
  fitFunc_RA->SetParLimits(0, reParam - 0.005, reParam + 0.005);
  fitFunc_RA->SetParLimits(1, imParam - 0.005, imParam + 0.005);
  fitFunc_RA->SetParLimits(2, 50.0, 310.0);

  TF1 *fitFunc_RB = new TF1("fitFunc_RB", func_pmpm_pm00_1D, 0.0, 20.0, 3);
  fitFunc_RB->SetParameters(reParam, imParam, 300.0);
  fitFunc_RB->SetParNames("Re", "Im", "Range");
  fitFunc_RB->SetParLimits(0, reParam - 0.005, reParam + 0.005);
  fitFunc_RB->SetParLimits(1, imParam - 0.005, imParam + 0.005);
  fitFunc_RB->SetParLimits(2, 50.0, 310.0);

  TF1 *fitFunc_RC = new TF1("fitFunc_RC", func_RA_RB_1D, 0.0, 20.0, 3);
  fitFunc_RC->SetParameters(reParam, imParam, 300.0);
  fitFunc_RC->SetParNames("Re", "Im", "Range");
  fitFunc_RC->SetParLimits(0, reParam - 0.005, reParam + 0.005);
  fitFunc_RC->SetParLimits(1, imParam - 0.005, imParam + 0.005);
  fitFunc_RC->SetParLimits(2, 50.0, 310.0);

  // Prepare plot for parameter w.r.t. number of events
  TGraph *graphRe_RA = new TGraph();
  graphRe_RA->SetTitle("Fitted Re parameter for R_A vs Number of Events;Number of Events;Fitted Re");
  TGraph *graphIm_RA = new TGraph();
  graphIm_RA->SetTitle("Fitted Im parameter for R_A vs Number of Events;Number of Events;Fitted Im");

  TCanvas *cR = new TCanvas("cR", "Ratios RA, RB, RC", 1800, 600);
  cR->Divide(3, 1);

  TFitResultPtr fitResult;
  for (const auto &entry : files)
  {
    cR->cd(1);
    hist_RA[entry.first]->SetTitle(Form("R_A - %lld events", entry.first));
    hist_RA[entry.first]->GetXaxis()->SetRangeUser(0, 20);
    fitResult = hist_RA[entry.first]->Fit(fitFunc_RA, "RSME");

    hist_RA[entry.first]->Draw();

    if (fitResult->IsValid())
    {
      graphRe_RA->SetPoint(graphRe_RA->GetN(), entry.first, fitResult->Error(0));
      graphIm_RA->SetPoint(graphIm_RA->GetN(), entry.first, fitResult->Error(1));
    }
    else
    {
      std::cout << "Fit for R_A (" << entry.first << " events) was not valid." << std::endl;
    }

    cR->cd(2);
    hist_RB[entry.first]->SetTitle(Form("R_B - %lld events", entry.first));
    hist_RB[entry.first]->GetXaxis()->SetRangeUser(0, 20);
    hist_RB[entry.first]->Fit(fitFunc_RB, "RMES");
    hist_RB[entry.first]->Draw();
    cR->cd(3);
    hist_RC[entry.first]->SetTitle(Form("R_C - %lld events", entry.first));
    hist_RC[entry.first]->GetXaxis()->SetRangeUser(0, 20);
    hist_RC[entry.first]->Fit(fitFunc_RC, "RMES");
    hist_RC[entry.first]->Draw();
    cR->SaveAs(Form("img/R_histograms_%lld.pdf", entry.first));

    fitFunc_RA->SetParameter(0, reParam);
    fitFunc_RA->SetParameter(1, imParam);
  }

  TCanvas *cGraphRe_RA = new TCanvas("cGraphRe_RA", "Fitted Re parameter for R_A vs Number of Events", 800, 600);
  graphRe_RA->SetMarkerColor(kBlue);
  graphRe_RA->SetMarkerSize(1.2);
  graphRe_RA->SetMarkerStyle(20);
  graphRe_RA->Draw("AP");
  cGraphRe_RA->SaveAs("img/Fitted_Re_RA_vs_Nevents.pdf");

  TCanvas *cGraphIm_RA = new TCanvas("cGraphIm_RA", "Fitted Im parameter for R_A vs Number of Events", 800, 600);
  graphIm_RA->SetMarkerColor(kRed);
  graphIm_RA->SetMarkerSize(1.2);
  graphIm_RA->SetMarkerStyle(20);
  graphIm_RA->Draw("AP");
  cGraphIm_RA->SaveAs("img/Fitted_Im_RA_vs_Nevents.pdf");

  return 0;
}