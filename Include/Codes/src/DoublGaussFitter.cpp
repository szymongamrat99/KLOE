#include "DoublGaussFitter.h"
#include <TStyle.h>
#include <TRandom.h>
#include <algorithm>
#include <ctime>

#include "double_gaus.h"

using namespace KLOE;

DoublGaussFitter::DoublGaussFitter()
{
  gStyle->SetOptFit(1);
}

DoublGaussFitter::~DoublGaussFitter()
{
  CleanupPreviousResults();
  if (fCanvas) {
    delete fCanvas;
    fCanvas = nullptr;
  }
}

Bool_t DoublGaussFitter::FitHistogram(TH1 *hist, const TString &setTitle)
{
  if (!hist) {
    std::cerr << "[DoublGaussFitter] Histogram is null!" << std::endl;
    return false;
  }

  if (hist->GetEntries() == 0) {
    std::cerr << "[DoublGaussFitter] Histogram is empty!" << std::endl;
    return false;
  }

  FitParams params = PrepareDefaultParams(hist);
  return FitHistogram(hist, params, setTitle);
}

Bool_t DoublGaussFitter::FitHistogram(TH1 *hist, const FitParams &params, const TString &setTitle)
{
  if (!hist) {
    std::cerr << "[DoublGaussFitter] Histogram is null!" << std::endl;
    return false;
  }

  CleanupPreviousResults();

  fLastHistogram = hist;
  fLastTitle = setTitle.Length() == 0 ? TString(hist->GetTitle()) : setTitle;

  if (fVerbose) {
    std::cout << "[DoublGaussFitter] Fitting histogram: " << hist->GetName() << std::endl;
  }

  return DoFit(hist, params);
}

DoublGaussFitter::FitParams DoublGaussFitter::PrepareDefaultParams(TH1 *hist)
{
  FitParams params;

  if (!hist || hist->GetEntries() == 0) {
    std::cerr << "[DoublGaussFitter] Cannot prepare params for null/empty histogram" << std::endl;
    return params;
  }

  Double_t xmin = hist->GetXaxis()->GetXmin();
  Double_t xmax = hist->GetXaxis()->GetXmax();
  Double_t mean = hist->GetMean();
  Double_t rms = hist->GetRMS();
  Double_t integral = hist->Integral();

  params.SetFitRange(xmin, xmax);

  // Parametry początkowe: [A1, sigma1, A2, sigma2, mean]
  // A1, A2: amplitudy (60% i 40%)
  // sigma1, sigma2: szerokości (różne - jeden węższy, jeden szerszy)
  // mean: wspólna średnia
  std::vector<Double_t> initial(5);
  
  initial[0] = integral * 0.6;      // A1 - główna amplituda
  initial[1] = rms * 0.3;           // sigma1 - węższy gauss (core)
  initial[2] = integral * 0.4;      // A2 - druga amplituda
  initial[3] = rms * 1.4;           // sigma2 - szerszy gauss (tail)
  initial[4] = mean;                // mean - wspólna średnia

  params.SetInitialValues(initial);

  // Granice parametrów
  std::vector<Double_t> lower(5), upper(5);
  
  lower[0] = 0.0;       // A1 dolna granica
  upper[0] = integral * 1000.0;        // A1 górna granica
  
  lower[1] = rms * 0.001;             // sigma1 dolna granica
  upper[1] = rms * 2.0;             // sigma1 górna granica (dla core)
  
  lower[2] = 0.0;       // A2 dolna granica
  upper[2] = integral * 2.0;       // A2 górna granica
  
  lower[3] = rms * 0.2;             // sigma2 dolna granica
  upper[3] = rms * 5.0;             // sigma2 górna granica (dla tail)
  
  lower[4] = mean - rms * 0.2;      // mean dolna granica
  upper[4] = mean + rms * 0.2;      // mean górna granica

  params.SetBounds(lower, upper);

  if (fVerbose) {
    std::cout << "[DoublGaussFitter] Default params prepared" << std::endl;
    std::cout << "  Mean: " << mean << ", RMS: " << rms << ", Integral: " << integral << std::endl;
  }

  return params;
}

Bool_t DoublGaussFitter::DoFit(TH1 *hist, const FitParams &params)
{
  Double_t fitMin = params.fitRangeMin > -999 ? params.fitRangeMin : hist->GetXaxis()->GetXmin();
  Double_t fitMax = params.fitRangeMax > -999 ? params.fitRangeMax : hist->GetXaxis()->GetXmax();

  TString funcName = Form("doubleGaus_%p_%ld", (void *)hist, (long)time(nullptr));
  TF1 *fitFunc = new TF1(funcName, double_gaus, fitMin, fitMax, 5);

  // Ustaw nazwy parametrów
  TString paramNames[5] = {"A1", "σ1", "A2", "σ2", "mean"};
  for (Int_t i = 0; i < 5; ++i) {
    fitFunc->SetParName(i, paramNames[i]);
    fitFunc->SetParameter(i, params.initialValues[i]);
    fitFunc->SetParLimits(i, params.lowerBounds[i], params.upperBounds[i]);
  }

  fitFunc->SetLineColor(fMainColor);
  fitFunc->SetLineStyle(fMainLineStyle);
  fitFunc->SetLineWidth(fMainLineWidth);

  if (fVerbose) {
    std::cout << "[DoublGaussFitter] Starting fit..." << std::endl;
  }

  // Wykonaj fit
  Int_t fitStatus = hist->Fit(fitFunc, params.fitOptions, "");

  // Sprawdź wyniki
  fLastResult.status = fitStatus;
  fLastResult.converged = (fitStatus == 0 && fitFunc->GetChisquare() > 0);

  if (fLastResult.converged) {
    fLastResult.chi2 = fitFunc->GetChisquare();
    fLastResult.ndf = fitFunc->GetNDF();
    fLastResult.probability = fitFunc->GetProb();

    // Pobierz parametry i błędy
    for (Int_t i = 0; i < 5; ++i) {
      fLastResult.parameters[i] = fitFunc->GetParameter(i);
      fLastResult.errors[i] = fitFunc->GetParError(i);
    }

    // Oblicz wartości kombinowane
    CalculateCombinedValues();

    // Stwórz TF1 dla rysowania
    fLastResult.fitFunction = (TF1*)fitFunc->Clone();

    if (fVerbose) {
      std::cout << "[DoublGaussFitter] Fit successful!" << std::endl;
      std::cout << "  Chi2/NDF: " << GetChi2NDF() << std::endl;
      std::cout << "  Core Sigma: " << fLastResult.coreSigma << " ± " << fLastResult.coreSigmaErr << std::endl;
    }

    delete fitFunc;
    return true;
  } else {
    std::cerr << "[DoublGaussFitter] Fit failed with status: " << fitStatus << std::endl;
    std::cerr << "  Chi2: " << fitFunc->GetChisquare() << std::endl;
    delete fitFunc;
    return false;
  }
}

void DoublGaussFitter::CalculateCombinedValues()
{
  if (!fLastResult.converged)
    return;

  // Użyj funkcji z double_gaus.h
  fLastResult.mean = comb_mean_double(fLastResult.parameters.data(), fLastResult.errors.data());
  fLastResult.meanErr = comb_mean_err_double(fLastResult.parameters.data(), fLastResult.errors.data());
  fLastResult.combinedSigma = comb_std_dev_double(fLastResult.parameters.data(), fLastResult.errors.data());
  fLastResult.combinedSigmaErr = comb_std_dev_err_double(fLastResult.parameters.data(), fLastResult.errors.data());
  
  // Pobierz core sigma (węższy rozkład)
  fLastResult.coreSigma = get_core_sigma(fLastResult.parameters.data());
  fLastResult.coreSigmaErr = get_core_sigma_err(fLastResult.parameters.data(), fLastResult.errors.data());
}

void DoublGaussFitter::PrintResults(Bool_t detailed) const
{
  if (!fLastResult.converged) {
    std::cout << "Last fit did not converge!" << std::endl;
    return;
  }

  std::cout << "\n==================== DOUBLE GAUSSIAN FIT RESULTS ====================" << std::endl;
  std::cout << "Title: " << fLastTitle.Data() << std::endl;
  std::cout << "Histogram: " << (fLastHistogram ? fLastHistogram->GetName() : "Unknown") << std::endl;
  std::cout << "Status: " << fLastResult.status << " (0 = success)" << std::endl;
  std::cout << "Chi2/NDF: " << fLastResult.chi2 << "/" << fLastResult.ndf;
  if (fLastResult.ndf > 0) {
    std::cout << " = " << GetChi2NDF();
  }
  std::cout << std::endl;
  std::cout << "Probability: " << fLastResult.probability << std::endl;

  if (detailed) {
    std::cout << "\nIndividual Parameters:" << std::endl;
    std::cout << "  A1:      " << fLastResult.parameters[0] << " ± " << fLastResult.errors[0] << std::endl;
    std::cout << "  Sigma1:  " << fLastResult.parameters[1] << " ± " << fLastResult.errors[1] << std::endl;
    std::cout << "  A2:      " << fLastResult.parameters[2] << " ± " << fLastResult.errors[2] << std::endl;
    std::cout << "  Sigma2:  " << fLastResult.parameters[3] << " ± " << fLastResult.errors[3] << std::endl;
    std::cout << "  Mean:    " << fLastResult.parameters[4] << " ± " << fLastResult.errors[4] << std::endl;
  }

  std::cout << "\nCombined Results:" << std::endl;
  std::cout << "  Mean:            " << fLastResult.mean << " ± " << fLastResult.meanErr << std::endl;
  std::cout << "  Combined Sigma:  " << fLastResult.combinedSigma << " ± " << fLastResult.combinedSigmaErr << std::endl;
  std::cout << "  Core Sigma:      " << fLastResult.coreSigma << " ± " << fLastResult.coreSigmaErr << std::endl;
  std::cout << "=======================================================================" << std::endl;
}

void DoublGaussFitter::PrintSummary() const
{
  if (!fLastResult.converged) {
    std::cout << "Last fit did not converge!" << std::endl;
    return;
  }

  std::cout << "Double Gauss Fit: Chi2/NDF=" << GetChi2NDF()
            << ", Mean=" << fLastResult.mean << "±" << fLastResult.meanErr
            << ", Sigma=" << fLastResult.combinedSigma << "±" << fLastResult.combinedSigmaErr
            << ", CoreSigma=" << fLastResult.coreSigma << "±" << fLastResult.coreSigmaErr << std::endl;
}

TString DoublGaussFitter::GetSummaryText() const
{
  if (!fLastResult.converged) {
    return "Fit not converged!";
  }

  return Form("Double Gauss: #chi^{2}/NDF=%.2f, #mu=%.3f#pm%.3f, #sigma=%.3f#pm%.3f, #sigma_{core}=%.3f#pm%.3f",
              GetChi2NDF(), fLastResult.mean, fLastResult.meanErr,
              fLastResult.combinedSigma, fLastResult.combinedSigmaErr,
              fLastResult.coreSigma, fLastResult.coreSigmaErr);
}

TCanvas *DoublGaussFitter::DrawFit(const TString &canvasName, const TString &canvasTitle,
                                    Int_t width, Int_t height,
                                    Bool_t drawComponents, Bool_t showStats, Bool_t logScale)
{
  if (!fLastResult.converged || !fLastHistogram) {
    std::cerr << "[DoublGaussFitter] Cannot draw: fit not converged or histogram null" << std::endl;
    return nullptr;
  }

  TString cName = canvasName.Length() == 0 ? TString(Form("c_doubleGauss_%p", (void *)fLastHistogram)) : canvasName;
  TString cTitle = canvasTitle.Length() == 0 ? TString(Form("Double Gaussian Fit - %s", fLastTitle.Data())) : canvasTitle;

  if (fCanvas) {
    delete fCanvas;
  }
  fCanvas = new TCanvas(cName, cTitle, width, height);
  fCanvas->SetLeftMargin(0.12);
  fCanvas->SetBottomMargin(0.12);
  fCanvas->SetRightMargin(showStats ? 0.3 : 0.05);
  fCanvas->SetTopMargin(0.08);

  if (logScale) {
    fCanvas->SetLogy();
  }

  fLastHistogram->Draw("HIST");
  SetOptimalYAxisRange(fLastHistogram, logScale);

  DrawFitOnCurrentPad(drawComponents, showStats);

  fCanvas->Update();
  return fCanvas;
}

void DoublGaussFitter::DrawFitOnCurrentPad(Bool_t drawComponents, Bool_t showStats)
{
  if (!fLastResult.converged || !fLastResult.fitFunction) {
    std::cerr << "[DoublGaussFitter] Cannot draw: fit not converged" << std::endl;
    return;
  }

  fLastResult.fitFunction->SetLineColor(fMainColor);
  fLastResult.fitFunction->SetLineStyle(fMainLineStyle);
  fLastResult.fitFunction->SetLineWidth(fMainLineWidth);
  fLastResult.fitFunction->Draw("HIST SAME");

  if (drawComponents) {
    std::vector<TF1*> components = CreateComponentFunctions();
    for (size_t i = 0; i < components.size(); ++i) {
      if (components[i]) {
        Color_t colors[] = {fComp1Color, fComp2Color};
        components[i]->SetLineColor(colors[i]);
        components[i]->SetLineStyle(fComponentLineStyle);
        components[i]->SetLineWidth(fComponentLineWidth);
        components[i]->Draw("SAME");
      }
    }
  }

  if (showStats) {
    TPaveText *statsBox = CreateStatsBox();
    statsBox->Draw();
  }

  gPad->Update();
}

std::vector<TF1 *> DoublGaussFitter::CreateComponentFunctions()
{
  std::vector<TF1 *> components;

  if (!fLastResult.converged || !fLastResult.fitFunction) {
    return components;
  }

  Double_t xmin = fLastResult.fitFunction->GetXmin();
  Double_t xmax = fLastResult.fitFunction->GetXmax();

  for (Int_t i = 0; i < 2; ++i) {
    TString compName = Form("comp%d_%p", i, (void *)fLastHistogram);
    TF1 *comp = new TF1(compName, single_gaus_core, xmin, xmax, 3);
    
    comp->SetParameter(0, fLastResult.parameters[i*2]);     // Amplitude
    comp->SetParameter(1, fLastResult.parameters[i*2 + 1]); // Sigma
    comp->SetParameter(2, fLastResult.parameters[4]);       // Mean (wspólna)
    
    components.push_back(comp);
  }

  return components;
}

TPaveText *DoublGaussFitter::CreateStatsBox()
{
  TPaveText *stats = new TPaveText(0.6, 0.4, 0.9, 0.9, "NDC");
  stats->SetFillColor(kWhite);
  stats->SetBorderSize(1);
  stats->SetTextAlign(12);
  stats->SetTextSize(0.025);

  stats->AddText("Double Gaussian Fit");
  stats->AddText("");
  stats->AddText("Combined Results:");
  stats->AddText(Form("#mu = %.3f #pm %.3f", fLastResult.mean, fLastResult.meanErr));
  stats->AddText(Form("#sigma = %.3f #pm %.3f", fLastResult.combinedSigma, fLastResult.combinedSigmaErr));
  stats->AddText(Form("#sigma_{core} = %.3f #pm %.3f", fLastResult.coreSigma, fLastResult.coreSigmaErr));
  stats->AddText("");

  // Krótko o komponentach
  stats->AddText(Form("Gauss 1: #sigma=%.3f", fLastResult.parameters[1]));
  stats->AddText(Form("Gauss 2: #sigma=%.3f", fLastResult.parameters[3]));

  return stats;
}

Bool_t DoublGaussFitter::SaveFitToFile(const TString &filename)
{
  if (!fCanvas) {
    std::cerr << "[DoublGaussFitter] No canvas to save!" << std::endl;
    return false;
  }

  fCanvas->SaveAs(filename);
  if (fVerbose) {
    std::cout << "[DoublGaussFitter] Canvas saved to: " << filename << std::endl;
  }
  return true;
}

void DoublGaussFitter::SetComponentColors(Color_t mainColor, Color_t comp1Color, Color_t comp2Color)
{
  fMainColor = mainColor;
  fComp1Color = comp1Color;
  fComp2Color = comp2Color;
}

void DoublGaussFitter::SetLineStyles(Style_t mainStyle, Style_t componentStyle)
{
  fMainLineStyle = mainStyle;
  fComponentLineStyle = componentStyle;
}

void DoublGaussFitter::SetLineWidths(Width_t mainWidth, Width_t componentWidth)
{
  fMainLineWidth = mainWidth;
  fComponentLineWidth = componentWidth;
}

void DoublGaussFitter::CleanupPreviousResults()
{
  fLastResult = FitResult();
  fLastHistogram = nullptr;
  fLastTitle = "";
}

void DoublGaussFitter::SetOptimalYAxisRange(TH1 *hist, Bool_t isLogScale)
{
  if (!hist)
    return;

  Double_t maxValue = hist->GetMaximum();
  Double_t minValue = hist->GetMinimum();

  if (maxValue <= 0)
    return;

  if (isLogScale) {
    // Dla skali logarytmicznej: od 0.5 do 10x max
    hist->SetMinimum(0.5);
    hist->SetMaximum(maxValue * 10.0);
  } else {
    // Dla skali liniowej: od -5% do 120% max
    Double_t range = maxValue - minValue;
    hist->SetMinimum(minValue - 0.05 * range);
    hist->SetMaximum(maxValue * 1.2);
  }
}
