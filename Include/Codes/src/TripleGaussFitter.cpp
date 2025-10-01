#include "TripleGaussFitter.h"
#include <TStyle.h>
#include <TRandom.h>
#include <algorithm>

using namespace KLOE;

TripleGaussFitter::TripleGaussFitter() {
    // Inicjalizacja domyślnych ustawień
    gStyle->SetOptFit(1);
}

TripleGaussFitter::~TripleGaussFitter() {
    CleanupPreviousResults();
    if(fCanvas) {
        delete fCanvas;
    }
}

Bool_t TripleGaussFitter::FitHistogram(TH1* hist, const TString& setTitle) {
    if(!hist) {
        std::cerr << "ERROR: TripleGaussFitter::FitHistogram - null histogram" << std::endl;
        return false;
    }
    
    if(hist->GetEntries() == 0) {
        std::cerr << "ERROR: TripleGaussFitter::FitHistogram - empty histogram" << std::endl;
        return false;
    }
    
    // Przygotuj domyślne parametry
    FitParams params = PrepareDefaultParams(hist);
    
    // Wykonaj fit
    return FitHistogram(hist, params, setTitle);
}

Bool_t TripleGaussFitter::FitHistogram(TH1* hist, const FitParams& params, const TString& setTitle) {
    if(!hist) {
        std::cerr << "ERROR: TripleGaussFitter::FitHistogram - null histogram" << std::endl;
        return false;
    }
    
    CleanupPreviousResults();
    
    fLastHistogram = hist;
    fLastTitle = setTitle.Length() == 0 ? TString(hist->GetTitle()) : setTitle;
    
    if(fVerbose) {
        std::cout << "INFO: Starting Triple Gaussian fit for: " << fLastTitle.Data() << std::endl;
    }
    
    return DoFit(hist, params);
}

TripleGaussFitter::FitParams TripleGaussFitter::PrepareDefaultParams(TH1* hist) {
    FitParams params;
    
    if(!hist || hist->GetEntries() == 0) {
        std::cerr << "WARNING: Cannot prepare params for null/empty histogram" << std::endl;
        return params;
    }
    
    // Podstawowe właściwości histogramu
    Double_t xmin = hist->GetXaxis()->GetXmin();
    Double_t xmax = hist->GetXaxis()->GetXmax();
    Double_t range = xmax - xmin;
    Double_t mean = hist->GetMean();
    Double_t rms = hist->GetRMS();
    Double_t integral = hist->Integral();
    Double_t maxValue = hist->GetMaximum();
    
    // Ustaw zakres fitu na cały histogram
    params.SetFitRange(xmin, xmax);
    
    // Parametry początkowe dla trzech gaussów
    std::vector<Double_t> initial(9);
    
    // Gauss 1 (główny) - 60% amplitudy, wyśrodkowany
    initial[0] = integral * 0.6;    // A1
    initial[1] = mean;              // μ1
    initial[2] = rms * 0.8;         // σ1
    
    // Gauss 2 (lewy) - 25% amplitudy, przesunięty w lewo
    initial[3] = integral * 0.25;   // A2
    initial[4] = mean - rms * 0.7;  // μ2
    initial[5] = rms * 1.2;         // σ2
    
    // Gauss 3 (prawy) - 15% amplitudy, przesunięty w prawo
    initial[6] = integral * 0.15;   // A3
    initial[7] = mean + rms * 0.7;  // μ3
    initial[8] = rms * 1.5;         // σ3
    
    params.SetInitialValues(initial);
    
    // Granice parametrów
    std::vector<Double_t> lower(9), upper(9);
    
    for(Int_t i = 0; i < 3; ++i) {
        Int_t idx = i * 3;
        // Amplitudy - od 0 do 2x całkowita wartość
        lower[idx] = 0.0;
        upper[idx] = integral * 2.0;
        
        // Średnie - w zakresie histogramu ± 20%
        lower[idx + 1] = xmin - range * 0.2;
        upper[idx + 1] = xmax + range * 0.2;
        
        // Szerokości - od bardzo małej do całego zakresu
        lower[idx + 2] = range * 0.01;
        upper[idx + 2] = range;
    }
    
    params.SetBounds(lower, upper);
    
    if(fVerbose) {
        std::cout << "INFO: Default parameters prepared for histogram with:" << std::endl;
        std::cout << "  Mean: " << mean << ", RMS: " << rms << ", Integral: " << integral << std::endl;
    }
    
    return params;
}

Bool_t TripleGaussFitter::DoFit(TH1* hist, const FitParams& params) {
    // Przygotuj zakres fitu
    Double_t fitMin = params.fitRangeMin > -999 ? params.fitRangeMin : hist->GetXaxis()->GetXmin();
    Double_t fitMax = params.fitRangeMax > -999 ? params.fitRangeMax : hist->GetXaxis()->GetXmax();
    
    // Stwórz funkcję fitu
    TString funcName = Form("tripleGaus_%p_%ld", (void*)hist, (long)time(nullptr));
    TF1* fitFunc = new TF1(funcName, triple_gaus, fitMin, fitMax, 9);
    
    // Ustaw nazwy parametrów
    TString paramNames[9] = {"A1", "μ1", "σ1", "A2", "μ2", "σ2", "A3", "μ3", "σ3"};
    for(Int_t i = 0; i < 9; ++i) {
        fitFunc->SetParName(i, paramNames[i]);
        fitFunc->SetParameter(i, params.initialValues[i]);
        fitFunc->SetParLimits(i, params.lowerBounds[i], params.upperBounds[i]);
    }
    
    // Ustawienia funkcji
    fitFunc->SetLineColor(fMainColor);
    fitFunc->SetLineStyle(fMainLineStyle);
    fitFunc->SetLineWidth(fMainLineWidth);
    
    if(fVerbose) {
        std::cout << "INFO: Performing fit in range [" << fitMin << ", " << fitMax << "]" << std::endl;
    }
    
    // Wykonaj fit
    Int_t fitStatus = hist->Fit(fitFunc, params.fitOptions, "");
    
    // Sprawdź wyniki
    fLastResult.status = fitStatus;
    fLastResult.converged = (fitStatus == 0 && fitFunc->GetChisquare() > 0);
    
    if(fLastResult.converged) {
        // Pobierz parametry
        for(Int_t i = 0; i < 9; ++i) {
            fLastResult.parameters[i] = fitFunc->GetParameter(i);
            fLastResult.errors[i] = fitFunc->GetParError(i);
        }
        
        fLastResult.chi2 = fitFunc->GetChisquare();
        fLastResult.ndf = fitFunc->GetNDF();
        fLastResult.probability = TMath::Prob(fLastResult.chi2, fLastResult.ndf);
        
        // Oblicz wartości kombinowane
        CalculateCombinedValues();
        
        // Przechowaj funkcję
        fLastResult.fitFunction = (TF1*)fitFunc->Clone();
        
        if(fVerbose) {
            std::cout << "INFO: Fit successful! Chi2/NDF = " << fLastResult.chi2 
                      << "/" << fLastResult.ndf << " = " << GetChi2NDF() << std::endl;
        }
        
        delete fitFunc;
        return true;
    } else {
        if(fVerbose) {
            std::cerr << "WARNING: Fit failed or did not converge (status = " << fitStatus << ")" << std::endl;
        }
        delete fitFunc;
        return false;
    }
}

void TripleGaussFitter::CalculateCombinedValues() {
    if(!fLastResult.converged) return;
    
    // Użyj funkcji z triple_gaus.h
    fLastResult.combinedMean = comb_mean(fLastResult.parameters.data(), fLastResult.errors.data());
    fLastResult.combinedMeanErr = comb_mean_err(fLastResult.parameters.data(), fLastResult.errors.data());
    fLastResult.combinedSigma = comb_std_dev(fLastResult.parameters.data(), fLastResult.errors.data());
    fLastResult.combinedSigmaErr = comb_std_dev_err(fLastResult.parameters.data(), fLastResult.errors.data());
}

void TripleGaussFitter::PrintResults(Bool_t detailed) const {
    if(!fLastResult.converged) {
        std::cout << "No successful fit results to display." << std::endl;
        return;
    }
    
    std::cout << "\n==================== TRIPLE GAUSSIAN FIT RESULTS ====================" << std::endl;
    std::cout << "Title: " << fLastTitle.Data() << std::endl;
    std::cout << "Histogram: " << (fLastHistogram ? fLastHistogram->GetName() : "Unknown") << std::endl;
    std::cout << "Status: " << fLastResult.status << " (0 = success)" << std::endl;
    std::cout << "Chi2/NDF: " << fLastResult.chi2 << "/" << fLastResult.ndf;
    if(fLastResult.ndf > 0) {
        std::cout << " = " << GetChi2NDF();
    }
    std::cout << std::endl;
    std::cout << "Probability: " << fLastResult.probability << std::endl;
    
    if(detailed) {
        std::cout << "\nIndividual Gaussian Components:" << std::endl;
        for(Int_t i = 0; i < 3; ++i) {
            Int_t idx = i * 3;
            std::cout << "Gaussian " << (i+1) << ":" << std::endl;
            std::cout << "  Amplitude: " << fLastResult.parameters[idx] << " ± " << fLastResult.errors[idx] << std::endl;
            std::cout << "  Mean:      " << fLastResult.parameters[idx+1] << " ± " << fLastResult.errors[idx+1] << std::endl;
            std::cout << "  Sigma:     " << fLastResult.parameters[idx+2] << " ± " << fLastResult.errors[idx+2] << std::endl;
            std::cout << std::endl;
        }
    }
    
    std::cout << "Combined Results:" << std::endl;
    std::cout << "  Combined Mean:  " << fLastResult.combinedMean << " ± " << fLastResult.combinedMeanErr << std::endl;
    std::cout << "  Combined Sigma: " << fLastResult.combinedSigma << " ± " << fLastResult.combinedSigmaErr << std::endl;
    std::cout << "=======================================================================" << std::endl;
}

void TripleGaussFitter::PrintSummary() const {
    if(!fLastResult.converged) {
        std::cout << "Triple Gauss Fit: FAILED" << std::endl;
        return;
    }
    
    std::cout << "Triple Gauss Fit: Chi2/NDF=" << GetChi2NDF() 
              << ", Mean=" << fLastResult.combinedMean << "±" << fLastResult.combinedMeanErr
              << ", Sigma=" << fLastResult.combinedSigma << "±" << fLastResult.combinedSigmaErr << std::endl;
}

TString TripleGaussFitter::GetSummaryText() const {
    if(!fLastResult.converged) {
        return "Triple Gauss Fit: FAILED";
    }
    
    return Form("Triple Gauss: #chi^{2}/NDF=%.2f, #mu=%.3f#pm%.3f, #sigma=%.3f#pm%.3f",
                GetChi2NDF(), fLastResult.combinedMean, fLastResult.combinedMeanErr,
                fLastResult.combinedSigma, fLastResult.combinedSigmaErr);
}

TCanvas* TripleGaussFitter::DrawFit(const TString& canvasName, const TString& canvasTitle,
                                   Int_t width, Int_t height,
                                   Bool_t drawComponents, Bool_t showStats, Bool_t logScale) {
    if(!fLastResult.converged || !fLastHistogram) {
        std::cerr << "ERROR: No successful fit results to draw" << std::endl;
        return nullptr;
    }
    
    // Nazwa i tytuł canvasu
    TString cName = canvasName.Length() == 0 ? TString(Form("c_tripleGauss_%p", (void*)fLastHistogram)) : canvasName;
    TString cTitle = canvasTitle.Length() == 0 ? TString(Form("Triple Gaussian Fit - %s", fLastTitle.Data())) : canvasTitle;
    
    // Stwórz canvas
    if(fCanvas) {
        delete fCanvas;
    }
    fCanvas = new TCanvas(cName, cTitle, width, height);
    fCanvas->SetLeftMargin(0.12);
    fCanvas->SetBottomMargin(0.12);
    fCanvas->SetRightMargin(showStats ? 0.3 : 0.05);
    fCanvas->SetTopMargin(0.08);
    
    // Ustaw skalę logarytmiczną jeśli wymagana
    if(logScale) {
        fCanvas->SetLogy(1);
    }
    
    // Narysuj histogram
    fLastHistogram->Draw("HIST");
    
    // Ustaw optymalny zakres osi Y
    SetOptimalYAxisRange(fLastHistogram, logScale);
    
    // Narysuj fit na tym canvasie
    DrawFitOnCurrentPad(drawComponents, showStats);
    
    fCanvas->Update();
    return fCanvas;
}

void TripleGaussFitter::DrawFitOnCurrentPad(Bool_t drawComponents, Bool_t showStats) {
    if(!fLastResult.converged || !fLastResult.fitFunction) {
        std::cerr << "ERROR: No successful fit results to draw" << std::endl;
        return;
    }
    
    // Narysuj główną funkcję fitu
    fLastResult.fitFunction->SetLineColor(fMainColor);
    fLastResult.fitFunction->SetLineStyle(fMainLineStyle);
    fLastResult.fitFunction->SetLineWidth(fMainLineWidth);
    fLastResult.fitFunction->Draw("SAME");
    
    // Narysuj komponenty jeśli wymagane
    if(drawComponents) {
        std::vector<TF1*> components = CreateComponentFunctions();
        
        // Kolory i style dla komponentów
        Color_t colors[3] = {fComp1Color, fComp2Color, fComp3Color};
        
        for(Int_t i = 0; i < 3; ++i) {
            components[i]->SetLineColor(colors[i]);
            components[i]->SetLineStyle(fComponentLineStyle);
            components[i]->SetLineWidth(fComponentLineWidth);
            components[i]->Draw("SAME");
        }
    }
    
    // Dodaj statystyki jeśli wymagane
    if(showStats) {
        TPaveText* statsBox = CreateStatsBox();
        statsBox->Draw();
    }
    
    gPad->Update();
}

std::vector<TF1*> TripleGaussFitter::CreateComponentFunctions() {
    std::vector<TF1*> components;
    
    if(!fLastResult.converged || !fLastResult.fitFunction) {
        return components;
    }
    
    Double_t xmin = fLastResult.fitFunction->GetXmin();
    Double_t xmax = fLastResult.fitFunction->GetXmax();
    
    for(Int_t i = 0; i < 3; ++i) {
        TString compName = Form("comp%d_%p", i+1, (void*)fLastHistogram);
        TF1* comp = new TF1(compName, single_gaus, xmin, xmax, 3);
        
        Int_t idx = i * 3;
        comp->SetParameter(0, fLastResult.parameters[idx]);     // A
        comp->SetParameter(1, fLastResult.parameters[idx+1]);   // μ
        comp->SetParameter(2, fLastResult.parameters[idx+2]);   // σ
        
        components.push_back(comp);
    }
    
    return components;
}

TPaveText* TripleGaussFitter::CreateStatsBox() {
    TPaveText* stats = new TPaveText(0.72, 0.4, 0.98, 0.9, "NDC");
    stats->SetFillColor(kWhite);
    stats->SetBorderSize(1);
    stats->SetTextAlign(12);
    stats->SetTextSize(0.025);
    
    stats->AddText("Triple Gaussian Fit");
    stats->AddText("");
    stats->AddText(Form("#chi^{2}/NDF = %.2f/%d", fLastResult.chi2, fLastResult.ndf));
    stats->AddText(Form("= %.3f", GetChi2NDF()));
    stats->AddText(Form("Prob = %.4f", fLastResult.probability));
    stats->AddText("");
    stats->AddText("Combined Results:");
    stats->AddText(Form("#mu = %.3f #pm %.3f", fLastResult.combinedMean, fLastResult.combinedMeanErr));
    stats->AddText(Form("#sigma = %.3f #pm %.3f", fLastResult.combinedSigma, fLastResult.combinedSigmaErr));
    stats->AddText("");
    
    // Dodaj skrócone info o komponentach
    for(Int_t i = 0; i < 3; ++i) {
        Int_t idx = i * 3;
        stats->AddText(Form("G%d: A=%.1f, #mu=%.2f, #sigma=%.2f", 
                           i+1, fLastResult.parameters[idx], 
                           fLastResult.parameters[idx+1], fLastResult.parameters[idx+2]));
    }
    
    return stats;
}

Bool_t TripleGaussFitter::SaveFitToFile(const TString& filename) {
    if(!fCanvas) {
        std::cerr << "ERROR: No canvas to save" << std::endl;
        return false;
    }
    
    fCanvas->SaveAs(filename);
    if(fVerbose) {
        std::cout << "INFO: Fit canvas saved to " << filename.Data() << std::endl;
    }
    return true;
}

void TripleGaussFitter::SetComponentColors(Color_t mainColor, Color_t comp1Color, 
                                          Color_t comp2Color, Color_t comp3Color) {
    fMainColor = mainColor;
    fComp1Color = comp1Color;
    fComp2Color = comp2Color;
    fComp3Color = comp3Color;
}

void TripleGaussFitter::SetLineStyles(Style_t mainStyle, Style_t componentStyle) {
    fMainLineStyle = mainStyle;
    fComponentLineStyle = componentStyle;
}

void TripleGaussFitter::SetLineWidths(Width_t mainWidth, Width_t componentWidth) {
    fMainLineWidth = mainWidth;
    fComponentLineWidth = componentWidth;
}

void TripleGaussFitter::CleanupPreviousResults() {
    // FitResult destruktor sam posprzata fitFunction
    fLastResult = FitResult();
    fLastHistogram = nullptr;
    fLastTitle = "";
}

void TripleGaussFitter::SetOptimalYAxisRange(TH1* hist, Bool_t isLogScale) {
    if(!hist) return;
    
    // Znajdź maksymalną wartość w histogramie
    Double_t maxValue = hist->GetMaximum();
    Double_t minValue = hist->GetMinimum();
    
    if(maxValue <= 0) return;
    
    if(isLogScale) {
        // Dla skali logarytmicznej
        Double_t yMin = minValue > 0 ? minValue * 0.1 : 0.1;
        Double_t yMax = maxValue * 10.0;
        
        hist->SetMinimum(yMin);
        hist->SetMaximum(yMax);
        
        if(fVerbose) {
            std::cout << "INFO: Y-axis range set to [" << yMin << ", " << yMax << "] (log scale)" << std::endl;
        }
    } else {
        // Dla skali liniowej - +40% nad maksimum
        Double_t yMin = TMath::Min(0.0, minValue);
        Double_t yMax = maxValue * 1.4;
        
        hist->SetMinimum(yMin);
        hist->SetMaximum(yMax);
        
        if(fVerbose) {
            std::cout << "INFO: Y-axis range set to [" << yMin << ", " << yMax << "] (linear scale)" << std::endl;
        }
    }
}