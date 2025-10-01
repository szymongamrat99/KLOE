/**
 * @file TripleGaussUtils.h
 * @brief Funkcje pomocnicze dla łatwego użycia KLOE::TripleGaussFitter z HistManager
 */

#pragma once

#include "TripleGaussFitter.h"
#include "HistManager.h"
#include <vector>
#include <TString.h>

namespace KLOE {

/**
 * @brief Fituje potrójnego Gaussa do wybranego histogramu i kanału
 * @param histMgr Wskaźnik na HistManager
 * @param setName Nazwa zestawu histogramów
 * @param mctruth Numer kanału MC (0=suma, 1-N=kanały, -1=dane)
 * @param verbose Czy wypisać szczegółowe informacje
 * @return true jeśli fit się udał
 */
inline Bool_t FitTripleGaussToChannel(HistManager* histMgr, 
                                     const TString& setName, 
                                     Int_t mctruth, 
                                     Bool_t verbose = true) {
    if(!histMgr) {
        std::cerr << "ERROR: HistManager is null" << std::endl;
        return false;
    }
    
    // Pobierz histogram
    TH1* hist = histMgr->GetHistogram1D(setName, mctruth);
    if(!hist) {
        if(verbose) {
            std::cerr << "WARNING: No histogram found for set '" << setName.Data() 
                      << "' mctruth=" << mctruth << std::endl;
        }
        return false;
    }
    
    if(hist->GetEntries() < 50) {
        if(verbose) {
            std::cerr << "WARNING: Too few entries in histogram (" 
                      << hist->GetEntries() << ")" << std::endl;
        }
        return false;
    }
    
    // Stwórz fitter
    TripleGaussFitter fitter;
    fitter.SetVerbose(verbose);
    fitter.SetComponentColors(kRed, kBlue, kGreen, kMagenta);
    
    // Sprawdź czy histogram ma skalę logarytmiczną
    Bool_t isLogScale = histMgr->IsLogScale(setName);
    
    // Przygotuj tytuł
    TString title;
    if(mctruth == -1) {
        title = Form("%s Data", setName.Data());
    } else if(mctruth == 0) {
        title = Form("%s Sum MC", setName.Data());
    } else {
        title = Form("%s Channel %d", setName.Data(), mctruth);
    }
    
    // Wykonaj fit
    Bool_t success = fitter.FitHistogram(hist, title);
    
    if(success) {
        if(verbose) {
            fitter.PrintResults(true);
        }
        
        // Ustaw odpowiednią skalę i zakres Y na canvasie fitu
        TCanvas* canvas = fitter.DrawFit("", "", 1000, 700, true, true, isLogScale);
        if(canvas) {
            fitter.SetOptimalYAxisRange(hist, isLogScale);
            canvas->Update();
            
                // Zapisz plot
            TString plotName;
            if(mctruth == -1) {
                plotName = Form("img/%s_data_triple_gauss.png", setName.Data());
            } else if(mctruth == 0) {
                plotName = Form("img/%s_sum_triple_gauss.png", setName.Data());
            } else {
                plotName = Form("img/%s_ch%d_triple_gauss.png", setName.Data(), mctruth);
            }
            
            canvas->SaveAs(plotName);
            if(verbose) {
                std::cout << "Plot saved to: " << plotName.Data() << std::endl;
            }
        }
    }
    
    return success;
}

/**
 * @brief Fituje potrójnego Gaussa do wszystkich kanałów w zestawie
 * @param histMgr Wskaźnik na HistManager  
 * @param setName Nazwa zestawu histogramów
 * @param fitSum Czy fitować sumę MC (mctruth=0)
 * @param fitChannels Czy fitować poszczególne kanały MC
 * @param fitData Czy fitować dane (mctruth=-1)
 * @param minEntries Minimalna liczba wpisów do fitu
 * @param verbose Czy wypisać szczegółowe informacje
 * @return Liczba udanych fitów
 */
inline Int_t FitTripleGaussToAllChannels(HistManager* histMgr,
                                        const TString& setName,
                                        Bool_t fitSum = true,
                                        Bool_t fitChannels = true, 
                                        Bool_t fitData = true,
                                        Int_t minEntries = 100,
                                        Bool_t verbose = true) {
    if(!histMgr) {
        std::cerr << "ERROR: HistManager is null" << std::endl;
        return 0;
    }
    
    Int_t successfulFits = 0;
    
    if(verbose) {
        std::cout << "\n=== Triple Gaussian Fitting for " << setName.Data() << " ===" << std::endl;
    }
    
    // Fituj sumę MC
    if(fitSum) {
        TH1* sumHist = histMgr->GetHistogram1D(setName, 0);
        if(sumHist && sumHist->GetEntries() >= minEntries) {
            if(verbose) {
                std::cout << "\nFitting Sum MC (entries: " << sumHist->GetEntries() << ")..." << std::endl;
            }
            if(FitTripleGaussToChannel(histMgr, setName, 0, verbose)) {
                successfulFits++;
            }
        }
    }
    
    // Fituj poszczególne kanały MC
    if(fitChannels) {
        // Sprawdź ile kanałów MC mamy
        std::vector<TH1*> allHists = histMgr->GetAllHistograms1D(setName);
        Int_t numChannels = allHists.size() > 0 ? allHists.size() - 1 : 0; // -1 bo pierwszy to suma
        
        for(Int_t ch = 1; ch <= numChannels; ++ch) {
            TH1* channelHist = histMgr->GetHistogram1D(setName, ch);
            if(channelHist && channelHist->GetEntries() >= minEntries) {
                if(verbose) {
                    std::cout << "\nFitting Channel " << ch << " (entries: " 
                              << channelHist->GetEntries() << ")..." << std::endl;
                }
                if(FitTripleGaussToChannel(histMgr, setName, ch, verbose)) {
                    successfulFits++;
                }
            }
        }
    }
    
    // Fituj dane
    if(fitData) {
        TH1* dataHist = histMgr->GetDataHistogram1D(setName);
        if(dataHist && dataHist->GetEntries() >= minEntries) {
            if(verbose) {
                std::cout << "\nFitting Data (entries: " << dataHist->GetEntries() << ")..." << std::endl;
            }
            if(FitTripleGaussToChannel(histMgr, setName, -1, verbose)) {
                successfulFits++;
            }
        }
    }
    
    if(verbose) {
        std::cout << "\n=== Triple Gaussian Fitting Summary ===" << std::endl;
        std::cout << "Set: " << setName.Data() << std::endl;
        std::cout << "Successful fits: " << successfulFits << std::endl;
        std::cout << "========================================" << std::endl;
    }
    
    return successfulFits;
}

/**
 * @brief Porównuje wyniki Triple Gaussian fitów na jednym canvasie
 * @param histMgr Wskaźnik na HistManager
 * @param setName Nazwa zestawu histogramów
 * @param canvasName Nazwa canvasu (opcjonalne)
 * @param canvasTitle Tytuł canvasu (opcjonalne)
 * @return Wskaźnik na canvas z porównaniem
 */
inline TCanvas* CompareTripleGaussFits(HistManager* histMgr,
                                      const TString& setName,
                                      const TString& canvasName = "",
                                      const TString& canvasTitle = "") {
    if(!histMgr) {
        std::cerr << "ERROR: HistManager is null" << std::endl;
        return nullptr;
    }
    
    // Sprawdź dostępne histogramy
    std::vector<TH1*> allHists = histMgr->GetAllHistograms1D(setName);
    TH1* dataHist = histMgr->GetDataHistogram1D(setName);
    
    Int_t numPads = 0;
    if(allHists.size() > 0 && allHists[0] && allHists[0]->GetEntries() > 50) numPads++; // Suma
    if(dataHist && dataHist->GetEntries() > 50) numPads++; // Dane
    
    // Dodaj kanały z wystarczającą statystyką
    Int_t numChannels = 0;
    for(size_t i = 1; i < allHists.size(); ++i) {
        if(allHists[i] && allHists[i]->GetEntries() > 50) {
            numChannels++;
        }
    }
    numPads += numChannels;
    
    if(numPads == 0) {
        std::cerr << "ERROR: No histograms with sufficient statistics found" << std::endl;
        return nullptr;
    }
    
    // Określ układ canvasu
    Int_t nx = (numPads == 1) ? 1 : (numPads <= 4) ? 2 : 3;
    Int_t ny = (numPads + nx - 1) / nx;
    
    // Stwórz canvas
    TString cName = canvasName.Length() == 0 ? TString(Form("c_%s_comparison", setName.Data())) : canvasName;
    TString cTitle = canvasTitle.Length() == 0 ? TString(Form("Triple Gaussian Fits - %s", setName.Data())) : canvasTitle;
    
    TCanvas* canvas = new TCanvas(cName, cTitle, nx * 400, ny * 300);
    canvas->Divide(nx, ny);
    
    TripleGaussFitter fitter;
    fitter.SetVerbose(false);
    fitter.SetComponentColors(kRed, kBlue, kGreen, kMagenta);
    
    // Sprawdź czy histogram ma skalę logarytmiczną
    Bool_t isLogScale = histMgr->IsLogScale(setName);
    
    Int_t pad = 1;
    
    // Rysuj sumę MC
    if(allHists.size() > 0 && allHists[0] && allHists[0]->GetEntries() > 50) {
        canvas->cd(pad++);
        if(isLogScale) {
            gPad->SetLogy(1);
        }
        allHists[0]->SetTitle("Sum MC");
        
        // Ustaw optymalny zakres osi Y
        Double_t maxVal = allHists[0]->GetMaximum();
        if(isLogScale) {
            allHists[0]->SetMinimum(0.1);
            allHists[0]->SetMaximum(maxVal * 10.0);
        } else {
            allHists[0]->SetMinimum(0.0);
            allHists[0]->SetMaximum(maxVal * 1.4);
        }
        
        allHists[0]->Draw("HIST");
        
        if(fitter.FitHistogram(allHists[0], "Sum MC")) {
            fitter.DrawFitOnCurrentPad(true, false);
        }
    }
    
    // Rysuj kanały MC
    for(size_t i = 1; i < allHists.size(); ++i) {
        if(allHists[i] && allHists[i]->GetEntries() > 50) {
            canvas->cd(pad++);
            if(isLogScale) {
                gPad->SetLogy(1);
            }
            allHists[i]->SetTitle(Form("Channel %d", (Int_t)i));
            
            // Ustaw optymalny zakres osi Y
            Double_t maxVal = allHists[i]->GetMaximum();
            if(isLogScale) {
                allHists[i]->SetMinimum(0.1);
                allHists[i]->SetMaximum(maxVal * 10.0);
            } else {
                allHists[i]->SetMinimum(0.0);
                allHists[i]->SetMaximum(maxVal * 1.4);
            }
            
            allHists[i]->Draw("HIST");
            
            if(fitter.FitHistogram(allHists[i], Form("Channel %d", (Int_t)i))) {
                fitter.DrawFitOnCurrentPad(true, false);
            }
        }
    }
    
    // Rysuj dane
    if(dataHist && dataHist->GetEntries() > 50) {
        canvas->cd(pad++);
        if(isLogScale) {
            gPad->SetLogy(1);
        }
        dataHist->SetTitle("Data");
        
        // Ustaw optymalny zakres osi Y
        Double_t maxVal = dataHist->GetMaximum();
        if(isLogScale) {
            dataHist->SetMinimum(0.1);
            dataHist->SetMaximum(maxVal * 10.0);
        } else {
            dataHist->SetMinimum(0.0);
            dataHist->SetMaximum(maxVal * 1.4);
        }
        
        dataHist->Draw("HIST");
        
        if(fitter.FitHistogram(dataHist, "Data")) {
            fitter.DrawFitOnCurrentPad(true, false);
        }
    }
    
    canvas->Update();
    return canvas;
}

/**
 * @brief Funkcja pomocnicza do szybkiego fitu wybranego kanału z terminala
 * Przykład użycia: FitChannelQuick(histMgr, "timeDiff", 0)  // fit suma MC
 */
inline void FitChannelQuick(HistManager* histMgr, const TString& setName, Int_t mctruth) {
    std::cout << "\n=== Quick Triple Gaussian Fit ===" << std::endl;
    std::cout << "Set: " << setName.Data() << ", mctruth: " << mctruth << std::endl;
    
    Bool_t success = FitTripleGaussToChannel(histMgr, setName, mctruth, true);
    
    if(success) {
        std::cout << "✓ Fit completed successfully!" << std::endl;
    } else {
        std::cout << "✗ Fit failed!" << std::endl;
    }
    std::cout << "=================================" << std::endl;
}

} // namespace KLOE
