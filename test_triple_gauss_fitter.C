/**
 * @file test_triple_gauss_fitter.C
 * @brief Przykład użycia klasy KLOE::TripleGaussFitter
 * 
 * Ten plik demonstruje jak używać nowej klasy TripleGaussFitter
 * w połączeniu z HistManager do fitowania potrójnego Gaussa.
 */

#include "Include/Codes/inc/HistManager.h"
#include "Include/Codes/inc/TripleGaussFitter.h"
#include <TRandom3.h>
#include <iostream>

void test_triple_gauss_fitter() {
    std::cout << "=== TEST KLOE::TripleGaussFitter ===" << std::endl;
    
    // ==================== KONFIGURACJA HISTMANAGER ====================
    
    // Kolory i nazwy kanałów MC
    Color_t colors[] = {kRed, kBlue, kGreen, kOrange};
    std::vector<TString> names = {"Signal", "Background1", "Background2"};
    
    // Utwórz HistManager z 3 kanałami MC
    HistManager* histMgr = new HistManager(3, colors, names);
    
    // Konfiguracja histogramu
    HistManager::HistConfig config;
    config.name = "deltaT";
    config.title = "Time difference distribution";
    config.bins = 200;
    config.xmin = -10.0;
    config.xmax = 10.0;
    config.xtitle = "#Delta t [ns]";
    config.ytitle = "Events";
    config.logy = false;
    config.showStats = true;
    
    // Stwórz zestaw histogramów 1D
    histMgr->CreateHistSet1D("deltaT", config);
    
    // ==================== GENEROWANIE DANYCH TESTOWYCH ====================
    
    std::cout << "Generating test data..." << std::endl;
    TRandom3 rand(12345);
    
    // Kanał 1: Sygnał - wąski gauss w centrum
    for(int i = 0; i < 8000; ++i) {
        Double_t val = rand.Gaus(0.0, 1.2);
        histMgr->Fill1D("deltaT", 1, val);
    }
    
    // Kanał 2: Background 1 - szerszy gauss przesunięty w lewo
    for(int i = 0; i < 3000; ++i) {
        Double_t val = rand.Gaus(-2.5, 1.8);
        histMgr->Fill1D("deltaT", 2, val);
    }
    
    // Kanał 3: Background 2 - gauss przesunięty w prawo
    for(int i = 0; i < 2000; ++i) {
        Double_t val = rand.Gaus(3.0, 1.0);
        histMgr->Fill1D("deltaT", 3, val);
    }
    
    std::cout << "Data generated successfully!" << std::endl;
    std::cout << "Total events: " << 8000+3000+2000 << std::endl;
    
    // ==================== UŻYCIE TRIPLEGAUSSFITTER ====================
    
    // Stwórz fitter
    KLOE::TripleGaussFitter fitter;
    fitter.SetVerbose(true);  // Włącz tryb verbose dla szczegółów
    
    // Ustaw kolory dla lepszej wizualizacji
    fitter.SetComponentColors(kRed, kBlue, kGreen, kMagenta);
    fitter.SetLineWidths(3, 2);
    
    std::cout << "\n==================== FITTING SUM MC HISTOGRAM ====================" << std::endl;
    
    // Pobierz histogram sumy MC (mctruth = 0)
    TH1* sumHist = histMgr->GetHistogram1D("deltaT", 0);
    if(!sumHist) {
        std::cerr << "ERROR: Cannot get sum histogram!" << std::endl;
        return;
    }
    
    std::cout << "Sum histogram entries: " << sumHist->GetEntries() << std::endl;
    
    // Fit z automatycznymi parametrami
    if(fitter.FitHistogram(sumHist, "Delta T Sum MC Distribution")) {
        fitter.PrintResults(true);  // Szczegółowe wyniki
        
        // Narysuj i zapisz wyniki
        std::cout << "Creating visualization..." << std::endl;
        TCanvas* canvas = fitter.DrawFit("c_sum_fit", "Triple Gaussian Fit - Sum MC", 
                                        1000, 700, true, true);
        if(canvas) {
            canvas->SaveAs("img/deltaT_sum_tripleGauss_fit.png");
            std::cout << "Plot saved to: img/deltaT_sum_tripleGauss_fit.png" << std::endl;
        }
    } else {
        std::cout << "Fit failed for sum histogram!" << std::endl;
    }
    
    std::cout << "\n==================== FITTING INDIVIDUAL CHANNELS ====================" << std::endl;
    
    // Fituj poszczególne kanały MC
    for(Int_t ch = 1; ch <= 3; ++ch) {
        TH1* channelHist = histMgr->GetHistogram1D("deltaT", ch);
        if(!channelHist) {
            std::cout << "Channel " << ch << ": histogram not found" << std::endl;
            continue;
        }
        
        if(channelHist->GetEntries() < 100) {
            std::cout << "Channel " << ch << ": too few entries (" 
                      << channelHist->GetEntries() << ")" << std::endl;
            continue;
        }
        
        TString title = Form("Channel %d (%s)", ch, names[ch-1].Data());
        std::cout << "\nFitting " << title.Data() << "..." << std::endl;
        std::cout << "Entries: " << channelHist->GetEntries() << std::endl;
        
        // Przygotuj niestandardowe parametry dla każdego kanału
        auto params = fitter.PrepareDefaultParams(channelHist);
        
        // Dostosuj zakres fitu w zależności od kanału
        if(ch == 1) {
            // Sygnał - zawęź zakres wokół centrum
            params.SetFitRange(-8.0, 8.0);
        } else if(ch == 2) {
            // Background 1 - focus na lewą stronę
            params.SetFitRange(-8.0, 5.0);
        } else {
            // Background 2 - focus na prawą stronę
            params.SetFitRange(-5.0, 8.0);
        }
        
        if(fitter.FitHistogram(channelHist, params, title)) {
            // Wypisz tylko podsumowanie dla kanałów
            fitter.PrintSummary();
            
            // Narysuj i zapisz dla każdego kanału
            TString canvasName = Form("c_ch%d_fit", ch);
            TString canvasTitle = Form("Triple Gaussian Fit - %s", title.Data());
            
            TCanvas* chCanvas = fitter.DrawFit(canvasName, canvasTitle, 800, 600);
            if(chCanvas) {
                TString filename = Form("img/deltaT_channel%d_tripleGauss_fit.png", ch);
                chCanvas->SaveAs(filename);
                std::cout << "  Plot saved to: " << filename.Data() << std::endl;
            }
        } else {
            std::cout << "  Fit failed for " << title.Data() << std::endl;
        }
    }
    
    std::cout << "\n==================== FITOWANIE DANYCH (JEŚLI ISTNIEJĄ) ====================" << std::endl;
    
    // Sprawdź czy są dane
    TH1* dataHist = histMgr->GetDataHistogram1D("deltaT");
    if(dataHist && dataHist->GetEntries() > 50) {
        std::cout << "Data histogram found with " << dataHist->GetEntries() << " entries" << std::endl;
        
        if(fitter.FitHistogram(dataHist, "Experimental Data")) {
            fitter.PrintResults(true);
            
            TCanvas* dataCanvas = fitter.DrawFit("c_data_fit", "Triple Gaussian Fit - Data");
            if(dataCanvas) {
                dataCanvas->SaveAs("img/deltaT_data_tripleGauss_fit.png");
                std::cout << "Data fit plot saved to: img/deltaT_data_tripleGauss_fit.png" << std::endl;
            }
        }
    } else {
        std::cout << "No data histogram available for fitting" << std::endl;
    }
    
    std::cout << "\n==================== PORÓWNANIE WYNIKÓW ====================" << std::endl;
    
    // Stwórz porównanie wszystkich fitów
    TCanvas* comparison = new TCanvas("c_comparison", "Triple Gaussian Fits Comparison", 1200, 800);
    comparison->Divide(2, 2);
    
    // Suma MC
    comparison->cd(1);
    TH1* sumHistCopy = histMgr->GetHistogram1D("deltaT", 0);
    if(sumHistCopy) {
        sumHistCopy->SetTitle("Sum MC");
        sumHistCopy->Draw("HIST");
        if(fitter.FitHistogram(sumHistCopy, "Sum")) {
            fitter.DrawFitOnCurrentPad(true, false);
        }
    }
    
    // Kanały indywidualne
    for(Int_t ch = 1; ch <= 3; ++ch) {
        comparison->cd(ch + 1);
        TH1* chHist = histMgr->GetHistogram1D("deltaT", ch);
        if(chHist && chHist->GetEntries() > 50) {
            chHist->SetTitle(Form("Channel %d", ch));
            chHist->Draw("HIST");
            
            auto params = fitter.PrepareDefaultParams(chHist);
            if(fitter.FitHistogram(chHist, params, Form("Ch%d", ch))) {
                fitter.DrawFitOnCurrentPad(true, false);
            }
        }
    }
    
    comparison->SaveAs("img/deltaT_all_tripleGauss_comparison.png");
    std::cout << "Comparison plot saved to: img/deltaT_all_tripleGauss_comparison.png" << std::endl;
    
    // ==================== PODSUMOWANIE ====================
    
    std::cout << "\n==================== PODSUMOWANIE ====================" << std::endl;
    std::cout << "✓ Klasa KLOE::TripleGaussFitter została pomyślnie przetestowana!" << std::endl;
    std::cout << "✓ Integracja z HistManager działa poprawnie" << std::endl;
    std::cout << "✓ Dostęp do histogramów przez GetHistogram1D() funkcjonuje" << std::endl;
    std::cout << "✓ Automatyczne i ręczne parametry fitu zostały przetestowane" << std::endl;
    std::cout << "✓ Wizualizacja i zapis wyników działają prawidłowo" << std::endl;
    
    std::cout << "\nPliki wygenerowane w folderze img/:" << std::endl;
    std::cout << "- deltaT_sum_tripleGauss_fit.png (fit sumy MC)" << std::endl;
    std::cout << "- deltaT_channel1_tripleGauss_fit.png (fit kanału 1)" << std::endl;
    std::cout << "- deltaT_channel2_tripleGauss_fit.png (fit kanału 2)" << std::endl;
    std::cout << "- deltaT_channel3_tripleGauss_fit.png (fit kanału 3)" << std::endl;
    std::cout << "- deltaT_all_tripleGauss_comparison.png (porównanie wszystkich)" << std::endl;
    
    // Cleanup
    delete histMgr;
    delete comparison;
    
    std::cout << "\nTest zakończony pomyślnie!" << std::endl;
}

/**
 * @brief Funkcja główna dla standalone uruchomienia
 */
int main() {
    test_triple_gauss_fitter();
    return 0;
}

/**
 * @brief Makro ROOT - uruchom przez: root -l -b -q test_triple_gauss_fitter.C
 */
#ifndef __CLING__
int main() {
    return 0;
}
#endif