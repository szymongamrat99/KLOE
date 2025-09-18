#include "Include/Codes/inc/HistManager.h"
#include <TRandom3.h>
#include <TApplication.h>
#include <iostream>

void test_triple_gauss() {
    // Inicjalizuj HistManager z prostymi kolorami
    Color_t colors[] = {kRed, kBlue, kGreen};
    std::vector<TString> names = {"Signal", "Background1", "Background2"};
    HistManager* histMgr = new HistManager(3, colors, names);
    
    // Konfiguracja histogramu do testowania
    HistManager::HistConfig config;
    config.name = "test_gauss";
    config.title = "Test Triple Gaussian";
    config.bins = 100;
    config.xmin = -10.0;
    config.xmax = 10.0;
    config.xtitle = "x [units]";
    config.ytitle = "Events";
    
    histMgr->CreateHistSet1D("test", config);
    
    // Generuj dane symulujące triple gaussian
    TRandom3 rand(12345);
    
    // Wypełnij histogram danymi simulującymi triple gaussian
    for(int i = 0; i < 10000; ++i) {
        Double_t val = 0;
        Double_t r = rand.Rndm();
        
        if(r < 0.5) {
            // Pierwszy gauss: μ=-2, σ=0.8
            val = rand.Gaus(-2.0, 0.8);
        } else if(r < 0.8) {
            // Drugi gauss: μ=0, σ=1.2  
            val = rand.Gaus(0.0, 1.2);
        } else {
            // Trzeci gauss: μ=3, σ=0.6
            val = rand.Gaus(3.0, 0.6);
        }
        
        // Wypełnij jako dane (mctruth=-1)
        histMgr->FillData1D("test", val);
    }
    
    // Przygotuj parametry Triple Gaussian
    HistManager::FitParams3Gauss fitParams;
    
    // Ustaw początkowe parametry (A1, μ1, σ1, A2, μ2, σ2, A3, μ3, σ3)
    std::vector<Double_t> initial = {
        2000.0, -2.0, 1.0,   // Pierwszy gauss
        3000.0,  0.0, 1.5,   // Drugi gauss  
        1500.0,  3.0, 0.8    // Trzeci gauss
    };
    fitParams.SetInitialParams(initial);
    
    // Ustaw granice parametrów
    std::vector<Double_t> lower = {
        100.0, -5.0, 0.1,    // Minimum dla pierwszego gaussa
        100.0, -3.0, 0.1,    // Minimum dla drugiego gaussa
        100.0,  1.0, 0.1     // Minimum dla trzeciego gaussa
    };
    std::vector<Double_t> upper = {
        10000.0, 0.0, 3.0,   // Maximum dla pierwszego gaussa
        10000.0, 3.0, 3.0,   // Maximum dla drugiego gaussa
        10000.0, 6.0, 3.0    // Maximum dla trzeciego gaussa
    };
    fitParams.SetParamBounds(lower, upper);
    
    // Ustaw zakres fitu
    fitParams.SetFitRange(-8.0, 8.0);
    
    std::cout << "=== Testing Triple Gaussian Fitting ===" << std::endl;
    
    // Wykonaj fit i narysuj histogram z fitem
    TCanvas* canvas = histMgr->DrawHistogram1DWithFit("test", -1, "HIST", true, &fitParams, false);
    
    if(canvas) {
        canvas->SaveAs("/data/ssd/gamrat/KLOE/build/test_triple_gauss.png");
        std::cout << "Canvas saved to build/test_triple_gauss.png" << std::endl;
    }
    
    // Sprawdź wyniki fitu
    const HistManager::FitParams3Gauss* results = histMgr->GetTripleGaussResults1D("test", -1);
    
    if(results && results->converged) {
        std::cout << "\n=== Fit Results ===" << std::endl;
        std::cout << "Fit converged: " << (results->converged ? "YES" : "NO") << std::endl;
        std::cout << "Chi2/NDF: " << results->chi2 << "/" << results->ndf << " = " 
                  << (results->ndf > 0 ? results->chi2/results->ndf : 0) << std::endl;
        std::cout << "Combined mean: " << results->combinedMean << " ± " << results->combinedMeanErr << std::endl;
        std::cout << "Combined sigma: " << results->combinedStdDev << " ± " << results->combinedStdDevErr << std::endl;
    } else {
        std::cout << "Fit failed or results not available!" << std::endl;
    }
    
    delete histMgr;
}

int main(int argc, char** argv) {
    TApplication app("TripleGaussTest", &argc, argv);
    
    test_triple_gauss();
    
    std::cout << "Test completed. Press Enter to exit..." << std::endl;
    getchar();
    
    return 0;
}
