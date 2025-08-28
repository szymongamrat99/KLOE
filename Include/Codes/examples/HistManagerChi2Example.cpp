#include "HistManager.h"
#include <iostream>
#include <TRandom.h>
#include <TApplication.h>

/**
 * @brief Example demonstrating the chi-squared calculation feature in HistManager
 * 
 * This example shows how to:
 * 1. Create histograms with MC and data
 * 2. Enable automatic chi-squared calculation
 * 3. Display results with chi-squared information
 */

void HistManagerChi2Example() {
    // Initialize with 3 MC channels
    Color_t colors[] = {kRed, kBlue, kGreen};
    std::vector<TString> names = {"Signal", "Background1", "Background2"};
    HistManager* histMgr = new HistManager(3, colors, names);
    
    // Configure histogram
    HistManager::HistConfig config;
    config.name = "mass";
    config.title = "Invariant Mass Distribution";
    config.bins = 100;
    config.xmin = 0.4;
    config.xmax = 0.6;
    config.xtitle = "Mass [GeV/c^{2}]";
    config.ytitle = "Events / 2 MeV";
    config.logy = false;
    config.showStats = false; // Wyłącz statBox
    
    histMgr->CreateHistSet1D("mass", config);
    
    // Enable automatic chi-squared calculation
    histMgr->SetCalculateDataMCChi2(true);
    
    // Generate some sample data
    TRandom* rand = new TRandom(12345);
    
    // Fill MC histograms with different distributions
    for(int i = 0; i < 10000; ++i) {
        // Signal: Gaussian around 0.498 GeV (K0 mass)
        if(rand->Uniform() < 0.4) {
            Double_t mass = rand->Gaus(0.498, 0.005);
            if(mass > 0.4 && mass < 0.6) {
                histMgr->Fill1D("mass", 1, mass, 1.0);
            }
        }
        
        // Background 1: Exponential
        if(rand->Uniform() < 0.3) {
            Double_t mass = 0.4 + rand->Exp(0.05);
            if(mass > 0.4 && mass < 0.6) {
                histMgr->Fill1D("mass", 2, mass, 1.0);
            }
        }
        
        // Background 2: Flat
        if(rand->Uniform() < 0.2) {
            Double_t mass = rand->Uniform(0.4, 0.6);
            histMgr->Fill1D("mass", 3, mass, 1.0);
        }
    }
    
    // Generate data that's similar but not identical to MC sum
    for(int i = 0; i < 8500; ++i) {
        Double_t mass;
        Double_t type = rand->Uniform();
        
        if(type < 0.45) {
            // Signal component (slightly different from MC)
            mass = rand->Gaus(0.4985, 0.0048);
        } else if(type < 0.75) {
            // Background 1 component
            mass = 0.4 + rand->Exp(0.048);
        } else {
            // Background 2 component
            mass = rand->Uniform(0.4, 0.6);
        }
        
        if(mass > 0.4 && mass < 0.6) {
            histMgr->FillData1D("mass", mass);
        }
    }
    
    std::cout << "\n=== HistManager Chi2 Example ===" << std::endl;
    std::cout << "Generated MC and data histograms" << std::endl;
    
    // Example 1: Draw without normalization (chi2 calculation only)
    std::cout << "\n1. Drawing without normalization..." << std::endl;
    histMgr->SetNormalizationType(HistManager::NormalizationType::NONE);
    histMgr->DrawSet1D("mass", "HIST", true);
    
    // Example 2: Draw with simple scaling and chi2
    std::cout << "\n2. Drawing with simple scaling..." << std::endl;
    histMgr->SetNormalizationType(HistManager::NormalizationType::SIMPLE_SCALE);
    histMgr->DrawSet1D("mass", "HIST", true);
    
    // Example 3: Draw with fraction fit and chi2
    std::cout << "\n3. Drawing with fraction fit..." << std::endl;
    histMgr->SetNormalizationType(HistManager::NormalizationType::FRACTION_FIT);
    
    // Set some fit constraints
    HistManager::FitConstraints constraints(3);
    constraints.lowerBounds = {0.0, 0.0, 0.0};
    constraints.upperBounds = {2.0, 2.0, 2.0};
    constraints.fitRangeMin = 10;  // Use bins 10-90 for fit
    constraints.fitRangeMax = 90;
    histMgr->SetFitConstraints("mass", constraints);
    
    histMgr->DrawSet1D("mass", "HIST", true);
    
    // Example 4: Manual chi2 calculation
    std::cout << "\n4. Manual chi2 calculation..." << std::endl;
    auto chi2Result = histMgr->CalculateDataMCChi2("mass", 5.0); // min 5 events per bin
    
    // Example 5: Get stored results
    std::cout << "\n5. Accessing stored results..." << std::endl;
    if(histMgr->HasFitResult("mass")) {
        auto fitResult = histMgr->GetLastFitResult();
        std::cout << "Last fit result:" << std::endl;
        std::cout << "  Converged: " << (fitResult.converged ? "Yes" : "No") << std::endl;
        std::cout << "  Chi2/NDF: " << fitResult.chi2 << "/" << fitResult.ndf << std::endl;
        for(size_t i = 0; i < fitResult.fractions.size(); ++i) {
            std::cout << "  " << names[i] << ": " << fitResult.fractions[i] 
                      << " ± " << fitResult.errors[i] << std::endl;
        }
    }
    
    auto lastChi2 = histMgr->GetLastChi2Result();
    std::cout << "Last chi2 result:" << std::endl;
    std::cout << "  Chi2/NDF: " << lastChi2.chi2 << "/" << lastChi2.ndf 
              << " = " << lastChi2.chi2_ndf << std::endl;
    std::cout << "  p-value: " << lastChi2.pValue << std::endl;
    std::cout << "  Quality: " << lastChi2.comparisonInfo << std::endl;
    
    // Save results
    histMgr->SaveSet("mass", "chi2_example_{channel}.png");
    histMgr->SaveToRoot("chi2_example_results.root");
    
    std::cout << "\nResults saved to chi2_example_results.root and PNG files" << std::endl;
    
    delete rand;
    delete histMgr;
}

/**
 * @brief Advanced example showing different chi2 scenarios
 */
void AdvancedChi2Example() {
    std::cout << "\n=== Advanced Chi2 Scenarios ===" << std::endl;
    
    Color_t colors[] = {kRed, kBlue};
    std::vector<TString> names = {"Signal", "Background"};
    HistManager* histMgr = new HistManager(2, colors, names);
    
    HistManager::HistConfig config;
    config.name = "scenario";
    config.title = "Chi2 Test Scenarios";
    config.bins = 50;
    config.xmin = 0.0;
    config.xmax = 10.0;
    config.xtitle = "Variable";
    config.ytitle = "Events";
    config.showStats = false; // Wyłącz statBox
    
    TRandom* rand = new TRandom(54321);
    
    // Scenario 1: Perfect agreement
    std::cout << "\nScenario 1: Perfect agreement between MC and data" << std::endl;
    histMgr->CreateHistSet1D("perfect", config);
    
    for(int i = 0; i < 5000; ++i) {
        Double_t x = rand->Gaus(5.0, 1.5);
        if(x > 0 && x < 10) {
            histMgr->Fill1D("perfect", 1, x);
            histMgr->FillData1D("perfect", x); // Identical to MC
        }
    }
    
    histMgr->SetNormalizationType(HistManager::NormalizationType::SIMPLE_SCALE);
    histMgr->DrawSet1D("perfect", "HIST", true);
    
    // Scenario 2: Poor agreement
    std::cout << "\nScenario 2: Poor agreement between MC and data" << std::endl;
    histMgr->CreateHistSet1D("poor", config);
    
    // MC: Gaussian at 5.0
    for(int i = 0; i < 3000; ++i) {
        Double_t x = rand->Gaus(5.0, 1.0);
        if(x > 0 && x < 10) {
            histMgr->Fill1D("poor", 1, x);
        }
    }
    
    // Data: Gaussian at 6.0 (shifted)
    for(int i = 0; i < 3000; ++i) {
        Double_t x = rand->Gaus(6.0, 1.0);
        if(x > 0 && x < 10) {
            histMgr->FillData1D("poor", x);
        }
    }
    
    histMgr->DrawSet1D("poor", "HIST", true);
    
    // Scenario 3: Low statistics
    std::cout << "\nScenario 3: Low statistics" << std::endl;
    histMgr->CreateHistSet1D("lowstat", config);
    
    for(int i = 0; i < 50; ++i) {
        Double_t x = rand->Uniform(0, 10);
        histMgr->Fill1D("lowstat", 1, x);
        histMgr->FillData1D("lowstat", x + rand->Gaus(0, 0.1));
    }
    
    histMgr->DrawSet1D("lowstat", "HIST", true);
    
    delete rand;
    delete histMgr;
}

int main(int argc, char** argv) {
    TApplication app("HistManagerChi2Example", &argc, argv);
    
    HistManagerChi2Example();
    AdvancedChi2Example();
    
    std::cout << "\n=== Chi2 Calculation Features Summary ===" << std::endl;
    std::cout << "1. Automatic chi2 calculation when drawing data" << std::endl;
    std::cout << "2. Manual chi2 calculation with CalculateDataMCChi2()" << std::endl;
    std::cout << "3. Customizable minimum bin content threshold" << std::endl;
    std::cout << "4. Chi2 text display on canvas" << std::endl;
    std::cout << "5. Statistical quality assessment" << std::endl;
    std::cout << "6. p-value calculation" << std::endl;
    std::cout << "7. Integration with normalization methods" << std::endl;
    
    return 0;
}
