#include "HistManager.h"
#include <TMath.h>
#include <TStyle.h>
#include <Math/MinimizerOptions.h>
#include <algorithm>
#include <stdexcept>
#include <iostream>

HistManager::HistManager(Int_t channNum, const Color_t* channColors, 
                                 const std::vector<TString>& channelNames,
                                 Int_t dataStyle, Int_t dataColor, Int_t sumColor) 
    : fChannNum(channNum), fDataStyle(dataStyle), fDataColor(dataColor), fSumColor(sumColor) {
    
    for(Int_t i = 0; i < channNum; ++i) {
        fChannColors.push_back(channColors[i]);
    }
    
    // Initialize channel names
    if(channelNames.empty()) {
        // Use default names if none provided
        for(Int_t i = 0; i < channNum; ++i) {
            fChannelNames.push_back(Form("Channel %d", i+1));
        }
    } else {
        // Use provided names
        if(channelNames.size() != static_cast<size_t>(channNum)) {
            throw std::runtime_error("Number of channel names doesn't match number of channels");
        }
        fChannelNames = channelNames;
    }
}

HistManager::~HistManager() {
    // Cleanup all sets
    std::vector<TString> setsToClean;
    for(const auto& pair : fHists1D) {
        setsToClean.push_back(pair.first);
    }
    // Cleanup data histograms
    for(const auto& pair : fData1D) {
        delete pair.second;
    }
    for(const auto& pair : fData2D) {
        delete pair.second;
    }
    for(const auto& setName : setsToClean) {
        CleanupSet(setName);
    }
    
    // Cleanup array histograms
    std::vector<TString> arrayToClean;
    for(const auto& pair : fArrayHists1D) {
        arrayToClean.push_back(pair.first);
    }
    for(const auto& baseName : arrayToClean) {
        CleanupArraySet(baseName);
    }
}

void HistManager::CreateHistSet1D(const TString& setName, const HistConfig& config) {
    if(fHists1D.find(setName) != fHists1D.end()) {
        CleanupSet(setName);
    }

    std::vector<TH1*> hists;
    // Tworzenie histogramu sumy jako pierwszy element wektora
    TString sumHistName = Form("%s_sum", config.name.Data());
    TH1D* sumHist = new TH1D(sumHistName, config.title,
                            config.bins, config.xmin, config.xmax);
    sumHist->GetXaxis()->SetTitle(config.xtitle);
    sumHist->GetYaxis()->SetTitle(config.ytitle);
    ConfigureHistogram(sumHist, fSumColor, config.showStats);
    hists.push_back(sumHist);

    // Tworzenie pozostałych histogramów
    for(Int_t i = 0; i < fChannNum; ++i) {
        TString histName = Form("%s_%d", config.name.Data(), i+1);
        TH1D* hist = new TH1D(histName, config.title, 
                             config.bins, config.xmin, config.xmax);
        hist->GetXaxis()->SetTitle(config.xtitle);
        hist->GetYaxis()->SetTitle(config.ytitle);
        
        ConfigureHistogram(hist, fChannColors[i], config.showStats);
        hists.push_back(hist);
    }

    fHists1D[setName] = hists;
    fConfigs1D[setName] = config;
}

void HistManager::CreateHistSet2D(const TString& setName, const Hist2DConfig& config) {
    if(fHists2D.find(setName) != fHists2D.end()) {
        CleanupSet(setName);
    }

    std::vector<TH2*> hists;
    for(Int_t i = 0; i < fChannNum; ++i) {
        TString histName = Form("%s_%d", config.name.Data(), i+1);
        TH2D* hist = new TH2D(histName, config.title,
                             config.bins, config.xmin, config.xmax,
                             config.binsy, config.ymin, config.ymax);
        hist->GetXaxis()->SetTitle(config.xtitle);
        hist->GetYaxis()->SetTitle(config.ytitle);
        
        ConfigureHistogram(hist, fChannColors[i], config.showStats);
        hists.push_back(hist);
    }

    fHists2D[setName] = hists;
    fConfigs2D[setName] = config;
}

void HistManager::Fill1D(const TString& setName, Int_t mctruth, Double_t value, Double_t weight) {
    if(mctruth < 1 || mctruth > fChannNum) {
        return; // Invalid mctruth, silently ignore
    }
    
    auto it = fHists1D.find(setName);
    if(it == fHists1D.end()) {
        throw std::runtime_error("Histogram set not found: " + std::string(setName.Data()));
    }

    // Fill suma MC (indeks 0)
    it->second[0]->Fill(value, weight);
    // Fill pojedyncza składowa (przesunięte o 1 przez sumę)
    it->second[mctruth]->Fill(value, weight);
}

void HistManager::Fill2D(const TString& setName, Int_t mctruth, Double_t x, Double_t y, Double_t weight) {
    if(mctruth < 1 || mctruth > fChannNum) {
        return; // Invalid mctruth, silently ignore
    }
    
    auto it = fHists2D.find(setName);
    if(it == fHists2D.end()) {
        throw std::runtime_error("Histogram set not found: " + std::string(setName.Data()));
    }

    it->second[mctruth-1]->Fill(x, y, weight);
}

void HistManager::FillData1D(const TString& setName, Double_t value, Double_t weight) {
    auto it = fHists1D.find(setName);
    if(it == fHists1D.end()) {
        throw std::runtime_error("Histogram set not found: " + std::string(setName.Data()));
    }

    // Create data histogram if it doesn't exist
    if(fData1D.find(setName) == fData1D.end()) {
        auto configIt = fConfigs1D.find(setName);
        if(configIt == fConfigs1D.end()) {
            throw std::runtime_error("Configuration not found for set: " + std::string(setName.Data()));
        }

        TH1D* hist = new TH1D(Form("%s_data", configIt->second.name.Data()),
                             configIt->second.title,
                             configIt->second.bins,
                             configIt->second.xmin,
                             configIt->second.xmax);
        hist->GetXaxis()->SetTitle(configIt->second.xtitle);
        hist->GetYaxis()->SetTitle(configIt->second.ytitle);
        fData1D[setName] = hist;
    }

    fData1D[setName]->Fill(value, weight);
}

void HistManager::FillData2D(const TString& setName, Double_t x, Double_t y, Double_t weight) {
    auto it = fHists2D.find(setName);
    if(it == fHists2D.end()) {
        throw std::runtime_error("Histogram set not found: " + std::string(setName.Data()));
    }

    // Create data histogram if it doesn't exist
    if(fData2D.find(setName) == fData2D.end()) {
        auto configIt = fConfigs2D.find(setName);
        if(configIt == fConfigs2D.end()) {
            throw std::runtime_error("Configuration not found for set: " + std::string(setName.Data()));
        }

        TH2D* hist = new TH2D(Form("%s_data", configIt->second.name.Data()),
                             configIt->second.title,
                             configIt->second.bins,
                             configIt->second.xmin,
                             configIt->second.xmax,
                             configIt->second.binsy,
                             configIt->second.ymin,
                             configIt->second.ymax);
        hist->GetXaxis()->SetTitle(configIt->second.xtitle);
        hist->GetYaxis()->SetTitle(configIt->second.ytitle);
        fData2D[setName] = hist;
    }

    fData2D[setName]->Fill(x, y, weight);
}

void HistManager::DrawSet1D(const TString& setName, const TString& drawOpt, Bool_t drawData) {
    auto histIt = fHists1D.find(setName);
    if(histIt == fHists1D.end()) {
        throw std::runtime_error("Histogram set not found: " + std::string(setName.Data()));
    }

    auto configIt = fConfigs1D.find(setName);
    if(configIt == fConfigs1D.end()) {
        throw std::runtime_error("Configuration not found for set: " + std::string(setName.Data()));
    }

    TCanvas* canvas = new TCanvas(Form("c_%s", setName.Data()), 
                                 configIt->second.title, 790, 790);
    canvas->SetLeftMargin(0.15);
    canvas->SetBottomMargin(0.15);
    
    if(configIt->second.logy) {
        canvas->SetLogy();
    }

    // Find maximum Y value across all histograms
    Double_t maxY = CalculateYRange(histIt->second, setName, configIt->second.logy);

    // Draw histograms
    bool first = true;
    for(size_t i = 1; i < histIt->second.size(); ++i) {
        histIt->second[i]->GetYaxis()->SetRangeUser(configIt->second.logy ? 0.1 : 0, maxY);
        histIt->second[i]->Draw(first ? drawOpt : drawOpt + "SAME");
        first = false;
    }
    
    // Draw sum histogram last, always as HIST
    histIt->second[0]->GetYaxis()->SetRangeUser(configIt->second.logy ? 0.1 : 0, maxY);
    histIt->second[0]->Draw("HIST SAME");

    // Draw data if requested and available
    if(drawData) {
        auto dataIt = fData1D.find(setName);
        if(dataIt != fData1D.end()) {
            // Wykonaj normalizację w zależności od ustawionego typu
            if(fNormalizationType == NormalizationType::FRACTION_FIT) {
                try {
                    FitResult fitResult = PerformFractionFit(setName);
                    
                    if(fitResult.converged && fitResult.status == 0) {
                        // Po ficie, ponownie oblicz zakres Y uwzględniając przeskalowane histogramy
                        maxY = CalculateYRange(histIt->second, setName, configIt->second.logy);
                        
                        // Zaktualizuj zakres Y dla wszystkich histogramów
                        for(size_t i = 0; i < histIt->second.size(); ++i) {
                            histIt->second[i]->GetYaxis()->SetRangeUser(configIt->second.logy ? 0.1 : 0, maxY);
                        }
                        
                        // Wyświetl informacje o ficie
                        TPaveText* fitInfo = new TPaveText(0.15, 0.7, 0.5, 0.9, "NDC");
                        fitInfo->SetFillColor(kWhite);
                        fitInfo->SetBorderSize(1);
                        fitInfo->AddText("Fraction Fit Results:");
                        fitInfo->AddText(Form("#chi^{2}/NDF = %.2f/%d", fitResult.chi2, fitResult.ndf));
                        
                        for(Int_t i = 0; i < fChannNum && i < static_cast<Int_t>(fitResult.fractions.size()); ++i) {
                            if(fitResult.fractions[i] > 0) {
                                fitInfo->AddText(Form("%s: %.3f #pm %.3f", 
                                                    fChannelNames[i].Data(), 
                                                    fitResult.fractions[i], 
                                                    fitResult.errors[i]));
                            } else {
                                fitInfo->AddText(Form("%s: not used (empty)", 
                                                    fChannelNames[i].Data()));
                            }
                        }
                        
                        if(!fitResult.fitInfo.IsNull()) {
                            fitInfo->AddText(fitResult.fitInfo);
                        }
                        
                        fitInfo->Draw();
                    } else {
                        // Wyświetl informację o nieudanym ficie
                        TPaveText* fitInfo = new TPaveText(0.15, 0.8, 0.5, 0.9, "NDC");
                        fitInfo->SetFillColor(kYellow);
                        fitInfo->SetBorderSize(1);
                        fitInfo->AddText("Fraction Fit Failed");
                        if(!fitResult.fitInfo.IsNull()) {
                            fitInfo->AddText(fitResult.fitInfo);
                        }
                        fitInfo->Draw();
                    }
                } catch(const std::exception& e) {
                    std::cerr << "Error during FractionFit: " << e.what() << std::endl;
                }
            } else if(fNormalizationType == NormalizationType::SIMPLE_SCALE) {
                try {
                    Double_t scaleFactor = PerformSimpleScaling(setName);
                    
                    if(scaleFactor > 0) {
                        // Po skalowaniu, ponownie oblicz zakres Y
                        maxY = CalculateYRange(histIt->second, setName, configIt->second.logy);
                        
                        // Zaktualizuj zakres Y dla wszystkich histogramów
                        for(size_t i = 0; i < histIt->second.size(); ++i) {
                            histIt->second[i]->GetYaxis()->SetRangeUser(configIt->second.logy ? 0.1 : 0, maxY);
                        }
                        
                        // Wyświetl informacje o skalowaniu
                        TPaveText* scaleInfo = new TPaveText(0.15, 0.8, 0.5, 0.9, "NDC");
                        scaleInfo->SetFillColor(kCyan-10);
                        scaleInfo->SetBorderSize(1);
                        scaleInfo->AddText("Simple Scaling Applied:");
                        scaleInfo->AddText(Form("Scale Factor = %.3f", scaleFactor));
                        scaleInfo->AddText(Form("Data Events: %.0f", dataIt->second->Integral()));
                        scaleInfo->AddText(Form("MC Events: %.0f", histIt->second[0]->Integral()));
                        scaleInfo->Draw();
                    }
                } catch(const std::exception& e) {
                    std::cerr << "Error during Simple Scaling: " << e.what() << std::endl;
                }
            }
            
            dataIt->second->SetMarkerStyle(fDataStyle);
            dataIt->second->SetMarkerColor(fDataColor);
            dataIt->second->SetLineColor(fDataColor);
            // Upewnij się, że dane mają poprawny zakres Y (może być zaktualizowany po normalizacji)
            if(fNormalizationType != NormalizationType::NONE) {
                dataIt->second->GetYaxis()->SetRangeUser(configIt->second.logy ? 0.1 : 0, maxY);
            }
            dataIt->second->Draw("SAME PE1");
        }
    }

    // Add legend with data if present
    TH1* dataHist = nullptr;
    if(drawData && fData1D.find(setName) != fData1D.end()) {
        dataHist = fData1D[setName];
    }
    TLegend* legend = CreateLegend(histIt->second, dataHist);
    legend->Draw();

    canvas->Update();
    
    if(fCanvases.find(setName) == fCanvases.end()) {
        fCanvases[setName] = std::vector<TCanvas*>();
    }
    fCanvases[setName].push_back(canvas);
}

void HistManager::DrawSet2D(const TString& setName, const TString& drawOpt, Bool_t drawData) {
    auto histIt = fHists2D.find(setName);
    if(histIt == fHists2D.end()) {
        throw std::runtime_error("Histogram set not found: " + std::string(setName.Data()));
    }

    std::vector<TCanvas*> canvases;
    for(Int_t i = 0; i < fChannNum; ++i) {
        TCanvas* canvas = new TCanvas(Form("c_%s_%d", setName.Data(), i+1),
                                    Form("%s (Channel %d)", setName.Data(), i+1),
                                    790, 790);
        canvas->SetLeftMargin(0.15);
        canvas->SetBottomMargin(0.15);
        canvas->SetRightMargin(0.15); // For z-axis color scale

        histIt->second[i]->Draw(drawOpt);
        canvas->Update();
        canvases.push_back(canvas);
    }

    fCanvases[setName] = canvases;
}

void HistManager::SaveSet(const TString& setName, const TString& filePattern) {
    ExportSet(setName, filePattern, ImageFormat::PNG);
}

void HistManager::SaveToRoot(const TString& filename) {
    TFile* outFile = TFile::Open(filename, "RECREATE");
    if(!outFile || outFile->IsZombie()) {
        throw std::runtime_error("Cannot create ROOT file: " + std::string(filename.Data()));
    }

    // Save all 1D histograms
    for(const auto& pair : fHists1D) {
        TDirectory* dir = outFile->mkdir(pair.first);
        dir->cd();
        for(size_t i = 0; i < pair.second.size(); ++i) {
            pair.second[i]->Write(Form("channel_%lu", i+1));
        }
        // Save data histogram if exists
        auto dataIt = fData1D.find(pair.first);
        if(dataIt != fData1D.end()) {
            dataIt->second->Write("data");
        }
    }

    // Save all 2D histograms
    for(const auto& pair : fHists2D) {
        TDirectory* dir = outFile->mkdir(pair.first + "_2D");
        dir->cd();
        for(size_t i = 0; i < pair.second.size(); ++i) {
            pair.second[i]->Write(Form("channel_%lu", i+1));
        }
        // Save data histogram if exists
        auto dataIt = fData2D.find(pair.first);
        if(dataIt != fData2D.end()) {
            dataIt->second->Write("data");
        }
    }

    outFile->Close();
    delete outFile;
}

void HistManager::SaveSetToRoot(const TString& setName, const TString& filename) {
    TFile* outFile = TFile::Open(filename, "RECREATE");
    if(!outFile || outFile->IsZombie()) {
        throw std::runtime_error("Cannot create ROOT file: " + std::string(filename.Data()));
    }

    // Save 1D histograms if they exist
    auto hist1DIt = fHists1D.find(setName);
    if(hist1DIt != fHists1D.end()) {
        for(size_t i = 0; i < hist1DIt->second.size(); ++i) {
            hist1DIt->second[i]->Write(Form("channel_%lu", i+1));
        }
        // Save data histogram if exists
        auto dataIt = fData1D.find(setName);
        if(dataIt != fData1D.end()) {
            dataIt->second->Write("data");
        }
    }

    // Save 2D histograms if they exist
    auto hist2DIt = fHists2D.find(setName);
    if(hist2DIt != fHists2D.end()) {
        for(size_t i = 0; i < hist2DIt->second.size(); ++i) {
            hist2DIt->second[i]->Write(Form("channel_%lu", i+1));
        }
        // Save data histogram if exists
        auto dataIt = fData2D.find(setName);
        if(dataIt != fData2D.end()) {
            dataIt->second->Write("data");
        }
    }

    outFile->Close();
    delete outFile;
}

void HistManager::ExportSet(const TString& setName, const TString& filePattern, ImageFormat format) {
    auto canvasIt = fCanvases.find(setName);
    if(canvasIt == fCanvases.end()) {
        throw std::runtime_error("No canvases found for set: " + std::string(setName.Data()));
    }

    TString extension;
    switch(format) {
        case ImageFormat::PNG:
            extension = ".png";
            break;
        case ImageFormat::SVG:
            extension = ".svg";
            break;
        case ImageFormat::PDF:
            extension = ".pdf";
            break;
    }

    for(size_t i = 0; i < canvasIt->second.size(); ++i) {
        TString filename = filePattern;
        filename.ReplaceAll("{channel}", Form("%lu", i+1));
        // Add extension if not present
        if(!filename.EndsWith(extension)) {
            filename += extension;
        }
        canvasIt->second[i]->SaveAs(filename);
    }
}

void HistManager::ConfigureHistogram(TH1* hist, Int_t color, Bool_t showStats) {
    hist->SetLineColor(color);
    hist->SetMarkerColor(color);
    hist->GetYaxis()->CenterTitle(true);
    hist->GetXaxis()->CenterTitle(true);
    hist->GetXaxis()->SetMaxDigits(3);
    hist->GetXaxis()->SetTitleOffset(1.2);  // Adjust title position
    hist->GetYaxis()->SetTitleOffset(1.4);  // Adjust title position
    hist->SetStats(showStats);
    gStyle->SetTitleFont(62, "XYZ");       // Use Times New Roman for axis titles
    hist->GetXaxis()->SetTitleFont(62);    // Set font for X axis title
    hist->GetYaxis()->SetTitleFont(62);    // Set font for Y axis title
    hist->SetTitleFont(62);                // Set font for main title
    gStyle->SetLabelFont(62, "XYZ");       // Use Times New Roman for axis labels
}

Double_t HistManager::CalculateYRange(const std::vector<TH1*>& hists, const TString& setName, Bool_t logy) {
    Double_t maxY = 0;
    for(auto hist : hists) {
        maxY = TMath::Max(maxY, hist->GetMaximum());
    }

    // Check data histogram if it exists
    auto dataIt = fData1D.find(setName);
    if(dataIt != fData1D.end()) {
        maxY = TMath::Max(maxY, dataIt->second->GetMaximum());
    }

    return logy ? maxY * 5 : maxY * 1.1;
}

void HistManager::CleanupSet(const TString& setName) {
    // Clean up 1D histograms
    auto it1D = fHists1D.find(setName);
    if(it1D != fHists1D.end()) {
        for(auto hist : it1D->second) {
            delete hist;
        }
        fHists1D.erase(it1D);
    }

    // Clean up 2D histograms
    auto it2D = fHists2D.find(setName);
    if(it2D != fHists2D.end()) {
        for(auto hist : it2D->second) {
            delete hist;
        }
        fHists2D.erase(it2D);
    }

    // Clean up canvases
    auto itCanvas = fCanvases.find(setName);
    if(itCanvas != fCanvases.end()) {
        for(auto canvas : itCanvas->second) {
            delete canvas;
        }
        fCanvases.erase(itCanvas);
    }
}

TLegend* HistManager::CreateLegend(const std::vector<TH1*>& hists, const TH1* dataHist) {
    TLegend* legend = new TLegend(0.7, 0.7, 0.9, 0.9);
    // legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    
    // Dodaj najpierw sumę MC
    legend->AddEntry(hists[0], "MC sum", "l");
    
    // Dodaj pozostałe składowe
    for(size_t i = 1; i < hists.size(); ++i) {
        legend->AddEntry(hists[i], fChannelNames[i-1], "l");
    }
    
    if(dataHist) {
        legend->AddEntry(dataHist, "Data", "pe");
    }
    
    return legend;
}

// ==================== FRACTIONFIT IMPLEMENTATION ====================

void HistManager::SetFitConstraints(const TString& setName, const FitConstraints& constraints) {
    // Sprawdź czy zestaw istnieje
    if(fHists1D.find(setName) == fHists1D.end()) {
        throw std::runtime_error("Histogram set not found: " + std::string(setName.Data()));
    }
    
    // Sprawdź poprawność ograniczeń
    if(constraints.lowerBounds.size() != static_cast<size_t>(fChannNum) || 
       constraints.upperBounds.size() != static_cast<size_t>(fChannNum) ||
       constraints.initialValues.size() != static_cast<size_t>(fChannNum)) {
        throw std::runtime_error("Constraint vectors must have size equal to number of channels");
    }
    
    // Sprawdź zakresy
    for(Int_t i = 0; i < fChannNum; ++i) {
        if(constraints.lowerBounds[i] < 0.0 || 
           constraints.upperBounds[i] < 0.0 ||
           constraints.lowerBounds[i] >= constraints.upperBounds[i]) {
            throw std::runtime_error(Form("Invalid constraints for channel %d", i));
        }
    }
    
    fFitConstraints[setName] = constraints;
}

HistManager::FitResult HistManager::PerformFractionFit(const TString& setName, Bool_t useStoredConstraints) {
    auto histIt = fHists1D.find(setName);
    if(histIt == fHists1D.end()) {
        throw std::runtime_error("Histogram set not found: " + std::string(setName.Data()));
    }
    
    auto dataIt = fData1D.find(setName);
    if(dataIt == fData1D.end()) {
        throw std::runtime_error("No data histogram found for set: " + std::string(setName.Data()));
    }
    
    // Przygotuj wektor histogramów MC (bez sumy - indeks 0)
    std::vector<TH1*> mcHists;
    for(size_t i = 1; i < histIt->second.size(); ++i) {
        mcHists.push_back(histIt->second[i]);
    }
    
    // Użyj zapisanych ograniczeń lub domyślnych
    FitConstraints constraints(fChannNum);
    if(useStoredConstraints) {
        auto constraintIt = fFitConstraints.find(setName);
        if(constraintIt != fFitConstraints.end()) {
            constraints = constraintIt->second;
        }
    }
    
    // Wykonaj fit
    FitResult result = DoFractionFit(mcHists, dataIt->second, constraints);
    
    // Zapisz wyniki
    fFitResults[setName] = result;
    fLastFitResult = result;
    
    // Jeśli fit się udał, zaktualizuj histogramy
    if(result.converged && result.status == 0) {
        UpdateMCHistograms(mcHists, histIt->second[0], result);
    }
    
    return result;
}

HistManager::FitResult HistManager::DoFractionFit(const std::vector<TH1*>& mcHists, TH1* dataHist, 
                                                 const FitConstraints& constraints) {
    FitResult result;
    result.fractions.resize(fChannNum);
    result.errors.resize(fChannNum);
    
    try {
        // Sprawdź czy histogramy mają wystarczającą statystykę
        Double_t dataIntegral = dataHist->Integral();
        if(dataIntegral < 10) {
            result.fitInfo = "Insufficient statistics in data histogram";
            return result;
        }
        
        // Sprawdź które histogramy MC mają wystarczającą statystykę
        std::vector<Bool_t> useInFit(mcHists.size(), false);
        std::vector<Int_t> fitIndexMap; // Mapowanie z indeksu w ficie na indeks w mcHists
        Int_t nHistsForFit = 0;
        
        for(size_t i = 0; i < mcHists.size(); ++i) {
            if(mcHists[i]->Integral() >= 1) {
                useInFit[i] = true;
                fitIndexMap.push_back(i);
                nHistsForFit++;
            }
        }
        
        if(nHistsForFit == 0) {
            result.fitInfo = "No MC histograms with sufficient statistics for fit";
            return result;
        }
        
        // Stwórz kopie histogramów MC do fitowania (tylko te z wystarczającą statystyką)
        TObjArray mcArray;
        std::vector<TH1*> mcCopies;
        Double_t totalMC = 0;
        
        for(size_t i = 0; i < mcHists.size(); ++i) {
            if(useInFit[i]) {
                TH1* copy = (TH1*)mcHists[i]->Clone(Form("%s_fitcopy_%zu", mcHists[i]->GetName(), i));
                mcCopies.push_back(copy);
                totalMC += copy->Integral();
            }
        }
        
        // Normalizuj MC do danych jako punkt startowy
        if(totalMC > 0) {
            Double_t scale = dataIntegral / totalMC;
            for(auto copy : mcCopies) {
                copy->Scale(scale);
                mcArray.Add(copy);
            }
        } else {
            result.fitInfo = "No events in MC histograms";
            return result;
        }
        
        // Inicjalizuj TFractionFitter
        TFractionFitter fitter(dataHist, &mcArray);
        
        // Ustaw zakres fitu jeśli określony
        if(constraints.fitRangeMin > 0 && constraints.fitRangeMax > constraints.fitRangeMin) {
            fitter.SetRangeX(constraints.fitRangeMin, constraints.fitRangeMax);
        }
        
        // Ustaw ograniczenia tylko dla histogramów używanych w ficie
        for(Int_t i = 0; i < nHistsForFit; ++i) {
            Int_t originalIndex = fitIndexMap[i];
            fitter.Constrain(i, constraints.lowerBounds[originalIndex], constraints.upperBounds[originalIndex]);
        }
        
        // Zwiększ liczbę iteracji
        ROOT::Math::MinimizerOptions::SetDefaultMaxIterations(500000);
        ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(500000);
        
        // Wykonaj fit z zabezpieczeniem stabilności
        result.status = EnsureFitStability(&fitter, constraints, nHistsForFit);
        
        if(result.status == 0) {
            result.converged = true;
            
            // Pobierz wyniki - tylko dla histogramów używanych w ficie
            for(Int_t i = 0; i < nHistsForFit; ++i) {
                Int_t originalIndex = fitIndexMap[i];
                fitter.GetResult(i, result.fractions[originalIndex], result.errors[originalIndex]);
            }
            
            // Ustaw frakcje na 0 dla histogramów nie używanych w ficie
            for(size_t i = 0; i < mcHists.size(); ++i) {
                if(!useInFit[i]) {
                    result.fractions[i] = 0.0;
                    result.errors[i] = 0.0;
                }
            }
            
            result.chi2 = fitter.GetChisquare();
            result.ndf = fitter.GetNDF();
            
            // Sprawdź jakość fitu
            if(result.ndf > 0) {
                Double_t chi2_ndf = result.chi2 / result.ndf;
                if(chi2_ndf > 10.0) {
                    result.fitInfo = Form("Poor fit quality: chi2/ndf = %.2f", chi2_ndf);
                } else if(chi2_ndf > 3.0) {
                    result.fitInfo = Form("Moderate fit quality: chi2/ndf = %.2f", chi2_ndf);
                } else {
                    result.fitInfo = Form("Good fit quality: chi2/ndf = %.2f", chi2_ndf);
                }
            }
            
        } else {
            result.fitInfo = Form("Fit failed with status: %d", result.status);
        }
        
        // Cleanup kopii
        for(auto copy : mcCopies) {
            delete copy;
        }
        
    } catch(const std::exception& e) {
        result.fitInfo = Form("Exception during fit: %s", e.what());
        result.status = -999;
    }
    
    return result;
}

Int_t HistManager::EnsureFitStability(TFractionFitter* fitter, const FitConstraints& constraints, 
                                     Int_t maxRetries) {
    Int_t status = -1;
    
    for(Int_t attempt = 0; attempt < maxRetries; ++attempt) {
        status = fitter->Fit();
        
        if(status == 0) {
            // Sprawdź czy wyniki są sensowne
            Bool_t resultsValid = true;
            Double_t sumFractions = 0.0;
            
            // Sprawdź wszystkie parametry fitu używając liczby histogramów w ficie
            // Liczba parametrów to liczba histogramów MC używanych w ficie
            Int_t nPars = maxRetries; // maxRetries jest teraz przekazywane jako nHistsForFit
            for(Int_t i = 0; i < nPars; ++i) {
                Double_t fraction, error;
                fitter->GetResult(i, fraction, error);
                
                // Sprawdź czy frakcje są w rozsądnych granicach
                if(fraction < -0.1 || fraction > 10.0 || error <= 0) {
                    resultsValid = false;
                    break;
                }
                sumFractions += fraction;
            }
            
            // Sprawdź czy suma frakcji jest rozsądna
            if(resultsValid && sumFractions > 0.1 && sumFractions < 10.0) {
                break; // Fit OK
            }
        }
        
        // Jeśli fit nie udał się, spróbuj ponownie (TFractionFitter ma wewnętrzną randomizację)
        if(attempt < 3) { // Maksymalnie 3 próby
            std::cout << "Fit attempt " << attempt + 1 << " failed, retrying..." << std::endl;
        }
    }
    
    return status;
}

void HistManager::UpdateMCHistograms(const std::vector<TH1*>& mcHists, TH1* sumHist, 
                                    const FitResult& fitResult) {
    if(!fitResult.converged || fitResult.status != 0) {
        return;
    }
    
    // Wyczyść histogram sumy
    sumHist->Reset();
    
    // Znajdź histogram danych dla tego zestawu
    Double_t dataIntegral = 0;
    for(const auto& pair : fData1D) {
        dataIntegral = pair.second->Integral();
        break; // Używamy pierwszego znalezionego (powinien być tylko jeden dla tego zestawu)
    }
    
    if(dataIntegral <= 0) {
        return;
    }
    
    // Przeskaluj histogramy MC i zaktualizuj sumę
    for(Int_t i = 0; i < fChannNum && i < static_cast<Int_t>(mcHists.size()) && i < static_cast<Int_t>(fitResult.fractions.size()); ++i) {
        if(fitResult.fractions[i] > 0) {
            // Histogram został użyty w ficie
            Double_t currentIntegral = mcHists[i]->Integral();
            if(currentIntegral > 0) {
                Double_t scaleFactor = (fitResult.fractions[i] * dataIntegral) / currentIntegral;
                mcHists[i]->Scale(scaleFactor);
            }
            
            // Dodaj do sumy
            sumHist->Add(mcHists[i]);
        } else {
            // Histogram nie został użyty w ficie (frakcja = 0)
            // Wyzeruj go żeby nie był widoczny na płótnie
            mcHists[i]->Scale(0);
        }
    }
}

HistManager::FitResult HistManager::GetFitResult(const TString& setName) const {
    auto it = fFitResults.find(setName);
    if(it != fFitResults.end()) {
        return it->second;
    }
    return FitResult(); // Pusty wynik jeśli nie znaleziono
}

Bool_t HistManager::HasFitResult(const TString& setName) const {
    return fFitResults.find(setName) != fFitResults.end();
}

Double_t HistManager::PerformSimpleScaling(const TString& setName) {
    auto histIt = fHists1D.find(setName);
    if(histIt == fHists1D.end()) {
        throw std::runtime_error("Histogram set not found: " + std::string(setName.Data()));
    }
    
    auto dataIt = fData1D.find(setName);
    if(dataIt == fData1D.end()) {
        throw std::runtime_error("No data histogram found for set: " + std::string(setName.Data()));
    }
    
    // Pobierz liczbę zdarzeń w danych i w sumie MC
    Double_t dataEntries = dataIt->second->Integral();
    Double_t mcSumEntries = histIt->second[0]->Integral(); // Indeks 0 to suma MC
    
    if(mcSumEntries <= 0) {
        std::cerr << "Warning: MC sum histogram is empty for set " << setName << std::endl;
        return 0.0;
    }
    
    // Oblicz czynnik skalujący
    Double_t scaleFactor = dataEntries / mcSumEntries;
    
    // Przygotuj wektor histogramów MC (bez sumy - indeks 0)
    std::vector<TH1*> mcHists;
    for(size_t i = 1; i < histIt->second.size(); ++i) {
        mcHists.push_back(histIt->second[i]);
    }
    
    // Zastosuj skalowanie
    ApplySimpleScaling(mcHists, histIt->second[0], scaleFactor);
    
    std::cout << "Simple scaling applied to set " << setName 
              << ": scale factor = " << scaleFactor 
              << " (Data: " << dataEntries << ", MC: " << mcSumEntries << ")" << std::endl;
    
    return scaleFactor;
}

void HistManager::ApplySimpleScaling(const std::vector<TH1*>& mcHists, TH1* sumHist, Double_t scaleFactor) {
    if(scaleFactor <= 0) {
        return;
    }
    
    // Przeskaluj wszystkie składowe MC
    for(auto hist : mcHists) {
        hist->Scale(scaleFactor);
    }
    
    // Przeskaluj sumę MC
    sumHist->Scale(scaleFactor);
}

// ===== IMPLEMENTACJA METOD ARRAY HISTOGRAMÓW =====

Bool_t HistManager::CreateHistArray1D(const ArrayConfig& config) {
    if(config.arraySize <= 0) {
        std::cerr << "ERROR: Invalid array size: " << config.arraySize << std::endl;
        return false;
    }
    
    // Sprawdź czy array już istnieje i usuń go
    if(fArrayHists1D.find(config.baseName) != fArrayHists1D.end()) {
        CleanupArraySet(config.baseName);
    }
    
    // Stwórz strukturę [index][mctruth+1] gdzie mctruth+1: 0=suma MC, 1-N=kanały MC  
    std::vector<std::vector<TH1D*>> arraySet(config.arraySize);
    std::vector<TH1D*> dataSet(config.arraySize);
    
    for(Int_t index = 0; index < config.arraySize; ++index) {
        // Określ nazwę i tytuł dla tego indeksu
        TString varName = config.baseName + Form("[%d]", index);
        TString varTitle = config.baseTitle + Form(" [%d]", index);
        
        // Użyj custom names/titles jeśli dostępne
        if(index < static_cast<Int_t>(config.varNames.size())) {
            varName = config.varNames[index];
        }
        if(index < static_cast<Int_t>(config.varTitles.size())) {
            varTitle = config.varTitles[index];
        }
        
        // Stwórz histogramy dla tego indeksu (suma MC + wszystkie kanały MC)
        std::vector<TH1D*> indexHists(fChannNum + 1); // +1 dla sumy MC (index 0)
        
        // Histogram sumy MC (index 0)
        TString sumHistName = Form("%s_sum_%d", config.baseName.Data(), index);
        indexHists[0] = new TH1D(sumHistName, varTitle, 
                                config.commonConfig.bins, 
                                config.commonConfig.xmin, 
                                config.commonConfig.xmax);
        indexHists[0]->GetXaxis()->SetTitle(config.commonConfig.xtitle);
        indexHists[0]->GetYaxis()->SetTitle(config.commonConfig.ytitle);
        ConfigureHistogram(indexHists[0], fSumColor, config.commonConfig.showStats);
        
        // Histogramy MC (indeksy 1 do fChannNum)
        for(Int_t ch = 0; ch < fChannNum; ++ch) {
            TString mcHistName = Form("%s_mc%d_%d", config.baseName.Data(), ch+1, index);
            indexHists[ch + 1] = new TH1D(mcHistName, varTitle,
                                         config.commonConfig.bins,
                                         config.commonConfig.xmin,
                                         config.commonConfig.xmax);
            indexHists[ch + 1]->GetXaxis()->SetTitle(config.commonConfig.xtitle);
            indexHists[ch + 1]->GetYaxis()->SetTitle(config.commonConfig.ytitle);
            ConfigureHistogram(indexHists[ch + 1], fChannColors[ch], config.commonConfig.showStats);
        }
        
        arraySet[index] = indexHists;
        
        // Histogram danych
        TString dataHistName = Form("%s_data_%d", config.baseName.Data(), index);
        dataSet[index] = new TH1D(dataHistName, varTitle, 
                                 config.commonConfig.bins, 
                                 config.commonConfig.xmin, 
                                 config.commonConfig.xmax);
        dataSet[index]->GetXaxis()->SetTitle(config.commonConfig.xtitle);
        dataSet[index]->GetYaxis()->SetTitle(config.commonConfig.ytitle);
        dataSet[index]->SetMarkerStyle(fDataStyle);
        dataSet[index]->SetMarkerColor(fDataColor);
        dataSet[index]->SetLineColor(fDataColor);
        dataSet[index]->SetStats(config.commonConfig.showStats);
    }
    
    fArrayHists1D[config.baseName] = arraySet;
    fArrayData1D[config.baseName] = dataSet;
    fArrayConfigs[config.baseName] = config;
    
    return true;
}

Bool_t HistManager::FillArray1D(const TString& baseName, Int_t index, Int_t mctruth, Double_t value, Double_t weight) {
    auto it = fArrayHists1D.find(baseName);
    if(it == fArrayHists1D.end()) {
        std::cerr << "ERROR: Array histogram set '" << baseName << "' not found!" << std::endl;
        return false;
    }
    
    if(index < 0 || index >= static_cast<Int_t>(it->second.size())) {
        std::cerr << "ERROR: Array index " << index << " out of range [0," << it->second.size()-1 << "]!" << std::endl;
        return false;
    }
    
    if(mctruth < 1 || mctruth > fChannNum) {
        return false; // Invalid mctruth, silently ignore
    }
    
    // Wypełnij sumę MC (index 0)
    it->second[index][0]->Fill(value, weight);
    // Wypełnij konkretny kanał MC (index mctruth)
    it->second[index][mctruth]->Fill(value, weight);
    
    return true;
}

Bool_t HistManager::FillArrayData1D(const TString& baseName, Int_t index, Double_t value, Double_t weight) {
    auto it = fArrayData1D.find(baseName);
    if(it == fArrayData1D.end()) {
        std::cerr << "ERROR: Array data histogram set '" << baseName << "' not found!" << std::endl;
        return false;
    }
    
    if(index < 0 || index >= static_cast<Int_t>(it->second.size())) {
        std::cerr << "ERROR: Array index " << index << " out of range [0," << it->second.size()-1 << "]!" << std::endl;
        return false;
    }
    
    it->second[index]->Fill(value, weight);
    return true;
}

Bool_t HistManager::FillArrayAll1D(const TString& baseName, Int_t mctruth, const std::vector<Double_t>& values, Double_t weight) {
    auto it = fArrayHists1D.find(baseName);
    if(it == fArrayHists1D.end()) {
        std::cerr << "ERROR: Array histogram set '" << baseName << "' not found!" << std::endl;
        return false;
    }
    
    if(static_cast<Int_t>(values.size()) != static_cast<Int_t>(it->second.size())) {
        std::cerr << "ERROR: Values vector size (" << values.size() << ") doesn't match array size (" << it->second.size() << ")!" << std::endl;
        return false;
    }
    
    for(Int_t index = 0; index < static_cast<Int_t>(values.size()); ++index) {
        FillArray1D(baseName, index, mctruth, values[index], weight);
    }
    
    return true;
}

Bool_t HistManager::FillArrayAllData1D(const TString& baseName, const std::vector<Double_t>& values, Double_t weight) {
    auto it = fArrayData1D.find(baseName);
    if(it == fArrayData1D.end()) {
        std::cerr << "ERROR: Array data histogram set '" << baseName << "' not found!" << std::endl;
        return false;
    }
    
    if(static_cast<Int_t>(values.size()) != static_cast<Int_t>(it->second.size())) {
        std::cerr << "ERROR: Values vector size (" << values.size() << ") doesn't match array size (" << it->second.size() << ")!" << std::endl;
        return false;
    }
    
    for(Int_t index = 0; index < static_cast<Int_t>(values.size()); ++index) {
        FillArrayData1D(baseName, index, values[index], weight);
    }
    
    return true;
}

TCanvas* HistManager::DrawArray1D(const TString& baseName, Bool_t drawData, const TString& canvasName, 
                                  Int_t nCols, Int_t nRows) {
    auto it = fArrayHists1D.find(baseName);
    if(it == fArrayHists1D.end()) {
        std::cerr << "ERROR: Array histogram set '" << baseName << "' not found!" << std::endl;
        return nullptr;
    }
    
    Int_t nPanels = it->second.size();
    if(nPanels == 0) {
        std::cerr << "ERROR: Empty array histogram set!" << std::endl;
        return nullptr;
    }
    
    // Oblicz layout kanwy
    CalculateCanvasLayout(nPanels, nCols, nRows);
    
    // Stwórz kanwę
    TString finalCanvasName = canvasName.IsNull() ? Form("c_%s_array", baseName.Data()) : canvasName;
    TCanvas* canvas = new TCanvas(finalCanvasName, baseName + " Array", 
                                 200 * nCols, 150 * nRows);
    canvas->Divide(nCols, nRows);
    
    // Pobierz konfigurację dla log scale
    auto configIt = fArrayConfigs.find(baseName);
    Bool_t logY = (configIt != fArrayConfigs.end()) ? configIt->second.commonConfig.logy : false;
    
    // Rysuj każdy panel
    for(Int_t index = 0; index < nPanels; ++index) {
        canvas->cd(index + 1);
        
        if(logY) {
            gPad->SetLogy();
        }
        
        // Znajdź maksimum dla skalowania Y
        Double_t maxY = 0;
        for(auto* hist : it->second[index]) {
            if(hist && hist->GetEntries() > 0) {
                maxY = TMath::Max(maxY, hist->GetMaximum());
            }
        }
        
        // Sprawdź maksimum danych jeśli rysujemy
        if(drawData) {
            auto dataIt = fArrayData1D.find(baseName);
            if(dataIt != fArrayData1D.end() && dataIt->second[index] && dataIt->second[index]->GetEntries() > 0) {
                maxY = TMath::Max(maxY, dataIt->second[index]->GetMaximum());
            }
        }
        
        maxY = logY ? maxY * 5 : maxY * 1.1;
        
        bool first = true;
        
        // Rysuj histogramy MC (indeksy 1 do fChannNum)
        for(Int_t ch = 1; ch <= fChannNum; ++ch) {
            TH1D* mcHist = it->second[index][ch];
            if(mcHist && mcHist->GetEntries() > 0) {
                mcHist->GetYaxis()->SetRangeUser(logY ? 0.1 : 0, maxY);
                mcHist->Draw(first ? "HIST" : "HIST SAME");
                first = false;
            }
        }
        
        // Rysuj sumę MC (index 0)
        TH1D* sumHist = it->second[index][0];
        if(sumHist && sumHist->GetEntries() > 0) {
            sumHist->GetYaxis()->SetRangeUser(logY ? 0.1 : 0, maxY);
            sumHist->Draw(first ? "HIST" : "HIST SAME");
            first = false;
        }
        
        // Rysuj dane na końcu
        if(drawData) {
            auto dataIt = fArrayData1D.find(baseName);
            if(dataIt != fArrayData1D.end()) {
                TH1D* dataHist = dataIt->second[index];
                if(dataHist && dataHist->GetEntries() > 0) {
                    if(first) {
                        dataHist->GetYaxis()->SetRangeUser(logY ? 0.1 : 0, maxY);
                        dataHist->Draw("PE1");
                    } else {
                        dataHist->Draw("PE1 SAME");
                    }
                }
            }
        }
        
        // Dodaj tytuł panelu
        if(configIt != fArrayConfigs.end() && 
           index < static_cast<Int_t>(configIt->second.varTitles.size())) {
            gPad->SetTitle(configIt->second.varTitles[index]);
        } else {
            gPad->SetTitle(Form("%s[%d]", baseName.Data(), index));
        }
    }
    
    // Zapisz kanwę
    if(fCanvases.find(baseName) == fCanvases.end()) {
        fCanvases[baseName] = std::vector<TCanvas*>();
    }
    fCanvases[baseName].push_back(canvas);
    
    return canvas;
}

TH1D* HistManager::GetArrayHist1D(const TString& baseName, Int_t index, Int_t mctruth) {
    if(mctruth == -1) {
        // Dane
        auto dataIt = fArrayData1D.find(baseName);
        if(dataIt == fArrayData1D.end() || index < 0 || index >= static_cast<Int_t>(dataIt->second.size())) {
            return nullptr;
        }
        return dataIt->second[index];
    } else {
        // MC
        auto it = fArrayHists1D.find(baseName);
        if(it == fArrayHists1D.end() || index < 0 || index >= static_cast<Int_t>(it->second.size())) {
            return nullptr;
        }
        
        // mctruth: 0=suma MC, 1-N=kanały MC
        if(mctruth < 0 || mctruth > fChannNum) {
            return nullptr;
        }
        
        return it->second[index][mctruth];
    }
}

void HistManager::CalculateCanvasLayout(Int_t nPanels, Int_t& nCols, Int_t& nRows) {
    if(nCols > 0 && nRows > 0) {
        return; // Użytkownik podał wymiary
    }
    
    if(nPanels <= 1) {
        nCols = nRows = 1;
    } else if(nPanels <= 2) {
        nCols = 2; nRows = 1;
    } else if(nPanels <= 4) {
        nCols = nRows = 2;
    } else if(nPanels <= 6) {
        nCols = 3; nRows = 2;
    } else if(nPanels <= 9) {
        nCols = nRows = 3;
    } else if(nPanels <= 12) {
        nCols = 4; nRows = 3;
    } else if(nPanels <= 16) {
        nCols = nRows = 4;
    } else {
        // Dla większej liczby paneli, oblicz optymalny rozkład
        nCols = static_cast<Int_t>(TMath::Ceil(TMath::Sqrt(nPanels)));
        nRows = static_cast<Int_t>(TMath::Ceil(static_cast<Double_t>(nPanels) / nCols));
    }
}

void HistManager::CleanupArraySet(const TString& baseName) {
    // Cleanup MC histograms
    auto itMC = fArrayHists1D.find(baseName);
    if(itMC != fArrayHists1D.end()) {
        for(auto& indexVec : itMC->second) {
            for(auto* hist : indexVec) {
                delete hist;
            }
        }
        fArrayHists1D.erase(itMC);
    }
    
    // Cleanup data histograms  
    auto itData = fArrayData1D.find(baseName);
    if(itData != fArrayData1D.end()) {
        for(auto* hist : itData->second) {
            delete hist;
        }
        fArrayData1D.erase(itData);
    }
    
    // Remove config
    fArrayConfigs.erase(baseName);
}
