#include "HistManager.h"
#include <TMath.h>
#include <TStyle.h>
#include <TPaveText.h>
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

    // Calculate chi2 if enabled and data is available
    Chi2Result chi2Result;
    Bool_t hasChi2 = false;
    if(drawData && fCalculateDataMCChi2) {
        auto dataIt = fData1D.find(setName);
        if(dataIt != fData1D.end()) {
            try {
                chi2Result = CalculateDataMCChi2(setName);
                hasChi2 = true;
            } catch(const std::exception& e) {
                std::cerr << "Error during chi2 calculation: " << e.what() << std::endl;
            }
        }
    }

    // Add legend with data and chi2 if present
    TH1* dataHist = nullptr;
    if(drawData && fData1D.find(setName) != fData1D.end()) {
        dataHist = fData1D[setName];
    }
    TLegend* legend = CreateLegend(histIt->second, dataHist, hasChi2 ? &chi2Result : nullptr);
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
                                    800, 600);
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

TLegend* HistManager::CreateLegend(const std::vector<TH1*>& hists, const TH1* dataHist, const Chi2Result* chi2Result) {
    // Oblicz rozmiar legendy w zależności od zawartości
    Int_t nEntries = hists.size(); // MC sum + składowe
    if(dataHist) nEntries++; // Data
    if(chi2Result && chi2Result->nBinsUsed > 0) nEntries += 3; // Chi2 info (3 linie)
    
    // Dostosuj rozmiar legendy
    Double_t legendHeight = TMath::Max(0.15, TMath::Min(0.5, nEntries * 0.04));
    Double_t legendTop = 0.9;
    Double_t legendBottom = legendTop - legendHeight;
    
    TLegend* legend = new TLegend(0.65, legendBottom, 0.9, legendTop);
    legend->SetBorderSize(1);
    legend->SetFillStyle(0); // Przezroczyste tło
    legend->SetTextSize(0.03);
    
    // Dodaj najpierw sumę MC
    legend->AddEntry(hists[0], "MC sum", "l");
    
    // Dodaj pozostałe składowe
    for(size_t i = 1; i < hists.size(); ++i) {
        legend->AddEntry(hists[i], fChannelNames[i-1], "l");
    }
    
    // Dodaj dane jeśli są dostępne
    if(dataHist) {
        legend->AddEntry(dataHist, "Data", "pe");
    }
    
    // Dodaj informacje o chi-squared jeśli są dostępne
    if(chi2Result && chi2Result->nBinsUsed > 0) {
        // Dodaj separator
        legend->AddEntry((TObject*)0, "", "");
        
        // Dodaj chi2 informacje
        legend->AddEntry((TObject*)0, Form("#chi^{2}/NDF = %.1f/%d", chi2Result->chi2, chi2Result->ndf), "");
        legend->AddEntry((TObject*)0, Form("= %.2f", chi2Result->chi2_ndf), "");
        legend->AddEntry((TObject*)0, Form("p-value = %.3f", chi2Result->pValue), "");
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

// ==================== CHI2 DATA vs MC SUM IMPLEMENTATION ====================

HistManager::Chi2Result HistManager::CalculateDataMCChi2(const TString& setName, Double_t minBinContent) {
    auto histIt = fHists1D.find(setName);
    if(histIt == fHists1D.end()) {
        throw std::runtime_error("Histogram set not found: " + std::string(setName.Data()));
    }
    
    auto dataIt = fData1D.find(setName);
    if(dataIt == fData1D.end()) {
        throw std::runtime_error("No data histogram found for set: " + std::string(setName.Data()));
    }
    
    // Histogram sumy MC jest pod indeksem 0
    TH1* mcSumHist = histIt->second[0];
    TH1* dataHist = dataIt->second;
    
    // Wykonaj obliczenie chi-kwadrat
    Chi2Result result = DoDataMCChi2Calculation(dataHist, mcSumHist, minBinContent);
    
    // Zapisz wyniki
    fChi2Results[setName] = result;
    fLastChi2Result = result;
    
    // Wypisz wyniki
    std::cout << "=== Chi2 Data vs MC Sum for set: " << setName << " ===" << std::endl;
    std::cout << "Chi2: " << result.chi2 << std::endl;
    std::cout << "NDF: " << result.ndf << std::endl;
    std::cout << "Chi2/NDF: " << result.chi2_ndf << std::endl;
    std::cout << "p-value: " << result.pValue << std::endl;
    std::cout << "Bins used: " << result.nBinsUsed << std::endl;
    if(!result.comparisonInfo.IsNull()) {
        std::cout << "Info: " << result.comparisonInfo << std::endl;
    }
    std::cout << "=================================================" << std::endl;
    
    return result;
}

HistManager::Chi2Result HistManager::DoDataMCChi2Calculation(TH1* dataHist, TH1* mcSumHist, Double_t minBinContent) {
    Chi2Result result;
    
    if(!dataHist || !mcSumHist) {
        result.comparisonInfo = "Invalid histogram pointers";
        return result;
    }
    
    // Sprawdź czy histogramy mają ten sam binning
    if(dataHist->GetNbinsX() != mcSumHist->GetNbinsX()) {
        result.comparisonInfo = "Histograms have different number of bins";
        return result;
    }
    
    Int_t nBins = dataHist->GetNbinsX();
    Double_t chi2 = 0.0;
    Int_t nValidBins = 0;
    
    // Oblicz chi-kwadrat bin po binie
    for(Int_t bin = 1; bin <= nBins; ++bin) {
        Double_t dataContent = dataHist->GetBinContent(bin);
        Double_t mcContent = mcSumHist->GetBinContent(bin);
        
        // Pomiń biny z małą statystyką
        if(dataContent < minBinContent && mcContent < minBinContent) {
            continue;
        }
        
        // Użyj błędu Poissona dla danych
        Double_t dataError = dataHist->GetBinError(bin);
        if(dataError <= 0) {
            dataError = TMath::Sqrt(TMath::Max(1.0, dataContent)); // Błąd Poissona
        }
        
        // Oblicz wkład do chi-kwadrat
        if(dataError > 0) {
            Double_t diff = dataContent - mcContent;
            chi2 += (diff * diff) / (dataError * dataError);
            nValidBins++;
        }
    }
    
    // Wypełnij wyniki
    result.chi2 = chi2;
    result.nBinsUsed = nValidBins;
    result.ndf = TMath::Max(1, nValidBins - 1); // NDF = liczba binów - 1 (brak parametrów fitu)
    result.chi2_ndf = (result.ndf > 0) ? result.chi2 / result.ndf : 0.0;
    
    // Oblicz p-value używając funkcji gamma niekompletnej
    if(result.ndf > 0 && result.chi2 >= 0) {
        result.pValue = TMath::Prob(result.chi2, result.ndf);
    }
    
    // Dodaj informacje o jakości porównania
    if(nValidBins < 3) {
        result.comparisonInfo = "Too few bins for reliable chi2 calculation";
    } else if(result.chi2_ndf > 5.0) {
        result.comparisonInfo = "Poor agreement between data and MC";
    } else if(result.chi2_ndf > 2.0) {
        result.comparisonInfo = "Moderate agreement between data and MC";
    } else {
        result.comparisonInfo = "Good agreement between data and MC";
    }
    
    return result;
}
