#include "HistManager.h"
#include <TMath.h>
#include <TStyle.h>
#include <Math/MinimizerOptions.h>
#include <TPaletteAxis.h>
#include <TFitResult.h>
#include <algorithm>
#include <stdexcept>
#include <iostream>

HistManager::HistManager(Int_t channNum, const Color_t* channColors, 
                                 const std::vector<TString>& channelNames,
                                 Int_t dataStyle, Int_t dataColor, Float_t dataSize,Int_t sumColor) 
    : fChannNum(channNum), fDataStyle(dataStyle), fDataColor(dataColor), fDataSize(dataSize), fSumColor(sumColor) {
    
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
    
    // Tworzenie histogramu sumy jako pierwszy element wektora (konsystentnie z 1D)
    TString sumHistName = Form("%s_sum", config.name.Data());
    TH2D* sumHist = new TH2D(sumHistName, config.title,
                            config.bins, config.xmin, config.xmax,
                            config.binsy, config.ymin, config.ymax);
    sumHist->GetXaxis()->SetTitle(config.xtitle);
    sumHist->GetYaxis()->SetTitle(config.ytitle);
    ConfigureHistogram(sumHist, fSumColor, config.showStats);
    hists.push_back(sumHist);
    
    // Tworzenie histogramów dla poszczególnych kanałów
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

    // Fill suma MC (indeks 0) - konsystentnie z Fill1D
    it->second[0]->Fill(x, y, weight);
    // Fill pojedyncza składowa (przesunięte o 1 przez sumę)
    it->second[mctruth]->Fill(x, y, weight);
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

    // Check for Triple Gaussian fit results and draw them if available
    for(size_t i = 0; i < histIt->second.size(); ++i) {
        Int_t mctruth = (i == 0) ? 0 : static_cast<Int_t>(i);  // 0 for sum, i for channel
        TString key = GenerateTripleGaussKey(setName, mctruth, -1);
        auto fitIt = f3GaussFitResults.find(key);
        if(fitIt != f3GaussFitResults.end() && fitIt->second.converged) {
            DrawTripleGaussFit(histIt->second[i], fitIt->second, false);
        }
    }

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
            dataIt->second->SetMarkerSize(fDataSize);
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

    auto configIt = fConfigs2D.find(setName);
    if(configIt == fConfigs2D.end()) {
        throw std::runtime_error("Configuration not found for set: " + std::string(setName.Data()));
    }

    std::vector<TCanvas*> canvases;
    
    // Stwórz canvas z MC Sum i danymi obok siebie (jeśli dane są dostępne)
    Bool_t hasData = false;
    auto dataIt = fData2D.find(setName);
    if(drawData && dataIt != fData2D.end() && dataIt->second) {
        hasData = true;
    }
    
    if(hasData) {
        // Canvas z MC Sum i danymi obok siebie
        TString summaryCanvasName = Form("c_%s_summary", setName.Data());
        TString summaryCanvasTitle = Form("%s - MC Sum & Data", configIt->second.title.Data());
        
        TCanvas* summaryCanvas = new TCanvas(summaryCanvasName, summaryCanvasTitle, 1600, 600);
        summaryCanvas->Divide(2, 1);
        
        // Panel 1: MC Sum
        summaryCanvas->cd(1);
        gPad->SetLeftMargin(0.15);
        gPad->SetBottomMargin(0.15);
        gPad->SetRightMargin(0.18);
        gPad->SetTopMargin(0.12);
        
        if(histIt->second[0] && histIt->second[0]->GetEntries() > 0) {
            histIt->second[0]->Draw(drawOpt);
            
            // Aktualizuj pad aby uzyskać dostęp do primitives
            gPad->Update();
            
            // Znajdź i dostosuj pozycję color bar (palette)
            TPaletteAxis* palette = (TPaletteAxis*)histIt->second[0]->GetListOfFunctions()->FindObject("palette");
            if(palette) {
                palette->SetX1NDC(0.83);  // Lewa krawędź color bar - obok prawej osi Y
                palette->SetX2NDC(0.88);  // Prawa krawędź color bar
                palette->SetY1NDC(0.15);  // Dolna krawędź = dolna oś X
                palette->SetY2NDC(0.88);  // Górna krawędź = górna oś X
            }
            
            TPaveText* sumInfo = new TPaveText(0.02, 0.85, 0.35, 0.98, "NDC");
            sumInfo->SetFillColor(kOrange-10);
            sumInfo->SetBorderSize(1);
            sumInfo->AddText("MC Sum (All Channels)");
            sumInfo->AddText(Form("mctruth: All (1-%d)", fChannNum));
            sumInfo->AddText(Form("Entries: %.0f", histIt->second[0]->GetEntries()));
            sumInfo->Draw();
        } else {
            histIt->second[0]->Draw(drawOpt);
            
            // Aktualizuj pad aby uzyskać dostęp do primitives
            gPad->Update();
            
            // Znajdź i dostosuj pozycję color bar (palette) nawet dla pustego histogramu
            TPaletteAxis* palette = (TPaletteAxis*)histIt->second[0]->GetListOfFunctions()->FindObject("palette");
            if(palette) {
                palette->SetX1NDC(0.83);  // Lewa krawędź color bar - obok prawej osi Y
                palette->SetX2NDC(0.88);  // Prawa krawędź color bar
                palette->SetY1NDC(0.15);  // Dolna krawędź = dolna oś X
                palette->SetY2NDC(0.88);  // Górna krawędź = górna oś X
            }
            
            TPaveText* emptyInfo = new TPaveText(0.3, 0.4, 0.7, 0.6, "NDC");
            emptyInfo->SetFillColor(kYellow-10);
            emptyInfo->SetBorderSize(1);
            emptyInfo->AddText("Empty MC Sum");
            emptyInfo->Draw();
        }
        
        // Panel 2: Data
        summaryCanvas->cd(2);
        gPad->SetLeftMargin(0.15);
        gPad->SetBottomMargin(0.15);
        gPad->SetRightMargin(0.18);
        gPad->SetTopMargin(0.12);
        
        dataIt->second->Draw(drawOpt);
        
        // Aktualizuj pad aby uzyskać dostęp do primitives
        gPad->Update();
        
        // Znajdź i dostosuj pozycję color bar (palette) dla danych
        TPaletteAxis* dataPalette = (TPaletteAxis*)dataIt->second->GetListOfFunctions()->FindObject("palette");
        if(dataPalette) {
            dataPalette->SetX1NDC(0.83);  // Lewa krawędź color bar - obok prawej osi Y
            dataPalette->SetX2NDC(0.88);  // Prawa krawędź color bar
            dataPalette->SetY1NDC(0.15);  // Dolna krawędź = dolna oś X
            dataPalette->SetY2NDC(0.88);  // Górna krawędź = górna oś X
        }
        
        TPaveText* dataInfo = new TPaveText(0.02, 0.85, 0.35, 0.98, "NDC");
        dataInfo->SetFillColor(kCyan-10);
        dataInfo->SetBorderSize(1);
        dataInfo->AddText("Experimental Data");
        dataInfo->AddText("mctruth: Real Data");
        dataInfo->AddText(Form("Entries: %.0f", dataIt->second->GetEntries()));
        dataInfo->Draw();
        
        summaryCanvas->Update();
        canvases.push_back(summaryCanvas);
    } else {
        // Tylko MC Sum jeśli brak danych
        TString sumCanvasName = Form("c_%s_sum", setName.Data());
        TString sumCanvasTitle = Form("%s - MC Sum", configIt->second.title.Data());
        
        TCanvas* sumCanvas = new TCanvas(sumCanvasName, sumCanvasTitle, 800, 600);
        sumCanvas->SetLeftMargin(0.12);
        sumCanvas->SetBottomMargin(0.12);
        sumCanvas->SetRightMargin(0.15);
        sumCanvas->SetTopMargin(0.1);

        if(histIt->second[0] && histIt->second[0]->GetEntries() > 0) {
            histIt->second[0]->Draw(drawOpt);
            
            // Aktualizuj pad aby uzyskać dostęp do primitives
            gPad->Update();
            
            // Znajdź i dostosuj pozycję color bar (palette)
            TPaletteAxis* palette = (TPaletteAxis*)histIt->second[0]->GetListOfFunctions()->FindObject("palette");
            if(palette) {
                palette->SetX1NDC(0.83);  // Lewa krawędź color bar - obok prawej osi Y
                palette->SetX2NDC(0.88);  // Prawa krawędź color bar
                palette->SetY1NDC(0.12);  // Dolna krawędź = dolna oś X
                palette->SetY2NDC(0.90);  // Górna krawędź = górna oś X
            }
            
            TPaveText* sumInfo = new TPaveText(0.02, 0.85, 0.35, 0.98, "NDC");
            sumInfo->SetFillColor(kOrange-10);
            sumInfo->SetBorderSize(1);
            sumInfo->AddText("MC Sum (All Channels)");
            sumInfo->AddText(Form("mctruth: All (1-%d)", fChannNum));
            sumInfo->AddText(Form("Entries: %.0f", histIt->second[0]->GetEntries()));
            sumInfo->Draw();
        } else {
            histIt->second[0]->Draw(drawOpt);
            
            // Aktualizuj pad aby uzyskać dostęp do primitives
            gPad->Update();
            
            // Znajdź i dostosuj pozycję color bar (palette) nawet dla pustego histogramu
            TPaletteAxis* palette = (TPaletteAxis*)histIt->second[0]->GetListOfFunctions()->FindObject("palette");
            if(palette) {
                palette->SetX1NDC(0.83);  // Lewa krawędź color bar - obok prawej osi Y
                palette->SetX2NDC(0.88);  // Prawa krawędź color bar
                palette->SetY1NDC(0.12);  // Dolna krawędź = dolna oś X
                palette->SetY2NDC(0.90);  // Górna krawędź = górna oś X
            }
            
            TPaveText* emptyInfo = new TPaveText(0.3, 0.4, 0.7, 0.6, "NDC");
            emptyInfo->SetFillColor(kYellow-10);
            emptyInfo->SetBorderSize(1);
            emptyInfo->AddText("Empty MC Sum");
            emptyInfo->Draw();
        }
        
        sumCanvas->Update();
        canvases.push_back(sumCanvas);
    }
    
    // Rysuj histogramy dla każdego kanału MC (indeksy 1 do fChannNum)
    for(Int_t i = 0; i < fChannNum; ++i) {
        TString canvasName = Form("c_%s_ch%d", setName.Data(), i+1);
        TString canvasTitle = Form("%s - %s", configIt->second.title.Data(), fChannelNames[i].Data());
        
        TCanvas* canvas = new TCanvas(canvasName, canvasTitle, 800, 600);
        canvas->SetLeftMargin(0.15);
        canvas->SetBottomMargin(0.15);
        canvas->SetRightMargin(0.18);
        canvas->SetTopMargin(0.12);

        // Rysuj histogram MC (indeks i+1 przez histogram sumy na indeksie 0)
        if(histIt->second[i+1] && histIt->second[i+1]->GetEntries() > 0) {
            histIt->second[i+1]->Draw(drawOpt);
            
            // Aktualizuj pad aby uzyskać dostęp do primitives
            gPad->Update();
            
            // Znajdź i dostosuj pozycję color bar (palette)
            TPaletteAxis* palette = (TPaletteAxis*)histIt->second[i+1]->GetListOfFunctions()->FindObject("palette");
            if(palette) {
                palette->SetX1NDC(0.83);  // Lewa krawędź color bar - obok prawej osi Y
                palette->SetX2NDC(0.88);  // Prawa krawędź color bar
                palette->SetY1NDC(0.15);  // Dolna krawędź = dolna oś X
                palette->SetY2NDC(0.88);  // Górna krawędź = górna oś X
            }
            
            TPaveText* channelInfo = new TPaveText(0.02, 0.85, 0.35, 0.98, "NDC");
            channelInfo->SetFillColor(kWhite);
            channelInfo->SetBorderSize(1);
            channelInfo->AddText(Form("Channel: %s", fChannelNames[i].Data()));
            channelInfo->AddText(Form("mctruth: %d", i+1));
            channelInfo->AddText(Form("Entries: %.0f", histIt->second[i+1]->GetEntries()));
            channelInfo->Draw();
        } else {
            histIt->second[i+1]->Draw(drawOpt);
            
            // Aktualizuj pad aby uzyskać dostęp do primitives
            gPad->Update();
            
            // Znajdź i dostosuj pozycję color bar (palette) nawet dla pustego histogramu
            TPaletteAxis* palette = (TPaletteAxis*)histIt->second[i+1]->GetListOfFunctions()->FindObject("palette");
            if(palette) {
                palette->SetX1NDC(0.83);  // Lewa krawędź color bar - obok prawej osi Y
                palette->SetX2NDC(0.88);  // Prawa krawędź color bar
                palette->SetY1NDC(0.15);  // Dolna krawędź = dolna oś X
                palette->SetY2NDC(0.88);  // Górna krawędź = górna oś X
            }
            
            TPaveText* emptyInfo = new TPaveText(0.3, 0.4, 0.7, 0.6, "NDC");
            emptyInfo->SetFillColor(kYellow-10);
            emptyInfo->SetBorderSize(1);
            emptyInfo->AddText("Empty Histogram");
            emptyInfo->AddText(Form("Channel: %s", fChannelNames[i].Data()));
            emptyInfo->AddText(Form("mctruth: %d", i+1));
            emptyInfo->Draw();
        }
        
        canvas->Update();
        canvases.push_back(canvas);
    }
    
    // Stwórz canvas porównawczy ze wszystkimi kanałami (włącznie z sumą MC i danymi)
    Int_t totalHists = fChannNum + 1; // +1 dla sumy MC
    if(hasData) totalHists += 1; // +1 dla danych
    
    Int_t nCols, nRows;
    
    if(totalHists <= 4) {
        nCols = 2; nRows = 2;
    } else if(totalHists <= 6) {
        nCols = 3; nRows = 2;
    } else if(totalHists <= 9) {
        nCols = 3; nRows = 3;
    } else {
        nCols = 4; nRows = TMath::Ceil(static_cast<Double_t>(totalHists) / 4.0);
    }
    
    TString compareCanvasName = Form("c_%s_compare", setName.Data());
    TString compareCanvasTitle = Form("%s - All Histograms", configIt->second.title.Data());
    
    TCanvas* compareCanvas = new TCanvas(compareCanvasName, compareCanvasTitle, 
                                       200 * nCols, 160 * nRows);
    compareCanvas->Divide(nCols, nRows);
    
    Int_t padIndex = 1;
    
    // Panel 1: Suma MC
    compareCanvas->cd(padIndex++);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);
    gPad->SetRightMargin(0.18);
    gPad->SetTopMargin(0.15);
    
    if(histIt->second[0] && histIt->second[0]->GetEntries() > 0) {
        histIt->second[0]->Draw("COLZ");
        
        // Aktualizuj pad aby uzyskać dostęp do primitives
        gPad->Update();
        
        // Znajdź i dostosuj pozycję color bar (palette)
        TPaletteAxis* palette = (TPaletteAxis*)histIt->second[0]->GetListOfFunctions()->FindObject("palette");
        if(palette) {
            palette->SetX1NDC(0.83);  // Lewa krawędź color bar - obok prawej osi Y
            palette->SetX2NDC(0.88);  // Prawa krawędź color bar
            palette->SetY1NDC(0.15);  // Dolna krawędź = dolna oś X
            palette->SetY2NDC(0.85);  // Górna krawędź = górna oś X
        }
        
        // Dodaj tekst z informacją o MC Sum - wycentrowany względem górnej osi X
        TPaveText* sumLabel = new TPaveText(0.3, 0.88, 0.7, 0.96, "NDC");
        sumLabel->SetFillColor(kOrange-10);
        sumLabel->SetBorderSize(1);
        sumLabel->SetTextSize(0.06);
        sumLabel->AddText("MC Sum");
        sumLabel->Draw();
    } else {
        histIt->second[0]->Draw("COLZ");
        
        // Aktualizuj pad aby uzyskać dostęp do primitives
        gPad->Update();
        
        // Znajdź i dostosuj pozycję color bar (palette) nawet dla pustego histogramu
        TPaletteAxis* palette = (TPaletteAxis*)histIt->second[0]->GetListOfFunctions()->FindObject("palette");
        if(palette) {
            palette->SetX1NDC(0.83);  // Lewa krawędź color bar - obok prawej osi Y
            palette->SetX2NDC(0.88);  // Prawa krawędź color bar
            palette->SetY1NDC(0.15);  // Dolna krawędź = dolna oś X
            palette->SetY2NDC(0.85);  // Górna krawędź = górna oś X
        }
        
        // Dodaj tekst z informacją o MC Sum nawet gdy pusty - wycentrowany względem górnej osi X
        TPaveText* sumLabel = new TPaveText(0.3, 0.88, 0.7, 0.96, "NDC");
        sumLabel->SetFillColor(kOrange-10);
        sumLabel->SetBorderSize(1);
        sumLabel->SetTextSize(0.06);
        sumLabel->AddText("MC Sum");
        sumLabel->Draw();
    }
    
    // Panel dla danych (jeśli dostępne)
    if(hasData) {
        compareCanvas->cd(padIndex++);
        gPad->SetLeftMargin(0.15);
        gPad->SetBottomMargin(0.15);
        gPad->SetRightMargin(0.18);
        gPad->SetTopMargin(0.15);
        
        dataIt->second->Draw("COLZ");
        
        // Aktualizuj pad aby uzyskać dostęp do primitives
        gPad->Update();
        
        // Znajdź i dostosuj pozycję color bar (palette) dla danych
        TPaletteAxis* dataPalette = (TPaletteAxis*)dataIt->second->GetListOfFunctions()->FindObject("palette");
        if(dataPalette) {
            dataPalette->SetX1NDC(0.83);  // Lewa krawędź color bar - obok prawej osi Y
            dataPalette->SetX2NDC(0.88);  // Prawa krawędź color bar
            dataPalette->SetY1NDC(0.15);  // Dolna krawędź = dolna oś X
            dataPalette->SetY2NDC(0.85);  // Górna krawędź = górna oś X
        }
        
        TPaveText* dataLabel = new TPaveText(0.3, 0.88, 0.7, 0.96, "NDC");
        dataLabel->SetFillColor(kCyan-10);
        dataLabel->SetBorderSize(1);
        dataLabel->SetTextSize(0.06);
        dataLabel->AddText("Data");
        dataLabel->Draw();
    }
    
    // Panele dla kanałów MC
    for(Int_t i = 0; i < fChannNum; ++i) {
        compareCanvas->cd(padIndex++);
        gPad->SetLeftMargin(0.15);
        gPad->SetBottomMargin(0.15);
        gPad->SetRightMargin(0.18);
        gPad->SetTopMargin(0.15);
        
        if(histIt->second[i+1] && histIt->second[i+1]->GetEntries() > 0) {
            histIt->second[i+1]->Draw("COLZ");
            
            // Aktualizuj pad aby uzyskać dostęp do primitives
            gPad->Update();
            
            // Znajdź i dostosuj pozycję color bar (palette)
            TPaletteAxis* palette = (TPaletteAxis*)histIt->second[i+1]->GetListOfFunctions()->FindObject("palette");
            if(palette) {
                palette->SetX1NDC(0.83);  // Lewa krawędź color bar - obok prawej osi Y
                palette->SetX2NDC(0.88);  // Prawa krawędź color bar
                palette->SetY1NDC(0.15);  // Dolna krawędź = dolna oś X
                palette->SetY2NDC(0.85);  // Górna krawędź = górna oś X
            }
            
            // Dodaj tekst z informacją o kanale MC - wycentrowany względem górnej osi X
            TPaveText* channelLabel = new TPaveText(0.3, 0.88, 0.7, 0.96, "NDC");
            channelLabel->SetFillColor(kWhite);
            channelLabel->SetBorderSize(1);
            channelLabel->SetTextSize(0.06);
            if(i < static_cast<Int_t>(fChannelNames.size()) && !fChannelNames[i].IsNull()) {
                channelLabel->AddText(fChannelNames[i]);
            } else {
                channelLabel->AddText(Form("Channel %d", i+1));
            }
            channelLabel->Draw();
        } else {
            histIt->second[i+1]->Draw("COLZ");
            
            // Aktualizuj pad aby uzyskać dostęp do primitives
            gPad->Update();
            
            // Znajdź i dostosuj pozycję color bar (palette) nawet dla pustego histogramu
            TPaletteAxis* palette = (TPaletteAxis*)histIt->second[i+1]->GetListOfFunctions()->FindObject("palette");
            if(palette) {
                palette->SetX1NDC(0.83);  // Lewa krawędź color bar - obok prawej osi Y
                palette->SetX2NDC(0.88);  // Prawa krawędź color bar
                palette->SetY1NDC(0.15);  // Dolna krawędź = dolna oś X
                palette->SetY2NDC(0.85);  // Górna krawędź = górna oś X
            }
            
            // Dodaj tekst z informacją o kanale MC nawet dla pustego histogramu - wycentrowany względem górnej osi X
            TPaveText* channelLabel = new TPaveText(0.3, 0.88, 0.7, 0.96, "NDC");
            channelLabel->SetFillColor(kWhite);
            channelLabel->SetBorderSize(1);
            channelLabel->SetTextSize(0.06);
            if(i < static_cast<Int_t>(fChannelNames.size()) && !fChannelNames[i].IsNull()) {
                channelLabel->AddText(fChannelNames[i]);
            } else {
                channelLabel->AddText(Form("Channel %d", i+1));
            }
            channelLabel->Draw();
        }
    }
    
    compareCanvas->Update();
    canvases.push_back(compareCanvas);

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

Double_t HistManager::ScaleChannelByEntries(const TString& setName, Int_t mctruth) {
    // Sprawdź poprawność parametrów
    if(mctruth < 1 || mctruth > fChannNum) {
        std::cerr << "ERROR in ScaleChannelByEntries: Invalid mctruth " << mctruth 
                  << " (valid range: 1-" << fChannNum << ")" << std::endl;
        return -1.0;
    }
    
    // Znajdź zestaw histogramów 1D
    auto histIt = fHists1D.find(setName);
    if(histIt == fHists1D.end()) {
        std::cerr << "ERROR in ScaleChannelByEntries: Histogram set not found: " 
                  << setName.Data() << std::endl;
        return -1.0;
    }
    
    // Sprawdź czy histogram kanału istnieje (indeks mctruth, bo indeks 0 to suma MC)
    if(mctruth >= static_cast<Int_t>(histIt->second.size())) {
        std::cerr << "ERROR in ScaleChannelByEntries: Channel histogram index out of range" << std::endl;
        return -1.0;
    }
    
    TH1* channelHist = histIt->second[mctruth];
    if(!channelHist) {
        std::cerr << "ERROR in ScaleChannelByEntries: Channel histogram is null" << std::endl;
        return -1.0;
    }
    
    // Pobierz liczbę zliczeń i obecną całkę
    Double_t entries = channelHist->GetEntries();
    Double_t integral = channelHist->Integral();
    
    if(entries <= 0) {
        std::cout << "WARNING in ScaleChannelByEntries: Channel " << mctruth 
                  << " in set " << setName.Data() << " has no entries. No scaling applied." << std::endl;
        return 0.0;
    }
    
    if(integral <= 0) {
        std::cout << "WARNING in ScaleChannelByEntries: Channel " << mctruth 
                  << " in set " << setName.Data() << " has zero integral. No scaling applied." << std::endl;
        return 0.0;
    }
    
    // Oblicz czynnik skalujący aby całka = entries
    Double_t scaleFactor = entries / integral;
    
    // Stwórz kopię histogramu przed skalowaniem do aktualizacji sumy
    TH1* oldChannelHist = (TH1*)channelHist->Clone("temp_for_sum_update");
    
    // Przeskaluj histogram kanału
    channelHist->Scale(scaleFactor);
    
    // Zaktualizuj histogram sumy MC (indeks 0)
    TH1* sumHist = histIt->second[0];
    if(sumHist) {
        // Odejmij stary wkład ze sumy
        sumHist->Add(oldChannelHist, -1.0);
        
        // Dodaj nowy przeskalowany wkład do sumy
        sumHist->Add(channelHist, 1.0);
    }
    
    delete oldChannelHist;
    
    std::cout << "INFO: Channel " << mctruth << " in set " << setName.Data() 
              << " normalized: integral " << integral << " -> " << entries 
              << " (scale factor = " << scaleFactor << ")" << std::endl;
    
    return scaleFactor;
}

Double_t HistManager::ScaleChannel2DByEntries(const TString& setName, Int_t mctruth) {
    // Sprawdź poprawność parametrów
    if(mctruth < 1 || mctruth > fChannNum) {
        std::cerr << "ERROR in ScaleChannel2DByEntries: Invalid mctruth " << mctruth 
                  << " (valid range: 1-" << fChannNum << ")" << std::endl;
        return -1.0;
    }
    
    // Znajdź zestaw histogramów 2D
    auto histIt = fHists2D.find(setName);
    if(histIt == fHists2D.end()) {
        std::cerr << "ERROR in ScaleChannel2DByEntries: Histogram set not found: " 
                  << setName.Data() << std::endl;
        return -1.0;
    }
    
    // Sprawdź czy histogram kanału istnieje (indeks mctruth, bo indeks 0 to suma MC)
    if(mctruth >= static_cast<Int_t>(histIt->second.size())) {
        std::cerr << "ERROR in ScaleChannel2DByEntries: Channel histogram index out of range" << std::endl;
        return -1.0;
    }
    
    TH2* channelHist = histIt->second[mctruth];
    if(!channelHist) {
        std::cerr << "ERROR in ScaleChannel2DByEntries: Channel histogram is null" << std::endl;
        return -1.0;
    }
    
    // Pobierz liczbę zliczeń i obecną całkę
    Double_t entries = channelHist->GetEntries();
    Double_t integral = channelHist->Integral();
    
    if(entries <= 0) {
        std::cout << "WARNING in ScaleChannel2DByEntries: Channel " << mctruth 
                  << " in set " << setName.Data() << " has no entries. No scaling applied." << std::endl;
        return 0.0;
    }
    
    if(integral <= 0) {
        std::cout << "WARNING in ScaleChannel2DByEntries: Channel " << mctruth 
                  << " in set " << setName.Data() << " has zero integral. No scaling applied." << std::endl;
        return 0.0;
    }
    
    // Oblicz czynnik skalujący aby całka = entries
    Double_t scaleFactor = entries / integral;
    
    // Stwórz kopię histogramu przed skalowaniem do aktualizacji sumy
    TH2* oldChannelHist = (TH2*)channelHist->Clone("temp_for_sum_update_2d");
    
    // Przeskaluj histogram kanału
    channelHist->Scale(scaleFactor);
    
    // Zaktualizuj histogram sumy MC (indeks 0)
    TH2* sumHist = histIt->second[0];
    if(sumHist) {
        // Odejmij stary wkład ze sumy
        sumHist->Add(oldChannelHist, -1.0);
        
        // Dodaj nowy przeskalowany wkład do sumy
        sumHist->Add(channelHist, 1.0);
    }
    
    delete oldChannelHist;
    
    std::cout << "INFO: 2D Channel " << mctruth << " in set " << setName.Data() 
              << " normalized: integral " << integral << " -> " << entries 
              << " (scale factor = " << scaleFactor << ")" << std::endl;
    
    return scaleFactor;
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

    return logy ? maxY * 10 : maxY * 1.4;
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
        
        // Pobierz parametry binowania i tytuły osi dla danego indeksu
        Int_t bins;
        Double_t xmin, xmax;
        TString xtitle, ytitle;
        config.GetBinning(index, bins, xmin, xmax);
        config.GetAxisTitles(index, xtitle, ytitle);
        
        // Stwórz histogramy dla tego indeksu (suma MC + wszystkie kanały MC)
        std::vector<TH1D*> indexHists(fChannNum + 1); // +1 dla sumy MC (index 0)
        
        // Histogram sumy MC (index 0)
        TString sumHistName = Form("%s_sum_%d", config.baseName.Data(), index);
        indexHists[0] = new TH1D(sumHistName, varTitle, bins, xmin, xmax);
        indexHists[0]->GetXaxis()->SetTitle(xtitle);
        indexHists[0]->GetYaxis()->SetTitle(ytitle);
        ConfigureHistogram(indexHists[0], fSumColor, config.commonConfig.showStats);
        
        // Histogramy MC (indeksy 1 do fChannNum)
        for(Int_t ch = 0; ch < fChannNum; ++ch) {
            TString mcHistName = Form("%s_mc%d_%d", config.baseName.Data(), ch+1, index);
            indexHists[ch + 1] = new TH1D(mcHistName, varTitle, bins, xmin, xmax);
            indexHists[ch + 1]->GetXaxis()->SetTitle(xtitle);
            indexHists[ch + 1]->GetYaxis()->SetTitle(ytitle);
            ConfigureHistogram(indexHists[ch + 1], fChannColors[ch], config.commonConfig.showStats);
        }
        
        arraySet[index] = indexHists;
        
        // Histogram danych
        TString dataHistName = Form("%s_data_%d", config.baseName.Data(), index);
        dataSet[index] = new TH1D(dataHistName, varTitle, bins, xmin, xmax);
        dataSet[index]->GetXaxis()->SetTitle(xtitle);
        dataSet[index]->GetYaxis()->SetTitle(ytitle);
        dataSet[index]->SetMarkerStyle(fDataStyle);
        dataSet[index]->SetMarkerSize(fDataSize);
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
        
        // Wykonaj simple scaling jeśli włączone i dane są dostępne
        if(drawData && fNormalizationType == NormalizationType::SIMPLE_SCALE) {
            try {
                Double_t scaleFactor = PerformArraySimpleScaling(baseName, index);
                
                if(scaleFactor > 0) {
                    // Po skalowaniu, ponownie oblicz zakres Y
                    Double_t newMaxY = 0;
                    for(auto* hist : it->second[index]) {
                        if(hist && hist->GetEntries() > 0) {
                            newMaxY = TMath::Max(newMaxY, hist->GetMaximum());
                        }
                    }
                    
                    // Sprawdź maksimum danych ponownie
                    auto dataIt = fArrayData1D.find(baseName);
                    if(dataIt != fArrayData1D.end() && dataIt->second[index] && dataIt->second[index]->GetEntries() > 0) {
                        newMaxY = TMath::Max(newMaxY, dataIt->second[index]->GetMaximum());
                    }
                    
                    maxY = logY ? newMaxY * 5 : newMaxY * 1.1;
                    
                    std::cout << "Simple scaling applied to array " << baseName.Data() 
                              << "[" << index << "] with factor " << scaleFactor << std::endl;
                }
            } catch(const std::exception& e) {
                std::cerr << "Error during Array Simple Scaling: " << e.what() << std::endl;
            }
        }
        
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

        // Check for Triple Gaussian fit results and draw them if available
        // Check MC histograms (mctruth 1 to fChannNum)
        for(Int_t ch = 1; ch <= fChannNum; ++ch) {
            TString key = GenerateTripleGaussKey(baseName, ch, index);
            auto fitIt = f3GaussFitResults.find(key);
            if(fitIt != f3GaussFitResults.end() && fitIt->second.converged) {
                DrawTripleGaussFit(it->second[index][ch], fitIt->second, false);
            }
        }
        // Check sum MC (mctruth 0)
        TString keySumMC = GenerateTripleGaussKey(baseName, 0, index);
        auto fitItSumMC = f3GaussFitResults.find(keySumMC);
        if(fitItSumMC != f3GaussFitResults.end() && fitItSumMC->second.converged) {
            DrawTripleGaussFit(it->second[index][0], fitItSumMC->second, false);
        }
        // Check data histogram (mctruth -1)
        if(drawData) {
            TString keyData = GenerateTripleGaussKey(baseName, -1, index);
            auto fitItData = f3GaussFitResults.find(keyData);
            if(fitItData != f3GaussFitResults.end() && fitItData->second.converged) {
                auto dataIt = fArrayData1D.find(baseName);
                if(dataIt != fArrayData1D.end()) {
                    DrawTripleGaussFit(dataIt->second[index], fitItData->second, false);
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

Double_t HistManager::ScaleArrayChannelByEntries(const TString& baseName, Int_t index, Int_t mctruth) {
    // Sprawdź poprawność parametrów
    if(mctruth < 1 || mctruth > fChannNum) {
        std::cerr << "ERROR in ScaleArrayChannelByEntries: Invalid mctruth " << mctruth 
                  << " (valid range: 1-" << fChannNum << ")" << std::endl;
        return -1.0;
    }
    
    // Znajdź zestaw array histogramów
    auto histIt = fArrayHists1D.find(baseName);
    if(histIt == fArrayHists1D.end()) {
        std::cerr << "ERROR in ScaleArrayChannelByEntries: Array histogram set not found: " 
                  << baseName.Data() << std::endl;
        return -1.0;
    }
    
    // Sprawdź poprawność indeksu
    if(index < 0 || index >= static_cast<Int_t>(histIt->second.size())) {
        std::cerr << "ERROR in ScaleArrayChannelByEntries: Invalid index " << index 
                  << " (valid range: 0-" << (histIt->second.size()-1) << ")" << std::endl;
        return -1.0;
    }
    
    // Sprawdź czy histogram kanału istnieje
    if(mctruth >= static_cast<Int_t>(histIt->second[index].size())) {
        std::cerr << "ERROR in ScaleArrayChannelByEntries: Channel histogram index out of range" << std::endl;
        return -1.0;
    }
    
    TH1D* channelHist = histIt->second[index][mctruth];
    if(!channelHist) {
        std::cerr << "ERROR in ScaleArrayChannelByEntries: Channel histogram is null" << std::endl;
        return -1.0;
    }
    
    // Pobierz liczbę zliczeń i obecną całkę
    Double_t entries = channelHist->GetEntries();
    Double_t integral = channelHist->Integral();
    
    if(entries <= 0) {
        std::cout << "WARNING in ScaleArrayChannelByEntries: Channel " << mctruth 
                  << " at index " << index << " in array " << baseName.Data() 
                  << " has no entries. No scaling applied." << std::endl;
        return 0.0;
    }
    
    if(integral <= 0) {
        std::cout << "WARNING in ScaleArrayChannelByEntries: Channel " << mctruth 
                  << " at index " << index << " in array " << baseName.Data() 
                  << " has zero integral. No scaling applied." << std::endl;
        return 0.0;
    }
    
    // Oblicz czynnik skalujący aby całka = entries
    Double_t scaleFactor = entries / integral;
    
    // Stwórz kopię histogramu przed skalowaniem do aktualizacji sumy
    TH1D* oldChannelHist = (TH1D*)channelHist->Clone("temp_for_array_sum_update");
    
    // Przeskaluj histogram kanału
    channelHist->Scale(scaleFactor);
    
    // Zaktualizuj histogram sumy MC (indeks 0)
    TH1D* sumHist = histIt->second[index][0];
    if(sumHist) {
        // Odejmij stary wkład ze sumy
        sumHist->Add(oldChannelHist, -1.0);
        
        // Dodaj nowy przeskalowany wkład do sumy
        sumHist->Add(channelHist, 1.0);
    }
    
    delete oldChannelHist;
    
    std::cout << "INFO: Array channel " << mctruth << " at index " << index 
              << " in array " << baseName.Data() 
              << " normalized: integral " << integral << " -> " << entries 
              << " (scale factor = " << scaleFactor << ")" << std::endl;
    
    return scaleFactor;
}

Double_t HistManager::PerformArraySimpleScaling(const TString& baseName, Int_t index) {
    // Znajdź zestaw array histogramów MC
    auto histIt = fArrayHists1D.find(baseName);
    if(histIt == fArrayHists1D.end()) {
        std::cerr << "ERROR: Array histogram set not found: " << baseName.Data() << std::endl;
        return -1.0;
    }
    
    // Sprawdź poprawność indeksu
    if(index < 0 || index >= static_cast<Int_t>(histIt->second.size())) {
        std::cerr << "ERROR: Invalid index " << index << " for array " << baseName.Data() 
                  << " (size: " << histIt->second.size() << ")" << std::endl;
        return -1.0;
    }
    
    // Znajdź histogram danych dla tego indeksu
    auto dataIt = fArrayData1D.find(baseName);
    if(dataIt == fArrayData1D.end()) {
        std::cerr << "ERROR: No data histogram array found for set: " << baseName.Data() << std::endl;
        return -1.0;
    }
    
    if(index >= static_cast<Int_t>(dataIt->second.size()) || !dataIt->second[index]) {
        std::cerr << "ERROR: No data histogram found at index " << index << " for set: " << baseName.Data() << std::endl;
        return -1.0;
    }
    
    // Pobierz liczbę zdarzeń w danych i w sumie MC
    Double_t dataEntries = dataIt->second[index]->Integral();
    Double_t mcSumEntries = histIt->second[index][0]->Integral(); // Indeks 0 to suma MC
    
    if(mcSumEntries <= 0) {
        std::cerr << "WARNING: MC sum histogram is empty for array " << baseName.Data() 
                  << " at index " << index << std::endl;
        return -1.0;
    }
    
    // Oblicz czynnik skalujący
    Double_t scaleFactor = dataEntries / mcSumEntries;
    
    // Przygotuj wektor histogramów MC (bez sumy - indeks 0)
    std::vector<TH1*> mcHists;
    for(Int_t ch = 1; ch <= fChannNum; ++ch) {
        if(ch < static_cast<Int_t>(histIt->second[index].size()) && histIt->second[index][ch]) {
            mcHists.push_back(histIt->second[index][ch]);
        }
    }
    
    // Zastosuj skalowanie używając istniejącej funkcji
    ApplySimpleScaling(mcHists, histIt->second[index][0], scaleFactor);
    
    std::cout << "Array simple scaling applied to " << baseName.Data() << "[" << index << "]"
              << ": scale factor = " << scaleFactor 
              << " (Data: " << dataEntries << ", MC: " << mcSumEntries << ")" << std::endl;
    
    return scaleFactor;
}

// ==================== TRIPLE GAUSSIAN FITTING IMPLEMENTATION ====================

Bool_t HistManager::PrepareDefaultTripleGaussParams(TH1* hist, FitParams3Gauss& params) {
    if(!hist || hist->GetEntries() == 0) {
        std::cerr << "ERROR: Cannot prepare parameters for null or empty histogram" << std::endl;
        return false;
    }
    
    // Reinicjalizuj parametry
    params.InitializeDefaults();
    
    // Pobierz podstawowe właściwości histogramu
    Double_t xmin = hist->GetXaxis()->GetXmin();
    Double_t xmax = hist->GetXaxis()->GetXmax();
    Double_t range = xmax - xmin;
    
    Double_t mean = hist->GetMean();
    Double_t rms = hist->GetRMS();
    Double_t integral = hist->Integral();
    Double_t maxValue = hist->GetMaximum();
    
    // Ustawienia początkowe dla trzech składowych Gaussa
    // Składowa 1: główny pik (najwyższa)
    params.initialParams[0] = integral * 0.6;  // A1 - 60% całkowitej amplitudy
    params.initialParams[1] = mean;             // μ1 - średnia histogramu
    params.initialParams[2] = rms * 0.8;       // σ1 - nieco mniejsza od RMS
    
    // Składowa 2: lewa strona
    params.initialParams[3] = integral * 0.25; // A2 - 25% amplitudy
    params.initialParams[4] = mean - rms * 0.7; // μ2 - przesunięta w lewo
    params.initialParams[5] = rms * 1.2;       // σ2 - szersza
    
    // Składowa 3: prawa strona
    params.initialParams[6] = integral * 0.15; // A3 - 15% amplitudy
    params.initialParams[7] = mean + rms * 0.7; // μ3 - przesunięta w prawo
    params.initialParams[8] = rms * 1.5;       // σ3 - najszersza
    
    // Ustawienia granic parametrów
    for(Int_t i = 0; i < 3; ++i) {
        Int_t idx = i * 3;
        // Amplitudy - od 0 do 2x całkowita wartość
        params.lowerBounds[idx] = 0.0;
        params.upperBounds[idx] = integral * 2.0;
        
        // Średnie - w zakresie histogramu +/- 20%
        params.lowerBounds[idx + 1] = xmin - range * 0.2;
        params.upperBounds[idx + 1] = xmax + range * 0.2;
        
        // Szerokości - od bardzo małej do całego zakresu
        params.lowerBounds[idx + 2] = range * 0.01;
        params.upperBounds[idx + 2] = range;
    }
    
    // Ustaw domyślny zakres fitu na cały histogram
    params.SetFitRange(xmin, xmax);
    
    std::cout << "INFO: Default Triple Gaussian parameters prepared:" << std::endl;
    std::cout << "  Histogram: mean=" << mean << ", RMS=" << rms << ", integral=" << integral << std::endl;
    std::cout << "  Component 1: A=" << params.initialParams[0] << ", μ=" << params.initialParams[1] << ", σ=" << params.initialParams[2] << std::endl;
    std::cout << "  Component 2: A=" << params.initialParams[3] << ", μ=" << params.initialParams[4] << ", σ=" << params.initialParams[5] << std::endl;
    std::cout << "  Component 3: A=" << params.initialParams[6] << ", μ=" << params.initialParams[7] << ", σ=" << params.initialParams[8] << std::endl;
    
    return true;
}

Bool_t HistManager::FitTripleGauss1D(const TString& setName, Int_t mctruth, FitParams3Gauss& params) {
    TH1* hist = GetHistogramForTripleGaussFit(setName, mctruth, -1);
    if(!hist) {
        std::cerr << "ERROR: Cannot find histogram for Triple Gaussian fit: " << setName.Data() 
                  << ", mctruth=" << mctruth << std::endl;
        return false;
    }
    
    // Wykonaj fit
    if(!DoTripleGaussFit(hist, params)) {
        return false;
    }
    
    // Przechowaj wyniki
    TString key = GenerateTripleGaussKey(setName, mctruth, -1);
    f3GaussFitResults[key] = params;
    
    std::cout << "INFO: Triple Gaussian fit completed for " << setName.Data() << ", mctruth=" << mctruth << std::endl;
    return true;
}

Bool_t HistManager::FitTripleGaussArray1D(const TString& baseName, Int_t index, Int_t mctruth, FitParams3Gauss& params) {
    TH1* hist = GetHistogramForTripleGaussFit(baseName, mctruth, index);
    if(!hist) {
        std::cerr << "ERROR: Cannot find array histogram for Triple Gaussian fit: " << baseName.Data() 
                  << ", index=" << index << ", mctruth=" << mctruth << std::endl;
        return false;
    }
    
    // Wykonaj fit
    if(!DoTripleGaussFit(hist, params)) {
        return false;
    }
    
    // Przechowaj wyniki
    TString key = GenerateTripleGaussKey(baseName, mctruth, index);
    f3GaussArrayResults[key] = params;
    
    std::cout << "INFO: Triple Gaussian fit completed for array " << baseName.Data() 
              << ", index=" << index << ", mctruth=" << mctruth << std::endl;
    return true;
}

TH1* HistManager::GetHistogramForTripleGaussFit(const TString& setName, Int_t mctruth, Int_t arrayIndex) {
    if(arrayIndex >= 0) {
        // Array histogram
        auto arrayIt = fArrayHists1D.find(setName);
        if(arrayIt == fArrayHists1D.end()) {
            return nullptr;
        }
        
        if(arrayIndex >= static_cast<Int_t>(arrayIt->second.size())) {
            return nullptr;
        }
        
        if(mctruth == -1) {
            // Data histogram
            auto dataIt = fArrayData1D.find(setName);
            if(dataIt == fArrayData1D.end() || arrayIndex >= static_cast<Int_t>(dataIt->second.size())) {
                return nullptr;
            }
            return dataIt->second[arrayIndex];
        } else if(mctruth >= 0 && mctruth <= fChannNum) {
            // MC histogram (0 = suma, 1-N = kanały)
            if(mctruth >= static_cast<Int_t>(arrayIt->second[arrayIndex].size())) {
                return nullptr;
            }
            return arrayIt->second[arrayIndex][mctruth];
        }
    } else {
        // Regular set histogram
        if(mctruth == -1) {
            // Data histogram
            auto dataIt = fData1D.find(setName);
            if(dataIt == fData1D.end()) {
                return nullptr;
            }
            return dataIt->second;
        } else if(mctruth >= 0 && mctruth <= fChannNum) {
            // MC histogram (0 = suma, 1-N = kanały)
            auto histIt = fHists1D.find(setName);
            if(histIt == fHists1D.end() || mctruth >= static_cast<Int_t>(histIt->second.size())) {
                return nullptr;
            }
            return histIt->second[mctruth];
        }
    }
    
    return nullptr;
}

Bool_t HistManager::DoTripleGaussFit(TH1* hist, FitParams3Gauss& params) {
    if(!hist || hist->GetEntries() == 0) {
        std::cerr << "ERROR: Cannot fit null or empty histogram" << std::endl;
        return false;
    }
    
    // Przygotuj funkcję fitu
    Double_t fitRangeMin = params.fitRangeMin;
    Double_t fitRangeMax = params.fitRangeMax;
    
    if(fitRangeMin <= -999 || fitRangeMax <= -999) {
        fitRangeMin = hist->GetXaxis()->GetXmin();
        fitRangeMax = hist->GetXaxis()->GetXmax();
        params.SetFitRange(fitRangeMin, fitRangeMax);
    }
    
    // Stwórz funkcję TF1 używającą triple_gaus
    TString funcName = Form("triple_gaus_%p_%ld", (void*)hist, (long)time(nullptr));
    TF1* fitFunc = new TF1(funcName, triple_gaus, fitRangeMin, fitRangeMax, 9);
    
    // Ustaw nazwy i początkowe wartości parametrów
    for(Int_t i = 0; i < 9; ++i) {
        fitFunc->SetParName(i, params.paramNames[i]);
        fitFunc->SetParameter(i, params.initialParams[i]);
        fitFunc->SetParLimits(i, params.lowerBounds[i], params.upperBounds[i]);
    }
    
    // Wykonaj fit
    Int_t fitStatus = hist->Fit(fitFunc, "RBQ"); // R=range, B=bounds, Q=quiet (bez S)
    
    // Sprawdź wyniki na podstawie statusu fitu i jakości funkcji
    params.converged = (fitStatus == 0 && fitFunc->GetChisquare() > 0);
    params.status = fitStatus;
    
    if(params.converged) {
        // Pobierz wyniki fitu
        for(Int_t i = 0; i < 9; ++i) {
            params.fittedParams[i] = fitFunc->GetParameter(i);
            params.paramErrors[i] = fitFunc->GetParError(i);
        }
        
        params.chi2 = fitFunc->GetChisquare();
        params.ndf = fitFunc->GetNDF();
        
        // Oblicz wartości kombinowane
        CalculateCombinedValues(params);
        
        // Przechowaj funkcję
        params.fitFunction = (TF1*)fitFunc->Clone();
        
        std::cout << "INFO: Triple Gaussian fit successful. Chi2/NDF = " 
                  << params.chi2 << "/" << params.ndf << " = " << (params.ndf > 0 ? params.chi2/params.ndf : 0) 
                  << std::endl;
    } else {
        std::cerr << "WARNING: Triple Gaussian fit failed or did not converge (status = " 
                  << fitStatus << ")" << std::endl;
    }
    
    delete fitFunc; // Usuwamy oryginalną funkcję, kopia jest przechowana w params
    return params.converged;
}

void HistManager::CalculateCombinedValues(FitParams3Gauss& params) {
    if(!params.converged) {
        return;
    }
    
    // Użyj funkcji z triple_gaus.h do obliczenia wartości kombinowanych
    params.combinedMean = comb_mean(params.fittedParams.data(), params.paramErrors.data());
    params.combinedMeanErr = comb_mean_err(params.fittedParams.data(), params.paramErrors.data());
    params.combinedStdDev = comb_std_dev(params.fittedParams.data(), params.paramErrors.data());
    params.combinedStdDevErr = comb_std_dev_err(params.fittedParams.data(), params.paramErrors.data());
}

TString HistManager::GenerateTripleGaussKey(const TString& setName, Int_t mctruth, Int_t arrayIndex) {
    if(arrayIndex >= 0) {
        return Form("%s_%d_%d", setName.Data(), arrayIndex, mctruth);
    } else {
        return Form("%s_%d", setName.Data(), mctruth);
    }
}

void HistManager::DrawTripleGaussFit(TH1* hist, const FitParams3Gauss& params, Bool_t drawComponents) {
    if(!hist || !params.converged || !params.fitFunction) {
        std::cerr << "ERROR: Cannot draw fit - histogram or fit results invalid" << std::endl;
        return;
    }
    
    // Narysuj główną funkcję fitu (przerywaną czarną linią)
    params.fitFunction->SetLineColor(kBlack);
    params.fitFunction->SetLineStyle(2); // przerywana linia
    params.fitFunction->SetLineWidth(2);
    params.fitFunction->Draw("SAME");
    
    if(drawComponents) {
        DrawTripleGaussComponents(hist, params);
    }
}

void HistManager::DrawTripleGaussComponents(TH1* hist, const FitParams3Gauss& params) {
    if(!hist || !params.converged) {
        std::cerr << "ERROR: Cannot draw components - histogram or fit results invalid" << std::endl;
        return;
    }
    
    // Narysuj poszczególne komponenty Gaussa
    Double_t xmin = params.fitRangeMin;
    Double_t xmax = params.fitRangeMax;
    
    // Komponent 1 (główny) - czerwona linia punktowa
    TF1* comp1 = new TF1(Form("comp1_%p", (void*)hist), single_gaus, xmin, xmax, 3);
    comp1->SetParameters(params.fittedParams[0], params.fittedParams[1], params.fittedParams[2]);
    comp1->SetLineColor(kRed);
    comp1->SetLineStyle(3); // punktowa
    comp1->SetLineWidth(1);
    comp1->Draw("SAME");
    
    // Komponent 2 - niebieska linia punktowa
    TF1* comp2 = new TF1(Form("comp2_%p", (void*)hist), single_gaus, xmin, xmax, 3);
    comp2->SetParameters(params.fittedParams[3], params.fittedParams[4], params.fittedParams[5]);
    comp2->SetLineColor(kBlue);
    comp2->SetLineStyle(3);
    comp2->SetLineWidth(1);
    comp2->Draw("SAME");
    
    // Komponent 3 - zielona linia punktowa
    TF1* comp3 = new TF1(Form("comp3_%p", (void*)hist), single_gaus, xmin, xmax, 3);
    comp3->SetParameters(params.fittedParams[6], params.fittedParams[7], params.fittedParams[8]);
    comp3->SetLineColor(kGreen);
    comp3->SetLineStyle(3);
    comp3->SetLineWidth(1);
    comp3->Draw("SAME");
}

void HistManager::DisplayTripleGaussResults(const FitParams3Gauss& params, const TH1* hist, Bool_t addToCanvas) {
    if(!params.converged) {
        std::cout << "WARNING: Triple Gaussian fit did not converge" << std::endl;
        return;
    }
    
    // Wypisz wyniki na terminalu
    std::cout << "================== TRIPLE GAUSSIAN FIT RESULTS ==================" << std::endl;
    std::cout << "Histogram: " << (hist ? hist->GetName() : "Unknown") << std::endl;
    std::cout << "Fit Status: " << params.status << " (0 = success)" << std::endl;
    std::cout << "Chi2/NDF: " << params.chi2 << "/" << params.ndf;
    if(params.ndf > 0) {
        std::cout << " = " << params.chi2/params.ndf;
    }
    std::cout << std::endl;
    std::cout << "Fit Range: [" << params.fitRangeMin << ", " << params.fitRangeMax << "]" << std::endl;
    std::cout << std::endl;
    
    // Wypisz parametry poszczególnych komponentów
    for(Int_t i = 0; i < 3; ++i) {
        Int_t idx = i * 3;
        std::cout << "Component " << (i+1) << ":" << std::endl;
        std::cout << "  Amplitude: " << params.fittedParams[idx] << " ± " << params.paramErrors[idx] << std::endl;
        std::cout << "  Mean:      " << params.fittedParams[idx+1] << " ± " << params.paramErrors[idx+1] << std::endl;
        std::cout << "  Sigma:     " << params.fittedParams[idx+2] << " ± " << params.paramErrors[idx+2] << std::endl;
        std::cout << std::endl;
    }
    
    // Wypisz wartości kombinowane
    std::cout << "Combined Values (using triple_gaus functions):" << std::endl;
    std::cout << "  Combined Mean:   " << params.combinedMean << " ± " << params.combinedMeanErr << std::endl;
    std::cout << "  Combined Sigma:  " << params.combinedStdDev << " ± " << params.combinedStdDevErr << std::endl;
    std::cout << "=========================================================" << std::endl;
    
    if(addToCanvas && gPad) {
        // Dodaj tekst z wynikami na canvas
        Double_t textX = 0.02;
        Double_t textY = 0.98;
        Double_t lineHeight = 0.05;
        
        // Główne wyniki
        TLatex* latex = new TLatex();
        latex->SetNDC();
        latex->SetTextSize(0.025);
        latex->SetTextFont(42);
        
        latex->DrawLatex(textX, textY, Form("Triple Gaussian Fit"));
        textY -= lineHeight;
        
        latex->DrawLatex(textX, textY, Form("#chi^{2}/NDF = %.1f/%d = %.2f", 
                                           params.chi2, params.ndf, 
                                           params.ndf > 0 ? params.chi2/params.ndf : 0));
        textY -= lineHeight;
        
        latex->DrawLatex(textX, textY, Form("Combined Mean: %.3f #pm %.3f", 
                                           params.combinedMean, params.combinedMeanErr));
        textY -= lineHeight;
        
        latex->DrawLatex(textX, textY, Form("Combined #sigma: %.3f #pm %.3f", 
                                           params.combinedStdDev, params.combinedStdDevErr));
        textY -= lineHeight;
        
        // Poszczególne komponenty (w skrócie)
        latex->SetTextSize(0.02);
        for(Int_t i = 0; i < 3; ++i) {
            Int_t idx = i * 3;
            latex->DrawLatex(textX, textY, 
                           Form("G%d: A=%.1f, #mu=%.2f, #sigma=%.2f", 
                                i+1, params.fittedParams[idx], 
                                params.fittedParams[idx+1], 
                                params.fittedParams[idx+2]));
            textY -= lineHeight * 0.8;
        }
        
        gPad->Update();
    }
}

TCanvas* HistManager::DrawHistogram1DWithFit(const TString& setName, Int_t mctruth, 
                                            const TString& drawOpt,
                                            Bool_t performFit,
                                            FitParams3Gauss* fitParams,
                                            Bool_t drawData) {
    // Znajdź histogram
    TH1* hist = GetHistogramForTripleGaussFit(setName, mctruth, -1);
    if(!hist) {
        std::cerr << "ERROR: Cannot find histogram " << setName.Data() << ", mctruth=" << mctruth << std::endl;
        return nullptr;
    }
    
    // Stwórz canvas
    TString canvasName = Form("c_%s_mc%d_fit", setName.Data(), mctruth);
    TString canvasTitle = Form("%s - MC %d with Triple Gaussian Fit", setName.Data(), mctruth);
    if(mctruth == -1) {
        canvasName = Form("c_%s_data_fit", setName.Data());
        canvasTitle = Form("%s - Data with Triple Gaussian Fit", setName.Data());
    } else if(mctruth == 0) {
        canvasName = Form("c_%s_sum_fit", setName.Data());
        canvasTitle = Form("%s - MC Sum with Triple Gaussian Fit", setName.Data());
    }
    
    TCanvas* canvas = new TCanvas(canvasName, canvasTitle, 800, 600);
    canvas->SetLeftMargin(0.15);
    canvas->SetBottomMargin(0.15);
    canvas->SetRightMargin(0.35); // Więcej miejsca na tekst wyników
    canvas->SetTopMargin(0.10);
    
    // Narysuj histogram
    hist->Draw(drawOpt);
    
    // Narysuj dane jeśli wymagane
    if(drawData && mctruth != -1) {
        auto dataIt = fData1D.find(setName);
        if(dataIt != fData1D.end()) {
            dataIt->second->SetMarkerStyle(fDataStyle);
            dataIt->second->SetMarkerSize(fDataSize);
            dataIt->second->SetMarkerColor(fDataColor);
            dataIt->second->SetLineColor(fDataColor);
            dataIt->second->Draw("PE SAME");
        }
    }
    
    // Wykonaj fit jeśli wymagany
    if(performFit) {
        FitParams3Gauss localParams;
        FitParams3Gauss* paramsToUse = &localParams;
        
        if(fitParams) {
            paramsToUse = fitParams;
        } else {
            // Przygotuj domyślne parametry
            if(!PrepareDefaultTripleGaussParams(hist, localParams)) {
                std::cerr << "ERROR: Cannot prepare default fit parameters" << std::endl;
                return canvas;
            }
        }
        
        // Wykonaj fit
        Bool_t fitSuccess = false;
        if(FitTripleGauss1D(setName, mctruth, *paramsToUse)) {
            fitSuccess = true;
            // Narysuj wyniki fitu
            DrawTripleGaussFit(hist, *paramsToUse, false);
            
            // Wyświetl wyniki
            DisplayTripleGaussResults(*paramsToUse, hist, true);
        }
        
        if(!fitSuccess) {
            std::cerr << "WARNING: Triple Gaussian fit failed" << std::endl;
        }
    }
    
    canvas->Update();
    
    // Zapisz canvas w mapie dla możliwości eksportu
    if(fCanvases.find(setName) == fCanvases.end()) {
        fCanvases[setName] = std::vector<TCanvas*>();
    }
    fCanvases[setName].push_back(canvas);
    
    return canvas;
}

TCanvas* HistManager::DrawArrayHistogramWithFit(const TString& baseName, Int_t index, Int_t mctruth,
                                               const TString& drawOpt, 
                                               Bool_t performFit,
                                               FitParams3Gauss* fitParams) {
    // Znajdź histogram
    TH1* hist = GetHistogramForTripleGaussFit(baseName, mctruth, index);
    if(!hist) {
        std::cerr << "ERROR: Cannot find array histogram " << baseName.Data() 
                  << ", index=" << index << ", mctruth=" << mctruth << std::endl;
        return nullptr;
    }
    
    // Stwórz canvas
    TString canvasName = Form("c_%s_%d_mc%d_fit", baseName.Data(), index, mctruth);
    TString canvasTitle = Form("%s[%d] - MC %d with Triple Gaussian Fit", baseName.Data(), index, mctruth);
    if(mctruth == -1) {
        canvasName = Form("c_%s_%d_data_fit", baseName.Data(), index);
        canvasTitle = Form("%s[%d] - Data with Triple Gaussian Fit", baseName.Data(), index);
    } else if(mctruth == 0) {
        canvasName = Form("c_%s_%d_sum_fit", baseName.Data(), index);
        canvasTitle = Form("%s[%d] - MC Sum with Triple Gaussian Fit", baseName.Data(), index);
    }
    
    TCanvas* canvas = new TCanvas(canvasName, canvasTitle, 800, 600);
    canvas->SetLeftMargin(0.15);
    canvas->SetBottomMargin(0.15);
    canvas->SetRightMargin(0.35); // Więcej miejsca na tekst wyników
    canvas->SetTopMargin(0.10);
    
    // Narysuj histogram
    hist->Draw(drawOpt);
    
    // Wykonaj fit jeśli wymagany
    if(performFit) {
        FitParams3Gauss localParams;
        FitParams3Gauss* paramsToUse = &localParams;
        
        if(fitParams) {
            paramsToUse = fitParams;
        } else {
            // Przygotuj domyślne parametry
            if(!PrepareDefaultTripleGaussParams(hist, localParams)) {
                std::cerr << "ERROR: Cannot prepare default fit parameters" << std::endl;
                return canvas;
            }
        }
        
        // Wykonaj fit
        Bool_t fitSuccess = false;
        if(FitTripleGaussArray1D(baseName, index, mctruth, *paramsToUse)) {
            fitSuccess = true;
            // Narysuj wyniki fitu
            DrawTripleGaussFit(hist, *paramsToUse, false);
            
            // Wyświetl wyniki
            DisplayTripleGaussResults(*paramsToUse, hist, true);
        }
        
        if(!fitSuccess) {
            std::cerr << "WARNING: Triple Gaussian fit failed" << std::endl;
        }
    }
    
    canvas->Update();
    
    // Zapisz canvas w mapie dla możliwości eksportu
    if(fCanvases.find(baseName) == fCanvases.end()) {
        fCanvases[baseName] = std::vector<TCanvas*>();
    }
    fCanvases[baseName].push_back(canvas);
    
    return canvas;
}

// ==================== TRIPLE GAUSSIAN UTILITY METHODS ====================

const HistManager::FitParams3Gauss* HistManager::GetTripleGaussResults1D(const TString& setName, Int_t mctruth) {
    TString key = GenerateTripleGaussKey(setName, mctruth, -1);
    auto it = f3GaussFitResults.find(key);
    if(it != f3GaussFitResults.end()) {
        return &(it->second);
    }
    return nullptr;
}

const HistManager::FitParams3Gauss* HistManager::GetTripleGaussResultsArray(const TString& baseName, Int_t index, Int_t mctruth) {
    TString key = GenerateTripleGaussKey(baseName, mctruth, index);
    auto it = f3GaussArrayResults.find(key);
    if(it != f3GaussArrayResults.end()) {
        return &(it->second);
    }
    return nullptr;
}

Bool_t HistManager::HasTripleGaussResults(const TString& setName, Int_t mctruth, Int_t arrayIndex) {
    if(arrayIndex >= 0) {
        return GetTripleGaussResultsArray(setName, arrayIndex, mctruth) != nullptr;
    } else {
        return GetTripleGaussResults1D(setName, mctruth) != nullptr;
    }
}

// Implementacje metod kontroli rozmiaru markerów danych
void HistManager::SetDataMarkerSize(Float_t size) {
    fDataSize = size;
    // Opcjonalnie można od razu zaktualizować istniejące histogramy
    UpdateExistingDataHistograms(); // Odkomentować jeśli potrzeba automatycznej aktualizacji
}

void HistManager::UpdateExistingDataHistograms() {
    // Aktualizuj histogramy 1D
    for(auto& pair : fData1D) {
        TH1* hist = pair.second;
        if(hist) {
            hist->SetMarkerSize(fDataSize);
        }
    }
    
    // Aktualizuj histogramy Array1D
    for(auto& setPair : fArrayData1D) {
        for(auto& hist : setPair.second) {
            if(hist) {
                hist->SetMarkerSize(fDataSize);
            }
        }
    }
    
    // Aktualizuj histogramy 2D (jeśli mają markery)
    for(auto& pair : fData2D) {
        TH2* hist = pair.second;
        if(hist) {
            hist->SetMarkerSize(fDataSize);
        }
    }
}
