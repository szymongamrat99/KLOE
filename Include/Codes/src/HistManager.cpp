#include "HistManager.h"
#include <TMath.h>
#include <TStyle.h>
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
            dataIt->second->SetMarkerStyle(fDataStyle);
            dataIt->second->SetMarkerColor(fDataColor);
            dataIt->second->SetLineColor(fDataColor);
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
