#pragma once
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TString.h>
#include <TLegend.h>
#include <TPad.h>
#include <TFile.h>
#include <TDirectory.h>
#include <vector>
#include <map>

class HistManager {
public:
    struct HistConfig {
        TString name;
        TString title;
        Int_t bins;
        Double_t xmin;
        Double_t xmax;
        TString xtitle;
        TString ytitle = "Counts";
        Bool_t logy = false;
        Bool_t showStats = true;
    };

    struct Hist2DConfig : public HistConfig {
        Int_t binsy;
        Double_t ymin;
        Double_t ymax;
    };

    // Constructor that takes channel number, array of channel colors, optional channel names and data parameters
    HistManager(Int_t channNum, const Color_t* channColors, 
                const std::vector<TString>& channelNames = std::vector<TString>(),
                Int_t dataStyle = kFullCircle, Int_t dataColor = kBlack,
                Int_t sumColor = kOrange);
    ~HistManager();

    // Create a new set of 1D histograms
    void CreateHistSet1D(const TString& setName, const HistConfig& config);
    
    // Create a new set of 2D histograms  
    void CreateHistSet2D(const TString& setName, const Hist2DConfig& config);

    // Fill histogram in a set for given channel (mctruth 1-7)
    void Fill1D(const TString& setName, Int_t mctruth, Double_t value, Double_t weight = 1.0);
    
    // Fill 2D histogram
    void Fill2D(const TString& setName, Int_t mctruth, Double_t x, Double_t y, Double_t weight = 1.0);

    // Fill data points
    void FillData1D(const TString& setName, Double_t value, Double_t weight = 1.0);
    void FillData2D(const TString& setName, Double_t x, Double_t y, Double_t weight = 1.0);

    // Draw all histograms in a 1D set on one canvas with optional data points
    void DrawSet1D(const TString& setName, const TString& drawOpt = "", Bool_t drawData = false);
    
    // Draw each 2D histogram in set on separate canvas with optional data points
    void DrawSet2D(const TString& setName, const TString& drawOpt = "COLZ", Bool_t drawData = false);

    // Save all canvases for a set to image files
    void SaveSet(const TString& setName, const TString& filePattern);
    
    // Save all histograms to ROOT file
    void SaveToRoot(const TString& filename);
    
    // Save specific set to ROOT file
    void SaveSetToRoot(const TString& setName, const TString& filename);
    
    // Export plots to various formats
    enum class ImageFormat {
        PNG,
        SVG,
        PDF
    };
    void ExportSet(const TString& setName, const TString& filePattern, ImageFormat format);

private:
    Int_t fChannNum;
    std::vector<Int_t> fChannColors;
    std::vector<TString> fChannelNames;
    Int_t fDataStyle;
    Int_t fDataColor;
    Int_t fSumColor;  // Kolor dla sumy MC
    
    // Maps to store histograms and canvases
    std::map<TString, std::vector<TH1*>> fHists1D;
    std::map<TString, std::vector<TH2*>> fHists2D;
    std::map<TString, TH1*> fData1D;  // Histogramy dla danych
    std::map<TString, TH2*> fData2D;  // Histogramy 2D dla danych
    std::map<TString, std::vector<TCanvas*>> fCanvases;
    std::map<TString, HistConfig> fConfigs1D;
    std::map<TString, Hist2DConfig> fConfigs2D;

    void ConfigureHistogram(TH1* hist, Int_t color, Bool_t showStats);
    Double_t CalculateYRange(const std::vector<TH1*>& hists, const TString& setName, Bool_t logy);
    void CleanupSet(const TString& setName);
    TLegend* CreateLegend(const std::vector<TH1*>& hists, const TH1* dataHist = nullptr);
};
