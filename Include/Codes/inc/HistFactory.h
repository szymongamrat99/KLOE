#pragma once

#include <map>

#include "TH1.h"
#include "TH2.h"
#include "TString.h"
#include "TFractionFitter.h"

struct HistConfig
{
    std::string id;
    std::string name;
    std::string title;
    std::string xAxisTitle;
    std::string yAxisTitle;
    Bool_t logScale;
    Bool_t byChannel;
    Bool_t is2D;
    Int_t nBinsX;
    Double_t xMin;
    Double_t xMax;
    Int_t nBinsY;
    Double_t yMin;
    Double_t yMax;
};

struct FitResult
{
  Double_t scaleFactor;
  Double_t scaleFactorError;
};

class HistFactory
{
  public:
      enum class Normalization
      {
        None,
        ToOne,
        MCToData
      };

      std::map<TString, FitResult> fitComponentsToData(
        const TString &id, 
        const std::map<TString, std::pair<Double_t, Double_t>>& customLimits = {}
      );

      HistFactory(const HistConfig& config);
      HistFactory(const std::vector<HistConfig>& configs);

      void fill(const TString& id, Double_t x, Double_t weight = 1.0);
      void fill(const TString& id, const Int_t mctruth, Double_t x, Double_t weight = 1.0);

      void fill(const TString& id, Double_t x,  Double_t y, Double_t weight = 1.0);
      void fill(const TString& id, const Int_t mctruth, Double_t x, Double_t y, Double_t weight = 1.0);

      ~HistFactory() = default;

  private:
      std::map<TString, TH1*> simple1DHistograms;
      std::map<TString, TH2*> simple2DHistograms;

      std::map<TString, std::map<TString, TH1*>> histograms1DByChannel;
      std::map<TString, std::map<TString, TH2*>> histograms2DByChannel;

      void addHistogram(const HistConfig& config);
      void setChannColors(const TString &id, const TString &channel);
      void buildMCSum(const TString &id);
      void buildMCSum2D(const TString &id);
      void normalize(Normalization norm);

      TH1* create1DHistogram(const HistConfig& config, const TString& channel = "");
      TH2* create2DHistogram(const HistConfig& config, const TString& channel = "");

      TH1* get1DHistogram(const TString& id, const TString& channel = "");
      TH2* get2DHistogram(const TString& id, const TString& channel = "");
};