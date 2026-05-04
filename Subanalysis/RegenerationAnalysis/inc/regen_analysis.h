#ifndef REGEN_ANALYSIS_H
#define REGEN_ANALYSIS_H

#include <TChain.h>
#include <TTreeReader.h>
#include <TH1D.h>
#include <map>

#include <TVector3.h>

class RegenAnalysis
{
public:
  RegenAnalysis(TChain *chain);
  ~RegenAnalysis();

  void Loop();
  void FitResults();
  void SaveToConfig();

private:
  TTreeReader fReader;
  // Definicje czytników (TTreeReaderValue/Array)

  // Histogramy pogrupowane w mapy dla łatwiejszej iteracji
  std::map<std::string, TH1D *> fHistos;
  std::map<std::string, std::pair<double, double>> histRanges, fitRanges; // Zakresy histogramów (min, max)
  std::map<std::string, Double_t> resolutions;

  std::vector<std::string> regions = {"CylindricalBeamPipe", "SphericalBeamPipe"};
  std::vector<std::string> types = {"Charged", "Neutral"};

  // Bounds
  Double_t
      beamPipeBound,
      driftChamberBound;

  // Ip position
  TVector3 ipVector;
  // Neutral vtx position
  TVector3 neutralVtxVector;
  // Charged vtx position
  TVector3 chargedVtxVector;

  // Config loaded
  Bool_t configLoaded = false;

  void LoadConfig();
  void InitHistograms();
};

#endif