#ifndef REGEN_ANALYSIS_H
#define REGEN_ANALYSIS_H

#include <map>
#include <nlohmann/json.hpp>

#include <TTreeReader.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TVector3.h>

#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

#include <TF1.h>

#include <const.h>
#include <interference.h>

enum class HistogramType
{
  RHO_CHARGED,
  RHO_NEUTRAL,
  R_CHARGED,
  R_NEUTRAL,
  UNKNOWN
};

struct Config
{
  double beamPipeBound;
  double driftChamberBound;
  std::map<HistogramType, std::pair<double, double>> histRanges;
  std::map<HistogramType, std::pair<double, double>> displayRanges;
  std::map<HistogramType, Double_t> resolutions;

  std::map<HistogramType, std::string> titles;
  std::map<HistogramType, std::string> xAxisTitles;
  std::map<HistogramType, std::string> yAxisTitles;

  std::map<HistogramType, bool> logarithmic;

  double u0;
  double v0;
  double su;
  double sv;
  double rho;
  double n_sigma_cut;

  TVector3 DCWallCenter;

  double DCWallBottomBoundCharged;
  double DCWallTopBoundCharged;
  double DCWallBottomBoundNeutral;
  double DCWallTopBoundNeutral;

  double sphericalBPBottomBoundCharged;
  double sphericalBPTopBoundCharged;
  double sphericalBPBottomBoundNeutral;
  double sphericalBPTopBoundNeutral;

  double cylindricalBPBottomBoundCharged;
  double cylindricalBPTopBoundCharged;
  double cylindricalBPBottomBoundNeutral;
  double cylindricalBPTopBoundNeutral;

  bool isSignalWeighted;
  bool isRegenerationWeighted;
  bool onlyKLong;

  double lumiFactor;
  double weightErrorLimit;
};

class RegenerationFractionFit
{
public:
  RegenerationFractionFit(TTreeReader *reader);
  RegenerationFractionFit(std::string weightsFilePath);
  ~RegenerationFractionFit();

  void LoadConfig();
  void SaveHistograms();
  void DrawHistograms() const;
  std::pair<double, double> GetRegenerationWeight(HistogramType histType, double radius) const;
  std::pair<double, double> GetRegenerationWeight(double radius) const;
  double GetContinuousRegenerationWeight(double dt, double radius);
  double GetContinuousRegenerationWeight(double dt, std::map<HistogramType, double> radius);

private:
  using json = nlohmann::json;

  Config fConfig;
  TTreeReader *fReader;

  KLOE::interference fInterf;

  std::map<HistogramType, std::unordered_map<std::string, TH1D *>> fHistos;
  std::map<HistogramType, TH1D *> fRegenerationWeights;
  std::map<std::string, TH1D *> fTimeDiffHist;
  std::map<std::string, TH2D *> fRhovsRChargedHist, fRhovsRNeutralHist;

  std::map<std::string, TTreeReaderValue<Double_t> *> fReaderFloat;
  std::map<std::string, TTreeReaderValue<Int_t> *> fReaderInteger;

  std::map<std::string, TTreeReaderArray<Double_t> *> fReaderFloatArray;

  std::map<HistogramType, std::vector<std::string>> _channInMCSum;

  std::vector<std::string> _channInTimeDiff;

  std::map<HistogramType, TF1 *> fcontWeightDC, fcontWeightCylBP, fcontWeightSphBP;

  void _FillHistograms();
  void _FillMCSumHistogram(bool addRegeneration);
  void _RenormalizeSignalHistogram();
  void _RenormalizeMCToLumi();
  void _ResetHistograms();
  void _SetupJson(json &inputConfig,json &outputConfig, std::string section) const;

  std::pair<TF1*, HistogramType> _GetFitFuncByDt(double dt);

  void _FitContinuousWeightFunction(HistogramType histType);

  // ---

  const std::vector<std::string> floatVars = {"Chi2SignalKinFit", "minv4gam", "Qmiss", "Bx", "By", "Bz", "KaonChTimeCMMC", "KaonNeTimeCMMC", "KaonChTimeCMSignalFit", "KaonNeTimeCMSignalFit", "bestErrorSixGamma"};
  const std::vector<std::string> integerVars = {"mcflag", "mctruth"};
  const std::vector<std::string> floatArrayVars = {"Knerec", "KnerecSix", "KchrecClosest", "KnerecFit", "KchrecFit", "ipFit", "Kchrec", "ip", "trk1Fit", "trk2Fit", "pi01Fit", "pi02Fit"};

  void _SetupReaderVariables();

  // ---

  std::array<HistogramType, 4> _histTypesIterative = {HistogramType::RHO_CHARGED, 
                                                      HistogramType::RHO_NEUTRAL, 
                                                      HistogramType::R_CHARGED, 
                                                      HistogramType::R_NEUTRAL};

  std::string _HistTypeToString(HistogramType type) const;
  HistogramType _StringToHistType(const std::string &str) const;

  // ---

  bool _isAcceptedEvent() const;
  bool _isOutsideEllipse(double u, double v) const;
  bool _applyGlobalCut() const;
  bool _checkRegenerationLimits(double xLeft, double xRight) const;
  HistogramType _checkRegenerationLimitsTypes(double radius) const;
  HistogramType _checkRegenerationLimitsTypesDeltaT(double dt) const;
  TFile *fFileWeights;

  TH1 *threePi0WeightsHist, *semileptonicWeightsHist;

  void _calculateRegenerationWeights();

  double _calculateRadius(HistogramType type) const;
};

#endif