#include <regen_frac_fit.h>

#include <nlohmann/json.hpp>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iomanip> // Wymagane dla std::setprecision i std::fixed
#include <limits>
#include <memory>
#include <sstream> // Wymagane dla std::ostringstream
#include <stdexcept>

#include <boost/filesystem.hpp>
#include <TH1D.h>
#include <TVector3.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TGraphAsymmErrors.h>
#include <TLine.h>
#include <TFitResultPtr.h>
#include <TFitResult.h>

#include <interference.h>
#include <BRCorrectionFactors.h>

using json = nlohmann::json;

namespace
{
double evaluateCorrectionSpline(const TSpline3 &spline, double x)
{
  if (!std::isfinite(x) || x < spline.GetXmin() || x > spline.GetXmax())
    return 1.0;
  const double value = spline.Eval(x);
  return std::isfinite(value) ? value : 1.0;
}
}

RegenerationFractionFit::RegenerationFractionFit(TTreeReader *reader) : fReader(reader), fInterf()
{
  LoadConfig();
  KLOE::setGlobalStyle();

  // TFile *fFileWeightsThreePi0 = TFile::Open("/data/ssd/gamrat/python-kloe-analysis/scripts/results/control_sample_corr_factors/control_sample_corr_factors.root", "READ");

  // TCanvas *c = (TCanvas*)fFileWeightsThreePi0->Get("three_pi0/correction_factors/cCorr");
  // threePi0WeightsHist = (TH1*)c->FindObject("correction_factors")->Clone();

  // TFile *fFileWeightsSemileptonic = TFile::Open("/data/ssd/gamrat/python-kloe-analysis/scripts/results/control_sample_corr_factors/control_sample_corr_factors.root", "READ");

  // c = (TCanvas*)fFileWeightsSemileptonic->Get("semileptonic/correction_factors/cCorr");
  // semileptonicWeightsHist = (TH1*)c->FindObject("correction_factors")->Clone();

  // TFitResultPtr threePi0WeightResult = threePi0WeightsHist->Fit("pol0", "QS");
  avgThreePi0Weight = 0.8868;//threePi0WeightResult->Parameter(0);

  std::cout << "Average Three Pi0 Weight: " << avgThreePi0Weight << std::endl;

  // TFitResultPtr semileptonicWeightResult = semileptonicWeightsHist->Fit("pol0", "QS");
  avgSemileptonicWeight = 0.9844; //semileptonicWeightResult->Parameter(0);

  std::cout << "Average Semileptonic Weight: " << avgSemileptonicWeight << std::endl;

  const auto resolveCPFitPath = [](const std::string &configuredPath)
  {
    const boost::filesystem::path path(configuredPath);
    return (path.is_absolute() ? path : boost::filesystem::path((std::string)Paths::cpfit_dir) / path).string();
  };
  const bool needsKsPmMean = fConfig.scaleThreePi0ByKsPmMean;
  const bool needsKs00Mean = fConfig.scaleSemileptonicByKs00Mean;
  for (const auto &component : fConfig.signalWeightComponents)
  {
    const bool neededForMean =
        (component.first == "ks_pm" && needsKsPmMean) ||
        (component.first == "ks_00" && needsKs00Mean);
    if (!component.second.enabled && !neededForMean)
      continue;

    const std::string correctionPath = resolveCPFitPath(component.second.correctionFile);
    {
      TDirectory::TContext directoryContext;
      std::unique_ptr<TFile> correctionFile(TFile::Open(correctionPath.c_str(), "READ"));
      if (!correctionFile || correctionFile->IsZombie())
        throw std::runtime_error("Could not open signal correction file for " + component.first + ": " +
                                 correctionPath);
      TH1 *source = dynamic_cast<TH1 *>(correctionFile->Get(component.second.correctionHistogram.c_str()));
      if (!source)
        throw std::runtime_error("Could not find signal correction histogram for " + component.first + ": " +
                                 component.second.correctionHistogram);
      std::unique_ptr<TH1> correction(dynamic_cast<TH1 *>(source->Clone(
          ("regen_" + component.first + "_correction").c_str())));
      if (!correction)
        throw std::runtime_error("Could not clone signal correction histogram for " + component.first);
      correction->SetDirectory(nullptr);
      signalWeightHistograms.emplace(component.first, std::move(correction));
      if (component.second.useCorrectionSpline)
      {
        TSpline3 *sourceSpline =
            dynamic_cast<TSpline3 *>(correctionFile->Get(component.second.correctionSpline.c_str()));
        if (!sourceSpline)
          throw std::runtime_error("Could not find signal correction spline for " + component.first + ": " +
                                   component.second.correctionSpline);
        std::unique_ptr<TSpline3> correctionSpline(dynamic_cast<TSpline3 *>(
            sourceSpline->Clone(("regen_" + component.first + "_correction_spline").c_str())));
        if (!correctionSpline)
          throw std::runtime_error("Could not clone signal correction spline for " + component.first);
        signalWeightCorrectionSplines.emplace(component.first, std::move(correctionSpline));
      }
    }

    if (!component.second.ratioFile.empty())
    {
      const std::string ratioPath = resolveCPFitPath(component.second.ratioFile);
      TDirectory::TContext directoryContext;
      std::unique_ptr<TFile> ratioFile(TFile::Open(ratioPath.c_str(), "READ"));
      if (!ratioFile || ratioFile->IsZombie())
        throw std::runtime_error("Could not open signal ratio file for " + component.first + ": " + ratioPath);
      TSpline3 *source = dynamic_cast<TSpline3 *>(ratioFile->Get(component.second.ratioSpline.c_str()));
      if (!source)
        throw std::runtime_error("Could not find signal ratio spline for " + component.first + ": " +
                                 component.second.ratioSpline);
      std::unique_ptr<TSpline3> ratio(dynamic_cast<TSpline3 *>(source->Clone(
          ("regen_" + component.first + "_ratio").c_str())));
      if (!ratio)
        throw std::runtime_error("Could not clone signal ratio spline for " + component.first);
      signalWeightRatios.emplace(component.first, std::move(ratio));
    }
  }

  for (const auto &histType : _histTypesIterative)
  {
    for (const auto &channel : KLOE::channName)
    {
      std::string channelName_str = (std::string)channel.second;
      std::ostringstream resStream;
      resStream << std::fixed << std::setprecision(2) << fConfig.resolutions[histType];

      std::string histUniqueName = _HistTypeToString(histType) + "_" + channelName_str;
      std::string histTotalTitle = fConfig.titles[histType] + ";" +
                                   fConfig.xAxisTitles[histType] + ";" +
                                   fConfig.yAxisTitles[histType] + " / " +
                                   resStream.str() + " cm";

      int nbins = floor((fConfig.histRanges[histType].second - fConfig.histRanges[histType].first) / fConfig.resolutions[histType]);

      // Create the histogram and store it in the map
      TH1D *h = new TH1D(histUniqueName.c_str(), histTotalTitle.c_str(), nbins, fConfig.histRanges[histType].first, fConfig.histRanges[histType].second);
      // h->SetDirectory(0);
      h->GetXaxis()->SetRangeUser(fConfig.displayRanges[histType].first, fConfig.displayRanges[histType].second);

      fHistos[histType][channelName_str] = h;
    }

    fRegenerationWeights[histType] = (TH1D *)fHistos[histType]["Regeneration"]->Clone(("RegenerationWeights" + _HistTypeToString(histType)).c_str());
  }

  for (const auto &channel : KLOE::channName)
  {
    std::string channelName_str = (std::string)channel.second;
    std::string histUniqueName = "TimeDiff_" + channelName_str;
    std::string histTotalTitle = "Time Difference (t_ch - t_ne) for " + channelName_str + ";t_ch - t_ne (ns);Entries";

    TH1D *h = new TH1D(histUniqueName.c_str(), histTotalTitle.c_str(), 150, -300, 300);
    fTimeDiffHist[channelName_str] = h;

    std::string histUniqueNameRhovsRCharged = "RhovsRCharged_" + channelName_str;
    std::string histTotalTitleRhovsRCharged = "Rho vs R Charged for " + channelName_str + "R vs. #rho charged;R [cm];#rho [cm]";
    fRhovsRChargedHist[channelName_str] = new TH2D(histUniqueNameRhovsRCharged.c_str(), histTotalTitleRhovsRCharged.c_str(), 100, 0, 50, 100, 0, 50);

    std::string histUniqueNameRhovsRNeutral = "RhovsRNeutral_" + channelName_str;
    std::string histTotalTitleRhovsRNeutral = "Rho vs R Neutral for " + channelName_str + "R vs. #rho neutral;R [cm];#rho [cm]";
    fRhovsRNeutralHist[channelName_str] = new TH2D(histUniqueNameRhovsRNeutral.c_str(), histTotalTitleRhovsRNeutral.c_str(), 100, 0, 50, 100, 0, 50);

  }

  _SetupReaderVariables();
}

RegenerationFractionFit::RegenerationFractionFit(std::string weightsFilePath) : fReader(nullptr), fInterf()
{
  LoadConfig();
  KLOE::setGlobalStyle();

  TDirectory::TContext directoryContext;
  fFileWeights = TFile::Open(weightsFilePath.c_str(), "READ");
  if (!fFileWeights || fFileWeights->IsZombie())
  {
    delete fFileWeights;
    fFileWeights = nullptr;
    throw std::runtime_error("Could not open regeneration weight file: " + weightsFilePath);
  }

  for (const auto &histType : _histTypesIterative)
  {
    std::string funcDCName = "FuncDC_" + _HistTypeToString(histType);
    fcontWeightDC[histType] = dynamic_cast<TF1*>(fFileWeights->Get(funcDCName.c_str()));

    std::string funcCylBPName = "FuncCylBP_" + _HistTypeToString(histType);
    fcontWeightCylBP[histType] = dynamic_cast<TF1*>(fFileWeights->Get(funcCylBPName.c_str()));

    std::string funcSphBPName = "FuncSphBP_" + _HistTypeToString(histType);
    fcontWeightSphBP[histType] = dynamic_cast<TF1*>(fFileWeights->Get(funcSphBPName.c_str()));

    if (!fcontWeightDC[histType] || !fcontWeightCylBP[histType] || !fcontWeightSphBP[histType])
    {
      delete fFileWeights;
      fFileWeights = nullptr;
      throw std::runtime_error("Missing regeneration weight function for " + _HistTypeToString(histType) +
                               " in " + weightsFilePath);
    }
  }
  for (const auto &histType : _histTypesIterative)
  {
    fRegenerationWeights[histType] = new TH1D(("RegenerationWeights" + _HistTypeToString(histType)).c_str(), ("Regeneration Weights for " + _HistTypeToString(histType)).c_str(), 100, fConfig.histRanges[histType].first, fConfig.histRanges[histType].second);
  }
}

RegenerationFractionFit::~RegenerationFractionFit()
{
  if (!fReaderFloat.empty())
    for (auto &pair : fReaderFloat)
    {
      delete pair.second;
    }

  if (!fReaderInteger.empty())
    for (auto &pair : fReaderInteger)
    {
      delete pair.second;
    }

  if (!fReaderFloatArray.empty())
    for (auto &pair : fReaderFloatArray)
    {
      delete pair.second;
    }

  // Clean up histograms

  if (!fHistos.empty())
  {
    for (auto &histTypePair : fHistos)
    {
      for (auto &channelPair : histTypePair.second)
      {
        delete channelPair.second;
      }
    }
  }

  if (!fRegenerationWeights.empty())
  {
    for (auto &pair : fRegenerationWeights)
    {
      delete pair.second;
    }
  }

  if (!fTimeDiffHist.empty())
  {
    for (auto &pair : fTimeDiffHist)
    {
      delete pair.second;
    }
  }

  if (!fRhovsRChargedHist.empty())
  {
    for (auto &pair : fRhovsRChargedHist)
    {
      delete pair.second;
    }
  }

  if (!fRhovsRNeutralHist.empty())
  {
    for (auto &pair : fRhovsRNeutralHist)
    {
      delete pair.second;
    }
  }
}

void RegenerationFractionFit::LoadConfig()
{
  // Load fit configuration from JSON file
  std::ifstream configFile(Paths::regen_analysis_dir + "config/fit_config.json");
  if (!configFile.is_open())
  {
    std::cerr << "Unable to open fit configuration file!" << std::endl;
    return;
  }

  json config, histConfig, constantCorrectionsConfig, generalSettingsConfig;
  configFile >> config;

  _SetupJson(config, histConfig, "histogramConfiguration");

  std::map<HistogramType, json> rhoRConfigs;

  for (const auto &histType : _histTypesIterative)
  {
    _SetupJson(histConfig, rhoRConfigs[histType], _HistTypeToString(histType));

    try
    {
      fConfig.histRanges[histType] = {rhoRConfigs.at(histType).at("range")[0],
                                      rhoRConfigs.at(histType).at("range")[1]};
      fConfig.displayRanges[histType] = {rhoRConfigs.at(histType).at("displayRange")[0],
                                         rhoRConfigs.at(histType).at("displayRange")[1]};
      fConfig.resolutions[histType] = rhoRConfigs.at(histType).at("resolution");
      fConfig.titles[histType] = rhoRConfigs.at(histType).at("title");
      fConfig.xAxisTitles[histType] = rhoRConfigs.at(histType).at("xAxisTitle");
      fConfig.yAxisTitles[histType] = rhoRConfigs.at(histType).at("yAxisTitle");
      fConfig.logarithmic[histType] = rhoRConfigs.at(histType).at("logarithmic");
    }
    catch (json::out_of_range &e)
    {
      std::cerr << "Missing parameter: " << e.what() << std::endl;
      return;
    }
  }

  _SetupJson(config, constantCorrectionsConfig, "constantCorrections");
  try
  {
    fConfig.DCWallCenter.SetXYZ(constantCorrectionsConfig.at("DCWallCenter")[0],
                                constantCorrectionsConfig.at("DCWallCenter")[1],
                                constantCorrectionsConfig.at("DCWallCenter")[2]);

    fConfig.DCWallBottomBoundCharged = constantCorrectionsConfig.at("DCWallBounds").at("charged").at("bottom");
    fConfig.DCWallTopBoundCharged = constantCorrectionsConfig.at("DCWallBounds").at("charged").at("top");
    fConfig.DCWallBottomBoundNeutral = constantCorrectionsConfig.at("DCWallBounds").at("neutral").at("bottom");
    fConfig.DCWallTopBoundNeutral = constantCorrectionsConfig.at("DCWallBounds").at("neutral").at("top");

    fConfig.sphericalBPBottomBoundCharged = constantCorrectionsConfig.at("sphericalBeamPipeBound").at("charged").at("bottom");
    fConfig.sphericalBPTopBoundCharged = constantCorrectionsConfig.at("sphericalBeamPipeBound").at("charged").at("top");
    fConfig.sphericalBPBottomBoundNeutral = constantCorrectionsConfig.at("sphericalBeamPipeBound").at("neutral").at("bottom");
    fConfig.sphericalBPTopBoundNeutral = constantCorrectionsConfig.at("sphericalBeamPipeBound").at("neutral").at("top");

    fConfig.cylindricalBPBottomBoundCharged = constantCorrectionsConfig.at("cylindricalBeamPipeBound").at("charged").at("bottom");
    fConfig.cylindricalBPTopBoundCharged = constantCorrectionsConfig.at("cylindricalBeamPipeBound").at("charged").at("top");
    fConfig.cylindricalBPBottomBoundNeutral = constantCorrectionsConfig.at("cylindricalBeamPipeBound").at("neutral").at("bottom");
    fConfig.cylindricalBPTopBoundNeutral = constantCorrectionsConfig.at("cylindricalBeamPipeBound").at("neutral").at("top");
  }
  catch (json::out_of_range &e)
  {
    std::cerr << "Missing constant correction parameter: " << e.what() << std::endl;
    return;
  }

  _SetupJson(config, generalSettingsConfig, "generalSettings");
  try
  {
    fConfig.isSignalWeighted = generalSettingsConfig.at("isSignalWeighted");
    fConfig.isRegenerationWeighted = generalSettingsConfig.at("isRegenerationWeighted");
    fConfig.onlyKLong = generalSettingsConfig.at("onlyKLong");
    fConfig.lumiFactor = generalSettingsConfig.at("lumiFactor");
    fConfig.weightErrorLimit = generalSettingsConfig.at("weightErrorLimit");
    fConfig.signalWeightRelativeErrorLimit = generalSettingsConfig.at("signalWeightRelativeErrorLimit");
    const auto meanWeights = generalSettingsConfig.value("crossChannelMeanWeights", json::object());
    fConfig.scaleThreePi0ByKsPmMean = meanWeights.value("threePi0FromKsPm", false);
    fConfig.scaleSemileptonicByKs00Mean = meanWeights.value("semileptonicFromKs00", false);

    const auto &signalComponents = generalSettingsConfig.at("signalWeightComponents");
    const std::array<std::string, 4> componentNames = {"ks_pm", "kl_pm", "ks_00", "kl_00"};
    for (const auto &name : componentNames)
    {
      const auto &entry = signalComponents.at(name);
      Config::SignalWeightComponent component;
      component.enabled = entry.at("enabled").get<bool>();
      component.correctionFile = entry.at("correctionFile").get<std::string>();
      component.correctionHistogram = entry.at("correctionHistogram").get<std::string>();
      component.correctionSpline = entry.value("correctionSpline", "");
      component.useCorrectionSpline = entry.value("useCorrectionSpline", false);
      component.ratioFile = entry.value("ratioFile", "");
      component.ratioSpline = entry.value("ratioSpline", "");
      fConfig.signalWeightComponents.emplace(name, std::move(component));
    }
    for (auto it = signalComponents.begin(); it != signalComponents.end(); ++it)
      if (std::find(componentNames.begin(), componentNames.end(), it.key()) == componentNames.end())
        throw std::runtime_error("Unknown signal weight component: " + it.key());
  }
  catch (json::out_of_range &e)
  {
    throw std::runtime_error("Missing general settings parameter: " + std::string(e.what()));
  }

  if (!std::isfinite(fConfig.weightErrorLimit) || fConfig.weightErrorLimit < 0.0 ||
      !std::isfinite(fConfig.signalWeightRelativeErrorLimit) || fConfig.signalWeightRelativeErrorLimit < 0.0)
    throw std::runtime_error("Signal and regeneration relative-error limits must be finite and non-negative");
  for (const auto &component : fConfig.signalWeightComponents)
  {
    const bool required = component.second.enabled ||
                          (component.first == "ks_pm" && fConfig.scaleThreePi0ByKsPmMean) ||
                          (component.first == "ks_00" && fConfig.scaleSemileptonicByKs00Mean);
    if (required &&
        (component.second.correctionFile.empty() || component.second.correctionHistogram.empty() ||
         (component.second.useCorrectionSpline && component.second.correctionSpline.empty())))
      throw std::runtime_error("Required signal correction component has incomplete ROOT object configuration: " +
                               component.first);
    if (required && !component.second.ratioFile.empty() && component.second.ratioSpline.empty())
      throw std::runtime_error("Signal ratio spline name is missing for component: " + component.first);
  }

  configFile.close();

  // ---
  // Load cut configuration from JSON file
  std::ifstream cutConfigFile(Paths::regen_analysis_dir + "config/cut_definition.json");
  if (!cutConfigFile.is_open())
  {
    std::cerr << "Unable to open cut configuration file!" << std::endl;
    return;
  }

  json cutConfig, ellipseCutConfig;
  cutConfigFile >> cutConfig;
  _SetupJson(cutConfig, ellipseCutConfig, "ellipseNeuInvMassCut");

  try
  {
    fConfig.u0 = ellipseCutConfig.at("u0");
    fConfig.v0 = ellipseCutConfig.at("v0");
    fConfig.su = ellipseCutConfig.at("su");
    fConfig.sv = ellipseCutConfig.at("sv");
    fConfig.rho = ellipseCutConfig.at("rho");
    fConfig.n_sigma_cut = ellipseCutConfig.at("n_sigma_cut");
  }
  catch (json::out_of_range &e)
  {
    std::cerr << "Missing cut parameter: " << e.what() << std::endl;
    return;
  }

  cutConfigFile.close();
}

void RegenerationFractionFit::_FillHistograms()
{
  fReader->SetEntry(0);
  _ResetHistograms();
  std::size_t invalidSignalWeightCount = 0;
  std::size_t excessiveSignalWeightErrorCount = 0;
  std::size_t nonFiniteFinalSignalWeightCount = 0;
  std::map<std::string, std::size_t> excessiveSignalWeightErrorCountByComponent;
  std::map<std::string, std::size_t> crossChannelMeanErrorCountByComponent;
  double threePi0KsPmMeanWeightSum = 0.0, semileptonicKs00MeanWeightSum = 0.0;
  std::size_t threePi0MeanSampleCount = 0, semileptonicMeanSampleCount = 0;
  const auto evaluateCrossChannelMeanWeight =
      [this, &crossChannelMeanErrorCountByComponent](const std::string &name, double dt, double radius)
  {
    const auto histogram = signalWeightHistograms.find(name);
    if (histogram == signalWeightHistograms.end())
      throw std::runtime_error("Signal weight component is not loaded for cross-channel mean: " + name);
    const auto spline = signalWeightCorrectionSplines.find(name);
    double correction = spline != signalWeightCorrectionSplines.end()
                            ? evaluateCorrectionSpline(*spline->second, dt)
                            : histogram->second->Interpolate(dt);
    const auto ratioSpline = signalWeightRatios.find(name);
    const double ratio = ratioSpline != signalWeightRatios.end() ? ratioSpline->second->Eval(radius) : 1.0;
    if (fConfig.signalWeightRelativeErrorLimit > 0.0)
    {
      const int bin = histogram->second->FindFixBin(dt);
      const double content = histogram->second->GetBinContent(bin);
      const double error = histogram->second->GetBinError(bin);
      const double relativeError = content > 0.0 && std::isfinite(content) &&
                                           std::isfinite(error) && error >= 0.0
                                       ? error / content
                                       : std::numeric_limits<double>::infinity();
      if (!std::isfinite(relativeError) || relativeError > fConfig.signalWeightRelativeErrorLimit)
      {
        correction = 1.0;
        ++crossChannelMeanErrorCountByComponent[name];
      }
    }
    const double weight = ratio * correction;
    return std::isfinite(weight) && weight > 0.0 ? weight : 1.0;
  };

  while (fReader->Next())
  {
    auto &mctruth = **fReaderInteger.at("mctruth");
    auto &tch_mc = **fReaderFloat.at("KaonChTimeCMMC");
    auto &tne_mc = **fReaderFloat.at("KaonNeTimeCMMC");

    auto &tch_fit = **fReaderFloat.at("KaonChTimeCMSignalFit");
    auto &tne_fit = **fReaderFloat.at("KaonNeTimeCMSignalFit");

    auto &tch_tri = **fReaderFloat.at("KaonChTimeCMBoostTriFit");
    auto &tne_tri = **fReaderFloat.at("KaonNeTimeCMBoostTriFit");

    auto &minv4gam = **fReaderFloat.at("minv4gam");
    auto &Kchrec = *fReaderFloatArray.at("Kchrec");
    auto &ip = *fReaderFloatArray.at("ip");

    if (Kchrec.GetSize() < 9 || ip.GetSize() < 3)
      throw std::runtime_error("Signal weight input requires Kchrec[9] and ip[3]");

    double R_ch = sqrt(pow(Kchrec[6] - ip[0],2) + pow(Kchrec[7] - ip[1],2) + pow(Kchrec[8] - ip[2],2));
    double dttri = tch_tri - tne_tri;
    double dtfit = tch_fit - tne_fit;

    double interpWeight = 1.0;
    double chargedModeWeight = 0.0, neutralModeWeight = 0.0;
    bool hasChargedMode = false, hasNeutralMode = false;
    bool excessiveSignalWeightError = false;
    std::vector<std::string> excessiveSignalWeightComponents;
    if (mctruth == 1)
    {
      interpWeight = 0.0;
      for (const auto &component : signalWeightHistograms)
      {
        if (!fConfig.signalWeightComponents.at(component.first).enabled)
          continue;
        const auto correctionSpline = signalWeightCorrectionSplines.find(component.first);
        double correction = correctionSpline != signalWeightCorrectionSplines.end()
                                ? evaluateCorrectionSpline(*correctionSpline->second, dtfit)
                                : component.second->Interpolate(dtfit);
        const auto ratioSpline = signalWeightRatios.find(component.first);
        const double ratio = ratioSpline != signalWeightRatios.end() ? ratioSpline->second->Eval(R_ch) : 1.0;

        if (fConfig.signalWeightRelativeErrorLimit > 0.0)
        {
          const int bin = component.second->FindFixBin(dtfit);
          const double content = component.second->GetBinContent(bin);
          const double error = component.second->GetBinError(bin);
          const double relativeError = content > 0.0 && std::isfinite(content) &&
                                                   std::isfinite(error) && error >= 0.0
                                           ? error / content
                                           : std::numeric_limits<double>::infinity();
          if (!std::isfinite(relativeError) || relativeError > fConfig.signalWeightRelativeErrorLimit)
          {
            correction = 1.0;
            excessiveSignalWeightError = true;
            excessiveSignalWeightComponents.push_back(component.first);
          }
        }

        const double componentWeight = ratio * correction;
        if (component.first == "ks_pm" || component.first == "kl_pm")
        {
          chargedModeWeight += componentWeight;
          hasChargedMode = true;
        }
        else
        {
          neutralModeWeight += componentWeight;
          hasNeutralMode = true;
        }
      }

      if (!hasChargedMode)
        chargedModeWeight = 1.0;
      if (!hasNeutralMode)
        neutralModeWeight = 1.0;
      interpWeight = chargedModeWeight * neutralModeWeight;
    }

    bool invalidSignalWeight = mctruth == 1 && (!std::isfinite(interpWeight) || interpWeight <= 0.0);
    if (invalidSignalWeight)
      interpWeight = 1.0;

    double signalWeight = interpWeight,
           regenerationWeight = 1.0,
           regenerationWeightError = 0.0;

    if (!_isAcceptedEvent())
      continue;

    if (!_applyGlobalCut())
      continue;

    std::string channelName_str = (std::string)KLOE::channName.at(mctruth);

    if (channelName_str == "Signal")
    {
      _sumSignalInterpWeights += interpWeight;
      if (invalidSignalWeight)
        ++invalidSignalWeightCount;
      if (excessiveSignalWeightError)
      {
        ++excessiveSignalWeightErrorCount;
        for (const auto &component : excessiveSignalWeightComponents)
          ++excessiveSignalWeightErrorCountByComponent[component];
      }
    }
    else if (channelName_str == "3pi0" && fConfig.scaleThreePi0ByKsPmMean)
    {
      double tmpWeight = evaluateCrossChannelMeanWeight("ks_pm", dtfit, R_ch) <= 1.0 ? evaluateCrossChannelMeanWeight("ks_pm", dtfit, R_ch) : 1.0;
      threePi0KsPmMeanWeightSum += tmpWeight;
      ++threePi0MeanSampleCount;
    }
    else if (channelName_str == "Semileptonic" && fConfig.scaleSemileptonicByKs00Mean)
    {
      semileptonicKs00MeanWeightSum += evaluateCrossChannelMeanWeight("ks_00", dtfit, R_ch);
      ++semileptonicMeanSampleCount;
    }

    if (fConfig.isSignalWeighted)
    {
      double dtmc = tch_mc - tne_mc;
      signalWeight *= fInterf.fit_function(dtmc, 1);
      if (channelName_str == "Signal" && !std::isfinite(signalWeight))
      {
        signalWeight = 1.0;
        ++nonFiniteFinalSignalWeightCount;
      }
    }

    for (const auto &histType : _histTypesIterative)
    {
      double radius = _calculateRadius(histType);

      double integralWeights = fRegenerationWeights[histType]->Integral();

      if (fConfig.isRegenerationWeighted && integralWeights > 0)
      {
        regenerationWeight = GetRegenerationWeight(histType, radius).first;
      }

      if (channelName_str == "Signal")
      {
        fHistos[histType][channelName_str]->Fill(radius, signalWeight);
      }
      else if (channelName_str == "Regeneration")
      {
        fHistos[histType][channelName_str]->Fill(radius, regenerationWeight);
      }
      else if (channelName_str == "3pi0")
      {
        fHistos[histType][channelName_str]->Fill(radius, avgThreePi0Weight);//threePi0WeightsHist->Interpolate(Kchrec[5]));
      }
      else if (channelName_str == "Semileptonic")
      {
        fHistos[histType][channelName_str]->Fill(radius, avgSemileptonicWeight);//semileptonicWeightsHist->Interpolate(minv4gam));
      }
      else
      {
        fHistos[histType][channelName_str]->Fill(radius);
      }
    }

    regenerationWeight = 1.0;

    // Calculate radius for time difference histogram
    double rhoCharged = _calculateRadius(HistogramType::RHO_CHARGED);
    double rhoNeutral = _calculateRadius(HistogramType::RHO_NEUTRAL);
    double RCharged = _calculateRadius(HistogramType::R_CHARGED);
    double RNeutral = _calculateRadius(HistogramType::R_NEUTRAL);

    if (fConfig.isRegenerationWeighted)
    {
      HistogramType histType = _checkRegenerationLimitsTypesDeltaT(tch_fit - tne_fit);

      if (histType == HistogramType::RHO_CHARGED)
        regenerationWeight = GetContinuousRegenerationWeight(tch_fit - tne_fit, rhoCharged);
      else if (histType == HistogramType::R_CHARGED)
        regenerationWeight = GetContinuousRegenerationWeight(tch_fit - tne_fit, RCharged);
      else if (histType == HistogramType::RHO_NEUTRAL)
        regenerationWeight = GetContinuousRegenerationWeight(tch_fit - tne_fit, rhoNeutral);
      else if (histType == HistogramType::R_NEUTRAL)
        regenerationWeight = GetContinuousRegenerationWeight(tch_fit - tne_fit, RNeutral);
    }

    if (channelName_str == "Signal")
    {
      fTimeDiffHist[channelName_str]->Fill(tch_fit - tne_fit, signalWeight);
      fRhovsRChargedHist[channelName_str]->Fill(RCharged, rhoCharged, signalWeight);
      fRhovsRNeutralHist[channelName_str]->Fill(RNeutral, rhoNeutral, signalWeight);
    }
    else if (channelName_str == "Regeneration")
    {
      fTimeDiffHist[channelName_str]->Fill(tch_fit - tne_fit, regenerationWeight);
      fRhovsRChargedHist[channelName_str]->Fill(RCharged, rhoCharged, regenerationWeight);
      fRhovsRNeutralHist[channelName_str]->Fill(RNeutral, rhoNeutral, regenerationWeight);
    }
    else if (channelName_str == "3pi0")
    {
      fTimeDiffHist[channelName_str]->Fill(tch_fit - tne_fit, avgThreePi0Weight);//threePi0WeightsHist->Interpolate(Kchrec[5]));
      fRhovsRChargedHist[channelName_str]->Fill(RCharged, rhoCharged, avgThreePi0Weight);//threePi0WeightsHist->Interpolate(Kchrec[5]));
      fRhovsRNeutralHist[channelName_str]->Fill(RNeutral, rhoNeutral, avgThreePi0Weight);//threePi0WeightsHist->Interpolate(Kchrec[5]));
    }
    else if (channelName_str == "Semileptonic")
    {
      fTimeDiffHist[channelName_str]->Fill(tch_fit - tne_fit,  avgSemileptonicWeight);//semileptonicWeightsHist->Interpolate(minv4gam));
      fRhovsRChargedHist[channelName_str]->Fill(RCharged, rhoCharged, avgSemileptonicWeight);//semileptonicWeightsHist->Interpolate(minv4gam));
      fRhovsRNeutralHist[channelName_str]->Fill(RNeutral, rhoNeutral, avgSemileptonicWeight);//semileptonicWeightsHist->Interpolate(minv4gam));
    }
    else
    {
      fTimeDiffHist[channelName_str]->Fill(tch_fit - tne_fit);
      fRhovsRChargedHist[channelName_str]->Fill(RCharged, rhoCharged);
      fRhovsRNeutralHist[channelName_str]->Fill(RNeutral, rhoNeutral);
    }
  }

  const auto scaleChannel = [this](const std::string &channel, double factor)
  {
    for (const auto &histType : _histTypesIterative)
      fHistos[histType][channel]->Scale(factor);
    fTimeDiffHist[channel]->Scale(factor);
    fRhovsRChargedHist[channel]->Scale(factor);
    fRhovsRNeutralHist[channel]->Scale(factor);
  };
  if (fConfig.scaleThreePi0ByKsPmMean || fConfig.scaleSemileptonicByKs00Mean)
  {
    if (fConfig.scaleThreePi0ByKsPmMean)
    {
      if (threePi0MeanSampleCount == 0)
        throw std::runtime_error("Cannot calculate ks_pm mean weight: no selected 3pi0 events");
      const double mean = threePi0KsPmMeanWeightSum / threePi0MeanSampleCount;
      if (!std::isfinite(mean) || mean <= 0.0)
        throw std::runtime_error("Invalid selected-3pi0 mean weight for ks_pm");
      scaleChannel("3pi0", mean);
      std::cout << "Applying selected-3pi0 ks_pm mean weight " << mean << std::endl;
    }
    if (fConfig.scaleSemileptonicByKs00Mean)
    {
      if (semileptonicMeanSampleCount == 0)
        throw std::runtime_error("Cannot calculate ks_00 mean weight: no selected Semileptonic events");
      const double mean = semileptonicKs00MeanWeightSum / semileptonicMeanSampleCount;
      if (!std::isfinite(mean) || mean <= 0.0)
        throw std::runtime_error("Invalid selected-Semileptonic mean weight for ks_00");
      scaleChannel("Semileptonic", mean);
      std::cout << "Applying selected-Semileptonic ks_00 mean weight " << mean << std::endl;
    }
    for (const auto &component : crossChannelMeanErrorCountByComponent)
      std::cout << "Cross-channel " << component.first << " error substitutions=" << component.second << std::endl;
  }

  std::cout << "Signal weight corrections substituted with 1: invalid composed weights="
            << invalidSignalWeightCount
            << ", events with high-error component substitutions=" << excessiveSignalWeightErrorCount;
  for (const auto &component : excessiveSignalWeightErrorCountByComponent)
    std::cout << ", " << component.first << "=" << component.second;
  std::cout
            << ", non-finite after interference=" << nonFiniteFinalSignalWeightCount
            << " (limit=" << fConfig.signalWeightRelativeErrorLimit << ")" << std::endl;

  if (fConfig.isSignalWeighted)
  {
    _RenormalizeSignalHistogram();
  }
}

void RegenerationFractionFit::_RenormalizeMCToLumi()
{
  KLOE::BRCorrectionFactors factors;

  for (const auto &histType : _histTypesIterative)
  {
    for (const auto &channel : KLOE::channName)
    {
      std::string channelName_str = (std::string)channel.second;

      if (channelName_str == "Data" || channelName_str == "MC sum")
        continue;

      fHistos[histType][channelName_str]->Scale(fConfig.lumiFactor * factors.BRcorrectionFactors.at(channelName_str));
    }
  }

  for (const auto &channel : KLOE::channName)
  {
    std::string channelName_str = (std::string)channel.second;

    if (channelName_str == "Data" || channelName_str == "MC sum")
      continue;

    fTimeDiffHist[channelName_str]->Scale(fConfig.lumiFactor * factors.BRcorrectionFactors.at(channelName_str));
    fRhovsRChargedHist[channelName_str]->Scale(fConfig.lumiFactor * factors.BRcorrectionFactors.at(channelName_str));
    fRhovsRNeutralHist[channelName_str]->Scale(fConfig.lumiFactor * factors.BRcorrectionFactors.at(channelName_str));
  }
}

void RegenerationFractionFit::_ResetHistograms()
{
  for (const auto &histType : _histTypesIterative)
  {
    for (const auto &channel : KLOE::channName)
    {
      std::string channelName_str = (std::string)channel.second;
      fHistos[histType][channelName_str]->Reset("ICESM");
    }

    _channInMCSum[histType].clear(); // Clear the list of channels included in MC sum for this histogram type
  }

  for (const auto &channel : KLOE::channName)
  {
    std::string channelName_str = (std::string)channel.second;
    fTimeDiffHist[channelName_str]->Reset("ICESM");
    fRhovsRChargedHist[channelName_str]->Reset("ICESM");
    fRhovsRNeutralHist[channelName_str]->Reset("ICESM");
  }

  _channInTimeDiff.clear(); // Clear the list of channels included in time difference sum
  _sumSignalInterpWeights = 0.0;
}

void RegenerationFractionFit::_FillMCSumHistogram(bool addRegeneration)
{

  for (const auto &histType : _histTypesIterative)
  {
    for (const auto &channel : KLOE::channName)
    {
      std::string channelName_str = (std::string)channel.second;

      if (channelName_str == "Data" || channelName_str == "MC sum")
        continue;

      bool channelInSum = std::find(_channInMCSum[histType].begin(),
                                    _channInMCSum[histType].end(), channelName_str) !=
                          _channInMCSum[histType].end();

      if (channelInSum) // Avoid adding the same channel multiple times if this function is called more than once
        continue;

      if (!addRegeneration && channelName_str == "Regeneration")
        continue;

      fHistos[histType][channelName_str]->Sumw2();
      fHistos[histType]["MC sum"]->Add(fHistos[histType][channelName_str]);

      _channInMCSum[histType].push_back(channelName_str);
    }
  }

  for (const auto &channel : KLOE::channName)
  {
    std::string channelName_str = (std::string)channel.second;

    if (channelName_str == "Data" || channelName_str == "MC sum")
      continue;

    bool channelInSum = std::find(_channInTimeDiff.begin(),
                                  _channInTimeDiff.end(), channelName_str) !=
                        _channInTimeDiff.end();

    if (channelInSum) // Avoid adding the same channel multiple times if this function is called more than once
      continue;

    if (!addRegeneration && channelName_str == "Regeneration")
      continue;

    fTimeDiffHist[channelName_str]->Sumw2();
    fTimeDiffHist["MC sum"]->Add(fTimeDiffHist[channelName_str]);

    _channInTimeDiff.push_back(channelName_str); // Using RHO_CHARGED just to track which channels have been added to the time difference sum

    fRhovsRChargedHist[channelName_str]->Sumw2();
    fRhovsRChargedHist["MC sum"]->Add(fRhovsRChargedHist[channelName_str]);

    fRhovsRNeutralHist[channelName_str]->Sumw2();
    fRhovsRNeutralHist["MC sum"]->Add(fRhovsRNeutralHist[channelName_str]);
  }
}

void RegenerationFractionFit::DrawHistograms() const
{
  std::string outputDir = (std::string)Paths::regen_analysis_dir + (std::string)Paths::img_dir;
  gStyle->SetOptStat(0);
  for (const auto &histType : _histTypesIterative)
  {
    TCanvas canvas;
    canvas.SetLogy(fConfig.logarithmic.at(histType));
    TLegend legend(0.65, 0.6, 0.85, 0.9);

    if (fConfig.logarithmic.at(histType))
    {
      fHistos.at(histType).at("Data")->SetMinimum(0.1);
    }
    else
    {
      fHistos.at(histType).at("Data")->SetMinimum(0.0);
    }

    fHistos.at(histType).at("Data")->SetMarkerColor(KLOE::channColor.at("Data"));
    fHistos.at(histType).at("Data")->SetLineColor(KLOE::channColor.at("Data"));
    fHistos.at(histType).at("Data")->Draw("E1");

    legend.AddEntry(fHistos.at(histType).at("Data"), "Data", "lep");

    for (const auto &channel : KLOE::channName)
    {
      std::string channelName_str = (std::string)channel.second;

      if (channelName_str == "Data")
        continue;

      fHistos.at(histType).at(channelName_str)->SetLineColor(KLOE::channColor.at(channelName_str));
      fHistos.at(histType).at(channelName_str)->Draw("HIST SAME");
      legend.AddEntry(fHistos.at(histType).at(channelName_str), channelName_str.c_str(), "l");
    }

    legend.Draw();
    std::string outputFilePath = outputDir + "histogram_" + _HistTypeToString(histType) + ".svg";
    canvas.SaveAs(outputFilePath.c_str());
  }

  // Draw time difference histograms
  TCanvas timeDiffCanvas;
  TLegend timeDiffLegend(0.65, 0.6, 0.85, 0.9);

  fTimeDiffHist.at("Data")->SetMinimum(0.0);
  fTimeDiffHist.at("Data")->SetMarkerColor(KLOE::channColor.at("Data"));
  fTimeDiffHist.at("Data")->SetLineColor(KLOE::channColor.at("Data"));
  fTimeDiffHist.at("Data")->SetTitle("Time Difference (t_{ch} - t_{ne});t_{ch} - t_{ne} [#tau_{S}];Entries / 1.5 #tau_{S}");
  fTimeDiffHist.at("Data")->GetYaxis()->SetRangeUser(0.0, fTimeDiffHist.at("Data")->GetMaximum() * 1.5);
  fTimeDiffHist.at("Data")->Draw("E1");

  timeDiffLegend.AddEntry(fTimeDiffHist.at("Data"), "Data", "lep");

  for (const auto &channel : KLOE::channName)
  {
    std::string channelName_str = (std::string)channel.second;

    if (channelName_str == "Data")
      continue;

    fTimeDiffHist.at(channelName_str)->SetLineColor(KLOE::channColor.at(channelName_str));
    fTimeDiffHist.at(channelName_str)->Draw("HIST SAME");
    timeDiffLegend.AddEntry(fTimeDiffHist.at(channelName_str), channelName_str.c_str(), "l");
  }

  timeDiffLegend.Draw();
  std::string timeDiffOutputFilePath = outputDir + "time_difference_histogram.svg";
  timeDiffCanvas.SaveAs(timeDiffOutputFilePath.c_str());

  // Draw regeneration weight histograms
  for (const auto &histType : _histTypesIterative)
  {
    TCanvas weightCanvas;

    fRegenerationWeights.at(histType)->SetLineColor(kBlack);
    fRegenerationWeights.at(histType)->SetMarkerColor(kBlack);
    fRegenerationWeights.at(histType)->SetMinimum(0.0);
    fRegenerationWeights.at(histType)->SetMaximum(7.0);
    fRegenerationWeights.at(histType)->SetTitle(("Regeneration Weights for " + fConfig.titles.at(histType)).c_str());
    fRegenerationWeights.at(histType)->GetXaxis()->SetTitle(fConfig.xAxisTitles.at(histType).c_str());
    fRegenerationWeights.at(histType)->GetYaxis()->SetTitle("Weight [-]");
    fRegenerationWeights.at(histType)->Draw("PE1");

    TLine lineDCLow(fConfig.DCWallBottomBoundCharged, 0, fConfig.DCWallBottomBoundCharged, fRegenerationWeights.at(histType)->GetMaximum());
    TLine lineDCHigh(fConfig.DCWallTopBoundCharged, 0, fConfig.DCWallTopBoundCharged, fRegenerationWeights.at(histType)->GetMaximum());
    lineDCLow.SetLineColor(kRed);
    lineDCHigh.SetLineColor(kRed);
    lineDCLow.SetLineStyle(2);

    TLine lineBPCylBottomCharged(fConfig.cylindricalBPBottomBoundCharged, 0, fConfig.cylindricalBPBottomBoundCharged, fRegenerationWeights.at(histType)->GetMaximum());
    TLine lineBPCylTopCharged(fConfig.cylindricalBPTopBoundCharged, 0, fConfig.cylindricalBPTopBoundCharged, fRegenerationWeights.at(histType)->GetMaximum());
    lineBPCylBottomCharged.SetLineColor(kBlue);
    lineBPCylTopCharged.SetLineColor(kBlue);
    lineBPCylBottomCharged.SetLineStyle(2);
    lineBPCylTopCharged.SetLineStyle(2);

    TLine lineBPSphBottomCharged(fConfig.sphericalBPBottomBoundCharged, 0, fConfig.sphericalBPBottomBoundCharged, fRegenerationWeights.at(histType)->GetMaximum());
    TLine lineBPSphTopCharged(fConfig.sphericalBPTopBoundCharged, 0, fConfig.sphericalBPTopBoundCharged, fRegenerationWeights.at(histType)->GetMaximum());
    lineBPSphBottomCharged.SetLineColor(kGreen);
    lineBPSphTopCharged.SetLineColor(kGreen);
    lineBPSphBottomCharged.SetLineStyle(2);
    lineBPSphTopCharged.SetLineStyle(2);

    lineDCLow.Draw("SAME");
    lineDCHigh.Draw("SAME");
    lineBPCylBottomCharged.Draw("SAME");
    lineBPCylTopCharged.Draw("SAME");
    lineBPSphBottomCharged.Draw("SAME");
    lineBPSphTopCharged.Draw("SAME");

    std::string outputFilePath = outputDir + "regeneration_weights_" + _HistTypeToString(histType) + ".svg";
    weightCanvas.SaveAs(outputFilePath.c_str());
  }

  // Draw time difference histograms

  for (const auto &channel : KLOE::channName)
  {
    std::string channelName_str = (std::string)channel.second;
    TCanvas *rhoVsRChargedCanvas = new TCanvas(("RhoVsRCharged_" + channelName_str).c_str(), ("#rho vs R Charged for " + channelName_str).c_str(), 1200, 1200);

    fRhovsRChargedHist.at(channelName_str)->Draw("COLZ");

    std::string rhovsRChargedPath = outputDir + "rho_vs_r_charged_" + channelName_str + ".svg";
    rhoVsRChargedCanvas->SaveAs(rhovsRChargedPath.c_str());
  }

  for (const auto &channel : KLOE::channName)
  {
    std::string channelName_str = (std::string)channel.second;
    TCanvas *rhoVsRNeutralCanvas = new TCanvas(("RhoVsRNeutral_" + channelName_str).c_str(), ("#rho vs R Neutral for " + channelName_str).c_str(), 1200, 1200);

    fRhovsRNeutralHist.at(channelName_str)->Draw("COLZ");

    std::string rhovsRNeutralPath = outputDir + "rho_vs_r_neutral_" + channelName_str + ".svg";
    rhoVsRNeutralCanvas->SaveAs(rhovsRNeutralPath.c_str());
  }

}

void RegenerationFractionFit::_RenormalizeSignalHistogram()
{
  for (const auto &histType : _histTypesIterative)
  {
    int signalBins = fHistos[histType]["Signal"]->GetNbinsX();
    double signalIntegral = fHistos[histType]["Signal"]->Integral(0, signalBins + 1);
    double signalEntries = fHistos[histType]["Signal"]->GetEntries();

    if (signalIntegral > 0)
    {
      double scaleFactor = (_sumSignalInterpWeights > 0)
                           ? _sumSignalInterpWeights / signalIntegral
                           : signalEntries / signalIntegral;
      fHistos[histType]["Signal"]->Scale(scaleFactor);
    }
  }

  // Renormalize the time difference histogram for the "Signal" channel
  int signalTimeDiffBins = fTimeDiffHist["Signal"]->GetNbinsX();
  double signalTimeDiffIntegral = fTimeDiffHist["Signal"]->Integral(0, signalTimeDiffBins + 1);
  double signalTimeDiffEntries = fTimeDiffHist["Signal"]->GetEntries();

  if (signalTimeDiffIntegral > 0)
  {
    double scaleFactor = (_sumSignalInterpWeights > 0)
                         ? _sumSignalInterpWeights / signalTimeDiffIntegral
                         : signalTimeDiffEntries / signalTimeDiffIntegral;
    fTimeDiffHist["Signal"]->Scale(scaleFactor);
  }
}

void RegenerationFractionFit::_calculateRegenerationWeights()
{
  // 1. Napełniamy surowe histogramy (bez wag i bez lumi)
  _FillHistograms();
  _RenormalizeMCToLumi();
  _FillMCSumHistogram(false); // Include Regeneration in the MC sum
  
  for (const auto &histType : _histTypesIterative)
  {
    TH1 *dataRegenerationHist = (TH1*)(fHistos[histType]["Data"]->Clone("DataRegeneration"));
    dataRegenerationHist->Add(fHistos[histType]["MC sum"], -1);

    fRegenerationWeights[histType]->Reset("ICESM");
    fRegenerationWeights[histType]->Divide(dataRegenerationHist, fHistos[histType]["Regeneration"], 1.0, 1.0);
      
    for (int i = 1; i <= fRegenerationWeights[histType]->GetNbinsX(); ++i)
    {
      double weight = fRegenerationWeights[histType]->GetBinContent(i);
      double weightError = fRegenerationWeights[histType]->GetBinError(i);
      double relativeError = (weight != 0) ? weightError / weight : 0.0;

      std::cout << "Histogram Type: " << _HistTypeToString(histType) << ", Bin: " << i
                << ", Weight: " << weight << ", Error: " << weightError
                << ", Relative Error: " << relativeError << std::endl;

      if (weight <= 0.0 || relativeError > fConfig.weightErrorLimit)
      {
        fRegenerationWeights[histType]->SetBinContent(i, 0.0);
        fRegenerationWeights[histType]->SetBinError(i, 0.0);
      }
    }

    dataRegenerationHist->Delete();
  }
}

std::pair<double, double> RegenerationFractionFit::GetRegenerationWeight(HistogramType histType, double radius) const
{
  int bin = fRegenerationWeights.at(histType)->FindBin(radius);
  double xLeft = fRegenerationWeights.at(histType)->GetBinLowEdge(bin);
  double xRight = fRegenerationWeights.at(histType)->GetBinLowEdge(bin + 1);

  if (!_checkRegenerationLimits(xLeft, xRight))
    return std::make_pair(1.0, 0.0); // Default weight if the radius is outside the defined limits

  double weight = fRegenerationWeights.at(histType)->GetBinContent(bin);
  double weightError = fRegenerationWeights.at(histType)->GetBinError(bin);
  double relativeError = (weight != 0) ? weightError / weight : 0.0;

  if (relativeError > fConfig.weightErrorLimit)
    return std::make_pair(1.0, 0.0); // Default weight if the relative error exceeds the configured limit

  if (weight <= 0.0)
    return std::make_pair(1.0, 0.0); // Default weight if the calculated weight is zero or negative

  return std::make_pair(weight, weightError);
}

std::pair<double, double> RegenerationFractionFit::GetRegenerationWeight(double radius) const
{
  HistogramType histType = _checkRegenerationLimitsTypes(radius);
  if (histType == HistogramType::UNKNOWN)
    return std::make_pair(1.0, 0.0); // Default weight if the radius does not fall within any defined histogram type

  return GetRegenerationWeight(histType, radius);
}

bool RegenerationFractionFit::_checkRegenerationLimits(double xLeft, double xRight) const
{
  bool isInsideDCWallBounds = (xLeft >= fConfig.DCWallBottomBoundCharged && xRight <= fConfig.DCWallTopBoundCharged) ||
                              (xLeft >= fConfig.DCWallBottomBoundNeutral && xRight <= fConfig.DCWallTopBoundNeutral);

  bool isInsideSphericalBPBounds = (xLeft >= fConfig.sphericalBPBottomBoundCharged &&
                                    xRight <= fConfig.sphericalBPTopBoundCharged) ||
                                   (xLeft >= fConfig.sphericalBPBottomBoundNeutral &&
                                    xRight <= fConfig.sphericalBPTopBoundNeutral);

  bool isInsideCylindricalBPBounds = (xLeft >= fConfig.cylindricalBPBottomBoundCharged &&
                                      xRight <= fConfig.cylindricalBPTopBoundCharged) ||
                                     (xLeft >= fConfig.cylindricalBPBottomBoundNeutral &&
                                      xRight <= fConfig.cylindricalBPTopBoundNeutral);

  return isInsideDCWallBounds || isInsideSphericalBPBounds || isInsideCylindricalBPBounds;
}

HistogramType RegenerationFractionFit::_checkRegenerationLimitsTypes(double radius) const
{
  bool isInsideChargedCylindricalBP = (radius >= fConfig.cylindricalBPBottomBoundCharged && radius <= fConfig.cylindricalBPTopBoundCharged);
  bool isInsideNeutralCylindricalBP = (radius >= fConfig.cylindricalBPBottomBoundNeutral && radius <= fConfig.cylindricalBPTopBoundNeutral);

  bool isInsideChargedDCWall = (radius >= fConfig.DCWallBottomBoundCharged && radius <= fConfig.DCWallTopBoundCharged);
  bool isInsideNeutralDCWall = (radius >= fConfig.DCWallBottomBoundNeutral && radius <= fConfig.DCWallTopBoundNeutral);

  bool isInsideChargedSphericalBP = (radius >= fConfig.sphericalBPBottomBoundCharged && radius <= fConfig.sphericalBPTopBoundCharged);
  bool isInsideNeutralSphericalBP = (radius >= fConfig.sphericalBPBottomBoundNeutral && radius <= fConfig.sphericalBPTopBoundNeutral);

  if ((isInsideChargedCylindricalBP || isInsideChargedDCWall))
    return HistogramType::RHO_CHARGED; // Assuming rhoCharged corresponds to charged regeneration

  if ((isInsideNeutralCylindricalBP || isInsideNeutralDCWall))
    return HistogramType::RHO_NEUTRAL; // Assuming rhoNeutral corresponds to neutral regeneration

  if (isInsideChargedSphericalBP)
    return HistogramType::R_CHARGED; // Assuming RCharged corresponds to charged regeneration

  if (isInsideNeutralSphericalBP)
    return HistogramType::R_NEUTRAL; // Assuming RNeutral corresponds to neutral regeneration

  return HistogramType::UNKNOWN;
}

HistogramType RegenerationFractionFit::_checkRegenerationLimitsTypesDeltaT(double dt) const
{
  if (dt > -300.0 && dt <= -30.0)
    return HistogramType::RHO_NEUTRAL;
  if (dt > -30.0 && dt <= -8.0)
    return HistogramType::R_NEUTRAL;
  if (dt > -8.0 && dt <= 0.0)
    return HistogramType::RHO_NEUTRAL;
  if (dt > 0.0 && dt <= 8.0)
    return HistogramType::RHO_CHARGED;
  if (dt > 8.0 && dt <= 30.0)
    return HistogramType::R_CHARGED;
  if (dt > 30.0 && dt <= 300.0)
    return HistogramType::RHO_CHARGED;

  return HistogramType::UNKNOWN;
}

std::pair<TF1*, HistogramType> RegenerationFractionFit::_GetFitFuncByDt(double dt)
{
  if (dt > -300.0 && dt <= -30.0)
    return {fcontWeightDC[HistogramType::RHO_NEUTRAL], HistogramType::RHO_NEUTRAL};
  if (dt > -30.0 && dt <= -12.0)
    return {fcontWeightSphBP[HistogramType::R_NEUTRAL], HistogramType::R_NEUTRAL};
  if (dt > -12.0 && dt <= 0.0)
    return {fcontWeightCylBP[HistogramType::RHO_NEUTRAL], HistogramType::RHO_NEUTRAL};
  if (dt > 0.0 && dt <= 12.0)
    return {fcontWeightCylBP[HistogramType::RHO_CHARGED], HistogramType::RHO_CHARGED};
  if (dt > 12.0 && dt <= 30.0)
    return {fcontWeightSphBP[HistogramType::R_CHARGED], HistogramType::R_CHARGED};
  if (dt > 30.0 && dt <= 300.0)
    return {fcontWeightDC[HistogramType::RHO_CHARGED], HistogramType::RHO_CHARGED};

  return {nullptr, HistogramType::UNKNOWN};
}

void RegenerationFractionFit::SaveHistograms()
{
  // 1. Uruchomienie procedury wyliczenia wag na bazie surowego MC pomnożonego w locie
  _calculateRegenerationWeights();

  // 2. Ponowne napełnienie pętli zdarzeń – tym razem z zaaplikowaniem świeżo wyliczonych wag
  if (fConfig.isRegenerationWeighted)
  {
    _FillHistograms();
  }

  _RenormalizeMCToLumi();
  // 4. Budujemy finalny sumaryczny histogram "MC sum" (zawierający już wagi i lumi) do pliku ROOT
  _FillMCSumHistogram(true);

  // 4b. Chi2 Data vs MC sum dla wszystkich rysowanych histogramów
  std::cout << "\n--- Chi2 (Data vs MC sum) after regeneration weighting ---" << std::endl;
  for (const auto &histType : _histTypesIterative)
  {
    TH1 *hData = fHistos.at(histType).at("Data");
    TH1 *hMC   = fHistos.at(histType).at("MC sum");
    if (hData->GetEntries() > 0 && hMC->GetEntries() > 0)
    {
      Double_t chi2 = 0.0;
      Int_t ndf = 0;
      Int_t igood = 0;
      hData->Chi2TestX(hMC, chi2, ndf, igood, "WW");
      double chi2ndf = (ndf > 0) ? chi2 / ndf : 0.0;
      std::cout << "  " << _HistTypeToString(histType)
                << ":  chi2 = " << chi2
                << ",  ndf = " << ndf
                << ",  chi2/ndf = " << chi2ndf << std::endl;
    }
    else
    {
      std::cout << "  " << _HistTypeToString(histType) << ":  skipped (empty histogram)" << std::endl;
    }
  }
  std::cout << "----------------------------------------------------------\n" << std::endl;

  // 5. Zapis struktur do pliku ROOT
  std::string outputFilePath = (std::string)Paths::regen_analysis_dir + "results/regeneration_analysis_results.root";
  TFile outputFile(outputFilePath.c_str(), "RECREATE");

  for (const auto &histType : _histTypesIterative)
  {
    for (const auto &channel : KLOE::channName)
    {
      std::string channelName_str = (std::string)channel.second;
      fHistos.at(histType).at(channelName_str)->Write();
    }
    fRegenerationWeights.at(histType)->Write();
  }

  // Zapisujemy również nowo dodane widma różnic czasu
  for (const auto &channel : KLOE::channName)
  {
    std::string channelName_str = (std::string)channel.second;
    fTimeDiffHist.at(channelName_str)->Write();
  }

  for (const auto &histType : _histTypesIterative)
  {
    if (fcontWeightDC[histType] != nullptr)
    {
      fcontWeightDC[histType]->Write();
    }
    if (fcontWeightCylBP[histType] != nullptr)
    {
      fcontWeightCylBP[histType]->Write();
    }
    if (fcontWeightSphBP[histType] != nullptr)
    {
      fcontWeightSphBP[histType]->Write();
    }
  }

  outputFile.Close();
}

// ----

void RegenerationFractionFit::_SetupJson(json &inputConfig, json &outputConfig, std::string section) const
{
  try
  {
    outputConfig = inputConfig.at(section);
  }
  catch (json::out_of_range &e)
  {
    std::cerr << "Missing section " << section << ": " << e.what() << std::endl;
    return;
  }
};

void RegenerationFractionFit::_SetupReaderVariables()
{
  for (const auto &var : floatVars)
  {
    fReaderFloat[var] = new TTreeReaderValue<Double_t>(*fReader, var.c_str());
  }
  for (const auto &var : integerVars)
  {
    fReaderInteger[var] = new TTreeReaderValue<Int_t>(*fReader, var.c_str());
  }
  for (const auto &var : floatArrayVars)
  {
    fReaderFloatArray[var] = new TTreeReaderArray<Double_t>(*fReader, var.c_str());
  }
}

std::string RegenerationFractionFit::_HistTypeToString(HistogramType type) const
{
  switch (type)
  {
  case HistogramType::RHO_CHARGED:
    return "rhoCharged";
  case HistogramType::RHO_NEUTRAL:
    return "rhoNeutral";
  case HistogramType::R_CHARGED:
    return "RCharged";
  case HistogramType::R_NEUTRAL:
    return "RNeutral";
  default:
    return "Unknown";
  }
}

HistogramType RegenerationFractionFit::_StringToHistType(const std::string &str) const
{
  if (str == "rhoCharged")
    return HistogramType::RHO_CHARGED;
  else if (str == "rhoNeutral")
    return HistogramType::RHO_NEUTRAL;
  else if (str == "RCharged")
    return HistogramType::R_CHARGED;
  else if (str == "RNeutral")
    return HistogramType::R_NEUTRAL;
  else
    return HistogramType::UNKNOWN; // Default value
}

// ---

bool RegenerationFractionFit::_isAcceptedEvent() const
{
  auto &mcflag = **fReaderInteger.at("mcflag");
  auto &mctruth = **fReaderInteger.at("mctruth");

  return (mcflag == 0 || (mcflag == 1 && mctruth != -1 && mctruth != 0));
}

bool RegenerationFractionFit::_isOutsideEllipse(double u, double v) const
{
  double du = u - fConfig.u0;
  double dv = v - fConfig.v0;
  double inv_rho2 = 1.0 / (1.0 - fConfig.rho * fConfig.rho);

  double dist2 = inv_rho2 * ((du * du) / (fConfig.su * fConfig.su) +
                             (dv * dv) / (fConfig.sv * fConfig.sv) -
                             2.0 * fConfig.rho * (du * dv) / (fConfig.su * fConfig.sv));

  return std::sqrt(std::max(0.0, dist2)) > fConfig.n_sigma_cut;
}

bool RegenerationFractionFit::_applyGlobalCut() const
{
  auto &KchrecFit = *fReaderFloatArray.at("KchrecFit");
  auto &KnerecFit = *fReaderFloatArray.at("KnerecFit");
  auto &Knerec = *fReaderFloatArray.at("Knerec");
  auto &KchrecClosest = *fReaderFloatArray.at("KchrecClosest");
  auto &KnerecSix = *fReaderFloatArray.at("KnerecSix");
  auto &trk1Fit = *fReaderFloatArray.at("trk1Fit");
  auto &trk2Fit = *fReaderFloatArray.at("trk2Fit");

  double Bx = **fReaderFloat.at("Bx");
  double By = **fReaderFloat.at("By");
  double Bz = **fReaderFloat.at("Bz");
  double minv4gam = **fReaderFloat.at("minv4gam");

  // 2. Geometria wierzchołków (Distances)
  double distNC_X = KchrecFit[6] - KnerecFit[6];
  double distNC_Y = KchrecFit[7] - KnerecFit[7];
  double distNC_Z = KchrecFit[8] - KnerecFit[8];

  double radius00 = std::hypot(Knerec[6] - Bx, Knerec[7] - By);
  double radiuspm = std::hypot(KchrecClosest[6] - Bx, KchrecClosest[7] - By);

  double zdist00 = std::abs(Knerec[8] - KchrecClosest[8]);
  double zdistpm = std::abs(KchrecClosest[8] - Bz);

  // 3. Wyliczanie zmiennej RHO (wielowymiarowa odległość od punktu interakcji)
  double rho_pm = std::hypot(KchrecClosest[6] - Bx, KchrecClosest[7] - By);
  double rho_00 = std::hypot(Knerec[6] - Bx, Knerec[7] - By);
  double rho = std::hypot(rho_pm, rho_00);

  // 4. Definicja Fiducial Volumes (Obszarów detektora)
  bool fiducialVolume = std::hypot(distNC_X, distNC_Y) < 2.05 && std::abs(distNC_Z) < 2.45;
  bool fiducialVolumeClose = radius00 < 1.5 && radiuspm < 2.0 && zdist00 < 1.5 && zdistpm < 1.5;

  // 5. Kinematyka i kąty (KLOE Track/Cluster matching)
  TVector3 KchrecVec(KchrecFit[0], KchrecFit[1], KchrecFit[2]);
  TVector3 trk1VecFit(trk1Fit[0], trk1Fit[1], trk1Fit[2]);
  TVector3 trk2VecFit(trk2Fit[0], trk2Fit[1], trk2Fit[2]);

  double phiTrk1Angle = std::cos(trk1VecFit.Angle(KchrecVec));
  double phiTrk2Angle = std::cos(trk2VecFit.Angle(KchrecVec));

  // 6. Warunki logiczne cięć (Cuts criteria)
  bool angleCut = !fiducialVolume || (std::abs(phiTrk1Angle) < 0.8 || std::abs(phiTrk2Angle) < 0.8);
  bool rhoCut = !fiducialVolumeClose || (rho > 1.5);

  // Wywołanie elipsy (zakładam, że parametry u0, v0 są w fConfig, a elipsa to metoda klasy)
  bool ellipseCut = _isOutsideEllipse(minv4gam - PhysicsConstants::mK0, KnerecSix[5] - PhysicsConstants::mK0);

  // Zwracamy finalny werdykt
  return angleCut && rhoCut && ellipseCut;
}

double RegenerationFractionFit::_calculateRadius(HistogramType type) const
{
  auto &Kchrec = *fReaderFloatArray.at("Kchrec");
  auto &Knerec = *fReaderFloatArray.at("Knerec");

  TVector3 KchrecVec(Kchrec[6], Kchrec[7], Kchrec[8]);
  TVector3 KnerecVec(Knerec[6], Knerec[7], Knerec[8]);

  double rho_ch_init = KchrecVec.Perp();
  double rho_ne_init = KnerecVec.Perp();

  if (rho_ch_init > fConfig.DCWallBottomBoundCharged && rho_ch_init < fConfig.DCWallTopBoundCharged)
    KchrecVec -= fConfig.DCWallCenter;

  if (rho_ne_init > fConfig.DCWallBottomBoundNeutral && rho_ne_init < fConfig.DCWallTopBoundNeutral)
    KnerecVec -= fConfig.DCWallCenter;

  switch (type)
  {
  case HistogramType::RHO_CHARGED:
    return KchrecVec.Perp();
  case HistogramType::RHO_NEUTRAL:
    return KnerecVec.Perp();
  case HistogramType::R_CHARGED:
    return KchrecVec.Mag();
  case HistogramType::R_NEUTRAL:
    return KnerecVec.Mag();
  default:
    return -1; // Invalid type
  }
}

void RegenerationFractionFit::_FitContinuousWeightFunction(HistogramType histType)
{
  if (fcontWeightDC.find(histType) != fcontWeightDC.end() && fcontWeightDC[histType] != nullptr)
  {
    delete fcontWeightDC[histType];
    delete fcontWeightSphBP[histType];
    delete fcontWeightCylBP[histType];
  }

  // Assuming fInterf.fit_function is a callable object that takes x and type as parameters
  std::string funcNameDC = "FuncDC_" + _HistTypeToString(histType),
              funcNameSphBP = "FuncSphBP_" + _HistTypeToString(histType),
              funcNameCylBP = "FuncCylBP_" + _HistTypeToString(histType);
  if (histType == HistogramType::RHO_CHARGED || histType == HistogramType::R_CHARGED)
  {
    fcontWeightDC[histType] = new TF1(funcNameDC.c_str(), "pol1", fConfig.DCWallBottomBoundCharged, fConfig.DCWallTopBoundCharged);
    fcontWeightSphBP[histType] = new TF1(funcNameSphBP.c_str(), "pol1", fConfig.sphericalBPBottomBoundCharged, fConfig.sphericalBPTopBoundCharged);
    fcontWeightCylBP[histType] = new TF1(funcNameCylBP.c_str(), "pol1", fConfig.cylindricalBPBottomBoundCharged, fConfig.cylindricalBPTopBoundCharged);
  }
  else if (histType == HistogramType::RHO_NEUTRAL || histType == HistogramType::R_NEUTRAL)
  {
    fcontWeightDC[histType] = new TF1(funcNameDC.c_str(), "pol1", fConfig.DCWallBottomBoundNeutral, fConfig.DCWallTopBoundNeutral);
    fcontWeightSphBP[histType] = new TF1(funcNameSphBP.c_str(), "pol1", fConfig.sphericalBPBottomBoundNeutral, fConfig.sphericalBPTopBoundNeutral);
    fcontWeightCylBP[histType] = new TF1(funcNameCylBP.c_str(), "pol1", fConfig.cylindricalBPBottomBoundNeutral, fConfig.cylindricalBPTopBoundNeutral);
  }
  else
  {
    std::cerr << "Unknown histogram type for fitting continuous weight function." << std::endl;
    return;
  }

  TFitResultPtr fitResultDC = fRegenerationWeights[histType]->Fit(fcontWeightDC[histType], "RES");
  TFitResultPtr fitResultSphBP = fRegenerationWeights[histType]->Fit(fcontWeightSphBP[histType], "RES+");
  TFitResultPtr fitResultCylBP = fRegenerationWeights[histType]->Fit(fcontWeightCylBP[histType], "RES+");

  TFitResult *resDC    = fitResultDC.Get();
  TFitResult *resSphBP = fitResultSphBP.Get();
  TFitResult *resCylBP = fitResultCylBP.Get();

  if (resDC != nullptr && resDC->NPar() > 0)
  {
    const double *p = resDC->GetParams();
    const auto &e = resDC->Errors();
    std::cout << "Fit successful for DC Wall bounds, histogram type: " << _HistTypeToString(histType) << std::endl;
    if (p && !e.empty())
      std::cout << "Fit parameters for DC Wall: b = (" << p[0] << " +- " << e[0] << "), a = (" << p[1] << " +- " << e[1] << ")" << std::endl << std::endl;
  }
  else
  {
    std::cerr << "Fit failed for DC Wall bounds, histogram type: " << _HistTypeToString(histType) << std::endl;
  }

  if (resSphBP != nullptr && resSphBP->NPar() > 0)
  {
    const double *p = resSphBP->GetParams();
    const auto &e = resSphBP->Errors();
    std::cout << "Fit successful for Spherical BP bounds, histogram type: " << _HistTypeToString(histType) << std::endl;
    if (p && !e.empty())
      std::cout << "Fit parameters for Spherical BP: b = (" << p[0] << " +- " << e[0] << "), a = (" << p[1] << " +- " << e[1] << ")" << std::endl;
  }
  else
  {
    std::cerr << "Fit failed for Spherical BP bounds, histogram type: " << _HistTypeToString(histType) << std::endl;
  }

  if (resCylBP != nullptr && resCylBP->NPar() > 0)
  {
    const double *p = resCylBP->GetParams();
    const auto &e = resCylBP->Errors();
    std::cout << "Fit successful for Cylindrical BP bounds, histogram type: " << _HistTypeToString(histType) << std::endl;
    if (p && !e.empty())
      std::cout << "Fit parameters for Cylindrical BP: b = (" << p[0] << " +- " << e[0] << "), a = (" << p[1] << " +- " << e[1] << ")" << std::endl;
  }
  else
  {
    std::cerr << "Fit failed for Cylindrical BP bounds, histogram type: " << _HistTypeToString(histType) << std::endl;
  }
}

double RegenerationFractionFit::GetContinuousRegenerationWeight(double dt, double radius)
{
  HistogramType histType = _checkRegenerationLimitsTypesDeltaT(dt);

  if (fcontWeightDC.find(histType) == fcontWeightDC.end() || fcontWeightDC.at(histType) == nullptr)
  {
    if (fRegenerationWeights.at(histType)->GetEntries() > 0)
      _FitContinuousWeightFunction(histType);
    else
      return 1.0; // Default weight if the histogram type is not found and cannot be fitted
  }

  TF1 *weightFunction = nullptr;

  if (radius >= fConfig.DCWallBottomBoundCharged && radius <= fConfig.DCWallTopBoundCharged)
  {
    weightFunction = fcontWeightDC[histType];
  }
  else if (radius >= fConfig.sphericalBPBottomBoundCharged && radius <= fConfig.sphericalBPTopBoundCharged)
  {
    weightFunction = fcontWeightSphBP[histType];
  }
  else if (radius >= fConfig.cylindricalBPBottomBoundCharged && radius <= fConfig.cylindricalBPTopBoundCharged)
  {
    weightFunction = fcontWeightCylBP[histType];
  }
  else if (radius >= fConfig.DCWallBottomBoundNeutral && radius <= fConfig.DCWallTopBoundNeutral)
  {
    weightFunction = fcontWeightDC[histType];
  }
  else if (radius >= fConfig.sphericalBPBottomBoundNeutral && radius <= fConfig.sphericalBPTopBoundNeutral)
  {
    weightFunction = fcontWeightSphBP[histType];
  }
  else if (radius >= fConfig.cylindricalBPBottomBoundNeutral && radius <= fConfig.cylindricalBPTopBoundNeutral)
  {
    weightFunction = fcontWeightCylBP[histType];
  }

  if (weightFunction == nullptr)
  {
    return 1.0; // Default weight if the radius is outside the defined bounds
  }

  double weight = weightFunction->Eval(radius);
  return weight;
}

double RegenerationFractionFit::GetContinuousRegenerationWeight(double dt, std::map<HistogramType, double> radius)
{
  std::pair<TF1*, HistogramType> fitFuncPair = _GetFitFuncByDt(dt);

  TF1 *weightFunction = nullptr;//fitFuncPair.first;
  HistogramType histTypeFromDt = fitFuncPair.second;

  if (radius[histTypeFromDt] >= fConfig.DCWallBottomBoundCharged && radius[histTypeFromDt] <= fConfig.DCWallTopBoundCharged)
  {
    weightFunction = fcontWeightDC[histTypeFromDt];
  }
  else if (radius[histTypeFromDt] >= fConfig.sphericalBPBottomBoundCharged && radius[histTypeFromDt] <= fConfig.sphericalBPTopBoundCharged)
  {
    weightFunction = fcontWeightSphBP[histTypeFromDt];
  }
  else if (radius[histTypeFromDt] >= fConfig.cylindricalBPBottomBoundCharged && radius[histTypeFromDt] <= fConfig.cylindricalBPTopBoundCharged)
  {
    weightFunction = fcontWeightCylBP[histTypeFromDt];
  }
  else if (radius[histTypeFromDt] >= fConfig.DCWallBottomBoundNeutral && radius[histTypeFromDt] <= fConfig.DCWallTopBoundNeutral)
  {
    weightFunction = fcontWeightDC[histTypeFromDt];
  }
  else if (radius[histTypeFromDt] >= fConfig.sphericalBPBottomBoundNeutral && radius[histTypeFromDt] <= fConfig.sphericalBPTopBoundNeutral)
  {
    weightFunction = fcontWeightSphBP[histTypeFromDt];
  }
  else if (radius[histTypeFromDt] >= fConfig.cylindricalBPBottomBoundNeutral && radius[histTypeFromDt] <= fConfig.cylindricalBPTopBoundNeutral)
  {
    weightFunction = fcontWeightCylBP[histTypeFromDt];
  }

  if (weightFunction == nullptr)
  {
    return 1.0; // Default weight if the radius is outside the defined bounds
  }

  double weight = weightFunction->Eval(radius[histTypeFromDt]);
  return weight;
}