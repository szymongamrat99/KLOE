#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cmath>

#include "TMath.h"
#include "nlohmann/json.hpp"

using json = nlohmann::json;

struct Parameter
{
  Double_t initial_value;
  Double_t step;
  Double_t limit_lower;
  Double_t limit_upper;
  std::vector<std::string> channels;
  bool enabled;
  bool fixed = false; // Indicates if the parameter is fixed (not to be fitted)
  std::array<Double_t, 2> limits() const
  {
    const Double_t scale = std::abs(initial_value);
    return {initial_value - limit_lower * scale, initial_value + limit_upper * scale};
  };
  std::vector<std::array<Double_t, 2>> parameter_ranges;
};

struct DeltaTConfig
{
  std::array<Double_t, 2> xRange;
  std::array<Double_t, 2> xRangeDisplay;
  Double_t resolution;
};

struct CutConfig
{
  bool ellipseEnabled;
  Double_t ellipseU0, ellipseV0, ellipseSigmaU, ellipseSigmaV, ellipseRho, ellipseNSigma;
  bool trackAngleEnabled;
  Double_t vertexDistanceXY, vertexDistanceZ, maxTrackCosine;
  bool closeVertexEnabled;
  Double_t neutralRadius, chargedRadius, neutralZ, chargedZ, minCombinedRadius;

  bool outsideEllipse(Double_t u, Double_t v) const
  {
    const Double_t du = (u - ellipseU0) / ellipseSigmaU;
    const Double_t dv = (v - ellipseV0) / ellipseSigmaV;
    const Double_t distance2 = (du * du + dv * dv - 2.0 * ellipseRho * du * dv) /
                               (1.0 - ellipseRho * ellipseRho);
    return std::sqrt(std::max(0.0, distance2)) > ellipseNSigma;
  }
};

#ifndef REGEN_KIN_RANGE_DEFINED
#define REGEN_KIN_RANGE_DEFINED
struct RegenKinRange
{
  std::string variable;   // "R_Charged" | "Rho_Charged" | "R_Neutral" | "Rho_Neutral"
  std::string histName;   // nazwa TH2 w pliku ROOT (klucz osi: deltaT x var)
  Double_t rangeMin;      // dolna granica zmiennej kinematycznej
  Double_t rangeMax;      // górna granica zmiennej kinematycznej
};
#endif

struct RegenShapeCorrectionConfig
{
  std::string correctionFile;
  std::string weightFunctionsFile;
  std::vector<RegenKinRange> corrections;
};

struct EventWeightConfig
{
  bool applySignal = true;
  bool applySemileptonic = true;
  bool applyThreePi0 = true;
  bool applyRegeneration = true;
  std::string signalCorrectionFile;
  std::string signalCorrectionHistogram;
  Double_t threePi0 = 0.9375;
  Double_t semileptonic = 0.9272;
  Double_t signalMaxChargedRadius = 20.0;
};

struct RegenerationExclusionConfig
{
  bool enabled = false;
  Double_t beamPipeBoundary = 5.0;
  Double_t meanCharged = 10.4941;
  Double_t sigmaCharged = 0.957544;
  Double_t meanNeutral = 10.3769;
  Double_t sigmaNeutral = 1.23898;
  Double_t sigmaMultiplier = 1.5;
};

struct FitConfig
{
  std::map<std::string, Parameter> parameters;
  std::map<std::string, std::vector<std::string>> channel_to_param;
  std::string minimizer_type;
  std::string minimizer_algorithm;
  Int_t max_function_calls;
  Double_t tolerance;
  Int_t print_level;
  Int_t strategy;
  DeltaTConfig deltaTConfig;
  CutConfig cuts;
  RegenerationExclusionConfig regenerationExclusion;
  EventWeightConfig eventWeights;
  std::array<Double_t, 3> regenerationSplit = {-30.0, 0.0, 30.0};
  RegenShapeCorrectionConfig regenShapeCorrection;

  Int_t getNumOfEnabledParameters() const
  {
    Int_t count = 0;
    for (const auto &entry : parameters)
    {
      if (entry.second.enabled)
        ++count;
    }
    return count;
  }
};

class FitSetter
{
private:
  json _config;

  Bool_t _range_overlaps(const std::array<Double_t, 2> &r1, const std::array<Double_t, 2> &r2) const
  {
    return r1[0] < r2[1] && r2[0] < r1[1];
  };

public:
  const json &getRawConfig() const { return _config; }

  FitSetter(std::string config_path)
  {
    try
    {
      std::ifstream configFile(config_path);

      if (!configFile.is_open())
      {
        throw std::runtime_error("ERROR: Could not open config file: " + config_path);
      }

      try
      {
        _config = json::parse(configFile);
      }
      catch (const json::parse_error &e)
      {
        throw std::runtime_error("ERROR: Failed to parse config file: " + std::string(e.what()));
      }
    }
    catch (const std::exception &e)
    {
      std::cerr << e.what() << '\n';
      throw; // Rethrow the exception after logging
    }
  }

  FitConfig getFitConfig() const
  {
    FitConfig fitConfig;

    // 1. Extract fit settings
    try
    {
      const auto &fs = _config.at("fitSettings");
      fitConfig.minimizer_type = fs.at("minimizerType").get<std::string>();
      fitConfig.minimizer_algorithm = fs.at("minimizerAlgorithm").get<std::string>();
      fitConfig.max_function_calls = fs.at("maxFunctionCalls").get<Int_t>();
      fitConfig.tolerance = fs.at("tolerance").get<Double_t>();
      fitConfig.print_level = fs.at("printLevel").get<Int_t>();
      fitConfig.strategy = fs.at("strategy").get<Int_t>();
    }
    catch (const json::out_of_range &e)
    {
      throw std::runtime_error("ERROR: Missing fit setting in config file: " + std::string(e.what()));
    }

    try
    {
      const auto &exclusion = _config.at("regenerationExclusion");
      fitConfig.regenerationExclusion.enabled = exclusion.at("enabled").get<bool>();
      fitConfig.regenerationExclusion.beamPipeBoundary = exclusion.at("beamPipeBoundary").get<Double_t>();
      fitConfig.regenerationExclusion.meanCharged = exclusion.at("spherical").at("meanCharged").get<Double_t>();
      fitConfig.regenerationExclusion.sigmaCharged = exclusion.at("spherical").at("sigmaCharged").get<Double_t>();
      fitConfig.regenerationExclusion.meanNeutral = exclusion.at("spherical").at("meanNeutral").get<Double_t>();
      fitConfig.regenerationExclusion.sigmaNeutral = exclusion.at("spherical").at("sigmaNeutral").get<Double_t>();
      fitConfig.regenerationExclusion.sigmaMultiplier = exclusion.at("sigmaMultiplier").get<Double_t>();
    }
    catch (const json::out_of_range &e)
    {
      throw std::runtime_error("ERROR: Missing regeneration exclusion setting in config file: " + std::string(e.what()));
    }

    try
    {
      const auto &weights = _config.at("eventWeights");
      const auto &applyWeights = weights.at("apply");
      fitConfig.eventWeights.applySignal = applyWeights.at("signal").get<bool>();
      fitConfig.eventWeights.applySemileptonic = applyWeights.at("semileptonic").get<bool>();
      fitConfig.eventWeights.applyThreePi0 = applyWeights.at("threePi0").get<bool>();
      fitConfig.eventWeights.applyRegeneration = applyWeights.at("regeneration").get<bool>();
      fitConfig.eventWeights.signalCorrectionFile = weights.at("signalCorrectionFile").get<std::string>();
      fitConfig.eventWeights.signalCorrectionHistogram = weights.at("signalCorrectionHistogram").get<std::string>();
      fitConfig.eventWeights.threePi0 = weights.at("threePi0").get<Double_t>();
      fitConfig.eventWeights.semileptonic = weights.at("semileptonic").get<Double_t>();
      fitConfig.eventWeights.signalMaxChargedRadius = weights.at("signalMaxChargedRadius").get<Double_t>();
    }
    catch (const json::out_of_range &e)
    {
      throw std::runtime_error("ERROR: Missing event weight setting in config file: " + std::string(e.what()));
    }

    const auto &cuts = _config.at("cuts");
    const auto &ellipse = cuts.at("ellipse");
    fitConfig.cuts.ellipseEnabled = ellipse.at("enabled").get<bool>();
    fitConfig.cuts.ellipseU0 = ellipse.at("u0").get<Double_t>();
    fitConfig.cuts.ellipseV0 = ellipse.at("v0").get<Double_t>();
    fitConfig.cuts.ellipseSigmaU = ellipse.at("sigmaU").get<Double_t>();
    fitConfig.cuts.ellipseSigmaV = ellipse.at("sigmaV").get<Double_t>();
    fitConfig.cuts.ellipseRho = ellipse.at("rho").get<Double_t>();
    fitConfig.cuts.ellipseNSigma = ellipse.at("nSigma").get<Double_t>();

    const auto &track = cuts.at("trackAngle");
    fitConfig.cuts.trackAngleEnabled = track.at("enabled").get<bool>();
    fitConfig.cuts.vertexDistanceXY = track.at("vertexDistanceXY").get<Double_t>();
    fitConfig.cuts.vertexDistanceZ = track.at("vertexDistanceZ").get<Double_t>();
    fitConfig.cuts.maxTrackCosine = track.at("maxCosine").get<Double_t>();

    const auto &close = cuts.at("closeVertex");
    fitConfig.cuts.closeVertexEnabled = close.at("enabled").get<bool>();
    fitConfig.cuts.neutralRadius = close.at("neutralRadius").get<Double_t>();
    fitConfig.cuts.chargedRadius = close.at("chargedRadius").get<Double_t>();
    fitConfig.cuts.neutralZ = close.at("neutralZ").get<Double_t>();
    fitConfig.cuts.chargedZ = close.at("chargedZ").get<Double_t>();
    fitConfig.cuts.minCombinedRadius = close.at("minCombinedRadius").get<Double_t>();
    if (!std::isfinite(fitConfig.cuts.ellipseSigmaU) || fitConfig.cuts.ellipseSigmaU <= 0 ||
        !std::isfinite(fitConfig.cuts.ellipseSigmaV) || fitConfig.cuts.ellipseSigmaV <= 0 ||
        !std::isfinite(fitConfig.cuts.ellipseU0) || !std::isfinite(fitConfig.cuts.ellipseV0) ||
        !std::isfinite(fitConfig.cuts.ellipseRho) || std::abs(fitConfig.cuts.ellipseRho) >= 1 ||
        !std::isfinite(fitConfig.cuts.ellipseNSigma) || fitConfig.cuts.ellipseNSigma < 0)
      throw std::runtime_error("ERROR: Invalid ellipse cut parameters");
    if (!std::isfinite(fitConfig.cuts.vertexDistanceXY) || fitConfig.cuts.vertexDistanceXY < 0 ||
        !std::isfinite(fitConfig.cuts.vertexDistanceZ) || fitConfig.cuts.vertexDistanceZ < 0 ||
        !std::isfinite(fitConfig.cuts.maxTrackCosine) || fitConfig.cuts.maxTrackCosine < 0 ||
        fitConfig.cuts.maxTrackCosine > 1 ||
        !std::isfinite(fitConfig.cuts.neutralRadius) || fitConfig.cuts.neutralRadius < 0 ||
        !std::isfinite(fitConfig.cuts.chargedRadius) || fitConfig.cuts.chargedRadius < 0 ||
        !std::isfinite(fitConfig.cuts.neutralZ) || fitConfig.cuts.neutralZ < 0 ||
        !std::isfinite(fitConfig.cuts.chargedZ) || fitConfig.cuts.chargedZ < 0 ||
        !std::isfinite(fitConfig.cuts.minCombinedRadius) || fitConfig.cuts.minCombinedRadius < 0)
      throw std::runtime_error("ERROR: Invalid track angle or close vertex cut parameters");
    if (!std::isfinite(fitConfig.eventWeights.threePi0) || fitConfig.eventWeights.threePi0 < 0 ||
        !std::isfinite(fitConfig.eventWeights.semileptonic) || fitConfig.eventWeights.semileptonic < 0 ||
        !std::isfinite(fitConfig.eventWeights.signalMaxChargedRadius) ||
        fitConfig.eventWeights.signalMaxChargedRadius < 0)
      throw std::runtime_error("ERROR: Invalid event weight settings");
    if (!std::isfinite(fitConfig.regenerationExclusion.beamPipeBoundary) ||
        fitConfig.regenerationExclusion.beamPipeBoundary < 0 ||
        !std::isfinite(fitConfig.regenerationExclusion.meanCharged) ||
        !std::isfinite(fitConfig.regenerationExclusion.meanNeutral) ||
        !std::isfinite(fitConfig.regenerationExclusion.sigmaCharged) ||
        fitConfig.regenerationExclusion.sigmaCharged <= 0 ||
        !std::isfinite(fitConfig.regenerationExclusion.sigmaNeutral) ||
        fitConfig.regenerationExclusion.sigmaNeutral <= 0 ||
        !std::isfinite(fitConfig.regenerationExclusion.sigmaMultiplier) ||
        fitConfig.regenerationExclusion.sigmaMultiplier < 0)
      throw std::runtime_error("ERROR: Invalid regeneration exclusion settings");

    // Regeneration split boundaries (used in single A_regen mode)
    if (_config.contains("regenerationSplit"))
    {
      fitConfig.regenerationSplit = _config.at("regenerationSplit").get<std::array<Double_t, 3>>();
    }

    // 2. Extract DeltaTConfig (DOPISANE)
    try
    {
      const auto &dtc = _config.at("deltaT");
      fitConfig.deltaTConfig.xRange = dtc.at("rangeX").get<std::array<Double_t, 2>>();
      fitConfig.deltaTConfig.xRangeDisplay = dtc.at("xRangeDisplay").get<std::array<Double_t, 2>>();
      fitConfig.deltaTConfig.resolution = dtc.at("resolution").get<Double_t>();
    }
    catch (const json::out_of_range &e)
    {
      throw std::runtime_error("ERROR: Missing deltaT setting: " + std::string(e.what()));
    }
    if (!std::isfinite(fitConfig.deltaTConfig.resolution) || fitConfig.deltaTConfig.resolution <= 0 ||
        !std::isfinite(fitConfig.deltaTConfig.xRange[0]) || !std::isfinite(fitConfig.deltaTConfig.xRange[1]) ||
        fitConfig.deltaTConfig.xRange[1] <= fitConfig.deltaTConfig.xRange[0] ||
        fitConfig.deltaTConfig.resolution > fitConfig.deltaTConfig.xRange[1] - fitConfig.deltaTConfig.xRange[0])
      throw std::runtime_error("ERROR: Invalid deltaT range or resolution");

    // 3. Extract parameter settings
    const auto &enabled_params = _config.at("parameterSettings");

    for (auto it = enabled_params.begin(); it != enabled_params.end(); ++it)
    {
      const std::string &paramName = it.key();
      const auto &val = it.value();

      Parameter p;
      p.enabled = val.at("enabled").get<bool>();
      p.fixed = val.at("fixed").get<bool>();
      p.initial_value = val.at("initialValue").get<Double_t>();
      p.step = val.at("step").get<Double_t>();
      p.limit_lower = val.at("limitLower").get<Double_t>();
      p.limit_upper = val.at("limitUpper").get<Double_t>();

      if (val.contains("channels"))
      {
        p.channels = val.at("channels").get<std::vector<std::string>>();
      }

      if (val.contains("parameter_ranges"))
      {
        p.parameter_ranges = val.at("parameter_ranges").get<std::vector<std::array<Double_t, 2>>>();
      }

      fitConfig.parameters[paramName] = p;

      if (p.enabled)
      {
        for (const auto &channel : p.channels)
        {
          fitConfig.channel_to_param[channel].push_back(paramName);
        }
      }
    }

    // 4. Walidacja rozłączności zakresów wewnątrz kanałów
    for (auto const &entry : fitConfig.channel_to_param)
    {
      const std::string &channelName = entry.first;
      const std::vector<std::string> &paramsInChannel = entry.second;

      std::vector<std::pair<std::array<Double_t, 2>, std::string>> allRangesInChannel;

      for (const auto &pName : paramsInChannel)
      {
        const Parameter &p = fitConfig.parameters.at(pName);
        if (!p.enabled)
          continue;

        for (const auto &range : p.parameter_ranges)
        {
          // OPCJONALNIE: Sprawdzenie czy zakres parametru mieści się w globalnym xRange
          if (range[0] < fitConfig.deltaTConfig.xRange[0] || range[1] > fitConfig.deltaTConfig.xRange[1])
          {
            std::cerr << "WARNING: Range (" << range[0] << "," << range[1] << ") for param ["
                      << pName << "] is outside global xRange.\n";
          }
          allRangesInChannel.push_back(std::make_pair(range, pName));
        }
      }

      for (size_t i = 0; i < allRangesInChannel.size(); ++i)
      {
        for (size_t j = i + 1; j < allRangesInChannel.size(); ++j)
        {
          const auto &r1 = allRangesInChannel[i].first;
          const auto &r2 = allRangesInChannel[j].first;

          if (_range_overlaps(r1, r2))
          {
            std::stringstream ss;
            ss << "FATAL: Overlap detected in channel [" << channelName << "].\n"
               << "Range (" << r1[0] << ", " << r1[1] << ") from [" << allRangesInChannel[i].second << "]\n"
               << "overlaps with (" << r2[0] << ", " << r2[1] << ") from [" << allRangesInChannel[j].second << "].";
            throw std::runtime_error(ss.str());
          }
        }
      }
    }

    if (_config.contains("regenerationShapeCorrection"))
    {
      const auto &rsc = _config.at("regenerationShapeCorrection");
      fitConfig.regenShapeCorrection.correctionFile = rsc.value("correctionFile", "");
      fitConfig.regenShapeCorrection.weightFunctionsFile = rsc.value("weightFunctionsFile", "regeneration_analysis_results.root");
      if (rsc.contains("corrections"))
      {
        for (const auto &c : rsc.at("corrections"))
        {
          RegenKinRange entry;
          entry.histName = c.at("histName").get<std::string>();
          entry.variable = c.at("variable").get<std::string>();
          entry.rangeMin = c.at("rangeMin").get<Double_t>();
          entry.rangeMax = c.at("rangeMax").get<Double_t>();
          fitConfig.regenShapeCorrection.corrections.push_back(entry);
        }
      }
    }

    return fitConfig;
  }
};