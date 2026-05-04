#include <iostream>
#include <fstream>
#include <sstream>

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
    return {initial_value - limit_lower * initial_value, initial_value + limit_upper * initial_value};
  };
  std::vector<std::array<Double_t, 2>> parameter_ranges;
};

struct DeltaTConfig
{
  std::array<Double_t, 2> xRange;
  std::array<Double_t, 2> xRangeDisplay;
  Double_t resolution;
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
  Bool_t regenerationExclusionFlag;

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

    // Regeneration exclusion flag    try
    try 
    {
      fitConfig.regenerationExclusionFlag = _config.at("regenerationExclusion").at("enabled").get<Bool_t>();
    }
    catch (const json::out_of_range &e)
    {
      throw std::runtime_error("ERROR: Missing regenerationExclusionFlag in config file: " + std::string(e.what()));
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
        if (!p.enabled) continue;

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

    return fitConfig;
  }
};