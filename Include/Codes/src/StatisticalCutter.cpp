#include "StatisticalCutter.h"

#include <iostream>

using json = nlohmann::json;

StatisticalCutter::StatisticalCutter(const std::string &jsonPath, int signalMctruth)
    : signalMctruth_(signalMctruth)
{
  LoadCuts(jsonPath);
  survivedSignal_.resize(cuts_.size(), 0);
  survivedBackground_.resize(cuts_.size(), 0);
}

void StatisticalCutter::RegisterVariableGetter(const std::string &varName, std::function<double()> getter)
{
  variableGetters_[varName] = getter;
  // Przypisz getter do odpowiednich cięć
  for (auto &cut : cuts_)
  {
    if (cut.cutId == varName)
    {
      cut.valueGetter = getter;
    }
  }
}

void StatisticalCutter::LoadCuts(const std::string &jsonPath)
{
  std::ifstream file(jsonPath);
  if (!file.is_open())
    throw std::runtime_error("Cannot open cut file: " + jsonPath);

  nlohmann::json j;

  try
  {
    j = nlohmann::json::parse(file);
  }
  catch (const std::exception &e)
  {
    throw std::runtime_error("JSON parse error: " + std::string(e.what()));
  }

  for (const auto &cutj : j["listOfCuts"])
  {
    Cut cut;
    cut.cutId = cutj["cutId"];
    cut.cutType = cutj["cutType"];
    cut.cutCondition = cutj["cutCondition"];
    cut.cutValue = cutj["cutValue"];
    cut.centralValue = cutj.value("centralValue", 0.0);
    cuts_.push_back(cut);
  }
}

bool StatisticalCutter::EvaluateCondition(double value, const Cut &cut) const
{
  if (cut.cutCondition == "<=")
    return value <= cut.centralValue + cut.cutValue;
  if (cut.cutCondition == ">=")
    return value >= cut.centralValue - cut.cutValue;
  if (cut.cutCondition == "<")
    return value < cut.centralValue + cut.cutValue;
  if (cut.cutCondition == ">")
    return value > cut.centralValue - cut.cutValue;
  if (cut.cutCondition == "==")
    return value == cut.centralValue;
  // Możesz dodać więcej warunków
  throw std::runtime_error("Unknown cut condition: " + cut.cutCondition);
}

bool StatisticalCutter::PassCut(size_t cutIndex)
{
  if (cutIndex >= cuts_.size())
    throw std::out_of_range("Cut index out of range");
  const auto &cut = cuts_[cutIndex];
  if (!cut.valueGetter)
    throw std::runtime_error("No getter registered for variable: " + cut.cutId);
  double value = cut.valueGetter();
  return EvaluateCondition(value, cut);
}

bool StatisticalCutter::PassAllCuts()
{
  for (size_t i = 0; i < cuts_.size(); ++i)
  {
    if (!PassCut(i))
      return false;
  }
  return true;
}

void StatisticalCutter::UpdateStats(int mctruth)
{
  // Zliczaj sygnał/tło przed cięciami
  if (mctruth == signalMctruth_)
    totalSignal_++;
  else
    totalBackground_++;

  // Sprawdź przeżycie po kolejnych cięciach
  bool survived = true;
  for (size_t i = 0; i < cuts_.size(); ++i)
  {
    if (survived && PassCut(i))
    {
      if (mctruth == signalMctruth_)
        survivedSignal_[i]++;
      else
        survivedBackground_[i]++;
    }
    else
    {
      survived = false;
    }
  }
}

double StatisticalCutter::GetEfficiency(size_t cutIndex) const
{
  if (totalSignal_ == 0)
    return 0.0;
  return static_cast<double>(survivedSignal_[cutIndex]) / totalSignal_;
}

double StatisticalCutter::GetPurity(size_t cutIndex) const
{
  size_t sig = survivedSignal_[cutIndex];
  size_t bkg = survivedBackground_[cutIndex];
  if (sig + bkg == 0)
    return 0.0;
  return static_cast<double>(sig) / (sig + bkg);
}

double StatisticalCutter::GetSignalToBackground(size_t cutIndex) const
{
  size_t sig = survivedSignal_[cutIndex];
  size_t bkg = survivedBackground_[cutIndex];
  if (bkg == 0)
    return sig > 0 ? 1e9 : 0.0;
  return static_cast<double>(sig) / bkg;
}

size_t StatisticalCutter::GetSurvivedSignal(size_t cutIndex) const
{
  return survivedSignal_[cutIndex];
}

size_t StatisticalCutter::GetSurvivedBackground(size_t cutIndex) const
{
  return survivedBackground_[cutIndex];
}