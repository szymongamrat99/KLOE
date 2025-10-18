#pragma once
#include <string>
#include <vector>
#include <map>
#include <functional>
#include <stdexcept>
#include <json.hpp>
#include <fstream>
#include <algorithm>
// #include "../../klspm00.hpp"  // dla KLOE::HypothesisCode
#include <HypothesisCodes.h>

struct Cut {
    int order = 0;
    std::string cutId;
    std::string cutType;
    std::string cutCondition;
    double cutValue;
    double centralValue;
    bool centralValueDynamic = false;
    std::function<double()> valueGetter;
    std::function<double()> centralValueGetter;
};

class StatisticalCutter {
public:
    // Konstruktor z jednym plikiem JSON (dla kompatybilności)
    StatisticalCutter(const std::string& jsonPath, int signalMctruth, KLOE::HypothesisCode hypoCode);

    // Konstruktor z dwoma plikami (Utils::properties i cuts)
    StatisticalCutter(const std::string& propertiesPath, const std::string& cutsPath, KLOE::HypothesisCode hypoCode);

    void RegisterVariableGetter(const std::string& varName, std::function<double()> getter);
    void RegisterCentralValueGetter(const std::string& cutId, std::function<double()> getter);

    bool PassCut(size_t cutIndex);
    bool PassAllCuts();

    void UpdateStats(int mctruth);

    double GetEfficiency(size_t cutIndex) const;
    double GetEfficiencyError(size_t cutIndex) const;
    
    double GetPurity(size_t cutIndex) const;
    double GetPurityError(size_t cutIndex) const;
    
    double GetSignalToBackground(size_t cutIndex) const;
    double GetSignalToBackgroundError(size_t cutIndex) const;

    size_t GetSurvivedSignal(size_t cutIndex) const;
    size_t GetSurvivedBackground(size_t cutIndex) const;

    // Dostęp do cięć (np. po nazwach)
    const std::vector<Cut>& GetCuts() const { return cuts_; }

private:
    void LoadCuts(const nlohmann::json& j);
    void LoadCuts(const std::string& jsonPath);
    void LoadCutsFromFiles(const std::string& propertiesPath, const std::string& cutsPath);
    bool EvaluateCondition(double value, const Cut& cut) const;

    KLOE::HypothesisCode hypoCode_;

    std::vector<Cut> cuts_;
    int signalMctruth_;
    std::map<std::string, std::function<double()>> variableGetters_;

    std::vector<size_t> survivedSignal_;
    std::vector<size_t> survivedBackground_;
    size_t totalSignal_ = 0;
    size_t totalBackground_ = 0;
};