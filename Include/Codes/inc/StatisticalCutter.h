#pragma once
#include <string>
#include <vector>
#include <map>
#include <functional>
#include <stdexcept>
#include <json.hpp>
#include <fstream>

struct Cut {
    std::string cutId;
    std::string cutType;
    std::string cutCondition; // "<=", ">=", "==", etc.
    double cutValue;
    double centralValue;
    std::function<double()> valueGetter; // Getter do zmiennej
};

class StatisticalCutter {
public:
    StatisticalCutter(const std::string& jsonPath, int signalMctruth);

    void RegisterVariableGetter(const std::string& varName, std::function<double()> getter);

    bool PassCut(size_t cutIndex);
    bool PassAllCuts();

    void UpdateStats(int mctruth);

    double GetEfficiency(size_t cutIndex) const;
    double GetPurity(size_t cutIndex) const;
    double GetSignalToBackground(size_t cutIndex) const;

    size_t GetSurvivedSignal(size_t cutIndex) const;
    size_t GetSurvivedBackground(size_t cutIndex) const;

private:
    void LoadCuts(const std::string& jsonPath);
    bool EvaluateCondition(double value, const Cut& cut) const;

    std::vector<Cut> cuts_;
    int signalMctruth_;
    std::map<std::string, std::function<double()>> variableGetters_;

    std::vector<size_t> survivedSignal_;
    std::vector<size_t> survivedBackground_;
    size_t totalSignal_ = 0;
    size_t totalBackground_ = 0;
};