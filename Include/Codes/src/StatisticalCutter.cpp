#include "StatisticalCutter.h"
#include <iostream>

using json = nlohmann::json;

StatisticalCutter::StatisticalCutter(const std::string& jsonPath, int signalMctruth, KLOE::HypothesisCode hypoCode)
    : signalMctruth_(signalMctruth), hypoCode_(hypoCode)
{
    LoadCuts(jsonPath);
    survivedSignal_.resize(cuts_.size(), 0);
    survivedBackground_.resize(cuts_.size(), 0);
}

StatisticalCutter::StatisticalCutter(const std::string& propertiesPath, const std::string& cutsPath, KLOE::HypothesisCode hypoCode)
    : hypoCode_(hypoCode)
{
    LoadCutsFromFiles(propertiesPath, cutsPath);
    survivedSignal_.resize(cuts_.size(), 0);
    survivedBackground_.resize(cuts_.size(), 0);
}

void StatisticalCutter::RegisterVariableGetter(const std::string& varName, std::function<double()> getter) {
    variableGetters_[varName] = getter;
    for (auto& cut : cuts_) {
        if (cut.cutId == varName) {
            cut.valueGetter = getter;
        }
    }
}

void StatisticalCutter::RegisterCentralValueGetter(const std::string& cutId, std::function<double()> getter) {
    for (auto& cut : cuts_) {
        if (cut.cutId == cutId) {
            cut.centralValueGetter = getter;
            cut.centralValueDynamic = true;
        }
    }
}

void StatisticalCutter::LoadCuts(const std::string& jsonPath) {
    std::ifstream file(jsonPath);
    if (!file.is_open())
        throw std::runtime_error("Cannot open cut file: " + jsonPath);

    json j;
    try {
        file >> j;
    } catch (const std::exception& e) {
        throw std::runtime_error("JSON parse error: " + std::string(e.what()));
    }
    LoadCuts(j);
}

void StatisticalCutter::LoadCuts(const json& j) {
    cuts_.clear();
    
    // Wybierz właściwą listę cięć na podstawie hypoCode_
    std::string cutListKey;
    switch(hypoCode_) {
        case KLOE::HypothesisCode::SIGNAL:
            cutListKey = "SIGNAL";
            break;
        case KLOE::HypothesisCode::FOUR_PI:
            cutListKey = "FOUR_PI";
            break;
        case KLOE::HypothesisCode::OMEGAPI:
            cutListKey = "OMEGA";
            break;
        default:
            std::cout << "Warning: Using default cuts for unknown hypothesis code" << std::endl;
            return;
    }

    // Sprawdź czy mamy odpowiednią sekcję w listOfCuts
    if (!j.contains("listOfCuts") || !j["listOfCuts"].contains(cutListKey)) {
        std::cout << "Warning: No cuts found for hypothesis " << cutListKey << std::endl;
        return;
    }

    for (const auto& cutj : j["listOfCuts"][cutListKey]) {
        Cut cut;
        cut.order = cutj.value("order", 0);
        cut.cutId = cutj["cutId"];
        cut.cutType = cutj["cutType"];
        cut.cutCondition = cutj["cutCondition"];
        cut.cutValue = cutj["cutValue"];
        cut.centralValue = cutj.value("centralValue", 0.0);
        cut.centralValueDynamic = false;
        cuts_.push_back(cut);
    }
    
    // Sortuj po order malejąco (zmień na rosnąco jeśli wolisz)
    std::sort(cuts_.begin(), cuts_.end(), [](const Cut& a, const Cut& b) {
        return a.order < b.order;
    });
}

void StatisticalCutter::LoadCutsFromFiles(const std::string& propertiesPath, const std::string& cutsPath) {
    json props, cuts;
    {
        std::ifstream f(propertiesPath);
        if (!f.is_open()) throw std::runtime_error("Cannot open properties file: " + propertiesPath);
        f >> props;
    }
    {
        std::ifstream f(cutsPath);
        if (!f.is_open()) throw std::runtime_error("Cannot open cuts file: " + cutsPath);
        f >> cuts;
    }
    int analysisCode = 0;
    if (props["flags"]["analysisCode"].is_string())
        analysisCode = std::stoi(props["flags"]["analysisCode"].get<std::string>());
    else
        analysisCode = props["flags"]["analysisCode"];
    std::string listName = "listOfCuts";
    std::string altList = "cutsFor" + std::to_string(analysisCode);
    if (cuts.contains(altList))
        listName = altList;
    LoadCuts(cuts[listName]);
}

bool StatisticalCutter::EvaluateCondition(double value, const Cut& cut) const {
    double central = cut.centralValueDynamic && cut.centralValueGetter
        ? cut.centralValueGetter()
        : cut.centralValue;

    // One-sided cuts ("O")
    if (cut.cutType == "O") {
        if (cut.cutCondition == "<=") return value <= cut.cutValue;
        if (cut.cutCondition == ">=") return value >= cut.cutValue;
        if (cut.cutCondition == "<")  return value < cut.cutValue;
        if (cut.cutCondition == ">")  return value > cut.cutValue;
    }
    // Two-sided cuts ("T")
    else if (cut.cutType == "T") {
        double deviation = std::abs(value - central);
        if (cut.cutCondition == "<=") return deviation <= cut.cutValue;
        if (cut.cutCondition == ">=") return deviation >= cut.cutValue;
        if (cut.cutCondition == "<")  return deviation < cut.cutValue;
        if (cut.cutCondition == ">")  return deviation > cut.cutValue;
        if (cut.cutCondition == "==") return deviation == cut.cutValue;
    }
    
    throw std::runtime_error("Invalid cut type '" + cut.cutType + "' or condition '" + cut.cutCondition + "'");
}

bool StatisticalCutter::PassCut(size_t cutIndex) {
    if (cutIndex >= cuts_.size())
        throw std::out_of_range("Cut index out of range");
    const auto& cut = cuts_[cutIndex];
    if (!cut.valueGetter)
        throw std::runtime_error("No getter registered for variable: " + cut.cutId);
    double value = cut.valueGetter();
    return EvaluateCondition(value, cut);
}

bool StatisticalCutter::PassAllCuts() {
    for (size_t i = 0; i < cuts_.size(); ++i) {
        if (!PassCut(i))
            return false;
    }
    return true;
}

void StatisticalCutter::UpdateStats(int mctruth) {
    if (mctruth == signalMctruth_ || mctruth == 0 || mctruth == -1)
        totalSignal_++;
    else
        totalBackground_++;

    bool survived = true;
    for (size_t i = 0; i < cuts_.size(); ++i) {
        if (survived && PassCut(i)) {
            if (mctruth == signalMctruth_)
                survivedSignal_[i]++;
            else
                survivedBackground_[i]++;
        } else {
            survived = false;
        }
    }
}

double StatisticalCutter::GetEfficiency(size_t cutIndex) const {
    if (totalSignal_ == 0) return 0.0;
    return static_cast<double>(survivedSignal_[cutIndex]) / totalSignal_;
}

double StatisticalCutter::GetPurity(size_t cutIndex) const {
    size_t sig = survivedSignal_[cutIndex];
    size_t bkg = survivedBackground_[cutIndex];
    if (sig + bkg == 0) return 0.0;
    return static_cast<double>(sig) / (sig + bkg);
}

double StatisticalCutter::GetSignalToBackground(size_t cutIndex) const {
    size_t sig = survivedSignal_[cutIndex];
    size_t bkg = survivedBackground_[cutIndex];
    if (bkg == 0) return sig > 0 ? 1e9 : 0.0;
    return static_cast<double>(sig) / bkg;
}

size_t StatisticalCutter::GetSurvivedSignal(size_t cutIndex) const {
    return survivedSignal_[cutIndex];
}

size_t StatisticalCutter::GetSurvivedBackground(size_t cutIndex) const {
    return survivedBackground_[cutIndex];
}

double StatisticalCutter::GetEfficiencyError(size_t cutIndex) const {
    if (totalSignal_ == 0) return 0.0;
    double eff = GetEfficiency(cutIndex);
    // Niepewność z rozkładu dwumianowego
    return std::sqrt(eff * (1.0 - eff) / totalSignal_);
}

double StatisticalCutter::GetPurityError(size_t cutIndex) const {
    size_t sig = survivedSignal_[cutIndex];
    size_t bkg = survivedBackground_[cutIndex];
    double total = sig + bkg;
    
    if (total == 0) return 0.0;
    
    double purity = static_cast<double>(sig) / total;
    // Propagacja niepewności dla stosunku, zakładając rozkład Poissona
    if (sig == 0) return 0.0;  // Unikamy dzielenia przez 0
    
    return purity * std::sqrt(
        (1.0 - purity) * (1.0 - purity) / sig + 
        purity * purity / total
    );
}

double StatisticalCutter::GetSignalToBackgroundError(size_t cutIndex) const {
    size_t sig = survivedSignal_[cutIndex];
    size_t bkg = survivedBackground_[cutIndex];
    
    if (bkg == 0) return 0.0;  // Unikamy dzielenia przez 0
    
    double sb = static_cast<double>(sig) / bkg;
    // Propagacja niepewności dla stosunku, zakładając rozkład Poissona
    return sb * std::sqrt(1.0/sig + 1.0/bkg);
}