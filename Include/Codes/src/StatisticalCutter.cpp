#include "StatisticalCutter.h"
#include <iostream>

using json = nlohmann::json;

StatisticalCutter::StatisticalCutter(const std::string& jsonPath, int signalMctruth, KLOE::HypothesisCode hypoCode)
    : signalMctruth_(signalMctruth), hypoCode_(hypoCode)
{
    LoadCuts(jsonPath);
    survivedSignal_.resize(cuts_.size(), 0);
    survivedBackground_.resize(cuts_.size(), 0);
    survivedSignalExcludingMinus1_.resize(cuts_.size(), 0);
    survivedBackgroundExcludingMinus1_.resize(cuts_.size(), 0);
}

StatisticalCutter::StatisticalCutter(const std::string& propertiesPath, const std::string& cutsPath, KLOE::HypothesisCode hypoCode)
    : hypoCode_(hypoCode)
{
    LoadCutsFromFiles(propertiesPath, cutsPath);
    survivedSignal_.resize(cuts_.size(), 0);
    survivedBackground_.resize(cuts_.size(), 0);
    survivedSignalExcludingMinus1_.resize(cuts_.size(), 0);
    survivedBackgroundExcludingMinus1_.resize(cuts_.size(), 0);
}

void StatisticalCutter::SetTree(TTree* tree) {
    tree_ = tree;
    
    // Rekompiluj TreeFormulas z nowym drzewem
    for (auto& cut : cuts_) {
        if (cut.isComplexCut && !cut.expression.empty()) {
            try {
                cut.treeFormula = std::make_unique<TTreeFormula>(
                    cut.cutId.c_str(),
                    cut.expression.c_str(),
                    tree_
                );
                if (cut.treeFormula->GetNdim() <= 0) {
                    std::cerr << "Warning: TreeFormula for cut " << cut.cutId 
                              << " has no dimensions. Expression: " << cut.expression << std::endl;
                }
            } catch (const std::exception& e) {
                std::cerr << "Error creating TTreeFormula for cut " << cut.cutId 
                          << ": " << e.what() << std::endl;
            }
        }
    }
}

void StatisticalCutter::RegisterVariableGetter(const std::string& varName, std::function<double()> getter) {
    variableGetters_[varName] = getter;
    // Aktualizuj gettery dla cięć które mogą używać tej zmiennej
    for (auto& cut : cuts_) {
        if (cut.cutId == varName && !cut.isComplexCut) {
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
        case KLOE::HypothesisCode::SIMONA_ANALYSIS:
            cutListKey = "SIMONA_ANALYSIS";
            break;
        default:
            std::cout << "Warning: Using default cuts for unknown hypothesis code" << std::endl;
            return;
    }

    const json* listRoot = &j;
    if (j.contains("listOfCuts")) {
        if (!j["listOfCuts"].contains(cutListKey)) {
            std::cout << "Warning: No cuts found for hypothesis " << cutListKey << std::endl;
            return;
        }
        listRoot = &j["listOfCuts"][cutListKey];
    }

    for (const auto& cutj : *listRoot) {
        Cut cut;
        cut.order = cutj.value("order", 0);
        cut.cutId = cutj.value("cutId", std::string());
        
        if (cutj.contains("expression")) {
            // Złożone cięcie - wyrażenie na zmiennych
            cut.isComplexCut = true;
            cut.expression = cutj["expression"].get<std::string>();
        } else {
            // Zwykłe cięcie na pojedynczej zmiennej
            cut.cutType = cutj.value("cutType", std::string());
            cut.cutCondition = cutj.value("cutCondition", std::string());
            cut.cutValue = cutj.value("cutValue", 0.0);
            cut.centralValue = cutj.value("centralValue", 0.0);
            cut.centralValueDynamic = false;
        }
        cuts_.push_back(cut);
    }
    
    std::sort(cuts_.begin(), cuts_.end(), [](const Cut& a, const Cut& b) {
        return a.order < b.order;
    });
}

void StatisticalCutter::LoadCutsFromFiles(const std::string& propertiesPath, const std::string& cutsPath) {
    json props, cuts;
    {
        std::ifstream f(propertiesPath);
        if (!f.is_open()) throw std::runtime_error("Cannot open Utils::properties file: " + propertiesPath);
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
    if (cut.isComplexCut && !cut.expression.empty()) {
        return EvaluateExpression(cut);
    }

    double central = cut.centralValueDynamic && cut.centralValueGetter
        ? cut.centralValueGetter()
        : cut.centralValue;

    if (cut.cutType == "O") {
        if (cut.cutCondition == "<=") return value <= cut.cutValue;
        if (cut.cutCondition == ">=") return value >= cut.cutValue;
        if (cut.cutCondition == "<")  return value < cut.cutValue;
        if (cut.cutCondition == ">")  return value > cut.cutValue;
    }
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

double StatisticalCutter::EvaluateExpressionToDouble(const Cut& cut) const {
    if (!tree_) {
        std::cerr << "Warning: No TTree set for cut " << cut.cutId << std::endl;
        return 0.0;
    }

    if (!cut.treeFormula) {
        std::cerr << "Warning: No compiled TTreeFormula for cut " << cut.cutId << std::endl;
        return 0.0;
    }

    try {
        return cut.treeFormula->EvalInstance();
    } catch (const std::exception& e) {
        std::cerr << "Error evaluating TTreeFormula for cut " << cut.cutId 
                  << ": " << e.what() << std::endl;
        return 0.0;
    }
}

bool StatisticalCutter::EvaluateExpression(const Cut& cut) const {
    double result = EvaluateExpressionToDouble(cut);
    return result != 0.0;
}

bool StatisticalCutter::PassCut(size_t cutIndex) {
    if (cutIndex >= cuts_.size())
        throw std::out_of_range("Cut index out of range");
    
    const auto& cut = cuts_[cutIndex];
    
    if (cut.isComplexCut) {
        return EvaluateExpression(cut);
    }
    
    // Dla zwykłych cięć - szukaj gettera
    if (cut.valueGetter) {
        double value = cut.valueGetter();
        return EvaluateCondition(value, cut);
    }
    
    // Spróbuj znaleźć w mapie getterów
    auto it = variableGetters_.find(cut.cutId);
    if (it != variableGetters_.end()) {
        double value = it->second();
        return EvaluateCondition(value, cut);
    }
    
    throw std::runtime_error("No getter registered for variable: " + cut.cutId);
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

    if (mctruth != -1) {
        if (mctruth == signalMctruth_ || mctruth == 0)
            totalSignalExcludingMinus1_++;
        else
            totalBackgroundExcludingMinus1_++;

        survived = true;
        for (size_t i = 0; i < cuts_.size(); ++i) {
            if (survived && PassCut(i)) {
                if (mctruth == signalMctruth_)
                    survivedSignalExcludingMinus1_[i]++;
                else
                    survivedBackgroundExcludingMinus1_[i]++;
            } else {
                survived = false;
            }
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
    return std::sqrt(eff * (1.0 - eff) / totalSignal_);
}

double StatisticalCutter::GetPurityError(size_t cutIndex) const {
    size_t sig = survivedSignal_[cutIndex];
    size_t bkg = survivedBackground_[cutIndex];
    double total = sig + bkg;
    
    if (total == 0) return 0.0;
    
    double purity = static_cast<double>(sig) / total;
    if (sig == 0) return 0.0;
    
    return purity * std::sqrt(
        (1.0 - purity) * (1.0 - purity) / sig + 
        purity * purity / total
    );
}

double StatisticalCutter::GetSignalToBackgroundError(size_t cutIndex) const {
    size_t sig = survivedSignal_[cutIndex];
    size_t bkg = survivedBackground_[cutIndex];
    
    if (bkg == 0) return 0.0;
    
    double sb = static_cast<double>(sig) / bkg;
    return sb * std::sqrt( (sig>0?1.0/sig:0.0) + 1.0/bkg );
}

bool StatisticalCutter::PassCutsAnd(const std::vector<size_t>& cutIndices) {
    for (size_t idx : cutIndices) {
        if (!PassCut(idx))
            return false;
    }
    return true;
}

bool StatisticalCutter::PassCutsOr(const std::vector<size_t>& cutIndices) {
    for (size_t idx : cutIndices) {
        if (PassCut(idx))
            return true;
    }
    return false;
}

bool StatisticalCutter::PassCutsCustom(const std::vector<size_t>& cutIndices, const std::string& logicOperator) {
    if (logicOperator == "&&") {
        return PassCutsAnd(cutIndices);
    } else if (logicOperator == "||") {
        return PassCutsOr(cutIndices);
    }
    throw std::runtime_error("Unknown logic operator: " + logicOperator);
}

std::function<bool()> StatisticalCutter::GetCutCombination(const std::vector<size_t>& cutIndices, const std::string& logicOperator) {
    // Zwróć lambdę która ewaluuje kombinację cięć
    return [this, cutIndices, logicOperator]() {
        return PassCutsCustom(cutIndices, logicOperator);
    };
}

double StatisticalCutter::GetEfficiencyBetweenCuts(size_t cutIndex1, size_t cutIndex2) const {
    if (cutIndex1 >= cuts_.size() || cutIndex2 >= cuts_.size())
        throw std::out_of_range("Cut index out of range");
    
    size_t first = std::min(cutIndex1, cutIndex2);
    size_t second = std::max(cutIndex1, cutIndex2);
    
    if (survivedSignal_[first] == 0) return 0.0;
    return static_cast<double>(survivedSignal_[second]) / survivedSignal_[first];
}

double StatisticalCutter::GetEfficiencyBetweenCutsExcludingMctruthMinus1(size_t cutIndex1, size_t cutIndex2) const {
    if (cutIndex1 >= cuts_.size() || cutIndex2 >= cuts_.size())
        throw std::out_of_range("Cut index out of range");
    
    size_t first = std::min(cutIndex1, cutIndex2);
    size_t second = std::max(cutIndex1, cutIndex2);
    
    if (survivedSignalExcludingMinus1_[first] == 0) return 0.0;
    return static_cast<double>(survivedSignalExcludingMinus1_[second]) / survivedSignalExcludingMinus1_[first];
}

double StatisticalCutter::GetEfficiencyExcludingMctruthMinus1(size_t cutIndex) const {
    if (totalSignalExcludingMinus1_ == 0) return 0.0;
    return static_cast<double>(survivedSignalExcludingMinus1_[cutIndex]) / totalSignalExcludingMinus1_;
}