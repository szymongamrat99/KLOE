#pragma once
#include <string>
#include <vector>
#include <map>
#include <functional>
#include <stdexcept>
#include <json.hpp>
#include <fstream>
#include <algorithm>
#include <memory>
#include <HypothesisCodes.h>

#include <TFormula.h>
#include <TTreeFormula.h>
#include <TTree.h>

struct Cut {
    int order = 0;
    std::string cutId;
    std::string cutType;
    std::string cutCondition;
    double cutValue = 0.0;
    double centralValue = 0.0;
    bool centralValueDynamic = false;
    std::function<double()> valueGetter;
    std::function<double()> centralValueGetter;

    bool isComplexCut = false;
    std::string expression;
    std::shared_ptr<TTreeFormula> treeFormula;
};

class StatisticalCutter {
public:
    StatisticalCutter(const std::string& jsonPath, int signalMctruth, KLOE::HypothesisCode hypoCode);
    StatisticalCutter(const std::string& propertiesPath, const std::string& cutsPath, KLOE::HypothesisCode hypoCode);

    void SetTree(TTree* tree);

    // Rejestracja zmiennych zwykłych (z TTree lub pochodnych)
    void RegisterVariableGetter(const std::string& varName, std::function<double()> getter);
    
    // Rejestracja wartości centralnej dla danego cięcia
    void RegisterCentralValueGetter(const std::string& cutId, std::function<double()> getter);

    // Nowe metody dla kombinacji cięć
    bool PassCutsAnd(const std::vector<size_t>& cutIndices);
    bool PassCutsOr(const std::vector<size_t>& cutIndices);
    bool PassCutsCustom(const std::vector<size_t>& cutIndices, const std::string& logicOperator);

    // Alternatywa - zwróć funkcję do ewaluacji
    std::function<bool()> GetCutCombination(const std::vector<size_t>& cutIndices, const std::string& logicOperator = "&&");


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

    const std::vector<Cut>& GetCuts() const { return cuts_; }

private:
    void LoadCuts(const nlohmann::json& j);
    void LoadCuts(const std::string& jsonPath);
    void LoadCutsFromFiles(const std::string& propertiesPath, const std::string& cutsPath);
    bool EvaluateCondition(double value, const Cut& cut) const;
    bool EvaluateExpression(const Cut& cut) const;
    double EvaluateExpressionToDouble(const Cut& cut) const;

    KLOE::HypothesisCode hypoCode_;

    std::vector<Cut> cuts_;
    int signalMctruth_;
    // Mapa getterów dla zmiennych pochodnych/pomocniczych
    std::map<std::string, std::function<double()>> variableGetters_;

    std::vector<size_t> survivedSignal_;
    std::vector<size_t> survivedBackground_;
    size_t totalSignal_ = 0;
    size_t totalBackground_ = 0;
    
    TTree* tree_ = nullptr;
};