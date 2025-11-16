#pragma once
#include <string>
#include <vector>
#include <map>
#include <set>
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
    bool cutValueDynamic = false;
    double centralValue = 0.0;
    bool centralValueDynamic = false;
    std::function<double()> valueGetter;
    std::function<double()> cutValueGetter;
    std::function<double()> centralValueGetter;

    bool isComplexCut = false;
    std::string expression;
    std::shared_ptr<TTreeFormula> treeFormula;
    
    bool isFiducialVolume = false;
    bool isBackgroundRejection = false;
    std::string channel = "default";
    
    // Specjalne cięcia warunkowe
    bool conditionCut = false;              // Czy to cięcie warunkowe (inactive -> skip w UpdateStats)
    std::string condition = "";             // Wyrażenie logiczne na innych conditionCut'ach
    
    // Pole dla syntetycznych cięć grupowych
    bool isSyntheticGroup = false;  // ← NOWE: czy to syntetyczne cięcie grupy
    std::vector<size_t> groupMembers;  // ← NOWE: indeksy cięć wchodzących w skład grupy
    std::string groupCondition = "";  // ← NOWE: warunek na całą grupę background rejection
};

// Enum do wyboru metody normalizacji
enum class NormalizationMode {
    TOTAL_EVENTS,          // Normalizacja do całkowitej liczby zdarzeń
    FIDUCIAL_VOLUME       // Normalizacja do liczby zdarzeń w fiducial volume
};

class StatisticalCutter {
public:
    StatisticalCutter(const std::string& jsonPath, int signalMctruth, KLOE::HypothesisCode hypoCode);
    StatisticalCutter(const std::string& propertiesPath, const std::string& cutsPath, KLOE::HypothesisCode hypoCode);

    void SetTree(TTree* tree);
    
    // Ustawienie trybu normalizacji
    void SetNormalizationMode(NormalizationMode mode) { normalizationMode_ = mode; }
    NormalizationMode GetNormalizationMode() const { return normalizationMode_; }

    // Rejestracja zmiennych zwykłych (z TTree lub pochodnych)
    void RegisterVariableGetter(const std::string& varName, std::function<double()> getter);
    
    // Rejestracja wartości centralnej dla danego cięcia (sprawdza czy nie syntetyczne)
    void RegisterCentralValueGetter(const std::string& cutId, std::function<double()> getter);
    
    // Rejestracja limitu cięcia (dynamiczna wartość cutValue) - sprawdza czy nie syntetyczne
    void RegisterCutValueGetter(const std::string& cutId, std::function<double()> getter);

    // Nowe metody dla kombinacji cięć
    bool PassCutsAnd(const std::vector<size_t>& cutIndices);
    bool PassCutsOr(const std::vector<size_t>& cutIndices);
    bool PassCutsCustom(const std::vector<size_t>& cutIndices, const std::string& logicOperator);

    // Alternatywa - zwróć funkcję do ewaluacji
    std::function<bool()> GetCutCombination(const std::vector<size_t>& cutIndices, const std::string& logicOperator = "&&");

    // Metody dla fiducial volume
    bool PassFiducialVolume() const;
    size_t GetTotalSignalInFV() const { return totalSignalInFV_; }
    size_t GetTotalBackgroundInFV() const { return totalBackgroundInFV_; }
    size_t GetTotalSignalInFVExcludingMinus1() const { return totalSignalInFVExcludingMinus1_; }
    size_t GetTotalBackgroundInFVExcludingMinus1() const { return totalBackgroundInFVExcludingMinus1_; }

    size_t GetTotalSignal() const { return totalSignal_; }
    size_t GetTotalBackground() const { return totalBackground_; }
    size_t GetTotalSignalExcludingMinus1() const { return totalSignalExcludingMinus1_; }
    size_t GetTotalBackgroundExcludingMinus1() const { return totalBackgroundExcludingMinus1_; }

    bool PassCut(size_t cutIndex);
    bool PassAllCuts();
    
    // Sprawdź czy wszystkie cięcia (z listy lub wszystkie jeśli lista pusta) przechodzą
    // Automatycznie obsługuje grupy background rejection z negacją całej grupy
    bool PassCuts(const std::vector<size_t>& cutIndices = std::vector<size_t>());

    void UpdateStats(int mctruth);

    // Metody z automatyczną normalizacją do FV (jeśli włączony)
    double GetEfficiency(size_t cutIndex) const;
    double GetEfficiencyError(size_t cutIndex) const;
    
    double GetPurity(size_t cutIndex) const;
    double GetPurityError(size_t cutIndex) const;
    
    double GetSignalToBackground(size_t cutIndex) const;
    double GetSignalToBackgroundError(size_t cutIndex) const;

    size_t GetSurvivedSignal(size_t cutIndex) const;
    size_t GetSurvivedBackground(size_t cutIndex) const;

    double GetEfficiencyBetweenCuts(size_t cutIndex1, size_t cutIndex2) const;
    double GetEfficiencyBetweenCutsExcludingMctruthMinus1(size_t cutIndex1, size_t cutIndex2) const;
    
    double GetEfficiencyExcludingMctruthMinus1(size_t cutIndex) const;

    // Wersje jawne z wybraną normalizacją (ignorują SetNormalizationMode)
    double GetEfficiencyWithMode(size_t cutIndex, NormalizationMode mode) const;
    double GetEfficiencyErrorWithMode(size_t cutIndex, NormalizationMode mode) const;
    double GetPurityWithMode(size_t cutIndex, NormalizationMode mode) const;
    double GetPurityErrorWithMode(size_t cutIndex, NormalizationMode mode) const;

    const std::vector<Cut>& GetCuts() const { return cuts_; }
    
    // Zwróć listę cięć pomijając członków grup (zwraca tylko widoczne cięcia i syntetyczne)
    std::vector<size_t> GetVisibleCuts() const;

    // Zarządzanie aktywnymi cięciami
    void SetActiveCuts(const std::vector<size_t>& indices) { activeCutIndices_ = indices; }
    void ClearActiveCuts() { activeCutIndices_.clear(); }
    const std::vector<size_t>& GetActiveCuts() const { return activeCutIndices_; }
    
    // Metody pomocnicze do grupowych cięć
    const std::map<std::string, size_t>& GetGroupChannelToIndexMap() const { return groupChannelToSyntheticIndex_; }
    bool IsCutInGroup(size_t cutIndex) const;

private:
    void LoadCuts(const nlohmann::json& j);
    void LoadCuts(const std::string& jsonPath);
    void LoadCutsFromFiles(const std::string& propertiesPath, const std::string& cutsPath);
    bool EvaluateCondition(double value, const Cut& cut) const;
    bool EvaluateExpression(const Cut& cut) const;
    double EvaluateExpressionToDouble(const Cut& cut) const;
    
    // Ewaluacja warunkowych wyrażeń logicznych
    bool EvaluateConditionExpression(const std::string& expression) const;
    
    // Budowanie syntetycznych cięć dla grup background rejection
    void BuildSyntheticGroupCuts();
    
    // Sprawdzenie warunku grupy (groupCondition)
    bool CheckGroupCondition(const Cut& groupCut) const;

    struct CutStats {
        size_t survivedSignal = 0;
        size_t survivedBackground = 0;
        size_t survivedSignalExcludingMinus1 = 0;
        size_t survivedBackgroundExcludingMinus1 = 0;
    };

    void UpdateStatsInternal(int mctruth, std::vector<CutStats>& stats);
    
    // Przetwarzanie cięć z obsługą per-channel grup background rejection
    bool ProcessCutsLoop(const std::vector<size_t>& cutsToProcess, 
                         std::vector<size_t>& survived_array, int mctruth);

    KLOE::HypothesisCode hypoCode_;
    NormalizationMode normalizationMode_ = NormalizationMode::TOTAL_EVENTS;

    std::vector<Cut> cuts_;
    int signalMctruth_;
    // Mapa getterów dla zmiennych pochodnych/pomocniczych
    std::map<std::string, std::function<double()>> variableGetters_;

    std::vector<size_t> survivedSignal_;
    std::vector<size_t> survivedBackground_;
    std::vector<size_t> survivedSignalExcludingMinus1_;
    std::vector<size_t> survivedBackgroundExcludingMinus1_;
    
    // Oddzielne liczniki dla zdarzeń w fiducial volume
    std::vector<size_t> survivedSignalInFV_;
    std::vector<size_t> survivedBackgroundInFV_;
    std::vector<size_t> survivedSignalInFVExcludingMinus1_;
    std::vector<size_t> survivedBackgroundInFVExcludingMinus1_;
    
    size_t totalSignal_ = 0;
    size_t totalBackground_ = 0;
    size_t totalSignalExcludingMinus1_ = 0;
    size_t totalBackgroundExcludingMinus1_ = 0;
    
    // Statystyki dla fiducial volume
    size_t totalSignalInFV_ = 0;
    size_t totalBackgroundInFV_ = 0;
    size_t totalSignalInFVExcludingMinus1_ = 0;
    size_t totalBackgroundInFVExcludingMinus1_ = 0;
    
    // Indeksy aktywnych cięć (jeśli puste - żadne cięcie nie jest aktywne)
    std::vector<size_t> activeCutIndices_;
    
    // Mapowanie channel na indeks syntetycznego cięcia grupy
    std::map<std::string, size_t> groupChannelToSyntheticIndex_;
    
    // Zbiór indeksów cięć które są członkami grup (ukrywamy je z GetCuts)
    std::set<size_t> groupMemberIndices_;
    
    // Mapowanie ID conditionCut na ich indeksy dla szybkiego lookup
    std::map<std::string, size_t> conditionCutIdToIndex_;
    
    TTree* tree_ = nullptr;
};