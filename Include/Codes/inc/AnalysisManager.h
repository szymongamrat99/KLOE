#ifndef ANALYSIS_CONFIG_H
#define ANALYSIS_CONFIG_H

#include <HypothesisCodes.h>
#include <string>
#include <vector>
#include <map>
#include <kloe_class.h>

namespace KLOE {

// ===== RUN MODE (tylko to dodajemy nowe) =====

enum class RunMode {
    FULL,
    QUICK,
    DEBUG,
    CUTS_ONLY,
    SYSTEMATICS
};

inline std::string RunModeToString(RunMode mode) {
    switch (mode) {
        case RunMode::FULL:        return "FULL";
        case RunMode::QUICK:       return "QUICK";
        case RunMode::DEBUG:       return "DEBUG";
        case RunMode::CUTS_ONLY:   return "CUTS_ONLY";
        case RunMode::SYSTEMATICS: return "SYSTEMATICS";
        default:                   return "UNKNOWN";
    }
}

inline RunMode StringToRunMode(const std::string& str) {
    if (str == "FULL")        return RunMode::FULL;
    if (str == "QUICK")       return RunMode::QUICK;
    if (str == "DEBUG")       return RunMode::DEBUG;
    if (str == "CUTS_ONLY")   return RunMode::CUTS_ONLY;
    if (str == "SYSTEMATICS") return RunMode::SYSTEMATICS;
    return RunMode::FULL;
}

// ===== KONFIGURACJA (proste struktury) =====

class AnalysisConfig
{
public:
    // Singleton
    static AnalysisConfig& getInstance() {
        static AnalysisConfig instance;
        return instance;
    }
    
    // Delete copy/move
    AnalysisConfig(const AnalysisConfig&) = delete;
    AnalysisConfig& operator=(const AnalysisConfig&) = delete;
    
    // ===== STRUKTURY =====
    
    struct Metadata {
        std::string version;
        std::string description;
        std::string lastModified;
        std::string author;
    };
    
    struct HypothesisModules {
        bool momentumSmearing = true;
        bool trilaterationKinFit = true;
        bool signalKinFit = true;
        bool omegaKinFit = false;
        bool triangleReconstruction = true;
        bool photonPairing = true;
        bool kaonProperTimes = true;
    };
    
    struct HypothesisCuts {
        double minMassWindow = 480.0;
        double maxMassWindow = 540.0;
        double maxChi2 = 50.0;
        bool applyPreselection = true;
        bool applyKinematicFit = true;
        bool applyStatisticalCuts = true;
    };
    
    struct HypothesisConfig {
        bool enabled = false;
        std::string description;
        HypothesisModules modules;
        HypothesisCuts cuts;
    };
    
    struct OutputConfig {
        bool savePulls = true;
        bool saveIntermediateResults = false;
        bool saveDebugInfo = false;
        bool saveMCTruthAlways = true;
        int verboseLevel = 1;
    };
    
    struct PerformanceConfig {
        bool parallelProcessing = false;
        int nThreads = 4;
    };
    
    struct ValidationConfig {
        bool checkPhysicalBounds = true;
        bool checkEnergyMomentum = true;
    };
    
    // ===== PUBLIC MEMBERS =====
    
    Metadata metadata;
    RunMode runMode = RunMode::FULL;
    HypothesisCode activeHypothesis = HypothesisCode::SIGNAL;
    std::vector<HypothesisCode> enabledHypotheses;
    std::map<HypothesisCode, HypothesisConfig> hypothesisConfigs;
    OutputConfig output;
    PerformanceConfig performance;
    ValidationConfig validation;
    
    // ===== METODY =====
    
    bool LoadFromFile(const std::string& filename);
    const HypothesisConfig& GetActiveHypothesisConfig() const;
    const HypothesisConfig& GetHypothesisConfig(HypothesisCode hyp) const;
    
    bool ShouldRunModule(const std::string& moduleName) const;
    bool ShouldRunKinFit() const;
    bool ShouldSaveDetails() const;
    bool IsDebugMode() const;
    
    void Print() const;
    
private:
    AnalysisConfig() = default;
};

}  // namespace KLOE

#endif  // ANALYSIS_CONFIG_H