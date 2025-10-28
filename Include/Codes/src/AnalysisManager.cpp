#include "AnalysisManager.h"
#include <nlohmann/json.hpp>
#include <fstream>
#include <iostream>
#include <iomanip>


namespace KLOE {

// ===== HELPER dla konwersji HypothesisCode <-> string =====
// (wykorzystuje Twoją istniejącą funkcję!)


HypothesisCode JsonToHypothesisCode(const std::string& str) {
    // Możesz użyć Twojej funkcji pm00::StringToHypothesisCode
    // lub zrobić prostą mapę:
    if (str == "SIGNAL")   return HypothesisCode::SIGNAL;
    if (str == "OMEGAPI")  return HypothesisCode::OMEGAPI;
    if (str == "FOUR_PI")  return HypothesisCode::FOUR_PI;
    
    std::cerr << "Unknown hypothesis." << std::endl;
    return HypothesisCode::INVALID_VALUE;
}

std::string HypothesisCodeToString(HypothesisCode code) {
    switch (code) {
        case HypothesisCode::SIGNAL:  return "SIGNAL";
        case HypothesisCode::OMEGAPI: return "OMEGAPI";
        case HypothesisCode::FOUR_PI: return "FOUR_PI";
        default: return "UNKNOWN";
    }
}

// ===== LOAD FROM JSON =====

bool AnalysisConfig::LoadFromFile(const std::string& filename) {
    try {
        std::ifstream file(filename);
        if (!file.is_open()) {
            std::cerr << "ERROR: Cannot open config file: " << filename << std::endl;
            return false;
        }
        
        nlohmann::json j;
        file >> j;
        
        std::cout << "Loading analysis config from: " << filename << std::endl;
        
        // 1. Metadata
        if (j.contains("metadata")) {
            metadata.version = j["metadata"].value("version", "1.0");
            metadata.description = j["metadata"].value("description", "");
            metadata.lastModified = j["metadata"].value("lastModified", "");
            metadata.author = j["metadata"].value("author", "");
        }
        
        // 2. Run mode
        if (j.contains("runMode")) {
            runMode = StringToRunMode(j["runMode"]);
        }
        
        // 3. Active hypothesis
        if (j.contains("hypothesis")) {
            if (j["hypothesis"].contains("active")) {
                activeHypothesis = JsonToHypothesisCode(j["hypothesis"]["active"]);
            }
            
            if (j["hypothesis"].contains("enabled")) {
                enabledHypotheses.clear();
                for (const auto& hypStr : j["hypothesis"]["enabled"]) {
                    enabledHypotheses.push_back(JsonToHypothesisCode(hypStr));
                }
            }
        }
        
        // 4. Hypothesis configs
        if (j.contains("hypothesisConfigs")) {
            for (auto item : j["hypothesisConfigs"].items()) {
                HypothesisCode hyp = JsonToHypothesisCode(item.key());
                HypothesisConfig config;
                
                config.enabled = item.value().value("enabled", false);
                config.description = item.value().value("description", "");
                config.signal = item.value().value("signal", 1);
                
                // Modules
                if (item.value().contains("modules")) {
                    auto& mod = item.value()["modules"];
                    config.modules.classifyMCVariables = mod.value("classifyMCVariables", true);
                    config.modules.signalOnly = mod.value("signalOnly", false);
                    config.modules.momentumSmearing = mod.value("momentumSmearing", true);
                    config.modules.trilaterationKinFit = mod.value("trilaterationKinFit", true);
                    config.modules.signalKinFit = mod.value("signalKinFit", true);
                    config.modules.omegaKinFit = mod.value("omegaKinFit", false);
                    config.modules.triangleReconstruction = mod.value("triangleReconstruction", true);
                    config.modules.photonPairing = mod.value("photonPairing", true);
                    config.modules.kaonProperTimes = mod.value("kaonProperTimes", true);
                }
                
                // Cuts
                if (item.value().contains("cuts")) {
                    auto& cuts = item.value()["cuts"];
                    config.cuts.minMassWindow = cuts.value("minMassWindow", 480.0);
                    config.cuts.maxMassWindow = cuts.value("maxMassWindow", 540.0);
                    config.cuts.maxChi2 = cuts.value("maxChi2", 50.0);
                    config.cuts.applyPreselection = cuts.value("applyPreselection", true);
                    config.cuts.applyKinematicFit = cuts.value("applyKinematicFit", true);
                    config.cuts.applyStatisticalCuts = cuts.value("applyStatisticalCuts", true);
                }
                
                hypothesisConfigs[hyp] = config;
            }
        }
        
        // 5. Output
        if (j.contains("output")) {
            auto& out = j["output"];
            output.savePulls = out.value("savePulls", true);
            output.saveIntermediateResults = out.value("saveIntermediateResults", false);
            output.saveDebugInfo = out.value("saveDebugInfo", false);
            output.saveMCTruthAlways = out.value("saveMCTruthAlways", true);
            output.verboseLevel = out.value("verboseLevel", 1);
        }
        
        // 6. Performance
        if (j.contains("performance")) {
            auto& perf = j["performance"];
            performance.parallelProcessing = perf.value("parallelProcessing", false);
            performance.nThreads = perf.value("nThreads", 4);
        }
        
        // 7. Validation
        if (j.contains("validation")) {
            auto& val = j["validation"];
            validation.checkPhysicalBounds = val.value("checkPhysicalBounds", true);
            validation.checkEnergyMomentum = val.value("checkEnergyMomentum", true);
        }
        
        std::cout << "Config loaded: " << HypothesisCodeToString(activeHypothesis) 
                  << " (" << RunModeToString(runMode) << ")" << std::endl;
        return true;
        
    } catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << std::endl;
        return false;
    }
}

// ===== GETTERS =====

const AnalysisConfig::HypothesisConfig& 
AnalysisConfig::GetActiveHypothesisConfig() const {
    return GetHypothesisConfig(activeHypothesis);
}

const AnalysisConfig::HypothesisConfig& 
AnalysisConfig::GetHypothesisConfig(HypothesisCode hyp) const {
    auto it = hypothesisConfigs.find(hyp);
    if (it != hypothesisConfigs.end()) {
        return it->second;
    }
    
    static HypothesisConfig defaultConfig;
    return defaultConfig;
}

// ===== HELPERS =====

bool AnalysisConfig::ShouldRunModule(const std::string& moduleName) const {
    const auto& hypConfig = GetActiveHypothesisConfig();
    
    if (moduleName == "momentumSmearing")        return hypConfig.modules.momentumSmearing;
    if (moduleName == "trilaterationKinFit")     return hypConfig.modules.trilaterationKinFit;
    if (moduleName == "signalKinFit")           return hypConfig.modules.signalKinFit;
    if (moduleName == "omegaKinFit")            return hypConfig.modules.omegaKinFit;
    if (moduleName == "triangleReconstruction")  return hypConfig.modules.triangleReconstruction;
    if (moduleName == "photonPairing")          return hypConfig.modules.photonPairing;
    if (moduleName == "kaonProperTimes")        return hypConfig.modules.kaonProperTimes;
    
    return false;
}

bool AnalysisConfig::ShouldRunKinFit() const {
    if (runMode == RunMode::QUICK) return false;
    
    const auto& hypConfig = GetActiveHypothesisConfig();
    return hypConfig.modules.trilaterationKinFit || 
           hypConfig.modules.signalKinFit || 
           hypConfig.modules.omegaKinFit;
}

bool AnalysisConfig::ShouldSaveDetails() const {
    return runMode == RunMode::DEBUG || output.saveIntermediateResults;
}

bool AnalysisConfig::IsDebugMode() const {
    return runMode == RunMode::DEBUG || output.verboseLevel >= 2;
}

// ===== PRINT =====

void AnalysisConfig::Print() const {
    std::cout << "\n╔═══════════════════════════════════════════╗" << std::endl;
    std::cout << "║     KLOE Analysis Configuration           ║" << std::endl;
    std::cout << "╚═══════════════════════════════════════════╝" << std::endl;
    
    if (!metadata.author.empty()) {
        std::cout << "Metadata:" << std::endl;
        std::cout << "   Author: " << metadata.author << std::endl;
        std::cout << "   Version: " << metadata.version << std::endl;
    }
    
    std::cout << "Analysis Setup:" << std::endl;
    std::cout << "   Run Mode: " << RunModeToString(runMode) << std::endl;
    std::cout << "   Hypothesis: " << HypothesisCodeToString(activeHypothesis) << std::endl;
    
    const auto& hypConfig = GetActiveHypothesisConfig();
    if (!hypConfig.description.empty()) {
        std::cout << "   (" << hypConfig.description << ")" << std::endl;
    }
    
    std::cout << "Modules:" << std::endl;
    std::cout << "   Classify MC variables:       " << (hypConfig.modules.classifyMCVariables ? "true" : "false") << std::endl;
    std::cout << "   Signal only:       " << (hypConfig.modules.signalOnly ? "true" : "false") << std::endl;
    std::cout << "   Momentum Smearing:       " << (hypConfig.modules.momentumSmearing ? "true" : "false") << std::endl;
    std::cout << "   Trilateration KinFit:    " << (hypConfig.modules.trilaterationKinFit ? "true" : "false") << std::endl;
    std::cout << "   Signal KinFit:           " << (hypConfig.modules.signalKinFit ? "true" : "false") << std::endl;
    std::cout << "   Omega KinFit:            " << (hypConfig.modules.omegaKinFit ? "true" : "false") << std::endl;
    std::cout << "   Triangle Reconstruction: " << (hypConfig.modules.triangleReconstruction ? "true" : "false") << std::endl;
    std::cout << "   Photon Pairing:          " << (hypConfig.modules.photonPairing ? "true" : "false") << std::endl;
    std::cout << "   Kaon Proper Times:       " << (hypConfig.modules.kaonProperTimes ? "true" : "false") << std::endl;
    
    std::cout << "Cuts:" << std::endl;
    std::cout << "   Mass window: [" << hypConfig.cuts.minMassWindow 
              << ", " << hypConfig.cuts.maxMassWindow << "] MeV" << std::endl;
    std::cout << "   Max χ²: " << hypConfig.cuts.maxChi2 << std::endl;
    
    std::cout << "Output:" << std::endl;
    std::cout << "   Save Pulls:          " << (output.savePulls ? "true" : "false") << std::endl;
    std::cout << "   Save MC Truth:       " << (output.saveMCTruthAlways ? "true" : "false") << std::endl;
    std::cout << "   Verbose Level:       " << output.verboseLevel << std::endl;
    
    std::cout << "\n══════════════════════════════════════════════\n" << std::endl;
}

}  // namespace KLOE