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
    if (str == "THREE_PI0") return HypothesisCode::THREE_PI0;
    
    std::cerr << "Unknown hypothesis." << std::endl;
    return HypothesisCode::INVALID_VALUE;
}

std::string HypothesisCodeToString(HypothesisCode code) {
    switch (code) {
        case HypothesisCode::SIGNAL:  return "SIGNAL";
        case HypothesisCode::OMEGAPI: return "OMEGAPI";
        case HypothesisCode::FOUR_PI: return "FOUR_PI";
        case HypothesisCode::THREE_PI0: return "THREE_PI0";
        default: return "UNKNOWN";
    }
}

// ===== KONSTRUKTOR I DESTRUKTOR =====
void AnalysisConfig::SetupLogger(ErrorHandling::ErrorLogs* logger) {
    _logger = logger;
}

// ===== LOAD FROM JSON =====
bool AnalysisConfig::LoadFromFile(const std::string& filename) {
    try {
        std::ifstream file(filename);
        if (!file.is_open()) {
          _lastErrorCode = ErrorHandling::ErrorCodes::FILE_NOT_EXIST;
          _logger->getErrLog(_lastErrorCode, "Cannot open config file: " + filename);
          return false;
        }
        
        nlohmann::json j;
        file >> j;
        
        _lastInfoCode = ErrorHandling::InfoCodes::FILE_ADDED;
        _logger->getLog(_lastInfoCode, "Config file loaded: " + filename);
        
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
        
        _lastInfoCode = ErrorHandling::InfoCodes::CONFIG_LOADED;
        _logger->getLog(_lastInfoCode, Form("Hypothesis %s loaded with run mode %s", 
                  HypothesisCodeToString(activeHypothesis).c_str(), 
                  RunModeToString(runMode).c_str()));
        return true;
        
    } catch (const std::exception& e) {
        _lastErrorCode = ErrorHandling::ErrorCodes::FILE_NOT_EXIST; // You might want a different code for JSON parsing errors
        _logger->getErrLog(_lastErrorCode, std::string("Exception while loading config: ") + e.what());
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

void AnalysisConfig::PrintToScreen() const {
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

void AnalysisConfig::Print() const {

    std::string currentConfiguration = "";

    currentConfiguration += "\n╔═══════════════════════════════════════════╗\n";
    currentConfiguration += "║     KLOE Analysis Configuration           ║\n";
    currentConfiguration += "╚═══════════════════════════════════════════╝\n";

    if (!metadata.author.empty()) {
        currentConfiguration += "Metadata:\n";
        currentConfiguration += "   Author: " + metadata.author + "\n";
        currentConfiguration += "   Version: " + metadata.version + "\n";
    }

    currentConfiguration += "Analysis Setup:\n";
    currentConfiguration += "   Run Mode: " + RunModeToString(runMode) + "\n";
    currentConfiguration += "   Hypothesis: " + HypothesisCodeToString(activeHypothesis) + "\n";
    
    const auto& hypConfig = GetActiveHypothesisConfig();
    if (!hypConfig.description.empty()) {
        currentConfiguration += "   (" + hypConfig.description + ")\n";
    }
    
    currentConfiguration += "Modules:\n";
    currentConfiguration += "   Classify MC variables:       " + std::string(hypConfig.modules.classifyMCVariables ? "true" : "false") + "\n";
    currentConfiguration += "   Signal only:       " + std::string(hypConfig.modules.signalOnly ? "true" : "false") + "\n";
    currentConfiguration += "   Momentum Smearing:       " + std::string(hypConfig.modules.momentumSmearing ? "true" : "false") + "\n";
    currentConfiguration += "   Trilateration KinFit:    " + std::string(hypConfig.modules.trilaterationKinFit ? "true" : "false") + "\n";
    currentConfiguration += "   Signal KinFit:           " + std::string(hypConfig.modules.signalKinFit ? "true" : "false") + "\n";
    currentConfiguration += "   Omega KinFit:            " + std::string(hypConfig.modules.omegaKinFit ? "true" : "false") + "\n";
    currentConfiguration += "   Triangle Reconstruction: " + std::string(hypConfig.modules.triangleReconstruction ? "true" : "false") + "\n";
    currentConfiguration += "   Photon Pairing:          " + std::string(hypConfig.modules.photonPairing ? "true" : "false") + "\n";
    currentConfiguration += "   Kaon Proper Times:       " + std::string(hypConfig.modules.kaonProperTimes ? "true" : "false") + "\n";

    currentConfiguration += "Cuts:\n";
    currentConfiguration += "   Mass window: [" + std::to_string(hypConfig.cuts.minMassWindow) + ", " + std::to_string(hypConfig.cuts.maxMassWindow) + "] MeV\n";
    currentConfiguration += "   Max χ²: " + std::to_string(hypConfig.cuts.maxChi2) + "\n";

    currentConfiguration += "Output:\n";
    currentConfiguration += "   Save Pulls:          " + std::string(output.savePulls ? "true" : "false") + "\n";
    currentConfiguration += "   Save MC Truth:       " + std::string(output.saveMCTruthAlways ? "true" : "false") + "\n";
    currentConfiguration += "   Verbose Level:       " + std::to_string(output.verboseLevel) + "\n";

    currentConfiguration += "\n══════════════════════════════════════════════\n\n";

    std::string beginMsg = "Current analysis configuration.";
    std::string endMsg = "End of configuration details.";
    _logger->prettyPrint(currentConfiguration, beginMsg, endMsg, ErrorHandling::LogFiles::LogType::ANALYSIS_CONFIG);
}

}  // namespace KLOE