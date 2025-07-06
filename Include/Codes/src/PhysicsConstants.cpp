// PhysicsConstants.cpp
#include <PhysicsConstants.h>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <ConstantsKeys.h>

// --- Template Helper to safely get value with detailed error messages ---
template<typename T>
T PhysicsConstants::getValue(const std::string& jsonPath, const std::string& paramName) const {
    try {
        return m_rawConfigData.at(nlohmann::json::json_pointer(jsonPath)).get<T>();
    } catch (const nlohmann::json::exception& e) {
        throw std::runtime_error("ERROR: Critical config parameter '" + paramName + "' (JSON path: '" + jsonPath + "') missing or invalid type: " + e.what());
    }
}

// --- Template Helper to get value with default, with detailed warnings ---
template<typename T>
T PhysicsConstants::getValueOrDefault(const std::string& jsonPath, const std::string& paramName, T defaultValue) const {
    try {
        return m_rawConfigData.at(nlohmann::json::json_pointer(jsonPath)).get<T>();
    } catch (const nlohmann::json::exception& e) {
        std::cerr << "WARNING: Optional config parameter '" << paramName << "' (JSON path: '" << jsonPath << "') not found or invalid type. Using default value." << std::endl;
        return defaultValue;
    }
}


// --- Private Method to Load JSON File ---
void PhysicsConstants::loadJsonFile(const std::string& filePath) {
    std::cout << "INFO: Attempting to load experiment configuration from: " << filePath << std::endl;
    std::ifstream fileStream(filePath);

    if (!fileStream.is_open()) {
        throw std::runtime_error("ERROR: Could not open experiment configuration file: " + filePath);
    }

    try {
        fileStream >> m_rawConfigData;
        std::cout << "INFO: Successfully parsed JSON from " << filePath << std::endl;
    } catch (const nlohmann::json::parse_error& e) {
        throw std::runtime_error("ERROR: JSON parsing error in experiment configuration file (" + filePath + "): " + e.what());
    } catch (const std::exception& e) {
        throw std::runtime_error("ERROR: Failed to load experiment configuration from " + filePath + ": " + e.what());
    }
}

// --- Private Method to Populate Member Variables from JSON ---
void PhysicsConstants::populateMembersFromJson() {
    std::cout << "INFO: Populating PhysicsConstants member variables from loaded JSON." << std::endl;

    // Detector Parameters
    mK0        = getValue<double>("/values/" + ConstantsKeys::K0mass, "Neutral Kaon Mass");
    tau_S_nonCPT   = getValue<double>("/values/" + ConstantsKeys::TauS, "K-short lifetime (non-CPT)") * 1E9; // Convert to ns
    tau_L   = getValue<double>("/values/" + ConstantsKeys::TauL, "K-long lifetime (non-CPT)") * 1E9; // Convert to ns
    delta_mass_nonCPT = getValue<double>("/values/" + ConstantsKeys::deltaM, "Delta Mass (non-CPT)");

    mod_epsilon = getValue<double>("/values/" + ConstantsKeys::modEps, "Modulus of Epsilon");

    Re = getValue<double>("/values/" + ConstantsKeys::Re, "Real part of Epsilon");

    Im_nonCPT = getValue<double>("/values/" + ConstantsKeys::Im, "Imaginary part of Epsilon (non-CPT)") * (M_PI / 180.); // Convert degrees to radians

    phi_pm_nonCPT = getValue<double>("/values/" + ConstantsKeys::phiPM, "Phi PM (non-CPT)");

    phi_00_nonCPT = getValue<double>("/values/" + ConstantsKeys::phi00, "Phi 00 (non-CPT)");

    std::cout << "INFO: All PhysicsConstants member variables populated." << std::endl;

    // Clear the raw JSON data if no longer needed to save memory
    m_rawConfigData.clear();
}

// --- Constructor ---
PhysicsConstants::PhysicsConstants(const std::string& jsonFilePath) {
    loadJsonFile(jsonFilePath); // Load the raw JSON first
    populateMembersFromJson();  // Then populate the member variables
}