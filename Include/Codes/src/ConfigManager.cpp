#include "../inc/ConfigManager.h"
#include <iostream>
#include <const.h>
#include <kloe_class.h>

// Static member definitions - moved here to avoid multiple definition errors
ConfigManager* ConfigManager::instance = nullptr;
std::mutex ConfigManager::mutex_;

// Private constructor implementation
ConfigManager::ConfigManager() : propertiesLoaded(false), constantsLoaded(false), errorLogger(nullptr) {
    // Initialize error logger with config-specific log file
    KLOE::pm00 Obj;

    std::string logFileName = Paths::logsCNAFDir + "config_manager_" + Obj.getCurrentDate() + ".log";
    errorLogger = new ErrorHandling::ErrorLogs(logFileName);
    errorLogger->setPrintToScreen(true);
    
    loadConfigurations();
}

// Destructor implementation
ConfigManager::~ConfigManager() {
    if (errorLogger != nullptr) {
        delete errorLogger;
        errorLogger = nullptr;
    }
}

// Load both properties and constants JSON files
void ConfigManager::loadConfigurations() {
    loadProperties();
    loadConstants();
}

// Load properties.json file
void ConfigManager::loadProperties() {
    try {
        const char* propertiesPath = getenv("PROPERTIESKLOE");
        if (!propertiesPath) {
            ErrorHandling::ErrorCodes errCode = ErrorHandling::ErrorCodes::FILE_NOT_EXIST;
            errorLogger->getErrLog(errCode, "PROPERTIESKLOE environment variable not set");
            return;
        }
        
        propertiesFilePath = std::string(propertiesPath) + "/properties.json";
        std::ifstream propertyFile(propertiesFilePath.c_str());
        
        if (!propertyFile.is_open()) {
            ErrorHandling::ErrorCodes errCode = ErrorHandling::ErrorCodes::FILE_NOT_EXIST;
            errorLogger->getErrLog(errCode, "Cannot open properties file: " + propertiesFilePath);
            return;
        }
        
        properties = json::parse(propertyFile);
        propertiesLoaded = true;
        
        ErrorHandling::InfoCodes infoCode = ErrorHandling::InfoCodes::FILE_ADDED;
        errorLogger->getLog(infoCode, "Properties loaded successfully from: " + propertiesFilePath);
        
    } catch (const json::parse_error& e) {
        ErrorHandling::ErrorCodes errCode = ErrorHandling::ErrorCodes::DATA_TYPE;
        errorLogger->getErrLog(errCode, "JSON parse error in properties: " + std::string(e.what()));
    } catch (const std::exception& e) {
        ErrorHandling::ErrorCodes errCode = ErrorHandling::ErrorCodes::NOT_RECOGNIZED;
        errorLogger->getErrLog(errCode, "Unexpected error loading properties: " + std::string(e.what()));
    }
}

// Load pdg_const.json file
void ConfigManager::loadConstants() {
    try {
        const char* pdgApiPath = getenv("PDGAPI");
        if (!pdgApiPath) {
            ErrorHandling::ErrorCodes errCode = ErrorHandling::ErrorCodes::FILE_NOT_EXIST;
            errorLogger->getErrLog(errCode, "PDGAPI environment variable not set");
            return;
        }
        
        std::string pdgConstFilePath = std::string(pdgApiPath) + "/pdg_const.json";
        std::ifstream fconst(pdgConstFilePath);
        
        if (!fconst.is_open()) {
            ErrorHandling::ErrorCodes errCode = ErrorHandling::ErrorCodes::FILE_NOT_EXIST;
            errorLogger->getErrLog(errCode, "Cannot open constants file: " + pdgConstFilePath);
            return;
        }
        
        constants = json::parse(fconst);
        constantsLoaded = true;
        
        ErrorHandling::InfoCodes infoCode = ErrorHandling::InfoCodes::FILE_ADDED;
        errorLogger->getLog(infoCode, "Constants loaded successfully from: " + pdgConstFilePath);
        
    } catch (const json::parse_error& e) {
        ErrorHandling::ErrorCodes errCode = ErrorHandling::ErrorCodes::DATA_TYPE;
        errorLogger->getErrLog(errCode, "JSON parse error in constants: " + std::string(e.what()));
    } catch (const std::exception& e) {
        ErrorHandling::ErrorCodes errCode = ErrorHandling::ErrorCodes::NOT_RECOGNIZED;
        errorLogger->getErrLog(errCode, "Unexpected error loading constants: " + std::string(e.what()));
    }
}

// Navigate through nested JSON using dot notation
json ConfigManager::navigateJsonPath(const json& jsonObj, const std::string& path) const {
    if (path.empty()) return jsonObj;
    
    json current = jsonObj;
    std::string remainingPath = path;
    
    while (!remainingPath.empty()) {
        size_t dotPos = remainingPath.find('.');
        std::string key = (dotPos == std::string::npos) ? remainingPath : remainingPath.substr(0, dotPos);
        
        if (!current.contains(key)) {
            ErrorHandling::ErrorCodes errCode = ErrorHandling::ErrorCodes::RANGE;
            errorLogger->getErrLog(errCode, "JSON path not found: " + path + " (missing key: " + key + ")");
            return json();
        }
        
        current = current[key];
        remainingPath = (dotPos == std::string::npos) ? "" : remainingPath.substr(dotPos + 1);
    }
    
    return current;
}

// Set nested JSON value using dot notation - adds to existing tree structure
void ConfigManager::setNestedJsonValue(json& jsonObj, const std::string& path, const json& value) {
    if (path.empty()) return;
    
    json* current = &jsonObj;
    std::string remainingPath = path;
    
    while (!remainingPath.empty()) {
        size_t dotPos = remainingPath.find('.');
        std::string key = (dotPos == std::string::npos) ? remainingPath : remainingPath.substr(0, dotPos);
        
        if (dotPos == std::string::npos) {
            // Last key - set the value (adds to existing structure without overwriting)
            if (current->contains(key)) {
                ErrorHandling::InfoCodes infoCode = ErrorHandling::InfoCodes::FUNC_EXECUTED;
                errorLogger->getLog(infoCode, "Updated existing property in tree: " + path);
            } else {
                ErrorHandling::InfoCodes infoCode = ErrorHandling::InfoCodes::FUNC_EXECUTED;
                errorLogger->getLog(infoCode, "Added new property to tree: " + path);
            }
            (*current)[key] = value;
            break;
        } else {
            // Navigate deeper, create nested objects if they don't exist
            if (!current->contains(key)) {
                (*current)[key] = json::object();
                ErrorHandling::InfoCodes infoCode = ErrorHandling::InfoCodes::FUNC_EXECUTED;
                errorLogger->getLog(infoCode, "Created new branch in config tree: " + key);
            } else if (!(*current)[key].is_object()) {
                // If existing value is not an object, we need to decide what to do
                // For safety, we'll create a new object and preserve the old value under "_value" key
                json oldValue = (*current)[key];
                (*current)[key] = json::object();
                (*current)[key]["_value"] = oldValue;
                ErrorHandling::InfoCodes infoCode = ErrorHandling::InfoCodes::FUNC_EXECUTED;
                errorLogger->getLog(infoCode, "Converted existing value to object for nesting, preserved old value: " + key);
            }
            current = &((*current)[key]);
            remainingPath = remainingPath.substr(dotPos + 1);
        }
    }
}

// Get singleton instance with thread-safe lazy initialization
ConfigManager& ConfigManager::getInstance() {
    std::lock_guard<std::mutex> lock(mutex_);
    if (instance == nullptr) {
        instance = new ConfigManager();
    }
    return *instance;
}

// Cleanup singleton instance (optional, for explicit cleanup)
void ConfigManager::cleanup() {
    std::lock_guard<std::mutex> lock(mutex_);
    if (instance != nullptr) {
        delete instance;
        instance = nullptr;
    }
}

// Check if property exists
bool ConfigManager::hasProperty(const std::string& path) const {
    if (!propertiesLoaded) return false;
    json value = navigateJsonPath(properties, path);
    return !value.is_null();
}

// Check if constant exists
bool ConfigManager::hasConstant(const std::string& path) const {
    if (!constantsLoaded) return false;
    return constants.contains("values") && constants["values"].contains(path);
}

// Get property as TString (ROOT-specific convenience method)
TString ConfigManager::getPropertyTString(const std::string& path, const TString& defaultValue) const {
    std::string value = getProperty<std::string>(path, std::string(defaultValue.Data()));
    return TString(value.c_str());
}

// Save properties to the original file
bool ConfigManager::saveProperties() const {
    if (!propertiesLoaded || propertiesFilePath.empty()) {
        ErrorHandling::ErrorCodes errCode = ErrorHandling::ErrorCodes::FILE_NOT_EXIST;
        errorLogger->getErrLog(errCode, "Cannot save properties: not loaded or file path empty");
        return false;
    }
    
    return savePropertiesToFile(propertiesFilePath);
}

// Save properties to specific file
bool ConfigManager::savePropertiesToFile(const std::string& filePath) const {
    if (!propertiesLoaded) {
        ErrorHandling::ErrorCodes errCode = ErrorHandling::ErrorCodes::FILE_NOT_EXIST;
        errorLogger->getErrLog(errCode, "Cannot save properties: not loaded");
        return false;
    }
    
    try {
        std::ofstream propertyFile(filePath);
        if (!propertyFile.is_open()) {
            ErrorHandling::ErrorCodes errCode = ErrorHandling::ErrorCodes::FILE_NOT_EXIST;
            errorLogger->getErrLog(errCode, "Cannot open properties file for writing: " + filePath);
            return false;
        }
        
        // Write JSON with nice formatting (indentation = 2)
        propertyFile << properties.dump(4) << std::endl;
        propertyFile.close();
        
        ErrorHandling::InfoCodes infoCode = ErrorHandling::InfoCodes::FUNC_EXECUTED;
        errorLogger->getLog(infoCode, "Properties saved successfully to: " + filePath);
        
        return true;
        
    } catch (const std::exception& e) {
        ErrorHandling::ErrorCodes errCode = ErrorHandling::ErrorCodes::NOT_RECOGNIZED;
        errorLogger->getErrLog(errCode, "Error saving properties to " + filePath + ": " + std::string(e.what()));
        return false;
    }
}

// Reload configurations from files
bool ConfigManager::reloadConfigurations() {
    std::lock_guard<std::mutex> lock(mutex_);
    propertiesLoaded = false;
    constantsLoaded = false;
    properties.clear();
    constants.clear();
    
    loadConfigurations();
    
    bool success = propertiesLoaded && constantsLoaded;
    if (success) {
        ErrorHandling::InfoCodes infoCode = ErrorHandling::InfoCodes::FUNC_EXECUTED;
        errorLogger->getLog(infoCode, "Configuration files reloaded successfully");
    }
    
    return success;
}

// Check if configurations are loaded successfully
bool ConfigManager::isConfigurationLoaded() const {
    return propertiesLoaded && constantsLoaded;
}

// Get access to error logger for external use
ErrorHandling::ErrorLogs& ConfigManager::getErrorLogger() const {
    return *errorLogger;
}
