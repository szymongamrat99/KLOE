#ifndef CONFIG_MANAGER_H
#define CONFIG_MANAGER_H

#include <json.hpp>
#include <fstream>
#include <string>
#include <mutex>
#include <TString.h>
#include <ErrorLogs.h>

using json = nlohmann::json;

/**
 * @class ConfigManager
 * @brief Thread-safe singleton for managing properties and constants configuration
 * 
 * This class provides centralized access to JSON configuration files (properties.json and pdg_const.json)
 * with proper error handling using the ErrorHandling system and thread-safe lazy initialization.
 */
class ConfigManager {
private:
    static ConfigManager* instance;
    static std::mutex mutex_;
    
    json properties;
    json constants;
    bool propertiesLoaded;
    bool constantsLoaded;
    
    // Error logging system
    mutable ErrorHandling::ErrorLogs* errorLogger;
    
    /**
     * @brief Private constructor implementing lazy loading
     */
    ConfigManager() : propertiesLoaded(false), constantsLoaded(false), errorLogger(nullptr) {
        // Initialize error logger with config-specific log file
        std::string logFileName = "log/config_manager.log";
        errorLogger = new ErrorHandling::ErrorLogs(logFileName);
        errorLogger->setPrintToScreen(true);
        
        loadConfigurations();
    }
    
    /**
     * @brief Load both properties and constants JSON files
     */
    void loadConfigurations() {
        loadProperties();
        loadConstants();
    }
    
    /**
     * @brief Load properties.json file
     */
    void loadProperties() {
        try {
            const char* propertiesPath = getenv("PROPERTIESKLOE");
            if (!propertiesPath) {
                ErrorHandling::ErrorCodes errCode = ErrorHandling::ErrorCodes::FILE_NOT_EXIST;
                errorLogger->getErrLog(errCode, "PROPERTIESKLOE environment variable not set");
                return;
            }
            
            std::string propName = std::string(propertiesPath) + "/properties.json";
            std::ifstream propertyFile(propName.c_str());
            
            if (!propertyFile.is_open()) {
                ErrorHandling::ErrorCodes errCode = ErrorHandling::ErrorCodes::FILE_NOT_EXIST;
                errorLogger->getErrLog(errCode, "Cannot open properties file: " + propName);
                return;
            }
            
            properties = json::parse(propertyFile);
            propertiesLoaded = true;
            
            ErrorHandling::InfoCodes infoCode = ErrorHandling::InfoCodes::FILE_ADDED;
            errorLogger->getLog(infoCode, "Properties loaded successfully from: " + propName);
            
        } catch (const json::parse_error& e) {
            ErrorHandling::ErrorCodes errCode = ErrorHandling::ErrorCodes::DATA_TYPE;
            errorLogger->getErrLog(errCode, "JSON parse error in properties: " + std::string(e.what()));
        } catch (const std::exception& e) {
            ErrorHandling::ErrorCodes errCode = ErrorHandling::ErrorCodes::NOT_RECOGNIZED;
            errorLogger->getErrLog(errCode, "Unexpected error loading properties: " + std::string(e.what()));
        }
    }
    
    /**
     * @brief Load pdg_const.json file
     */
    void loadConstants() {
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
    
    /**
     * @brief Navigate through nested JSON using dot notation
     * @param jsonObj JSON object to navigate
     * @param path Dot-separated path (e.g., "variables.rootFiles.firstFileMax")
     * @return JSON value at the specified path
     */
    json navigateJsonPath(const json& jsonObj, const std::string& path) const {
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
    
public:
    // Delete copy constructor and assignment operator (C++14 compatible)
    ConfigManager(const ConfigManager&) = delete;
    ConfigManager& operator=(const ConfigManager&) = delete;
    
    /**
     * @brief Destructor - cleanup error logger
     */
    ~ConfigManager() {
        if (errorLogger != nullptr) {
            delete errorLogger;
            errorLogger = nullptr;
        }
    }
    
    /**
     * @brief Get singleton instance with thread-safe lazy initialization
     * @return Reference to the singleton instance
     */
    static ConfigManager& getInstance() {
        std::lock_guard<std::mutex> lock(mutex_);
        if (instance == nullptr) {
            instance = new ConfigManager();
        }
        return *instance;
    }
    
    /**
     * @brief Cleanup singleton instance (optional, for explicit cleanup)
     */
    static void cleanup() {
        std::lock_guard<std::mutex> lock(mutex_);
        if (instance != nullptr) {
            delete instance;
            instance = nullptr;
        }
    }
    
    /**
     * @brief Get property value with type conversion and error handling
     * @tparam T Type to convert the value to
     * @param path Dot-separated path to the property (e.g., "variables.rootFiles.firstFileMax")
     * @param defaultValue Default value to return if property not found
     * @return Property value or default value
     */
    template<typename T>
    T getProperty(const std::string& path, const T& defaultValue) const {
        if (!propertiesLoaded) {
            ErrorHandling::ErrorCodes errCode = ErrorHandling::ErrorCodes::FILE_NOT_EXIST;
            errorLogger->getErrLog(errCode, "Properties not loaded, using default for: " + path);
            return defaultValue;
        }
        
        try {
            json value = navigateJsonPath(properties, path);
            if (value.is_null()) {
                return defaultValue;
            }
            return value.get<T>();
        } catch (const json::type_error& e) {
            ErrorHandling::ErrorCodes errCode = ErrorHandling::ErrorCodes::DATA_TYPE;
            errorLogger->getErrLog(errCode, "Type conversion error for property " + path + ": " + std::string(e.what()));
            return defaultValue;
        } catch (const std::exception& e) {
            ErrorHandling::ErrorCodes errCode = ErrorHandling::ErrorCodes::NOT_RECOGNIZED;
            errorLogger->getErrLog(errCode, "Error getting property " + path + ": " + std::string(e.what()));
            return defaultValue;
        }
    }
    
    /**
     * @brief Get property value with default-constructed default value (C++14 compatible)
     * @tparam T Type to convert the value to
     * @param path Dot-separated path to the property
     * @return Property value or default-constructed T
     */
    template<typename T>
    T getProperty(const std::string& path) const {
        return getProperty<T>(path, T());
    }
    
    /**
     * @brief Get constant value with type conversion and error handling
     * @tparam T Type to convert the value to
     * @param path PDG path to the constant (e.g., "/S011M" for K0 mass)
     * @param defaultValue Default value to return if constant not found
     * @return Constant value or default value
     */
    template<typename T>
    T getConstant(const std::string& path, const T& defaultValue) const {
        if (!constantsLoaded) {
            ErrorHandling::ErrorCodes errCode = ErrorHandling::ErrorCodes::FILE_NOT_EXIST;
            errorLogger->getErrLog(errCode, "Constants not loaded, using default for: " + path);
            return defaultValue;
        }
        
        try {
            if (!constants.contains("values") || !constants["values"].contains(path)) {
                ErrorHandling::ErrorCodes errCode = ErrorHandling::ErrorCodes::RANGE;
                errorLogger->getErrLog(errCode, "Constant not found: " + path);
                return defaultValue;
            }
            
            return constants["values"][path].get<T>();
        } catch (const json::type_error& e) {
            ErrorHandling::ErrorCodes errCode = ErrorHandling::ErrorCodes::DATA_TYPE;
            errorLogger->getErrLog(errCode, "Type conversion error for constant " + path + ": " + std::string(e.what()));
            return defaultValue;
        } catch (const std::exception& e) {
            ErrorHandling::ErrorCodes errCode = ErrorHandling::ErrorCodes::NOT_RECOGNIZED;
            errorLogger->getErrLog(errCode, "Error getting constant " + path + ": " + std::string(e.what()));
            return defaultValue;
        }
    }
    
    /**
     * @brief Get constant value with default-constructed default value (C++14 compatible)
     * @tparam T Type to convert the value to
     * @param path PDG path to the constant
     * @return Constant value or default-constructed T
     */
    template<typename T>
    T getConstant(const std::string& path) const {
        return getConstant<T>(path, T());
    }
    
    /**
     * @brief Check if property exists
     * @param path Dot-separated path to the property
     * @return True if property exists, false otherwise
     */
    bool hasProperty(const std::string& path) const {
        if (!propertiesLoaded) return false;
        json value = navigateJsonPath(properties, path);
        return !value.is_null();
    }
    
    /**
     * @brief Check if constant exists
     * @param path PDG path to the constant
     * @return True if constant exists, false otherwise
     */
    bool hasConstant(const std::string& path) const {
        if (!constantsLoaded) return false;
        return constants.contains("values") && constants["values"].contains(path);
    }
    
    /**
     * @brief Get property as TString (ROOT-specific convenience method)
     * @param path Dot-separated path to the property
     * @param defaultValue Default value to return if property not found
     * @return Property value as TString or default value
     */
    TString getPropertyTString(const std::string& path, const TString& defaultValue = "") const {
        std::string value = getProperty<std::string>(path, std::string(defaultValue.Data()));
        return TString(value.c_str());
    }
    
    /**
     * @brief Reload configurations from files
     * @return True if reload was successful, false otherwise
     */
    bool reloadConfigurations() {
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
    
    /**
     * @brief Check if configurations are loaded successfully
     * @return True if both properties and constants are loaded
     */
    bool isConfigurationLoaded() const {
        return propertiesLoaded && constantsLoaded;
    }
    
    /**
     * @brief Get access to error logger for external use
     * @return Reference to the error logger
     */
    ErrorHandling::ErrorLogs& getErrorLogger() const {
        return *errorLogger;
    }
};

// Static member definitions
ConfigManager* ConfigManager::instance = nullptr;
std::mutex ConfigManager::mutex_;

#endif // CONFIG_MANAGER_H
