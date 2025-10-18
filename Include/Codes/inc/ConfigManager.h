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
    std::string propertiesFilePath; // Store path for saving
    
    // Error logging system
    mutable ErrorHandling::ErrorLogs* errorLogger;
    
    /**
     * @brief Private constructor implementing lazy loading
     */
    ConfigManager();
    
    /**
     * @brief Load both properties and constants JSON files
     */
    void loadConfigurations();
    
    /**
     * @brief Load properties.json file
     */
    void loadProperties();
    
    /**
     * @brief Load pdg_const.json file
     */
    void loadConstants();
    
    /**
     * @brief Navigate through nested JSON using dot notation
     * @param jsonObj JSON object to navigate
     * @param path Dot-separated path (e.g., "variables.rootFiles.firstFileMax")
     * @return JSON value at the specified path
     */
    json navigateJsonPath(const json& jsonObj, const std::string& path) const;
    
    /**
     * @brief Set nested JSON value using dot notation
     * @param jsonObj JSON object to modify
     * @param path Dot-separated path (e.g., "variables.rootFiles.firstFileMax")
     * @param value Value to set
     */
    void setNestedJsonValue(json& jsonObj, const std::string& path, const json& value);
    
public:
    // Delete copy constructor and assignment operator (C++14 compatible)
    ConfigManager(const ConfigManager&) = delete;
    ConfigManager& operator=(const ConfigManager&) = delete;
    
    /**
     * @brief Destructor - cleanup error logger
     */
    ~ConfigManager();
    
    /**
     * @brief Get singleton instance with thread-safe lazy initialization
     * @return Reference to the singleton instance
     */
    static ConfigManager& getInstance();
    
    /**
     * @brief Cleanup singleton instance (optional, for explicit cleanup)
     */
    static void cleanup();
    
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
    bool hasProperty(const std::string& path) const;
    
    /**
     * @brief Check if constant exists
     * @param path PDG path to the constant
     * @return True if constant exists, false otherwise
     */
    bool hasConstant(const std::string& path) const;
    
    /**
     * @brief Set property value in the configuration tree (adds to existing structure)
     * @tparam T Type of the value to set
     * @param path Dot-separated path to the property (e.g., "variables.rootFiles.newParam")
     * @param value Value to set at the specified path
     * @return True if property was set successfully, false otherwise
     */
    template<typename T>
    bool setProperty(const std::string& path, const T& value) {
        if (!propertiesLoaded) {
            ErrorHandling::ErrorCodes errCode = ErrorHandling::ErrorCodes::FILE_NOT_EXIST;
            errorLogger->getErrLog(errCode, "Cannot set property, properties not loaded: " + path);
            return false;
        }
        
        try {
            json jsonValue = value; // Convert to JSON
            setNestedJsonValue(properties, path, jsonValue);
            
            ErrorHandling::InfoCodes infoCode = ErrorHandling::InfoCodes::FUNC_EXECUTED;
            errorLogger->getLog(infoCode, "Property set successfully: " + path);
            
            return true;
        } catch (const std::exception& e) {
            ErrorHandling::ErrorCodes errCode = ErrorHandling::ErrorCodes::NOT_RECOGNIZED;
            errorLogger->getErrLog(errCode, "Error setting property " + path + ": " + std::string(e.what()));
            return false;
        }
    }
    
    /**
     * @brief Save current properties to the original file
     * @return True if save was successful, false otherwise
     */
    bool saveProperties() const;
    
    /**
     * @brief Save current properties to a specific file
     * @param filePath Path to the file where properties should be saved
     * @return True if save was successful, false otherwise
     */
    bool savePropertiesToFile(const std::string& filePath) const;
    
    /**
     * @brief Get property as TString (ROOT-specific convenience method)
     * @param path Dot-separated path to the property
     * @param defaultValue Default value to return if property not found
     * @return Property value as TString or default value
     */
    TString getPropertyTString(const std::string& path, const TString& defaultValue = "") const;
    
    /**
     * @brief Reload configurations from files
     * @return True if reload was successful, false otherwise
     */
    bool reloadConfigurations();
    
    /**
     * @brief Check if configurations are loaded successfully
     * @return True if both properties and constants are loaded
     */
    bool isConfigurationLoaded() const;
    
    /**
     * @brief Get access to error logger for external use
     * @return Reference to the error logger
     */
    ErrorHandling::ErrorLogs& getErrorLogger() const;
};

#endif // CONFIG_MANAGER_H
