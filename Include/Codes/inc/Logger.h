#pragma once

#include <ErrorLogs.h>
#include <Logs.h>

#include <const.h>

enum class LogFile
{
  GENERAL,
  ANALYSIS,
  EFFICIENCY,
  NEUTREC,
  CPFIT
};

using json = nlohmann::json;

class Logger
{
  public:
    static Logger& getInstance()
    {
        static Logger instance;
        return instance;
    } 

    void logToFile(LogFile fileType, const std::string& message)
    {
      if (logFileMap.find(fileType) == logFileMap.end())
      {
        logFileMap[fileType] = std::make_unique<std::ofstream>();
        std::string filename;

      }
    }

  private:
    Logger() {}
    ~Logger() {}

    // Zapobiegamy kopiowaniu
    Logger(const Logger&) = delete;
    void operator=(const Logger&) = delete;

    void CollectLogFilePaths()
    {
      json pathsFile;
      std::ifstream fpath(Paths::pathsExtensionsPath);
      fpath >> pathsFile;

      if (pathsFile.contains("logFilesBaseNames"))
      {
        auto logFilesBaseNames = pathsFile["logFilesBaseNames"];
        for (const auto& item : logFilesBaseNames.items())
        {
          
        }
      }
    }

    std::map<LogFile, std::unique_ptr<std::ofstream>> logFileMap;
};


