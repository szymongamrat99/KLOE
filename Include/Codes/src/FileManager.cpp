#include "FileManager.h"
#include <boost/filesystem.hpp>
#include <fstream>
#include <json.hpp>
#include <regex>
#include <set>
#include <vector>
#include <string>
#include <iostream>

using json = nlohmann::json;

namespace KLOE
{

  RunStats FileManager::getRunStats(const std::string &directory, const std::string &regex_pattern)
  {
    std::set<int> runs;
    size_t totalBytes = 0;
    std::regex run_regex(regex_pattern);

    for (auto &entry : boost::filesystem::directory_iterator(directory))
    {
      if (!boost::filesystem::is_regular_file(entry.status()))
        continue;
      std::string filename = entry.path().filename().string();
      std::smatch match;
      if (std::regex_search(filename, match, run_regex) && match.size() > 1)
      {
        try
        {
          int run = std::stoi(match[1].str());
          runs.insert(run);
          totalBytes += boost::filesystem::file_size(entry.path());
        }
        catch (...)
        {
        }
      }
    }

    RunStats stats;
    if (!runs.empty())
    {
      stats.minRun = *runs.begin();
      stats.maxRun = *runs.rbegin();
    }
    else
    {
      stats.minRun = stats.maxRun = -1;
    }
    stats.fileCount = runs.size();
    stats.totalGB = totalBytes / (1024.0 * 1024.0 * 1024.0);
    stats.runList.assign(runs.begin(), runs.end());
    return stats;
  }

  void FileManager::chainInit(TChain &chain_init, Controls::DataType &dataTypeOpt, UInt_t &firstData, UInt_t &lastData, UInt_t &firstMC, UInt_t &lastMC, ErrorHandling::ErrorLogs &logger, Int_t csFlag)
  {
    ErrorHandling::InfoCodes infoCode;

    std::ifstream rootFiles(rootfilesName);
    json filePaths = json::parse(rootFiles);

    Bool_t *all_phys;

    TString fullname = "",
            dirnamedata,
            filenamedata,
            extension = ext_root,
            path = (std::string)properties["variables"]["rootFiles"]["path"];

    std::vector<TString>
        dirnamemc,
        filenamemc;

    for (Int_t i = 0; i < 3; i++)
    {
      if (all_phys[i])
      {
        dirnamemc.push_back((std::string)properties["variables"]["rootFiles"]["directoryMC"][i]);
        filenamemc.push_back((std::string)properties["variables"]["rootFiles"]["filenameMC"][i]);
      }
    }

    dirnamedata = "DATA";
    filenamedata = "data_stream42_";

    switch (dataTypeOpt)
    {
    case Controls::DataType::MC_DATA:
    {
      for (Int_t i = firstData; i <= lastData; i++)
      {
        fullname = path + "/" + dirnamedata + "/" + filenamedata + std::to_string(i) + extension;

        // Check if file exists
        boost::filesystem::path pathExist(fullname);

        if (boost::filesystem::exists(pathExist))
        {
          infoCode = ErrorHandling::InfoCodes::FILE_ADDED;

          chain_init.Add(fullname);
          logger.getLog(infoCode, (std::string)filenamedata + std::to_string(i) + (std::string)extension);
        }
      }

      for (Int_t j = firstMC; j <= lastMC; j++)
      {
        // fullname = path + "/" + dirnamemc + "/" + filenamemc + std::to_string(j) + extension;

        // Check if file exists
        boost::filesystem::path pathExist(fullname);

        if (boost::filesystem::exists(pathExist))
        {
          infoCode = ErrorHandling::InfoCodes::FILE_ADDED;

          chain_init.Add(fullname);
          // logger.getLog(infoCode, (std::string)filenamemc + std::to_string(j) + (std::string)extension);
        }
      }

      break;
    }
    default:
    {
      for (Int_t j = firstMC; j <= lastMC; j++)
      {
        // fullname = path + "/" + dirnamemc + "/" + filenamemc + std::to_string(j) + extension;

        // Check if file exists
        boost::filesystem::path pathExist(fullname);

        if (boost::filesystem::exists(pathExist))
        {
          infoCode = ErrorHandling::InfoCodes::FILE_ADDED;

          chain_init.Add(fullname);
          // logger.getLog(infoCode, (std::string)filenamemc + std::to_string(j) + (std::string)extension);
        }
      }

      break;
    }
    }
  }

  void FileManager::chainInit(TChain &chain_init, ErrorHandling::ErrorLogs &logger)
  {
    ErrorHandling::InfoCodes infoCode;

    std::ifstream rootFiles(rootfilesName);
    json filePaths = json::parse(rootFiles);

    std::string
        DataPath = filePaths["MC"]["path"][2],
        DataFilenameBase = filePaths["MC"]["filenameBase"][2],
        DataExtension = filePaths["MC"]["extension"];

    std::string fullnameData = DataPath + DataFilenameBase + "41618" + DataExtension;

    // Check if file exists
    boost::filesystem::path pathExistData(fullnameData);
    if (boost::filesystem::exists(pathExistData))
    {
      infoCode = ErrorHandling::InfoCodes::FILE_ADDED;

      chain_init.Add(fullnameData.c_str());
      logger.getLog(infoCode, DataFilenameBase + "41618" + DataExtension);
    }
  }

  void FileManager::chainInit(TChain &chain_init, ErrorHandling::ErrorLogs &logger,
                              const std::string &directory, const std::string &regex_pattern,
                              int minRun, int maxRun)
  {
    ErrorHandling::InfoCodes infoCode;
    std::regex run_regex(regex_pattern);
    for (auto &entry : boost::filesystem::directory_iterator(directory))
    {
      if (!boost::filesystem::is_regular_file(entry.status()))
        continue;
      std::string filename = entry.path().filename().string();
      std::smatch match;
      if (std::regex_search(filename, match, run_regex) && match.size() > 1)
      {
        int run = -1;
        try
        {
          run = std::stoi(match[1].str());
        }
        catch (...)
        {
          continue;
        }
        if (run >= minRun && run <= maxRun)
        {
          std::string fullpath = entry.path().string();
          chain_init.Add(fullpath.c_str());
          infoCode = ErrorHandling::InfoCodes::FILE_ADDED;
          logger.getLog(infoCode, filename);
        }
      }
    }
  }
} // namespace KLOE
