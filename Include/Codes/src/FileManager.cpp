#include "FileManager.h"
#include <boost/filesystem.hpp>
#include <fstream>
#include <json.hpp>
#include <regex>
#include <set>
#include <vector>
#include <string>
#include <iostream>
#include <const.h>


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
            newestDateStamp = (std::string)properties["variables"]["rootFiles"]["newestFileDateStamp"],
            path = (std::string)properties["variables"]["rootFiles"]["path"];

    Int_t fileNumData[2] = {properties["variables"]["rootFiles"]["Data"]["firstFile"],
                            properties["variables"]["rootFiles"]["Data"]["lastFile"]};

    std::vector<Int_t> fileNumMC[3];

    std::vector<TString>
        dirnamemc,
        filenamemc;

    for (Int_t i = 0; i < 3; i++)
    {

      dirnamemc.push_back(path + newestDateStamp);
      filenamemc.push_back((std::string)filePaths["MC"]["filenameBase"][i]);
      fileNumMC[i].push_back(properties["variables"]["rootFiles"]["MC"][i]["firstFile"]);
      fileNumMC[i].push_back(properties["variables"]["rootFiles"]["MC"][i]["lastFile"]);
    }

    dirnamedata = path + newestDateStamp;
    filenamedata = (std::string)filePaths["Data"]["filenameBase"];

    switch (dataTypeOpt)
    {
    case Controls::DataType::MC_DATA:
    {
      for (Int_t i = firstData; i <= lastData; i++)
      {
        fullname = dirnamedata + "/" + filenamedata + "_" + std::to_string(i) + extension;

        // Check if file exists
        boost::filesystem::path pathExist(fullname);

        if (boost::filesystem::exists(pathExist))
        {
          infoCode = ErrorHandling::InfoCodes::FILE_ADDED;

          chain_init.Add(fullname);
          logger.getLog(infoCode, (std::string)filenamedata + "_" + std::to_string(i) + (std::string)extension);
        }
      }

      for (Int_t k = 0; k < dirnamemc.size(); k++)
      {
        for (Int_t j = fileNumMC[k][0]; j <= fileNumMC[k][1]; j++)
        {

          fullname = dirnamemc[k] + "/" + filenamemc[k] + "_" + std::to_string(j) + extension;

          // Check if file exists
          boost::filesystem::path pathExist(fullname);

          if (boost::filesystem::exists(pathExist))
          {
            infoCode = ErrorHandling::InfoCodes::FILE_ADDED;

            chain_init.Add(fullname);
            logger.getLog(infoCode, (std::string)filenamemc[k] + "_" + std::to_string(j) + (std::string)extension);
          }
        }
      }

      break;
    }
    case Controls::DataType::MC_ONLY:
    {
      for (Int_t k = 0; k < dirnamemc.size(); k++)
      {
        for (Int_t j = fileNumMC[k][0]; j <= fileNumMC[k][1]; j++)
        {

          fullname = dirnamemc[k] + "/" + filenamemc[k] + "_" + std::to_string(j) + extension;

          // Check if file exists
          boost::filesystem::path pathExist(fullname);

          if (boost::filesystem::exists(pathExist))
          {
            infoCode = ErrorHandling::InfoCodes::FILE_ADDED;

            chain_init.Add(fullname);
            logger.getLog(infoCode, (std::string)filenamemc[k] + "_" + std::to_string(j) + (std::string)extension);
          }
        }
      }

      break;
    }
    case Controls::DataType::DATA_ONLY:
    {
      for (Int_t i = firstData; i <= lastData; i++)
      {
        fullname = dirnamedata + "/" + filenamedata + "_" + std::to_string(i) + extension;

        // Check if file exists
        boost::filesystem::path pathExist(fullname);

        if (boost::filesystem::exists(pathExist))
        {
          infoCode = ErrorHandling::InfoCodes::FILE_ADDED;

          chain_init.Add(fullname);
          logger.getLog(infoCode, (std::string)filenamedata + "_" + std::to_string(i) + (std::string)extension);
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

  void FileManager::chainInit(TChain &chain, ErrorHandling::ErrorLogs &logger,
                              const std::string &DataPath, const std::string &runRegexPattern,
                              int minRun, int maxRun)
  {
    std::regex runRegex(runRegexPattern);
    std::map<int, std::string> runFiles;

    for (boost::filesystem::directory_iterator it(DataPath), end; it != end; ++it)
    {
      if (!boost::filesystem::is_regular_file(it->status()))
        continue;
      const std::string fname = it->path().filename().string();
      std::smatch match;
      if (std::regex_match(fname, match, runRegex))
      {
        int runNum = std::stoi(match[1]);
        if (runNum < minRun || runNum > maxRun)
          continue;

        // Priorytet dla v2.root
        if (fname.find("_v2.root") != std::string::npos)
        {
          runFiles[runNum] = it->path().string();
        }
        // } else {
        //     // Dodaj tylko jeśli nie ma już v2 dla tego runu
        //     if (runFiles.count(runNum) == 0)
        //         runFiles[runNum] = it->path().string();
        // }
      }
    }

    // Dodaj wybrane pliki do TChain
    for (const auto &kv : runFiles)
    {
      chain.Add(kv.second.c_str());
    }
  }

  void FileManager::chainInit(TChain &chain, ErrorHandling::ErrorLogs &logger,
                              const std::vector<std::string> &fileList, const std::string &runRegexPattern,
                              int minRun, int maxRun)
  {
    std::regex runRegex(runRegexPattern);
    std::map<int, std::string> runFiles;

    for (const auto &filePath : fileList)
    {
      boost::filesystem::path pathObj(filePath);
      if (!boost::filesystem::is_regular_file(pathObj))
        continue;
      const std::string fname = pathObj.filename().string();
      std::smatch match;
      if (std::regex_match(fname, match, runRegex))
      {
        int runNum = std::stoi(match[1]);
        if (runNum < minRun || runNum > maxRun)
          continue;

        // Priorytet dla v2.root
        if (fname.find("_v2.root") != std::string::npos)
        {
          runFiles[runNum] = filePath;
        }
        // } else {
        //     if (runFiles.count(runNum) == 0)
        //         runFiles[runNum] = filePath;
        // }
      }
    }

    // Dodaj wybrane pliki do TChain
    for (const auto &kv : runFiles)
    {
      chain.Add(kv.second.c_str());
    }
  }
} // namespace KLOE
