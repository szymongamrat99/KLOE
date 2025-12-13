#include "FileManager.h"
#include <boost/filesystem.hpp>
#include <fstream>
#include <json.hpp>
#include <regex>
#include <set>
#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include <ctime>
#include <const.h>
#include <TFile.h>
#include <TTree.h>

using json = nlohmann::json;

namespace KLOE
{

  double FileManager::EventsToLuminosity(Long64_t nEvents)
  {
    // Wzór do przeliczenia liczby zdarzeń na luminozność
    //
    // Dla KLOE:
    // - Typowy plik z danymi ma ~N wydarzeń
    // - Odpowiada to ~X nb^-1 luminozności
    //
    // Wzór empiryczny bazujący na liczbie zdarzeń:
    // Luminosity [nb^-1] = nEvents * conversion_factor + displacement

    const double conversion_factor = 0.000908; // [nb^-1 / event]
    const double displacement = 1.38;          // [nb^-1]
    return nEvents * conversion_factor + displacement;
  }

  void FileManager::LogChainLuminosity(TChain &chain, const std::string &logFile)
  {
    std::ofstream log(logFile, std::ios::out | std::ios::trunc);

    if (!log.is_open())
    {
      std::cerr << "ERROR: Cannot open log file: " << logFile << std::endl;
      return;
    }

    // Nagłówek
    auto t = std::time(nullptr);
    auto tm = *std::localtime(&t);
    log << "=== Input Files Luminosity Log ===" << std::endl;
    log << "Timestamp: " << std::put_time(&tm, "%Y-%m-%d %H:%M:%S") << std::endl;
    log << "Chain name: " << chain.GetName() << std::endl;
    log << "Total entries: " << chain.GetEntries() << std::endl;
    log << "==================================\n"
        << std::endl;

    double totalLuminosity = 0.0;
    Long64_t totalEvents = 0;

    // Iteruj po plikach w TChain
    TObjArray *fileElements = chain.GetListOfFiles();
    TIter next(fileElements);
    TChainElement *chEl = nullptr;

    log << std::left << std::setw(50) << "File"
        << std::right << std::setw(15) << "Events"
        << std::setw(15) << "Lumi [nb^-1]" << std::endl;
    log << std::string(80, '-') << std::endl;

    while ((chEl = (TChainElement *)next()))
    {
      std::string filename = chEl->GetTitle();
      boost::filesystem::path p(filename);

      // Pobierz liczbę zdarzeń z tego pliku
      Long64_t entries = chEl->GetEntries();

      // Przelicz liczbę zdarzeń na luminozność
      double lumi = EventsToLuminosity(entries);

      // Loguj
      log << std::left << std::setw(50) << p.filename().string()
          << std::right << std::fixed
          << std::setw(15) << entries
          << std::setw(15) << std::setprecision(4) << lumi << std::endl;

      totalLuminosity += lumi;
      totalEvents += entries;
    }

    log << std::string(80, '=') << std::endl;
    log << std::left << std::setw(50) << "TOTAL"
        << std::right << std::fixed
        << std::setw(15) << totalEvents
        << std::setw(15) << std::setprecision(4) << totalLuminosity << std::endl;

    log.close();

    std::cout << "Input files luminosity log saved to: " << logFile << std::endl;
    std::cout << "Total integrated luminosity: " << std::setprecision(4)
              << totalLuminosity << " nb^-1 (" << totalEvents << " events)" << std::endl;
  }

  RunStats FileManager::getRunStats(const std::string &directory, const std::string &regex_pattern)
  {
    std::set<int> runs;
    std::regex run_regex(regex_pattern);

    // Tylko przeskanuj nazwy plików - NIE otwieraj plików ROOT
    // To dramatycznie przyspiesza działanie
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
    stats.totalEvents = 0; // Będzie policzone przez TChain później
    stats.runList.assign(runs.begin(), runs.end());
    stats.totalLuminosity = 0.0; // Będzie obliczone później
    return stats;
  }

  void FileManager::UpdateRunStatsFromChain(RunStats &stats, TChain &chain)
  {
    // Oblicz całkowitą liczbę zdarzeń z TChain (już załadowanego)
    stats.totalEvents = chain.GetEntries();

    // Oblicz luminozność
    stats.totalLuminosity = EventsToLuminosity(stats.totalEvents);
  }

  void FileManager::chainInit(TChain &chain_init, Controls::DataType &dataTypeOpt, UInt_t &firstData, UInt_t &lastData, UInt_t &firstMC, UInt_t &lastMC, ErrorHandling::ErrorLogs &logger, Int_t csFlag, Int_t oldAnaFlag)
  {
    ErrorHandling::InfoCodes infoCode;

    std::ifstream rootFiles(Paths::rootfilesName);
    json filePaths = json::parse(rootFiles);

    Bool_t *all_phys;

    TString fullname = "",
            dirnamedata,
            filenamedata,
            extension = Paths::ext_root,
            newestDateStampData = (std::string)Utils::properties["variables"]["rootFiles"]["newestFileDateStampData"],
            newestDateStampMC = (std::string)Utils::properties["variables"]["rootFiles"]["newestFileDateStampMC"],
            path = (std::string)Utils::properties["variables"]["rootFiles"]["pathOldAna"],
            mc_dir = "MONTE_CARLO",
            data_dir = "DATA";

    Int_t fileNumData[2] = {Utils::properties["variables"]["rootFiles"]["Data"]["firstFile"],
                            Utils::properties["variables"]["rootFiles"]["Data"]["lastFile"]};

    std::vector<Int_t> fileNumMC[3];

    std::vector<TString>
        dirnamemc,
        filenamemc;

    // Przygotuj nazwy plików MC - nowa analiza
    for (Int_t i = 0; i < 3; i++)
    {
      dirnamemc.push_back(path + mc_dir);
      // filenamemc.push_back((std::string)filePaths["MC"]["filenameBase"][i]);
      // fileNumMC[i].push_back(Utils::properties["variables"]["rootFiles"]["MC"][i]["firstFile"]);
      // fileNumMC[i].push_back(Utils::properties["variables"]["rootFiles"]["MC"][i]["lastFile"]);
    }
    //////////////////////////////////////////////
    // Przygotuj nazwy plików dla starej analizy
    TString dirnamemcOld = path + mc_dir;
    TString filenamemcOld = (std::string)filePaths["MC"]["filenameBaseOldAna"];
    Int_t fileNumMCOld[2] = {Utils::properties["variables"]["rootFiles"]["firstFileOld"],
                            Utils::properties["variables"]["rootFiles"]["lastFileOld"]};
    //////////////////////////////////////////////

    dirnamedata = path + data_dir;
    filenamedata = (std::string)filePaths["Data"]["filenameBaseOldAna"];

    switch (dataTypeOpt)
    {
    case Controls::DataType::MC_DATA:
    {
      for (Int_t i = firstData; i <= lastData; i++)
      {
        fullname = dirnamedata + "/" + filenamedata + "_" + std::to_string(i) + extension;

        // Check if file exists
        boost::filesystem::path pathExist(fullname);

        if (1) // boost::filesystem::exists(pathExist))
        {
          infoCode = ErrorHandling::InfoCodes::FILE_ADDED;

          chain_init.Add(fullname);
          logger.getLog(infoCode, (std::string)filenamedata + "_" + std::to_string(i) + (std::string)extension);
        }
      }

      if (oldAnaFlag == false)
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
      }
      else
      {


        for (Int_t j = fileNumMCOld[0]; j <= fileNumMCOld[1]; j++)
        {

          fullname = dirnamemcOld + "/" + filenamemcOld + "_" + std::to_string(j) + extension;

          // Check if file exists
          boost::filesystem::path pathExist(fullname);

          if (boost::filesystem::exists(pathExist))
          {
            infoCode = ErrorHandling::InfoCodes::FILE_ADDED;

            chain_init.Add(fullname);
            logger.getLog(infoCode, (std::string)filenamemcOld + "_" + std::to_string(j) + (std::string)extension);

            std::cout << "Added MC file: " << fullname << std::endl;
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

    std::ifstream rootFiles(Paths::rootfilesName);
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

        // Priorytet dla v2.root i v1.root
        if (fname.find("_v2.root") != std::string::npos)
        {
          runFiles[runNum] = it->path().string();
        }
        else
        {
          // Dodaj tylko jeśli nie ma już v2 dla tego runu
          if (runFiles.count(runNum) == 0)
            runFiles[runNum] = it->path().string();
        }
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
        else
        {
          if (runFiles.count(runNum) == 0)
            runFiles[runNum] = filePath;
        }
      }
    }

    // Dodaj wybrane pliki do TChain
    for (const auto &kv : runFiles)
    {
      chain.Add(kv.second.c_str());
    }
  }

  bool FileManager::ValidateJobListFilename(const std::string &filename)
  {
    // Format: job_v{wersja}_{typ}_{luminosity}_inv_pb_{numer}.txt
    // Przykład: job_v1_data_5000_inv_pb_001.txt

    std::regex jobFileRegex(R"(^job_v\d+_(data|all_phys|all_phys2|all_phys3)_[\d.]+_inv_pb_\d+\.txt$)");
    return std::regex_match(filename, jobFileRegex);
  }

  std::vector<std::string> FileManager::LoadFileListFromFile(const std::string &filePath)
  {
    std::vector<std::string> fileList;
    std::ifstream file(filePath);

    if (!file.is_open())
    {
      throw std::runtime_error("Cannot open file: " + filePath);
    }

    std::string line;
    while (std::getline(file, line))
    {
      // Usuń białe znaki z początku i końca linii
      line.erase(0, line.find_first_not_of(" \t\r\n"));
      line.erase(line.find_last_not_of(" \t\r\n") + 1);

      // Pomiń puste linie i komentarze
      if (line.empty() || line[0] == '#')
        continue;

      fileList.push_back(line);
    }

    file.close();

    if (fileList.empty())
    {
      throw std::runtime_error("No files found in job list file: " + filePath);
    }

    return fileList;
  }

  void FileManager::chainInit(TChain &chain, ErrorHandling::ErrorLogs &logger,
                              const std::vector<std::string> &fileList)
  {
    ErrorHandling::ErrorCodes errorCode;

    for (const auto &filepath : fileList)
    {
      // Sprawdź czy plik istnieje
      boost::filesystem::path pathObj(filepath);
      if (!boost::filesystem::exists(pathObj))
      {
        std::cerr << "WARNING: File does not exist: " << filepath << std::endl;
        errorCode = ErrorHandling::ErrorCodes::FILE_NOT_EXIST;
        std::cout << "Dupsko" << std::endl;
        logger.getErrLog(errorCode, "File not found: " + filepath);
        continue;
      }

      // Sprawdź czy to plik ROOT
      if (pathObj.extension() != ".root")
      {
        std::cerr << "WARNING: File is not ROOT file: " << filepath << std::endl;
        continue;
      }

      // Dodaj plik do TChain
      chain.Add(filepath.c_str());
      errorCode = ErrorHandling::ErrorCodes::FILE_NOT_EXIST;
      logger.getErrLog(errorCode, filepath);
    }
  }
} // namespace KLOE
