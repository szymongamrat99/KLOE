#include "../Include/klspm00.hpp"

#include "../Include/MainMenuHandler/MainMenuHandler.h"
#include "../Include/MainMenuHandler/InputParamsHandler.h"
#include <FileManager.h>
#include <ConfigManager.h>
#include <event_data.h>
#include <AnalysisManager.h>
#include <TROOT.h>

#include <const.h>
#include <boost/filesystem.hpp>

using namespace std;
using namespace std::chrono;

int main(int argc, char *argv[])
{
  ROOT::EnableImplicitMT(16); // Użyj do 16 wątków

  bool
      firstFileRangeErr,
      lastFileRangeErr,
      dataTypeErr,
      menuRangeErr;

  UInt_t firstFile, lastFile, csFlag;

  // -------------------------------------------------------------------
  // Check if a job list file was provided as argument
  std::string jobListFile = "";
  bool useJobListFile = false;

  if (argc > 1)
  {
    jobListFile = argv[1];

    // Sprawdzenie czy plik istnieje i ma prawidłową nazwę
    boost::filesystem::path jobListPath(jobListFile);
    if (!boost::filesystem::exists(jobListPath))
    {
      std::cerr << "ERROR: Job list file does not exist: " << jobListFile << std::endl;
      return 1;
    }

    // Walidacja nazwy pliku
    if (!KLOE::FileManager::ValidateJobListFilename(jobListPath.filename().string()))
    {
      std::cerr << "ERROR: Job list file has invalid name format." << std::endl;
      std::cerr << "Expected format: job_v{version}_{type}_{luminosity}_inv_pb_{number}.txt" << std::endl;
      std::cerr << "Valid types: data, all_phys, all_phys2, all_phys3" << std::endl;
      return 1;
    }

    useJobListFile = true;
  }

  // Set KLOE class instance
  KLOE::pm00 eventAnalysis;
  // Set logger for error logging
  TString logDirectory = Paths::base_path + Paths::logs_dir;
  ErrorHandling::ErrorLogs logger((std::string)logDirectory);
  ErrorHandling::InfoCodes infoCode;
  // -------------------------------------------------------------------
  // Initialize utility variables
  Utils::InitializeVariables(logger);
  // -------------------------------------------------------------------
  // Analysis flags and settings
  ConfigManager &config = ConfigManager::getInstance();
  config.setupLogger(&logger);

  KLOE::AnalysisConfig &analysisConfig = KLOE::AnalysisConfig::getInstance();
  analysisConfig.SetupLogger(&logger);
  analysisConfig.LoadFromFile(Paths::analysisConfigPath);
  // -------------------------------------------------------------------
  // Set tree name
  const std::string generalTreeName = Utils::properties["variables"]["tree"]["treename"]["general"];
  // -------------------------------------------------------------------
  // Set Menu instance
  Controls::Menu mainMenu(10); // For analysis options
  Controls::MainMenu mainMenuOpt;

  Controls::Menu mainMenuControlSample(12); // For control sample options
  Controls::MainMenuControlSample mainMenuControlSampleOpt;

  Controls::Menu *dataType = new Controls::Menu(2); // For data type
  // -------------------------------------------------------------------
  // Set global style for histograms
  KLOE::setGlobalStyle();
  // -------------------------------------------------------------------
  // Set config file watcher
  // ConfigWatcher cfgWatcher(Paths::propName);
  // cfgWatcher.start();
  // -------------------------------------------------------------------
  Controls::DataType dataTypeOpt;
  Controls::FileType fileTypeOpt;

  // Set flag for initial analysis
  Bool_t initialAnalysisExecution = Utils::properties["flags"]["initialAnalysisExec"]["flag"];

  // Initialize and fill the TChain object
  TChain chain(generalTreeName.c_str());

  KLOE::FileManager initObj(logger);
  MainMenuHandler mainMenuHandler;

  // -------------------------------------------------------------------
  // PROCESS JOB LIST FILE IF PROVIDED
  // -------------------------------------------------------------------
  if (useJobListFile)
  {
    try
    {
      // Wczytaj listę plików z pliku
      std::vector<std::string> fileList = KLOE::FileManager::LoadFileListFromFile(jobListFile);

      std::cout << "Loaded " << fileList.size() << " files from job list: " << jobListFile << std::endl;

      // Inicjalizuj TChain z listy plików
      initObj.chainInit(chain, fileList);

      infoCode = ErrorHandling::InfoCodes::FILE_ADDED;
      std::string infoMsg = "Initialized TChain with " + std::to_string(chain.GetEntries()) + " entries from job list.";
      logger.getLog(infoCode, infoMsg);

      // Pobierz informacje o typie analizy z nazwy pliku
      boost::filesystem::path jobListPath(jobListFile);
      std::string filename = jobListPath.filename().string();

      // Wyznacz fileTypeOpt z nazwy pliku
      // Format: job_v{version}_{type}_{luminosity}_inv_pb_{number}.txt

      const std::regex re(R"(^job_v[^_]+_[^_]+_[^_]+_inv_pb_(\d+)\.txt$)");
      std::smatch match;
      int jobNumber = -1;

      if (std::regex_match(filename, match, re) && match.size() > 1)
      {
        jobNumber = std::stoi(match[1].str());
        std::cout << "Job number: " << jobNumber << std::endl;
      }
      else
      {
        std::cerr << "ERROR: Cannot extract job number from filename: " << filename << std::endl;
        return 1;
      }

      if (filename.find("_data_") != std::string::npos)
      {
        fileTypeOpt = Controls::FileType::DATA;
        std::cout << "Analysis type: DATA" << std::endl;
      }
      else if (filename.find("_all_phys3_") != std::string::npos)
      {
        fileTypeOpt = Controls::FileType::ALL_PHYS3;
        std::cout << "Analysis type: ALL_PHYS3" << std::endl;
      }
      else if (filename.find("_all_phys2_") != std::string::npos)
      {
        fileTypeOpt = Controls::FileType::ALL_PHYS2;
        std::cout << "Analysis type: ALL_PHYS2" << std::endl;
      }
      else if (filename.find("_all_phys_") != std::string::npos)
      {
        fileTypeOpt = Controls::FileType::ALL_PHYS;
        std::cout << "Analysis type: ALL_PHYS" << std::endl;
      }
      else
      {
        std::cerr << "ERROR: Cannot determine analysis type from filename." << std::endl;
        return 1;
      }

      // Wykonaj początkową analizę
      infoCode = ErrorHandling::InfoCodes::FUNC_EXECUTED;
      logger.getLog(infoCode, "Initial analysis execution from job list file.");
      InitAnalysis_main(chain, fileTypeOpt, eventAnalysis, true, logger, jobNumber);

      logger.printErrStats();
      return 0;
    }
    catch (const std::exception &e)
    {
      std::cerr << "ERROR: " << e.what() << std::endl;
      ErrorHandling::ErrorCodes errCode = ErrorHandling::ErrorCodes::TREE_NOT_EXIST;
      logger.getErrLog(errCode, std::string("Failed to load job list file: ") + e.what());
      logger.printErrStats();
      return 1;
    }
  }

  if (!initialAnalysisExecution)
  {
    // Set the data type options
    auto callParams = InputParamsHandler::getParams(
        Utils::properties, KLOE::lastFileMax, Paths::propName, logger, dataType);

    initObj.chainInit(chain, callParams.dataTypeOpt, callParams.firstFile, callParams.lastFile, callParams.firstFile, callParams.lastFile, callParams.csFlag, 1);

    mainMenuHandler.runMenuLoop(
        callParams.csFlag,
        mainMenu,
        mainMenuOpt,
        mainMenuControlSample,
        mainMenuControlSampleOpt,
        chain,
        eventAnalysis,
        callParams.dataTypeOpt,
        logger,
        infoCode);
    // cfgWatcher);
  }
  else
  {
    std::ifstream rootFiles(Paths::rootfilesName);
    json filePaths = json::parse(rootFiles);

    std::cout << "Choose the file type to analyze: " << std::endl;
    std::cin >> fileTypeOpt;

    std::string
        DataPath = "",
        runRegexPattern = R"(.*_(\d{5})_v\d\.root$)";

    std::vector<std::string> DataPathList(std::begin(filePaths["MC"]["path"][2]), std::end(filePaths["MC"]["path"][2]));

    KLOE::RunStats runs;

    switch (fileTypeOpt)
    {
    case Controls::FileType::DATA:
    {
      DataPath = filePaths["Data"]["path"];
      runs = initObj.getRunStats(DataPath, runRegexPattern);
      initObj.chainInit(chain, DataPath, runRegexPattern,
                        runs.minRun, runs.maxRun);

      break;
    }
    case Controls::FileType::ALL_PHYS:
    {
      DataPath = filePaths["MC"]["path"][0];
      runs = initObj.getRunStats(DataPath, runRegexPattern);
      initObj.chainInit(chain, DataPath, runRegexPattern,
                        runs.minRun, runs.maxRun);
      break;
    }
    case Controls::FileType::ALL_PHYS2:
    {
      DataPath = filePaths["MC"]["path"][1];
      runs = initObj.getRunStats(DataPath, runRegexPattern);
      initObj.chainInit(chain, DataPath, runRegexPattern,
                        runs.minRun, runs.maxRun);
      break;
    }
    case Controls::FileType::ALL_PHYS3:
    {
      for (const auto &path : DataPathList)
      {
        runs = initObj.getRunStats(path, runRegexPattern);
        initObj.chainInit(chain, path, runRegexPattern,
                          runs.minRun, runs.minRun);
      }

      break;
    }
    default:
    {
      std::cerr << "Invalid file type selected." << std::endl;
      return 1;
    }
    }

    // -------------------------------------------------------
    infoCode = ErrorHandling::InfoCodes::FILE_ADDED;

    std::string infoMsg = "Initialized TChain with " + std::to_string(chain.GetEntries()) + " entries.";
    logger.getLog(infoCode, infoMsg);

    // -------------------------------------------------------
    infoCode = ErrorHandling::InfoCodes::FUNC_EXECUTED;

    logger.getLog(infoCode, "Initial analysis execution.");
    InitAnalysis_main(chain, fileTypeOpt, eventAnalysis, false, logger);
  }
  // -------------------------------------------------------------------
  // cfgWatcher.stop();
  delete dataType;
  return 0;
}
