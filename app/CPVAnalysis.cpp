#include "../Include/klspm00.hpp"

#include "../Include/MainMenuHandler/MainMenuHandler.h"
#include "../Include/MainMenuHandler/InputParamsHandler.h"
#include <FileManager.h>
#include <event_data.h>

using namespace std;
using namespace std::chrono;

int main(int argc, char *argv[])
{
  bool
      firstFileRangeErr,
      lastFileRangeErr,
      dataTypeErr,
      menuRangeErr;

  UInt_t firstFile, lastFile, csFlag;

  // Set KLOE class instance
  KLOE::pm00 eventAnalysis;
  // -------------------------------------------------------------------
  // Set logger for error logging
  std::string logFilename = (std::string)base_path + (std::string)logs_dir + "general.prog_" + eventAnalysis.getCurrentDate() + ".log";
  ErrorHandling::ErrorLogs logger(logFilename);
  ErrorHandling::InfoCodes infoCode;
  // -------------------------------------------------------------------
  // Set tree name
  const std::string generalTreeName = properties["variables"]["tree"]["treename"]["general"];
  // -------------------------------------------------------------------
  // Set Menu instance
  Controls::Menu mainMenu(10); // For analysis options
  Controls::MainMenu mainMenuOpt;

  Controls::Menu mainMenuControlSample(12); // For control sample options
  Controls::MainMenuControlSample mainMenuControlSampleOpt;

  Controls::Menu *dataType = new Controls::Menu(2); // For data type
  // -------------------------------------------------------------------
  // Set global style for histograms
  setGlobalStyle();
  // -------------------------------------------------------------------
  // Set config file watcher
  // ConfigWatcher cfgWatcher(propName);
  // cfgWatcher.start();
  // -------------------------------------------------------------------
  Controls::DataType dataTypeOpt;
  Controls::FileType fileTypeOpt;

  // Set flag for initial analysis
  Bool_t initialAnalysisExecution = properties["flags"]["initialAnalysisExec"]["flag"];

  // Initialize and fill the TChain object
  TChain chain(generalTreeName.c_str());

  KLOE::FileManager initObj;
  MainMenuHandler mainMenuHandler;

  if (!initialAnalysisExecution)
  {
    // Set the data type options
    auto callParams = InputParamsHandler::getParams(
        properties, lastFileMax, propName, logger, dataType);

    initObj.chainInit(chain, callParams.dataTypeOpt, callParams.firstFile, callParams.lastFile, callParams.firstFile, callParams.lastFile, logger, callParams.csFlag);

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
    std::ifstream rootFiles(rootfilesName);
    json filePaths = json::parse(rootFiles);

    std::cout << "Choose the file type to analyze: " << std::endl;
    std::cin >> fileTypeOpt;

    std::string
        DataPath = "",
        runRegexPattern = R"(.*_(\d{5})_v2\.root$)";

    std::vector<std::string> DataPathList(std::begin(filePaths["MC"]["path"][2]), std::end(filePaths["MC"]["path"][2]));

    KLOE::RunStats runs;

    switch (fileTypeOpt)
    {
    case Controls::FileType::DATA:
    {
      DataPath = filePaths["Data"]["path"];
      runs = initObj.getRunStats(DataPath, runRegexPattern);
      initObj.chainInit(chain, logger, DataPath, runRegexPattern,
                        runs.minRun, runs.maxRun);

      break;
    }
    case Controls::FileType::ALL_PHYS:
    {
      DataPath = filePaths["MC"]["path"][0];
      runs = initObj.getRunStats(DataPath, runRegexPattern);
      initObj.chainInit(chain, logger, DataPath, runRegexPattern,
                        runs.minRun, runs.minRun);
      break;
    }
    case Controls::FileType::ALL_PHYS2:
    {
      DataPath = filePaths["MC"]["path"][1];
      runs = initObj.getRunStats(DataPath, runRegexPattern);
      initObj.chainInit(chain, logger, DataPath, runRegexPattern,
                        runs.minRun, runs.minRun);
      break;
    }
    case Controls::FileType::ALL_PHYS3:
    {
      for (const auto &path : DataPathList)
      {
        runs = initObj.getRunStats(path, runRegexPattern);
        initObj.chainInit(chain, logger, path, runRegexPattern,
                          runs.minRun, runs.minRun + 1);
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
    InitAnalysis_main(chain, fileTypeOpt, eventAnalysis);
  }
  // -------------------------------------------------------------------
  logger.printErrStats();

  // cfgWatcher.stop();
  return 0;
}
