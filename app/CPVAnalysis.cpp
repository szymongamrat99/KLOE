#include "../Include/klspm00.hpp"

#include "../Include/MainMenuHandler/MainMenuHandler.h"

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
  ConfigWatcher cfgWatcher(propName);
  cfgWatcher.start();
  // -------------------------------------------------------------------

  try
  {
    cout << "Measurement (1) / Control Sample (2)?";
    cin >> csFlag;
    cout << endl;

    dataTypeErr = !cin;

    if (dataTypeErr)
    {
      throw ErrorHandling::ErrorCodes::DATA_TYPE;
    }
    else if (csFlag != 1 && csFlag != 2)
    {
      throw ErrorHandling::ErrorCodes::RANGE;
    }
  }
  catch (ErrorHandling::ErrorCodes err)
  {
    cin.clear();
    cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    logger.getErrLog(err);
    return int(err);
  }

  try
  {
    cout << "Choose the first file: ";
    cin >> firstFile;
    cout << endl;

    dataTypeErr = !cin;
    firstFileRangeErr = firstFile < 1 || firstFile > lastFileMax;

    if (dataTypeErr)
    {
      throw ErrorHandling::ErrorCodes::DATA_TYPE;
    }
    else if (firstFileRangeErr)
    {
      throw ErrorHandling::ErrorCodes::RANGE;
    }
  }
  catch (ErrorHandling::ErrorCodes err)
  {
    cin.clear();
    cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    logger.getErrLog(err);
    return int(err);
  }

  try
  {
    cout << "Choose the last file: ";
    cin >> lastFile;
    cout << endl;

    dataTypeErr = !cin;
    lastFileRangeErr = lastFile < firstFile || lastFile > lastFileMax;

    if (dataTypeErr)
    {
      throw ErrorHandling::ErrorCodes::DATA_TYPE;
    }
    else if (lastFileRangeErr)
    {
      throw ErrorHandling::ErrorCodes::RANGE;
    }
  }
  catch (ErrorHandling::ErrorCodes err)
  {
    cin.clear();
    cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    logger.getErrLog(err);
    return int(err);
  }

  properties["variables"]["rootFiles"]["Data"]["firstFile"] = firstFile;
  properties["variables"]["rootFiles"]["Data"]["lastFile"] = lastFile;
  properties["variables"]["rootFiles"]["MC"]["firstFile"] = firstFile;
  properties["variables"]["rootFiles"]["MC"]["lastFile"] = lastFile;

  std::ofstream outfile(propName);
  outfile << properties.dump(4);
  outfile.close();

  Controls::DataType dataTypeOpt;

  try
  {
    dataType->InitMenu();
    dataType->ShowOpt();
    dataType->EndMenu();

    cin >> dataTypeOpt;

    if (!cin)
    {
      throw ErrorHandling::ErrorCodes::DATA_TYPE;
    }
    else if (dataTypeOpt < Controls::DataType::SIGNAL_TOT || dataTypeOpt > Controls::DataType::SIGNAL_MAX)
    {
      throw ErrorHandling::ErrorCodes::MENU_RANGE;
    }
  }
  catch (ErrorHandling::ErrorCodes err)
  {
    logger.getErrLog(err);
    return int(err);
  }

  // Initialize and fill the TChain object
  TChain chain(generalTreeName.c_str());
  eventAnalysis.chainInit(chain, dataTypeOpt, firstFile, lastFile, firstFile, lastFile, logger, csFlag);
  // -------------------------------------------------------------------

  MainMenuHandler mainMenuHandler;

  mainMenuHandler.runMenuLoop(
      csFlag,
      mainMenu,
      mainMenuOpt,
      mainMenuControlSample,
      mainMenuControlSampleOpt,
      chain,
      eventAnalysis,
      dataTypeOpt,
      logger,
      infoCode,
      cfgWatcher
  );
  // -------------------------------------------------------------------
  logger.printErrStats();

  cfgWatcher.stop();
  return 0;
}
