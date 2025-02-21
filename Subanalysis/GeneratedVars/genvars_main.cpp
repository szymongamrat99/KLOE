#include <chrono>

#include <TMath.h>

#include "inc/genvars.hpp"

using namespace std;
using std::chrono::duration;
using std::chrono::duration_cast;
using std::chrono::high_resolution_clock;
using std::chrono::minutes;

int GenVars_main(TChain &chain, KLOE::pm00 &Obj, Controls::DataType &dataTypeOpt)
{
  // Set logger for error logging
  std::string logFilename = (std::string)gen_vars_dir + (std::string)logs_dir + "genVars_" + Obj.getCurrentDate() + ".log";
  ErrorHandling::ErrorLogs logger(logFilename);
  ErrorHandling::InfoCodes infoCode;
  // -------------------------------------------------------------------
  // Set Menu instance
  Controls::Menu *menu = new Controls::Menu(5);
  Controls::GenVars menuOpt;
  // -------------------------------------------------------------------

  bool
      dataTypeErr,
      menuRangeErr;

  do
  {
    menu->InitMenu();
    menu->ShowOpt();
    menu->EndMenu();

    try
    {
      cin >> menuOpt;

      dataTypeErr = !cin;
      menuRangeErr = menuOpt < Controls::GenVars::GEN_VARS || menuOpt > Controls::GenVars::EXIT;

      if (dataTypeErr)
      {
        throw ErrorHandling::ErrorCodes::DATA_TYPE;
      }
      else if (menuRangeErr)
      {
        throw ErrorHandling::ErrorCodes::MENU_RANGE;
      }
    }
    catch (ErrorHandling::ErrorCodes err)
    {
      cin.clear();
      cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

      logger.getErrLog(err);
    }

    switch (menuOpt)
    {
    case Controls::GenVars::GEN_VARS:
    {
      // genvars(firstFile, lastFile, 4);
      break;
    }
    case Controls::GenVars::SPLIT_CHANN:
    {
      infoCode = ErrorHandling::InfoCodes::FUNC_EXECUTED;
      logger.getLog(infoCode, "Split of files into MCTruth variable.");

      Obj.startTimer();
      split_channels(chain, dataTypeOpt, logger, Obj);

      infoCode = ErrorHandling::InfoCodes::FUNC_EXEC_TIME;
      logger.getLog(infoCode, Obj.endTimer());

      break;
    }
    case Controls::GenVars::EXIT:
    {
      break;
    }
    default:
      break;
    }

  } while (menuOpt != Controls::GenVars::EXIT);

  return 0;
}