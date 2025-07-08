#include <chrono>
#include <TMath.h>

#include "../inc/initialanalysis.hpp"

using namespace std;

int InitAnalysis_main(TChain &chain, KLOE::pm00 &Obj, Controls::DataType &dataTypeOpt)
{
  // Set logger for error logging
  std::string logFilename = (std::string)initialanalysis_dir + (std::string)logs_dir + "InitialAnalysis.log";
  ErrorHandling::ErrorLogs logger(logFilename);
  ErrorHandling::InfoCodes infoCode;
  // -------------------------------------------------------------------
  // Set Menu instance
  Controls::Menu *menu = new Controls::Menu(11);
  Controls::CovMatrix menuOpt;
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
      menuRangeErr = menuOpt < Controls::CovMatrix::USING_MC_DATA || menuOpt > Controls::CovMatrix::EXIT;

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
    case Controls::CovMatrix::USING_CONTROL_SAMPLE:
    {
      infoCode = ErrorHandling::InfoCodes::FUNC_EXECUTED;
      logger.getLog(infoCode, "Full initial analysis");

      Obj.startTimer();
      InitialAnalysis_full(chain, dataTypeOpt, logger, Obj);
      
      infoCode = ErrorHandling::InfoCodes::FUNC_EXEC_TIME;
      logger.getLog(infoCode, Obj.endTimer());
      break;
    }
    case Controls::CovMatrix::EXIT:
    {
      break;
    }
    default:
      break;
    }

  } while (menuOpt != Controls::CovMatrix::EXIT);

  logger.printErrStats();

  return 0;
}