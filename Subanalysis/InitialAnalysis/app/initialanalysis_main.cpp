#include <chrono>
#include <TMath.h>

#include "../inc/initialanalysis.hpp"

using namespace std;

int InitAnalysis_main(TChain &chain, Controls::FileType &fileTypeOpt, KLOE::pm00 &Obj, bool singleFile, ErrorHandling::ErrorLogs &logger, Int_t jobNumber = -1)
{
  // Set logger for error logging
  ErrorHandling::InfoCodes infoCode;
  // -------------------------------------------------------------------
  // Set Menu instance
  std::unique_ptr<Controls::Menu> menu = std::make_unique<Controls::Menu>(13);
  Controls::InitialAnalysisMenu menuOpt;
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
      menuRangeErr = menuOpt < Controls::InitialAnalysisMenu::INIT_ANALYSIS_FULL || menuOpt > Controls::InitialAnalysisMenu::EXIT;

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
    case Controls::InitialAnalysisMenu::INIT_ANALYSIS_FULL:
    {
      infoCode = ErrorHandling::InfoCodes::FUNC_EXECUTED;
      logger.getLog(infoCode, "Full initial analysis");

      Obj.startTimer();
      InitialAnalysis_full(chain, fileTypeOpt, logger, Obj, singleFile, jobNumber);
      
      infoCode = ErrorHandling::InfoCodes::FUNC_EXEC_TIME;
      logger.getLog(infoCode, Obj.endTimer());
      break;
    }
    case Controls::InitialAnalysisMenu::EXIT:
    {
      break;
    }
    default:
      break;
    }

  } while (menuOpt != Controls::InitialAnalysisMenu::EXIT);
  
  return 0;
}