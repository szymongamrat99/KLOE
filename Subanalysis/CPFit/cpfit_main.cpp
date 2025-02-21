#include <chrono>

#include <TMath.h>
#include <TStyle.h>

#include "inc/cpfit.hpp"

using namespace std;

int CPFit_main(TChain &chain, KLOE::pm00 &Obj, Controls::DataType &dataTypeOpt)
{
  Short_t 
          loopcount = properties["variables"]["KinFit"]["Trilateration"]["loopCount"],
          numOfConstraints = properties["variables"]["KinFit"]["Trilateration"]["numOfConstraints"],
          jmin = properties["variables"]["KinFit"]["Trilateration"]["bunchMin"],
          jmax = properties["variables"]["KinFit"]["Trilateration"]["bunchMax"];

  // Set logger for error logging
  std::string logFilename = (std::string)cpfit_dir + (std::string)logs_dir + "cpfit.log";
  ErrorHandling::ErrorLogs logger(logFilename);
  ErrorHandling::InfoCodes infoCode;
  // -------------------------------------------------------------------
  // Set Menu instance
  Controls::Menu *menu = new Controls::Menu(4);
  Controls::CPFitMenu menuOpt;
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
      menuRangeErr = menuOpt < Controls::CPFitMenu::HALF_SIGNAL_MC || menuOpt > Controls::CPFitMenu::EXIT;

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
    case Controls::CPFitMenu::HALF_SIGNAL_MC:
    {
      break;
    }
    case Controls::CPFitMenu::HALF_SIG_BCG_MC:
    {

      break;
    }
    case Controls::CPFitMenu::MC_DATA:
    {
      infoCode = ErrorHandling::InfoCodes::FUNC_EXECUTED;
      logger.getLog(infoCode, "CP Final Fit");

      Obj.startTimer();
      cp_fit_mc_data(chain, "split", false, dataTypeOpt, logger, Obj);
      
      infoCode = ErrorHandling::InfoCodes::FUNC_EXEC_TIME;
      logger.getLog(infoCode, Obj.endTimer());

      break;
    }
    case Controls::CPFitMenu::EXIT:
    {
      break;
    }
    default:
      break;
    }

  } while (menuOpt != Controls::CPFitMenu::EXIT);

  logger.printErrStats();

  return 0;
}