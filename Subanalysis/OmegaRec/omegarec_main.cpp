#include <chrono>

#include <TMath.h>

#include "inc/omegarec.hpp"

using namespace std;

int OmegaRec_main(TChain &chain, KLOE::pm00 &Obj, Controls::DataType &dataTypeOpt)
{
  // Set logger for error logging
  std::string logFilename = (std::string)omegarec_dir + (std::string)logs_dir + "omegaRec_" + Obj.getCurrentDate() + ".log";
  ErrorHandling::ErrorLogs logger(logFilename);
  // -------------------------------------------------------------------
  // Set Menu instance
  Controls::Menu *menu = new Controls::Menu(6);
  Controls::OmegaRec menuOpt;
  // -------------------------------------------------------------------

  Short_t good_clus = properties["variables"]["KinFit"]["Trilateration"]["goodClus"],
          loopcount = properties["variables"]["KinFit"]["Trilateration"]["loopCount"],
          numOfConstraints = properties["variables"]["KinFit"]["Trilateration"]["numOfConstraints"],
          jmin = properties["variables"]["KinFit"]["Trilateration"]["bunchMin"],
          jmax = properties["variables"]["KinFit"]["Trilateration"]["bunchMax"];

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
      menuRangeErr = menuOpt < Controls::OmegaRec::OMEGA_REC || menuOpt > Controls::OmegaRec::EXIT;

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
    case Controls::OmegaRec::OMEGA_REC:
    {
      Obj.startTimer();
      omegarec(chain, dataTypeOpt, logger, Obj);
      std::cout << Obj.endTimer() << std::endl;

      break;
    }
    case Controls::OmegaRec::OMEGA_CUTS:
    {
      // omegacuts(firstFile, lastFile);
      break;
    }
    case Controls::OmegaRec::PLOTS:
    {
      Obj.startTimer();
      plots(chain, loopcount, numOfConstraints, jmin, jmax, dataTypeOpt, logger);
      std::cout << Obj.endTimer() << std::endl;

      break;
    }
    case Controls::OmegaRec::EXIT:
    {
      break;
    }
    default:
      break;
    }

  } while (menuOpt != Controls::OmegaRec::EXIT);

  logger.printErrStats();

  return 0;
}