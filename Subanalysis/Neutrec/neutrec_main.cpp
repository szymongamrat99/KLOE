#include <chrono>

#include <TMath.h>
#include <TStyle.h>

#include "inc/trilateration.hpp"


int Neutrec_main(TChain &chain, KLOE::pm00 &Obj, Controls::DataType &dataTypeOpt)
{

  Short_t 
          loopcount = properties["variables"]["KinFit"]["Trilateration"]["loopCount"],
          numOfConstraints = properties["variables"]["KinFit"]["Trilateration"]["numOfConstraints"],
          jmin = properties["variables"]["KinFit"]["Trilateration"]["bunchMin"],
          jmax = properties["variables"]["KinFit"]["Trilateration"]["bunchMax"];

  // Set logger for error logging
  std::string logFilename = (std::string)neutrec_dir + (std::string)logs_dir + "neutRec.log";
  ErrorHandling::ErrorLogs logger(logFilename);
  ErrorHandling::InfoCodes infoCode;
  // -------------------------------------------------------------------
  // Set Menu instance
  Controls::Menu *menu = new Controls::Menu(1);
  Controls::NeutRecMenu menuOpt;
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
      std::cin >> menuOpt;

      dataTypeErr = !std::cin;
      menuRangeErr = menuOpt < Controls::NeutRecMenu::BARE_TRILATERATION || menuOpt > Controls::NeutRecMenu::EXIT;

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
      std::cin.clear();
      std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

      logger.getErrLog(err);
    }

    switch (menuOpt)
    {
    case Controls::NeutRecMenu::BARE_TRILATERATION:
    {
      // auto start = std::chrono::system_clock::now();
      // tri_neurec(chain, dataTypeOpt, logger, Obj);
      // auto end = std::chrono::system_clock::now();
      // std::chrono::duration<double> elapsed_seconds = end - start;

      // std::cout << elapsedTimeHMS(elapsed_seconds.count()) << std::endl;

      break;
    }
    case Controls::NeutRecMenu::TRILATERATION_KIN_FIT:
    {
      auto start = std::chrono::system_clock::now();
      TrilaterationNeurecKinfit(chain, dataTypeOpt, logger, Obj);
      auto end = std::chrono::system_clock::now();
      std::chrono::duration<double> elapsed_seconds = end - start;

      std::cout << elapsedTimeHMS(elapsed_seconds.count()) << std::endl;

      break;
    }
    case Controls::NeutRecMenu::TRIANGLE:
    {
      infoCode = ErrorHandling::InfoCodes::FUNC_EXECUTED;
      logger.getLog(infoCode, "Neutral particles reconstruction: Triangle");

      Obj.startTimer();
      TriangleNeurec(chain, dataTypeOpt, logger, Obj);
      
      infoCode = ErrorHandling::InfoCodes::FUNC_EXEC_TIME;
      logger.getLog(infoCode, Obj.endTimer());

      break;
    }
    case Controls::NeutRecMenu::COMP_OF_MET:
    {
      infoCode = ErrorHandling::InfoCodes::FUNC_EXECUTED;
      logger.getLog(infoCode, "Comparison of Neutral Reconstruction Methods");

      Obj.startTimer();
      CompOfMethods(chain, dataTypeOpt, logger, Obj);
      
      infoCode = ErrorHandling::InfoCodes::FUNC_EXEC_TIME;
      logger.getLog(infoCode, Obj.endTimer());

      break;
    }
    case Controls::NeutRecMenu::EXIT:
    {
      break;
    }
    default:
      break;
    }

  } while (menuOpt != Controls::NeutRecMenu::EXIT);

  logger.printErrStats();

  return 0;
}