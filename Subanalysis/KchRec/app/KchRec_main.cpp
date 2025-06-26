#include "../inc/kchrec.hpp"

using namespace std;

int KchRec_main(TChain &chain, KLOE::pm00 &Obj, Controls::DataType &dataTypeOpt)
{
  // Set logger for error logging
  std::string logFilename = (std::string)charged_dir + (std::string)logs_dir + "KchRec_" + Obj.getCurrentDate() + ".log";
  ErrorHandling::ErrorLogs logger(logFilename);
  // -----------------------------------------------------------------------------------
  // Set Menu instance
  Controls::Menu *menu = new Controls::Menu(0); // For analysis options
  Controls::KchRecMenu menuOpt;
  // -----------------------------------------------------------------------------------

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
      menuRangeErr = menuOpt < Controls::KchRecMenu::K_MASS || menuOpt > Controls::KchRecMenu::EXIT;

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

      logger.getErrLog(err, "KchRec analysis option choice");
    }

    switch (menuOpt)
    {
    case Controls::KchRecMenu::K_MASS:
    {
      Obj.startTimer();
      kchrec_Kmass(chain, dataTypeOpt, logger, Obj);
      Obj.endTimer();

      break;
    }
    case Controls::KchRecMenu::KSKL:
    {
      Obj.startTimer();
      // kchrec_KSKL(chain, dataTypeOpt, logger, Obj);
      Obj.endTimer();

      break;
    }
    case Controls::KchRecMenu::CLOSEST:
    {
      Obj.startTimer();
      // kchrec_Closest(chain, dataTypeOpt, logger, Obj);
      Obj.endTimer();

      break;
    }
    case Controls::KchRecMenu::BOOST:
    {

      break;
    }
    case Controls::KchRecMenu::EXIT:
    {

      break;
    }
    default:
      break;
    }

  } while (menuOpt != Controls::KchRecMenu::EXIT);

  logger.printErrStats();

  return 0;
}