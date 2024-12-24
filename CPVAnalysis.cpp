#include <chrono>
#include <iostream>
#include <string>

#include <const.h>

#include <MainMenu.h>
#include <ErrorLogs.h>

#include "inc/klspm00.hpp"

bool 
    firstFileRangeErr, 
    lastFileRangeErr, 
    dataTypeErr,
    menuRangeErr;

ErrorHandling::ErrorLogs logger;
Controls::Menu mainMenu(0);

using namespace std;
using namespace std::chrono;

int main(int argc, char *argv[])
{
  int first, last;
  Controls::MainMenu mainMenuOpt;

  setGlobalStyle();


  try
  {
    cout << "Choose the first file: ";
    cin >> first;
    cout << endl;

    dataTypeErr = !cin;
    firstFileRangeErr = first < 1 || first > lastFileMax;

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
    cin >> last;
    cout << endl;

    dataTypeErr = !cin;
    firstFileRangeErr = last < first || last > lastFileMax;

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

  do
  {
    mainMenu.InitMenu();
    cout << "1. Generated variables." << endl;
    cout << "2. K->pi+pi- reconstruction - no boost." << endl;
    cout << "3. K->pi+pi- reconstruction - boost." << endl;
    cout << "4. IP event by event - using K->pi+pi-." << endl;
    cout << "5. Omega-pi0 reconstruction." << endl;
    cout << "6. Regeneration rejection." << endl;
    cout << "7. K->pi0pi0 - using trilateration." << endl;
    cout << "8. K->pi0pi0 - using triangle, photons with trilateration." << endl;
    cout << "9. Efficiency and signal-to-background plots." << endl;
    cout << "10. Kinematic fits." << endl;
    cout << "11. Recalculation of variables to CM frame." << endl;
    cout << "12. CPV normalization." << endl;
    cout << "13. Exit." << endl;
    mainMenu.EndMenu();

    try
    {
      cin >> mainMenuOpt;

      dataTypeErr = !cin;
      menuRangeErr = mainMenuOpt < Controls::MainMenu::GEN_VARS || mainMenuOpt > Controls::MainMenu::EXIT;

      if (dataTypeErr)
      {
        throw ErrorHandling::ErrorCodes::DATA_TYPE;
      }
      else if (menuRangeErr)
      {
        throw ErrorHandling::ErrorCodes::RANGE;
      }
    }
    catch (ErrorHandling::ErrorCodes err)
    {
      cin.clear();
      cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

      logger.getErrLog(err);
    }

    switch (mainMenuOpt)
    {
    case Controls::MainMenu::GEN_VARS:
    {
      GenVars_main();
    }
    case Controls::MainMenu::KCHREC_NO_BOOST:
      break;
    case Controls::MainMenu::KCHREC_BOOST:
      break;
    case Controls::MainMenu::IP_EV_BY_EV:
      break;
    case Controls::MainMenu::OMEGA_REC:
    {
      OmegaRec_main();
    }
    case Controls::MainMenu::REGEN_REJ:
    {
      Regen_main();
    }
    case Controls::MainMenu::KNEREC_TRILAT:
      break;
    case Controls::MainMenu::KNEREC_TRIANGLE:
      break;
    case Controls::MainMenu::EFF_SIG_TO_BCG:
      break;
    case Controls::MainMenu::KIN_FITS:
      break;
    case Controls::MainMenu::TRANSF_TO_CM:
      break;
    case Controls::MainMenu::CPV_NORM:
      CPFit_main();
    }

  } while (mainMenuOpt != Controls::MainMenu::EXIT);

  return 0;
}
