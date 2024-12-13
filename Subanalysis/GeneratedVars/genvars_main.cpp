#include <chrono>

#include <TMath.h>

#include "inc/genvars.hpp"

using namespace std;
using std::chrono::duration;
using std::chrono::duration_cast;
using std::chrono::high_resolution_clock;
using std::chrono::minutes;

int GenVars_main()
{
  Controls::Menu *menu = new Controls::Menu(5);
  ErrorHandling::ErrorLogs logger;
  ofstream LogFile;
  LogFile.open(gen_vars_dir + logs_dir + "GenVars.log");

  Int_t firstFile, lastFile, ind_data_mc;

  try
  {
    cout << "Input the first file number to be analyzed: ";
    cin >> firstFile;

    if (!cin)
    {
      throw ErrorHandling::ErrorCodes::DATA_TYPE;
    }
    else if (firstFile < firstFileMax || firstFile > lastFileMax)
    {
      throw ErrorHandling::ErrorCodes::RANGE;
    }
  }
  catch (ErrorHandling::ErrorCodes err)
  {
    logger.getErrLog(err);
    logger.getErrLog(err, LogFile);
    return int(err);
  }

  try
  {
    cout << "Input the last file number to be analyzed: ";
    cin >> lastFile;

    if (!cin)
    {
      throw ErrorHandling::ErrorCodes::DATA_TYPE;
    }
    else if (lastFile < firstFile || lastFile > lastFileMax)
    {
      throw ErrorHandling::ErrorCodes::RANGE;
    }
  }
  catch (ErrorHandling::ErrorCodes err)
  {
    logger.getErrLog(err);
    logger.getErrLog(err, LogFile);
    return int(err);
  }

  Controls::GenVars menuOpt;
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
      genvars(firstFile, lastFile, 4);
      break;
    }
    case Controls::GenVars::SPLIT_CHANN:
    {
      split_channels(firstFile, lastFile);
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

  LogFile.close();

  return 0;
}