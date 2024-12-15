#include <chrono>

#include <TMath.h>

#include "inc/regenrejec.hpp"

using namespace std;
using std::chrono::duration;
using std::chrono::duration_cast;
using std::chrono::high_resolution_clock;
using std::chrono::minutes;

int Regen_main()
{
  Controls::Menu *menu = new Controls::Menu(8);
  ErrorHandling::ErrorLogs logger;
  ofstream LogFile;
  LogFile.open(gen_vars_dir + logs_dir + "Regen.log");
  setGlobalStyle();

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

  Controls::Regen menuOpt;
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
      menuRangeErr = menuOpt < Controls::Regen::REGEN_REJEC_TEST || menuOpt > Controls::Regen::EXIT;

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
    case Controls::Regen::REGEN_REJEC_TEST:
    {
      auto start = std::chrono::system_clock::now();
      regenrejec(firstFile, lastFile);
      auto end = std::chrono::system_clock::now();
      std::chrono::duration<double> elapsed_seconds = end - start;
      break;
    }
    case Controls::Regen::PLOTS:
    {
      auto start = std::chrono::system_clock::now();
      plots(firstFile, lastFile);
      auto end = std::chrono::system_clock::now();
      std::chrono::duration<double> elapsed_seconds = end - start;
      break;
    }
    case Controls::Regen::EXIT:
    {
      break;
    }
    default:
      break;
    }

  } while (menuOpt != Controls::Regen::EXIT);

  LogFile.close();

  return 0;
}