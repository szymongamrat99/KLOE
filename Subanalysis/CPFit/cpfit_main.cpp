#include <chrono>

#include <TMath.h>
#include <TStyle.h>

#include "inc/cpfit.hpp"

#define __FILENAME__ (__builtin_strrchr(__FILE__, '/') ? __builtin_strrchr(__FILE__, '/') + 1 : __FILE__)

using namespace std;
using std::chrono::duration;
using std::chrono::duration_cast;
using std::chrono::high_resolution_clock;
using std::chrono::minutes;

int CPFit_main()
{
  Short_t good_clus = properties["variables"]["KinFit"]["Trilateration"]["goodClus"],
          loopcount = properties["variables"]["KinFit"]["Trilateration"]["loopCount"], 
          jmin = properties["variables"]["KinFit"]["Trilateration"]["bunchMin"], 
          jmax = properties["variables"]["KinFit"]["Trilateration"]["bunchMax"];

  ErrorHandling::ErrorLogs errLogger;
  LogsHandling::Logs logger;
  std::ofstream LogFileMain, LogFileTriangle, LogFileTri, LogFileTriKinFit;
  std::ofstream ErrFileMain, ErrFileTriangle, ErrFileTri, ErrFileTriKinFit;

  Controls::Menu *menu = new Controls::Menu(4);
  Controls::Menu *dataType = new Controls::Menu(2);

  logger.getErrLog(LogsHandling::GeneralLogs::TIME_STAMP);

  ErrFileMain.open(neutrec_dir + logs_dir + __FILENAME__ + ".err");

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
    errLogger.getErrLog(err);
    errLogger.getErrLog(err, ErrFileMain);
    errLogger.setErrCount(err);
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
    errLogger.getErrLog(err);
    errLogger.getErrLog(err, ErrFileMain);
    errLogger.setErrCount(err);
    return int(err);
  }

  Controls::DataType dataTypeOpt;

  try
  {
    dataType->InitMenu();
    dataType->ShowOpt();
    dataType->EndMenu();

    cin >> dataTypeOpt;

    if (!cin)
    {
      throw ErrorHandling::ErrorCodes::DATA_TYPE;
    }
    else if (dataTypeOpt < Controls::DataType::SIGNAL_TOT || dataTypeOpt > Controls::DataType::MC_DATA)
    {
      throw ErrorHandling::ErrorCodes::MENU_RANGE;
    }
  }
  catch (ErrorHandling::ErrorCodes err)
  {
    errLogger.getErrLog(err);
    errLogger.getErrLog(err, ErrFileMain);
    errLogger.setErrCount(err);
    return int(err);
  }

  Controls::CPFitMenu menuOpt;
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

      errLogger.getErrLog(err);
      errLogger.getErrLog(err, ErrFileMain);
      errLogger.setErrCount(err);
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
      auto start = std::chrono::system_clock::now();
      cp_fit_mc_data(firstFile, lastFile, "split", false, 10, 5, 1);
      auto end = std::chrono::system_clock::now();
      std::chrono::duration<double> elapsed_seconds = end - start;

      std::cout << elapsedTimeHMS(elapsed_seconds.count()) << std::endl;

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

  errLogger.errLogStats();
  errLogger.errLogStats(ErrFileMain);

  ErrFileMain.close();
  ErrFileMain.close();

  return 0;
}