#include <chrono>

#include <TMath.h>
#include <TStyle.h>

#include "inc/trilateration.hpp"

using namespace std;
using std::chrono::duration;
using std::chrono::duration_cast;
using std::chrono::high_resolution_clock;
using std::chrono::minutes;

int main(int argc, char *argv[])
{
  setGlobalStyle();

  ErrorHandling::ErrorLogs errLogger;
  LogsHandling::Logs logger;
  std::ofstream LogFileMain, LogFileTriangle, LogFileTri, LogFileTriKinFit;
  std::ofstream ErrFileMain, ErrFileTriangle, ErrFileTri, ErrFileTriKinFit;

  Controls::Menu *menu = new Controls::Menu(1);
  Controls::Menu *dataType = new Controls::Menu(2);

  logger.getErrLog(LogsHandling::GeneralLogs::TIME_STAMP);

  ErrFileMain.open(neutrec_dir + logs_dir + "NeutRecMain.err");

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

  Short_t good_clus = atoi(argv[3]), loopcount = atoi(argv[4]), jmin = atoi(argv[5]), jmax = atoi(argv[6]);

  Controls::NeutRecMenu menuOpt;
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
      cin.clear();
      cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

      errLogger.getErrLog(err);
      errLogger.getErrLog(err, ErrFileMain);
      errLogger.setErrCount(err);
    }

    switch (menuOpt)
    {
    case Controls::NeutRecMenu::BARE_TRILATERATION:
    {
      tri_neurec(ind_data_mc, firstFile, lastFile, good_clus);
      break;
    }
    case Controls::NeutRecMenu::TRILATERATION_KIN_FIT:
    {
      auto start = std::chrono::system_clock::now();
      tri_neurec_kinfit_corr(firstFile, lastFile, 10, 0, 0, dataTypeOpt);
      auto end = std::chrono::system_clock::now();
      std::chrono::duration<double> elapsed_seconds = end - start;

      break;
    }
    case Controls::NeutRecMenu::TRIANGLE:
    {
      auto start = std::chrono::system_clock::now();
      triangle_neurec(firstFile, lastFile, 10, 5, 1, dataTypeOpt, 0);
      auto end = std::chrono::system_clock::now();
      std::chrono::duration<double> elapsed_seconds = end - start;

      std::cout << elapsedTimeHMS(elapsed_seconds.count()) << std::endl;

      break;
    }
    case Controls::NeutRecMenu::COMP_OF_MET:
    {
      auto start = std::chrono::system_clock::now();
      comp_of_methods(firstFile, lastFile, 10, 5, 1, dataTypeOpt);
      auto end = std::chrono::system_clock::now();
      std::chrono::duration<double> elapsed_seconds = end - start;

      std::cout << elapsedTimeHMS(elapsed_seconds.count()) << std::endl;

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

  errLogger.errLogStats();
  errLogger.errLogStats(ErrFileMain);

  ErrFileMain.close();

  return 0;
}