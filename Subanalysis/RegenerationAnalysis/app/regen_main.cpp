#include <chrono>

#include <TMath.h>

#include "../inc/regenrejec.hpp"

using namespace std;
using std::chrono::duration;
using std::chrono::duration_cast;
using std::chrono::high_resolution_clock;
using std::chrono::minutes;

int Regen_main(TChain &chain, KLOE::pm00 &Obj, Controls::DataType &dataTypeOpt, ErrorHandling::ErrorLogs &logger)
{
  // Set logger for error logging
  ErrorHandling::InfoCodes infoCode;
  // -------------------------------------------------------------------
  // Set Menu instance
  std::unique_ptr<Controls::Menu> menu = std::make_unique<Controls::Menu>(8);
  Controls::Regen menuOpt;
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
      menuRangeErr = int(menuOpt) < int(Controls::Regen::REGEN_REJEC_TEST) || int(menuOpt) > int(Controls::Regen::EXIT);

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
      regenrejec(chain, Obj, dataTypeOpt, logger);
      auto end = std::chrono::system_clock::now();
      std::chrono::duration<double> elapsed_seconds = end - start;
      break;
    }
    case Controls::Regen::PLOTS:
    {
      auto start = std::chrono::system_clock::now();
      // plots(chain, Obj, dataTypeOpt, logger);
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

  return 0;
}