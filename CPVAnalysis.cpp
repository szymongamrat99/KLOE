#include "inc/klspm00.hpp"

using namespace std;
using namespace std::chrono;

int main(int argc, char *argv[])
{
  bool
      firstFileRangeErr,
      lastFileRangeErr,
      dataTypeErr,
      menuRangeErr;

  UInt_t firstFile, lastFile, csFlag;

  // Set KLOE class instance
  KLOE::pm00 eventAnalysis;
  // -------------------------------------------------------------------
  // Set logger for error logging
  std::string logFilename = (std::string)logs_dir + "general.prog_" + eventAnalysis.getCurrentDate() + ".log";
  ErrorHandling::ErrorLogs logger(logFilename);
  ErrorHandling::InfoCodes infoCode;
  // -------------------------------------------------------------------
  // Set Menu instance
  Controls::Menu mainMenu(10); // For analysis options
  Controls::MainMenu mainMenuOpt;

  Controls::Menu mainMenuControlSample(12); // For control sample options
  Controls::MainMenuControlSample mainMenuControlSampleOpt;

  Controls::Menu *dataType = new Controls::Menu(2); // For data type
  // -------------------------------------------------------------------
  // Set global style for histograms
  setGlobalStyle();
  // -------------------------------------------------------------------

  try
  {
    cout << "Measurement (1) / Control Sample (2)?";
    cin >> csFlag;
    cout << endl;

    dataTypeErr = !cin;

    if (dataTypeErr)
    {
      throw ErrorHandling::ErrorCodes::DATA_TYPE;
    }
    else if (csFlag != 1 && csFlag != 2)
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
    cout << "Choose the first file: ";
    cin >> firstFile;
    cout << endl;

    dataTypeErr = !cin;
    firstFileRangeErr = firstFile < 1 || firstFile > lastFileMax;

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
    cin >> lastFile;
    cout << endl;

    dataTypeErr = !cin;
    lastFileRangeErr = lastFile < firstFile || lastFile > lastFileMax;

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

  properties["variables"]["rootFiles"]["Data"]["firstFile"] = firstFile;
  properties["variables"]["rootFiles"]["Data"]["lastFile"] = lastFile;
  properties["variables"]["rootFiles"]["MC"]["firstFile"] = firstFile;
  properties["variables"]["rootFiles"]["MC"]["lastFile"] = lastFile;

  std::ofstream outfile(propName);
  outfile << properties.dump(4);
  outfile.close();

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
    else if (dataTypeOpt < Controls::DataType::SIGNAL_TOT || dataTypeOpt > Controls::DataType::SIGNAL_MAX)
    {
      throw ErrorHandling::ErrorCodes::MENU_RANGE;
    }
  }
  catch (ErrorHandling::ErrorCodes err)
  {
    logger.getErrLog(err);
    return int(err);
  }

  // Initialize and fill the TChain object
  TChain chain("INTERF/h1");
  eventAnalysis.chainInit(chain, dataTypeOpt, firstFile, lastFile, firstFile, lastFile, logger, csFlag);
  // -------------------------------------------------------------------

  std::cout << "Number of entries in the batch: " << chain.GetEntries() << std::endl;
  std::cout << "Number of entries - signal with improperly reconstructed Kch: " << chain.GetEntries("mctruth == 0") << std::endl;
  std::cout << "Number of entries - signal: " << chain.GetEntries("mctruth == 1") << std::endl;
  std::cout << "Number of entries - regeneration: " << chain.GetEntries("mctruth == 3") << std::endl;
  std::cout << "Number of entries - omega: " << chain.GetEntries("mctruth == 4") << std::endl;
  std::cout << "Number of entries - threepi0: " << chain.GetEntries("mctruth == 5") << std::endl;
  std::cout << "Number of entries - semileptonic: " << chain.GetEntries("mctruth == 6") << std::endl;
  std::cout << "Number of entries - pi+pi-pi+pi-: " << chain.GetEntries("mctruth == 8") << std::endl;
  std::cout << "Number of entries - other bcg: " << chain.GetEntries("mctruth == 7") << std::endl;

  do
  {
    mainMenu.InitMenu();
    switch (csFlag)
    {
    case 1:
    {
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
      cout << "13. Plots." << endl;
      cout << "14. Covariance Matrix Determination." << endl;
      cout << "15. Exit." << endl;

      break;
    }
    case 2:
    {
      mainMenuControlSample.ShowOpt();
      break;
    }
    }
    mainMenu.EndMenu();

    try
    {
      if (csFlag == 1)
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
      else if (csFlag == 2)
      {
        cin >> mainMenuControlSampleOpt;

        dataTypeErr = !cin;
        menuRangeErr = mainMenuControlSampleOpt < Controls::MainMenuControlSample::COV_MATRIX || mainMenuControlSampleOpt > Controls::MainMenuControlSample::EXIT;

        if (dataTypeErr)
        {
          throw ErrorHandling::ErrorCodes::DATA_TYPE;
        }
        else if (menuRangeErr)
        {
          throw ErrorHandling::ErrorCodes::RANGE;
        }
      }
    }
    catch (ErrorHandling::ErrorCodes err)
    {
      cin.clear();
      cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

      logger.getErrLog(err);
    }

    if (csFlag == 1)
    {
      switch (mainMenuOpt)
      {
      case Controls::MainMenu::GEN_VARS:
      {
        infoCode = ErrorHandling::InfoCodes::FUNC_EXECUTED;
        logger.getLog(infoCode, "MC generated variables for analysis");

        GenVars_main(chain, eventAnalysis, dataTypeOpt);
        break;
      }
      case Controls::MainMenu::KCHREC_NO_BOOST:
      {
        infoCode = ErrorHandling::InfoCodes::FUNC_EXECUTED;
        logger.getLog(infoCode, "Kchrec reconstruction");

        KchRec_main(chain, eventAnalysis, dataTypeOpt);
        break;
      }
      case Controls::MainMenu::KCHREC_BOOST:
        break;
      case Controls::MainMenu::IP_EV_BY_EV:
        break;
      case Controls::MainMenu::OMEGA_REC:
      {
        infoCode = ErrorHandling::InfoCodes::FUNC_EXECUTED;
        logger.getLog(infoCode, "Omega-pi0 Reconstruction");

        OmegaRec_main(chain, eventAnalysis, dataTypeOpt);
        break;
      }
      case Controls::MainMenu::REGEN_REJ:
      {
        // Regen_main(chain);
        break;
      }
      case Controls::MainMenu::KNEREC_TRILAT:
      {
        Neutrec_main(chain, eventAnalysis, dataTypeOpt);
        break;
      }
      case Controls::MainMenu::KNEREC_TRIANGLE:
        break;
      case Controls::MainMenu::EFF_SIG_TO_BCG:
        break;
      case Controls::MainMenu::KIN_FITS:
        break;
      case Controls::MainMenu::TRANSF_TO_CM:
        break;
      case Controls::MainMenu::CPV_NORM:
      {
        infoCode = ErrorHandling::InfoCodes::FUNC_EXECUTED;
        logger.getLog(infoCode, "CP Final Fit");

        CPFit_main(chain, eventAnalysis, dataTypeOpt);
        break;
      }
      case Controls::MainMenu::PLOTS:
        break;
      case Controls::MainMenu::COVMATRIX:
      {
        infoCode = ErrorHandling::InfoCodes::FUNC_EXECUTED;
        logger.getLog(infoCode, "Covariance Matrix Determination");

        CovMatrix_main(chain, eventAnalysis, dataTypeOpt);
        break;
      }
        // Plots_main(chain);
      }
    }
    else if (csFlag == 2)
    {
      switch (mainMenuControlSampleOpt)
      {
      case Controls::MainMenuControlSample::KCHREC:
      {
        infoCode = ErrorHandling::InfoCodes::FUNC_EXECUTED;
        logger.getLog(infoCode, "Kchrec reconstruction");

        KchRec_main(chain, eventAnalysis, dataTypeOpt);
        break;
      }
      case Controls::MainMenuControlSample::COV_MATRIX:
      {
        infoCode = ErrorHandling::InfoCodes::FUNC_EXECUTED;
        logger.getLog(infoCode, "Covariance Matrix Determination");

        CovMatrix_main(chain, eventAnalysis, dataTypeOpt);
        break;
      }
      case Controls::MainMenuControlSample::CORR_FACTOR:
        break;
      }
    }

  } while (mainMenuOpt != Controls::MainMenu::EXIT && mainMenuControlSampleOpt != Controls::MainMenuControlSample::EXIT);

  logger.printErrStats();

  return 0;
}
