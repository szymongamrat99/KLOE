#include "MainMenuHandler.h"
#include <iostream>
#include <limits>

void MainMenuHandler::runMenuLoop(
    int csFlag,
    Controls::Menu& mainMenu,
    Controls::MainMenu& mainMenuOpt,
    Controls::Menu& mainMenuControlSample,
    Controls::MainMenuControlSample& mainMenuControlSampleOpt,
    TChain& chain,
    KLOE::pm00& eventAnalysis,
    Controls::DataType dataTypeOpt,
    ErrorHandling::ErrorLogs& logger,
    ErrorHandling::InfoCodes& infoCode
    // ConfigWatcher& cfgWatcher
) {
    ConfigWatcher cfgWatcher(Paths::propName);

    bool dataTypeErr = false, menuRangeErr = false;
    do
    {
        mainMenu.InitMenu();
        switch (csFlag)
        {
        case 1:
            std::cout << "1. Generated variables." << std::endl;
            std::cout << "2. K->pi+pi- reconstruction - no boost." << std::endl;
            std::cout << "3. K->pi+pi- reconstruction - boost." << std::endl;
            std::cout << "4. IP event by event - using K->pi+pi-." << std::endl;
            std::cout << "5. Omega-pi0 reconstruction." << std::endl;
            std::cout << "6. Regeneration rejection." << std::endl;
            std::cout << "7. K->pi0pi0 - using trilateration." << std::endl;
            std::cout << "8. K->pi0pi0 - using triangle, photons with trilateration." << std::endl;
            std::cout << "9. Efficiency and signal-to-background plots." << std::endl;
            std::cout << "10. Kinematic fits." << std::endl;
            std::cout << "11. Recalculation of variables to CM frame." << std::endl;
            std::cout << "12. CPV normalization." << std::endl;
            std::cout << "13. Plots." << std::endl;
            std::cout << "14. Covariance Matrix Determination." << std::endl;
            std::cout << "15. Exit." << std::endl;
            break;
        case 2:
            mainMenuControlSample.ShowOpt();
            break;
        }
        mainMenu.EndMenu();

        try
        {
            if (csFlag == 1)
            {
                std::cin >> mainMenuOpt;
                dataTypeErr = !std::cin;
                menuRangeErr = mainMenuOpt < Controls::MainMenu::GEN_VARS || mainMenuOpt > Controls::MainMenu::EXIT;
                if (dataTypeErr)
                    throw ErrorHandling::ErrorCodes::DATA_TYPE;
                else if (menuRangeErr)
                    throw ErrorHandling::ErrorCodes::RANGE;
            }
            else if (csFlag == 2)
            {
                std::cin >> mainMenuControlSampleOpt;
                dataTypeErr = !std::cin;
                menuRangeErr = mainMenuControlSampleOpt < Controls::MainMenuControlSample::COV_MATRIX || mainMenuControlSampleOpt > Controls::MainMenuControlSample::EXIT;
                if (dataTypeErr)
                    throw ErrorHandling::ErrorCodes::DATA_TYPE;
                else if (menuRangeErr)
                    throw ErrorHandling::ErrorCodes::RANGE;
            }
        }
        catch (ErrorHandling::ErrorCodes err)
        {
            std::cin.clear();
            std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            logger.getErrLog(err);
        }

        if (csFlag == 1)
        {
            switch (mainMenuOpt)
            {
            case Controls::MainMenu::GEN_VARS:
                infoCode = ErrorHandling::InfoCodes::FUNC_EXECUTED;
                logger.getLog(infoCode, "MC generated variables for analysis");
                GenVars_main(chain, eventAnalysis, dataTypeOpt);
                break;
            case Controls::MainMenu::KCHREC_NO_BOOST:
                infoCode = ErrorHandling::InfoCodes::FUNC_EXECUTED;
                logger.getLog(infoCode, "Kchrec reconstruction");
                KchRec_main(chain, eventAnalysis, dataTypeOpt);
                break;
            case Controls::MainMenu::KCHREC_BOOST:
                break;
            case Controls::MainMenu::IP_EV_BY_EV:
                break;
            case Controls::MainMenu::OMEGA_REC:
                infoCode = ErrorHandling::InfoCodes::FUNC_EXECUTED;
                logger.getLog(infoCode, "Omega-pi0 Reconstruction");
                OmegaRec_main(chain, eventAnalysis, dataTypeOpt);
                break;
            case Controls::MainMenu::REGEN_REJ:
                // Regen_main(chain);
                break;
            case Controls::MainMenu::KNEREC_TRILAT:
                Neutrec_main(chain, eventAnalysis, dataTypeOpt);
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
                infoCode = ErrorHandling::InfoCodes::FUNC_EXECUTED;
                logger.getLog(infoCode, "CP Final Fit");
                CPFit_main(chain, eventAnalysis, cfgWatcher, dataTypeOpt);
                break;
            case Controls::MainMenu::PLOTS:
                break;
            case Controls::MainMenu::COVMATRIX:
                infoCode = ErrorHandling::InfoCodes::FUNC_EXECUTED;
                logger.getLog(infoCode, "Covariance Matrix Determination");
                CovMatrix_main(chain, eventAnalysis, dataTypeOpt);
                break;
            }
        }
        else if (csFlag == 2)
        {
            switch (mainMenuControlSampleOpt)
            {
            case Controls::MainMenuControlSample::KCHREC:
                infoCode = ErrorHandling::InfoCodes::FUNC_EXECUTED;
                logger.getLog(infoCode, "Kchrec reconstruction");
                KchRec_main(chain, eventAnalysis, dataTypeOpt);
                break;
            case Controls::MainMenuControlSample::COV_MATRIX:
                infoCode = ErrorHandling::InfoCodes::FUNC_EXECUTED;
                logger.getLog(infoCode, "Covariance Matrix Determination");
                CovMatrix_main(chain, eventAnalysis, dataTypeOpt);
                break;
            case Controls::MainMenuControlSample::CORR_FACTOR:
                break;
            }
        }
    } while (mainMenuOpt != Controls::MainMenu::EXIT && mainMenuControlSampleOpt != Controls::MainMenuControlSample::EXIT);
}
