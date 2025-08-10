#ifndef MAIN_MENU_HANDLER_H
#define MAIN_MENU_HANDLER_H

#include <TChain.h>
#include "../Include/klspm00.hpp"
#include "ErrorLogs.h"

class MainMenuHandler {
public:
    static void runMenuLoop(
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
    );
};

#endif // MAIN_MENU_HANDLER_H
