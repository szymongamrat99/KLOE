#ifndef INPUT_PARAMS_HANDLER_H
#define INPUT_PARAMS_HANDLER_H

#include <string>
#include "ErrorLogs.h"
#include "MainMenu.h"
#include <nlohmann/json.hpp>

class InputParamsHandler {
public:
    struct Params {
        UInt_t firstFile;
        UInt_t lastFile;
        UInt_t csFlag;
        Controls::DataType dataTypeOpt;
    };

    static Params getParams(
        nlohmann::json& properties,
        int lastFileMax,
        const std::string& propName,
        ErrorHandling::ErrorLogs& logger,
        Controls::Menu* dataType
    );
};

#endif // INPUT_PARAMS_HANDLER_H
