#ifndef INPUT_PARAMS_HANDLER_H
#define INPUT_PARAMS_HANDLER_H

#include <string>
#include "ErrorLogs.h"
#include "MainMenu.h"
#include <json.hpp>

class InputParamsHandler {
public:
    struct Params {
        UInt_t firstFile;
        UInt_t lastFile;
        UInt_t csFlag;
        Controls::DataType dataTypeOpt;
    };

    static Params getParams(
        nlohmann::json& Utils::properties,
        int lastFileMax,
        const std::string& Paths::propName,
        ErrorHandling::ErrorLogs& logger,
        Controls::Menu* dataType
    );
};

#endif // INPUT_PARAMS_HANDLER_H
