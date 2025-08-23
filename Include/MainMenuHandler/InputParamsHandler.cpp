#include "InputParamsHandler.h"
#include <iostream>
#include <fstream>
#include <limits>
#include <const.h>
#include <MainMenu.h>

InputParamsHandler::Params InputParamsHandler::getParams(
    nlohmann::json& properties,
    int lastFileMax,
    const std::string& propName,
    ErrorHandling::ErrorLogs& logger,
    Controls::Menu* dataType
) {
    Params params;
    bool firstFileRangeErr = false, lastFileRangeErr = false, dataTypeErr = false;

    // Tryb pomiar/kontrola
    while (true) {
        std::cout << "Measurement (1) / Control Sample (2)? ";
        std::cin >> params.csFlag;
        std::cout << std::endl;
        dataTypeErr = !std::cin;
        if (dataTypeErr) {
            auto err = ErrorHandling::ErrorCodes::DATA_TYPE;
            logger.getErrLog(err);
            std::cin.clear();
            std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            continue;
        }
        if (params.csFlag != 1 && params.csFlag != 2) {
            auto err = ErrorHandling::ErrorCodes::RANGE;
            logger.getErrLog(err);
            continue;
        }
        break;
    }

    // Pierwszy plik
    while (true) {
        std::cout << "Choose the first file: ";
        std::cin >> params.firstFile;
        std::cout << std::endl;
        dataTypeErr = !std::cin;
        firstFileRangeErr = params.firstFile < 1 || params.firstFile > static_cast<UInt_t>(lastFileMax);
        if (dataTypeErr) {
            auto err = ErrorHandling::ErrorCodes::DATA_TYPE;
            logger.getErrLog(err);
            std::cin.clear();
            std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            continue;
        }
        if (firstFileRangeErr) {
            auto err = ErrorHandling::ErrorCodes::RANGE;
            logger.getErrLog(err);
            continue;
        }
        break;
    }

    // Ostatni plik
    while (true) {
        std::cout << "Choose the last file: ";
        std::cin >> params.lastFile;
        std::cout << std::endl;
        dataTypeErr = !std::cin;
        lastFileRangeErr = params.lastFile < params.firstFile || params.lastFile > static_cast<UInt_t>(lastFileMax);
        if (dataTypeErr) {
            auto err = ErrorHandling::ErrorCodes::DATA_TYPE;
            logger.getErrLog(err);
            std::cin.clear();
            std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            continue;
        }
        if (lastFileRangeErr) {
            auto err = ErrorHandling::ErrorCodes::RANGE;
            logger.getErrLog(err);
            continue;
        }
        break;
    }

    // Zapis do properties
    properties["variables"]["rootFiles"]["Data"]["firstFile"] = params.firstFile;
    properties["variables"]["rootFiles"]["Data"]["lastFile"] = params.lastFile;
    // properties["variables"]["rootFiles"]["MC"]["firstFile"] = params.firstFile;
    // properties["variables"]["rootFiles"]["MC"]["lastFile"] = params.lastFile;

    std::ofstream outfile(propName);
    if (outfile.is_open()) {
        outfile << properties.dump(4);
        outfile.close();
    } else {
        auto err = ErrorHandling::ErrorCodes::FILE_NOT_EXIST;
        logger.getErrLog(err);
    }

    // Wybór typu danych
    while (true) {
        dataType->InitMenu();
        dataType->ShowOpt();
        dataType->EndMenu();
        std::cin >> params.dataTypeOpt;
        if (!std::cin) {
            auto err = ErrorHandling::ErrorCodes::DATA_TYPE;
            logger.getErrLog(err);
            std::cin.clear();
            std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            continue;
        }
        if (params.dataTypeOpt < Controls::DataType::SIGNAL_TOT || params.dataTypeOpt > Controls::DataType::DATA_ONLY) {
            auto err = ErrorHandling::ErrorCodes::MENU_RANGE;
            logger.getErrLog(err);
            continue;
        }
        break;
    }

    return params;
}
