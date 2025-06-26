#include <ConfigManager.h>
#include <fstream>
#include <cstdlib>
#include <stdexcept>
#include <iostream>

nlohmann::json Config::properties_;
nlohmann::json Config::constants_;

void Config::load()
{
    std::string propPath = std::getenv("PROPERTIESKLOE") + std::string("/properties.json");
    std::ifstream pf(propPath);
    if (!pf)
        throw std::runtime_error("Brak pliku properties.json: " + propPath);
    pf >> properties_;

    std::string constPath = std::getenv("WORKDIR") + std::string("/scripts/Scripts/Subanalysis/Properties/pdg_api/pdg_const.json");
    std::ifstream cf(constPath);
    if (!cf)
        throw std::runtime_error("Brak pliku pdg_const.json: " + constPath);

    cf >> constants_;
}

nlohmann::json &Config::properties()
{
    static bool loaded = false;
    if (!loaded)
    {
        load();
        loaded = true;
    }
    return properties_;
}
nlohmann::json &Config::constants()
{
    static bool loaded = false;
    if (!loaded)
    {
        load();
        loaded = true;
    }
    return constants_;
}
void Config::reload()
{
    load();
}