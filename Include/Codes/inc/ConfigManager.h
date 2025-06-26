#pragma once

#include <json.hpp>
#include <string>

class Config {
public:
    static nlohmann::json& properties();
    static nlohmann::json& constants();
    static void reload();
private:
    static nlohmann::json properties_;
    static nlohmann::json constants_;
    static void load();
};