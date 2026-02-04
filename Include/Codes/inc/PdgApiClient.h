#pragma once
#include <string>
#include <curl/curl.h>

class PdgApiClient {
public:
    PdgApiClient();
    ~PdgApiClient();

    // Metoda budująca URL i wykonująca zapytanie
    std::string fetchParticleJson(const std::string& pdgid, int year);

private:
    static size_t WriteCallback(void* contents, size_t size, size_t nmemb, std::string* s);
};