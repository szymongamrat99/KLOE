#pragma once
#include <PdgSummaryEditionDataModel.h> // Quicktype
#include <PdgApiClient.h>
#include <map>
#include <memory>

class PdgManager {
public:
    // Singleton Access
    static PdgManager& getInstance() {
        static PdgManager instance;
        return instance;
    }

    // Główny punkt styku dla Twojej analizy
    PDGProvider::SummaryEdition& getParticleData(const std::string& pdgid, int year = 2025);

    // Nowa metoda logiczna
    static double getBestValue(const PDGProvider::Property& prop, double hardcodedFallback);

    // Zapobiegamy kopiowaniu
    PdgManager(const PdgManager&) = delete;
    void operator=(const PdgManager&) = delete;

private:
    PdgManager();
    ~PdgManager();

    PdgApiClient client;
    std::map<std::string, PDGProvider::SummaryEdition> cache;
};