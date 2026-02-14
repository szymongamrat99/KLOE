#pragma once
#include <PdgSummaryEditionDataModel.h> // Quicktype
#include <PdgApiClient.h>
#include <map>
#include <memory>
#include <TString.h>

enum class CPTStatus {
    CPT,
    NON_CPT,
    UNDEFINED
};

class PdgManager {
public:
    // Singleton Access
    static PdgManager& getInstance() {
        static PdgManager instance;
        return instance;
    }

    // Główny punkt styku dla Twojej analizy
    PDGProvider::SummaryEdition& getParticleData(const std::string& pdgid, int year = 2025, const std::string &cachePath = "/tmp");

    // Nowa metoda logiczna
    static double getBestValue(const PDGProvider::Property& prop, double hardcodedFallback, CPTStatus CPTOrNotCPT = CPTStatus::UNDEFINED);

    // Zapobiegamy kopiowaniu
    PdgManager(const PdgManager&) = delete;
    void operator=(const PdgManager&) = delete;

private:
    PdgManager();
    ~PdgManager();

    PdgApiClient client;
    std::map<std::string, PDGProvider::SummaryEdition> cache;
};