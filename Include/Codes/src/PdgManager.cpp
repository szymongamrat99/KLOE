#include <PdgManager.h>
#include <iostream>
#include <boost/optional/optional_io.hpp>
#include <fstream>
#include <boost/filesystem.hpp>

PdgManager::PdgManager()
{
  curl_global_init(CURL_GLOBAL_ALL);
}

PdgManager::~PdgManager()
{
  curl_global_cleanup();
}

PDGProvider::SummaryEdition &PdgManager::getParticleData(const std::string &pdgid, int year, const std::string &cachePath)
{
  std::string key = pdgid + "_" + std::to_string(year),
              filename = cachePath + "/pdg_cache_" + key + ".json";

  if(!boost::filesystem::exists(filename))
  {
    std::ofstream cacheFileOut(filename);

    if(!cacheFileOut.is_open())
    {
      std::cout << "Cache file cannot be created: " << filename << std::endl;
      return cache[key]; // Zwracamy pusty obiekt, który później i tak będzie uzupełniony
    }

    std::string jsonStr = client.fetchParticleJson(pdgid, year);

    cacheFileOut << jsonStr;
    cacheFileOut.flush();
    cacheFileOut.close();
  }

  // Jeśli nie ma w cache, wczytujemy z pliku (jeśli jest) i parsujemy do obiektu SummaryEdition
  std::ifstream cacheFileIn(filename);
  if (cache.find(key) == cache.end())
  {
    nlohmann::json jsonData;
    cacheFileIn >> jsonData;
    cache[key] = jsonData.get<PDGProvider::SummaryEdition>();
  }
  cacheFileIn.close();

  return cache[key];
}

double PdgManager::getBestValue(const PDGProvider::Property &prop, double hardcodedFallback, CPTStatus CPTOrNotCPT)
{
  static std::map<CPTStatus, TString> CPTStatusMap = {
      {CPTStatus::CPT, "Assuming CPT"},
      {CPTStatus::NON_CPT, "Not assuming CPT"},
      {CPTStatus::UNDEFINED, "CPT status undefined"}};

  if (!prop.get_pdg_values() || prop.get_pdg_values()->empty())
  {
    return hardcodedFallback;
  }
  const auto values = *prop.get_pdg_values();

  // 1. Szukamy OUR AVERAGE
  for (const auto &v : values)
  {
    if (v.get_type() == PDGProvider::Type::OUR_AVERAGE)
    {
      std::cout << v.get_comment().value_or((std::string)CPTStatusMap[CPTStatus::UNDEFINED]) << std::endl;

      if (v.get_comment().value_or((std::string)CPTStatusMap[CPTStatus::UNDEFINED]) == CPTStatusMap[CPTOrNotCPT])
      {
        return v.get_value().value();
      }
    }
  }

  // 2. Jeśli nie ma, szukamy OUR FIT
  for (const auto &v : values)
  {
    if (v.get_type() == PDGProvider::Type::OUR_FIT)
    {
      if (v.get_comment().value_or((std::string)CPTStatusMap[CPTStatus::UNDEFINED]) == CPTStatusMap[CPTOrNotCPT])
      {
        return v.get_value().value();
      }
    }
  }

  // 3. Jeśli nie znaleźliśmy nic pasującego w JSON, leci fallback
  return hardcodedFallback;
}