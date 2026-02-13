#include <PdgManager.h>
#include <iostream>
#include <boost/optional/optional_io.hpp>

PdgManager::PdgManager()
{
  curl_global_init(CURL_GLOBAL_ALL);
}

PdgManager::~PdgManager()
{
  curl_global_cleanup();
}

PDGProvider::SummaryEdition &PdgManager::getParticleData(const std::string &pdgid, int year)
{
  std::string key = pdgid + "_" + std::to_string(year);

  // Jeśli nie ma w cache, pobierz i sparsuj
  if (cache.find(key) == cache.end())
  {
    std::string jsonStr = client.fetchParticleJson(pdgid, year);
    cache[key] = nlohmann::json::parse(jsonStr).get<PDGProvider::SummaryEdition>();
  }

  return cache[key];
}

double PdgManager::getBestValue(const PDGProvider::Property &prop, double hardcodedFallback, CPTStatus CPTOrNotCPT)
{
  std::map<CPTStatus, TString> CPTStatusMap = {
      {CPTStatus::CPT, "Assuming CPT"},
      {CPTStatus::NON_CPT, "Not assuming CPT"},
      {CPTStatus::UNDEFINED, "CPT status undefined"}};

  if (!prop.get_pdg_values() || prop.get_pdg_values()->empty())
  {
    return hardcodedFallback;
  }
  const auto &values = *prop.get_pdg_values();

  // 1. Szukamy OUR AVERAGE
  for (const auto &v : values)
  {
    if (v.get_type() == PDGProvider::Type::OUR_AVERAGE)
    {
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