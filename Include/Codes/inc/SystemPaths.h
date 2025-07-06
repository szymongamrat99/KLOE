// SystemPath.h
#pragma once

#include <string>

/**
 * @brief Provides string constants for accessing physical constants in JSON.
 * Using these constants prevents typos and enables easier refactoring.
 * All keys are defined using JSON Pointer syntax for clarity in nesting.
 */
namespace SystemPath
{
  static const std::string kloedataPath = getenv("KLOE_DBV26_DK0");
  static const std::string kloeMCPath = getenv("KLOE_DBV26_MK0");
  static const std::string workdirPath = getenv("WORKDIR");
  static const std::string pdgConstFilePath = workdirPath + "/KLOE/Subanalysis/Properties/pdg_api/pdg_const.json";
  static const std::string propertiesPath = getenv("PROPERTIESKLOE");
  static const std::string generalPropertiesPath = propertiesPath + "/properties.json";

  static const std::string basePath = workdirPath + "/KLOE/";
  
  static const std::string gen_vars_dir = basePath + "Subanalysis/GeneratedVars/";
  static const std::string neutrec_dir = basePath + "Subanalysis/Neutrec/";
  static const std::string cpfit_dir = basePath + "Subanalysis/CPFit/";
  static const std::string covmatrix_dir = basePath + "Subanalysis/CovarianceMatrix/";
  static const std::string omegarec_dir = basePath + "Subanalysis/OmegaRec/";
  static const std::string efficiency_dir = basePath + "Subanalysis/EfficiencyAnalysis/";
  static const std::string charged_dir = basePath + "Subanalysis/KchRec/";
  static const std::string plots_dir = basePath + "Subanalysis/Plots/";
  static const std::string root_files_dir = "root_files/";
  static const std::string input_dir = "input/";
  static const std::string logs_dir = "log/";
  static const std::string result_dir = "results/";
  static const std::string img_dir = "img/";
}