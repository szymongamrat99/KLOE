#include "../inc/bttf_utils.h"

TChain *LoadPM00Chain(const TString &date1 = "2025-11-24", 
                      const TString &date2 = "2025-11-26")
{
  TString
      pm00_dir = Paths::initialanalysis_dir + Paths::root_files_dir + date1 + "/" + "mk0_initial_analysis_*_Signal*.root",
      pm00_dir_all_phys3 = Paths::initialanalysis_dir + Paths::root_files_dir + date2 + "/" + "mk0_initial_analysis_*_Signal*.root";

  TChain *chain_pm00 = new TChain("h1");
  chain_pm00->Add(pm00_dir);
  chain_pm00->Add(pm00_dir_all_phys3);

  return chain_pm00;
}

TChain *LoadPMPMChain(const TString &date = "2026-01-14")
{
  TString pmpm_dir = Paths::initialanalysis_dir + Paths::root_files_dir + date + "/" + "mk0_initial_analysis_*.root";

  TChain *chain_pmpm = new TChain("h1");
  chain_pmpm->Add(pmpm_dir);

  return chain_pmpm;
}

