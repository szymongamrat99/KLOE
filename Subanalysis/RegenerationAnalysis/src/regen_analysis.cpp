#include <iostream>

#include <TChain.h>

#include <const.h>
#include <ErrorLogs.h>
#include <kloe_class.h>
#include "../inc/regen_frac_fit.h"

int regeneration_correction_factors(TChain &chain, KLOE::pm00 &Obj, Controls::DataType &dataTypeOpt, ErrorHandling::ErrorLogs &logger)
{
  gErrorIgnoreLevel = kBreak;

  TTreeReader *reader = new TTreeReader(&chain);
  RegenerationFractionFit fit(reader);

  fit.SaveHistograms();
  fit.DrawHistograms();

  delete reader;

  std::string weightsFilePath = (std::string)Paths::regen_analysis_dir + (std::string)Paths::result_dir + "regeneration_analysis_results.root";
  RegenerationFractionFit fit_test(weightsFilePath);

  std::cout << "Test regeneration weight (DC, radius=50 cm): " << fit_test.GetContinuousRegenerationWeight(30.0, 40.0) << std::endl;
  return 0;
}