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

  return 0;
}