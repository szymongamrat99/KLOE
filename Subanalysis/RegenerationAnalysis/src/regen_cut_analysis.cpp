#include <iostream>
#include <fstream>

#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TH2.h>
#include <TFitResult.h>
#include <TGraphAsymmErrors.h>

#include "../inc/regen_analysis.h"

#include <fort_common.h>
#include <const.h>

#include "../inc/regenrejec.hpp"

using namespace std;

int regenrejec(TChain &chain, KLOE::pm00 &Obj, Controls::DataType &dataTypeOpt, ErrorHandling::ErrorLogs &logger)
{
  RegenAnalysis analysis(&chain);

  analysis.Loop();
  analysis.FitResults();

  return 0;
}