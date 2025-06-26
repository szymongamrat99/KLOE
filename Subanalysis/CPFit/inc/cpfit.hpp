#ifndef CPFIT_H
#define CPFIT_H

#include <iostream>
#include <fstream>
#include <vector>

#include <TCut.h>

#include <const.h>
#include "ErrorLogs.h"
#include "Logs.h"
#include "MainMenu.h"
#include <kloe_class.h>
#include <MeasQualityGraph.h>
#include <charged_mom.h>
#include <MomentumSmearing.h>
#include <CustomGraph.h>
#include <interference.h>
#include <fort_common.h>
#include <ConfigWatcher.h>

/**
* @file cpfit.hpp
* @brief Header file for CP Fit functions with plots generation
* @param chain Address to the TChain used to get the data
* @param mode Mode of operation ("split" by default)
* @param check_corr Flag, if the efficiency correction should be used
* @param data_type Parameter, which type of data should be taken into account: 
                  * - 1: MC signal only 
                  * - 2: MC signal + background 
                  * - 3: MC + Data 
                  * - 4: Entire MC signal (before cuts and with all events with errors)
* @param logger Address to the logger for error and info.
* @param Obj Address to the object from KLOE::pm00 class
* @return Integer flag: 
          * - 0: no error
          * - 1: error
*/

int cp_fit_mc_data(TChain &chain, TString mode, bool check_corr, Controls::DataType &data_type, ErrorHandling::ErrorLogs &logger, KLOE::pm00 &Obj, ConfigWatcher &cfgWatcher);

int cp_fit_func(KLOE::interference &event, std::vector<std::vector<Double_t>> &relativeErr, std::vector<std::vector<Double_t>> &real, std::vector<std::vector<Double_t>> &imaginary, ErrorHandling::ErrorLogs &logger, KLOE::pm00 &Obj);

int cut_search(TChain &chain, TString mode, bool check_corr, Controls::DataType &data_type, ErrorHandling::ErrorLogs &logger, KLOE::pm00 &Obj);

struct Event
{
  Int_t mctruth;
  Int_t mcflag;

  Double_t DeltaTRec;
  Double_t DeltaTGen;

  // Regeneration peaks cuts
  Double_t beamPipeCutNeg;
  Double_t beamPipeCutPos;

  // Omega variables for cuts
  Double_t KinEnePi0;
  Double_t MinvOmega;
  Double_t LineDown;
  Double_t LineUp;

  // Omega fiducial volume
  Bool_t rho00;
  Bool_t rhopm;
  Bool_t z00;
  Bool_t zpm;
  Bool_t totalOmegaFidVol;

  // General cut variables
  Double_t trcv;
  Double_t minvKch;
  Double_t minvKne;
  Double_t Qmiss;
};

#endif //! CPFIT_H