#ifndef CPFIT_H
#define CPFIT_H

#include <iostream>
#include <fstream>
#include <vector>

#include "const.h"
#include "ErrorLogs.h"
#include "Logs.h"
#include "MainMenu.h"
#include <kloe_class.h>
#include <MeasQualityGraph.h>
#include <interference.h>
#include <fort_common.h>

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

int cp_fit_mc_data(TChain &chain, TString mode, bool check_corr, Controls::DataType &data_type, ErrorHandling::ErrorLogs &logger, KLOE::pm00 &Obj);

#endif //! CPFIT_H