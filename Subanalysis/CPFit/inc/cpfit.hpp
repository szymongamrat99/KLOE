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

const TString cpfit_res_dir = cpfit_dir + result_dir;

int cp_fit_mc_data(TChain &chain, TString mode, Bool_t check_corr, Int_t loopcount, Int_t M, Int_t jmin, Int_t jmax, ErrorHandling::ErrorLogs &logger, Controls::DataType &dataType);

#endif //! CPFIT_H