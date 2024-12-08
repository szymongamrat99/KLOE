#ifndef CPFIT_H
#define CPFIT_H

#include <iostream>
#include <fstream>
#include <vector>

#include "const.h"
#include "ErrorLogs.h"
#include "Logs.h"
#include "MainMenu.h"

const TString cpfit_res_dir = cpfit_dir + result_dir;

int cp_fit_mc_data(Int_t firstFile, Int_t lastFile, TString mode, Bool_t check_corr, Int_t loopcount, Int_t M, Int_t range);

#endif //! CPFIT_H