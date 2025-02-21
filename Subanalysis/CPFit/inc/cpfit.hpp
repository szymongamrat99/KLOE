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
#include <interference.h>
#include <fort_common.h>
#include <lorentz_transf.h>

int cp_fit_mc_data(TChain &chain, TString mode, bool check_corr, Controls::DataType &data_type, ErrorHandling::ErrorLogs &logger, KLOE::pm00 &Obj);

#endif //! CPFIT_H