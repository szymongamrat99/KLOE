#ifndef OMEGAREC_HPP
#define OMEGAREC_HPP

#include <iostream>
#include <fstream>
#include <vector>
#include <boost/optional.hpp>

#include <TTree.h>
#include <TFile.h>

#include <const.h>
#include <ErrorLogs.h>
#include <MainMenu.h>
#include <clear_variables.h>
#include <kloe_class.h>
#include <charged_mom.h>

int omegarec(TChain &chain, Controls::DataType &dataType, ErrorHandling::ErrorLogs &logger, KLOE::pm00 &Obj);
int omegarec_kin_fit(TChain &chain, Controls::DataType &dataType, ErrorHandling::ErrorLogs &logger, KLOE::pm00 &Obj);
int plots(TChain &chain, Short_t &loopcount, Short_t &numOfConstraints, Short_t &jmin, Short_t &jmax, Controls::DataType &dataType, KLOE::pm00 &Obj, ErrorHandling::ErrorLogs &logger);

#endif //! OMEGAREC_HPP