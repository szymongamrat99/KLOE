#ifndef REGENREJEC_HPP
#define REGENREJEC_HPP

#include <iostream>
#include <fstream>
#include <vector>

#include <TTree.h>
#include <TFile.h>
#include <TChain.h>

#include <const.h>
#include <ErrorLogs.h>
#include <MainMenu.h>
#include <clear_variables.h>
#include <kloe_class.h>

int regenrejec(TChain &, KLOE::pm00 &, Controls::DataType &, ErrorHandling::ErrorLogs &);
// int plots(TChain &, KLOE::pm00 &, Controls::DataType &, ErrorHandling::ErrorLogs &);

#endif //! REGENREJEC_HPP