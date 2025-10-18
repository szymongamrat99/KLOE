#ifndef GENVARS_H
#define GENVARS_H

#include <iostream>
#include <fstream>
#include <vector>

#include <TTree.h>
#include <TFile.h>
#include <TChain.h>

#include <const.h>
#include <ErrorLogs.h>
#include <MainMenu.h>
#include <fort_common.h>
#include <kloe_class.h>

Int_t GenVars(TChain &chain, Controls::DataType &data_type, ErrorHandling::ErrorLogs &logger, KLOE::pm00 &Obj);
Int_t split_channels(TChain &chain, Controls::DataType &data_type, ErrorHandling::ErrorLogs &logger, KLOE::pm00 &Obj);

#endif //! GENVARS_H