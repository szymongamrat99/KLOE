#ifndef KLSPM00_H
#define KLSPM00_H

#include <chrono>
#include <iostream>
#include <string>

#include <const.h>

#include <MainMenu.h>
#include <ErrorLogs.h>
#include <kloe_class.h>

int CPFit_main(TChain &chain, KLOE::pm00 &Obj, Controls::DataType &dataTypeOpt);
int GenVars_main(TChain &chain, KLOE::pm00 &Obj, Controls::DataType &dataTypeOpt);
// int KchRec_main(TChain &chain, KLOE::pm00 &Obj, Controls::DataType &dataTypeOpt);
int OmegaRec_main(TChain &chain, KLOE::pm00 &Obj, Controls::DataType &dataTypeOpt);
// int Regen_main(TChain &chain);
int Neutrec_main(TChain &chain, KLOE::pm00 &Obj, Controls::DataType &dataTypeOpt);
// int Plots_main(TChain &chain);

#endif //! KLSPM00_H