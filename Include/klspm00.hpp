#ifndef KLSPM00_H
#define KLSPM00_H

#include <chrono>
#include <iostream>
#include <string>

#include <const.h>

#include <MainMenu.h>
#include <ErrorLogs.h>
#include <kloe_class.h>
#include <ConfigWatcher.h>
#include <PhysicsConstants.h>


int CPFit_main(TChain &chain, KLOE::pm00 &Obj, ConfigWatcher &cfgWatcher, Controls::DataType &dataTypeOpt, PhysicsConstants &physConst);
int GenVars_main(TChain &chain, KLOE::pm00 &Obj, Controls::DataType &dataTypeOpt, PhysicsConstants &physConst);
int KchRec_main(TChain &chain, KLOE::pm00 &Obj, Controls::DataType &dataTypeOpt, PhysicsConstants &physConst);
int OmegaRec_main(TChain &chain, KLOE::pm00 &Obj, Controls::DataType &dataTypeOpt, PhysicsConstants &physConst);
// int Regen_main(TChain &chain);
int Neutrec_main(TChain &chain, KLOE::pm00 &Obj, Controls::DataType &dataTypeOpt, PhysicsConstants &physConst);
int CovMatrix_main(TChain &chain, KLOE::pm00 &Obj, Controls::DataType &dataTypeOpt, PhysicsConstants &physConst);
// int Plots_main(TChain &chain);

#endif //! KLSPM00_H