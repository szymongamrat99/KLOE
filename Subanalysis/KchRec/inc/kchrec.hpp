#ifndef KCHREC_H
#define KCHREC_H

#include <iostream>
#include <fstream>
#include <vector>

#include <TBranch.h>
#include <TFile.h>
#include <TTree.h>

#include <charged_mom.h>
#include <ErrorLogs.h>
#include <MainMenu.h>
#include <kloe_class.h>

int kchrec_Kmass(TChain &chain, Controls::DataType &dataType, ErrorHandling::ErrorLogs &logger, KLOE::pm00 &Obj);
int kchrec_KSKL(TChain &chain, Controls::DataType &dataType, ErrorHandling::ErrorLogs &logger, KLOE::pm00 &Obj);
int kchrec_Closest(TChain &chain, Controls::DataType &dataType, ErrorHandling::ErrorLogs &logger, KLOE::pm00 &Obj);


#endif //! KCHREC_H