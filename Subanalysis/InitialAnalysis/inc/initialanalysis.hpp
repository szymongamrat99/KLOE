#ifndef COVARIANCEMATRIX_H
#define COVARIANCEMATRIX_H

#include <iostream>
#include <fstream>
#include <vector>

#include <TVectorT.h>
#include <TMatrixT.h>
#include <TChain.h>

#include <const.h>
#include <ErrorLogs.h>
#include <Logs.h>
#include <MainMenu.h>
#include <kloe_class.h>
#include <MeasQualityGraph.h>
#include <MomentumSmearing.h>
#include <CustomGraph.h>
#include <interference.h>
#include <fort_common.h>

int InitialAnalysis_full(TChain &chain, Controls::FileType &fileTypeOpt, ErrorHandling::ErrorLogs &logger, KLOE::pm00 &Obj, bool singleFile = false, std::string jobNumber = "");

void MctruthCounter(Int_t mctruth, UInt_t mctruth_num[8]);

#endif //! COVARIANCEMATRIX_H