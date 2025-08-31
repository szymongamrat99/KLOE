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
// int kchrec_KSKL(TChain &chain, Controls::DataType &dataType, ErrorHandling::ErrorLogs &logger, KLOE::pm00 &Obj);
// int kchrec_Closest(TChain &chain, Controls::DataType &dataType, ErrorHandling::ErrorLogs &logger, KLOE::pm00 &Obj);

ErrorHandling::ErrorCodes TwoBodyReconstruction(std::vector<Float_t> *Kchboost, std::vector<Float_t> *ip, std::vector<Float_t> *trk[2], KLOE::pm00 &Obj, std::vector<Float_t> &KchrecTwoBody, Float_t &gamma, TLorentzVector PiKaon4VecLAB[2], TLorentzVector trk4VecLAB[2]);


#endif //! KCHREC_H