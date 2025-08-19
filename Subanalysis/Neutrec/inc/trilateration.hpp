#ifndef TRILATERATION_H
#define TRILATERATION_H

#include <iostream>
#include <fstream>
#include <vector>
#include <boost/optional.hpp>

#include <const.h>
#include <ErrorLogs.h>
#include <MainMenu.h>
#include <kloe_class.h>

void tri_neurec(TChain &chain, Controls::DataType &dataType, ErrorHandling::ErrorLogs &logger, KLOE::pm00 &Obj);
Int_t TriangleNeurec(TChain &chain, Controls::DataType &dataType, ErrorHandling::ErrorLogs &logger, KLOE::pm00 &Obj);
Int_t CompOfMethods(TChain &chain, Controls::DataType &dataType, ErrorHandling::ErrorLogs &logger, KLOE::pm00 &Obj);
Int_t TrilaterationNeurecKinfit(TChain &chain, Controls::DataType &dataType, ErrorHandling::ErrorLogs &logger, KLOE::pm00 &Obj);
Double_t trilateration_chi_square(const Double_t *);

ErrorHandling::ErrorCodes TrilaterationKinFit(Int_t N_free, Int_t N_const, Int_t M, Int_t loopcount, Float_t chiSqrStep, Int_t jmin, Int_t jmax, Int_t nclu, std::vector<Float_t> cluster[5], std::vector<Int_t> Asscl, std::vector<Float_t> bhabha_mom, std::vector<Float_t> bhabha_mom_err, std::vector<Float_t> bhabha_vtx, Int_t &bunchnum, std::vector<Float_t> &iptri_kinfit, std::vector<Int_t> &g4takentri_kinfit, std::vector<Float_t> gamma_mom_final[4], std::vector<Float_t> &fourKnetri_kinfit, std::vector<Float_t> &neu_vtx_min, ErrorHandling::ErrorLogs &logger);

#endif //! TRILATERATION_H