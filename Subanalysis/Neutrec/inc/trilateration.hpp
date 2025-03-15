#ifndef TRILATERATION_H
#define TRILATERATION_H

#include <iostream>
#include <fstream>
#include <vector>

#include <const.h>
#include <ErrorLogs.h>
#include <MainMenu.h>
#include <kloe_class.h>

void tri_neurec(TChain &chain, Controls::DataType &dataType, ErrorHandling::ErrorLogs &logger, KLOE::pm00 &Obj);
Int_t TriangleNeurec(TChain &chain, Controls::DataType &dataType, ErrorHandling::ErrorLogs &logger, KLOE::pm00 &Obj);
Int_t CompOfMethods(TChain &chain, Controls::DataType &dataType, ErrorHandling::ErrorLogs &logger, KLOE::pm00 &Obj);
Int_t TrilaterationNeurecKinfit(TChain &chain, Controls::DataType &dataType, ErrorHandling::ErrorLogs &logger, KLOE::pm00 &Obj);
Double_t trilateration_chi_square(const Double_t *);

#endif //! TRILATERATION_H