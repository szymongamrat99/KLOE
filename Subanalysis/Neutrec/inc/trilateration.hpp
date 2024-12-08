#ifndef TRILATERATION_H
#define TRILATERATION_H

#include <iostream>
#include <fstream>
#include <vector>

#include "../../../Include/const.h"
#include "../../../Include/Codes/ErrorLogs.h"
#include "../../../Include/Codes/Logs.h"
#include "../../../Include/Codes/MainMenu.h"

void tri_neurec(int, int, int, int);
int triangle_neurec(int, int, Controls::DataType);
int comp_of_methods(int, int, Controls::DataType);
void tri_neurec_kinfit(int, int, int);
void tri_neurec_kinfit_corr(Int_t, Int_t, Short_t, Short_t, Short_t, Controls::DataType);
Double_t trilateration_chi_square(const Double_t *);

#endif //! TRILATERATION_H