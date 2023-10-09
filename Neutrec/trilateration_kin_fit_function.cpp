#include "trilateration_kin_fit_function.h"

Double_t trilateration_chi_square(const Double_t *x)
{
  Double_t clusters[5][4], clusters_err[5][4];

  for(Int_t i = 0; i < 4; i++)
    for(Int_t j = 0; j < 5; j++)
    {
      clusters[j][i] = x[i * 5 + j];
      clusters_err[j][i] = x[i * 5 + j + 2*23];
    }

  return 0;
}