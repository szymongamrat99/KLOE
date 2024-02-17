#ifndef TRILATERATION_H
#define TRILATERATION_H

#include "../src/chain_init.C"

void tri_neurec(int, int, int);
void tri_neurec_kinfit(int, int);
void tri_neurec_kinfit_corr(int, int);
Double_t trilateration_chi_square(const Double_t *);

#endif //! TRILATERATION_H