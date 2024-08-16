#ifndef TRILATERATION_H
#define TRILATERATION_H

void tri_neurec(int, int, int);
void triangle_neurec(int, int, int);
void tri_neurec_kinfit(int, int);
void tri_neurec_kinfit_corr(Short_t, Int_t, Int_t, Short_t, Short_t, Short_t);
Double_t trilateration_chi_square(const Double_t *);

#endif //! TRILATERATION_H