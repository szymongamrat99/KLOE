#if !defined(CONSTRAINTS_TRI_H)
#define CONSTRAINTS_TRI_H

#include <neutral_mom.h>
#include <lorentz_transf.h>
#include <plane_intersection.h>
#include <closest_approach.h>

#include <const.h>

Double_t ene_consv(Double_t *, Double_t *);
Double_t minv_consv(Double_t *, Double_t *);
Double_t x_consv(Double_t *, Double_t *);
Double_t y_consv(Double_t *, Double_t *);
Double_t z_consv(Double_t *, Double_t *);
Double_t gamma1_consv(Double_t *, Double_t *);
Double_t gamma2_consv(Double_t *, Double_t *);
Double_t gamma3_consv(Double_t *, Double_t *);
Double_t gamma4_consv(Double_t *, Double_t *);


#endif // CONSTRAINTS_TRI_H
