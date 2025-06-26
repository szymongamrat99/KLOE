#ifndef UNCERTAINTIES_H
#define UNCERTAINTIES_H

#include "TMath.h"

Double_t clu_ene_error(Float_t);
Double_t clu_time_error(Float_t);
Double_t clu_x_error(Float_t, Float_t, Float_t, Float_t);
Double_t clu_y_error(Float_t, Float_t, Float_t, Float_t);
Double_t clu_z_error(Float_t, Float_t, Float_t, Float_t);

#endif