#ifndef UNCERTAINTIES_H
#define UNCERTAINTIES_H

#include "TMath.h"

Double_t clu_ene_error(Double_t);
Double_t clu_time_error(Double_t);
Double_t clu_x_error(Double_t, Double_t, Double_t, Double_t);
Double_t clu_y_error(Double_t, Double_t, Double_t, Double_t);
Double_t clu_z_error(Double_t, Double_t, Double_t, Double_t);

#endif