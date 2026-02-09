#ifndef INTERF_FUNCTION_H
#define INTERF_FUNCTION_H

#include "TMath.h"

Double_t interf_function(const Double_t *x, const Double_t *par);

Double_t interf_function_00pm(const Double_t *x, const Double_t *par);
Double_t interf_function_pm00(const Double_t *x, const Double_t *par);
Double_t interf_function_pmpm(const Double_t *x, const Double_t *par);

Double_t interf_function_00pm_to_fit_mock(const Double_t *x, const Double_t *par);
Double_t interf_function_pm00_to_fit_mock(const Double_t *x, const Double_t *par);
Double_t interf_function_pmpm_to_fit_mock(const Double_t *x, const Double_t *par);


#endif