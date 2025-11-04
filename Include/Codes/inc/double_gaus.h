#ifndef DOUBLE_GAUS_H
#define DOUBLE_GAUS_H

#include <TMath.h>

/// Fit funkcja: 2 gaussy ze wspólną średnią
/// Parametry: [A1, sigma1, A2, sigma2, mean]
Double_t double_gaus(Double_t *x, Double_t *p);

/// Pojedynczy gauss (helper)
Double_t single_gaus_core(Double_t *x, Double_t *p);

/// Średnia kombinowana dla 2 gaussów (ważona amplitudami)
Double_t comb_mean_double(const Double_t *p, const Double_t *err);

/// Błąd średniej kombinowanej dla 2 gaussów
Double_t comb_mean_err_double(const Double_t *p, const Double_t *err);

/// Sigma kombinowana dla 2 gaussów (ważona amplitudami)
Double_t comb_std_dev_double(const Double_t *p, const Double_t *err);

/// Błąd sigmy kombinowanej dla 2 gaussów
Double_t comb_std_dev_err_double(const Double_t *p, const Double_t *err);

/// Sigma węższego rozkładu (core gaussian)
Double_t get_core_sigma(const Double_t *p);

/// Błąd sigmy core gaussu (propagacja błędu)
Double_t get_core_sigma_err(const Double_t *p, const Double_t *err);

#endif // !DOUBLE_GAUS_H
