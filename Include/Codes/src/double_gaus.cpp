#include "double_gaus.h"
#include <iostream>
#include <cmath>

/// Fit funkcja: 2 gaussy ze wspólną średnią
/// Parametry: [A1, sigma1, A2, sigma2, mean]
/// - A1, A2: amplitudy
/// - sigma1, sigma2: szerokości
/// - mean: wspólna średnia
Double_t double_gaus(Double_t *x, Double_t *p)
{
  Double_t norm1, norm2, gaus1, gaus2, value = 0.;
  
  // Parametry: A1, sigma1, A2, sigma2, mean
  Double_t A1 = p[0];
  Double_t sigma1 = p[1];
  Double_t A2 = p[2];
  Double_t sigma2 = p[3];
  Double_t mean = p[4];
  
  // Normalizacja każdego gaussa (z czynnikiem amplitudy)
  norm1 = A1 / (sigma1 * std::sqrt(2.0 * M_PI));
  norm2 = A2 / (sigma2 * std::sqrt(2.0 * M_PI));
  
  // Gaussy ze wspólną średnią
  Double_t dx1 = (x[0] - mean) / sigma1;
  Double_t dx2 = (x[0] - mean) / sigma2;
  
  gaus1 = norm1 * std::exp(-0.5 * dx1 * dx1);
  gaus2 = norm2 * std::exp(-0.5 * dx2 * dx2);
  
  value = gaus1 + gaus2;
  
  return value;
}

/// Pojedynczy gauss - helper
Double_t single_gaus_core(Double_t *x, Double_t *p)
{
  Double_t norm, gaus;
  
  // Parametry: A, sigma, mean
  Double_t A = p[0];
  Double_t sigma = p[1];
  Double_t mean = p[2];
  
  norm = A / (sigma * std::sqrt(2.0 * M_PI));
  Double_t dx = (x[0] - mean) / sigma;
  gaus = norm * std::exp(-0.5 * dx * dx);
  
  return gaus;
}

/// Średnia kombinowana dla 2 gaussów
/// p = [A1, sigma1, A2, sigma2, mean]
Double_t comb_mean_double(const Double_t *p, const Double_t *err)
{
  // Dla wspólnej średniej - zwracamy tę średnią
  return p[4];
}

/// Błąd średniej kombinowanej dla 2 gaussów
Double_t comb_mean_err_double(const Double_t *p, const Double_t *err)
{
  // Błąd wspólnej średniej
  // err[4] to błąd parametru mean
  return err[4];
}

/// Sigma kombinowana dla 2 gaussów (ważona amplitudami)
/// p = [A1, sigma1, A2, sigma2, mean]
Double_t comb_std_dev_double(const Double_t *p, const Double_t *err)
{
  Double_t A1 = p[0];
  Double_t sigma1 = p[1];
  Double_t A2 = p[2];
  Double_t sigma2 = p[3];
  
  Double_t Atot = A1 + A2;
  if (Atot <= 0.0)
    return 0.0;
  
  // Ważona średnia sigma (weighted average of sigmas)
  // sigma_comb = (A1*sigma1 + A2*sigma2) / (A1 + A2)
  Double_t sigma_comb = (A1 * sigma1 + A2 * sigma2) / Atot;
  
  return sigma_comb;
}

/// Błąd sigmy kombinowanej dla 2 gaussów
/// Propagacja błędu: d(sigma_comb)/d(A1), d(sigma_comb)/d(sigma1) itd.
Double_t comb_std_dev_err_double(const Double_t *p, const Double_t *err)
{
  Double_t A1 = p[0];
  Double_t sigma1 = p[1];
  Double_t A2 = p[2];
  Double_t sigma2 = p[3];
  
  Double_t A1_err = err[0];
  Double_t sigma1_err = err[1];
  Double_t A2_err = err[2];
  Double_t sigma2_err = err[3];
  
  Double_t Atot = A1 + A2;
  if (Atot <= 0.0)
    return 0.0;
  
  // sigma_comb = (A1*sigma1 + A2*sigma2) / Atot
  Double_t sigma_comb = (A1 * sigma1 + A2 * sigma2) / Atot;
  
  // Pochodne cząstkowe:
  // d(sigma_comb)/d(A1) = (sigma1 - sigma_comb) / Atot
  // d(sigma_comb)/d(sigma1) = A1 / Atot
  // d(sigma_comb)/d(A2) = (sigma2 - sigma_comb) / Atot
  // d(sigma_comb)/d(sigma2) = A2 / Atot
  
  Double_t d_dA1 = (sigma1 - sigma_comb) / Atot;
  Double_t d_dsigma1 = A1 / Atot;
  Double_t d_dA2 = (sigma2 - sigma_comb) / Atot;
  Double_t d_dsigma2 = A2 / Atot;
  
  Double_t err_squared = 0.0;
  err_squared += d_dA1 * d_dA1 * A1_err * A1_err;
  err_squared += d_dsigma1 * d_dsigma1 * sigma1_err * sigma1_err;
  err_squared += d_dA2 * d_dA2 * A2_err * A2_err;
  err_squared += d_dsigma2 * d_dsigma2 * sigma2_err * sigma2_err;
  
  return std::sqrt(err_squared);
}

/// Sigma węższego rozkładu (core gaussian)
/// Zwraca min(sigma1, sigma2)
Double_t get_core_sigma(const Double_t *p)
{
  Double_t sigma1 = p[1];
  Double_t sigma2 = p[3];
  
  return std::min(sigma1, sigma2);
}

/// Błąd sigmy core gaussu
/// Wybieramy błąd parametru sigma1 lub sigma2 (ten węższy)
Double_t get_core_sigma_err(const Double_t *p, const Double_t *err)
{
  Double_t sigma1 = p[1];
  Double_t sigma2 = p[3];
  Double_t sigma1_err = err[1];
  Double_t sigma2_err = err[3];
  
  // Jeśli sigma1 < sigma2, to core sigma to sigma1
  if (sigma1 < sigma2)
    return sigma1_err;
  else
    return sigma2_err;
}
