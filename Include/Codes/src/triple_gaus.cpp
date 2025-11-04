#include "triple_gaus.h"
#include <iostream>
#include <vector>

Double_t triple_gaus(Double_t *x, Double_t *p)
{
  Double_t norm[3], gaus[3];
  Double_t value = 0.;

  for (Int_t i = 0; i < 3; i++)
  {
    norm[i] = p[i * 3] / (p[2 + i * 3] * sqrt(2 * M_PI));
    gaus[i] = norm[i] * exp(-pow((x[0] - p[1 + i * 3]) / p[2 + i * 3], 2) / 2.);

    value += gaus[i];
  }

  return value;
}

Double_t single_gaus(Double_t *x, Double_t *p)
{
  Double_t norm, gaus;

  norm = p[0] / (p[2] * sqrt(2 * M_PI));
  gaus = norm * exp(-pow((x[0] - p[1]) / p[2], 2) / 2.);

  return gaus;
}

Double_t comb_mean(const Double_t *p, const Double_t *err)
{
  Double_t
      N[3],
      mean[3],
      sigma[3],
      rel_err_N[3],
      rel_err_mean[3],
      rel_err_sigma[3],
      tot_rel_err = 0.,
      NTot = 0.,
      nominator = 0.,
      std_dev = 0.;

  for (Int_t i = 0; i < 3; i++)
  {
    N[i] = p[i * 3];
    mean[i] = p[1 + i * 3];
    sigma[i] = p[2 + i * 3];

    NTot += N[i];
    nominator += N[i] * mean[i];
  }

  return nominator / NTot;
}

Double_t comb_mean_err(const Double_t *p, const Double_t *err)
{
  Double_t
      N[3],
      mean[3],
      sigma[3],
      N_err[3],
      mean_err[3],
      sigma_err[3],
      rel_err_N[3],
      rel_err_mean[3],
      rel_err_sigma[3],
      tot_rel_err = 0.,
      NTot = 0.,
      nominator = 0.,
      err_nominator = 0.,
      err_denominator = 0.,
      err_sigma_tot = 0.,
      err_N_tot = 0.,
      std_dev = 0.;
  Int_t
      perm[3][3] = {{0, 1, 2},
                    {1, 0, 2},
                    {2, 0, 1}};

  for (Int_t i = 0; i < 3; i++)
  {
    N[i] = p[i * 3];
    mean[i] = p[1 + i * 3];
    sigma[i] = p[2 + i * 3];

    N_err[i] = err[i * 3];
    mean_err[i] = err[1 + i * 3];
    sigma_err[i] = err[2 + i * 3];

    rel_err_N[i] = N_err[i] / N[i];
    rel_err_mean[i] = mean_err[i] / mean[i];
    rel_err_sigma[i] = sigma_err[i] / sigma[i];

    tot_rel_err = abs(rel_err_sigma[i]);

    // if (tot_rel_err > 0.9 || rel_err_N[i] > 0.9)
    //   N[i] = 0.;

    NTot += N[i];
    nominator += N[i] * mean[i];
  };

  // Weighted mean: mu_comb = sum(N_i * mu_i) / N_total
  Double_t mean_comb = nominator / NTot;

  // Propagacja błędu przez pochodne cząstkowe:
  // d(mu_comb)/d(mu_i) = N_i / N_tot
  // d(mu_comb)/d(N_i) = (mu_i - mu_comb) / N_tot

  Double_t err_squared = 0.;
  for (Int_t i = 0; i < 3; i++)
  {
    // Składowa od błędu mean_i
    Double_t d_dmean = N[i] / NTot;
    err_squared += pow(d_dmean * mean_err[i], 2);

    // Składowa od błędu N_i
    Double_t d_dN = (mean[i] - mean_comb) / NTot;
    err_squared += pow(d_dN * N_err[i], 2);
  }

  return sqrt(err_squared);
}

Double_t comb_std_dev(const Double_t *p, const Double_t *err)
{
  Double_t
      N[3],
      mean[3],
      sigma[3],
      rel_err_N[3],
      rel_err_mean[3],
      rel_err_sigma[3],
      tot_rel_err = 0.,
      NTot = 0.,
      nominator1 = 0.,
      nominator2 = 0.,
      nominator3 = 0.,
      nominator = 0.,
      std_dev = 0.,
      weight[3] = {0.};

  NTot = p[0] + p[3] + p[6];

  weight[0] = p[0] / NTot;
  weight[1] = p[3] / NTot;
  weight[2] = p[6] / NTot;

  for (Int_t i = 0; i < 3; i++)
  {
    N[i] = p[i * 3];
    mean[i] = p[1 + i * 3];
    sigma[i] = p[2 + i * 3];

    // ✅ POPRAWKA: weighted variance, nie (N*sigma)^2!
    // Stara formuła: pow(N[i] * sigma[i], 2) dawała zawyżone wartości
    nominator1 += weight[i] * sigma[i] * sigma[i];
    nominator2 += weight[i] * mean[i] * mean[i];
    nominator3 += weight[i] * mean[i];
  }

  nominator = nominator1 + nominator2 - nominator3 * nominator3;

  return sqrt(nominator);
}

Double_t comb_std_dev_err(const Double_t *p, const Double_t *err)
{
  Double_t
      N[3],
      mean[3],
      sigma[3],
      N_err[3],
      mean_err[3],
      sigma_err[3],
      rel_err_N[3],
      rel_err_mean[3],
      rel_err_sigma[3],
      tot_rel_err = 0.,
      NTot = 0.,
      nominator = 0.,
      nominator1 = 0.,
      nominator2 = 0.,
      nominator3 = 0.,
      err_nominator = 0.,
      err_denominator = 0.,
      err_sigma_tot = 0.,
      err_N_tot = 0.,
      std_dev = 0.,
      weight[3] = {0.};
  Int_t
      perm[3][3] = {{0, 1, 2},
                    {1, 0, 2},
                    {2, 0, 1}};

  NTot = p[0] + p[3] + p[6];

  weight[0] = p[0] / NTot;
  weight[1] = p[3] / NTot;
  weight[2] = p[6] / NTot;

  for (Int_t i = 0; i < 3; i++)
  {
    N[i] = p[i * 3];
    mean[i] = p[1 + i * 3];
    sigma[i] = p[2 + i * 3];

    N_err[i] = err[i * 3];
    mean_err[i] = err[1 + i * 3];
    sigma_err[i] = err[2 + i * 3];

    // ✅ POPRAWKA: zgodnie z nową formułą weighted variance
    nominator1 += weight[i] * sigma[i] * sigma[i];
    nominator2 += weight[i] * mean[i] * mean[i];
    nominator3 += weight[i] * mean[i];
  };

  nominator = nominator1 + nominator2 - nominator3 * nominator3;

  // Weighted RMS: sigma_comb = sqrt(sum(N_i * sigma_i^2) / N_total)
  std_dev = sqrt(nominator);

  // Propagacja błędu przez pochodne cząstkowe:
  // d(sigma_comb)/d(sigma_i) = (N_i * sigma_i) / (N_tot * sigma_comb)
  // d(sigma_comb)/d(N_i) = (sigma_i^2 - sigma_comb^2) / (2 * N_tot * sigma_comb)

  Double_t err_squared = 0.;
  for (Int_t i = 0; i < 3; i++)
  {
    // Składowa od błędu sigma_i
    Double_t d_dsigma = (2 * N[i] * sigma[i]) / (NTot);
    err_squared += pow(d_dsigma * sigma_err[i], 2);

    // Składowa od błędu mean_i
    Double_t d_dmean = (2 * N[i] / NTot ) * ( mean[i] - nominator3 );
    err_squared += pow(d_dmean * mean_err[i], 2);

    // Składowa od błędu N_i
    Double_t d_dN = pow(sigma[i], 2) + pow(mean[i], 2) - nominator1 - nominator2 - 2 * mean[i] * nominator3 + 2 * nominator3 * nominator3;

    d_dN /= (NTot);

    err_squared += pow(d_dN * N_err[i], 2);
  }

  return sqrt(err_squared);
}