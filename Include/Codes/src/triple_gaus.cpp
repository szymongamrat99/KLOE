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

    rel_err_N[i] = err[i * 3] / N[i];
    rel_err_mean[i] = err[1 + i * 3] / mean[i];
    rel_err_sigma[i] = err[2 + i * 3] / sigma[i];

    tot_rel_err = abs(rel_err_mean[i]);

    // if (tot_rel_err > 0.9 || rel_err_N[i] > 0.9)
    //   N[i] = 0.;

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

    nominator += pow(N[i] * sigma[i], 2);
    NTot += N[i];

    err_sigma_tot += pow(mean_err[i], 2) * pow(N[i], 2) * pow(NTot, 2);
    err_N_tot += pow(N_err[perm[i][0]], 2) * pow(N[perm[i][1]] * (mean[perm[i][0]] - mean[perm[i][1]]) + N[perm[i][2]] * (mean[perm[i][0]] - mean[perm[i][2]]), 2);
  };

  err_denominator = pow(NTot, 2);
  err_nominator = err_sigma_tot + err_N_tot;

  return sqrt(err_nominator) / err_denominator;
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
      nominator = 0.,
      std_dev = 0.;

  for (Int_t i = 0; i < 3; i++)
  {
    N[i] = p[i * 3];
    mean[i] = p[1 + i * 3];
    sigma[i] = p[2 + i * 3];

    rel_err_N[i] = err[i * 3] / N[i];
    rel_err_mean[i] = err[1 + i * 3] / mean[i];
    rel_err_sigma[i] = err[2 + i * 3] / sigma[i];

    tot_rel_err = abs(rel_err_sigma[i]);

    // if (tot_rel_err > 0.9 || rel_err_N[i] > 0.9)
    //   N[i] = 0.;

    NTot += N[i];
    nominator += pow(N[i] * sigma[i], 2);
  }

  return sqrt(nominator) / NTot;
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

    nominator += pow(N[i] * sigma[i], 2);
    NTot += N[i];

    err_sigma_tot += pow(sigma_err[i], 2) * pow(N[i], 4) * pow(NTot, 2) * pow(sigma[i], 2);
    err_N_tot += pow(N_err[perm[i][0]], 2) * pow(pow(N[perm[i][1]] * sigma[perm[i][1]], 2) + pow(N[perm[i][2]] * sigma[perm[i][2]], 2) - N[perm[i][0]] * (N[perm[i][1]] + N[perm[i][2]]) * pow(sigma[perm[i][0]], 2), 2);
  };

  err_denominator = pow(NTot, 2) * sqrt(nominator);
  err_nominator = err_sigma_tot + err_N_tot;

  return sqrt(err_nominator) / err_denominator;
}