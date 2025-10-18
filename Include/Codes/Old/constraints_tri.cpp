#include <string.h>

#include <reconstructor.h>
#include <constraints_tri.h>

const Double_t Trf = 2.715; // ns
Bool_t cond_detector = kFALSE;

Int_t selected[4] = {1, 2, 3, 4};

Double_t ene_consv(Double_t *x, Double_t *p)
{
  Float_t boost_vec[3] = {-p[20] / p[23], -p[21] / p[23], -p[22] / p[23]};
  Float_t gamma_mom[4][4], vec_init[4], vec_end[4], neu_vtx[4];

  Double_t value[2], value_min, T0;

  T0 = p[27] * Trf;

  Reconstructor R;
  Solution S;

  // Setting clusters for a solution
  for (Int_t k = 0; k < 4; k++)
  {
    R.SetClu(k, p[k * 5],
             p[k * 5 + 1],
             p[k * 5 + 2],
             p[k * 5 + 3] + T0,
             p[k * 5 + 4]);

    R.SetClu(4, 0., 0., 0., 0., 0.);
    R.SetClu(5, 0., 0., 0., 0., 0.);
  }

  S = R.MySolve(selected);

  for (Int_t i = 0; i < 2; i++)
  {

    if (S.error[i] == 0)
    {
      neu_vtx[0] = S.sol[i][0];
      neu_vtx[1] = S.sol[i][1];
      neu_vtx[2] = S.sol[i][2];
      neu_vtx[3] = S.sol[i][3];

      neutral_mom(p[0], p[1], p[2], p[4], neu_vtx, gamma_mom[0]);
      neutral_mom(p[5], p[6], p[7], p[9], neu_vtx, gamma_mom[1]);
      neutral_mom(p[10], p[11], p[12], p[14], neu_vtx, gamma_mom[2]);
      neutral_mom(p[15], p[16], p[17], p[19], neu_vtx, gamma_mom[3]);

      vec_init[0] = gamma_mom[0][0] + gamma_mom[1][0] + gamma_mom[2][0] + gamma_mom[3][0];
      vec_init[1] = gamma_mom[0][1] + gamma_mom[1][1] + gamma_mom[2][1] + gamma_mom[3][1];
      vec_init[2] = gamma_mom[0][2] + gamma_mom[1][2] + gamma_mom[2][2] + gamma_mom[3][2];
      vec_init[3] = gamma_mom[0][3] + gamma_mom[1][3] + gamma_mom[2][3] + gamma_mom[3][3];

      lorentz_transf(boost_vec, vec_init, vec_end);

      value[i] = vec_end[3] - (p[23] / 2.);
    }
    else
    {
      value[i] = 999999.;
    }
  }

  if (abs(value[0]) < abs(value[1]))
    value_min = value[0];
  else
    value_min = value[1];

  return value_min;
}

Double_t minv_consv(Double_t *x, Double_t *p)
{
  Float_t gamma_mom[4][4], kaon_mom[4], neu_vtx[4], inv_mass_kaon;

  Double_t value[2], value_min, T0;

  T0 = p[27] * Trf;

  Reconstructor R;
  Solution S;

  // Setting clusters for a solution
  for (Int_t k = 0; k < 4; k++)
  {
    R.SetClu(k, p[k * 5],
             p[k * 5 + 1],
             p[k * 5 + 2],
             p[k * 5 + 3] + T0,
             p[k * 5 + 4]);

    R.SetClu(4, 0., 0., 0., 0., 0.);
    R.SetClu(5, 0., 0., 0., 0., 0.);
  }

  S = R.MySolve(selected);

  for (Int_t i = 0; i < 2; i++)
  {

    if (!S.error[i])
    {
      neu_vtx[0] = S.sol[i][0];
      neu_vtx[1] = S.sol[i][1];
      neu_vtx[2] = S.sol[i][2];
      neu_vtx[3] = S.sol[i][3];

      neutral_mom(p[0], p[1], p[2], p[4], neu_vtx, gamma_mom[0]);
      neutral_mom(p[5], p[6], p[7], p[9], neu_vtx, gamma_mom[1]);
      neutral_mom(p[10], p[11], p[12], p[14], neu_vtx, gamma_mom[2]);
      neutral_mom(p[15], p[16], p[17], p[19], neu_vtx, gamma_mom[3]);

      kaon_mom[0] = gamma_mom[0][0] + gamma_mom[1][0] + gamma_mom[2][0] + gamma_mom[3][0];
      kaon_mom[1] = gamma_mom[0][1] + gamma_mom[1][1] + gamma_mom[2][1] + gamma_mom[3][1];
      kaon_mom[2] = gamma_mom[0][2] + gamma_mom[1][2] + gamma_mom[2][2] + gamma_mom[3][2];
      kaon_mom[3] = gamma_mom[0][3] + gamma_mom[1][3] + gamma_mom[2][3] + gamma_mom[3][3];

      inv_mass_kaon = sqrt(pow(kaon_mom[3], 2) - pow(kaon_mom[0], 2) - pow(kaon_mom[1], 2) - pow(kaon_mom[2], 2));

      value[i] = inv_mass_kaon - PhysicsConstants::mK0;
    }
    else
    {
      value[i] = 999999.;
    }
  }

  if (abs(value[0]) < abs(value[1]))
    value_min = value[0];
  else
    value_min = value[1];

  return value_min;
}

Double_t x_consv(Double_t *x, Double_t *p)
{
  Float_t gamma_mom[4][4], kaon_mom[4], neu_vtx[4], bhabha_vtx[3] = {p[24], p[25], p[26]},
                                                    y_axis[3] = {0., p[21], 0.}, ip[3], dist[3], kaon_vel[3], tot_length, tot_vel;

  Double_t value[2], value_min, T0;

  T0 = p[27] * Trf;

  Reconstructor R;
  Solution S;

  // Setting clusters for a solution
  for (Int_t k = 0; k < 4; k++)
  {
    R.SetClu(k, p[k * 5],
             p[k * 5 + 1],
             p[k * 5 + 2],
             p[k * 5 + 3] + T0,
             p[k * 5 + 4]);

    R.SetClu(4, 0., 0., 0., 0., 0.);
    R.SetClu(5, 0., 0., 0., 0., 0.);
  }

  S = R.MySolve(selected);

  for (Int_t i = 0; i < 2; i++)
  {

    if (S.error[i] == 0)
    {
      neu_vtx[0] = S.sol[i][0];
      neu_vtx[1] = S.sol[i][1];
      neu_vtx[2] = S.sol[i][2];
      neu_vtx[3] = S.sol[i][3];

      neutral_mom(p[0], p[1], p[2], p[4], neu_vtx, gamma_mom[0]);
      neutral_mom(p[5], p[6], p[7], p[9], neu_vtx, gamma_mom[1]);
      neutral_mom(p[10], p[11], p[12], p[14], neu_vtx, gamma_mom[2]);
      neutral_mom(p[15], p[16], p[17], p[19], neu_vtx, gamma_mom[3]);

      kaon_mom[0] = gamma_mom[0][0] + gamma_mom[1][0] + gamma_mom[2][0] + gamma_mom[3][0];
      kaon_mom[1] = gamma_mom[0][1] + gamma_mom[1][1] + gamma_mom[2][1] + gamma_mom[3][1];
      kaon_mom[2] = gamma_mom[0][2] + gamma_mom[1][2] + gamma_mom[2][2] + gamma_mom[3][2];
      kaon_mom[3] = gamma_mom[0][3] + gamma_mom[1][3] + gamma_mom[2][3] + gamma_mom[3][3];

      plane_intersection(bhabha_vtx, y_axis, neu_vtx, kaon_mom, ip);

      ip[0] = p[24];
      ip[1] = p[25];
      if (abs(p[26] - ip[2]) > 2)
        ip[2] = p[26];

      dist[0] = neu_vtx[0] - ip[0];
      dist[1] = neu_vtx[1] - ip[1];
      dist[2] = neu_vtx[2] - ip[2];

      kaon_vel[0] = PhysicsConstants::cVel * kaon_mom[0] / kaon_mom[3];
      kaon_vel[1] = PhysicsConstants::cVel * kaon_mom[1] / kaon_mom[3];
      kaon_vel[2] = PhysicsConstants::cVel * kaon_mom[2] / kaon_mom[3];

      tot_length = sqrt(pow(dist[0], 2) + pow(dist[1], 2) + pow(dist[2], 2));
      tot_vel = sqrt(pow(kaon_vel[0], 2) + pow(kaon_vel[1], 2) + pow(kaon_vel[2], 2));

      value[i] = neu_vtx[3] - (dist[0] / kaon_vel[0]);
    }
    else
    {
      value[i] = 999999.;
    }
  }

  if (abs(value[0]) < abs(value[1]))
    value_min = value[0];
  else
    value_min = value[1];

  return value_min;
}

Double_t y_consv(Double_t *x, Double_t *p)
{
  Float_t gamma_mom[4][4], kaon_mom[4], neu_vtx[4], bhabha_vtx[3] = {p[24], p[25], p[26]},
                                                    y_axis[3] = {0., p[21], 0.}, ip[3], dist[3], kaon_vel[3];

  Double_t value[2], value_min, T0;

  T0 = p[27] * Trf;

  Reconstructor R;
  Solution S;

  // Setting clusters for a solution
  for (Int_t k = 0; k < 4; k++)
  {
    R.SetClu(k, p[k * 5],
             p[k * 5 + 1],
             p[k * 5 + 2],
             p[k * 5 + 3] + T0,
             p[k * 5 + 4]);

    R.SetClu(4, 0., 0., 0., 0., 0.);
    R.SetClu(5, 0., 0., 0., 0., 0.);
  }

  S = R.MySolve(selected);

  for (Int_t i = 0; i < 2; i++)
  {

    if (S.error[i] == 0)
    {
      neu_vtx[0] = S.sol[i][0];
      neu_vtx[1] = S.sol[i][1];
      neu_vtx[2] = S.sol[i][2];
      neu_vtx[3] = S.sol[i][3];

      neutral_mom(p[0], p[1], p[2], p[4], neu_vtx, gamma_mom[0]);
      neutral_mom(p[5], p[6], p[7], p[9], neu_vtx, gamma_mom[1]);
      neutral_mom(p[10], p[11], p[12], p[14], neu_vtx, gamma_mom[2]);
      neutral_mom(p[15], p[16], p[17], p[19], neu_vtx, gamma_mom[3]);

      kaon_mom[0] = gamma_mom[0][0] + gamma_mom[1][0] + gamma_mom[2][0] + gamma_mom[3][0];
      kaon_mom[1] = gamma_mom[0][1] + gamma_mom[1][1] + gamma_mom[2][1] + gamma_mom[3][1];
      kaon_mom[2] = gamma_mom[0][2] + gamma_mom[1][2] + gamma_mom[2][2] + gamma_mom[3][2];
      kaon_mom[3] = gamma_mom[0][3] + gamma_mom[1][3] + gamma_mom[2][3] + gamma_mom[3][3];

      plane_intersection(bhabha_vtx, y_axis, neu_vtx, kaon_mom, ip);

      ip[0] = p[24];
      ip[1] = p[25];
      if (abs(p[26] - ip[2]) > 2)
        ip[2] = p[26];

      dist[0] = neu_vtx[0] - ip[0];
      dist[1] = neu_vtx[1] - ip[1];
      dist[2] = neu_vtx[2] - ip[2];

      kaon_vel[0] = PhysicsConstants::cVel * kaon_mom[0] / kaon_mom[3];
      kaon_vel[1] = PhysicsConstants::cVel * kaon_mom[1] / kaon_mom[3];
      kaon_vel[2] = PhysicsConstants::cVel * kaon_mom[2] / kaon_mom[3];

      value[i] = neu_vtx[3] - (dist[1] / kaon_vel[1]);
    }
    else
    {
      value[i] = 999999.;
    }
  }

  if (abs(value[0]) < abs(value[1]))
    value_min = value[0];
  else
    value_min = value[1];

  return value_min;
}

Double_t z_consv(Double_t *x, Double_t *p)
{
  Float_t gamma_mom[4][4], kaon_mom[4], neu_vtx[4], bhabha_vtx[3] = {p[24], p[25], p[26]},
                                                    y_axis[3] = {0., p[21], 0.}, ip[3], dist[3], kaon_vel[3];

  Double_t value[2], value_min, T0;

  T0 = p[27] * Trf;

  Reconstructor R;
  Solution S;

  // Setting clusters for a solution
  for (Int_t k = 0; k < 4; k++)
  {
    R.SetClu(k, p[k * 5],
             p[k * 5 + 1],
             p[k * 5 + 2],
             p[k * 5 + 3] + T0,
             p[k * 5 + 4]);

    R.SetClu(4, 0., 0., 0., 0., 0.);
    R.SetClu(5, 0., 0., 0., 0., 0.);
  }

  S = R.MySolve(selected);

  for (Int_t i = 0; i < 2; i++)
  {

    if (S.error[i] == 0)
    {
      neu_vtx[0] = S.sol[i][0];
      neu_vtx[1] = S.sol[i][1];
      neu_vtx[2] = S.sol[i][2];
      neu_vtx[3] = S.sol[i][3];

      neutral_mom(p[0], p[1], p[2], p[4], neu_vtx, gamma_mom[0]);
      neutral_mom(p[5], p[6], p[7], p[9], neu_vtx, gamma_mom[1]);
      neutral_mom(p[10], p[11], p[12], p[14], neu_vtx, gamma_mom[2]);
      neutral_mom(p[15], p[16], p[17], p[19], neu_vtx, gamma_mom[3]);

      kaon_mom[0] = gamma_mom[0][0] + gamma_mom[1][0] + gamma_mom[2][0] + gamma_mom[3][0];
      kaon_mom[1] = gamma_mom[0][1] + gamma_mom[1][1] + gamma_mom[2][1] + gamma_mom[3][1];
      kaon_mom[2] = gamma_mom[0][2] + gamma_mom[1][2] + gamma_mom[2][2] + gamma_mom[3][2];
      kaon_mom[3] = gamma_mom[0][3] + gamma_mom[1][3] + gamma_mom[2][3] + gamma_mom[3][3];

      plane_intersection(bhabha_vtx, y_axis, neu_vtx, kaon_mom, ip);

      ip[0] = p[24];
      ip[1] = p[25];
      if (abs(p[26] - ip[2]) > 2)
        ip[2] = p[26];

      dist[0] = neu_vtx[0] - ip[0];
      dist[1] = neu_vtx[1] - ip[1];
      dist[2] = neu_vtx[2] - ip[2];

      kaon_vel[0] = PhysicsConstants::cVel * kaon_mom[0] / kaon_mom[3];
      kaon_vel[1] = PhysicsConstants::cVel * kaon_mom[1] / kaon_mom[3];
      kaon_vel[2] = PhysicsConstants::cVel * kaon_mom[2] / kaon_mom[3];

      value[i] = neu_vtx[3] - (dist[2] / kaon_vel[2]);
    }
    else
    {
      value[i] = 999999.;
    }
  }

  if (abs(value[0]) < abs(value[1]))
    value_min = value[0];
  else
    value_min = value[1];

  return value_min;
}

Double_t gamma1_consv(Double_t *x, Double_t *p)
{
  Float_t gamma_mom[4][4], kaon_mom[4], inv_mass_kaon, neu_vtx[4],
      bhabha_vtx[3] = {p[24], p[25], p[26]}, y_axis[3] = {0., p[21], 0.}, ip[3];

  Double_t value[2], value_min, T0;

  T0 = p[27] * Trf;

  Reconstructor R;
  Solution S;

  // Setting clusters for a solution
  for (Int_t k = 0; k < 4; k++)
  {
    R.SetClu(k, p[k * 5],
             p[k * 5 + 1],
             p[k * 5 + 2],
             p[k * 5 + 3] + T0,
             p[k * 5 + 4]);

    R.SetClu(4, 0., 0., 0., 0., 0.);
    R.SetClu(5, 0., 0., 0., 0., 0.);
  }

  S = R.MySolve(selected);

  for (Int_t i = 0; i < 2; i++)
  {

    if (S.error[i] == 0)
    {
      neu_vtx[0] = S.sol[i][0];
      neu_vtx[1] = S.sol[i][1];
      neu_vtx[2] = S.sol[i][2];
      neu_vtx[3] = S.sol[i][3];

      neutral_mom(p[0], p[1], p[2], p[4], neu_vtx, gamma_mom[0]);
      neutral_mom(p[5], p[6], p[7], p[9], neu_vtx, gamma_mom[1]);
      neutral_mom(p[10], p[11], p[12], p[14], neu_vtx, gamma_mom[2]);
      neutral_mom(p[15], p[16], p[17], p[19], neu_vtx, gamma_mom[3]);

      kaon_mom[0] = gamma_mom[0][0] + gamma_mom[1][0] + gamma_mom[2][0] + gamma_mom[3][0];
      kaon_mom[1] = gamma_mom[0][1] + gamma_mom[1][1] + gamma_mom[2][1] + gamma_mom[3][1];
      kaon_mom[2] = gamma_mom[0][2] + gamma_mom[1][2] + gamma_mom[2][2] + gamma_mom[3][2];
      kaon_mom[3] = gamma_mom[0][3] + gamma_mom[1][3] + gamma_mom[2][3] + gamma_mom[3][3];

      plane_intersection(bhabha_vtx, y_axis, neu_vtx, kaon_mom, ip);

      ip[0] = p[24];
      ip[1] = p[25];
      if (abs(p[26] - ip[2]) > 2)
        ip[2] = p[26];

      Float_t gamma_path = sqrt(pow(p[0] - neu_vtx[0], 2) + pow(p[1] - neu_vtx[1], 2) + pow(p[2] - neu_vtx[2], 2));

      value[i] = p[3] + T0 - (gamma_path / PhysicsConstants::cVel) - neu_vtx[3];
    }
    else
    {
      value[i] = 999999.;
    }
  }

  if (abs(value[0]) < abs(value[1]))
    value_min = value[0];
  else
    value_min = value[1];

  return value_min;
}

Double_t gamma2_consv(Double_t *x, Double_t *p)
{
  Float_t gamma_mom[4][4], kaon_mom[4], inv_mass_kaon, neu_vtx[4],
      bhabha_vtx[3] = {p[24], p[25], p[26]}, y_axis[3] = {0., p[21], 0.}, ip[3];

  Double_t value[2], value_min, T0;

  T0 = p[27] * Trf;

  Reconstructor R;
  Solution S;

  // Setting clusters for a solution
  for (Int_t k = 0; k < 4; k++)
  {
    R.SetClu(k, p[k * 5],
             p[k * 5 + 1],
             p[k * 5 + 2],
             p[k * 5 + 3] + T0,
             p[k * 5 + 4]);

    R.SetClu(4, 0., 0., 0., 0., 0.);
    R.SetClu(5, 0., 0., 0., 0., 0.);
  }

  S = R.MySolve(selected);

  for (Int_t i = 0; i < 2; i++)
  {

    if (S.error[i] == 0)
    {
      neu_vtx[0] = S.sol[i][0];
      neu_vtx[1] = S.sol[i][1];
      neu_vtx[2] = S.sol[i][2];
      neu_vtx[3] = S.sol[i][3];

      neutral_mom(p[0], p[1], p[2], p[4], neu_vtx, gamma_mom[0]);
      neutral_mom(p[5], p[6], p[7], p[9], neu_vtx, gamma_mom[1]);
      neutral_mom(p[10], p[11], p[12], p[14], neu_vtx, gamma_mom[2]);
      neutral_mom(p[15], p[16], p[17], p[19], neu_vtx, gamma_mom[3]);

      kaon_mom[0] = gamma_mom[0][0] + gamma_mom[1][0] + gamma_mom[2][0] + gamma_mom[3][0];
      kaon_mom[1] = gamma_mom[0][1] + gamma_mom[1][1] + gamma_mom[2][1] + gamma_mom[3][1];
      kaon_mom[2] = gamma_mom[0][2] + gamma_mom[1][2] + gamma_mom[2][2] + gamma_mom[3][2];
      kaon_mom[3] = gamma_mom[0][3] + gamma_mom[1][3] + gamma_mom[2][3] + gamma_mom[3][3];

      plane_intersection(bhabha_vtx, y_axis, neu_vtx, kaon_mom, ip);

      ip[0] = p[24];
      ip[1] = p[25];
      if (abs(p[26] - ip[2]) > 2)
        ip[2] = p[26];

      Float_t gamma_path = sqrt(pow(p[5] - neu_vtx[0], 2) + pow(p[6] - neu_vtx[1], 2) + pow(p[7] - neu_vtx[2], 2));

      value[i] = p[8] + T0 - (gamma_path / PhysicsConstants::cVel) - neu_vtx[3];
    }
    else
    {
      value[i] = 999999.;
    }
  }

  if (abs(value[0]) < abs(value[1]))
    value_min = value[0];
  else
    value_min = value[1];

  return value_min;
}

Double_t gamma3_consv(Double_t *x, Double_t *p)
{
  Float_t gamma_mom[4][4], kaon_mom[4], inv_mass_kaon, neu_vtx[4],
      bhabha_vtx[3] = {p[24], p[25], p[26]}, y_axis[3] = {0., p[21], 0.}, ip[3];

  Double_t value[2], value_min, T0;

  T0 = p[27] * Trf;

  Reconstructor R;
  Solution S;

  // Setting clusters for a solution
  for (Int_t k = 0; k < 4; k++)
  {
    R.SetClu(k, p[k * 5],
             p[k * 5 + 1],
             p[k * 5 + 2],
             p[k * 5 + 3] + T0,
             p[k * 5 + 4]);

    R.SetClu(4, 0., 0., 0., 0., 0.);
    R.SetClu(5, 0., 0., 0., 0., 0.);
  }

  S = R.MySolve(selected);

  for (Int_t i = 0; i < 2; i++)
  {

    if (S.error[i] == 0)
    {
      neu_vtx[0] = S.sol[i][0];
      neu_vtx[1] = S.sol[i][1];
      neu_vtx[2] = S.sol[i][2];
      neu_vtx[3] = S.sol[i][3];

      neutral_mom(p[0], p[1], p[2], p[4], neu_vtx, gamma_mom[0]);
      neutral_mom(p[5], p[6], p[7], p[9], neu_vtx, gamma_mom[1]);
      neutral_mom(p[10], p[11], p[12], p[14], neu_vtx, gamma_mom[2]);
      neutral_mom(p[15], p[16], p[17], p[19], neu_vtx, gamma_mom[3]);

      kaon_mom[0] = gamma_mom[0][0] + gamma_mom[1][0] + gamma_mom[2][0] + gamma_mom[3][0];
      kaon_mom[1] = gamma_mom[0][1] + gamma_mom[1][1] + gamma_mom[2][1] + gamma_mom[3][1];
      kaon_mom[2] = gamma_mom[0][2] + gamma_mom[1][2] + gamma_mom[2][2] + gamma_mom[3][2];
      kaon_mom[3] = gamma_mom[0][3] + gamma_mom[1][3] + gamma_mom[2][3] + gamma_mom[3][3];

      plane_intersection(bhabha_vtx, y_axis, neu_vtx, kaon_mom, ip);

      ip[0] = p[24];
      ip[1] = p[25];
      if (abs(p[26] - ip[2]) > 2)
        ip[2] = p[26];

      Float_t gamma_path = sqrt(pow(p[10] - neu_vtx[0], 2) + pow(p[11] - neu_vtx[1], 2) + pow(p[12] - neu_vtx[2], 2));

      value[i] = p[13] + T0 - (gamma_path / PhysicsConstants::cVel) - neu_vtx[3];
    }
    else
    {
      value[i] = 999999.;
    }
  }

  if (abs(value[0]) < abs(value[1]))
    value_min = value[0];
  else
    value_min = value[1];

  return value_min;
}

Double_t gamma4_consv(Double_t *x, Double_t *p)
{
  Float_t gamma_mom[4][4], kaon_mom[4], inv_mass_kaon, neu_vtx[4],
      bhabha_vtx[3] = {p[24], p[25], p[26]}, y_axis[3] = {0., p[21], 0.}, ip[3];

  Double_t value[2], value_min, T0;

  T0 = p[27] * Trf;

  Reconstructor R;
  Solution S;

  // Setting clusters for a solution
  for (Int_t k = 0; k < 4; k++)
  {
    R.SetClu(k, p[k * 5],
             p[k * 5 + 1],
             p[k * 5 + 2],
             p[k * 5 + 3] + T0,
             p[k * 5 + 4]);

    R.SetClu(4, 0., 0., 0., 0., 0.);
    R.SetClu(5, 0., 0., 0., 0., 0.);
  }

  S = R.MySolve(selected);

  for (Int_t i = 0; i < 2; i++)
  {

    if (S.error[i] == 0)
    {
      neu_vtx[0] = S.sol[i][0];
      neu_vtx[1] = S.sol[i][1];
      neu_vtx[2] = S.sol[i][2];
      neu_vtx[3] = S.sol[i][3];

      neutral_mom(p[0], p[1], p[2], p[4], neu_vtx, gamma_mom[0]);
      neutral_mom(p[5], p[6], p[7], p[9], neu_vtx, gamma_mom[1]);
      neutral_mom(p[10], p[11], p[12], p[14], neu_vtx, gamma_mom[2]);
      neutral_mom(p[15], p[16], p[17], p[19], neu_vtx, gamma_mom[3]);

      kaon_mom[0] = gamma_mom[0][0] + gamma_mom[1][0] + gamma_mom[2][0] + gamma_mom[3][0];
      kaon_mom[1] = gamma_mom[0][1] + gamma_mom[1][1] + gamma_mom[2][1] + gamma_mom[3][1];
      kaon_mom[2] = gamma_mom[0][2] + gamma_mom[1][2] + gamma_mom[2][2] + gamma_mom[3][2];
      kaon_mom[3] = gamma_mom[0][3] + gamma_mom[1][3] + gamma_mom[2][3] + gamma_mom[3][3];

      plane_intersection(bhabha_vtx, y_axis, neu_vtx, kaon_mom, ip);

      ip[0] = p[24];
      ip[1] = p[25];
      if (abs(p[26] - ip[2]) > 2)
        ip[2] = p[26];

      Float_t gamma_path = sqrt(pow(p[15] - neu_vtx[0], 2) + pow(p[16] - neu_vtx[1], 2) + pow(p[17] - neu_vtx[2], 2));

      value[i] = p[18] + T0 - (gamma_path / PhysicsConstants::cVel) - neu_vtx[3];
    }
    else
    {
      value[i] = 999999.;
    }
  }

  if (abs(value[0]) < abs(value[1]))
    value_min = value[0];
  else
    value_min = value[1];

  return value_min;
}