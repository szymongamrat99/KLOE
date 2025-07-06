#include <ConstraintsTrilateration.h>

namespace KLOE
{
  Double_t ConstraintsTrilateration::EnergyConsvCM(Double_t *x, Double_t *p)
  {
    Double_t boost_vec[3] = {-p[20] / p[23], -p[21] / p[23], -p[22] / p[23]},
             gamma_mom[4][4],
             neu_vtx[4],
             value[2],
             value_min;

    Reconstructor R; // Reconstructor object
    Solution S;      // Solution struct

    // Setting clusters for a solution
    for (Int_t k = 0; k < 4; k++)
    {
      R.SetClu(k, p[k * 5],
               p[k * 5 + 1],
               p[k * 5 + 2],
               p[k * 5 + 3],
               p[k * 5 + 4]);

      R.SetClu(4, 0., 0., 0., 0., 0.);
      R.SetClu(5, 0., 0., 0., 0., 0.);
    }
    // -------------------------------

    S = R.MySolve(_selected); // Filling up the structure

    for (Int_t i = 0; i < 2; i++)
    {
      if (S.error[i] == 0)
      {
        neu_vtx[0] = S.sol[i][0];
        neu_vtx[1] = S.sol[i][1];
        neu_vtx[2] = S.sol[i][2];
        neu_vtx[3] = S.sol[i][3];

        for (Int_t j = 0; j < 4; j++)
          neutral_mom(p[j * 5],
                      p[j * 5 + 1],
                      p[j * 5 + 2],
                      p[j * 5 + 4],
                      neu_vtx,
                      gamma_mom[j]);

        Double_t
            vec_init[4] = {0.},
            vec_end[4] = {0.};

        for (Int_t j = 0; j < 4; j++)
          for (Int_t k = 0; k < 4; k++)
            vec_init[j] += gamma_mom[k][j];

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

  Double_t ConstraintsTrilateration::MinvConsv(Double_t *x, Double_t *p)
  {
    Double_t gamma_mom[4][4],
        neu_vtx[4],
        inv_mass_kaon,
        value[2],
        value_min;

    Reconstructor R; // Reconstructor object
    Solution S;      // Solution struct

    // Setting clusters for a solution
    for (Int_t k = 0; k < 4; k++)
    {
      R.SetClu(k, p[k * 5],
               p[k * 5 + 1],
               p[k * 5 + 2],
               p[k * 5 + 3],
               p[k * 5 + 4]);

      R.SetClu(4, 0., 0., 0., 0., 0.);
      R.SetClu(5, 0., 0., 0., 0., 0.);
    }
    // -------------------------------

    S = R.MySolve(_selected); // Filling up the structure

    for (Int_t i = 0; i < 2; i++)
    {
      if (!S.error[i])
      {
        neu_vtx[0] = S.sol[i][0];
        neu_vtx[1] = S.sol[i][1];
        neu_vtx[2] = S.sol[i][2];
        neu_vtx[3] = S.sol[i][3];

        for (Int_t j = 0; j < 4; j++)
          neutral_mom(p[j * 5],
                      p[j * 5 + 1],
                      p[j * 5 + 2],
                      p[j * 5 + 4],
                      neu_vtx,
                      gamma_mom[j]);

        Double_t
            kaon_mom[4] = {0.};

        for (Int_t j = 0; j < 4; j++)
          for (Int_t k = 0; k < 4; k++)
            kaon_mom[j] += gamma_mom[k][j];

        inv_mass_kaon = sqrt(pow(kaon_mom[3], 2) -
                             pow(kaon_mom[0], 2) -
                             pow(kaon_mom[1], 2) -
                             pow(kaon_mom[2], 2));

        value[i] = inv_mass_kaon - constants.getKaonMass();
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

  Double_t ConstraintsTrilateration::NeutralPathConsvLAB(Double_t *x, Double_t *p)
  {
    Double_t
        gamma_mom[4][4],
        neu_vtx[4],
        bhabha_vtx[3] = {p[24], p[25], p[26]},
        y_axis[3] = {0., p[21], 0.},
        ip[3],
        dist,
        kaon_vel,
        tot_length,
        tot_vel,
        value[2],
        value_min;

    Reconstructor R;
    Solution S;

    // Setting clusters for a solution
    for (Int_t k = 0; k < 4; k++)
    {
      R.SetClu(k, p[k * 5],
               p[k * 5 + 1],
               p[k * 5 + 2],
               p[k * 5 + 3],
               p[k * 5 + 4]);

      R.SetClu(4, 0., 0., 0., 0., 0.);
      R.SetClu(5, 0., 0., 0., 0., 0.);
    }

    S = R.MySolve(_selected);

    for (Int_t i = 0; i < 2; i++)
    {
      if (S.error[i] == 0)
      {
        neu_vtx[0] = S.sol[i][0];
        neu_vtx[1] = S.sol[i][1];
        neu_vtx[2] = S.sol[i][2];
        neu_vtx[3] = S.sol[i][3];

        for (Int_t j = 0; j < 4; j++)
          neutral_mom(p[j * 5],
                      p[j * 5 + 1],
                      p[j * 5 + 2],
                      p[j * 5 + 4],
                      neu_vtx,
                      gamma_mom[j]);

        Double_t
            kaon_mom[4] = {0.};

        for (Int_t j = 0; j < 4; j++)
          for (Int_t k = 0; k < 4; k++)
            kaon_mom[j] += gamma_mom[k][j];

        plane_intersection(bhabha_vtx, y_axis, neu_vtx, kaon_mom, ip);

        ip[0] = p[24];
        ip[1] = p[25];
        if (abs(p[26] - ip[2]) > 2)
          ip[2] = p[26];

        dist = neu_vtx[_chosenComponent] - ip[_chosenComponent];
        kaon_vel = constants.cVel * kaon_mom[_chosenComponent] / kaon_mom[3];

        value[i] = neu_vtx[3] - (dist / kaon_vel);
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

}
