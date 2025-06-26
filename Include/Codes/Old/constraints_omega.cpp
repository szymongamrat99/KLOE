#include <constraints_omega.h>

namespace OmegaConstraints
{
  Double_t ene_consv(Double_t *x, Double_t *p)
  {
    Float_t
        cluster[4][5],
        gamma_mom[4][4],
        neu_vtx[3],
        bhabha_mom[4],
        Kchrec[4];

    Double_t
        value_min;

    for (Int_t i = 0; i < 3; i++)
      neu_vtx[i] = p[20 + i];
      
    for (Int_t i = 0; i < 4; i++)
    {
      bhabha_mom[i] = p[27 + i];
      Kchrec[i] = p[23 + i];

      for (Int_t j = 0; j < 5; j++)
        cluster[i][j] = p[i * 5 + j];

      neutral_mom(cluster[i][0], cluster[i][1], cluster[i][2], cluster[i][4], neu_vtx, gamma_mom[i]);
    }

    value_min = (bhabha_mom[3] - Kchrec[3] - gamma_mom[0][3] - gamma_mom[1][3] - gamma_mom[2][3] - gamma_mom[3][3]);

    return value_min;
  }

  Double_t px_consv(Double_t *x, Double_t *p)
  {
    Float_t
        cluster[4][5],
        gamma_mom[4][4],
        neu_vtx[3],
        bhabha_mom[4],
        Kchrec[4];

    Double_t
        value_min;

    for (Int_t i = 0; i < 3; i++)
      neu_vtx[i] = p[20 + i];

    for (Int_t i = 0; i < 4; i++)
    {
      bhabha_mom[i] = p[27 + i];
      Kchrec[i] = p[23 + i];

      for (Int_t j = 0; j < 5; j++)
        cluster[i][j] = p[i * 5 + j];

      neutral_mom(cluster[i][0], cluster[i][1], cluster[i][2], cluster[i][4], neu_vtx, gamma_mom[i]);
    }

    value_min = (bhabha_mom[0] - Kchrec[0] - gamma_mom[0][0] - gamma_mom[1][0] - gamma_mom[2][0] - gamma_mom[3][0]);

    return value_min;
  }

  Double_t py_consv(Double_t *x, Double_t *p)
  {
    Float_t
        cluster[4][5],
        gamma_mom[4][4],
        neu_vtx[3],
        bhabha_mom[4],
        Kchrec[4];

    Double_t
        value_min;

    for(Int_t i = 0; i < 3; i++)
      neu_vtx[i] = p[20 + i];

    for(Int_t i = 0; i < 4; i++)
    {
      bhabha_mom[i] = p[27 + i];
      Kchrec[i] = p[23 + i];

      for(Int_t j = 0; j < 5; j++)
        cluster[i][j] = p[i * 5 + j];

      neutral_mom(cluster[i][0], cluster[i][1], cluster[i][2], cluster[i][4], neu_vtx, gamma_mom[i]);
    }

    value_min = (bhabha_mom[1] - Kchrec[1] - gamma_mom[0][1] - gamma_mom[1][1] - gamma_mom[2][1] - gamma_mom[3][1]);

    return value_min;
  }

  Double_t pz_consv(Double_t *x, Double_t *p)
  {
    Float_t
        cluster[4][5],
        gamma_mom[4][4],
        neu_vtx[3],
        bhabha_mom[4],
        Kchrec[4];

    Double_t
        value_min;

    for(Int_t i = 0; i < 3; i++)
      neu_vtx[i] = p[20 + i];

    for(Int_t i = 0; i < 4; i++)
    {
      bhabha_mom[i] = p[27 + i];
      Kchrec[i] = p[23 + i];

      for(Int_t j = 0; j < 5; j++)
        cluster[i][j] = p[i * 5 + j];

      neutral_mom(cluster[i][0], cluster[i][1], cluster[i][2], cluster[i][4], neu_vtx, gamma_mom[i]);
    }

    value_min = (bhabha_mom[2] - Kchrec[2] - gamma_mom[0][2] - gamma_mom[1][2] - gamma_mom[2][2] - gamma_mom[3][2]);

    return value_min;
  }

  Double_t photon1_consv(Double_t *x, Double_t *p)
  {
    Float_t
        cluster[4][5],
        neu_vtx[3],
        R_gamma = 0.;

    Double_t
        value_min;

    for(Int_t i = 0; i < 3; i++)
      neu_vtx[i] = p[20 + i];

    for(Int_t i = 0; i < 4; i++)
    {
      for(Int_t j = 0; j < 5; j++)
        cluster[i][j] = p[i * 5 + j];
    }

    R_gamma = sqrt(pow(cluster[0][0] - neu_vtx[0], 2) + 
                   pow(cluster[0][1] - neu_vtx[1], 2) +
                   pow(cluster[0][2] - neu_vtx[2], 2) );
    

    value_min = cVel * cluster[0][3] - R_gamma;

    return value_min;
  }

  Double_t photon2_consv(Double_t *x, Double_t *p)
  {
    Float_t
        cluster[4][5],
        neu_vtx[3],
        R_gamma = 0.;

    Double_t
        value_min;

    for(Int_t i = 0; i < 3; i++)
      neu_vtx[i] = p[20 + i];

    for(Int_t i = 0; i < 4; i++)
    {
      for(Int_t j = 0; j < 5; j++)
        cluster[i][j] = p[i * 5 + j];
    }

    R_gamma = sqrt(pow(cluster[1][0] - neu_vtx[0], 2) + 
                   pow(cluster[1][1] - neu_vtx[1], 2) +
                   pow(cluster[1][2] - neu_vtx[2], 2) );
    

    value_min = cVel * cluster[1][3] - R_gamma;

    return value_min;
  }

  Double_t photon3_consv(Double_t *x, Double_t *p)
  {
    Float_t
        cluster[4][5],
        neu_vtx[3],
        R_gamma = 0.;

    Double_t
        value_min;

    for(Int_t i = 0; i < 3; i++)
      neu_vtx[i] = p[20 + i];

    for(Int_t i = 0; i < 4; i++)
    {
      for(Int_t j = 0; j < 5; j++)
        cluster[i][j] = p[i * 5 + j];
    }

    R_gamma = sqrt(pow(cluster[2][0] - neu_vtx[0], 2) + 
                   pow(cluster[2][1] - neu_vtx[1], 2) +
                   pow(cluster[2][2] - neu_vtx[2], 2) );
    

    value_min = cVel * cluster[2][3] - R_gamma;

    return value_min;
  };
  
  Double_t photon4_consv(Double_t *x, Double_t *p)
  {
        Float_t
        cluster[4][5],
        neu_vtx[3],
        R_gamma = 0.;

    Double_t
        value_min;

    for(Int_t i = 0; i < 3; i++)
      neu_vtx[i] = p[20 + i];

    for(Int_t i = 0; i < 4; i++)
    {
      for(Int_t j = 0; j < 5; j++)
        cluster[i][j] = p[i * 5 + j];
    }

    R_gamma = sqrt(pow(cluster[3][0] - neu_vtx[0], 2) + 
                   pow(cluster[3][1] - neu_vtx[1], 2) +
                   pow(cluster[3][2] - neu_vtx[2], 2) );
    

    value_min = cVel * cluster[3][3] - R_gamma;

    return value_min;
  }

}