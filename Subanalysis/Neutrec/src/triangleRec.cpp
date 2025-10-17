
#include <ErrorLogs.h>
#include <KinFitter.h>
#include <uncertainties.h>
#include <reconstructor.h>
#include <boost/optional.hpp>

#include "../inc/trilateration.hpp"

ErrorHandling::ErrorCodes TriangleRec(std::vector<Int_t> g4taken_kinfit, std::vector<Float_t> cluster[5], std::vector<Int_t> Asscl, std::vector<Float_t> bhabha_mom, std::vector<Float_t> Kchboost, std::vector<Float_t> ip, std::vector<Float_t> &Knetriangle, std::vector<Float_t> gammatriangle[4], Float_t &minv4gam, std::vector<Float_t> &trcfinal, ErrorHandling::ErrorLogs &logger)
{
  Int_t ind_gam[4];

  Bool_t cond_ene;
  Bool_t cond_clus[4];
  Float_t Clu5Vec[4][5];

  ind_gam[0] = g4taken_kinfit[0];
  ind_gam[1] = g4taken_kinfit[1];
  ind_gam[2] = g4taken_kinfit[2];
  ind_gam[3] = g4taken_kinfit[3];

  cond_ene = cluster[4][Asscl[ind_gam[0]] - 1] > 0 && cluster[4][Asscl[ind_gam[1]] - 1] > 0 &&
             cluster[4][Asscl[ind_gam[2]] - 1] > 0 && cluster[4][Asscl[ind_gam[3]] - 1] > 0;

  cond_clus[0] = cluster[0][Asscl[ind_gam[0]] - 1] != 0 && cluster[1][Asscl[ind_gam[0]] - 1] != 0 && cluster[2][Asscl[ind_gam[0]] - 1] != 0;
  cond_clus[1] = cluster[0][Asscl[ind_gam[1]] - 1] != 0 && cluster[1][Asscl[ind_gam[1]] - 1] != 0 && cluster[2][Asscl[ind_gam[1]] - 1] != 0;
  cond_clus[2] = cluster[0][Asscl[ind_gam[2]] - 1] != 0 && cluster[1][Asscl[ind_gam[2]] - 1] != 0 && cluster[2][Asscl[ind_gam[2]] - 1] != 0;
  cond_clus[3] = cluster[0][Asscl[ind_gam[3]] - 1] != 0 && cluster[1][Asscl[ind_gam[3]] - 1] != 0 && cluster[2][Asscl[ind_gam[3]] - 1] != 0;
  
  if (cond_ene == true && cond_clus[0] && cond_clus[1] && cond_clus[2] && cond_clus[3])
  {
    for (Int_t k = 0; k < 4; k++)
    {
      Clu5Vec[k][0] = cluster[0][Asscl[ind_gam[k]] - 1];
      Clu5Vec[k][1] = cluster[1][Asscl[ind_gam[k]] - 1];
      Clu5Vec[k][2] = cluster[2][Asscl[ind_gam[k]] - 1];
      Clu5Vec[k][3] = cluster[3][Asscl[ind_gam[k]] - 1];
      Clu5Vec[k][4] = cluster[4][Asscl[ind_gam[k]] - 1];
    }

    //! Using the charged part of the decay

    Float_t Knerec[9] = {};

    Knerec[0] = bhabha_mom[0] - Kchboost[0];
    Knerec[1] = bhabha_mom[1] - Kchboost[1];
    Knerec[2] = bhabha_mom[2] - Kchboost[2];
    Knerec[3] = bhabha_mom[3] - Kchboost[3];

    //!

    Float_t TrcSum = 0., vtxSigma = 0., vtxSigmaMin = 1.e6, TrcSumMin = 1.e6, trcsum = 0.;

    Float_t neu_vtx[4], trc[4];

    Int_t done = 0;

    neu_triangle(&TrcSum, &vtxSigma, Clu5Vec, ip.data(), bhabha_mom.data(), Knerec, neu_vtx, trc, logger);

    if (sqrt(pow(vtxSigma, 2) + pow(TrcSum, 2)) < sqrt(pow(vtxSigmaMin, 2) + pow(TrcSumMin, 2)))
    {
      vtxSigmaMin = vtxSigma;
      TrcSumMin = TrcSum;

      trcsum = TrcSumMin;

      done = 1;

      for (Int_t l = 0; l < 4; l++)
      {
        neutral_mom(cluster[0][Asscl[ind_gam[l]] - 1], cluster[1][Asscl[ind_gam[l]] - 1], cluster[2][Asscl[ind_gam[l]] - 1], cluster[4][Asscl[ind_gam[l]] - 1], neu_vtx, gammatriangle[l].data());

        Knetriangle[l] = gammatriangle[0][l] + gammatriangle[1][l] + gammatriangle[2][l] + gammatriangle[3][l];

        Knetriangle[6 + l] = neu_vtx[l];

        trcfinal[l] = trc[l];
      }

      minv4gam = sqrt(pow(gammatriangle[0][3] + gammatriangle[1][3] + gammatriangle[2][3] + gammatriangle[3][3], 2) -
                      pow(gammatriangle[0][0] + gammatriangle[1][0] + gammatriangle[2][0] + gammatriangle[3][0], 2) -
                      pow(gammatriangle[0][1] + gammatriangle[1][1] + gammatriangle[2][1] + gammatriangle[3][1], 2) -
                      pow(gammatriangle[0][2] + gammatriangle[1][2] + gammatriangle[2][2] + gammatriangle[3][2], 2));

      Knetriangle[4] = 0.;
      for (Int_t l = 0; l < 3; l++)
      {
        Knetriangle[4] += pow(Knetriangle[l], 2);
      }

      Knetriangle[5] = sqrt(pow(Knetriangle[3], 2) - Knetriangle[4]);
      Knetriangle[4] = sqrt(Knetriangle[4]);
    }
  }

  return ErrorHandling::ErrorCodes::NO_ERROR;
}