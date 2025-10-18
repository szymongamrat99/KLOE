// Author: Szymon Gamrat
// Last update: 02.06.2024

#include <neu_triangle.h>

using namespace std;

int neu_triangle(Float_t *TrcSumFinal, Float_t *vtxSigmaFinal, Float_t Clu5Vec[4][5], Float_t *ip, Float_t *Phi4Mom, Float_t *Kne4Mom, Float_t *Kne4Vec, Float_t *trc, ErrorHandling::ErrorLogs &logger)
{
  Float_t KneTotMom = 0., BetaK, CosPkD[4], DTot[4];
  Float_t A, B[4], C[4], Delta[4], lK[4][2], lGamma[4][2], lGammaFinal[4], trctmp[4][2], lKtrue[4];

  Float_t NeuVtxTmp[4][2][4], NeuVtxTrueClu[4][4], NeuVtxAvg[4] = {0.}, EneTot = 0., TrcSum = 0., vtxSigma = 0.;

  TVector3 D[4], pK(Kne4Mom[0], Kne4Mom[1], Kne4Mom[2]);

  TrcSum = 0.;
  vtxSigma = 0.;

  KneTotMom = pK.Mag();
  BetaK = KneTotMom / Kne4Mom[3];

  A = (1 - pow(BetaK, 2)) / pow(BetaK, 2);

  for (Int_t i = 0; i < 4; i++)
  {
    for (Int_t j = 0; j < 3; j++)
    {
      D[i](j) = Clu5Vec[i][j] - ip[j];
    }

    DTot[i] = D[i].Mag();

    CosPkD[i] = (D[i].Dot(pK)) / (DTot[i] * KneTotMom);

    B[i] = 2 * (DTot[i] * CosPkD[i] - (PhysicsConstants::cVel * Clu5Vec[i][3] / BetaK));
    C[i] = pow(PhysicsConstants::cVel * Clu5Vec[i][3], 2) - pow(DTot[i], 2);

    Delta[i] = pow(B[i], 2) - 4 * A * C[i];

    try
    {
      if (Delta[i] < 0.)
      {
        throw ErrorHandling::ErrorCodes::DELTA_LT_ZERO;
      }
      else if (A == 0.)
      {
        throw ErrorHandling::ErrorCodes::DENOM_EQ_ZERO;
      }
      else
      {
        lK[i][0] = (-B[i] - sqrt(Delta[i])) / (2. * A);
        lK[i][1] = (-B[i] + sqrt(Delta[i])) / (2. * A);

        for (Int_t j = 0; j < 2; j++)
        {
          for (Int_t k = 0; k < 3; k++)
            NeuVtxTmp[i][j][k] = ip[k] + (lK[i][j] * (pK(k) / KneTotMom));

          NeuVtxTmp[i][j][3] = (lK[i][j] / (PhysicsConstants::cVel * BetaK));

          lGamma[i][j] = sqrt(pow(Clu5Vec[i][0] - NeuVtxTmp[i][j][0], 2) +
                              pow(Clu5Vec[i][1] - NeuVtxTmp[i][j][1], 2) +
                              pow(Clu5Vec[i][2] - NeuVtxTmp[i][j][2], 2));

          trctmp[i][j] = Clu5Vec[i][3] - NeuVtxTmp[i][j][3] - (lGamma[i][j] / PhysicsConstants::cVel);
        }

        if (abs(trctmp[i][0]) < abs(trctmp[i][1]))
        {
          for (Int_t j = 0; j < 4; j++)
            NeuVtxTrueClu[i][j] = NeuVtxTmp[i][0][j];
        }
        else
        {
          for (Int_t j = 0; j < 4; j++)
            NeuVtxTrueClu[i][j] = NeuVtxTmp[i][1][j];
        }

        for (Int_t j = 0; j < 3; j++)
          NeuVtxAvg[j] += Clu5Vec[i][4] * NeuVtxTrueClu[i][j];

        EneTot += Clu5Vec[i][4];
      }
    }
    catch(ErrorHandling::ErrorCodes err)
    {
      logger.getErrLog(err);
      return int(err);
    }
  }

  for (Int_t j = 0; j < 3; j++)
    NeuVtxAvg[j] = NeuVtxAvg[j] / EneTot;

  NeuVtxAvg[3] = sqrt(pow(NeuVtxAvg[0] - ip[0], 2) +
                      pow(NeuVtxAvg[1] - ip[1], 2) +
                      pow(NeuVtxAvg[2] - ip[2], 2)) /
                 (BetaK * PhysicsConstants::cVel);

  for (Int_t i = 0; i < 4; i++)
  {
    lGammaFinal[i] = sqrt(pow(Clu5Vec[i][0] - NeuVtxAvg[0], 2) +
                          pow(Clu5Vec[i][1] - NeuVtxAvg[1], 2) +
                          pow(Clu5Vec[i][2] - NeuVtxAvg[2], 2));
    trc[i] = Clu5Vec[i][3] - NeuVtxAvg[3] - (lGammaFinal[i] / PhysicsConstants::cVel);

    vtxSigma += Clu5Vec[i][4] * sqrt(pow(NeuVtxTrueClu[i][0] - NeuVtxAvg[0], 2) + pow(NeuVtxTrueClu[i][1] - NeuVtxAvg[1], 2) + pow(NeuVtxTrueClu[i][2] - NeuVtxAvg[2], 2)) / EneTot;
  }

  Kne4Vec[0] = NeuVtxAvg[0];
  Kne4Vec[1] = NeuVtxAvg[1];
  Kne4Vec[2] = NeuVtxAvg[2];
  Kne4Vec[3] = NeuVtxAvg[3];

  TrcSum = trc[0] + trc[1] + trc[2] + trc[3];

  *TrcSumFinal = TrcSum;
  *vtxSigmaFinal = vtxSigma;

  return 0;
};
