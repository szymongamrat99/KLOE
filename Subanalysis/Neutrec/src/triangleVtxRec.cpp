
#include <ErrorLogs.h>
#include <KinFitter.h>
#include <uncertainties.h>
#include <reconstructor.h>
#include <boost/optional.hpp>

#include "../inc/trilateration.hpp"

ErrorHandling::ErrorCodes TriangleVtxRec(std::vector<Double_t> cluster[5], std::vector<Int_t> &neu_clu_list, std::vector<Double_t> &bhabha_mom, std::vector<Double_t> &ip, std::vector<Double_t> &Kchboost, std::vector<Int_t> &g4taken_triangle, std::vector<Double_t> photon_triangle[4], std::vector<Double_t> &Knetriangle, Double_t &minv4gamTriangle, Double_t &TrcFinalTriangle, ErrorHandling::ErrorLogs &logger)
{

  const Int_t
      photonNum = 4,
      neuCluSize = neu_clu_list.size();

  Double_t MIN_CLU_ENE = 20.0; // Minimum energy of a cluster to be considered in the reconstruction

  std::vector<std::array<Double_t, 4>>
      partialSolutions;

  std::vector<Double_t>
      partialSolutionEnergy;

  std::vector<std::array<Int_t, 4>>
      indices;

  // Implementation of the four-gamma vertex reconstruction algorithm
  if (neuCluSize < photonNum)
    return ErrorHandling::ErrorCodes::LESS_THAN_FOUR_NEUTRAL_CLUSTERS;

  // Vector of combinations of photonNum clusters
  std::vector<std::array<Int_t, photonNum>> combinations;
  for (Int_t i1 = 0; i1 < neuCluSize - 3; i1++)
    for (Int_t i2 = i1 + 1; i2 < neuCluSize - 2; i2++)
      for (Int_t i3 = i2 + 1; i3 < neuCluSize - 1; i3++)
        for (Int_t i4 = i3 + 1; i4 < neuCluSize; i4++)
        {
          combinations.push_back({neu_clu_list[i1],
                                  neu_clu_list[i2],
                                  neu_clu_list[i3],
                                  neu_clu_list[i4]});
          
          indices.push_back({i1, i2, i3, i4});
        }

  //

  std::array<Double_t, 9> Knetriangleboost;

  // Setting up the Knetriangleboost vector based on bhabha_mom and Kchboost
  Knetriangleboost[0] = bhabha_mom[0] - Kchboost[0];
  Knetriangleboost[1] = bhabha_mom[1] - Kchboost[1];
  Knetriangleboost[2] = bhabha_mom[2] - Kchboost[2];
  Knetriangleboost[3] = bhabha_mom[3] - Kchboost[3];
  Knetriangleboost[4] = std::sqrt(pow(Knetriangleboost[0], 2) + pow(Knetriangleboost[1], 2) + pow(Knetriangleboost[2], 2));
  Knetriangleboost[5] = PhysicsConstants::mK0;

  std::array<Double_t, 3> boostKneu = {Knetriangleboost[0] / Knetriangleboost[3], Knetriangleboost[1] / Knetriangleboost[3], Knetriangleboost[2] / Knetriangleboost[3]};
  Double_t boostKneuMag = Knetriangleboost[4] / Knetriangleboost[3];

  ////////////////////////////////////////////////////////////////////////////

  std::array<Double_t, 3> ipCluster[photonNum];
  std::array<Double_t, 3> neuvtxTmp[photonNum], neuvtx[photonNum];
  Double_t ipClusterMag[photonNum];
  Double_t cosTheta[photonNum];
  Double_t a, b[photonNum], c[photonNum], delta[photonNum], kpath[photonNum], kpathTmp[photonNum], trc[photonNum], trcTmp[photonNum];

  std::array<Double_t, 3> vcoor;

  Double_t chigam = 0.0, chigamMin = 1e6;

  a = 1.0 - (1.0 / pow(boostKneuMag, 2));

  bool warning = false, success = false;

  // ---

  for (auto &combo : combinations)
  {
    warning = false;
    chigam = 0.0;

    Double_t totalEnergy = cluster[4][combo[0] - 1] + cluster[4][combo[1] - 1] +
                           cluster[4][combo[2] - 1] + cluster[4][combo[3] - 1];

    Bool_t energyLimitPerCluster = cluster[4][combo[0] - 1] > MIN_CLU_ENE &&
                                   cluster[4][combo[1] - 1] > MIN_CLU_ENE &&
                                   cluster[4][combo[2] - 1] > MIN_CLU_ENE &&
                                   cluster[4][combo[3] - 1] > MIN_CLU_ENE,
           condTotal = energyLimitPerCluster;

    // Reject clusters, which do not meet energy conditions
    if (!condTotal)
      continue;

    //////////////////////////////////////////////////////////////

    for (Int_t i = 0; i < photonNum; i++)
    {
      ipCluster[i] = {cluster[0][combo[i] - 1] - ip[0], 
                      cluster[1][combo[i] - 1] - ip[1], 
                      cluster[2][combo[i] - 1] - ip[2]};

      ipClusterMag[i] = std::sqrt(std::pow(ipCluster[i][0], 2) + std::pow(ipCluster[i][1], 2) + std::pow(ipCluster[i][2], 2));

      cosTheta[i] = (ipCluster[i][0] * boostKneu[0] + ipCluster[i][1] * boostKneu[1] + ipCluster[i][2] * boostKneu[2]) / (ipClusterMag[i] * boostKneuMag);

      b[i] = 2.0 * ((PhysicsConstants::cVel * cluster[3][combo[i] - 1] / boostKneuMag) - ipClusterMag[i] * cosTheta[i]);
      c[i] = std::pow(ipClusterMag[i], 2) - std::pow(PhysicsConstants::cVel * cluster[3][combo[i] - 1], 2);

      delta[i] = std::pow(b[i], 2) - 4.0 * a * c[i];

      if (delta[i] < 0)
      {
        warning = true;
        break;
      }
      else if (delta[i] == 0)
      {
        kpath[i] = -b[i] / (2.0 * a);
        kpathTmp[i] = kpath[i];
      }
      else
      {
        kpath[i] = (-b[i] + std::sqrt(delta[i])) / (2.0 * a);
        kpathTmp[i] = (-b[i] - std::sqrt(delta[i])) / (2.0 * a);        
      }

      neuvtx[i][0] = ip[0] + (kpath[i] * boostKneu[0] / boostKneuMag);
      neuvtx[i][1] = ip[1] + (kpath[i] * boostKneu[1] / boostKneuMag);
      neuvtx[i][2] = ip[2] + (kpath[i] * boostKneu[2] / boostKneuMag);

      neuvtxTmp[i][0] = ip[0] + (kpathTmp[i] * boostKneu[0] / boostKneuMag);
      neuvtxTmp[i][1] = ip[1] + (kpathTmp[i] * boostKneu[1] / boostKneuMag);
      neuvtxTmp[i][2] = ip[2] + (kpathTmp[i] * boostKneu[2] / boostKneuMag);

      trc[i] = cluster[3][combo[i] - 1] - (std::sqrt(std::pow(cluster[0][combo[i] - 1] - neuvtx[i][0], 2) + std::pow(cluster[1][combo[i] - 1] - neuvtx[i][1], 2) + std::pow(cluster[2][combo[i] - 1] - neuvtx[i][2], 2)) / PhysicsConstants::cVel) - (kpath[i] / (boostKneuMag * PhysicsConstants::cVel));

      trcTmp[i] = cluster[3][combo[i] - 1] - (std::sqrt(std::pow(cluster[0][combo[i] - 1] - neuvtxTmp[i][0], 2) + std::pow(cluster[1][combo[i] - 1] - neuvtxTmp[i][1], 2) + std::pow(cluster[2][combo[i] - 1] - neuvtxTmp[i][2], 2)) / PhysicsConstants::cVel) - (kpathTmp[i] / (boostKneuMag * PhysicsConstants::cVel));

      if (abs(trc[i]) > abs(trcTmp[i]))
      {
        kpath[i] = kpathTmp[i];
        neuvtx[i][0] = neuvtxTmp[i][0];
        neuvtx[i][1] = neuvtxTmp[i][1];
        neuvtx[i][2] = neuvtxTmp[i][2];
        trc[i] = trcTmp[i];
      }
    }

    if (warning)
      continue;

    for (Int_t i = 0; i < 3; i++)
    {
      vcoor[i] = (neuvtx[0][i] * cluster[4][combo[0] - 1] + neuvtx[1][i] * cluster[4][combo[1] - 1] + neuvtx[2][i] * cluster[4][combo[2] - 1] + neuvtx[3][i] * cluster[4][combo[3] - 1]) / totalEnergy;
    }

    for (Int_t i = 0; i < photonNum; i++)
    {
      chigam += std::sqrt(std::pow(neuvtx[i][0] - vcoor[0], 2) + std::pow(neuvtx[i][1] - vcoor[1], 2) + std::pow(neuvtx[i][2] - vcoor[2], 2)) * cluster[4][combo[i] - 1];
    }

    chigam /= totalEnergy;

    if (chigam < chigamMin)
    {
      success = true;
      chigamMin = chigam;

      Knetriangle[0] = 0.0;
      Knetriangle[1] = 0.0;
      Knetriangle[2] = 0.0;
      Knetriangle[3] = 0.0;

      for (Int_t i = 0; i < photonNum; i++)
      {
        neutral_mom(cluster[0][combo[i] - 1], cluster[1][combo[i] - 1], cluster[2][combo[i] - 1], cluster[4][combo[i] - 1], vcoor.data(), photon_triangle[i].data());
        Knetriangle[0] += photon_triangle[i][0];
        Knetriangle[1] += photon_triangle[i][1];
        Knetriangle[2] += photon_triangle[i][2];
        Knetriangle[3] += photon_triangle[i][3];
      }

      Knetriangle[4] = std::sqrt(std::pow(Knetriangle[0], 2) + std::pow(Knetriangle[1], 2) + std::pow(Knetriangle[2], 2));
      Knetriangle[5] = std::sqrt(std::pow(Knetriangle[3], 2) - std::pow(Knetriangle[4], 2));

      Knetriangle[6] = vcoor[0];
      Knetriangle[7] = vcoor[1];
      Knetriangle[8] = vcoor[2];

      minv4gamTriangle = Knetriangle[5];

      g4taken_triangle = {indices[&combo - &combinations[0]][0], indices[&combo - &combinations[0]][1], indices[&combo - &combinations[0]][2], indices[&combo - &combinations[0]][3]};

      TrcFinalTriangle = (trc[0] + trc[1] + trc[2] + trc[3]);
    }
  }

  if (!success)
    return ErrorHandling::ErrorCodes::TRIANGLE_REC;

  return ErrorHandling::ErrorCodes::NO_ERROR;
}