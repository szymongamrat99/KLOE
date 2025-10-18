#include <GeneratedVariables.h>
#include <const.h>

void GeneratedVariables::classifyChannel(
    Int_t ntmc,
    Int_t nvtxmc,
    const Int_t *pidmcOld,
    const Int_t *vtxmcOld,
    const Int_t *motherOld,
    UInt_t mcflag,
    Int_t &mctruth_int)
{
  UInt_t Ks = 0, Kl = 0, Ksregen = 0, piplusks = 0, pipluskl = 0, piminusks = 0, piminuskl = 0,
         muonplusks = 0, muonpluskl = 0, muonminusks = 0, muonminuskl = 0, electronks = 0, electronkl = 0,
         positronks = 0, positronkl = 0, pi0ks = 0, pi0kl = 0, pi0phi = 0, piplusphi = 0, piminusphi = 0,
         otherphi = 0, otherkl = 0, otherks = 0, gammaphi = 0;

  if (mcflag == 1)
  {
    for (Int_t j = 0; j < ntmc; ++j)
    {
      if (motherOld[vtxmcOld[j] - 1] == 50)
      {
        switch (pidmcOld[j])
        {
        case 10:
          Kl++;
          break;
        case 16:
          Ks++;
          break;
        case 7:
          pi0phi++;
          break;
        case 8:
          piplusphi++;
          break;
        case 9:
          piminusphi++;
          break;
        case 1:
          gammaphi++;
          break;
        default:
          otherphi++;
          break;
        }
      }
      else if (motherOld[vtxmcOld[j] - 1] == 10)
      {
        switch (pidmcOld[j])
        {
        case 16:
          Ksregen++;
          break;
        case 7:
          pi0kl++;
          break;
        case 8:
          pipluskl++;
          break;
        case 9:
          piminuskl++;
          break;
        case 5:
          muonpluskl++;
          break;
        case 6:
          muonminuskl++;
          break;
        case 2:
          positronkl++;
          break;
        case 3:
          electronkl++;
          break;
        default:
          otherkl++;
          break;
        }
      }
      else if (motherOld[vtxmcOld[j] - 1] == 16)
      {
        switch (pidmcOld[j])
        {
        case 7:
          pi0ks++;
          break;
        case 8:
          piplusks++;
          break;
        case 9:
          piminusks++;
          break;
        case 5:
          muonplusks++;
          break;
        case 6:
          muonminusks++;
          break;
        case 2:
          positronks++;
          break;
        case 3:
          electronks++;
          break;
        default:
          otherks++;
          break;
        }
      }
    }

    Bool_t signal_cond = (pi0phi == 0 && piplusphi == 0 && piminusphi == 0 && otherphi == 0 && otherks == 0 && otherkl == 0 &&
                          positronkl + positronks == 0 && electronkl + electronks == 0 && muonminuskl + muonminusks == 0 &&
                          muonpluskl + muonplusks == 0 && Ksregen == 0 &&
                          Ks == 1 && Kl == 1 &&
                          ((pi0ks == 2 && pipluskl == 1 && piminuskl == 1 && pi0kl == 0 && piplusks == 0 && piminusks == 0) ||
                           (pi0kl == 2 && piplusks == 1 && piminusks == 1 && pi0ks == 0 && pipluskl == 0 && piminuskl == 0)));

    Bool_t pipi_cond = (pi0phi == 0 && piplusphi == 0 && piminusphi == 0 && otherphi == 0 && otherks == 0 && otherkl == 0 &&
                        positronkl + positronks == 0 && electronkl + electronks == 0 && muonminuskl + muonminusks == 0 &&
                        muonpluskl + muonplusks == 0 && Ksregen == 0 &&
                        Ks == 1 && Kl == 1 &&
                        ((piminusks == 1 && piplusks == 1 && pipluskl == 1 && piminuskl == 1 && pi0kl == 0 && pi0ks == 0)));

    Bool_t regen_cond = (Ksregen == 1 && Ks == 1 && Kl == 1);

    Bool_t omega_cond = (pi0phi == 2 && piplusphi == 1 && piminusphi == 1 && otherphi == 0 && otherks == 0 && otherkl == 0 &&
                         positronkl + positronks == 0 && electronkl + electronks == 0 && muonminuskl + muonminusks == 0 &&
                         muonpluskl + muonplusks == 0 && Ksregen == 0 &&
                         Ks == 0 && Kl == 0 && pi0ks == 0 && pi0kl == 0 && pipluskl + piplusks == 0 && piminuskl + piminusks == 0);

    Bool_t three_cond = (pi0phi == 0 && piplusphi == 0 && piminusphi == 0 && otherphi == 0 && otherks == 0 && otherkl == 0 &&
                         positronkl + positronks == 0 && electronkl + electronks == 0 && muonminuskl + muonminusks == 0 &&
                         muonpluskl + muonplusks == 0 && Ksregen == 0 &&
                         Ks == 1 && Kl == 1 && (pi0kl == 3 && piplusks == 1 && piminusks == 1 && pi0ks == 0 && pipluskl == 0 && piminuskl == 0));

    Bool_t semi_cond = (pi0phi == 0 && piplusphi == 0 && piminusphi == 0 && otherphi == 0 && otherks == 0 && otherkl == 0 &&
                        Ksregen == 0 && Ks == 1 && Kl == 1 &&
                        ((pi0ks == 2 && positronkl == 1 && piminuskl == 1 && pi0kl == 0) ||
                         (pi0ks == 2 && pipluskl == 1 && electronkl == 1 && pi0kl == 0) ||
                         (pi0ks == 2 && pipluskl == 1 && muonminuskl == 1 && pi0kl == 0) ||
                         (pi0ks == 2 && piminuskl == 1 && muonpluskl == 1 && pi0kl == 0) ||
                         (pi0kl == 2 && positronks == 1 && piminusks == 1 && pi0ks == 0) ||
                         (pi0kl == 2 && piplusks == 1 && electronks == 1 && pi0ks == 0) ||
                         (pi0kl == 2 && piplusks == 1 && muonminusks == 1 && pi0ks == 0) ||
                         (pi0kl == 2 && piminusks == 1 && muonplusks == 1 && pi0ks == 0)));

    if (signal_cond)
      mctruth_int = 1; // Signal channel: KSKL -> pi+pi-pi0pi0
    else if (regen_cond)
      mctruth_int = 2; // Regeneration
    else if (omega_cond)
      mctruth_int = 3; // omega
    else if (three_cond)
      mctruth_int = 4; // 3pi0
    else if (semi_cond)
      mctruth_int = 5; // semi
    else if (pipi_cond)
      mctruth_int = 7; // pipi
    else
      mctruth_int = 6; // Other background
  }
  else if (mcflag == 0)
  {
    mctruth_int = 0; // Data event
  }
}

ErrorHandling::ErrorCodes GeneratedVariables::FindNeutralCluster(
    Int_t nclu,
    Int_t ntcl,
    const Int_t *asscl,
    Int_t NCLMIN,
    ErrorHandling::ErrorLogs &logger,
    std::vector<Int_t> &neuclulist)
{
  neuclulist.clear();
  Int_t neucluind = 0;
  for (Int_t i = 1; i <= nclu; ++i)
  { // Fortran: 1-based
    Bool_t neuclu = true;
    for (Int_t j = 1; j <= ntcl; ++j)
    {
      if (asscl[j - 1] == i)
      { // C++: 0-based
        neuclu = false;
        break;
      }
    }
    if (neuclu)
    {
      ++neucluind;
      neuclulist.push_back(i);
    }
  }
  if (neucluind < NCLMIN)
  {
    auto err = ErrorHandling::ErrorCodes::NOT_RECOGNIZED;

    if (NCLMIN == 4)
      err = ErrorHandling::ErrorCodes::LESS_THAN_FOUR_NEUTRAL_CLUSTERS;
    else if (NCLMIN == 6)
      err = ErrorHandling::ErrorCodes::LESS_THAN_SIX_NEUTRAL_CLUSTERS;

    return err;
  }
  return ErrorHandling::ErrorCodes::NO_ERROR;
}

ErrorHandling::ErrorCodes GeneratedVariables::genVars(
    Int_t ntmc,
    Int_t nvtxmc,
    Int_t nclu,
    const Int_t *pidmc,
    const Int_t *vtxmc,
    const Int_t *mother,
    const Float_t *xvmc,
    const Float_t *yvmc,
    const Float_t *zvmc,
    const Float_t *pxmc,
    const Float_t *pymc,
    const Float_t *pzmc,
    Int_t mcflag,
    Int_t mctruth,
    std::vector<Float_t> &ipmc,
    std::vector<Float_t> &Knemc,
    std::vector<Float_t> &Kchmc,
    std::vector<std::vector<Float_t>> &trkMC,
    const Int_t numberOfClusters,
    std::vector<std::vector<Float_t>> &pgammaMC,
    std::vector<Float_t> &CurvMC,
    std::vector<Float_t> &PhivMC,
    std::vector<Float_t> &CotvMC,
    std::vector<Int_t> &good_clus_ind,
    std::vector<std::vector<Float_t>> cluster_rec)
{
  Float_t
      Ks[9] = {},
      Kl[9] = {},
      neu_vtx[3] = {},
      clus_diff_min;
  Int_t
      count = 0,
      ind_gam[4] = {},
      min_ind[4] = {};

  std::vector<Bool_t>
      clus_time;

  const Int_t
      max_count = TMath::Factorial(numberOfClusters);

  std::vector<Float_t>
      clus_diff,
      region,
      cluster;

  std::vector<Int_t>
      mc_ind(numberOfClusters);

  KLOE::CylinderIntersection CylIndObj;

  ipmc.clear();
  Knemc.clear();
  Kchmc.clear();
  trkMC.clear();

  cluster.clear();
  region.clear();
  clus_diff.clear();

  pgammaMC.clear();
  CurvMC.clear();
  PhivMC.clear();
  CotvMC.clear();

  ipmc.resize(3);
  Knemc.resize(9);
  Kchmc.resize(9);

  cluster.resize(3);

  IPGenerated(nvtxmc, mother, ipmc, xvmc, yvmc, zvmc);

  if (mctruth == 1 || mctruth == 3 || mctruth == 4 || mctruth == 5 || mctruth == 7)
  {
    KSLGenerated(nvtxmc, mother, Kl, xvmc, yvmc, zvmc, Ks, ntmc, pidmc, pxmc, pymc, pzmc);

    twoTracksFinder(ntmc, mother, vtxmc, pidmc, Knemc, Kl, Kchmc, Ks, trkMC, pxmc, pymc, pzmc, mctruth, CurvMC, PhivMC, CotvMC);

    if (mctruth == 1 || mctruth == 3 || mctruth == 4 || mctruth == 5)
    {
      ClusterVariableFinder(ntmc, mother, vtxmc, pidmc, pgammaMC, count, pxmc, pymc, pzmc, neu_vtx, Knemc, region, CylIndObj, cluster, ipmc);
    }
  };

  count = 0;

  return ErrorHandling::ErrorCodes::NO_ERROR;
}
void GeneratedVariables::ClusterVariableFinder(Int_t ntmc, const Int_t *mother, const Int_t *vtxmc, const Int_t *pidmc, std::vector<std::vector<Float_t>> &pgammaMC, Int_t &count, const Float_t *pxmc, const Float_t *pymc, const Float_t *pzmc, Float_t neu_vtx[3], std::vector<Float_t> &Knemc, std::vector<Float_t> &region, KLOE::CylinderIntersection &CylIndObj, std::vector<Float_t> &cluster, std::vector<Float_t> &ipmc)
{
  for (Int_t j = 0; j < ntmc; j++)
  {
    if ((mother[vtxmc[j] - 1] == 7) && pidmc[j] == 1)
    {
      Float_t auxEne = sqrt(pow(pxmc[j], 2) + pow(pymc[j], 2) + pow(pzmc[j], 2));
      std::vector<Float_t> pgammaAux = {pxmc[j], pymc[j], pzmc[j], auxEne};

      neu_vtx[0] = Knemc[6];
      neu_vtx[1] = Knemc[7];
      neu_vtx[2] = Knemc[8];

      region.push_back(CylIndObj.inter_point(pgammaAux.data(), neu_vtx, cluster.data()));

      Float_t
          beta_c = PhysicsConstants::cVel * Knemc[4] / Knemc[3],
          length = sqrt(pow(Knemc[6] - ipmc[0], 2) +
                        pow(Knemc[7] - ipmc[1], 2) +
                        pow(Knemc[8] - ipmc[2], 2)),
          time_K = length / beta_c,
          length_clus = sqrt(pow(cluster[0] - Knemc[6], 2) +
                             pow(cluster[1] - Knemc[7], 2) +
                             pow(cluster[2] - Knemc[8], 2));

      Float_t auxTim = time_K + (length_clus / PhysicsConstants::cVel);
      std::vector<Float_t> auxiliaryVec = {pxmc[j], pymc[j], pzmc[j], auxEne, cluster[0], cluster[1], cluster[2], auxTim};

      pgammaMC.push_back(auxiliaryVec);
    }
  }
}

void GeneratedVariables::GeneratedClusterFinder(Int_t nclu, Int_t ind_gam[4], const Int_t max_count, std::vector<Float_t> &clus_diff, std::vector<std::vector<Float_t>> &cluster_rec, std::vector<std::vector<Float_t>> &pgammaMC, Int_t mc_ind[4], std::vector<bool> &clus_time, Int_t min_ind[4], Float_t &clus_diff_min, std::vector<Int_t> &good_clus_ind)
{
  for (Int_t j1 = 0; j1 < nclu - 3; j1++)
    for (Int_t j2 = j1 + 1; j2 < nclu - 2; j2++)
      for (Int_t j3 = j2 + 1; j3 < nclu - 1; j3++)
        for (Int_t j4 = j3 + 1; j4 < nclu; j4++)
        {
          ind_gam[0] = j1;
          ind_gam[1] = j2;
          ind_gam[2] = j3;
          ind_gam[3] = j4;

          for (Int_t k = 0; k < max_count; k++)
          {

            clus_diff.push_back(sqrt(pow(cluster_rec[0][ind_gam[0]] - pgammaMC[mc_ind[0]][4], 2) +
                                     pow(cluster_rec[1][ind_gam[0]] - pgammaMC[mc_ind[0]][5], 2) +
                                     pow(cluster_rec[2][ind_gam[0]] - pgammaMC[mc_ind[0]][6], 2)) +
                                sqrt(pow(cluster_rec[0][ind_gam[1]] - pgammaMC[mc_ind[1]][4], 2) +
                                     pow(cluster_rec[1][ind_gam[1]] - pgammaMC[mc_ind[1]][5], 2) +
                                     pow(cluster_rec[2][ind_gam[1]] - pgammaMC[mc_ind[1]][6], 2)) +
                                sqrt(pow(cluster_rec[0][ind_gam[2]] - pgammaMC[mc_ind[2]][4], 2) +
                                     pow(cluster_rec[1][ind_gam[2]] - pgammaMC[mc_ind[2]][5], 2) +
                                     pow(cluster_rec[2][ind_gam[2]] - pgammaMC[mc_ind[2]][6], 2)) +
                                sqrt(pow(cluster_rec[0][ind_gam[3]] - pgammaMC[mc_ind[3]][4], 2) +
                                     pow(cluster_rec[1][ind_gam[3]] - pgammaMC[mc_ind[3]][5], 2) +
                                     pow(cluster_rec[2][ind_gam[3]] - pgammaMC[mc_ind[3]][6], 2)));

            clus_time.push_back(pgammaMC[mc_ind[0]][7] > 0. && pgammaMC[mc_ind[1]][7] > 0. && pgammaMC[mc_ind[2]][7] > 0. && pgammaMC[mc_ind[3]][7] > 0.);

            std::next_permutation(mc_ind, mc_ind + 4);
          }

          TMath::Sort(max_count, clus_diff.data(), min_ind, kFALSE);

          if (clus_diff_min > clus_diff[min_ind[0]])
          {
            clus_diff_min = clus_diff[min_ind[0]];

            good_clus_ind[0] = ind_gam[0];
            good_clus_ind[1] = ind_gam[1];
            good_clus_ind[2] = ind_gam[2];
            good_clus_ind[3] = ind_gam[3];
          }
        }
}

void GeneratedVariables::MCvsReconstructedClustersComparator(const std::vector<Int_t> neuclulist, const std::vector<Int_t> gtaken, const std::vector<Int_t> Pnum1, const Int_t ntmc, const std::vector<Int_t> mother, const std::vector<Int_t> vtxmc, const std::vector<Int_t> pidmc, const std::vector<Int_t> kine, const std::vector<Int_t> kinmom, std::vector<Int_t> &goodCluster)
{
  goodCluster.clear();
  goodCluster.resize(0);

  for (Int_t i = 0; i < ntmc; ++i)
  {
    for (Int_t j = 0; j < gtaken.size(); ++j)
    {
      if (kine[i] == Pnum1[neuclulist[gtaken[j]] - 1])
      {
        if (mother[vtxmc[i] - 1] == 7 && pidmc[i] == 1)
        {
          Int_t kinPi0 = kinmom[vtxmc[i] - 1];
          for (Int_t k = 0; k < ntmc; ++k)
          {
            if(pidmc[k] == 7 && kine[k] == kinPi0 && (mother[vtxmc[k] - 1] == 10 || mother[vtxmc[k] - 1] == 16))
            {
              goodCluster.push_back(neuclulist[gtaken[j]]);

              // std::cout << "DEBUG: Found good cluster " << neuclulist[gtaken[j]] << " for pi0 with kine " << kinPi0 << " from mother " << mother[vtxmc[k] - 1] << " and mother's kine " << kinmom[vtxmc[k] - 1] << std::endl;
              break;
            }
          }
        }
      }
    }
  }
}

void GeneratedVariables::twoTracksFinder(Int_t ntmc, const Int_t *mother, const Int_t *vtxmc, const Int_t *pidmc, std::vector<Float_t> &Knemc, Float_t Kl[9], std::vector<Float_t> &Kchmc, Float_t Ks[9], std::vector<std::vector<Float_t>> &trkMC, const Float_t *pxmc, const Float_t *pymc, const Float_t *pzmc, Int_t mctruth, std::vector<Float_t> &CurvMC, std::vector<Float_t> &PhivMC, std::vector<Float_t> &CotvMC)
{
  for (Int_t j = 0; j < ntmc; j++)
  {
    if (mother[vtxmc[j] - 1] == 10)
    {
      if (mctruth != 7)
      {
        if (pidmc[j] == 7)
        {
          std::copy(Ks, Ks + 9, Kchmc.begin());
          std::copy(Kl, Kl + 9, Knemc.begin());
        }
        else if (pidmc[j] == 8 || pidmc[j] == 9)
        {
          std::copy(Ks, Ks + 9, Knemc.begin());
          std::copy(Kl, Kl + 9, Kchmc.begin());
        }
      }
      else
      {
        std::copy(Ks, Ks + 9, Kchmc.begin());
        std::copy(Kl, Kl + 9, Knemc.begin());
      }
    }

    if ((mother[vtxmc[j] - 1] == 10) && (pidmc[j] == 8 || pidmc[j] == 9))
    {
      Float_t auxEne = sqrt(pow(pxmc[j], 2) +
                            pow(pymc[j], 2) +
                            pow(pzmc[j], 2) +
                            pow(PhysicsConstants::mPiCh, 2));
      std::vector<Float_t> auxiliaryVec = {pxmc[j], pymc[j], pzmc[j], auxEne, 10};

      trkMC.push_back(auxiliaryVec);

      CurvMC.push_back(1000. / sqrt(pow(pxmc[j], 2) + pow(pymc[j], 2)));
      PhivMC.push_back(atan2(pymc[j], pxmc[j]));
      CotvMC.push_back(pzmc[j] / sqrt(pow(pxmc[j], 2) + pow(pymc[j], 2)));

      if (pidmc[j] == 9)
      {
        CurvMC.back() = -CurvMC.back();
      }
    }

    if ((mother[vtxmc[j] - 1] == 16) && (pidmc[j] == 8 || pidmc[j] == 9))
    {
      Float_t auxEne = sqrt(pow(pxmc[j], 2) +
                            pow(pymc[j], 2) +
                            pow(pzmc[j], 2) +
                            pow(PhysicsConstants::mPiCh, 2));
      std::vector<Float_t> auxiliaryVec = {pxmc[j], pymc[j], pzmc[j], auxEne, 16};

      trkMC.push_back(auxiliaryVec);

      CurvMC.push_back(1000. / sqrt(pow(pxmc[j], 2) + pow(pymc[j], 2)));
      PhivMC.push_back(atan2(pymc[j], pxmc[j]));
      CotvMC.push_back(pzmc[j] / sqrt(pow(pxmc[j], 2) + pow(pymc[j], 2)));

      if (pidmc[j] == 9)
      {
        CurvMC.back() = -CurvMC.back();
      }
    }
  }
}

void GeneratedVariables::KSLGenerated(Int_t nvtxmc, const Int_t *mother, Float_t Kl[9], const Float_t *xvmc, const Float_t *yvmc, const Float_t *zvmc, Float_t Ks[9], Int_t ntmc, const Int_t *pidmc, const Float_t *pxmc, const Float_t *pymc, const Float_t *pzmc)
{
  for (Int_t j = 0; j < nvtxmc; j++)
  {
    if (mother[j] == 10)
    {
      Kl[6] = xvmc[j];
      Kl[7] = yvmc[j];
      Kl[8] = zvmc[j];
    }
  }

  for (Int_t j = 0; j < nvtxmc; j++)
  {
    if (mother[j] == 16)
    {
      Ks[6] = xvmc[j];
      Ks[7] = yvmc[j];
      Ks[8] = zvmc[j];
    }
  }

  for (Int_t j = 0; j < ntmc; j++)
  {
    if (pidmc[j] == 10)
    {
      Kl[0] = pxmc[j];
      Kl[1] = pymc[j];
      Kl[2] = pzmc[j];
      Kl[5] = PhysicsConstants::mK0;
      Kl[4] = pow(Kl[0], 2) + pow(Kl[1], 2) + pow(Kl[2], 2);
      Kl[3] = sqrt(Kl[4] + pow(Kl[5], 2));
      Kl[4] = sqrt(Kl[4]);
    }
  }

  for (Int_t j = 0; j < ntmc; j++)
  {
    if (pidmc[j] == 16)
    {
      Ks[0] = pxmc[j];
      Ks[1] = pymc[j];
      Ks[2] = pzmc[j];
      Ks[5] = PhysicsConstants::mK0;
      Ks[4] = pow(Ks[0], 2) + pow(Ks[1], 2) + pow(Ks[2], 2);
      Ks[3] = sqrt(Ks[4] + pow(Ks[5], 2));
      Ks[4] = sqrt(Ks[4]);
    }
  }
}
void GeneratedVariables::IPGenerated(Int_t nvtxmc, const Int_t *mother, std::vector<Float_t> &ipmc, const Float_t *xvmc, const Float_t *yvmc, const Float_t *zvmc)
{
  for (Int_t j = 0; j < nvtxmc; j++)
  {
    if (mother[j] == 50)
    {
      ipmc[0] = xvmc[j];
      ipmc[1] = yvmc[j];
      ipmc[2] = zvmc[j];
    }
  }
};