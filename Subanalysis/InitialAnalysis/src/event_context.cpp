#include "../inc/event_context.h"

int EventContext::GetNumberOfVerticesWithTwoTracks()
{
  if (fIv.empty())
    fIv.assign(fdataAccess.GetIv().begin(), fdataAccess.GetIv().end());

  fVtxTracksAssignedNum = fobj.CountRepeatingElements(fIv);

  int verticesWithTwoTracksCount = 0;
  for (const auto &pair : fVtxTracksAssignedNum)
  {
    if (pair.second == 2)
    {
      verticesWithTwoTracksCount++;
    }
  }

  return verticesWithTwoTracksCount;
}

// Zwraca kod błędu lub listę przefiltrowanych klastrów
ErrorHandling::ErrorCodes EventContext::FilterNeutralClusters(int minClusters, double minClusterEnergy)
{
  // 1. Odpalasz bazowe szukanie klastrów
  auto errorCode = fgenVarClassifier.FindNeutralCluster(
      fdataAccess.GetNClu(),
      fdataAccess.GetNTCl(),
      fdataAccess.GetAssCl().data(),
      minClusters,
      flogger,
      fNeuCluList);

  if (errorCode != ErrorHandling::ErrorCodes::NO_ERROR)
    return errorCode;

  // 2. Filtrowanie po energii (jeśli minEnergy > 0)
  if (minClusterEnergy > 0.0)
  {
    const auto &enecl = fdataAccess.GetEneCl();
    fNeuCluList.erase(std::remove_if(fNeuCluList.begin(), fNeuCluList.end(),
                                     [&](Int_t idx)
                                 { return enecl[idx - 1] < minClusterEnergy; }),
                  fNeuCluList.end());

    if (static_cast<Int_t>(fNeuCluList.size()) < minClusters)
    {
      if (minClusters == 4)
      {
        return ErrorHandling::ErrorCodes::LESS_THAN_FOUR_CLUSTERS_WITH_GOOD_ENERGY;
      }
      else if (minClusters == 6)
      {
        return ErrorHandling::ErrorCodes::LESS_THAN_SIX_NEUTRAL_CLUSTERS;
      }
      else
      {
        return ErrorHandling::ErrorCodes::LESS_THAN_FOUR_NEUTRAL_CLUSTERS;
      }
    }
  }

  return ErrorHandling::ErrorCodes::NO_ERROR;
}

// ---

void EventContext::ProcessMonteCarloTruth(bool isMcEnabled, BaseKinematics &baseKin,
                                          int &mcflag, int &mctruth, int mctruthSignal,
                                          KLOE::KaonProperTimes &kaonTimesMC,
                                          std::array<UInt_t, 8> &mctruth_num)
{
  if (!isMcEnabled)
  {
    // -------------------------------------------------------------------
    // REAL DATA: Resetowanie struktur do wartości domyślnych
    // ----------------------------------------------------------------===
    mcflag = 0;
    mctruth = 0;
    baseKin.clearMC(); // Warto przenieść ten wielki reset do metody w BaseKinematics!
    kaonTimesMC = KLOE::KaonProperTimes();
    return;
  }

  // -------------------------------------------------------------------
  // MONTE CARLO: Pełna rekonstrukcja prawdy generatora
  // -------------------------------------------------------------------
  mcflag = 1;

  // Definiujemy tymczasowe kontenery (zamiast trzymać je globalnie w makrze)
  std::vector<std::vector<Double_t>> trkMC;
  std::vector<std::vector<Double_t>> pgammaMC;
  std::vector<std::vector<Double_t>> clusterMC;

  // 1. Klasyfikacja kanału MC
  fgenVarClassifier.classifyChannel(
      fdataAccess.GetNTMC(), fdataAccess.GetNVtxMC(), fdataAccess.GetPidMC().data(),
      fdataAccess.GetVtxMC().data(), fdataAccess.GetMother().data(),
      mcflag, mctruth, baseKin.semileptonic_flag);

  try
  {
    if (mctruth >= 0 && mctruth < 8)
    {
      mctruth_num[mctruth]++;
    }
    else
    {
      throw ErrorHandling::ErrorCodes::RANGE;
    }
  }
  catch (ErrorHandling::ErrorCodes &e)
  {
    std::string message = "MC truth value out of range: " + std::to_string(mctruth);
    LOG_CRITICAL(flogger, e, message, ErrorHandling::LogFiles::LogType::ERROR);
  }
  KLOE::channEventCount[mctruth]++;

  // 2. Generowanie zmiennych MC
  fgenVarClassifier.genVars(
      baseKin.ntmc, baseKin.nvtxmc, baseKin.nclu, baseKin.pidmc.data(), baseKin.vtxmc.data(),
      baseKin.mother.data(), baseKin.xvmc.data(), baseKin.yvmc.data(), baseKin.zvmc.data(),
      baseKin.pxmc.data(), baseKin.pymc.data(), baseKin.pzmc.data(), mcflag, mctruth,
      baseKin.ipmc, baseKin.Knemc, baseKin.Kchmc, trkMC, 4, pgammaMC,
      baseKin.CurvMC, baseKin.PhivMC, baseKin.CotvMC, baseKin.goodClusIndex, clusterMC,
      baseKin.muonAlertPlus, baseKin.muonAlertMinus);

  // 3. Obliczanie czasów życia kaonów
  std::vector<Double_t> kaonChMom = {baseKin.Kchmc[0], baseKin.Kchmc[1], baseKin.Kchmc[2], baseKin.Kchmc[3]};
  std::vector<Double_t> kaonChPos = {baseKin.Kchmc[6], baseKin.Kchmc[7], baseKin.Kchmc[8]};
  std::vector<Double_t> kaonNeMom = {baseKin.Knemc[0], baseKin.Knemc[1], baseKin.Knemc[2], baseKin.Knemc[3]};
  std::vector<Double_t> kaonNePos = {baseKin.Knemc[6], baseKin.Knemc[7], baseKin.Knemc[8]};
  std::vector<Double_t> ipPos = {baseKin.ipmc[0], baseKin.ipmc[1], baseKin.ipmc[2]};

  kaonTimesMC = fobj.CalculateKaonProperTimes(kaonChMom, kaonChPos, kaonNeMom, kaonNePos, ipPos);

  // 4. Specyficzna logika per-kanał (Rozbicie na mctruth)
  if (mctruth == 1 && trkMC.size() >= 2)
  {
    std::copy_n(trkMC[0].begin(), 4, baseKin.trk1MC.begin());
    std::copy_n(trkMC[1].begin(), 4, baseKin.trk2MC.begin());
  }
  else if (mctruth == 7)
  {
    for (const auto &singleTrk : trkMC)
    {
      if (singleTrk.size() < 5)
        continue;
      int pdgType = singleTrk[4];

      if (pdgType == 10) // KL
      {
        auto targetIdx = (baseKin.trkKLmc[0].empty()) ? 0 : 1;
        baseKin.trkKLmc[targetIdx].assign(singleTrk.begin(), singleTrk.end() - 1);
      }
      else if (pdgType == 16) // KS
      {
        auto targetIdx = (baseKin.trkKSmc[0].empty()) ? 0 : 1;
        baseKin.trkKSmc[targetIdx].assign(singleTrk.begin(), singleTrk.end() - 1);
      }
    }
  }
}