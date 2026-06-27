#include "../inc/event_context.h"

int EventContext::GetNumberOfVerticesWithTwoTracks()
{
  fVtxTracksAssignedNum = fobj.CountRepeatingElements(fBaseKin.iv);

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
      fBaseKin.neuclulist);

  if (errorCode != ErrorHandling::ErrorCodes::NO_ERROR)
    return errorCode;

  // 2. Filtrowanie po energii (jeśli minEnergy > 0)
  if (minClusterEnergy > 0.0)
  {
    const auto &enecl = fdataAccess.GetEneCl();
    fBaseKin.neuclulist.erase(std::remove_if(fBaseKin.neuclulist.begin(), fBaseKin.neuclulist.end(),
                                     [&](Int_t idx)
                                 { return enecl[idx - 1] < minClusterEnergy; }),
                  fBaseKin.neuclulist.end());

    if (static_cast<Int_t>(fBaseKin.neuclulist.size()) < minClusters)
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

void EventContext::ProcessMonteCarloTruth(bool isMcEnabled,
                                          int mctruthSignal,
                                          std::array<UInt_t, 8> &mctruth_num)
{
  if (!isMcEnabled)
  {
    // -------------------------------------------------------------------
    // REAL DATA: Resetowanie struktur do wartości domyślnych
    // -------------------------------------------------------------------
    fBaseKin.clearMC(); // Warto przenieść ten wielki reset do metody w BaseKinematics!
    return;
  }

  // -------------------------------------------------------------------
  // MONTE CARLO: Pełna rekonstrukcja prawdy generatora
  // -------------------------------------------------------------------

  // Definiujemy tymczasowe kontenery (zamiast trzymać je globalnie w makrze)
  std::vector<std::vector<Double_t>> trkMC;
  std::vector<std::vector<Double_t>> pgammaMC;
  std::vector<std::vector<Double_t>> clusterMC;

  // 1. Klasyfikacja kanału MC
  fgenVarClassifier.classifyChannel(
      fdataAccess.GetNTMC(), fdataAccess.GetNVtxMC(), fdataAccess.GetPidMC().data(),
      fdataAccess.GetVtxMC().data(), fdataAccess.GetMother().data(),
      fBaseKin.mcflag, fBaseKin.mctruth, fBaseKin.semileptonic_flag);

  try
  {
    if (fBaseKin.mctruth >= 0 && fBaseKin.mctruth < 8)
    {
      mctruth_num[fBaseKin.mctruth]++;
    }
    else
    {
      throw ErrorHandling::ErrorCodes::RANGE;
    }
  }
  catch (ErrorHandling::ErrorCodes &e)
  {
    std::string message = "MC truth value out of range: " + std::to_string(fBaseKin.mctruth);
    LOG_CRITICAL(flogger, e, message, ErrorHandling::LogFiles::LogType::ERROR);
  }
  KLOE::channEventCount[fBaseKin.mctruth]++;

  // 2. Generowanie zmiennych MC
  fgenVarClassifier.genVars(
      fBaseKin.ntmc, fBaseKin.nvtxmc, fBaseKin.nclu, fBaseKin.pidmc.data(), fBaseKin.vtxmc.data(),
      fBaseKin.mother.data(), fBaseKin.xvmc.data(), fBaseKin.yvmc.data(), fBaseKin.zvmc.data(),
      fBaseKin.pxmc.data(), fBaseKin.pymc.data(), fBaseKin.pzmc.data(), fBaseKin.mcflag, fBaseKin.mctruth,
      fBaseKin.ipmc, fBaseKin.Knemc, fBaseKin.Kchmc, trkMC, 4, pgammaMC,
      fBaseKin.CurvMC, fBaseKin.PhivMC, fBaseKin.CotvMC, fBaseKin.goodClusIndex, clusterMC,
      fBaseKin.muonAlertPlus, fBaseKin.muonAlertMinus);

  // 3. Obliczanie czasów życia kaonów
  std::vector<Double_t> kaonChMom = {fBaseKin.Kchmc[0], fBaseKin.Kchmc[1], fBaseKin.Kchmc[2], fBaseKin.Kchmc[3]};
  std::vector<Double_t> kaonChPos = {fBaseKin.Kchmc[6], fBaseKin.Kchmc[7], fBaseKin.Kchmc[8]};
  std::vector<Double_t> kaonNeMom = {fBaseKin.Knemc[0], fBaseKin.Knemc[1], fBaseKin.Knemc[2], fBaseKin.Knemc[3]};
  std::vector<Double_t> kaonNePos = {fBaseKin.Knemc[6], fBaseKin.Knemc[7], fBaseKin.Knemc[8]};
  std::vector<Double_t> ipPos = {fBaseKin.ipmc[0], fBaseKin.ipmc[1], fBaseKin.ipmc[2]};

  fBaseKin.kaonTimesMC = fobj.CalculateKaonProperTimes(kaonChMom, kaonChPos, kaonNeMom, kaonNePos, ipPos);

  // 4. Specyficzna logika per-kanał (Rozbicie na fmctruth)
  if (fBaseKin.mctruth == 1 && trkMC.size() >= 2)
  {
    std::copy_n(trkMC[0].begin(), 4, fBaseKin.trk1MC.begin());
    std::copy_n(trkMC[1].begin(), 4, fBaseKin.trk2MC.begin());
  }
  else if (fBaseKin.mctruth == 7)
  {
    for (const auto &singleTrk : trkMC)
    {
      if (singleTrk.size() < 5)
        continue;
      int pdgType = singleTrk[4];

      if (pdgType == 10) // KL
      {
        auto targetIdx = (fBaseKin.trkKLmc[0].empty()) ? 0 : 1;
        fBaseKin.trkKLmc[targetIdx].assign(singleTrk.begin(), singleTrk.end() - 1);
      }
      else if (pdgType == 16) // KS
      {
        auto targetIdx = (fBaseKin.trkKSmc[0].empty()) ? 0 : 1;
        fBaseKin.trkKSmc[targetIdx].assign(singleTrk.begin(), singleTrk.end() - 1);
      }
    }
  }
}

ErrorHandling::ErrorCodes EventContext::ReconstructKaonIntoChargedPions()
{
  TMatrixT<Double_t> covMatrixTot(6, 6);
  covMatrixTot.Zero();

  return fCurrentEventAnalysis->findKchRec(fBaseKin.mcflag, 0, covMatrixTot, fBaseKin.Kchrecnew, fBaseKin.trknew[0], fBaseKin.trknew[1], fBaseKin.vtaken, flogger, 1);
}

ErrorHandling::ErrorCodes EventContext::ReconstructKSIntoChargedParticles()
{
  TMatrixT<Double_t> covMatrixTot(6, 6);
  covMatrixTot.Zero();

  return fCurrentEventAnalysis->findKSLRec(16, -1, fBaseKin.KchrecKS, fBaseKin.trkKS[0], fBaseKin.trkKS[1], fBaseKin.vtakenKS, flogger);
}

ErrorHandling::ErrorCodes EventContext::ReconstructKLIntoChargedParticles()
{
  TMatrixT<Double_t> covMatrixTot(6, 6);
  covMatrixTot.Zero();

  return fCurrentEventAnalysis->findKSLRec(10, fBaseKin.vtakenKS[0], fBaseKin.KchrecKS, fBaseKin.trkKS[0], fBaseKin.trkKS[1], fBaseKin.vtakenKS, flogger);
}

ErrorHandling::ErrorCodes EventContext::ReconstructKaonClosestToIP()
{
  TMatrixT<Double_t> covMatrixTot(6, 6);
  covMatrixTot.Zero();

  return fCurrentEventAnalysis->findKClosestRec(fBaseKin.KchrecClosest, fBaseKin.trkClosest[0], fBaseKin.trkClosest[1], fBaseKin.vtakenClosest, flogger);
}

