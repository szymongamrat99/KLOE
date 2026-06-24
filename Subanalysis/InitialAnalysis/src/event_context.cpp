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
ErrorHandling::ErrorCodes EventContext::FilterNeutralClusters(int minClusters, double minEnergy, std::vector<Int_t> &outList)
{
  // 1. Odpalasz bazowe szukanie klastrów
  auto errorCode = fgenVarClassifier.FindNeutralCluster(
      fdataAccess.GetNClu(),
      fdataAccess.GetNTCl(),
      fdataAccess.GetAssCl().data(),
      minClusters,
      flogger,
      outList);

  if (errorCode != ErrorHandling::ErrorCodes::NO_ERROR)
    return errorCode;

  // 2. Filtrowanie po energii (jeśli minEnergy > 0)
  if (minEnergy > 0.0)
  {
    const auto &enecl = fdataAccess.GetEneCl();
    outList.erase(std::remove_if(outList.begin(), outList.end(),
                                 [&](Int_t idx)
                                 { return enecl[idx - 1] < minEnergy; }),
                  outList.end());

    if (static_cast<Int_t>(outList.size()) < minClusters)
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