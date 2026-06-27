#include <three_pi0_hypothesis.h>

ErrorHandling::ErrorCodes ThreePi0Hypothesis::Process(EventContext &ctx)
{
  // --- Get the number of vertices with two tracks ---
  int numVerticesWithTwoTracks = ctx.GetNumberOfVerticesWithTwoTracks();
  if (numVerticesWithTwoTracks < 1)
    errorCode = ErrorHandling::ErrorCodes::NO_VTX_WITH_TWO_TRACKS;

  noError = CheckErrorCode(errorCode);
  // --------------------------------------------------

  // --- Filter neutral clusters for the 3pi0 hypothesis ---
  if (noError)
  {
    errorCode = ctx.FilterNeutralClusters(minClusters, minClusterEnergy);
    noError = CheckErrorCode(errorCode);
  }
  // -------------------------------------------------------

  // --- Reconstruct kaon into charged pions with invariant mass condition ---
  if (noError)
  {
    errorCode = ctx.ReconstructKaonIntoChargedPions();
    noError = CheckErrorCode(errorCode);
  }
  // -------------------------------------------------------------------------

  // --- Reconstruct Kaon into neutral particles -----------------------------
  if (noError)
  {
    errorCode = ctx.ReconstructKaonClosestToIP();
    noError = CheckErrorCode(errorCode);
  }



  return errorCode;
}