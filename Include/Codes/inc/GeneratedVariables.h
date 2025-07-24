#ifndef GENERATED_VARIABLES_H
#define GENERATED_VARIABLES_H

#include <TMath.h>
#include <ErrorLogs.h>
#include <cylinder_intersection.h>

class GeneratedVariables
{
private:
	void IPGenerated(Int_t nvtxmc, Int_t *mother, std::vector<Float_t> &ipmc, Float_t *xvmc, Float_t *yvmc, Float_t *zvmc);
	void KSLGenerated(Int_t nvtxmc, Int_t *mother, Float_t Kl[9], Float_t *xvmc, Float_t *yvmc, Float_t *zvmc, Float_t Ks[9], Int_t ntmc, Int_t *pidmc, Float_t *pxmc, Float_t *pymc, Float_t *pzmc);
	void twoTracksFinder(Int_t ntmc, Int_t *mother, Int_t *vtxmc, Int_t *pidmc, std::vector<Float_t> &Knemc, Float_t Kl[9], std::vector<Float_t> &Kchmc, Float_t Ks[9], std::vector<std::vector<Float_t>> &trkMC, Float_t *pxmc, Float_t *pymc, Float_t *pzmc, Int_t mctruth);
	void GeneratedClusterFinder(Int_t nclu, Int_t ind_gam[4], const Int_t max_count, std::vector<Float_t> &clus_diff, std::vector<std::vector<Float_t>> &cluster_rec, std::vector<std::vector<Float_t>> &pgammaMC, Int_t mc_ind[4], std::vector<bool> &clus_time, Int_t min_ind[4], Float_t &clus_diff_min, std::vector<Int_t> &good_clus_ind);
	void ClusterVariableFinder(Int_t ntmc, Int_t *mother, Int_t *vtxmc, Int_t *pidmc, std::vector<std::vector<Float_t>> &pgammaMC, Int_t &count, Float_t *pxmc, Float_t *pymc, Float_t *pzmc, Float_t neu_vtx[3], std::vector<Float_t> &Knemc, std::vector<Float_t> &region, KLOE::CylinderIntersection &CylIndObj, std::vector<Float_t> &cluster, std::vector<Float_t> &ipmc);

public:
	/**
	 * @brief Classifies the MC channel and sets mctruth_int accordingly.
	 * @param ntmc Number of MC particles
	 * @param nvtxmc Number of MC vertices (not used in logic, but kept for compatibility)
	 * @param pidmcOld Array of MC particle IDs (size >= ntmc)
	 * @param vtxmcOld Array of MC vertex indices (size >= ntmc)
	 * @param motherOld Array of MC motherOld IDs (size >= ntmc)
	 * @param mcflag MC flag (1 = MC, 0 = data)
	 * @param mctruth MC truth code (input from MC, used in logic)
	 * @param[out] mctruth_int Output: classified channel code
	 */
	static void classifyChannel(
			Int_t ntmc,
			Int_t nvtxmc,
			Int_t *pidmcOld,
			Int_t *vtxmcOld,
			Int_t *motherOld,
			UInt_t mcflag,
			Int_t &mctruth);

	/**
	 * @brief Find clusters not associated with any track.
	 * @param nclu Number of clusters (Fortran: 1-based, C++: 1..nclu)
	 * @param ntcl Number of tracks
	 * @param asscl Array of associations (size ntcl), asscl[j] = cluster index (Fortran: 1-based)
	 * @param NCLMIN Minimal number of clusters required
	 * @param logger Error logger (ErrorHandling::ErrorLogs)
	 * @param neuclulist Output: vector of cluster indices (1-based) not associated with any track
	 * @return true if error (less than NCLMIN clusters found), false otherwise
	 */
	static ErrorHandling::ErrorCodes FindNeutralCluster(
			Int_t nclu,
			Int_t ntcl,
			const Int_t *asscl,
			Int_t NCLMIN,
			ErrorHandling::ErrorLogs &logger,
			std::vector<Int_t> &neuclulist);

  ErrorHandling::ErrorCodes genVars(
      Int_t ntmc,
      Int_t nvtxmc,
      Int_t nclu,
      Int_t *pidmc,
      Int_t *vtxmc,
      Int_t *mother,
      Float_t *xvmc,
      Float_t *yvmc,
      Float_t *zvmc,
      Float_t *pxmc,
      Float_t *pymc,
      Float_t *pzmc,
      Int_t mcflag,
      Int_t mctruth,
      std::vector<Float_t> &ipmc,
      std::vector<Float_t> &Knemc,
      std::vector<Float_t> &Kchmc,
      std::vector<std::vector<Float_t>> &trkMC,
      const Int_t numberOfClusters,
      std::vector<std::vector<Float_t>> &pgammaMC,
      std::vector<Int_t> &good_clus_ind,
      std::vector<std::vector<Float_t>> cluster_rec);
};

#endif // GENERATED_VARIABLES_H
