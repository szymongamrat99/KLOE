#ifndef GENERATED_VARIABLES_H
#define GENERATED_VARIABLES_H

#include <TMath.h>
#include <ErrorLogs.h>

class GeneratedVariables
{
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
};

#endif // GENERATED_VARIABLES_H
