#ifndef GENERATED_VARIABLES_H
#define GENERATED_VARIABLES_H

#include <TMath.h>

class GeneratedVariables {
public:
    /**
     * @brief Classifies the MC channel and sets mctruth_int accordingly.
     * @param ntmc Number of MC particles
     * @param nvtxmc Number of MC vertices (not used in logic, but kept for compatibility)
     * @param pidmc Array of MC particle IDs (size >= ntmc)
     * @param vtxmc Array of MC vertex indices (size >= ntmc)
     * @param mother Array of MC mother IDs (size >= ntmc)
     * @param mcflag MC flag (1 = MC, 0 = data)
     * @param mctruth MC truth code (input from MC, used in logic)
     * @param[out] mctruth_int Output: classified channel code
     */
    static void classifyChannel(
        Int_t ntmc,
        Int_t nvtxmc,
        Int_t* pidmc,
        Int_t* vtxmc,
        Int_t* mother,
        UInt_t mcflag,
        Int_t& mctruth
    );
};

#endif // GENERATED_VARIABLES_H
