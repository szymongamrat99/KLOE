#ifndef KAON_MASS_RECONSTRUCTOR_H
#define KAON_MASS_RECONSTRUCTOR_H

#include <TVector3.h>
#include <TLorentzVector.h>
#include <vector>
#include "klspm00.hpp"

namespace KLOE {

struct KaonReconstructionResult {
    std::vector<Float_t> KaonTwoBody;     // 9 elements: px,py,pz,E,|p|,m,vx,vy,vz
    std::vector<Float_t> track1TwoBody;   // 4 elements: px,py,pz,E
    std::vector<Float_t> track2TwoBody;   // 4 elements: px,py,pz,E
};

class KaonMassReconstructor {
public:
    static KaonReconstructionResult reconstructKaonMass(
        const std::vector<Float_t>& kaonBoost,    // 9 elements: px,py,pz,E,unused,unused,unused,vx,vy,vz
        const std::vector<Float_t>& interactionPoint, // 3 elements: x,y,z
        const std::vector<std::vector<Float_t>*>& tracks, // 2 tracks (pointers), each 4 elements: px,py,pz,E
        const Float_t mPiCh,  // charged pion mass
        KLOE::pm00& obj       // KLOE object for utility functions
    );

private:
    static std::vector<Float_t> calculateKaonMomentum(
        const std::vector<Float_t>& kaonBoost,
        const std::vector<Float_t>& interactionPoint
    );

    static void transformToKaonCMFrame(
        TVector3& kaonMomLAB,
        const std::vector<std::vector<Float_t>*>& tracks,
        std::vector<TLorentzVector>& tracksCM,
        KLOE::pm00& obj
    );

    static void reconstructPionsInCMFrame(
        const std::vector<TLorentzVector>& tracksCM,
        const Float_t kaonMass,
        const Float_t mPiCh,
        std::vector<TVector3>& pionsCM,
        KLOE::pm00& obj
    );
};

} // namespace KLOE

#endif // KAON_MASS_RECONSTRUCTOR_H
