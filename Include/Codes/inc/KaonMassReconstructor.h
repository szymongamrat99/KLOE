#ifndef KAON_MASS_RECONSTRUCTOR_H
#define KAON_MASS_RECONSTRUCTOR_H

#include <TVector3.h>
#include <TLorentzVector.h>
#include <vector>
#include "klspm00.hpp"

namespace KLOE {

struct KaonReconstructionResult {
    std::vector<Double_t> KaonTwoBody;     // 9 elements: px,py,pz,E,|p|,m,vx,vy,vz
    std::vector<Double_t> track1TwoBody;   // 4 elements: px,py,pz,E
    std::vector<Double_t> track2TwoBody;   // 4 elements: px,py,pz,E
};

class KaonMassReconstructor {
public:
    static KaonReconstructionResult reconstructKaonMass(
        const std::vector<Double_t>& kaonBoost,    // 9 elements: px,py,pz,E,vx,vy,vz
        const std::vector<Double_t>& interactionPoint, // 3 elements: x,y,z
        const std::vector<std::vector<Double_t>>& tracks, // 2 tracks, each 4 elements: px,py,pz,E
        const Double_t PhysicsConstants::mPiCh,  // charged pion mass
        KLOE::pm00& obj       // KLOE object for utility functions
    );

private:
    static std::vector<Double_t> calculateKaonMomentum(
        const std::vector<Double_t>& kaonBoost,
        const std::vector<Double_t>& interactionPoint
    );

    static void transformToKaonCMFrame(
        const TVector3& kaonMomLAB,
        const std::vector<std::vector<Double_t>>& tracks,
        std::vector<TLorentzVector>& tracksCM,
        KLOE::pm00& obj
    );

    static void reconstructPionsInCMFrame(
        const std::vector<TLorentzVector>& tracksCM,
        const Double_t kaonMass,
        const Double_t PhysicsConstants::mPiCh,
        std::vector<TVector3>& pionsCM,
        KLOE::pm00& obj
    );
};

} // namespace KLOE

#endif // KAON_MASS_RECONSTRUCTOR_H
