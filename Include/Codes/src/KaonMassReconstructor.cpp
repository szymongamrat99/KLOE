#include "KaonMassReconstructor.h"
#include <cmath>

namespace KLOE {

KaonReconstructionResult KaonMassReconstructor::reconstructKaonMass(
    const std::vector<Double_t>& kaonBoost,
    const std::vector<Double_t>& interactionPoint,
    const std::vector<std::vector<Double_t>>& tracks,
    const Double_t mPiCh,
    KLOE::pm00& obj)
{
    KaonReconstructionResult result;
    result.KaonTwoBody.resize(9);
    result.track1TwoBody.resize(4);
    result.track2TwoBody.resize(4);

    // 1. Calculate kaon flight direction and momentum
    Double_t KLpath = sqrt(pow(kaonBoost[6] - interactionPoint[0], 2) +
                          pow(kaonBoost[7] - interactionPoint[1], 2) +
                          pow(kaonBoost[8] - interactionPoint[2], 2));
    
    TVector3 flightDirection(
        (kaonBoost[6] - interactionPoint[0]) / KLpath,
        (kaonBoost[7] - interactionPoint[1]) / KLpath,
        (kaonBoost[8] - interactionPoint[2]) / KLpath
    );

    Double_t momMag = sqrt(pow(kaonBoost[0], 2) +
                          pow(kaonBoost[1], 2) +
                          pow(kaonBoost[2], 2));

    // Fill KaonTwoBody vector
    for(Int_t j = 0; j < 3; j++) {
        result.KaonTwoBody[j] = flightDirection[j] * momMag;
    }
    result.KaonTwoBody[3] = kaonBoost[3];  // Energy
    result.KaonTwoBody[4] = momMag;        // Momentum magnitude
    result.KaonTwoBody[5] = sqrt(pow(kaonBoost[3], 2) - pow(momMag, 2));  // Mass
    // Copy vertex position
    result.KaonTwoBody[6] = kaonBoost[6];
    result.KaonTwoBody[7] = kaonBoost[7];
    result.KaonTwoBody[8] = kaonBoost[8];

    // 2. Transform to kaon rest frame
    TVector3 
        z_axis(0.0, 0.0, 1.0),
        x_axis(1.0, 0.0, 0.0),
        kaonMomLAB(-result.KaonTwoBody[0] / result.KaonTwoBody[3],
                   -result.KaonTwoBody[1] / result.KaonTwoBody[3],
                   -result.KaonTwoBody[2] / result.KaonTwoBody[3]);

    TVector3 Pi1MomLAB(tracks[0][0], tracks[0][1], tracks[0][2]);
    TVector3 Pi2MomLAB(tracks[1][0], tracks[1][1], tracks[1][2]);
    
    TVector3 nVec = Pi1MomLAB.Cross(Pi2MomLAB);
    TVector3 cross = nVec.Cross(z_axis);
    
    Double_t rotAngle = nVec.Angle(z_axis);
    kaonMomLAB.Rotate(rotAngle, cross);
    
    Double_t rotAngleKaonX = kaonMomLAB.Angle(x_axis);
    kaonMomLAB.Rotate(rotAngleKaonX, z_axis);

    // Transform tracks to CM frame
    std::vector<TLorentzVector> tracksCM(2);
    transformToKaonCMFrame(kaonMomLAB, tracks, tracksCM, obj);

    // 3. Reconstruct pions in CM frame
    std::vector<TVector3> pionsCM(2);
    reconstructPionsInCMFrame(tracksCM, result.KaonTwoBody[5], mPiCh, pionsCM, obj);

    // 4. Transform back to lab frame
    kaonMomLAB = -kaonMomLAB;  // Invert boost direction
    
    std::vector<TLorentzVector> pionsLAB(2);
    for(int i = 0; i < 2; i++) {
        TLorentzVector pionCM;
        Double_t E = sqrt(pionsCM[i].Mag2() + pow(mPiCh, 2));
        pionCM.SetPxPyPzE(pionsCM[i].X(), pionsCM[i].Y(), pionsCM[i].Z(), E);
        
        TLorentzVector pionLAB;
        obj.lorentz_transf(kaonMomLAB, pionCM, pionLAB);
        
        TVector3 pionMomLAB(pionLAB.Px(), pionLAB.Py(), pionLAB.Pz());
        pionMomLAB.Rotate(-rotAngleKaonX, z_axis);
        pionMomLAB.Rotate(-rotAngle, cross);
        
        result.track1TwoBody[i] = pionMomLAB.X();
        result.track1TwoBody[i+1] = pionMomLAB.Y();
        result.track1TwoBody[i+2] = pionMomLAB.Z();
        result.track1TwoBody[3] = sqrt(pionMomLAB.Mag2() + pow(mPiCh, 2));
    }

    return result;
}

std::vector<Double_t> KaonMassReconstructor::calculateKaonMomentum(
    const std::vector<Double_t>& kaonBoost,
    const std::vector<Double_t>& interactionPoint)
{
    std::vector<Double_t> result(9);
    // Implementation...
    return result;
}

void KaonMassReconstructor::transformToKaonCMFrame(
    const TVector3& kaonMomLAB,
    const std::vector<std::vector<Double_t>>& tracks,
    std::vector<TLorentzVector>& tracksCM,
    KLOE::pm00& obj)
{
    // Implementation...
}

void KaonMassReconstructor::reconstructPionsInCMFrame(
    const std::vector<TLorentzVector>& tracksCM,
    const Double_t kaonMass,
    const Double_t mPiCh,
    std::vector<TVector3>& pionsCM,
    KLOE::pm00& obj)
{
    // Implementation...
}

} // namespace KLOE
