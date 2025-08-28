#include "KaonMassReconstructor.h"
#include <cmath>
#include <const.h>

namespace KLOE
{

    KaonReconstructionResult KaonMassReconstructor::reconstructKaonMass(
        const std::vector<Float_t> &kaonBoost,
        const std::vector<Float_t> &interactionPoint,
        const std::vector<std::vector<Float_t> *> &tracks,
        const Float_t mPiCh,
        KLOE::pm00 &obj)
    {
        KaonReconstructionResult result;
        result.KaonTwoBody.resize(9);
        result.track1TwoBody.resize(4);
        result.track2TwoBody.resize(4);

        // 1. Calculate kaon flight direction and momentum
        Float_t KLpath = sqrt(pow(kaonBoost[6] - interactionPoint[0], 2) +
                              pow(kaonBoost[7] - interactionPoint[1], 2) +
                              pow(kaonBoost[8] - interactionPoint[2], 2));

        TVector3 KLflightDirection(
            (kaonBoost[6] - interactionPoint[0]) / KLpath,
            (kaonBoost[7] - interactionPoint[1]) / KLpath,
            (kaonBoost[8] - interactionPoint[2]) / KLpath);

        Float_t KLmomMag = sqrt(pow(kaonBoost[0], 2) +
                                pow(kaonBoost[1], 2) +
                                pow(kaonBoost[2], 2));

        // 2. Fill KaonTwoBody vector (reconstructed kaon from 2-body decay)
        for (Int_t j = 0; j < 3; j++)
        {
            result.KaonTwoBody[j] = KLflightDirection[j] * KLmomMag;
        }
        result.KaonTwoBody[3] = kaonBoost[3]; // Energy
        result.KaonTwoBody[4] = sqrt(pow(result.KaonTwoBody[0], 2) +
                                     pow(result.KaonTwoBody[1], 2) +
                                     pow(result.KaonTwoBody[2], 2)); // Momentum magnitude
        result.KaonTwoBody[5] = sqrt(pow(result.KaonTwoBody[3], 2) -
                                     pow(result.KaonTwoBody[4], 2)); // Invariant mass
        // Copy vertex position
        result.KaonTwoBody[6] = kaonBoost[6];
        result.KaonTwoBody[7] = kaonBoost[7];
        result.KaonTwoBody[8] = kaonBoost[8];

        // 3. Go to KL CM frame
        TVector3
            z_axis(0.0, 0.0, 1.0),
            x_axis(1.0, 0.0, 0.0),
            kaonMomLAB(-result.KaonTwoBody[0] / result.KaonTwoBody[3],
                       -result.KaonTwoBody[1] / result.KaonTwoBody[3],
                       -result.KaonTwoBody[2] / result.KaonTwoBody[3]);

        TVector3 Pi1MomLAB(tracks[0]->at(0), tracks[0]->at(1), tracks[0]->at(2));
        TVector3 Pi2MomLAB(tracks[1]->at(0), tracks[1]->at(1), tracks[1]->at(2));

        TVector3 nVec = Pi1MomLAB.Cross(Pi2MomLAB);
        TVector3 cross = nVec.Cross(z_axis);

        Float_t rotAngle = nVec.Angle(z_axis);
        kaonMomLAB.Rotate(rotAngle, cross);

        Float_t rotAngleKaonX = kaonMomLAB.Angle(x_axis);
        kaonMomLAB.Rotate(rotAngleKaonX, z_axis);

        // Transform tracks to CM frame
        std::vector<TLorentzVector> tracksCM(2);
        transformToKaonCMFrame(kaonMomLAB, tracks, tracksCM, obj);

        // 4. Reconstruct pions in CM frame
        std::vector<TVector3> pionsCM(2);
        reconstructPionsInCMFrame(tracksCM, result.KaonTwoBody[5], mPiCh, pionsCM, obj);

        // 5. Transform pions back to lab frame
        kaonMomLAB = -kaonMomLAB; // Invert boost direction

        std::vector<TLorentzVector> pionsLAB(2);
        for (int i = 0; i < 2; i++)
        {
            TLorentzVector pionCM;
            Float_t E = sqrt(pionsCM[i].Mag2() + pow(mPiCh, 2));
            pionCM.SetPxPyPzE(pionsCM[i].X(), pionsCM[i].Y(), pionsCM[i].Z(), E);

            TLorentzVector pionLAB;
            obj.lorentz_transf(kaonMomLAB, pionCM, pionLAB);

            TVector3 pionMomLAB(pionLAB.Px(), pionLAB.Py(), pionLAB.Pz());
            pionMomLAB.Rotate(-rotAngleKaonX, z_axis);
            pionMomLAB.Rotate(-rotAngle, cross);

            if (i == 0)
            {
                result.track1TwoBody[0] = pionMomLAB.X();
                result.track1TwoBody[1] = pionMomLAB.Y();
                result.track1TwoBody[2] = pionMomLAB.Z();
                result.track1TwoBody[3] = sqrt(pionMomLAB.Mag2() + pow(mPiCh, 2));
            }
            else
            {
                result.track2TwoBody[0] = pionMomLAB.X();
                result.track2TwoBody[1] = pionMomLAB.Y();
                result.track2TwoBody[2] = pionMomLAB.Z();
                result.track2TwoBody[3] = sqrt(pionMomLAB.Mag2() + pow(mPiCh, 2));
            }
        }

        std::cout << tracks[0]->at(0) << " " << tracks[0]->at(1) << " " << tracks[0]->at(2) << " " << tracks[0]->at(3) << std::endl;
        std::cout << tracks[1]->at(0) << " " << tracks[1]->at(1) << " " << tracks[1]->at(2) << " " << tracks[1]->at(3) << std::endl;

        return result;
    }

    std::vector<Float_t> KaonMassReconstructor::calculateKaonMomentum(
        const std::vector<Float_t> &kaonBoost,
        const std::vector<Float_t> &interactionPoint)
    {
        std::vector<Float_t> result(9);

        // Calculate kaon flight direction
        Float_t KLpath = sqrt(pow(kaonBoost[6] - interactionPoint[0], 2) +
                              pow(kaonBoost[7] - interactionPoint[1], 2) +
                              pow(kaonBoost[8] - interactionPoint[2], 2));

        TVector3 flightDirection(
            (kaonBoost[6] - interactionPoint[0]) / KLpath,
            (kaonBoost[7] - interactionPoint[1]) / KLpath,
            (kaonBoost[8] - interactionPoint[2]) / KLpath);

        Float_t momMag = sqrt(pow(kaonBoost[0], 2) +
                              pow(kaonBoost[1], 2) +
                              pow(kaonBoost[2], 2));

        // Fill momentum vector
        for (Int_t j = 0; j < 3; j++)
        {
            result[j] = flightDirection[j] * momMag;
        }
        result[3] = kaonBoost[3];                                // Energy
        result[4] = momMag;                                      // Momentum magnitude
        result[5] = sqrt(pow(kaonBoost[3], 2) - pow(momMag, 2)); // Mass
        result[6] = kaonBoost[6];                                // Vertex x
        result[7] = kaonBoost[7];                                // Vertex y
        result[8] = kaonBoost[8];                                // Vertex z

        return result;
    }

    void KaonMassReconstructor::transformToKaonCMFrame(
        TVector3 &kaonMomLAB,
        const std::vector<std::vector<Float_t> *> &tracks,
        std::vector<TLorentzVector> &tracksCM,
        KLOE::pm00 &obj)
    {
        TVector3
            z_axis(0.0, 0.0, 1.0),
            x_axis(1.0, 0.0, 0.0);

        TVector3 trkKLMomVecLAB[2];
        TLorentzVector trkKL4VecLAB[2];

        // Get rotation angles from the kaon momentum setup
        TVector3 Pi1MomLAB(tracks[0]->at(0), tracks[0]->at(1), tracks[0]->at(2));
        TVector3 Pi2MomLAB(tracks[1]->at(0), tracks[1]->at(1), tracks[1]->at(2));

        TVector3 nVec = Pi1MomLAB.Cross(Pi2MomLAB);
        TVector3 cross = nVec.Cross(z_axis);

        Double_t rotAngle = nVec.Angle(z_axis);

        kaonMomLAB.Rotate(rotAngle, cross);

        Double_t rotAngleKaonX = kaonMomLAB.Angle(x_axis);

        // Rotate track momenta to align with kaon direction
        for (Int_t j = 0; j < 2; j++)
        {
            for (Int_t k = 0; k < 3; k++)
            {
                trkKLMomVecLAB[j][k] = tracks[j]->at(k);
            }

            trkKLMomVecLAB[j].Rotate(rotAngle, cross);
            trkKLMomVecLAB[j].Rotate(rotAngleKaonX, z_axis);

            trkKL4VecLAB[j].SetPxPyPzE(trkKLMomVecLAB[j][0],
                                       trkKLMomVecLAB[j][1],
                                       trkKLMomVecLAB[j][2],
                                       tracks[j]->at(3));

            // Transform to kaon CM frame
            obj.lorentz_transf(kaonMomLAB, trkKL4VecLAB[j], tracksCM[j]);
        }
    }

    void KaonMassReconstructor::reconstructPionsInCMFrame(
        const std::vector<TLorentzVector> &tracksCM,
        const Float_t kaonMass,
        const Float_t mPiCh,
        std::vector<TVector3> &pionsCM,
        KLOE::pm00 &obj)
    {
        TVector3 trkKLMomVecKaonCM[2];

        // Extract 3-momentum from 4-momentum in CM frame
        trkKLMomVecKaonCM[0].SetXYZ(tracksCM[0].Px(), tracksCM[0].Py(), tracksCM[0].Pz());
        trkKLMomVecKaonCM[1].SetXYZ(tracksCM[1].Px(), tracksCM[1].Py(), tracksCM[1].Pz());

        // Calculate momentum magnitude for 2-body decay
        Float_t PiMomMagKaonCM = obj.TwoBodyDecayMass(kaonMass, mPiCh, mPiCh);

        // Get angles from the track directions in CM frame
        Float_t angle1 = trkKLMomVecKaonCM[0].Phi();
        Float_t angle2 = trkKLMomVecKaonCM[1].Phi();
        Float_t theta1 = trkKLMomVecKaonCM[0].Theta();
        Float_t theta2 = trkKLMomVecKaonCM[1].Theta();

        // Apply gamma correction (from original algorithm)
        Float_t gamma = (M_PI_2 - 0.5 * abs(angle1) - 0.5 * abs(angle2));

        // Reconstruct pion momenta in CM frame with correct magnitude
        if (angle1 < 0.0)
        {
            pionsCM[0].SetXYZ(PiMomMagKaonCM * sin(theta1) * cos(angle1 - gamma),
                              PiMomMagKaonCM * sin(theta1) * sin(angle1 - gamma),
                              PiMomMagKaonCM * cos(theta1));
        }
        else
        {
            pionsCM[0].SetXYZ(PiMomMagKaonCM * sin(theta1) * cos(angle1 + gamma),
                              PiMomMagKaonCM * sin(theta1) * sin(angle1 + gamma),
                              PiMomMagKaonCM * cos(theta1));
        }

        if (angle2 < 0.0)
        {
            pionsCM[1].SetXYZ(PiMomMagKaonCM * sin(theta2) * cos(angle2 - gamma),
                              PiMomMagKaonCM * sin(theta2) * sin(angle2 - gamma),
                              PiMomMagKaonCM * cos(theta2));
        }
        else
        {
            pionsCM[1].SetXYZ(PiMomMagKaonCM * sin(theta2) * cos(angle2 + gamma),
                              PiMomMagKaonCM * sin(theta2) * sin(angle2 + gamma),
                              PiMomMagKaonCM * cos(theta2));
        }
    }

} // namespace KLOE
