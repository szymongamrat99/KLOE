#include <pi0_photon_pair.h>

// PhotonMom: N photons with 4 mom components, unsorted
// PhotonMomPi0: photons paired to pi0s with 4 mom components

Int_t Pi0PhotonPair(Int_t ClusterIndex[4], Float_t PhotonMom[4][8], Int_t ClusterIndexPi0[2][2], Float_t PhotonMomPi0[2][2][4], Float_t Pi0Mom[2][4], Bool_t OmegaFlag, Float_t PichFourMom[2][4], Float_t OmegaMom[4])
{
  const Int_t
      nums = 6;
  Int_t
      UniqueSets[3][4] = {
          {0, 1, 2, 3},
          {0, 2, 1, 3},
          {0, 3, 1, 2}};
  Double_t
      MinPseudoChi2 = 9E9;

  std::vector<std::vector<Int_t>>
      TmpClusterIndex(nums, std::vector<Int_t>(4));

  std::vector<Double_t>
      PseudoChi2(nums),
      InvMassPi01(nums),
      InvMassPi02(nums),
      InvMassOmega(nums);

  std::vector<std::vector<Double_t>>
      Pi0FourMom1(nums, std::vector<Double_t>(4)),
      Pi0FourMom2(nums, std::vector<Double_t>(4)),
      OmegaFourMom(nums, std::vector<Double_t>(4));

  for (Int_t i = 0; i < 3; i++)
  {
    for (Int_t k = 0; k < 4; k++)
    {
      TmpClusterIndex[i * 2][k] = ClusterIndex[UniqueSets[i][k]];
      TmpClusterIndex[i * 2 + 1][k] = ClusterIndex[UniqueSets[i][k]];

      Pi0FourMom1[i * 2][k] = PhotonMom[UniqueSets[i][0]][k] + PhotonMom[UniqueSets[i][1]][k];
      Pi0FourMom1[i * 2 + 1][k] = PhotonMom[UniqueSets[i][0]][k] + PhotonMom[UniqueSets[i][1]][k];

      Pi0FourMom2[i * 2][k] = PhotonMom[UniqueSets[i][2]][k] + PhotonMom[UniqueSets[i][3]][k];
      Pi0FourMom2[i * 2 + 1][k] = PhotonMom[UniqueSets[i][2]][k] + PhotonMom[UniqueSets[i][3]][k];

      if (OmegaFlag == true)
      {
        OmegaFourMom[i * 2][k] = Pi0FourMom1[i * 2][k] + PichFourMom[0][k] + PichFourMom[1][k];
        OmegaFourMom[i * 2 + 1][k] = Pi0FourMom2[i * 2 + 1][k] + PichFourMom[0][k] + PichFourMom[1][k];
      }
    }

    InvMassPi01[i * 2] = sqrt(pow(Pi0FourMom1[i * 2][3], 2) -
                              pow(Pi0FourMom1[i * 2][0], 2) -
                              pow(Pi0FourMom1[i * 2][1], 2) -
                              pow(Pi0FourMom1[i * 2][2], 2));
    InvMassPi01[i * 2 + 1] = sqrt(pow(Pi0FourMom1[i * 2 + 1][3], 2) -
                                  pow(Pi0FourMom1[i * 2 + 1][0], 2) -
                                  pow(Pi0FourMom1[i * 2 + 1][1], 2) -
                                  pow(Pi0FourMom1[i * 2 + 1][2], 2));

    InvMassPi02[i * 2] = sqrt(pow(Pi0FourMom2[i * 2][3], 2) -
                              pow(Pi0FourMom2[i * 2][0], 2) -
                              pow(Pi0FourMom2[i * 2][1], 2) -
                              pow(Pi0FourMom2[i * 2][2], 2));
    InvMassPi02[i * 2 + 1] = sqrt(pow(Pi0FourMom2[i * 2 + 1][3], 2) -
                                  pow(Pi0FourMom2[i * 2 + 1][0], 2) -
                                  pow(Pi0FourMom2[i * 2 + 1][1], 2) -
                                  pow(Pi0FourMom2[i * 2 + 1][2], 2));

    if (OmegaFlag == true)
    {
      InvMassOmega[i * 2] = sqrt(pow(OmegaFourMom[i * 2][3], 2) -
                                 pow(OmegaFourMom[i * 2][0], 2) -
                                 pow(OmegaFourMom[i * 2][1], 2) -
                                 pow(OmegaFourMom[i * 2][2], 2));
      InvMassOmega[i * 2 + 1] = sqrt(pow(OmegaFourMom[i * 2 + 1][3], 2) -
                                     pow(OmegaFourMom[i * 2 + 1][0], 2) -
                                     pow(OmegaFourMom[i * 2 + 1][1], 2) -
                                     pow(OmegaFourMom[i * 2 + 1][2], 2));
    }
    else
    {
      InvMassOmega[i * 2] = PhysicsConstants::mOmega;
      InvMassOmega[i * 2 + 1] = PhysicsConstants::mOmega;
    }
  }

  for (Int_t i = 0; i < 3; i++)
  {
    PseudoChi2[i * 2] = pow((InvMassPi01[i * 2] - PhysicsConstants::mPi0) / 17.0, 2) + pow((InvMassPi02[i * 2] - PhysicsConstants::mPi0) / 17.0, 2) + pow((InvMassOmega[i * 2] - PhysicsConstants::mOmega) / 20.0, 2);
    PseudoChi2[i * 2 + 1] = pow((InvMassPi01[i * 2 + 1] - PhysicsConstants::mPi0) / 17.0, 2) + pow((InvMassPi02[i * 2 + 1] - PhysicsConstants::mPi0) / 17.0, 2) + pow((InvMassOmega[i * 2 + 1] - PhysicsConstants::mOmega) / 20.0, 2);
  }

  Int_t min_iter = std::distance(PseudoChi2.begin(), std::min_element(PseudoChi2.begin(), PseudoChi2.end()));

  if (PseudoChi2[min_iter] < MinPseudoChi2)
  {
    Float_t
        test = (min_iter / 2.) - floor(min_iter / 2.);
    Int_t i = floor(min_iter / 2.);

    for (Int_t k = 0; k < 4; k++)
    {
      PhotonMomPi0[0][0][k] = PhotonMom[UniqueSets[i][0]][k];
      PhotonMomPi0[0][1][k] = PhotonMom[UniqueSets[i][1]][k];
      PhotonMomPi0[1][0][k] = PhotonMom[UniqueSets[i][2]][k];
      PhotonMomPi0[1][1][k] = PhotonMom[UniqueSets[i][3]][k];
    }

    if (1)
    {
      if (min_iter % 2)
      {
        ClusterIndexPi0[0][0] = TmpClusterIndex[min_iter][2];
        ClusterIndexPi0[0][1] = TmpClusterIndex[min_iter][3];

        ClusterIndexPi0[1][0] = TmpClusterIndex[min_iter][0];
        ClusterIndexPi0[1][1] = TmpClusterIndex[min_iter][1];

        for (Int_t k = 0; k < 4; k++)
        {
          Pi0Mom[0][k] = Pi0FourMom2[min_iter][k];
          Pi0Mom[1][k] = Pi0FourMom1[min_iter][k];

          OmegaMom[k] = Pi0Mom[1][k] + PichFourMom[0][k] + PichFourMom[1][k];
        }
      }
      else
      {
        ClusterIndexPi0[0][0] = TmpClusterIndex[min_iter][0];
        ClusterIndexPi0[0][1] = TmpClusterIndex[min_iter][1];

        ClusterIndexPi0[1][0] = TmpClusterIndex[min_iter][2];
        ClusterIndexPi0[1][1] = TmpClusterIndex[min_iter][3];

        for (Int_t k = 0; k < 4; k++)
        {
          Pi0Mom[0][k] = Pi0FourMom1[min_iter][k];
          Pi0Mom[1][k] = Pi0FourMom2[min_iter][k];

          OmegaMom[k] = Pi0Mom[0][k] + PichFourMom[0][k] + PichFourMom[1][k];
        }
      }
    }
  }

  return 0;
}