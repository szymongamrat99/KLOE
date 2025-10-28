#include <ConstraintsTrilateration.h>

namespace KLOE
{
  void ConstraintsTrilateration::ResetParameters()
  {
    // Setting parameters from the array p to the member variables
    // of the base classes to be used in the constraint calculations.
    if (bhabha_mom.size() != 4)
    {
      bhabha_mom.clear();
      bhabha_mom.resize(4);
    }

    if (cluster.size() != 4)
    {
      for (Int_t i = 0; i < cluster.size(); i++)
        cluster[i].clear();
      cluster.clear();
      cluster.resize(4);
    }

    for (Int_t i = 0; i < cluster.size(); i++)
    {
      if (cluster[i].size() != 5)
      {
        cluster[i].clear();
        cluster[i].resize(5);
      }
    }

    if (pionCh.size() != 2)
    {
      pionCh.clear();
      pionCh.resize(2);
    }

    if (photon.size() != 4)
    {
      photon.clear();
      photon.resize(4);
    }

    if (ip.size() != 3)
    {
      ip.clear();
      ip.resize(3);
    }
  }

  void ConstraintsTrilateration::SetParameters(Double_t *p)
  {
    // Setting parameters from the array p to the member variables
    // of the base classes to be used in the constraint calculations.
    ResetParameters();

    for (Int_t i = 0; i < 4; i++)
    {
      for (Int_t j = 0; j < 5; j++)
      {
        fphoton[i].clusterParams[j] = p[i * 5 + j];
      }
    }

    for (Int_t i = 0; i < 4; i++)
    {
      fphi.fourMom[i] = p[20 + i];

      if (i < 3)
      {
        fphi.vtxPos[i] = p[24 + i];
      };
    }

    fphi.SetTotalVector();
  }

  void ConstraintsTrilateration::IntermediateReconstruction()
  {
    static kaonNeutral KnerecTmp[2];        // Temporary kaon neutral objects
    static neutralParticle photonTmp[2][4]; // Temporary photon objects

    static Float_t ipTmp[2][3]; // Temporary interaction point

    std::array<Double_t, 2> value = {0., 0.};

    for (Int_t i = 0; i < 2; i++)
      KnerecTmp[i] = Knerec;

    for (Int_t i = 0; i < 2; i++)
    {
      for (Int_t j = 0; j < 4; j++)
        photonTmp[i][j] = photon[j];
    }

    Reconstructor R; // Reconstructor object
    Solution S;      // Solution struct

    // Setting clusters for a solution
    for (Int_t k = 0; k < 4; k++)
    {
      R.SetClu(k, photon[k].clusterParams[0],
               photon[k].clusterParams[1],
               photon[k].clusterParams[2],
               photon[k].clusterParams[3],
               photon[k].clusterParams[4]);
    }
    // --------------------------------------------

    S = R.MySolve(_selected); // Filling up the structure

    for (Int_t i = 0; i < 2; i++)
    {
      if (!S.error[i])
      {
        KnerecTmp[i].fourPos[0] = S.sol[i][0];
        KnerecTmp[i].fourPos[1] = S.sol[i][1];
        KnerecTmp[i].fourPos[2] = S.sol[i][2];
        KnerecTmp[i].fourPos[3] = S.sol[i][3];

        for (Int_t j = 0; j < 4; j++)
        {
          neutral_mom(photonTmp[i][j].clusterParams[0],
                      photonTmp[i][j].clusterParams[1],
                      photonTmp[i][j].clusterParams[2],
                      photonTmp[i][j].clusterParams[4],
                      KnerecTmp[i].fourPos.data(),
                      photonTmp[i][j].fourMom.data());

          photonTmp[i][j].fourPos[0] = photonTmp[i][j].clusterParams[0];
          photonTmp[i][j].fourPos[1] = photonTmp[i][j].clusterParams[1];
          photonTmp[i][j].fourPos[2] = photonTmp[i][j].clusterParams[2];
          photonTmp[i][j].fourPos[3] = photonTmp[i][j].clusterParams[3];

          photonTmp[i][j].calculatePath(KnerecTmp[i].fourPos.data());
          photonTmp[i][j].calculateTimeOfFlightPhoton();
          photonTmp[i][j].SetTotalVectorPhoton();
        }

        for (Int_t j = 0; j < 4; j++)
          KnerecTmp[i].fourMom[j] = photonTmp[i][0].fourMom[j] +
                                    photonTmp[i][1].fourMom[j] +
                                    photonTmp[i][2].fourMom[j] +
                                    photonTmp[i][3].fourMom[j];

        Float_t X_line[3] = {KnerecTmp[i].fourPos[0],
                             KnerecTmp[i].fourPos[1],
                             KnerecTmp[i].fourPos[2]},
                mom[3] = {KnerecTmp[i].fourMom[0],
                          KnerecTmp[i].fourMom[1],
                          KnerecTmp[i].fourMom[2]},
                xB[3] = {phi.vtxPos[0],
                         phi.vtxPos[1],
                         phi.vtxPos[2]},
                plane_perp[3] = {0.,
                                 phi.fourMom[1],
                                 0.};

        // Corrected IP event by event
        IPBoostCorr(X_line, mom, xB, plane_perp, ipTmp[i]);

        ipTmp[i][0] = phi.vtxPos[0];
        ipTmp[i][1] = phi.vtxPos[1];
        // // ip[2] is fitted
        if (abs(ipTmp[i][2] - phi.vtxPos[2]) > 2.)
          ipTmp[i][2] = phi.vtxPos[2];

        KnerecTmp[i].calculatePath(ipTmp[i]);
        KnerecTmp[i].SetTotalVector();
        KnerecTmp[i].calculateLifetimeLAB();
        KnerecTmp[i].fourPos[3] = S.sol[i][3];
        KnerecTmp[i].total[9] = S.sol[i][3];

        value[i] = sqrt(pow(KnerecTmp[i].total[5] - PhysicsConstants::mK0, 2) +
                        pow(KnerecTmp[i].fourPos[3] - KnerecTmp[i].lifetimeLAB, 2));
      }
      else
      {
        value[i] = 999999.;
      }
    }

    if (value[0] < value[1])
    {
      Knerec = KnerecTmp[0];

      for (Int_t j = 0; j < 4; j++)
        photon[j] = photonTmp[0][j];

      for (Int_t j = 0; j < 3; j++)
        ip[j] = ipTmp[0][j];
    }
    else if (value[1] < value[0])
    {
      Knerec = KnerecTmp[1];

      for (Int_t j = 0; j < 4; j++)
        photon[j] = photonTmp[1][j];

      for (Int_t j = 0; j < 3; j++)
        ip[j] = ipTmp[1][j];
    }
    else // If strange values just put first solution in
    {
      Knerec = KnerecTmp[0];
      for (Int_t j = 0; j < 4; j++)
        photon[j] = photonTmp[0][j];

      for (Int_t j = 0; j < 3; j++)
        ip[j] = ipTmp[0][j];
    }

    Float_t boostVec[3] = {-phi.fourMom[0] / phi.fourMom[3],
                           -phi.fourMom[1] / phi.fourMom[3],
                           -phi.fourMom[2] / phi.fourMom[3]};

    lorentz_transf(boostVec,
                   Knerec.fourMom.data(),
                   KnerecCMPhi.fourMom.data());
  }

  void ConstraintsTrilateration::IntermediateReconstruction(Double_t *p)
  {
    if (fip.size() != 3)
    {
      fip.clear();
      fip.resize(3);
    }


    for (Int_t i = 0; i < 4; i++)
    {
      for (Int_t j = 0; j < 5; j++)
      {
        fphoton[i].clusterParams[j] = p[i * 5 + j];
      }
    }

    for (Int_t i = 0; i < 4; i++)
    {
      fphi.fourMom[i] = p[20 + i];

      if (i < 3)
      {
        fphi.vtxPos[i] = p[24 + i];
      };
    }

    fphi.SetTotalVector();

    kaonNeutral KnerecTmp[2];        // Temporary kaon neutral objects
    neutralParticle photonTmp[2][4]; // Temporary photon objects

    Float_t ipTmp[2][3]; // Temporary interaction point

    std::array<Double_t, 2> value = {0., 0.};

    for (Int_t i = 0; i < 2; i++)
    {
      for (Int_t j = 0; j < 4; j++)
        photonTmp[i][j] = fphoton[j];
    }

    Reconstructor R; // Reconstructor object
    Solution S;      // Solution struct

    // Setting clusters for a solution
    for (Int_t k = 0; k < 4; k++)
    {
      R.SetClu(k, fphoton[k].clusterParams[0],
               fphoton[k].clusterParams[1],
               fphoton[k].clusterParams[2],
               fphoton[k].clusterParams[3],
               fphoton[k].clusterParams[4]);
    }
    // --------------------------------------------

    S = R.MySolve(_selected); // Filling up the structure

    for (Int_t i = 0; i < 2; i++)
    {
      if (!S.error[i])
      {
        KnerecTmp[i].fourPos[0] = S.sol[i][0];
        KnerecTmp[i].fourPos[1] = S.sol[i][1];
        KnerecTmp[i].fourPos[2] = S.sol[i][2];
        KnerecTmp[i].fourPos[3] = S.sol[i][3];

        for (Int_t j = 0; j < 4; j++)
        {
          neutral_mom(photonTmp[i][j].clusterParams[0],
                      photonTmp[i][j].clusterParams[1],
                      photonTmp[i][j].clusterParams[2],
                      photonTmp[i][j].clusterParams[4],
                      KnerecTmp[i].fourPos.data(),
                      photonTmp[i][j].fourMom.data());

          photonTmp[i][j].fourPos[0] = photonTmp[i][j].clusterParams[0];
          photonTmp[i][j].fourPos[1] = photonTmp[i][j].clusterParams[1];
          photonTmp[i][j].fourPos[2] = photonTmp[i][j].clusterParams[2];
          photonTmp[i][j].fourPos[3] = photonTmp[i][j].clusterParams[3];

          photonTmp[i][j].calculatePath(KnerecTmp[i].fourPos.data());
          photonTmp[i][j].calculateTimeOfFlightPhoton();
          photonTmp[i][j].SetTotalVectorPhoton();
        }

        for (Int_t j = 0; j < 4; j++)
          KnerecTmp[i].fourMom[j] = photonTmp[i][0].fourMom[j] +
                                    photonTmp[i][1].fourMom[j] +
                                    photonTmp[i][2].fourMom[j] +
                                    photonTmp[i][3].fourMom[j];

        Float_t X_line[3] = {KnerecTmp[i].fourPos[0],
                             KnerecTmp[i].fourPos[1],
                             KnerecTmp[i].fourPos[2]},
                mom[3] = {KnerecTmp[i].fourMom[0],
                          KnerecTmp[i].fourMom[1],
                          KnerecTmp[i].fourMom[2]},
                xB[3] = {fphi.vtxPos[0],
                         fphi.vtxPos[1],
                         fphi.vtxPos[2]},
                plane_perp[3] = {0.,
                                 phi.fourMom[1],
                                 0.};

        // Corrected IP event by event
        IPBoostCorr(X_line, mom, xB, plane_perp, ipTmp[i]);

        ipTmp[i][0] = fphi.vtxPos[0];
        ipTmp[i][1] = fphi.vtxPos[1];
        // // ip[2] is fitted
        if (abs(ipTmp[i][2] - fphi.vtxPos[2]) > 2.)
          ipTmp[i][2] = fphi.vtxPos[2];

        KnerecTmp[i].calculatePath(ipTmp[i]);
        KnerecTmp[i].SetTotalVector();
        KnerecTmp[i].calculateLifetimeLAB();
        KnerecTmp[i].fourPos[3] = S.sol[i][3];
        KnerecTmp[i].total[9] = S.sol[i][3];

        value[i] = sqrt(pow(KnerecTmp[i].total[5] - PhysicsConstants::mK0, 2) +
                        pow(KnerecTmp[i].fourPos[3] - KnerecTmp[i].lifetimeLAB, 2));
      }
      else
      {
        value[i] = 999999.;
      }
    }

    if (value[0] < value[1])
    {
      fKnerec = KnerecTmp[0];

      for (Int_t j = 0; j < 4; j++)
        fphoton[j] = photonTmp[0][j];

      for (Int_t j = 0; j < 3; j++)
        fip[j] = ipTmp[0][j];
    }
    else if (value[1] < value[0])
    {
      fKnerec = KnerecTmp[1];

      for (Int_t j = 0; j < 4; j++)
        fphoton[j] = photonTmp[1][j];

      for (Int_t j = 0; j < 3; j++)
        fip[j] = ipTmp[1][j];
    }
    else // If strange values just put first solution in
    {
      fKnerec = KnerecTmp[0];
      for (Int_t j = 0; j < 4; j++)
        fphoton[j] = photonTmp[0][j];

      for (Int_t j = 0; j < 3; j++)
        fip[j] = ipTmp[0][j];
    }

    Float_t boostVec[3] = {-fphi.fourMom[0] / fphi.fourMom[3],
                           -fphi.fourMom[1] / fphi.fourMom[3],
                           -fphi.fourMom[2] / fphi.fourMom[3]};

    lorentz_transf(boostVec,
                   fKnerec.fourMom.data(),
                   fKnerecCMPhi.fourMom.data());
  }

  Double_t ConstraintsTrilateration::EnergyConsvCM(Double_t *x, Double_t *p)
  {
    IntermediateReconstruction(p);

    return fKnerecCMPhi.fourMom[3] - (fphi.fourMom[3] / 2.);
  }

  Double_t ConstraintsTrilateration::MinvConsv(Double_t *x, Double_t *p)
  {
    IntermediateReconstruction(p);

    return fKnerec.total[5] - PhysicsConstants::mK0;
  }

  Double_t ConstraintsTrilateration::NeutralPathConsvLAB(Double_t *x, Double_t *p)
  {
    IntermediateReconstruction(p);

    Double_t pathComp = abs(fKnerec.fourPos[_chosenComponent] -
                            fip[_chosenComponent]),
             kneVelComp = abs(fKnerec.fourMom[_chosenComponent] / fKnerec.fourMom[3] * PhysicsConstants::cVel);

    return fKnerec.fourPos[3] - pathComp / kneVelComp;
  }
}
