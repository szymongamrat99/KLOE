#include <ConstraintsOmega.h>

using namespace KLOE;

void ConstraintsOmega::ResetParameters()
{
  // Setting parameters from the array p to the member variables
  // of the base classes to be used in the constraint calculations.

  if (bhabha_mom.size() != 4)
  {
    bhabha_mom.clear();
    bhabha_mom.resize(4);
  }

  if (pionCh.size() != 2)
  {
    pionCh.clear();
    pionCh.resize(2);
  }

  if (pionNe.size() != 2)
  {
    pionNe.clear();
    pionNe.resize(2);
  }

  if (photon.size() != 4)
  {
    photon.clear();
    photon.resize(4);
  }
}

void ConstraintsOmega::SetParameters(Double_t *p)
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

  for (Int_t i = 0; i < 2; i++)
  {
    fpionCh[i].fourMom[0] = p[20 + i * 3];
    fpionCh[i].fourMom[1] = p[20 + i * 3 + 1];
    fpionCh[i].fourMom[2] = p[20 + i * 3 + 2];
    fpionCh[i].fourMom[3] = sqrt(pow(fpionCh[i].fourMom[0], 2) +
                                 pow(fpionCh[i].fourMom[1], 2) +
                                 pow(fpionCh[i].fourMom[2], 2) +
                                 pow(PhysicsConstants::mPiCh, 2));

    pionCh[i].fourMomFilled = true;
  }

  for (Int_t i = 0; i < 4; i++)
  {
    fphi.fourMom[i] = p[26 + i];

    if (i < 3)
    {
      fomega.fourPos[i] = p[30 + i];
      fphi.vtxPos[i] = p[33 + i];
    };
  }
}

void ConstraintsOmega::IntermediateReconstruction()
{
  // Intermediate reconstruction to be done after setting the parameters
  // and before calculating the constraints.

  if (fip.size() != 3)
  {
    fip.clear();
    fip.resize(3);
  }

  phi.SetTotalVector();

  for (Int_t i = 0; i < 4; i++)
  {
    neutral_mom(photon[i].clusterParams[0],
                photon[i].clusterParams[1],
                photon[i].clusterParams[2],
                photon[i].clusterParams[4],
                omega.fourPos.data(),
                photon[i].fourMom.data());

    photon[i].fourPos[0] = photon[i].clusterParams[0];
    photon[i].fourPos[1] = photon[i].clusterParams[1];
    photon[i].fourPos[2] = photon[i].clusterParams[2];
    photon[i].fourPos[3] = photon[i].clusterParams[3];

    photon[i].calculatePath(omega.fourPos.data());
    photon[i].calculateTimeOfFlightPhoton();
    photon[i].SetTotalVectorPhoton();

    photon[i].fourMomFilled = true;
  }

  std::vector<Int_t> bestPairingIndexNeutral, bestPairingIndexCharged;

  neutRec.PhotonPairingToPi0WithOmega(photon, pionCh, bestPairingIndexNeutral, bestPairingIndexCharged, omega, pionNe);

  neutRec.Pi0Reconstruction(pionNe);

  for (Int_t i = 0; i < 2; i++)
  {
    pionNe[i].SetTotalVector();
  }

  omega.total[6] = phi.vtxPos[0];
  omega.total[7] = phi.vtxPos[1];
  omega.total[8] = phi.vtxPos[2];
}

void ConstraintsOmega::IntermediateReconstruction(Double_t *p)
{
  // Intermediate reconstruction to be done after setting the parameters
  // and before calculating the constraints.

  neutRec.SetNumberOfPhotons(4);

  for (Int_t i = 0; i < 4; i++)
  {
    for (Int_t j = 0; j < 5; j++)
    {
      fphoton[i].clusterParams[j] = p[i * 5 + j];
    }
  }

  for (Int_t i = 0; i < 2; i++)
  {
    fpionCh[i].fourMom[0] = p[20 + i * 3];
    fpionCh[i].fourMom[1] = p[20 + i * 3 + 1];
    fpionCh[i].fourMom[2] = p[20 + i * 3 + 2];
    fpionCh[i].fourMom[3] = sqrt(pow(fpionCh[i].fourMom[0], 2) +
                                 pow(fpionCh[i].fourMom[1], 2) +
                                 pow(fpionCh[i].fourMom[2], 2) +
                                 pow(PhysicsConstants::mPiCh, 2));

    fpionCh[i].fourMomFilled = true;
  }

  for (Int_t i = 0; i < 4; i++)
  {
    fphi.fourMom[i] = p[26 + i];

    if (i < 3)
    {
      fomega.fourPos[i] = p[30 + i];
      fphi.vtxPos[i] = p[33 + i];
    };
  }

  fphi.SetTotalVector();

  for (Int_t i = 0; i < 4; i++)
  {
    neutral_mom(fphoton[i].clusterParams[0],
                fphoton[i].clusterParams[1],
                fphoton[i].clusterParams[2],
                fphoton[i].clusterParams[4],
                fomega.fourPos.data(),
                fphoton[i].fourMom.data());

    fphoton[i].fourPos[0] = fphoton[i].clusterParams[0];
    fphoton[i].fourPos[1] = fphoton[i].clusterParams[1];
    fphoton[i].fourPos[2] = fphoton[i].clusterParams[2];
    fphoton[i].fourPos[3] = fphoton[i].clusterParams[3];

    fphoton[i].calculatePath(fomega.fourPos.data());
    fphoton[i].calculateTimeOfFlightPhoton();
    fphoton[i].SetTotalVectorPhoton();

    fphoton[i].fourMomFilled = true;
  }

  std::vector<Int_t> bestPairingIndexNeutral, bestPairingIndexCharged;

  std::vector<chargedParticle> s_pionCh = {fpionCh[0], fpionCh[1]};
  std::vector<neutralParticle> s_photon = {fphoton[0], fphoton[1], fphoton[2], fphoton[3]},
                               s_pionNe(2, neutralParticle());

#pragma omp critical
{
  neutRec.PhotonPairingToPi0WithOmega(s_photon, s_pionCh, bestPairingIndexNeutral, bestPairingIndexCharged, fomega, s_pionNe);
}

  for (Int_t i = 0; i < 2; i++)
  {
    fpionNe[i] = s_pionNe[i];
  }

  fomega.total[6] = fphi.vtxPos[0];
  fomega.total[7] = fphi.vtxPos[1];
  fomega.total[8] = fphi.vtxPos[2];
}

Double_t ConstraintsOmega::FourMomConsvLAB(Double_t *x, Double_t *p)
{
  IntermediateReconstruction(p);

  // Conservation of the chosen component of the 4-momentum in the LAB system
  return (fphi.fourMom[_chosen4MomComponent] -
          fpionCh[0].fourMom[_chosen4MomComponent] -
          fpionCh[1].fourMom[_chosen4MomComponent] -
          fphoton[0].fourMom[_chosen4MomComponent] -
          fphoton[1].fourMom[_chosen4MomComponent] -
          fphoton[2].fourMom[_chosen4MomComponent] -
          fphoton[3].fourMom[_chosen4MomComponent]);
};

Double_t ConstraintsOmega::PhotonPathConsvLAB(Double_t *x, Double_t *p)
{
  IntermediateReconstruction(p);

  return fphoton[_chosenPhoton].fourPos[3] - fphoton[_chosenPhoton].timeOfFlight;
}

Double_t ConstraintsOmega::MinvConsv(Double_t *x, Double_t *p)
{
  IntermediateReconstruction(p);

  std::map<std::string, Float_t> minvModes = {
      {"omega", fomega.mass}};

  return minvModes[_minvMode] - PhysicsConstants::mOmega; // MeV/c^2
}