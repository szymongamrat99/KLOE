#include <ConstraintsOmega.h>

using namespace KLOE;

void ConstraintsOmega::ResetParameters()
{
  // Setting parameters from the array p to the member variables
  // of the base classes to be used in the constraint calculations.

  neutRec.SetNumberOfPhotons(4);

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
      photon[i].clusterParams[j] = p[i * 5 + j];
    }
  }

  for (Int_t i = 0; i < 2; i++)
  {
    pionCh[i].trackParams[0] = p[20 + i * 3];
    pionCh[i].trackParams[1] = p[20 + i * 3 + 1];
    pionCh[i].trackParams[2] = p[20 + i * 3 + 2];
  }

  for (Int_t i = 0; i < 4; i++)
  {
    phi.fourMom[i] = p[26 + i];

    if (i < 3)
    {
      omega.fourPos[i] = p[30 + i];
      phi.vtxPos[i] = p[33 + i];
    };
  }
}

void ConstraintsOmega::IntermediateReconstruction()
{
  // Intermediate reconstruction to be done after setting the parameters
  // and before calculating the constraints.

  phi.SetTotalVector();

  for (Int_t i = 0; i < 2; i++)
  {
    charged_mom(pionCh[i].trackParams[0], pionCh[i].trackParams[1], pionCh[i].trackParams[2], pionCh[i].fourMom.data(), 1);
  }

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
  }

  std::vector<Int_t> bestPairingIndexNeutral, bestPairingIndexCharged;

  neutRec.PhotonPairingToPi0WithOmega(photon, pionCh, bestPairingIndexNeutral, bestPairingIndexCharged, omega);
  neutRec.Pi0Reconstruction(pionNe);
}

Double_t ConstraintsOmega::FourMomConsvLAB(Double_t *x, Double_t *p)
{
  SetParameters(p);
  IntermediateReconstruction();

  // Conservation of the chosen component of the 4-momentum in the LAB system
  return (bhabha_mom[_chosen4MomComponent] -
          pionCh[0].fourMom[_chosen4MomComponent] -
          pionCh[1].fourMom[_chosen4MomComponent] -
          photon[0].fourMom[_chosen4MomComponent] -
          photon[1].fourMom[_chosen4MomComponent] -
          photon[2].fourMom[_chosen4MomComponent] -
          photon[3].fourMom[_chosen4MomComponent]);
};

Double_t ConstraintsOmega::PhotonPathConsvLAB(Double_t *x, Double_t *p)
{
  SetParameters(p);
  IntermediateReconstruction();

  return photon[_chosenPhoton].fourPos[3] - photon[_chosenPhoton].timeOfFlight;
}

Double_t ConstraintsOmega::MinvConsv(Double_t *x, Double_t *p)
{
  SetParameters(p);
  IntermediateReconstruction();

  std::map<std::string, Float_t> minvModes = {
      {"omega", omega.total[5]}};

  return minvModes[_minvMode] - mOmega; // MeV/c^2
}