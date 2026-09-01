#include <ConstraintsNeutral.h>

using namespace KLOE;

void ConstraintsNeutral::ResetParameters()
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

void ConstraintsNeutral::SetParameters(Double_t *p)
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
    fphi.fourMom[i] = p[29 + i];

    if (i < 3)
    {
      fKnereclor.fourPos[i] = p[33 + i];
      fphi.vtxPos[i] = p[36 + i];
    };
  }
}

void ConstraintsNeutral::IntermediateReconstruction(Double_t *p)
{
  if (fip.size() != 3)
  {
    fip.clear();
    fip.resize(3);
  }

  int _offset = 0;

  for (Int_t i = 0; i < 4; i++)
  {
    for (Int_t j = 0; j < 5; j++)
    {
      fphoton.at(i).clusterParams[j] = p[_offset + i * 5 + j];
    }
  }
  _offset += 4 * 5;

  for (Int_t i = 0; i < 3; i++)
    fKnereclor.fourPos[i] = p[_offset + i];
  _offset += 3;

  for (Int_t i = 0; i < 4; i++)
    fphi.fourMom[i] = p[_offset + i];
  _offset += 4;

  for (Int_t i = 0; i < 3; i++)
    fphi.vtxPos[i] = p[_offset + i];
  _offset += 3;

  // === Setting ended === 

  fip[0] = fphi.vtxPos[0];
  fip[1] = fphi.vtxPos[1];
  fip[2] = fphi.vtxPos[2];

  for (Int_t i = 0; i < 4; i++)
  {
    neutral_mom(fphoton.at(i).clusterParams[0],
                fphoton.at(i).clusterParams[1],
                fphoton.at(i).clusterParams[2],
                fphoton.at(i).clusterParams[4],
                fKnereclor.fourPos.data(),
                fphoton.at(i).fourMom.data());

    fphoton.at(i).fourPos[0] = fphoton.at(i).clusterParams[0];
    fphoton.at(i).fourPos[1] = fphoton.at(i).clusterParams[1];
    fphoton.at(i).fourPos[2] = fphoton.at(i).clusterParams[2];
    fphoton.at(i).fourPos[3] = fphoton.at(i).clusterParams[3];

    fphoton.at(i).calculatePath(fKnereclor.fourPos.data());
    fphoton.at(i).calculateTimeOfFlightPhoton();
    fphoton.at(i).SetTotalVectorPhoton();
  }

  fKnerec.fourMom[0] = fphoton.at(0).fourMom[0] + fphoton.at(1).fourMom[0] + fphoton.at(2).fourMom[0] + fphoton.at(3).fourMom[0];
  fKnerec.fourMom[1] = fphoton.at(0).fourMom[1] + fphoton.at(1).fourMom[1] + fphoton.at(2).fourMom[1] + fphoton.at(3).fourMom[1];
  fKnerec.fourMom[2] = fphoton.at(0).fourMom[2] + fphoton.at(1).fourMom[2] + fphoton.at(2).fourMom[2] + fphoton.at(3).fourMom[2];
  fKnerec.fourMom[3] = fphoton.at(0).fourMom[3] + fphoton.at(1).fourMom[3] + fphoton.at(2).fourMom[3] + fphoton.at(3).fourMom[3];

  fKnerec.fourPos[0] = fKnereclor.fourPos[0];
  fKnerec.fourPos[1] = fKnereclor.fourPos[1];
  fKnerec.fourPos[2] = fKnereclor.fourPos[2];
  fKnerec.fourPos[3] = fKnereclor.fourPos[3];

  fKnerec.calculatePath(fip.data());
  fKnerec.SetTotalVector();

  if (fphi.fourMom[3] != 0.)
  {

    Double_t boostVec[3] = {-fphi.fourMom[0] / fphi.fourMom[3],
                            -fphi.fourMom[1] / fphi.fourMom[3],
                            -fphi.fourMom[2] / fphi.fourMom[3]};

    lorentz_transf(boostVec,
                   fKnerec.fourMom.data(),
                   fKnerecCMPhi.fourMom.data());
  }
  else
  {
    for (Int_t j = 0; j < 4; j++)
      fKnerecCMPhi.fourMom[j] = 999.;
  }
}

Double_t ConstraintsNeutral::EnergyConsvCM(Double_t *x, Double_t *p)
{
  IntermediateReconstruction(p);

  return (fphi.fourMom[3] / 2. - fKnerecCMPhi.fourMom[3]);
}

Double_t ConstraintsNeutral::PhotonPathConsvLAB(Double_t *x, Double_t *p)
{
  IntermediateReconstruction(p);

  return fphoton.at(_chosenPhoton).fourPos[3] - fKnerec.lifetimeLAB - fphoton.at(_chosenPhoton).timeOfFlight;
}

Double_t ConstraintsNeutral::MinvConsv(Double_t *x, Double_t *p)
{
  IntermediateReconstruction(p);

  std::map<std::string, Double_t> minvModes = {
      {"neutral", fKnerec.total[5]}};

  return minvModes[_minvMode] - PhysicsConstants::mK0; // MeV/c^2
}