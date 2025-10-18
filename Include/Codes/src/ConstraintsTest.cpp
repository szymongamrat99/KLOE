#include <ConstraintsTest.h>

using namespace KLOE;

void ConstraintsTest::ResetParameters()
{

  photon.resize(2);
  photon.clear();
}

void ConstraintsTest::SetParameters(Double_t *p)
{
  // Setting parameters from the array p to the member variables
  // of the base classes to be used in the constraint calculations.
  ResetParameters();

  for (Int_t i = 0; i < 2; i++)
  {
    photon[i].fourMom[3] = p[i];
  }
}

void ConstraintsTest::IntermediateReconstruction()
{
}

Double_t ConstraintsTest::FourMomConsvLAB(Double_t *x, Double_t *p)
{
  Float_t p4[36];

  for (Int_t i = 0; i < 36; i++)
    p4[i] = p[i];

  SetParameters(p);
  IntermediateReconstruction();

  return (phi.fourMom[_chosen4MomComponent] -
          Kchboost.fourMom[_chosen4MomComponent] -
          photon[0].fourMom[_chosen4MomComponent] -
          photon[1].fourMom[_chosen4MomComponent] -
          photon[2].fourMom[_chosen4MomComponent] -
          photon[3].fourMom[_chosen4MomComponent]);
}

Double_t ConstraintsTest::PhotonPathConsvLAB(Double_t *x, Double_t *p)
{
  Float_t p4[36];

  for (Int_t i = 0; i < 36; i++)
    p4[i] = p[i];

  SetParameters(p);
  IntermediateReconstruction();

  return photon[_chosenPhoton].fourPos[3] - Knereclor.lifetimeLAB - photon[_chosenPhoton].timeOfFlight;
}

Double_t ConstraintsTest::MinvConsv(Double_t *x, Double_t *p)
{
  Float_t p4[3];

  for (Int_t i = 0; i < 3; i++)
    p4[i] = p[i];

  SetParameters(p);
  IntermediateReconstruction();

  std::map<std::string, Float_t> minvModes = {
      {"neutral", 2 * photon[0].fourMom[3] * photon[1].fourMom[3] * (1. - cos(p[2]))}};

  return minvModes[_minvMode] - pow(PhysicsConstants::mPi0, 2); // MeV/c^2
}