#include <ConstraintsSignal.h>

using namespace KLOE;

void ConstraintsSignal::ResetParameters()
{
  // Setting parameters from the array p to the member variables
  // of the base classes to be used in the constraint calculations.
  if (bhabha_mom.size() != 4) {
        bhabha_mom.clear();
        bhabha_mom.resize(4);
    }
    
    if (trackParameters.size() != 6) {
        trackParameters.clear();
        trackParameters.resize(6);
    }
    
    if (cluster.size() != 4) {
        for (Int_t i = 0; i < cluster.size(); i++)
            cluster[i].clear();
        cluster.clear();
        cluster.resize(4);
    }
    
    for (Int_t i = 0; i < cluster.size(); i++) {
        if (cluster[i].size() != 5) {
            cluster[i].clear();
            cluster[i].resize(5);
        }
    }

    if (pionCh.size() != 2) {
        pionCh.clear();
        pionCh.resize(2);
    }
    
    if (photon.size() != 4) {
        photon.clear();
        photon.resize(4);
    }

    if (ip.size() != 3) {
        ip.clear();
        ip.resize(3);
    }
}

void ConstraintsSignal::SetParameters(Float_t *p)
{
  // Setting parameters from the array p to the member variables
  // of the base classes to be used in the constraint calculations.
  ResetParameters();

  for (Int_t i = 0; i < 2; i++)
  {
    pionCh[i].trackParams[0] = p[i * 3];
    pionCh[i].trackParams[1] = p[i * 3 + 1];
    pionCh[i].trackParams[2] = p[i * 3 + 2];
  }

  for (Int_t i = 0; i < 4; i++)
  {
    for (Int_t j = 0; j < 5; j++)
      photon[i].clusterParams[j] = p[6 + i * 5 + j];
  }

  for (Int_t i = 0; i < 3; i++)
    Kchrec.fourPos[i] = p[26 + i];

  for (Int_t i = 0; i < 4; i++)
  {
    phi.fourMom[i] = p[29 + i];

    if (i < 3)
      phi.vtxPos[i] = p[33 + i];
  }
}

void ConstraintsSignal::IntermediateReconstruction()
{
  // Intermediate reconstruction to be done after setting the parameters
  // and before calculating the constraints.

  for (Int_t i = 0; i < 2; i++)
    charged_mom(pionCh[i].trackParams[0], pionCh[i].trackParams[1], pionCh[i].trackParams[2], pionCh[i].fourMom.data(), 1);

  // Setting four momentum for kaon charged
  for (Int_t i = 0; i < 4; i++)
    Kchrec.fourMom[i] = pionCh[0].fourMom[i] + pionCh[1].fourMom[i];

  // Setting total vector for kaon charged reconstructed from pions
  Kchrec.SetLorentzVectors();
  Kchrec.SetTotalVector();

  KaonMomFromBoost(Kchrec.total.data(), phi.fourMom.data(), Kchboost.total.data());
  Kchboost.SetPositionAndMomentumFromTotal();
  Kchboost.SetLorentzVectors();

  Float_t X_line[3] = {Kchboost.fourPos[0],
                       Kchboost.fourPos[1],
                       Kchboost.fourPos[2]}, // Vertex laying on the line
      mom[3] = {Kchboost.fourMom[0],
                Kchboost.fourMom[1],
                Kchboost.fourMom[2]}, // Direction of the line
      xB[3] = {phi.vtxPos[0],
               phi.vtxPos[1],
               phi.vtxPos[2]}, // Bhabha vertex - laying on the plane
      plane_perp[3] = {0.,
                       phi.fourMom[1],
                       0.}; // Vector perpendicular to the plane from Bhabha momentum

  // Corrected IP event by event
  IPBoostCorr(X_line, mom, xB, plane_perp, ip.data());

  Kchrec.calculatePath(ip.data());
  Kchrec.SetTotalVector();

  Kchboost.calculatePath(ip.data());
  Kchboost.SetTotalVector();

  triangleReconstruction(photon, phi, Kchboost, ip.data(), Knereclor);

  Knereclor.calculatePath(ip.data());
  Knereclor.SetTotalVector();

  for (Int_t i = 0; i < 4; i++)
  {
    photon[i].fourPos[0] = photon[i].clusterParams[0];
    photon[i].fourPos[1] = photon[i].clusterParams[1];
    photon[i].fourPos[2] = photon[i].clusterParams[2];
    photon[i].fourPos[3] = photon[i].clusterParams[3];
    photon[i].calculatePath(Knereclor.fourPos.data());
    photon[i].calculateTimeOfFlightPhoton();
    photon[i].SetTotalVectorPhoton();
  }

  Knerec.fourMom[0] = photon[0].fourMom[0] + photon[1].fourMom[0] + photon[2].fourMom[0] + photon[3].fourMom[0];
  Knerec.fourMom[1] = photon[0].fourMom[1] + photon[1].fourMom[1] + photon[2].fourMom[1] + photon[3].fourMom[1];
  Knerec.fourMom[2] = photon[0].fourMom[2] + photon[1].fourMom[2] + photon[2].fourMom[2] + photon[3].fourMom[2];
  Knerec.fourMom[3] = photon[0].fourMom[3] + photon[1].fourMom[3] + photon[2].fourMom[3] + photon[3].fourMom[3];

  Knerec.fourPos[0] = Knereclor.fourPos[0];
  Knerec.fourPos[1] = Knereclor.fourPos[1];
  Knerec.fourPos[2] = Knereclor.fourPos[2];
  Knerec.fourPos[3] = Knereclor.fourPos[3];

  Knerec.calculatePath(ip.data());
  Knerec.SetTotalVector();
}

Double_t ConstraintsSignal::FourMomConsvLAB(Double_t *x, Double_t *p)
{
  Float_t p4[36];

  for (Int_t i = 0; i < 36; i++)
    p4[i] = p[i];

  SetParameters(p4);
  IntermediateReconstruction();

  return (phi.fourMom[_chosen4MomComponent] -
          Kchboost.fourMom[_chosen4MomComponent] -
          photon[0].fourMom[_chosen4MomComponent] -
          photon[1].fourMom[_chosen4MomComponent] -
          photon[2].fourMom[_chosen4MomComponent] -
          photon[3].fourMom[_chosen4MomComponent]);
}

Double_t ConstraintsSignal::PhotonPathConsvLAB(Double_t *x, Double_t *p)
{
  Float_t p4[36];

  for (Int_t i = 0; i < 36; i++)
    p4[i] = p[i];

  SetParameters(p4);
  IntermediateReconstruction();

  return photon[_chosenPhoton].fourPos[3] - Knereclor.lifetimeLAB - photon[_chosenPhoton].timeOfFlight;
}

Double_t ConstraintsSignal::MinvConsv(Double_t *x, Double_t *p)
{
  Float_t p4[36];

  for (Int_t i = 0; i < 36; i++)
    p4[i] = p[i];

  SetParameters(p4);
  IntermediateReconstruction();

  std::map<std::string, Float_t> minvModes = {
      {"neutral", Knerec.total[5]},
      {"charged", Kchrec.total[5]}};

  return minvModes[_minvMode] - mK0; // MeV/c^2
}