#include <ConstraintsSignal.h>

using namespace KLOE;

void ConstraintsSignal::ResetParameters()
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

void ConstraintsSignal::SetParameters(Double_t *p)
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

  for (Int_t i = 0; i < 3; i++)
    fKchrec.fourPos[i] = p[20 + i];

  for (Int_t i = 0; i < 2; i++)
  {
    fpionCh[i].trackParams[0] = p[23 + i * 3];
    fpionCh[i].trackParams[1] = p[23 + i * 3 + 1];
    fpionCh[i].trackParams[2] = p[23 + i * 3 + 2];
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

void ConstraintsSignal::IntermediateReconstruction()
{
  // Intermediate reconstruction to be done after setting the parameters
  // and before calculating the constraints.
  for (Int_t i = 0; i < 2; i++)
  {
    charged_mom(pionCh[i].trackParams[0], pionCh[i].trackParams[1], pionCh[i].trackParams[2], pionCh[i].fourMom.data(), 1);
  }

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

  ip[0] = phi.vtxPos[0];
  ip[1] = phi.vtxPos[1];
  // // ip[2] is fitted
  if (abs(ip[2] - phi.vtxPos[2]) > 2.)
    ip[2] = phi.vtxPos[2];

  Kchrec.calculatePath(ip.data());
  Kchrec.SetTotalVector();

  Kchboost.calculatePath(ip.data());
  Kchboost.SetTotalVector();

  // triangleReconstruction(photon, phi, Kchboost, ip.data(), Knereclor);

  Knereclor.fourMom[0] = phi.fourMom[0] - Kchboost.fourMom[0];
  Knereclor.fourMom[1] = phi.fourMom[1] - Kchboost.fourMom[1];
  Knereclor.fourMom[2] = phi.fourMom[2] - Kchboost.fourMom[2];
  Knereclor.fourMom[3] = phi.fourMom[3] - Kchboost.fourMom[3];

  Knereclor.calculatePath(ip.data());
  Knereclor.SetTotalVector();

  for (Int_t i = 0; i < 4; i++)
  {
    neutral_mom(photon[i].clusterParams[0],
                photon[i].clusterParams[1],
                photon[i].clusterParams[2],
                photon[i].clusterParams[4],
                Knereclor.fourPos.data(),
                photon[i].fourMom.data());

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

void ConstraintsSignal::IntermediateReconstruction(Double_t *p)
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

  for (Int_t i = 0; i < 3; i++)
    fKchrec.fourPos[i] = p[20 + i];

  for (Int_t i = 0; i < 2; i++)
  {
    fpionCh[i].fourMom[0] = p[23 + i * 3];
    fpionCh[i].fourMom[1] = p[23 + i * 3 + 1];
    fpionCh[i].fourMom[2] = p[23 + i * 3 + 2];
    fpionCh[i].fourMom[3] = sqrt(pow(fpionCh[i].fourMom[0], 2) +
                                 pow(fpionCh[i].fourMom[1], 2) +
                                 pow(fpionCh[i].fourMom[2], 2) +
                                 pow(PhysicsConstants::mPiCh, 2));
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

  // Setting four momentum for kaon charged
  for (Int_t i = 0; i < 4; i++)
    fKchrec.fourMom[i] = fpionCh[0].fourMom[i] + fpionCh[1].fourMom[i];

  // Setting total vector for kaon charged reconstructed from pions
  fKchrec.SetLorentzVectors();
  fKchrec.SetTotalVector();

  KaonMomFromBoost(fKchrec.total.data(), fphi.fourMom.data(), fKchboost.total.data());
  fKchboost.SetPositionAndMomentumFromTotal();
  fKchboost.SetLorentzVectors();

  Float_t X_line[3] = {fKchboost.fourPos[0],
                       fKchboost.fourPos[1],
                       fKchboost.fourPos[2]}, // Vertex laying on the line
      mom[3] = {fKchboost.fourMom[0],
                fKchboost.fourMom[1],
                fKchboost.fourMom[2]}, // Direction of the line
      xB[3] = {fphi.vtxPos[0],
               fphi.vtxPos[1],
               fphi.vtxPos[2]}, // Bhabha vertex - laying on the plane
      plane_perp[3] = {fphi.fourMom[0],
                       fphi.fourMom[1],
                       0.}; // Vector perpendicular to the plane from Bhabha momentum

  // Corrected IP event by event
  IPBoostCorr(X_line, mom, xB, plane_perp, fip.data());

  fip[0] = fphi.vtxPos[0];
  fip[1] = fphi.vtxPos[1];
  // fip[2] is fitted
  if (abs(fip[2] - fphi.vtxPos[2]) > 2.8)
    fip[2] = fphi.vtxPos[2];

  fKchrec.calculatePath(fip.data());
  fKchrec.SetTotalVector();

  fKchboost.calculatePath(fip.data());
  fKchboost.SetTotalVector();

  fKnereclor.fourMom[0] = fphi.fourMom[0] - fKchboost.fourMom[0];
  fKnereclor.fourMom[1] = fphi.fourMom[1] - fKchboost.fourMom[1];
  fKnereclor.fourMom[2] = fphi.fourMom[2] - fKchboost.fourMom[2];
  fKnereclor.fourMom[3] = fphi.fourMom[3] - fKchboost.fourMom[3];

  fKnereclor.calculatePath(fip.data());
  fKnereclor.SetTotalVector();

  for (Int_t i = 0; i < 4; i++)
  {
    neutral_mom(fphoton[i].clusterParams[0],
                fphoton[i].clusterParams[1],
                fphoton[i].clusterParams[2],
                fphoton[i].clusterParams[4],
                fKnereclor.fourPos.data(),
                fphoton[i].fourMom.data());

    fphoton[i].fourPos[0] = fphoton[i].clusterParams[0];
    fphoton[i].fourPos[1] = fphoton[i].clusterParams[1];
    fphoton[i].fourPos[2] = fphoton[i].clusterParams[2];
    fphoton[i].fourPos[3] = fphoton[i].clusterParams[3];

    fphoton[i].calculatePath(fKnereclor.fourPos.data());
    fphoton[i].calculateTimeOfFlightPhoton();
    fphoton[i].SetTotalVectorPhoton();
  }

  fKnerec.fourMom[0] = fphoton[0].fourMom[0] + fphoton[1].fourMom[0] + fphoton[2].fourMom[0] + fphoton[3].fourMom[0];
  fKnerec.fourMom[1] = fphoton[0].fourMom[1] + fphoton[1].fourMom[1] + fphoton[2].fourMom[1] + fphoton[3].fourMom[1];
  fKnerec.fourMom[2] = fphoton[0].fourMom[2] + fphoton[1].fourMom[2] + fphoton[2].fourMom[2] + fphoton[3].fourMom[2];
  fKnerec.fourMom[3] = fphoton[0].fourMom[3] + fphoton[1].fourMom[3] + fphoton[2].fourMom[3] + fphoton[3].fourMom[3];

  fKnerec.fourPos[0] = fKnereclor.fourPos[0];
  fKnerec.fourPos[1] = fKnereclor.fourPos[1];
  fKnerec.fourPos[2] = fKnereclor.fourPos[2];
  fKnerec.fourPos[3] = fKnereclor.fourPos[3];

  fKnerec.calculatePath(fip.data());
  fKnerec.SetTotalVector();
}

Double_t ConstraintsSignal::FourMomConsvLAB(Double_t *x, Double_t *p)
{
  // SetParameters(p);
  IntermediateReconstruction(p);

  return (fphi.fourMom[_chosen4MomComponent] -
          fpionCh[0].fourMom[_chosen4MomComponent] -
          fpionCh[1].fourMom[_chosen4MomComponent] -
          fphoton[0].fourMom[_chosen4MomComponent] -
          fphoton[1].fourMom[_chosen4MomComponent] -
          fphoton[2].fourMom[_chosen4MomComponent] -
          fphoton[3].fourMom[_chosen4MomComponent]);
}

Double_t ConstraintsSignal::PhotonPathConsvLAB(Double_t *x, Double_t *p)
{
  // SetParameters(p);
  IntermediateReconstruction(p);

  return fphoton[_chosenPhoton].fourPos[3] - fKnerec.lifetimeLAB - fphoton[_chosenPhoton].timeOfFlight;
}

Double_t ConstraintsSignal::MinvConsv(Double_t *x, Double_t *p)
{
  // SetParameters(p);
  IntermediateReconstruction(p);

  std::map<std::string, Float_t> minvModes = {
      {"neutral", fKnerec.total[5]},
      {"charged", fKchrec.total[5]}};

  return minvModes[_minvMode] - PhysicsConstants::mK0; // MeV/c^2
}