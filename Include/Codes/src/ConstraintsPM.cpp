#include <ConstraintsPM.h>

using namespace KLOE;

void ConstraintsPM::ResetParameters()
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

  if (ip.size() != 3)
  {
    ip.clear();
    ip.resize(3);
  }
}

void ConstraintsPM::SetParameters(Double_t *p)
{
  // Setting parameters from the array p to the member variables
  // of the base classes to be used in the constraint calculations.
  ResetParameters();

  for (Int_t i = 0; i < 2; i++)
  {
    fpionCh[i].trackParams[0] = p[i * 3];
    fpionCh[i].trackParams[1] = p[i * 3 + 1];
    fpionCh[i].trackParams[2] = p[i * 3 + 2];
  }

  for (Int_t i = 0; i < 4; i++)
  {
    fphi.fourMom[i] = p[6 + i];

    if (i < 3)
    {
      fphi.vtxPos[i] = p[10 + i];
    };
  }
}

void ConstraintsPM::IntermediateReconstruction()
{
  // Intermediate reconstruction to be done after setting the parameters
  // and before calculating the constraints.
  for (Int_t i = 0; i < 2; i++)
  {
    charged_mom(pionCh[i].trackParams[0], pionCh[i].trackParams[1], pionCh[i].trackParams[2], pionCh[i].fourMom.data(), 1, _logger);
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

  Double_t X_line[3] = {Kchboost.fourPos[0],
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
  if (std::abs(ip[2] - phi.vtxPos[2]) > 2.)
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

void ConstraintsPM::IntermediateReconstruction(Double_t *p)
{
  if (fip.size() != 3)
  {
    fip.clear();
    fip.resize(3);
  }

  for (Int_t i = 0; i < 2; i++)
  {
    fpionCh.at(i).fourMom[0] = p[i * 3];
    fpionCh.at(i).fourMom[1] = p[i * 3 + 1];
    fpionCh.at(i).fourMom[2] = p[i * 3 + 2];
    fpionCh.at(i).fourMom[3] = std::sqrt(std::pow(fpionCh.at(i).fourMom[0], 2) +
                                         std::pow(fpionCh.at(i).fourMom[1], 2) +
                                         std::pow(fpionCh.at(i).fourMom[2], 2) +
                                         std::pow(PhysicsConstants::mPiCh, 2));
  }

  for (Int_t i = 0; i < 4; i++)
  {
    fphi.fourMom[i] = p[6 + i];

    if (i < 3)
    {
      fphi.vtxPos[i] = p[10 + i];
    };
  }

  // Setting four momentum for kaon charged
  for (Int_t i = 0; i < 4; i++)
    fKchrec.fourMom[i] = fpionCh.at(0).fourMom[i] + fpionCh.at(1).fourMom[i];

  // Setting total vector for kaon charged reconstructed from pions
  fKchrec.SetLorentzVectors();
  fKchrec.SetTotalVector();

  KaonMomFromBoost(fKchrec.total.data(), fphi.fourMom.data(), fKchboost.total.data());
  fKchboost.SetPositionAndMomentumFromTotal();
  fKchboost.SetLorentzVectors();

  Double_t X_line[3] = {fKchboost.fourPos[0],
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
  if (std::abs(fip[2] - fphi.vtxPos[2]) > 2.8)
    fip[2] = fphi.vtxPos[2];

  fKchrec.calculatePath(fip.data());
  fKchrec.SetTotalVector();

  fKchboost.calculatePath(fip.data());
  fKchboost.SetTotalVector();

  if (fphi.fourMom[3] != 0.)
  {

    Double_t boostVec[3] = {-fphi.fourMom[0] / fphi.fourMom[3],
                            -fphi.fourMom[1] / fphi.fourMom[3],
                            -fphi.fourMom[2] / fphi.fourMom[3]};

    lorentz_transf(boostVec,
                   fKchrec.fourMom.data(),
                   fKchrecCMPhi.fourMom.data());
  }
  else
  {
    for (Int_t j = 0; j < 4; j++)
      fKchrecCMPhi.fourMom[j] = 999.;
  }
}

Double_t ConstraintsPM::EnergyConsvCM(Double_t *x, Double_t *p)
{
  IntermediateReconstruction(p);

  return fKchrecCMPhi.fourMom[3] - (fphi.fourMom[3] / 2.);
}

Double_t ConstraintsPM::MinvConsv(Double_t *x, Double_t *p)
{
  IntermediateReconstruction(p);

  return fKchrec.total[5] - PhysicsConstants::mK0; // MeV/c^2
}