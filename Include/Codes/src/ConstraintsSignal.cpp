#include <ConstraintsSignal.h>

#include "trilateration.hpp"

using namespace KLOE;

Double_t ConstraintsSignal::FourMomConsvLAB(Double_t *x, Double_t *p) const
{
  // Local variables
  Float_t
      gamma_mom[4][4], // Photon momentum in LAB: Px, Py, Pz, E 
      Curv[2],
      Phiv[2],
      Cotv[2],
      trk[2][4],
      bhabha_vtx[3],
      neu_vtx[3]; // Charged Kaons' momentum in LAB: Px, Py, Pz, E

  std::vector<Float_t> 
            cluster[4],  // Cluster in LAB: Xcl, Ycl, Zcl, TclOld, EneCl
            bhabha_mom(4), // Phi momentum average per run in LAB: Px, Py, Pz, Sqrt(S)
            Kchrec(9),
            Kchboost(9),
            ip(3);


  for (Int_t i = 0; i < 2; i++)
  {
    Curv[i] = p[i * 3];
    Phiv[i] = p[i * 3 + 1];
    Cotv[i] = p[i * 3 + 2];
  }

  for (Int_t i = 0; i < 4; i++)
  {
    cluster[i].resize(5);

    for (Int_t j = 0; j < 5; j++)
      cluster[i][j] = p[6 + i * 5 + j];

    bhabha_mom[i] = p[27 + i];
  }

  bhabha_vtx[0] = p[31];
  bhabha_vtx[1] = p[32];
  bhabha_vtx[2] = p[33];

  this->charged_mom(Curv[0], Phiv[0], Cotv[0], trk[0], 1);
  this->charged_mom(Curv[1], Phiv[1], Cotv[1], trk[1], 1);

  Kchrec[0] = trk[0][0] + trk[1][0];
  Kchrec[1] = trk[0][1] + trk[1][1];
  Kchrec[2] = trk[0][2] + trk[1][2];
  Kchrec[3] = trk[0][3] + trk[1][3];

  this->KaonMomFromBoost(Kchrec, bhabha_mom.data(), Kchboost);

  Float_t X_line[3] = {Kchboost[6],
                       Kchboost[7],
                       Kchboost[8]}, // Vertex laying on the line
      mom[3] = {Kchboost[0],
              Kchboost[1],
              Kchboost[2]}, // Direction of the line
      xB[3] = {bhabha_vtx[0],
               bhabha_vtx[1],
               bhabha_vtx[2]}, // Bhabha vertex - laying on the plane
      plane_perp[3] = {0.,
                       bhabha_mom[1],
                       0.}; // Vector perpendicular to the plane from Bhabha momentum

  // Corrected IP event by event
  this->IPBoostCorr(X_line, mom, xB, plane_perp, ip);

  ip[0] = bhabha_vtx[0];
  ip[1] = bhabha_vtx[1];
  if(abs(ip[2] - bhabha_vtx[2]) > 2.0)
    ip[2] = bhabha_vtx[2];

  std::vector<Float_t> 
                Knetriangle, 
                gammatriangle[4],
                trcfinal;
  Float_t minv4gam;

  // this->triangleReconstruction(_g4taken, cluster, _neuclulist, bhabha_mom, Kchboost, ip, Knetriangle, gammatriangle, minv4gam, trcfinal, *_logger);

  // Reconstruction of neutral momentum for the photons
  for (Int_t i = 0; i < 4; i++)
  {
    neutral_mom(cluster[i][0], cluster[i][1], cluster[i][2], cluster[i][4], neu_vtx, gamma_mom[i]);
  }

  return (bhabha_mom[_chosen4MomComponent] -
          Kchrec[_chosen4MomComponent] -
          gamma_mom[0][_chosen4MomComponent] -
          gamma_mom[1][_chosen4MomComponent] -
          gamma_mom[2][_chosen4MomComponent] -
          gamma_mom[3][_chosen4MomComponent]);
}

Double_t ConstraintsSignal::PhotonPathConsvLAB(Double_t *x, Double_t *p) const
{
  Double_t
      cluster[4],   // Cluster in LAB: Xcl, Ycl, Zcl, TclOld, EneCl
      neu_vtx[3],   // Derived neutral vertex: Xneu, Yneu, Zneu
      R_gamma = 0.; // Path of Photon from the neutral vtx

  for (Int_t i = 0; i < 4; i++)
  {
    cluster[i] = p[_chosenPhoton * 5 + i];

    if (i < 3)
      neu_vtx[i] = p[20 + i];
  }

  R_gamma = sqrt(pow(cluster[0] - neu_vtx[0], 2) +
                 pow(cluster[1] - neu_vtx[1], 2) +
                 pow(cluster[2] - neu_vtx[2], 2));

  return cVel * cluster[3] - R_gamma;
}