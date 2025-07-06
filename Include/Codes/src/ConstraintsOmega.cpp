#include <ConstraintsOmega.h>

using namespace KLOE;

Double_t ConstraintsOmega::FourMomConsvLAB(Double_t *x, Double_t *p) const
{
  // Local variables
  Double_t
      cluster[4][5],   // Cluster in LAB: Xcl, Ycl, Zcl, Tcl, EneCl
      gamma_mom[4][4], // Photon momentum in LAB: Px, Py, Pz, E
      neu_vtx[3],      // Derived neutral vertex: Xneu, Yneu, Zneu
      bhabha_mom[4],   // Phi momentum average per run in LAB: Px, Py, Pz, Sqrt(S)
      Kchrec[4];       // Charged Kaons' momentum in LAB: Px, Py, Pz, E

  for (Int_t i = 0; i < 4; i++)
  {
    for (Int_t j = 0; j < 5; j++)
      cluster[i][j] = p[i * 5 + j];

    if (i < 3)
      neu_vtx[i] = p[20 + i];

    Kchrec[i] = p[23 + i];
    bhabha_mom[i] = p[27 + i];

    // Reconstruction of neutral momentum for the photons
    neutral_mom(cluster[i][0], cluster[i][1], cluster[i][2], cluster[i][4], neu_vtx, gamma_mom[i]);
  }

  return (bhabha_mom[_chosen4MomComponent] -
          Kchrec[_chosen4MomComponent] -
          gamma_mom[0][_chosen4MomComponent] -
          gamma_mom[1][_chosen4MomComponent] -
          gamma_mom[2][_chosen4MomComponent] -
          gamma_mom[3][_chosen4MomComponent]);
}

Double_t ConstraintsOmega::PhotonPathConsvLAB(Double_t *x, Double_t *p) const
{
  Double_t
      cluster[4],   // Cluster in LAB: Xcl, Ycl, Zcl, Tcl, EneCl
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

  return constants.cVel * cluster[3] - R_gamma;
}