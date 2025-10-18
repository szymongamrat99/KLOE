#include <TMath.h>
#include <TMatrixD.h>
#include <TVectorD.h>

#include "kinematic_fits.h"
#include "charged_mom.h"
#include "neutral_mom.h"
#include "../const.h"

using namespace KLOE;

void kin_fits::FillConstructors(Double_t *par, Double_t *sigmas)
{
  for(Int_t i = 0; i < _N; i++)
  {
    _P(i) = par[i];
    _P0(i) = par[i];

    for(Int_t j = 0; j < _N; j++)
    {
      if(i == j) _V(i,j) = pow(sigmas[i],2);
      else _V(i,j) = 0.;

      if(i == j) _V0(i,j) = pow(sigmas[i],2);
      else _V0(i,j) = 0.;
    }
  }
}

void kin_fits::cons_vect()
{
  for(Int_t i = 0; i < _M; i++)
  {
    _constraints[i]->SetParameters(_P.GetMatrixArray());
    _C(i) = _constraints[i]->Eval(0.);
  }
}

void kin_fits::cons_diff()
{
  for(Int_t i = 0; i < _M; i++)
  {
    _constraints[i]->SetParameters(_P.GetMatrixArray());

    for(Int_t j = 0; j < _N; j++)
    {
      _D(i,j) = _constraints[i]->GradientPar(j, 0, 0.001);
    }
  }

  _D_T.Transpose(_D);
}

void kin_fits::solution()
{
  for(Int_t i = 0; i < _iter_num; i++)
  {
    cons_vect();
    cons_diff();

    _Aux = _D * _V * _D_T;

    _det = _Aux.Determinant();

    if(_det != 0.)
    { 
      _Aux.Invert();
      _Lambda = _Aux * _C;

      _Corr = _V * _D_T * _Lambda;

      if( _Corr.Abs() < 1000.)
      {
        _P1 = _P - _V * _D_T * _Lambda;
        _V_updated = _V * _D_T * _Aux * _D * _V;

        _P = _P1;
        _V = _V_updated;

        _fail = 0;
      }
      else{
        _fail = 1;
        break;
      }
    }
    else
    {
      _fail = 1;
      break;
    };
  }
  _chi2 = Dot(_P - _P0, _V0.Invert() * (_P - _P0) );
}


////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////

Double_t kin_fits::ene_consv_pi(const Double_t *xx, const Double_t *pp)
{
  Float_t Phi4mom[4] = {pp[0], pp[1], pp[2], pp[3]};
  Float_t mom_ch[2][4], mom_neu[4][4];
  Float_t neu_vtx[3] = {pp[30], pp[31], pp[32]};

  charged_mom(pp[4], pp[5], pp[6], mom_ch[0], 1);
  charged_mom(pp[7], pp[8], pp[9], mom_ch[1], 1);

  neutral_mom(pp[10], pp[11], pp[12], pp[14], neu_vtx, mom_neu[0]);
  neutral_mom(pp[15], pp[16], pp[17], pp[19], neu_vtx, mom_neu[1]);
  neutral_mom(pp[20], pp[21], pp[22], pp[24], neu_vtx, mom_neu[2]);
  neutral_mom(pp[25], pp[26], pp[27], pp[29], neu_vtx, mom_neu[3]);

  Double_t value = 0.;

  value = mom_ch[0][3] + mom_ch[1][3] + mom_neu[0][3] + mom_neu[1][3] + mom_neu[2][3] + mom_neu[3][3] - Phi4mom[3];

  return value; 
}

Double_t kin_fits::mom_x_consv_pi(const Double_t *xx, const Double_t *pp)
{
  Float_t Phi4mom[4] = {pp[0], pp[1], pp[2], pp[3]};
  Float_t mom_ch[2][4], mom_neu[4][4];
  Float_t neu_vtx[3] = {pp[30], pp[31], pp[32]};

  charged_mom(pp[4], pp[5], pp[6], mom_ch[0], 1);
  charged_mom(pp[7], pp[8], pp[9], mom_ch[1], 1);

  neutral_mom(pp[10], pp[11], pp[12], pp[14], neu_vtx, mom_neu[0]);
  neutral_mom(pp[15], pp[16], pp[17], pp[19], neu_vtx, mom_neu[1]);
  neutral_mom(pp[20], pp[21], pp[22], pp[24], neu_vtx, mom_neu[2]);
  neutral_mom(pp[25], pp[26], pp[27], pp[29], neu_vtx, mom_neu[3]);

  Double_t value = 0.;

  value = mom_ch[0][0] + mom_ch[1][0] + mom_neu[0][0] + mom_neu[1][0] + mom_neu[2][0] + mom_neu[3][0] - Phi4mom[0];

  return value;
}

Double_t kin_fits::mom_y_consv_pi(const Double_t *xx, const Double_t *pp)
{
  Float_t Phi4mom[4] = {pp[0], pp[1], pp[2], pp[3]};
  Float_t mom_ch[2][4], mom_neu[4][4];
  Float_t neu_vtx[3] = {pp[30], pp[31], pp[32]};

  charged_mom(pp[4], pp[5], pp[6], mom_ch[0], 1);
  charged_mom(pp[7], pp[8], pp[9], mom_ch[1], 1);

  neutral_mom(pp[10], pp[11], pp[12], pp[14], neu_vtx, mom_neu[0]);
  neutral_mom(pp[15], pp[16], pp[17], pp[19], neu_vtx, mom_neu[1]);
  neutral_mom(pp[20], pp[21], pp[22], pp[24], neu_vtx, mom_neu[2]);
  neutral_mom(pp[25], pp[26], pp[27], pp[29], neu_vtx, mom_neu[3]);

  Double_t value = 0.;

  value = mom_ch[0][1] + mom_ch[1][1] + mom_neu[0][1] + mom_neu[1][1] + mom_neu[2][1] + mom_neu[3][1] - Phi4mom[1];

  return value;
}

Double_t kin_fits::mom_z_consv_pi(const Double_t *xx, const Double_t *pp)
{
  Float_t Phi4mom[4] = {pp[0], pp[1], pp[2], pp[3]};
  Float_t mom_ch[2][4], mom_neu[4][4];
  Float_t neu_vtx[3] = {pp[30], pp[31], pp[32]};

  charged_mom(pp[4], pp[5], pp[6], mom_ch[0], 1);
  charged_mom(pp[7], pp[8], pp[9], mom_ch[1], 1);

  neutral_mom(pp[10], pp[11], pp[12], pp[14], neu_vtx, mom_neu[0]);
  neutral_mom(pp[15], pp[16], pp[17], pp[19], neu_vtx, mom_neu[1]);
  neutral_mom(pp[20], pp[21], pp[22], pp[24], neu_vtx, mom_neu[2]);
  neutral_mom(pp[25], pp[26], pp[27], pp[29], neu_vtx, mom_neu[3]);

  Double_t value = 0.;

  value = mom_ch[0][2] + mom_ch[1][2] + mom_neu[0][2] + mom_neu[1][2] + mom_neu[2][2] + mom_neu[3][2] - Phi4mom[2];

  return value;
}

Double_t kin_fits::cluster_time_cons_first(const Double_t *xx, const Double_t *pp)
{
  Float_t Phi4mom[4] = {pp[0], pp[1], pp[2], pp[3]};
  Float_t mom_ch[2][4], mom_neu[4][4];
  Float_t neu_vtx[3] = {pp[30], pp[31], pp[32]};
  Float_t ip_vtx[3] = {pp[33], pp[34], pp[35]};

  charged_mom(pp[4], pp[5], pp[6], mom_ch[0], 1);
  charged_mom(pp[7], pp[8], pp[9], mom_ch[1], 1);

  neutral_mom(pp[10], pp[11], pp[12], pp[14], neu_vtx, mom_neu[0]);
  neutral_mom(pp[15], pp[16], pp[17], pp[19], neu_vtx, mom_neu[1]);
  neutral_mom(pp[20], pp[21], pp[22], pp[24], neu_vtx, mom_neu[2]);
  neutral_mom(pp[25], pp[26], pp[27], pp[29], neu_vtx, mom_neu[3]);

  Double_t path_gamma, time_gamma, path_kaon, vK, time_kaon, time_diff;

  path_kaon = sqrt(pow(neu_vtx[0] - ip_vtx[0],2) + pow(neu_vtx[1] - ip_vtx[1],2) + pow(neu_vtx[2] - ip_vtx[2],2));
  vK = PhysicsConstants::cVel*sqrt(pow(mom_neu[0][0] + mom_neu[1][0] + mom_neu[2][0] + mom_neu[3][0], 2) + 
                  pow(mom_neu[0][1] + mom_neu[1][1] + mom_neu[2][1] + mom_neu[3][1], 2) +
                  pow(mom_neu[0][2] + mom_neu[1][2] + mom_neu[2][2] + mom_neu[3][2], 2))/(mom_neu[0][3] + mom_neu[1][3] + mom_neu[2][3] + mom_neu[3][3]);

  time_kaon = path_kaon/vK;

  Double_t value = 0.;


  path_gamma = sqrt(pow(pp[10] - neu_vtx[0],2) + pow(pp[11] - neu_vtx[1],2) + pow(pp[12] - neu_vtx[2],2));
  time_gamma = path_gamma/PhysicsConstants::cVel;

  time_diff = time_gamma + time_kaon - pp[13];

  value = time_diff;

  return value;

}

Double_t kin_fits::cluster_time_cons_second(const Double_t *xx, const Double_t *pp)
{
  Float_t Phi4mom[4] = {pp[0], pp[1], pp[2], pp[3]};
  Float_t mom_ch[2][4], mom_neu[4][4];
  Float_t neu_vtx[3] = {pp[30], pp[31], pp[32]};
  Float_t ip_vtx[3] = {pp[33], pp[34], pp[35]};

  charged_mom(pp[4], pp[5], pp[6], mom_ch[0], 1);
  charged_mom(pp[7], pp[8], pp[9], mom_ch[1], 1);

  neutral_mom(pp[10], pp[11], pp[12], pp[14], neu_vtx, mom_neu[0]);
  neutral_mom(pp[15], pp[16], pp[17], pp[19], neu_vtx, mom_neu[1]);
  neutral_mom(pp[20], pp[21], pp[22], pp[24], neu_vtx, mom_neu[2]);
  neutral_mom(pp[25], pp[26], pp[27], pp[29], neu_vtx, mom_neu[3]);

  Double_t path_gamma, time_gamma, path_kaon, vK, time_kaon, time_diff;

  path_kaon = sqrt(pow(neu_vtx[0] - ip_vtx[0],2) + pow(neu_vtx[1] - ip_vtx[1],2) + pow(neu_vtx[2] - ip_vtx[2],2));
  vK = PhysicsConstants::cVel*sqrt(pow(mom_neu[0][0] + mom_neu[1][0] + mom_neu[2][0] + mom_neu[3][0], 2) + 
                  pow(mom_neu[0][1] + mom_neu[1][1] + mom_neu[2][1] + mom_neu[3][1], 2) +
                  pow(mom_neu[0][2] + mom_neu[1][2] + mom_neu[2][2] + mom_neu[3][2], 2))/(mom_neu[0][3] + mom_neu[1][3] + mom_neu[2][3] + mom_neu[3][3]);

  time_kaon = path_kaon/vK;

  Double_t value = 0.;

  path_gamma = sqrt(pow(pp[15] - neu_vtx[0],2) + pow(pp[16] - neu_vtx[1],2) + pow(pp[17] - neu_vtx[2],2));
  time_gamma = path_gamma/PhysicsConstants::cVel;

  time_diff = time_gamma + time_kaon - pp[18];

  value = time_diff;

  return value;

}

Double_t kin_fits::cluster_time_cons_third(const Double_t *xx, const Double_t *pp)
{
  Float_t Phi4mom[4] = {pp[0], pp[1], pp[2], pp[3]};
  Float_t mom_ch[2][4], mom_neu[4][4];
  Float_t neu_vtx[3] = {pp[30], pp[31], pp[32]};
  Float_t ip_vtx[3] = {pp[33], pp[34], pp[35]};

  charged_mom(pp[4], pp[5], pp[6], mom_ch[0], 1);
  charged_mom(pp[7], pp[8], pp[9], mom_ch[1], 1);

  neutral_mom(pp[10], pp[11], pp[12], pp[14], neu_vtx, mom_neu[0]);
  neutral_mom(pp[15], pp[16], pp[17], pp[19], neu_vtx, mom_neu[1]);
  neutral_mom(pp[20], pp[21], pp[22], pp[24], neu_vtx, mom_neu[2]);
  neutral_mom(pp[25], pp[26], pp[27], pp[29], neu_vtx, mom_neu[3]);

  Double_t path_gamma, time_gamma, path_kaon, vK, time_kaon, time_diff;

  path_kaon = sqrt(pow(neu_vtx[0] - ip_vtx[0],2) + pow(neu_vtx[1] - ip_vtx[1],2) + pow(neu_vtx[2] - ip_vtx[2],2));
  vK = PhysicsConstants::cVel*sqrt(pow(mom_neu[0][0] + mom_neu[1][0] + mom_neu[2][0] + mom_neu[3][0], 2) + 
                  pow(mom_neu[0][1] + mom_neu[1][1] + mom_neu[2][1] + mom_neu[3][1], 2) +
                  pow(mom_neu[0][2] + mom_neu[1][2] + mom_neu[2][2] + mom_neu[3][2], 2))/(mom_neu[0][3] + mom_neu[1][3] + mom_neu[2][3] + mom_neu[3][3]);

  time_kaon = path_kaon/vK;

  Double_t value = 0.;

  path_gamma = sqrt(pow(pp[20] - neu_vtx[0],2) + pow(pp[21] - neu_vtx[1],2) + pow(pp[22] - neu_vtx[2],2));
  time_gamma = path_gamma/PhysicsConstants::cVel;

  time_diff = time_gamma + time_kaon - pp[23];

  value = time_diff;

  return value;

}

Double_t kin_fits::cluster_time_cons_fourth(const Double_t *xx, const Double_t *pp)
{
  Float_t Phi4mom[4] = {pp[0], pp[1], pp[2], pp[3]};
  Float_t mom_ch[2][4], mom_neu[4][4];
  Float_t neu_vtx[3] = {pp[30], pp[31], pp[32]};
  Float_t ip_vtx[3] = {pp[33], pp[34], pp[35]};

  charged_mom(pp[4], pp[5], pp[6], mom_ch[0], 1);
  charged_mom(pp[7], pp[8], pp[9], mom_ch[1], 1);

  neutral_mom(pp[10], pp[11], pp[12], pp[14], neu_vtx, mom_neu[0]);
  neutral_mom(pp[15], pp[16], pp[17], pp[19], neu_vtx, mom_neu[1]);
  neutral_mom(pp[20], pp[21], pp[22], pp[24], neu_vtx, mom_neu[2]);
  neutral_mom(pp[25], pp[26], pp[27], pp[29], neu_vtx, mom_neu[3]);

  Double_t path_gamma, time_gamma, path_kaon, vK, time_kaon, time_diff;

  path_kaon = sqrt(pow(neu_vtx[0] - ip_vtx[0],2) + pow(neu_vtx[1] - ip_vtx[1],2) + pow(neu_vtx[2] - ip_vtx[2],2));
  vK = PhysicsConstants::cVel*sqrt(pow(mom_neu[0][0] + mom_neu[1][0] + mom_neu[2][0] + mom_neu[3][0], 2) + 
                  pow(mom_neu[0][1] + mom_neu[1][1] + mom_neu[2][1] + mom_neu[3][1], 2) +
                  pow(mom_neu[0][2] + mom_neu[1][2] + mom_neu[2][2] + mom_neu[3][2], 2))/(mom_neu[0][3] + mom_neu[1][3] + mom_neu[2][3] + mom_neu[3][3]);

  time_kaon = path_kaon/vK;

  Double_t value = 0.;

  path_gamma = sqrt(pow(pp[25] - neu_vtx[0],2) + pow(pp[26] - neu_vtx[1],2) + pow(pp[27] - neu_vtx[2],2));
  time_gamma = path_gamma/PhysicsConstants::cVel;

  time_diff = time_gamma + time_kaon - pp[28];

  value = time_diff;

  return value;

}

Double_t kin_fits::minv_neu_cons(const Double_t *xx, const Double_t *pp)
{
  Float_t Phi4mom[4] = {pp[0], pp[1], pp[2], pp[3]};
  Float_t mom_ch[2][4], mom_neu[4][4];
  Float_t neu_vtx[3] = {pp[30], pp[31], pp[32]};

  charged_mom(pp[4], pp[5], pp[6], mom_ch[0], 1);
  charged_mom(pp[7], pp[8], pp[9], mom_ch[1], 1);

  neutral_mom(pp[10], pp[11], pp[12], pp[14], neu_vtx, mom_neu[0]);
  neutral_mom(pp[15], pp[16], pp[17], pp[19], neu_vtx, mom_neu[1]);
  neutral_mom(pp[20], pp[21], pp[22], pp[24], neu_vtx, mom_neu[2]);
  neutral_mom(pp[25], pp[26], pp[27], pp[29], neu_vtx, mom_neu[3]);

  Double_t value = 0.;

  value = sqrt(pow(mom_neu[0][3] + mom_neu[1][3] + mom_neu[2][3] + mom_neu[3][3],2) - 
               pow(mom_neu[0][0] + mom_neu[1][0] + mom_neu[2][0] + mom_neu[3][0],2) -
               pow(mom_neu[0][1] + mom_neu[1][1] + mom_neu[2][1] + mom_neu[3][1],2) -
               pow(mom_neu[0][2] + mom_neu[1][2] + mom_neu[2][2] + mom_neu[3][2],2)) - PhysicsConstants::mK0;

  return value;

}

Double_t kin_fits::minv_ch_cons(const Double_t *xx, const Double_t *pp)
{
  Float_t Phi4mom[4] = {pp[0], pp[1], pp[2], pp[3]};
  Float_t mom_ch[2][4], mom_neu[4][4];
  Float_t neu_vtx[3] = {pp[30], pp[31], pp[32]};

  charged_mom(pp[4], pp[5], pp[6], mom_ch[0], 1);
  charged_mom(pp[7], pp[8], pp[9], mom_ch[1], 1);

  neutral_mom(pp[10], pp[11], pp[12], pp[14], neu_vtx, mom_neu[0]);
  neutral_mom(pp[15], pp[16], pp[17], pp[19], neu_vtx, mom_neu[1]);
  neutral_mom(pp[20], pp[21], pp[22], pp[24], neu_vtx, mom_neu[2]);
  neutral_mom(pp[25], pp[26], pp[27], pp[29], neu_vtx, mom_neu[3]);

  Double_t value = 0.;

  value = sqrt(pow(mom_ch[0][3] + mom_ch[1][3], 2) - 
               pow(mom_ch[0][0] + mom_ch[1][0], 2) -
               pow(mom_ch[0][1] + mom_ch[1][1], 2) -
               pow(mom_ch[0][2] + mom_ch[1][2], 2)) - PhysicsConstants::mK0;

  return value;

}

// Trilateration conditions

Double_t kin_fits::tri_kaon_path(const Double_t *xx, const Double_t *pp)
{
  Reconstructor R;
  Solution S;

  //! Trilateration for every iteration
	for (Int_t i = 0; i < 4; i++)
		R.SetClu(i, clusters[0][i],
						    clusters[1][i],
						    clusters[2][i],
						    clusters[3][i],
						    clusters[4][i]);

	R.SetClu(4, 0., 0., 0., 0., 0.);
	R.SetClu(5, 0., 0., 0., 0., 0.);

	S = R.MySolve(selected);

	for (Int_t i = 0; i < 2; i++)
		for (Int_t j = 0; j < 4; j++)
			neu_vtx[i][j] = S.sol[i][j];
	//!

  Float_t Phi_vtx[3] = {pp[20], pp[21], pp[22]};
  Float_t mom_neu[4][4];

  Float_t path_length, vK;

  neutral_mom(pp[0], pp[1], pp[2], pp[4], neu_vtx_tri, mom_neu[0]);
  neutral_mom(pp[5], pp[6], pp[7], pp[9], neu_vtx_tri, mom_neu[1]);
  neutral_mom(pp[10], pp[11], pp[12], pp[14], neu_vtx_tri, mom_neu[2]);
  neutral_mom(pp[15], pp[16], pp[17], pp[19], neu_vtx_tri, mom_neu[3]);

  vK = PhysicsConstants::cVel * sqrt(pow(mom_neu[0][0] + mom_neu[1][0] + mom_neu[2][0] + mom_neu[3][0],2) + pow(mom_neu[0][1] + mom_neu[1][1] + mom_neu[2][1] + mom_neu[3][1],2) + pow(mom_neu[0][2] + mom_neu[1][2] + mom_neu[2][2] + mom_neu[3][2],2))/(mom_neu[0][3] + mom_neu[1][3] + mom_neu[2][3] + mom_neu[3][3]);

  path_length = sqrt(pow(neu_vtx_tri[0] - Phi_vtx[0],2) + pow(neu_vtx_tri[1] - Phi_vtx[1],2) + pow(neu_vtx_tri[2] - Phi_vtx[2],2));

  Double_t value = 0.;

  value = (vK * neu_vtx_tri[3] - path_length);

  return value;

}