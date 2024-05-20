#include <TMath.h>

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "Math/Functor.h"
#include "Math/Factory.h"
#include "Math/Minimizer.h"

#include "../../../Include/Codes/reconstructor.h"
#include "../../../Include/const.h"
#include "../../../Include/Codes/uncertainties.h"
#include "../../../Include/Codes/charged_mom.h"
#include "../../../Include/Codes/neutral_mom.h"
#include "../../../Include/Codes/lorentz_transf.h"
#include "../../../Include/Codes/plane_intersection.h"
#include "../../../Include/Codes/closest_approach.h"
#include "../../../Include/Codes/kloe_class.h"
#include "../../../Include/Codes/chi2_dist.h"
#include "../inc/trilateration.hpp"

Double_t trilateration_chi_square(const Double_t *x)
{
  const Int_t N = 27, M = 9;
  Int_t selected[4] = {1, 2, 3, 4};

  Reconstructor R;
  Solution S;

  Float_t T0 = 0.;

	Float_t clusters[5][4], clusters_meas[5][4], clusters_err[5][4];
	Float_t kaon_velocity[2][3], kaon_path[2][3], kaon_inv_mass[2], kaon_mom[2], kaon_mom_vec_lor[2][4], kaon_path_tot[2], kaon_velocity_tot[2];

	Float_t boost_vec[3], bhabha_vtx[3], phi_mom[4], phi_vtx[3], y_axis[3] = {0., 1., 0.}, kaon_mom_vec[2][4], ip_rec[2][3];
	Float_t gamma_mom[2][4][4], neu_vtx[2][4], lambda[M], neu_vtx_one[4], gamma_mom_one[4][4], gamma_path[2][4], time_diff[2], test[3][3];

	//! Momenta of gammas reconstructed w/o the coordinates of Phi vtx

	for (Int_t i = 0; i < M; i++)
		lambda[i] = x[3 * N + i];

	bhabha_vtx[0] = x[24];
	bhabha_vtx[1] = x[25];
	bhabha_vtx[2] = x[26];

	phi_mom[0] = x[20];
	phi_mom[1] = x[21];
	phi_mom[2] = x[22];
	phi_mom[3] = x[23];

	boost_vec[0] = -phi_mom[0] / phi_mom[3];
	boost_vec[1] = -phi_mom[1] / phi_mom[3];
	boost_vec[2] = -phi_mom[2] / phi_mom[3];

	for (Int_t i = 0; i < 4; i++)
		for (Int_t j = 0; j < 5; j++)
		{
			clusters[j][i] = x[i * 5 + j];
		}

	//! Trilateration for every iteration
	for (Int_t i = 0; i < 4; i++)
		R.SetClu(i, clusters[0][i],
						 clusters[1][i],
						 clusters[2][i],
						 clusters[3][i] - T0,
						 clusters[4][i]);

	R.SetClu(4, 0., 0., 0., 0., 0.);
	R.SetClu(5, 0., 0., 0., 0., 0.);

	S = R.MySolve(selected);

	for (Int_t i = 0; i < 2; i++)
		for (Int_t j = 0; j < 4; j++)
			neu_vtx[i][j] = S.sol[i][j];
	//!

	//! Parameters for function builder
	Double_t value[2] = {0.}, chi2[2] = {0.}, constraints[2][M], min_value;
	//!

	//! Values without constraints
	for (Int_t j = 0; j < N; j++)
	{
		value[0] += pow((x[j] - x[j + N]) / x[j + 2 * N], 2);
	}
	value[1] = value[0];
	//!

	chi2[0] = value[0];
	chi2[1] = value[1];

	//! Gamma 4-momentum reconstruction
	for (Int_t i = 0; i < 2; i++)
		for (Int_t j = 0; j < 4; j++)
		{
			neutral_mom(clusters[0][j], clusters[1][j], clusters[2][j], clusters[4][j], neu_vtx[i], gamma_mom[i][j]);
			gamma_path[i][j] = sqrt(pow(clusters[0][j] - neu_vtx[i][0], 2) + pow(clusters[1][j] - neu_vtx[i][1], 2) + pow(clusters[2][j] - neu_vtx[i][2], 2));
		}
	//!

	for (Int_t i = 0; i < 2; i++)
	{
		kaon_mom_vec[i][0] = gamma_mom[i][0][0] + gamma_mom[i][1][0] + gamma_mom[i][2][0] + gamma_mom[i][3][0];
		kaon_mom_vec[i][1] = gamma_mom[i][0][1] + gamma_mom[i][1][1] + gamma_mom[i][2][1] + gamma_mom[i][3][1];
		kaon_mom_vec[i][2] = gamma_mom[i][0][2] + gamma_mom[i][1][2] + gamma_mom[i][2][2] + gamma_mom[i][3][2];
		kaon_mom_vec[i][3] = gamma_mom[i][0][3] + gamma_mom[i][1][3] + gamma_mom[i][2][3] + gamma_mom[i][3][3];

		kaon_velocity[i][0] = c_vel * kaon_mom_vec[i][0] / kaon_mom_vec[i][3];
		kaon_velocity[i][1] = c_vel * kaon_mom_vec[i][1] / kaon_mom_vec[i][3];
		kaon_velocity[i][2] = c_vel * kaon_mom_vec[i][2] / kaon_mom_vec[i][3];

		kaon_mom[i] = sqrt(pow(kaon_mom_vec[i][0], 2) + pow(kaon_mom_vec[i][1], 2) + pow(kaon_mom_vec[i][2], 2));

		kaon_inv_mass[i] = pow(kaon_mom_vec[i][3], 2) - pow(kaon_mom[i], 2);

		plane_intersection(bhabha_vtx, y_axis, neu_vtx[i], kaon_mom_vec[i], ip_rec[i]); //! Plane rec
		// closest_approach(bhabha_vtx, z_axis, neu_vtx[i], kaon_mom_vec[i], ip_rec[i]); //! Closest approach

    ip_rec[i][0] = bhabha_vtx[0];
    ip_rec[i][1] = bhabha_vtx[1];

		if (abs(ip_rec[i][2] - bhabha_vtx[2]) > 2)
		{
			ip_rec[i][2] = bhabha_vtx[2];
		}

		kaon_path[i][0] = neu_vtx[i][0] - ip_rec[i][0];
		kaon_path[i][1] = neu_vtx[i][1] - ip_rec[i][1];
		kaon_path[i][2] = neu_vtx[i][2] - ip_rec[i][2];

		lorentz_transf(boost_vec, kaon_mom_vec[i], kaon_mom_vec_lor[i]); //! Lorentz transformation

    constraints[i][0] = pow(kaon_velocity[i][0] * neu_vtx[i][3] - kaon_path[i][0], 2);
		constraints[i][1] = pow(kaon_velocity[i][1] * neu_vtx[i][3] - kaon_path[i][1], 2);
		constraints[i][2] = pow(kaon_velocity[i][2] * neu_vtx[i][3] - kaon_path[i][2], 2);

		constraints[i][3] = pow(pow(kaon_inv_mass[i],2) - pow(m_k0,2), 2);

		constraints[i][4] = pow(kaon_mom_vec_lor[i][3] - (phi_mom[3] / 2.), 2);

		constraints[i][5] = pow(clusters[3][0] - neu_vtx[i][3] - (gamma_path[i][0] / c_vel),2);
		constraints[i][6] = pow(clusters[3][1] - neu_vtx[i][3] - (gamma_path[i][1] / c_vel),2);
		constraints[i][7] = pow(clusters[3][2] - neu_vtx[i][3] - (gamma_path[i][2] / c_vel),2);
		constraints[i][8] = pow(clusters[3][3] - neu_vtx[i][3] - (gamma_path[i][3] / c_vel),2);

		value[i] += lambda[0] * constraints[i][0] + // X-axis path
								lambda[1] * constraints[i][1] + // Y-axis path
								lambda[2] * constraints[i][2] + // Z-axis path
								lambda[3] * constraints[i][3] + // Inv kaon mass
								lambda[4] * constraints[i][4] + // Half of phi meson's energy
								lambda[5] * constraints[i][5] + // 1st tri eqn
								lambda[6] * constraints[i][6] + // 2nd tri eqn
								lambda[7] * constraints[i][7] + // 3rd tri eqn
								lambda[8] * constraints[i][8];  // 4th tri eqn
	}

	if (!TMath::IsNaN(value[0]) && !TMath::IsNaN(value[1]))
	{
		if (value[0] < value[1])
		{
			min_value = value[0];
		}
		else if (value[1] < value[0])
		{
			min_value = value[1];
		}
	}
	else if (TMath::IsNaN(value[0]) && !TMath::IsNaN(value[1]))
	{
		min_value = value[0];
	}
	else if (!TMath::IsNaN(value[0]) && TMath::IsNaN(value[1]))
	{
		min_value = value[1];
	}
	else
	{
		min_value = 999999.;
	}

  std::cout << min_value << " " << min_value -  chi2[0] << " " << chi2[0]<< std::endl;

	return min_value;
}