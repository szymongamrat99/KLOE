#include <string.h>

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "Math/Functor.h"
#include "Math/Factory.h"
#include "Math/Minimizer.h"

#include "../../Include/Codes/reconstructor.h"
#include "../../Include/const.h"
#include "../../Include/Codes/uncertainties.h"
#include "../../Include/Codes/charged_mom.h"
#include "../../Include/Codes/neutral_mom.h"
#include "../../Include/Codes/lorentz_transf.h"
#include "../../Include/Codes/plane_intersection.h"
#include "../../Include/Codes/closest_approach.h"
#include "../../Include/Codes/kloe_class.h"
#include "../../Include/Codes/kinematic_fits.h"
#include "../../Include/Codes/chi2_dist.h"
#include "chain_init.C"

const Int_t N = 27, M = 9;
const Float_t Trf = 2.715; // ns - time of a bunch (correction)
Float_t T0[2], T0_mean;

Reconstructor R;
Solution S;

Int_t selected[4] = {1, 2, 3, 4};

Double_t trilateration_chi_square(const Double_t *x)
{
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
						 clusters[3][i],
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

		kaon_inv_mass[i] = sqrt(pow(kaon_mom_vec[i][3], 2) - pow(kaon_mom[i], 2));

		plane_intersection(bhabha_vtx, y_axis, neu_vtx[i], kaon_mom_vec[i], ip_rec[i]); //! Plane rec
		// closest_approach(bhabha_vtx, z_axis, neu_vtx[i], kaon_mom_vec[i], ip_rec[i]); //! Plane rec

		if (abs(ip_rec[i][2] - bhabha_vtx[2]) > 2)
		{
			ip_rec[i][2] = bhabha_vtx[2];
		}

		ip_rec[i][0] = bhabha_vtx[0];

		kaon_path[i][0] = neu_vtx[i][0] - ip_rec[i][0];
		kaon_path[i][1] = neu_vtx[i][1] - ip_rec[i][1];
		kaon_path[i][2] = neu_vtx[i][2] - ip_rec[i][2];

		kaon_path_tot[i] = sqrt(pow(kaon_path[i][0], 2) + pow(kaon_path[i][1], 2) + pow(kaon_path[i][2], 2));
		kaon_velocity_tot[i] = sqrt(pow(kaon_velocity[i][0], 2) + pow(kaon_velocity[i][1], 2) + pow(kaon_velocity[i][2], 2));

		lorentz_transf(boost_vec, kaon_mom_vec[i], kaon_mom_vec_lor[i]); //! Lorentz transformation

		constraints[i][3] = pow(kaon_inv_mass[i] - m_k0, 2);

		// std::cout << kaon_mom_vec_lor[i][3] << std::endl;

		constraints[i][4] = pow(kaon_mom_vec_lor[i][3] - (m_phi / 2.), 2);

		constraints[i][5] = clusters[3][0] - neu_vtx[i][3] - (gamma_path[i][0] / c_vel);
		constraints[i][6] = clusters[3][1] - neu_vtx[i][3] - (gamma_path[i][1] / c_vel);
		constraints[i][7] = clusters[3][2] - neu_vtx[i][3] - (gamma_path[i][2] / c_vel);
		constraints[i][8] = clusters[3][3] - neu_vtx[i][3] - (gamma_path[i][3] / c_vel);

		T0[i] = 0;//TMath::Nint(constraints[i][5]/Trf)*Trf;

		//std::cout << constraints[i][5]/Trf << std::endl;

		constraints[i][5] = pow(constraints[i][5], 2);
		constraints[i][6] = pow(constraints[i][6], 2);
		constraints[i][7] = pow(constraints[i][7], 2);
		constraints[i][8] = pow(constraints[i][8], 2);

		// std::cout << T0[i] << std::endl;

		constraints[i][0] = pow(kaon_velocity_tot[i] * (neu_vtx[i][3] + T0[i]) - kaon_path_tot[i], 2);
		constraints[i][1] = pow(kaon_velocity[i][1] * (neu_vtx[i][3] + T0[i]) - kaon_path[i][1], 2);
		constraints[i][2] = pow(kaon_velocity[i][2] * (neu_vtx[i][3] + T0[i]) - kaon_path[i][2], 2);

		value[i] += lambda[0] * constraints[i][0] +
								/*lambda[1] * constraints[i][1]
								lambda[2] * constraints[i][2] +*/
								lambda[3] * constraints[i][3] +
								//lambda[4] * constraints[i][4]
								lambda[5] * constraints[i][5] +
								lambda[6] * constraints[i][6] +
								lambda[7] * constraints[i][7] +
								lambda[8] * constraints[i][8];
	}

	// std::cout << std::endl;

	if (!TMath::IsNaN(value[0]) && !TMath::IsNaN(value[1]))
	{
		if (value[0] < value[1])
		{
			min_value = value[0];
			// std::cout << chi2[0] << " " << min_value << " " << constraints[0][0] << " " << constraints[0][1] << " " << constraints[0][2] << " " << constraints[0][3] << " " << constraints[0][4] << std::endl;
		}
		else if (value[1] < value[0])
		{
			// std::cout << chi2[1] << " " << min_value << " " << constraints[1][0] << " " << constraints[1][1] << " " << constraints[1][2] << " " << constraints[1][3] << " " << constraints[1][4] << std::endl;
			min_value = value[1];
		}
	}
	else if (TMath::IsNaN(value[0]) && !TMath::IsNaN(value[1]))
	{
		min_value = value[0];
		// std::cout << chi2[0] << " " << min_value << " " << constraints[0][0] << " " << constraints[0][1] << " " << constraints[0][2] << " " << constraints[0][3] << " " << constraints[0][4] << std::endl;
	}
	else if (!TMath::IsNaN(value[0]) && TMath::IsNaN(value[1]))
	{
		min_value = value[1];
		// std::cout << chi2[1] << " " << min_value << " " << constraints[1][0] << " " << constraints[1][1] << " " << constraints[1][2] << " " << constraints[1][3] << " " << constraints[1][4] << std::endl;
	}
	else
	{
		min_value = 999999.;
	}

	// std::cout << std::endl;

	return min_value;
}

int main(int argc, char *argv[]) //	arguments are: 1. Number of points
																 //								 2. First analyzed file
																 //                3. Last analyzed file
{
	UInt_t first = atoi(argv[2]), last = atoi(argv[3]);

	TString first_s = argv[2], last_s = argv[3];

	TChain *chain = new TChain("INTERF/h1");
	chain_init(chain, first, last);

	TString name = "neuvtx_tri_kin_fit_" + first_s + "_" + last_s + ".root";

	TFile *file = new TFile(name, "recreate");
	TTree *tree = new TTree("h_tri_kin_fit", "Neu vtx rec with trilateration kin fit");

	// Branches' addresses
	// Bhabha vars
	Float_t bhabha_mom[4], bhabha_mom_err[4], bhabha_vtx[3], bhabha_vtx_err[3];

	chain->SetBranchAddress("Bpx", &bhabha_mom[0]);
	chain->SetBranchAddress("Bpy", &bhabha_mom[1]);
	chain->SetBranchAddress("Bpz", &bhabha_mom[2]);
	chain->SetBranchAddress("Broots", &bhabha_mom[3]);

	chain->SetBranchAddress("Bwidpx", &bhabha_mom_err[0]);
	chain->SetBranchAddress("Bwidpy", &bhabha_mom_err[1]);
	chain->SetBranchAddress("Bwidpz", &bhabha_mom_err[2]);
	chain->SetBranchAddress("Brootserr", &bhabha_mom_err[3]);

	chain->SetBranchAddress("Bx", &bhabha_vtx[0]);
	chain->SetBranchAddress("By", &bhabha_vtx[1]);
	chain->SetBranchAddress("Bz", &bhabha_vtx[2]);
	chain->SetBranchAddress("Bsx", &bhabha_vtx_err[0]);
	chain->SetBranchAddress("Bsy", &bhabha_vtx_err[1]);
	chain->SetBranchAddress("Bsz", &bhabha_vtx_err[2]);

	// Cluster vars
	Int_t nclu;
	UChar_t mctruth;
	Float_t cluster[5][500], Kchboost[9], Knerec[9], Knemc[9], ipmc[3], ip[3], Dtmc;

	chain->SetBranchAddress("nclu", &nclu);
	chain->SetBranchAddress("Xcl", cluster[0]);
	chain->SetBranchAddress("Ycl", cluster[1]);
	chain->SetBranchAddress("Zcl", cluster[2]);
	chain->SetBranchAddress("Tcl", cluster[3]);
	chain->SetBranchAddress("Enecl", cluster[4]);

	chain->SetBranchAddress("mctruth", &mctruth);

	chain->SetBranchAddress("Knemc", Knemc);
	chain->SetBranchAddress("Knerec", Knerec);

	chain->SetBranchAddress("Kchboost", Kchboost);
	chain->SetBranchAddress("Dtmc", &Dtmc);

	chain->SetBranchAddress("ipmc", ipmc);
	chain->SetBranchAddress("ip", ip);

	//! Parameters for fitting

	Double_t P[N], DP[N];

	Double_t CHISQR;
	Float_t CHISQRMIN;



	//!

	Int_t nentries = (Int_t)chain->GetEntries();

	Bool_t clusterEnergy, solError;
	Int_t ind_gam[4], sort_index[2][3], sort_ind_gam[2][4], chosen_ind_gam[4], found_best, isConverged;
	Float_t neu_vtx[2][4], inv_m_pi0[2][3][2], gamma_mom[2][4][4], mass_pair[2][3];

	Float_t P1[3 * N + M + 3], min_value, min_value_def, neu_vtx_min[2][4], clusters_min[4][5], neu_vtx_min_final[4],
			mom_kaon[2][4], ene_kaon[2], v_kaon[2], length_kaon[2], diff_kaon[2], length, length_mc, length_rec;

	Float_t length_ch, time_ch, velocity_ch, gamma_mom_final[4][8], fourKnetri_kinfit[10], constraints[M], iptri_kinfit[3], y_axis[3] = {0., 1., 0.}, ip_tri[2][3], bhabha_mom_fit[4];
	Int_t g4takentri_kinfit[4];

	TBranch *b_gamma1tri = tree->Branch("fourgamma1tri_kinfit", gamma_mom_final[0], "fourgamma1tri_kinfit[8]/F");
	TBranch *b_gamma2tri = tree->Branch("fourgamma2tri_kinfit", gamma_mom_final[1], "fourgamma2tri_kinfit[8]/F");
	TBranch *b_gamma3tri = tree->Branch("fourgamma3tri_kinfit", gamma_mom_final[2], "fourgamma3tri_kinfit[8]/F");
	TBranch *b_gamma4tri = tree->Branch("fourgamma4tri_kinfit", gamma_mom_final[3], "fourgamma4tri_kinfit[8]/F");
	TBranch *b_iptri = tree->Branch("iptri_kinfit", iptri_kinfit, "iptri_kinfit[3]/F");
	TBranch *b_Knetri = tree->Branch("fourKnetri_kinfit", fourKnetri_kinfit, "fourKnetri_kinfit[10]/F");
	TBranch *b_done = tree->Branch("done4_kinfit", &found_best, "done4_kinfit/I");
	TBranch *b_fourg4taken = tree->Branch("g4takentri_kinfit", g4takentri_kinfit, "g4takentri_kinfit[4]/I");

	TBranch *b_chisqr = tree->Branch("chi2min", &CHISQRMIN, "chi2min/F");

	UInt_t number_of_ev = atoi(argv[1]);
	TString number_of_ev_s = argv[1];

	if (number_of_ev_s == "all")
		number_of_ev = nentries;
	else if (number_of_ev > nentries)
		number_of_ev = nentries;

	TH1 *chi2_hist = new TH1F("chi2", "", 100, 0, 100);
	TH2 *chi2_corr = new TH2F("chi2_corr", "", N, 1, N, N, 1, N);

	for (Int_t i = 0; i < number_of_ev; i++)
	{
		chain->GetEntry(i);

		min_value_def = 99999.;
		CHISQRMIN = 99999.;
		found_best = 0;

		if (nclu >= 4 && (mctruth == 1 || mctruth == 2))
		{

			for (Int_t j1 = 0; j1 < nclu - 3; j1++)
				for (Int_t j2 = j1 + 1; j2 < nclu - 2; j2++)
					for (Int_t j3 = j2 + 1; j3 < nclu - 1; j3++)
						for (Int_t j4 = j3 + 1; j4 < nclu; j4++)
						{
							CHISQR = 0.;

							ind_gam[0] = j1;
							ind_gam[1] = j2;
							ind_gam[2] = j3;
							ind_gam[3] = j4;

							// Setting clusters for a solution
							for (Int_t l = 0; l < 4; l++)
								R.SetClu(l, cluster[0][ind_gam[l]],
												 cluster[1][ind_gam[l]],
												 cluster[2][ind_gam[l]],
												 cluster[3][ind_gam[l]],
												 cluster[4][ind_gam[l]]);

							R.SetClu(4, 0., 0., 0., 0., 0.);
							R.SetClu(5, 0., 0., 0., 0., 0.);

							S = R.MySolve(selected);

							clusterEnergy = (cluster[4][ind_gam[0]] >= 20 && cluster[4][ind_gam[1]] >= 20 && cluster[4][ind_gam[2]] >= 20 && cluster[4][ind_gam[3]] >= 20);

							solError = S.error[0] == false && S.error[1] == false;

							if (clusterEnergy && solError)
							{
								// !
								// !ITERATING OVER SOLUTIONS OF THE RECONSTRUCTION

								ROOT::Math::Minimizer *minimum;

								minimum = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");

								// set tolerance , etc...
								minimum->SetMaxFunctionCalls(10000000); // for Minuit/Minuit2
								minimum->SetTolerance(0.01);
								minimum->SetPrintLevel(0);

								// create function wrapper for minimizer
								// a IMultiGenFunction type
								ROOT::Math::Functor f(&trilateration_chi_square, 3 * N + M);

								minimum->SetFunction(f);

								for (Int_t l = 0; l < 3 * N + M; l++)
								{
									minimum->SetVariable(l, std::to_string(l), 0., 0.001);
									if (l >= N)
										minimum->FixVariable(l);
								}

								// Building the parameters for a fit
								for (Int_t k = 0; k < 4; k++)
								{
									P[k * 5] = cluster[0][ind_gam[k]];		 // X_k
									P[k * 5 + 1] = cluster[1][ind_gam[k]]; // Y_k
									P[k * 5 + 2] = cluster[2][ind_gam[k]]; // Z_k
									P[k * 5 + 3] = cluster[3][ind_gam[k]]; // T_k
									P[k * 5 + 4] = cluster[4][ind_gam[k]]; // E_k

									DP[k * 5] = 1.2;		 // cm
									DP[k * 5 + 1] = 1.2; // cm
									DP[k * 5 + 2] = 1.2 / sqrt(cluster[4][ind_gam[k]]/1000.); // cm
									DP[k * 5 + 3] = clu_time_error(cluster[4][ind_gam[k]]);
									DP[k * 5 + 4] = clu_ene_error(cluster[4][ind_gam[k]]);
								}

								P[20] = bhabha_mom[0]; // Px_phi
								P[21] = bhabha_mom[1]; // Py_phi
								P[22] = bhabha_mom[2]; // Pz_phi
								P[23] = bhabha_mom[3]; // Pz_phi

								DP[20] = bhabha_mom_err[0];
								DP[21] = bhabha_mom_err[1];
								DP[22] = bhabha_mom_err[2];
								DP[23] = bhabha_mom_err[3];

								P[24] = bhabha_vtx[0]; // Px_phi
								P[25] = bhabha_vtx[1]; // Py_phi
								P[26] = bhabha_vtx[2]; // Pz_phi

								DP[24] = bhabha_vtx_err[0];
								DP[25] = bhabha_vtx_err[1];
								DP[26] = bhabha_vtx_err[2];

								for (Int_t k = 0; k < N; k++)
								{
									minimum->SetVariableValue(k, P[k]);
									minimum->SetVariableLimits(k, P[k] - 5 * DP[k], P[k] + 5 * DP[k]);
									minimum->SetVariableValue(k + N, P[k]);
									minimum->SetVariableValue(k + 2 * N, DP[k]);
								}

								for (Int_t k = 0; k < M; k++)
								{
									minimum->SetVariableValue(k + 3 * N, 1.0);
								}

								// minimum->SetVariableStepSize(20,0.);
								// minimum->SetVariableStepSize(21,0.);
								// minimum->SetVariableStepSize(22,0.);
								// minimum->SetVariableStepSize(23,0.);
								minimum->SetVariableStepSize(24,0.);
								minimum->SetVariableStepSize(25,0.);
								minimum->SetVariableStepSize(26,0.);
								// minimum->SetVariableStepSize(3*N,0.1);
								// minimum->SetVariableStepSize(3*N + 1,0.1);
								// minimum->SetVariableStepSize(3*N + 2,0.1);
								// minimum->SetVariableStepSize(3*N + 3,0.1);
								// minimum->SetVariableStepSize(3*N + 4,0.1);
								// minimum->SetVariableStepSize(3*N + 5,0.1);
								// minimum->SetVariableStepSize(3*N + 6,0.1);
								// minimum->SetVariableStepSize(3*N + 7,0.1);
								// minimum->SetVariableStepSize(3*N + 8,0.1);
								// do the minimization
								isConverged = minimum->Minimize();

								if (isConverged == 1)
								{
									min_value = minimum->MinValue();
									for (Int_t m = 0; m < N; m++)
									{
										P1[m] = minimum->X()[m];
										CHISQR += pow((P1[m] - P[m]) / DP[m], 2);
										for (Int_t n = 0; n < N; n++)
											chi2_corr->SetBinContent(m+1,n+1,minimum->Correlation(m,n));
									}
								}
								else
								{
									for (Int_t m = 0; m < N; m++)
									{
										CHISQR = 999999;
										min_value = 999999;
										P1[m] = -999.;
									}
								}

								if (min_value < min_value_def)
								{
									min_value_def = min_value;
									CHISQRMIN = CHISQR;
									found_best = 1;

									// Setting clusters for a solution
									for (Int_t l = 0; l < 4; l++)
									{
										clusters_min[l][0] = P1[l * 5];
										clusters_min[l][1] = P1[l * 5 + 1];
										clusters_min[l][2] = P1[l * 5 + 2];
										clusters_min[l][3] = P1[l * 5 + 3];
										clusters_min[l][4] = P1[l * 5 + 4];

										R.SetClu(l, clusters_min[l][0],
														 clusters_min[l][1],
														 clusters_min[l][2],
														 clusters_min[l][3],
														 clusters_min[l][4]);
									}

									R.SetClu(4, 0., 0., 0., 0., 0.);
									R.SetClu(5, 0., 0., 0., 0., 0.);

									S = R.MySolve(selected);

									neu_vtx_min[0][0] = S.sol[0][0];
									neu_vtx_min[0][1] = S.sol[0][1];
									neu_vtx_min[0][2] = S.sol[0][2];
									neu_vtx_min[0][3] = S.sol[0][3];

									neu_vtx_min[1][0] = S.sol[1][0];
									neu_vtx_min[1][1] = S.sol[1][1];
									neu_vtx_min[1][2] = S.sol[1][2];
									neu_vtx_min[1][3] = S.sol[1][3];

									for (Int_t k = 0; k < 2; k++)
									{
										for (Int_t l = 0; l < 4; l++)
										{
											neutral_mom(clusters_min[l][0], clusters_min[l][1], clusters_min[l][2], clusters_min[l][4], neu_vtx_min[k], gamma_mom[k][l]);
										}

										mom_kaon[k][0] = gamma_mom[k][0][0] + gamma_mom[k][1][0] + gamma_mom[k][2][0] + gamma_mom[k][3][0];
										mom_kaon[k][1] = gamma_mom[k][0][1] + gamma_mom[k][1][1] + gamma_mom[k][2][1] + gamma_mom[k][3][1];
										mom_kaon[k][2] = gamma_mom[k][0][2] + gamma_mom[k][1][2] + gamma_mom[k][2][2] + gamma_mom[k][3][2];
										mom_kaon[k][3] = gamma_mom[k][0][3] + gamma_mom[k][1][3] + gamma_mom[k][2][3] + gamma_mom[k][3][3];

										v_kaon[k] = c_vel * sqrt(pow(mom_kaon[k][0], 2) + pow(mom_kaon[k][1], 2) + pow(mom_kaon[k][2], 2)) / mom_kaon[k][3];

										bhabha_mom_fit[0] = P1[20];
										bhabha_mom_fit[1] = P1[21];
										bhabha_mom_fit[2] = P1[22];
										bhabha_mom_fit[3] = P1[23];

										plane_intersection(bhabha_vtx, y_axis, neu_vtx_min[k], mom_kaon[k], ip_tri[k]); //! Plane rec

										if (abs(ip_tri[k][2] - bhabha_vtx[2]) > 2)
											ip_tri[k][2] = bhabha_vtx[2];
										if (abs(ip_tri[k][0] - bhabha_vtx[0]) > 0.2)
											ip_tri[k][0] = bhabha_vtx[0];

										length_kaon[k] = sqrt(pow(neu_vtx_min[k][0] - ip_tri[k][0], 2) + pow(neu_vtx_min[k][1] - ip_tri[k][1], 2) + pow(neu_vtx_min[k][2] - ip_tri[k][2], 2));

										diff_kaon[k] = v_kaon[k] * neu_vtx_min[k][3] - length_kaon[k];
									}

									if (abs(diff_kaon[0]) < abs(diff_kaon[1]))
									{
										gamma_mom_final[0][0] = gamma_mom[0][0][0];
										gamma_mom_final[0][1] = gamma_mom[0][0][1];
										gamma_mom_final[0][2] = gamma_mom[0][0][2];
										gamma_mom_final[0][3] = gamma_mom[0][0][3];
										gamma_mom_final[0][4] = clusters_min[0][0];
										gamma_mom_final[0][5] = clusters_min[0][1];
										gamma_mom_final[0][6] = clusters_min[0][2];
										gamma_mom_final[0][7] = clusters_min[0][3];

										gamma_mom_final[1][0] = gamma_mom[0][1][0];
										gamma_mom_final[1][1] = gamma_mom[0][1][1];
										gamma_mom_final[1][2] = gamma_mom[0][1][2];
										gamma_mom_final[1][3] = gamma_mom[0][1][3];
										gamma_mom_final[1][4] = clusters_min[1][0];
										gamma_mom_final[1][5] = clusters_min[1][1];
										gamma_mom_final[1][6] = clusters_min[1][2];
										gamma_mom_final[1][7] = clusters_min[1][3];

										gamma_mom_final[2][0] = gamma_mom[0][2][0];
										gamma_mom_final[2][1] = gamma_mom[0][2][1];
										gamma_mom_final[2][2] = gamma_mom[0][2][2];
										gamma_mom_final[2][3] = gamma_mom[0][2][3];
										gamma_mom_final[2][4] = clusters_min[2][0];
										gamma_mom_final[2][5] = clusters_min[2][1];
										gamma_mom_final[2][6] = clusters_min[2][2];
										gamma_mom_final[2][7] = clusters_min[2][3];

										gamma_mom_final[3][0] = gamma_mom[0][3][0];
										gamma_mom_final[3][1] = gamma_mom[0][3][1];
										gamma_mom_final[3][2] = gamma_mom[0][3][2];
										gamma_mom_final[3][3] = gamma_mom[0][3][3];
										gamma_mom_final[3][4] = clusters_min[3][0];
										gamma_mom_final[3][5] = clusters_min[3][1];
										gamma_mom_final[3][6] = clusters_min[3][2];
										gamma_mom_final[3][7] = clusters_min[3][3];

										fourKnetri_kinfit[0] = gamma_mom_final[0][0] + gamma_mom_final[1][0] +
																					 gamma_mom_final[2][0] + gamma_mom_final[3][0];
										fourKnetri_kinfit[1] = gamma_mom_final[0][1] + gamma_mom_final[1][1] +
																					 gamma_mom_final[2][1] + gamma_mom_final[3][1];
										fourKnetri_kinfit[2] = gamma_mom_final[0][2] + gamma_mom_final[1][2] +
																					 gamma_mom_final[2][2] + gamma_mom_final[3][2];
										fourKnetri_kinfit[3] = gamma_mom_final[0][3] + gamma_mom_final[1][3] +
																					 gamma_mom_final[2][3] + gamma_mom_final[3][3];
										fourKnetri_kinfit[4] = sqrt(pow(fourKnetri_kinfit[0], 2) +
																								pow(fourKnetri_kinfit[1], 2) +
																								pow(fourKnetri_kinfit[2], 2));
										fourKnetri_kinfit[5] = sqrt(pow(fourKnetri_kinfit[3], 2) - pow(fourKnetri_kinfit[4], 2));
										fourKnetri_kinfit[6] = neu_vtx_min[0][0];
										fourKnetri_kinfit[7] = neu_vtx_min[0][1];
										fourKnetri_kinfit[8] = neu_vtx_min[0][2];
										fourKnetri_kinfit[9] = neu_vtx_min[0][3];

										g4takentri_kinfit[0] = ind_gam[0];
										g4takentri_kinfit[1] = ind_gam[1];
										g4takentri_kinfit[2] = ind_gam[2];
										g4takentri_kinfit[3] = ind_gam[3];

										plane_intersection(bhabha_vtx, y_axis, neu_vtx_min[0], fourKnetri_kinfit, iptri_kinfit); //! Plane rec

										if (abs(iptri_kinfit[2] - bhabha_vtx[2]) > 2)
											iptri_kinfit[2] = bhabha_vtx[2];
										if (abs(iptri_kinfit[0] - bhabha_vtx[0]) > 0.2)
											iptri_kinfit[0] = bhabha_vtx[0];
									}
									else if (abs(diff_kaon[0]) > abs(diff_kaon[1]))
									{
										gamma_mom_final[0][0] = gamma_mom[1][0][0];
										gamma_mom_final[0][1] = gamma_mom[1][0][1];
										gamma_mom_final[0][2] = gamma_mom[1][0][2];
										gamma_mom_final[0][3] = gamma_mom[1][0][3];
										gamma_mom_final[0][4] = clusters_min[0][0];
										gamma_mom_final[0][5] = clusters_min[0][1];
										gamma_mom_final[0][6] = clusters_min[0][2];
										gamma_mom_final[0][7] = clusters_min[0][3];

										gamma_mom_final[1][0] = gamma_mom[1][1][0];
										gamma_mom_final[1][1] = gamma_mom[1][1][1];
										gamma_mom_final[1][2] = gamma_mom[1][1][2];
										gamma_mom_final[1][3] = gamma_mom[1][1][3];
										gamma_mom_final[1][4] = clusters_min[1][0];
										gamma_mom_final[1][5] = clusters_min[1][1];
										gamma_mom_final[1][6] = clusters_min[1][2];
										gamma_mom_final[1][7] = clusters_min[1][3];

										gamma_mom_final[2][0] = gamma_mom[1][2][0];
										gamma_mom_final[2][1] = gamma_mom[1][2][1];
										gamma_mom_final[2][2] = gamma_mom[1][2][2];
										gamma_mom_final[2][3] = gamma_mom[1][2][3];
										gamma_mom_final[2][4] = clusters_min[2][0];
										gamma_mom_final[2][5] = clusters_min[2][1];
										gamma_mom_final[2][6] = clusters_min[2][2];
										gamma_mom_final[2][7] = clusters_min[2][3];

										gamma_mom_final[3][0] = gamma_mom[1][3][0];
										gamma_mom_final[3][1] = gamma_mom[1][3][1];
										gamma_mom_final[3][2] = gamma_mom[1][3][2];
										gamma_mom_final[3][3] = gamma_mom[1][3][3];
										gamma_mom_final[3][4] = clusters_min[3][0];
										gamma_mom_final[3][5] = clusters_min[3][1];
										gamma_mom_final[3][6] = clusters_min[3][2];
										gamma_mom_final[3][7] = clusters_min[3][3];

										fourKnetri_kinfit[0] = gamma_mom_final[0][0] + gamma_mom_final[1][0] +
																					 gamma_mom_final[2][0] + gamma_mom_final[3][0];
										fourKnetri_kinfit[1] = gamma_mom_final[0][1] + gamma_mom_final[1][1] +
																					 gamma_mom_final[2][1] + gamma_mom_final[3][1];
										fourKnetri_kinfit[2] = gamma_mom_final[0][2] + gamma_mom_final[1][2] +
																					 gamma_mom_final[2][2] + gamma_mom_final[3][2];
										fourKnetri_kinfit[3] = gamma_mom_final[0][3] + gamma_mom_final[1][3] +
																					 gamma_mom_final[2][3] + gamma_mom_final[3][3];
										fourKnetri_kinfit[4] = sqrt(pow(fourKnetri_kinfit[0], 2) +
																								pow(fourKnetri_kinfit[1], 2) +
																								pow(fourKnetri_kinfit[2], 2));
										fourKnetri_kinfit[5] = sqrt(pow(fourKnetri_kinfit[3], 2) - pow(fourKnetri_kinfit[4], 2));
										fourKnetri_kinfit[6] = neu_vtx_min[1][0];
										fourKnetri_kinfit[7] = neu_vtx_min[1][1];
										fourKnetri_kinfit[8] = neu_vtx_min[1][2];
										fourKnetri_kinfit[9] = neu_vtx_min[1][3];

										g4takentri_kinfit[0] = ind_gam[0];
										g4takentri_kinfit[1] = ind_gam[1];
										g4takentri_kinfit[2] = ind_gam[2];
										g4takentri_kinfit[3] = ind_gam[3];

										plane_intersection(bhabha_vtx, y_axis, neu_vtx_min[1], fourKnetri_kinfit, iptri_kinfit); //! Plane rec

										if (abs(iptri_kinfit[2] - bhabha_vtx[2]) > 2)
											iptri_kinfit[2] = bhabha_vtx[2];
										if (abs(iptri_kinfit[0] - bhabha_vtx[0]) > 0.2)
											iptri_kinfit[0] = bhabha_vtx[0];
									}
								}

								delete minimum;
							}
							///////////////////////////////////////////////////////////////////
						}

			if (found_best == 1)
			{
				chi2_hist->Fill(CHISQRMIN);

				std::cout << CHISQRMIN << std::endl;
			}
		}

		if (found_best == 0)
		{
			gamma_mom_final[0][0] = -999.;
			gamma_mom_final[0][1] = -999.;
			gamma_mom_final[0][2] = -999.;
			gamma_mom_final[0][3] = -999.;
			gamma_mom_final[0][4] = -999.;
			gamma_mom_final[0][5] = -999.;
			gamma_mom_final[0][6] = -999.;
			gamma_mom_final[0][7] = -999.;

			gamma_mom_final[1][0] = -999.;
			gamma_mom_final[1][1] = -999.;
			gamma_mom_final[1][2] = -999.;
			gamma_mom_final[1][3] = -999.;
			gamma_mom_final[1][4] = -999.;
			gamma_mom_final[1][5] = -999.;
			gamma_mom_final[1][6] = -999.;
			gamma_mom_final[1][7] = -999.;

			gamma_mom_final[2][0] = -999.;
			gamma_mom_final[2][1] = -999.;
			gamma_mom_final[2][2] = -999.;
			gamma_mom_final[2][3] = -999.;
			gamma_mom_final[2][4] = -999.;
			gamma_mom_final[2][5] = -999.;
			gamma_mom_final[2][6] = -999.;
			gamma_mom_final[2][7] = -999.;

			gamma_mom_final[3][0] = -999.;
			gamma_mom_final[3][1] = -999.;
			gamma_mom_final[3][2] = -999.;
			gamma_mom_final[3][3] = -999.;
			gamma_mom_final[3][4] = -999.;
			gamma_mom_final[3][5] = -999.;
			gamma_mom_final[3][6] = -999.;
			gamma_mom_final[3][7] = -999.;

			fourKnetri_kinfit[0] = -999.;
			fourKnetri_kinfit[1] = -999.;
			fourKnetri_kinfit[2] = -999.;
			fourKnetri_kinfit[3] = -999.;
			fourKnetri_kinfit[4] = -999.;
			fourKnetri_kinfit[5] = -999.;
			fourKnetri_kinfit[6] = -999.;
			fourKnetri_kinfit[7] = -999.;
			fourKnetri_kinfit[8] = -999.;
			fourKnetri_kinfit[9] = -999.;

			g4takentri_kinfit[0] = -999;
			g4takentri_kinfit[1] = -999;
			g4takentri_kinfit[2] = -999;
			g4takentri_kinfit[3] = -999;

			iptri_kinfit[0] = -999.;
			iptri_kinfit[1] = -999.;
			iptri_kinfit[2] = -999.;
		}

		tree->Fill();
	}

	TCanvas *c1 = new TCanvas("c1", "", 750, 750);

	name = "#chi^{2}(" + std::to_string(M-4) + ")";

	chi2_hist->GetXaxis()->SetTitle(name);
	chi2_hist->GetYaxis()->SetTitle("Counts");
	chi2_hist->GetYaxis()->SetRangeUser(0,1.2*chi2_hist->GetMaximum());
	chi2_hist->Draw();

	c1->Print("chi2_test.png");

	TCanvas *c2 = new TCanvas("c2", "", 750, 750);
	c1->SetRightMargin(0.15);
	gStyle->SetOptStat(0);

	c2->SetLogz(0);
	chi2_corr->GetXaxis()->SetTitle("");
	chi2_corr->GetYaxis()->SetTitle("");
	chi2_corr->Draw("COLZ");

	c2->Print("chi2_corr.png");

	tree->Print();

	file->Write();
	file->Close();
	delete file;
	return 0;
}
