#include <string.h>

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "Math/Functor.h"
#include "Math/Factory.h"
#include "Math/Minimizer.h"

#include "../../Include/Codes/reconstructor.h"
#include "../../Include/const.h"
#include "../../Include/Codes/uncertainties.h"
#include "../../Include/Codes/kinematic_fits.h"
#include "../../Include/Codes/charged_mom.h"
#include "../../Include/Codes/neutral_mom.h"
#include "chain_init.C"

const Int_t N = 23, M = 4;
Int_t loopcount = 4;

Double_t trilateration_chi_square(const Double_t *x)
{
	Double_t clusters[5][4], clusters_meas[5][4], clusters_err[5][4], bhabha_vtx[3];

	Float_t gamma_mom[4][4], solution[4], lambda[4];

	lambda[0] = x[3 * N];
	lambda[1] = x[3 * N + 1];
	lambda[2] = x[3 * N + 2];
	lambda[3] = x[3 * N + 3];

	Double_t kaon_velocity, kaon_energy, kaon_mom, kaon_path;
	Double_t mass_inv[2];

	for (Int_t i = 0; i < 4; i++)
		for (Int_t j = 0; j < 5; j++)
		{
			clusters[j][i] = x[i * 5 + j];
		}

	bhabha_vtx[0] = x[20];
	bhabha_vtx[1] = x[21];
	bhabha_vtx[2] = x[22];

	solution[0] = x[3 * N + M];
	solution[1] = x[3 * N + M + 1];
	solution[2] = x[3 * N + M + 2];
	solution[3] = x[3 * N + M + 3];

	Double_t value = 0, constraints[4] = {0.};

	for (Int_t k = 0; k < 4; k++)
	{
		neutral_mom(clusters[0][k], clusters[1][k], clusters[2][k], clusters[4][k], solution, gamma_mom[k]);

		neutral_mom(clusters[0][k], clusters[1][k], clusters[2][k], clusters[4][k], solution, gamma_mom[k]);
	}

	mass_inv[0] = sqrt(pow(gamma_mom[0][3] + gamma_mom[1][3], 2) -
										 pow(gamma_mom[0][0] + gamma_mom[1][0], 2) -
										 pow(gamma_mom[0][1] + gamma_mom[1][1], 2) -
										 pow(gamma_mom[0][2] + gamma_mom[1][2], 2));
	mass_inv[1] = sqrt(pow(gamma_mom[2][3] + gamma_mom[3][3], 2) -
										 pow(gamma_mom[2][0] + gamma_mom[3][0], 2) -
										 pow(gamma_mom[2][1] + gamma_mom[3][1], 2) -
										 pow(gamma_mom[2][2] + gamma_mom[3][2], 2));

	kaon_mom = sqrt(pow(gamma_mom[0][0] + gamma_mom[1][0] + gamma_mom[2][0] + gamma_mom[3][0], 2) +
									pow(gamma_mom[0][1] + gamma_mom[1][1] + gamma_mom[2][1] + gamma_mom[3][1], 2) +
									pow(gamma_mom[0][2] + gamma_mom[1][2] + gamma_mom[2][2] + gamma_mom[3][2], 2));
	kaon_energy = gamma_mom[0][3] + gamma_mom[1][3] + gamma_mom[2][3] + gamma_mom[3][3];
	kaon_velocity = c_vel * kaon_mom / kaon_energy;
	kaon_path = sqrt(pow(solution[0] - bhabha_vtx[0], 2) +
									 pow(solution[1] - bhabha_vtx[1], 2) +
									 pow(solution[2] - bhabha_vtx[2], 2));

	constraints[0] = pow(kaon_velocity * solution[3] - kaon_path, 2);
	// constraints[1] = pow(mass_inv[0] - m_pi0, 2);
	// constraints[2] = pow(mass_inv[1] - m_pi0, 2);
	// constraints[3] = pow(kaon_path - 320, 2);

	for (Int_t i = 0; i < N; i++)
	{
		value += pow((x[i] - x[i + N]) / x[i + 2 * N], 2);
	}

	value += lambda[0] * constraints[0];
	value += lambda[1] * constraints[1];
	value += lambda[2] * constraints[2];
	value += lambda[3] * constraints[3];

	return value;
}

int main()
{
	TChain *chain = new TChain("INTERF/h1");

	chain_init(chain);

	// TFile *file = new TFile("neuvtx_tri_rec.root", "recreate");
	// TTree *tree = new TTree("h_tri", "Neu vtx rec with trilateration");

	// Branches' addresses
	// Bhabha vars
	Float_t bhabha_mom[4], bhabha_vtx[3], bhabha_vtx_err[3];

	chain->SetBranchAddress("Bpx", &bhabha_mom[0]);
	chain->SetBranchAddress("Bpy", &bhabha_mom[1]);
	chain->SetBranchAddress("Bpz", &bhabha_mom[2]);
	chain->SetBranchAddress("Broots", &bhabha_mom[3]);

	chain->SetBranchAddress("Bx", &bhabha_vtx[0]);
	chain->SetBranchAddress("By", &bhabha_vtx[1]);
	chain->SetBranchAddress("Bz", &bhabha_vtx[2]);
	chain->SetBranchAddress("Bsx", &bhabha_vtx_err[0]);
	chain->SetBranchAddress("Bsy", &bhabha_vtx_err[1]);
	chain->SetBranchAddress("Bsz", &bhabha_vtx_err[2]);

	// Cluster vars
	Int_t nclu;
	UChar_t mctruth;
	Float_t cluster[5][500], Kchboost[9], Knereclor[9], Knemc[9];

	chain->SetBranchAddress("nclu", &nclu);
	chain->SetBranchAddress("Xcl", cluster[0]);
	chain->SetBranchAddress("Ycl", cluster[1]);
	chain->SetBranchAddress("Zcl", cluster[2]);
	chain->SetBranchAddress("Tcl", cluster[3]);
	chain->SetBranchAddress("Enecl", cluster[4]);

	chain->SetBranchAddress("mctruth", &mctruth);

	//! Parameters for fitting

	Double_t P[N], DP[N];

	TMatrixD V_final(N, N);
	TVectorD P_final(N);

	Double_t CHISQR[2] = {0.}, CHISQRMIN[2] = {0.};

	//!

	Int_t nentries = (Int_t)chain->GetEntries();

	Bool_t clusterEnergy, solError, isConverged[2];
	Int_t ind_gam[4], sort_index[2][3], sort_ind_gam[2][4], chosen_ind_gam[4];
	Float_t neu_vtx[2][4], inv_m_pi0[2][3][2], gamma_mom[2][4][4], mass_pair[2][3];

	// Variables
	int selected[4] = {1, 2, 3, 4};

	Reconstructor R;
	Solution S;

	ROOT::Math::Minimizer *minimum;

	for (Int_t i = 0; i < nentries; i++)
	{
		chain->GetEntry(i);

		if (nclu >= 4)
		{
			for (Int_t j1 = 0; j1 < nclu - 3; j1++)
				for (Int_t j2 = j1 + 1; j2 < nclu - 2; j2++)
					for (Int_t j3 = j2 + 1; j3 < nclu - 1; j3++)
						for (Int_t j4 = j3 + 1; j4 < nclu; j4++)
						{
							CHISQRMIN[0] = 999999.;
							CHISQRMIN[1] = 999999.;

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
								////////////////////////////////////////

								for (Int_t k = 0; k < 4; k++)
								{
									neu_vtx[0][k] = S.sol[0][k];
									neu_vtx[1][k] = S.sol[1][k];
								}

								// Pairing of photons into pions
								for (Int_t k = 0; k < 4; k++)
								{
									neutral_mom(cluster[0][ind_gam[k]], cluster[1][ind_gam[k]], cluster[2][ind_gam[k]], cluster[4][ind_gam[k]], neu_vtx[0], gamma_mom[0][k]);

									neutral_mom(cluster[0][ind_gam[k]], cluster[1][ind_gam[k]], cluster[2][ind_gam[k]], cluster[4][ind_gam[k]], neu_vtx[1], gamma_mom[1][k]);
								}

								// All possible combinations for 1st solution
								inv_m_pi0[0][0][0] = sqrt(pow(gamma_mom[0][0][3] + gamma_mom[0][1][3], 2) -
																					pow(gamma_mom[0][0][0] + gamma_mom[0][1][0], 2) -
																					pow(gamma_mom[0][0][1] + gamma_mom[0][1][1], 2) -
																					pow(gamma_mom[0][0][2] + gamma_mom[0][1][2], 2));
								inv_m_pi0[0][0][1] = sqrt(pow(gamma_mom[0][2][3] + gamma_mom[0][3][3], 2) -
																					pow(gamma_mom[0][2][0] + gamma_mom[0][3][0], 2) -
																					pow(gamma_mom[0][2][1] + gamma_mom[0][3][1], 2) -
																					pow(gamma_mom[0][2][2] + gamma_mom[0][3][2], 2));

								inv_m_pi0[0][1][0] = sqrt(pow(gamma_mom[0][0][3] + gamma_mom[0][2][3], 2) -
																					pow(gamma_mom[0][0][0] + gamma_mom[0][2][0], 2) -
																					pow(gamma_mom[0][0][1] + gamma_mom[0][2][1], 2) -
																					pow(gamma_mom[0][0][2] + gamma_mom[0][2][2], 2));
								inv_m_pi0[0][1][1] = sqrt(pow(gamma_mom[0][1][3] + gamma_mom[0][3][3], 2) -
																					pow(gamma_mom[0][1][0] + gamma_mom[0][3][0], 2) -
																					pow(gamma_mom[0][1][1] + gamma_mom[0][3][1], 2) -
																					pow(gamma_mom[0][1][2] + gamma_mom[0][3][2], 2));

								inv_m_pi0[0][2][0] = sqrt(pow(gamma_mom[0][0][3] + gamma_mom[0][3][3], 2) -
																					pow(gamma_mom[0][0][0] + gamma_mom[0][3][0], 2) -
																					pow(gamma_mom[0][0][1] + gamma_mom[0][3][1], 2) -
																					pow(gamma_mom[0][0][2] + gamma_mom[0][3][2], 2));
								inv_m_pi0[0][2][1] = sqrt(pow(gamma_mom[0][2][3] + gamma_mom[0][1][3], 2) -
																					pow(gamma_mom[0][2][0] + gamma_mom[0][1][0], 2) -
																					pow(gamma_mom[0][2][1] + gamma_mom[0][1][1], 2) -
																					pow(gamma_mom[0][2][2] + gamma_mom[0][1][2], 2));

								// All possible combinations for 2nd solution
								inv_m_pi0[1][0][0] = sqrt(pow(gamma_mom[1][0][3] + gamma_mom[1][1][3], 2) -
																					pow(gamma_mom[1][0][0] + gamma_mom[1][1][0], 2) -
																					pow(gamma_mom[1][0][1] + gamma_mom[1][1][1], 2) -
																					pow(gamma_mom[1][0][2] + gamma_mom[1][1][2], 2));
								inv_m_pi0[1][0][1] = sqrt(pow(gamma_mom[1][2][3] + gamma_mom[1][3][3], 2) -
																					pow(gamma_mom[1][2][0] + gamma_mom[1][3][0], 2) -
																					pow(gamma_mom[1][2][1] + gamma_mom[1][3][1], 2) -
																					pow(gamma_mom[1][2][2] + gamma_mom[1][3][2], 2));

								inv_m_pi0[1][1][0] = sqrt(pow(gamma_mom[1][0][3] + gamma_mom[1][2][3], 2) -
																					pow(gamma_mom[1][0][0] + gamma_mom[1][2][0], 2) -
																					pow(gamma_mom[1][0][1] + gamma_mom[1][2][1], 2) -
																					pow(gamma_mom[1][0][2] + gamma_mom[1][2][2], 2));
								inv_m_pi0[1][1][1] = sqrt(pow(gamma_mom[1][1][3] + gamma_mom[1][3][3], 2) -
																					pow(gamma_mom[1][1][0] + gamma_mom[1][3][0], 2) -
																					pow(gamma_mom[1][1][1] + gamma_mom[1][3][1], 2) -
																					pow(gamma_mom[1][1][2] + gamma_mom[1][3][2], 2));

								inv_m_pi0[1][2][0] = sqrt(pow(gamma_mom[1][0][3] + gamma_mom[1][3][3], 2) -
																					pow(gamma_mom[1][0][0] + gamma_mom[1][3][0], 2) -
																					pow(gamma_mom[1][0][1] + gamma_mom[1][3][1], 2) -
																					pow(gamma_mom[1][0][2] + gamma_mom[1][3][2], 2));
								inv_m_pi0[1][2][1] = sqrt(pow(gamma_mom[1][2][3] + gamma_mom[1][1][3], 2) -
																					pow(gamma_mom[1][2][0] + gamma_mom[1][1][0], 2) -
																					pow(gamma_mom[1][2][1] + gamma_mom[1][1][1], 2) -
																					pow(gamma_mom[1][2][2] + gamma_mom[1][1][2], 2));

								// Pairing the photons for pi0

								mass_pair[0][0] = sqrt(pow(inv_m_pi0[0][0][0] - m_pi0, 2) + pow(inv_m_pi0[0][0][1] - m_pi0, 2));
								mass_pair[0][1] = sqrt(pow(inv_m_pi0[0][1][0] - m_pi0, 2) + pow(inv_m_pi0[0][1][1] - m_pi0, 2));
								mass_pair[0][2] = sqrt(pow(inv_m_pi0[0][2][0] - m_pi0, 2) + pow(inv_m_pi0[0][2][1] - m_pi0, 2));

								mass_pair[1][0] = sqrt(pow(inv_m_pi0[1][0][0] - m_pi0, 2) + pow(inv_m_pi0[1][0][1] - m_pi0, 2));
								mass_pair[1][1] = sqrt(pow(inv_m_pi0[1][1][0] - m_pi0, 2) + pow(inv_m_pi0[1][1][1] - m_pi0, 2));
								mass_pair[1][2] = sqrt(pow(inv_m_pi0[1][2][0] - m_pi0, 2) + pow(inv_m_pi0[1][2][1] - m_pi0, 2));

								TMath::Sort(3, mass_pair[0], sort_index[0], kFALSE);
								TMath::Sort(3, mass_pair[1], sort_index[1], kFALSE);

								// !
								// !ITERATING OVER SOLUTIONS OF THE RECONSTRUCTION
								for (Int_t l = 0; l < 2; l++)
								{
									minimum = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");

									// set tolerance , etc...
									minimum->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2
									minimum->SetTolerance(0.01);
									minimum->SetPrintLevel(0);

									// create function wrapper for minimizer
									// a IMultiGenFunction type
									ROOT::Math::Functor f(&trilateration_chi_square, 3 * N + M + 4);

									minimum->SetFunction(f);

									for (Int_t l = 0; l < 3 * N + M + 4; l++)
									{
										minimum->SetVariable(l, std::to_string(l), 0., 0.01);
										if (l >= N)
											minimum->FixVariable(l);
									}

									if (sort_index[l][0] == 0)
									{
										sort_ind_gam[l][0] = ind_gam[0];
										sort_ind_gam[l][1] = ind_gam[1];

										sort_ind_gam[l][2] = ind_gam[2];
										sort_ind_gam[l][3] = ind_gam[3];
									}
									else if (sort_index[l][0] == 1)
									{
										sort_ind_gam[l][0] = ind_gam[0];
										sort_ind_gam[l][1] = ind_gam[2];

										sort_ind_gam[l][2] = ind_gam[1];
										sort_ind_gam[l][3] = ind_gam[3];
									}
									else if (sort_index[l][0] == 2)
									{
										sort_ind_gam[l][0] = ind_gam[0];
										sort_ind_gam[l][1] = ind_gam[3];

										sort_ind_gam[l][2] = ind_gam[1];
										sort_ind_gam[l][3] = ind_gam[2];
									}

									// Building the parameters for a fit
									for (Int_t k = 0; k < 4; k++)
									{
										P[k * 5] = cluster[0][sort_ind_gam[l][k]];		 // X_k
										P[k * 5 + 1] = cluster[1][sort_ind_gam[l][k]]; // Y_k
										P[k * 5 + 2] = cluster[2][sort_ind_gam[l][k]]; // Z_k
										P[k * 5 + 3] = cluster[3][sort_ind_gam[l][k]]; // T_k
										P[k * 5 + 4] = cluster[4][sort_ind_gam[l][k]]; // E_k

										DP[k * 5] = 1.2;		 // cm
										DP[k * 5 + 1] = 1.2; // cm
										DP[k * 5 + 2] = 1.2; // cm
										DP[k * 5 + 3] = clu_time_error(cluster[4][sort_ind_gam[l][k]]);
										DP[k * 5 + 4] = clu_ene_error(cluster[4][sort_ind_gam[l][k]]);
									}

									P[20] = bhabha_vtx[0]; // X_phi
									P[21] = bhabha_vtx[1]; // Y_phi
									P[22] = bhabha_vtx[2]; // Z_phi

									DP[20] = bhabha_vtx_err[0];
									DP[21] = bhabha_vtx_err[1];
									DP[22] = bhabha_vtx_err[2];

									for (Int_t k = 0; k < N; k++)
									{
										minimum->SetVariableValue(k, P[k]);
										minimum->SetVariableLimits(k, P[k] - 3 * DP[k], P[k] + 3 * DP[k]);
										minimum->SetVariableValue(k + N, P[k]);
										minimum->SetVariableValue(k + 2 * N, DP[k]);
									}

									for (Int_t k = 0; k < M; k++)
									{
										minimum->SetVariableValue(k + 3 * N, 1.0);
									}

									minimum->SetVariableValue(3 * N + M, neu_vtx[l][0]);
									minimum->SetVariableValue(3 * N + M + 1, neu_vtx[l][1]);
									minimum->SetVariableValue(3 * N + M + 2, neu_vtx[l][2]);
									minimum->SetVariableValue(3 * N + M + 3, neu_vtx[l][3]);

									// do the minimization
									isConverged[l] = minimum->Minimize();

									// if (isConverged[l])
										

									delete minimum;
								}

							}
							///////////////////////////////////////////////////////////////////
						}
		}
	}
	return 0;
}
