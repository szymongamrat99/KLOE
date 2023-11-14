#include <string.h>

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TF1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "Math/Functor.h"
#include "Math/Factory.h"
#include "Math/Minimizer.h"
#include "TMatrixD.h"
#include "TVectorD.h"

#include "../../Include/Codes/reconstructor.h"
#include "../../Include/const.h"
#include "../../Include/Codes/uncertainties.h"
#include "../../Include/Codes/charged_mom.h"
#include "../../Include/Codes/neutral_mom.h"
#include "../../Include/Codes/lorentz_transf.h"
#include "../../Include/Codes/plane_intersection.h"
#include "../../Include/Codes/closest_approach.h"
#include "../../Include/Codes/constraints_tri.h"
#include "../../Include/Codes/chi2_dist.h"
#include "chain_init.C"

using namespace std;

const Int_t N_free = 24, N_const = 3, M = 5;
const Float_t Trf = 2.715; // ns - time of a bunch (correction)

const Int_t loopcount = 5;

TF1 *constraints[M];

Float_t det, CHISQR, CHISQRTMP, FUNVAL, FUNVALTMP, FUNVALMIN;
Int_t fail;

Reconstructor R;
Solution S;

Int_t selected[4] = {1, 2, 3, 4};

Int_t main(int argc, char *argv[])
{
	UInt_t first = atoi(argv[2]), last = atoi(argv[3]);

	TString first_s = argv[2], last_s = argv[3];

	TChain *chain = new TChain("INTERF/h1");
	chain_init(chain, first, last);

	TString name = "neuvtx_tri_kin_fit_" + first_s + "_" + last_s + "_" + loopcount + "_" + M + ".root";

	TFile *file = new TFile(name, "recreate");
	TTree *tree = new TTree("h_tri_kin_fit", "Neu vtx rec with trilateration kin fit");

	// Branches' addresses
	// Bhabha vars
	Float_t bhabha_mom[4], bhabha_mom_err[4], bhabha_vtx[3];

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

	Int_t nentries = (Int_t)chain->GetEntries();

	UInt_t number_of_ev = atoi(argv[1]);
	TString number_of_ev_s = argv[1];

	if (number_of_ev_s == "all")
		number_of_ev = nentries;
	else if (number_of_ev > nentries)
		number_of_ev = nentries;

	Bool_t clusterEnergy, solError;
	Int_t ind_gam[4], sort_index[2][3], sort_ind_gam[2][4], chosen_ind_gam[4], found_best, isConverged;
	Float_t CHISQRMIN, min_value_def;
	Double_t reinitialize[(N_free + N_const) * (N_free + N_const)] = {0.};
	Float_t gamma_mom_min[4][4], neu_vtx_min[4], kaon_mom_min[4], kaon_vel[3], kaon_vel_tot, kaon_path_tot, bhabha_vtx_min[3], ip_min[3], time_diff[2][4], time_diff_fin, gamma_path[2][4], neu_vtx[2][4];

	TMatrixD V(N_free + N_const, N_free + N_const), D(M, N_free + N_const), D_T(N_free + N_const, M), V_final(N_free + N_const, N_free + N_const), V_aux(N_free + N_const, N_free + N_const), V_min(N_free + N_const, N_free + N_const), Aux(M, M), V_invert(N_free, N_free), V_init(N_free + N_const, N_free + N_const);
	TVectorD X(N_free + N_const), C(M), X_final(N_free + N_const), L(M), CORR(N_free + N_const), X_init(N_free + N_const), X_min(N_free + N_const), C_min(M), L_min(M), C_aux(M), L_aux(M), X_init_min(N_free + N_const), X_init_aux(N_free + N_const);

	Float_t length_ch, time_ch, velocity_ch, gamma_mom_final[4][8], fourKnetri_kinfit[10], iptri_kinfit[3], y_axis[3], ip_tri[2][3], bhabha_mom_fit[4], gamma_mom_tmp[2][4][4], fourKnetri_tmp[2][4], value[2], kaon_vel_tmp[2], dist_tmp[2], ip_tmp[2][3];
	Int_t g4takentri_kinfit[4];

	TBranch *b_gamma1tri = tree->Branch("fourgamma1tri_kinfit", gamma_mom_final[0], "fourgamma1tri_kinfit[8]/F");
	TBranch *b_gamma2tri = tree->Branch("fourgamma2tri_kinfit", gamma_mom_final[1], "fourgamma2tri_kinfit[8]/F");
	TBranch *b_gamma3tri = tree->Branch("fourgamma3tri_kinfit", gamma_mom_final[2], "fourgamma3tri_kinfit[8]/F");
	TBranch *b_gamma4tri = tree->Branch("fourgamma4tri_kinfit", gamma_mom_final[3], "fourgamma4tri_kinfit[8]/F");
	TBranch *b_iptri = tree->Branch("iptri_kinfit", iptri_kinfit, "iptri_kinfit[3]/F");
	TBranch *b_Knetri = tree->Branch("fourKnetri_kinfit", fourKnetri_kinfit, "fourKnetri_kinfit[10]/F");
	TBranch *b_done = tree->Branch("done4_kinfit", &isConverged, "done4_kinfit/I");
	TBranch *b_fourg4taken = tree->Branch("g4takentri_kinfit", g4takentri_kinfit, "g4takentri_kinfit[4]/I");

	TBranch *b_Lmin = tree->Branch("lag_mult", "TVectorD", &L_min);
	TBranch *b_Cmin = tree->Branch("const_min", "TVectorD", &C_min);
	TBranch *b_Xmin = tree->Branch("min_vars", "TVectorD", &X_min);
	TBranch *b_Vmin = tree->Branch("min_cov", "TMatrixD", &V_min);
	TBranch *b_Xinit = tree->Branch("init_vars", "TVectorD", &X_init_min);

	TBranch *b_chisqr = tree->Branch("chi2min", &CHISQRMIN, "chi2min/F");

	constraints[0] = new TF1("Ene consv", &ene_consv, 0, 1, N_free + N_const);
	constraints[4] = new TF1("Minv consv", &minv_consv, 0, 1, N_free + N_const);
	constraints[5] = new TF1("gamma1 consv", &gamma1_consv, 0, 1, N_free + N_const);
	constraints[6] = new TF1("gamma2 consv", &gamma2_consv, 0, 1, N_free + N_const);
	constraints[7] = new TF1("gamma3 consv", &gamma3_consv, 0, 1, N_free + N_const);
	constraints[8] = new TF1("gamma4 consv", &gamma4_consv, 0, 1, N_free + N_const);
	constraints[1] = new TF1("x consv", &x_consv, 0, 1, N_free + N_const);
	constraints[2] = new TF1("y consv", &y_consv, 0, 1, N_free + N_const);
	constraints[3] = new TF1("z consv", &z_consv, 0, 1, N_free + N_const);

	TH1 *chi2 = new TH1F("chi2", "", 50, 0.0, 1.0);

	for (Int_t i = 0; i < number_of_ev; i++)
	{
		chain->GetEntry(i);

		min_value_def = 99999.;
		FUNVALMIN = 999999.;
		CHISQRMIN = 999999.;

		isConverged = 0;

		if (nclu >= 4 && (mctruth == 1 || mctruth == 2))
		{

			for (Int_t j1 = 0; j1 < nclu - 3; j1++)
				for (Int_t j2 = j1 + 1; j2 < nclu - 2; j2++)
					for (Int_t j3 = j2 + 1; j3 < nclu - 1; j3++)
						for (Int_t j4 = j3 + 1; j4 < nclu; j4++)
						{
							V.SetMatrixArray(reinitialize);
							fail = 0;

							CHISQR = 999999.;
							FUNVALTMP = 999999.;

							ind_gam[0] = j1;
							ind_gam[1] = j2;
							ind_gam[2] = j3;
							ind_gam[3] = j4;

							clusterEnergy = (cluster[4][ind_gam[0]] >= 20 && cluster[4][ind_gam[1]] >= 20 && cluster[4][ind_gam[2]] >= 20 && cluster[4][ind_gam[3]] >= 20);

							if (clusterEnergy)
							{

								for (Int_t k = 0; k < 4; k++)
								{
									X(k * 5) = cluster[0][ind_gam[k]];
									X(k * 5 + 1) = cluster[1][ind_gam[k]];
									X(k * 5 + 2) = cluster[2][ind_gam[k]];
									X(k * 5 + 3) = cluster[3][ind_gam[k]];
									X(k * 5 + 4) = cluster[4][ind_gam[k]];

									V(k * 5, k * 5) = pow(1.2, 2);
									V(k * 5 + 1, k * 5 + 1) = pow(1.2, 2);
									V(k * 5 + 2, k * 5 + 2) = pow(1.2 / sqrt(X(k * 5 + 4) / 1000.), 2); // cm
									V(k * 5 + 3, k * 5 + 3) = pow(clu_time_error(X(k * 5 + 4)), 2);			// ns
									V(k * 5 + 4, k * 5 + 4) = pow(clu_ene_error(X(k * 5 + 4)), 2);			// MeV
								}

								X(20) = bhabha_mom[0];
								X(21) = bhabha_mom[1];
								X(22) = bhabha_mom[2];
								X(23) = bhabha_mom[3];

								V(20, 20) = pow(bhabha_mom_err[0], 2);
								V(21, 21) = pow(bhabha_mom_err[1], 2);
								V(22, 22) = pow(bhabha_mom_err[2], 2);
								V(23, 23) = pow(bhabha_mom_err[3], 2);

								X(24) = bhabha_vtx[0];
								X(25) = bhabha_vtx[1];
								X(26) = bhabha_vtx[2];

								V(24, 24) = 0.;
								V(25, 25) = 0.;
								V(26, 26) = 0.;

								X_init = X;

								for (Int_t j = 0; j < loopcount; j++)
								{
									for (Int_t l = 0; l < M; l++)
									{
										C(l) = constraints[l]->EvalPar(0, X.GetMatrixArray());
										for (Int_t m = 0; m < N_free + N_const; m++)
										{
											constraints[l]->SetParameters(X.GetMatrixArray());
											if (m < N_free)
												D(l, m) = constraints[l]->GradientPar(m, 0, 0.0001);
											else
												D(l, m) = 0;
										}
									}

									if (1)
									{
										D_T.Transpose(D);

										Aux = D * V * D_T;

										det = Aux.Determinant();

										if (!TMath::IsNaN(det) && det != 0)
										{
											Aux.Invert();

											L = Aux * C;

											CORR = V * D_T * L;

											X_final = X - CORR;
											V_final = V - V * D_T * Aux * D * V;

											// Avoiding small energy photons
											/*if (X_final(4) < 20)
												X_final(4) = 20;

											if (X_final(9) < 20)
												X_final(9) = 20;

											if (X_final(14) < 20)
												X_final(14) = 20;

											if (X_final(19) < 20)
												X_final(19) = 20;*/

											V_invert = V.GetSub(0, N_free - 1, 0, N_free - 1).Invert();

											FUNVAL = Dot((X_final - X_init).GetSub(0, N_free - 1), V_invert * (X_final - X_init).GetSub(0, N_free - 1));

											CHISQR = Dot((X_final - X_init).GetSub(0, N_free - 1), V_invert * (X_final - X_init).GetSub(0, N_free - 1));
										}
										else
										{
											fail = 1;
										}
									}
									else
									{
										break;
									}

									if (fail == 0)
									{
										X = X_final;
										X_init_aux = X_init;
										V_aux = V_final;
										L_aux = L;
										C_aux = C;
										FUNVALTMP = FUNVAL;
										CHISQRTMP = CHISQR;
									}
									else
										break;
								}

								if (abs(CHISQRTMP) < abs(CHISQRMIN))
								{
									isConverged = 1;
									FUNVALMIN = FUNVALTMP;
									CHISQRMIN = CHISQRTMP;
									X_min = X;
									X_init_min = X_init_aux;
									V_min = V_aux;
									V_init = V;
									C_min = C_aux;
									L_min = L_aux;

									g4takentri_kinfit[0] = ind_gam[0];
									g4takentri_kinfit[1] = ind_gam[1];
									g4takentri_kinfit[2] = ind_gam[2];
									g4takentri_kinfit[3] = ind_gam[3];
								}
							}
						}

			if (isConverged)
			{
				// Setting clusters for a solution
				for (Int_t k = 0; k < 4; k++)
				{
					R.SetClu(k, X_min.GetMatrixArray()[k * 5],
									 X_min.GetMatrixArray()[k * 5 + 1],
									 X_min.GetMatrixArray()[k * 5 + 2],
									 X_min.GetMatrixArray()[k * 5 + 3],
									 X_min.GetMatrixArray()[k * 5 + 4]);

					R.SetClu(4, 0., 0., 0., 0., 0.);
					R.SetClu(5, 0., 0., 0., 0., 0.);
				}

				S = R.MySolve(selected);

				for (Int_t k = 0; k < 2; k++)
				{
					if (!S.error[k])
					{
						neu_vtx[k][0] = S.sol[k][0];
						neu_vtx[k][1] = S.sol[k][1];
						neu_vtx[k][2] = S.sol[k][2];
						neu_vtx[k][3] = S.sol[k][3];
					}
					else
					{
						neu_vtx[k][0] = 999.;
						neu_vtx[k][1] = 999.;
						neu_vtx[k][2] = 999.;
						neu_vtx[k][3] = 999.;
					}

					neutral_mom(X_min.GetMatrixArray()[0], X_min.GetMatrixArray()[1], X_min.GetMatrixArray()[2], X_min.GetMatrixArray()[4], neu_vtx[k], gamma_mom_tmp[k][0]);
					neutral_mom(X_min.GetMatrixArray()[5], X_min.GetMatrixArray()[6], X_min.GetMatrixArray()[7], X_min.GetMatrixArray()[9], neu_vtx[k], gamma_mom_tmp[k][1]);
					neutral_mom(X_min.GetMatrixArray()[10], X_min.GetMatrixArray()[11], X_min.GetMatrixArray()[12], X_min.GetMatrixArray()[14], neu_vtx[k], gamma_mom_tmp[k][2]);
					neutral_mom(X_min.GetMatrixArray()[15], X_min.GetMatrixArray()[16], X_min.GetMatrixArray()[17], X_min.GetMatrixArray()[19], neu_vtx[k], gamma_mom_tmp[k][3]);

					fourKnetri_tmp[k][0] = gamma_mom_tmp[k][0][0] + gamma_mom_tmp[k][1][0] + gamma_mom_tmp[k][2][0] + gamma_mom_tmp[k][3][0];
					fourKnetri_tmp[k][1] = gamma_mom_tmp[k][0][1] + gamma_mom_tmp[k][1][1] + gamma_mom_tmp[k][2][1] + gamma_mom_tmp[k][3][1];
					fourKnetri_tmp[k][2] = gamma_mom_tmp[k][0][2] + gamma_mom_tmp[k][1][2] + gamma_mom_tmp[k][2][2] + gamma_mom_tmp[k][3][2];
					fourKnetri_tmp[k][3] = gamma_mom_tmp[k][0][3] + gamma_mom_tmp[k][1][3] + gamma_mom_tmp[k][2][3] + gamma_mom_tmp[k][3][3];
					fourKnetri_tmp[k][4] = sqrt(pow(fourKnetri_tmp[k][0], 2) + pow(fourKnetri_tmp[k][1], 2) + pow(fourKnetri_tmp[k][2], 2));
					fourKnetri_tmp[k][5] = sqrt(pow(fourKnetri_tmp[k][3], 2) - pow(fourKnetri_tmp[k][4], 2));

					kaon_vel_tmp[k] = c_vel * fourKnetri_tmp[k][4] / fourKnetri_tmp[k][3];

					y_axis[0] = 0.;
					y_axis[1] = X_min.GetMatrixArray()[21];
					y_axis[2] = 0.;

					plane_intersection(bhabha_vtx, y_axis, neu_vtx[k], fourKnetri_tmp[k], ip_tmp[k]);

					ip_tmp[k][0] = bhabha_vtx[0];
					ip_tmp[k][1] = bhabha_vtx[1];
					if (abs(ip_tmp[k][2] - bhabha_vtx[2]) > 2)
						ip_tmp[k][2] = bhabha_vtx[2];

					dist_tmp[k] = sqrt(pow(neu_vtx[k][0] - ip_tmp[k][0], 2) +
														 pow(neu_vtx[k][1] - ip_tmp[k][1], 2) +
														 pow(neu_vtx[k][2] - ip_tmp[k][2], 2));

					value[k] = sqrt(pow(neu_vtx[k][3] - (dist_tmp[k] / kaon_vel_tmp[k]), 2) + pow(fourKnetri_tmp[k][5] - m_k0, 2));
				}

				if (value[0] < value[1])
				{
					neu_vtx_min[0] = neu_vtx[0][0];
					neu_vtx_min[1] = neu_vtx[0][1];
					neu_vtx_min[2] = neu_vtx[0][2];
					neu_vtx_min[3] = neu_vtx[0][3];

					for (Int_t k = 0; k < 4; k++)
					{
						for (Int_t l = 0; l < 4; l++)
						{
							gamma_mom_final[k][l] = gamma_mom_tmp[0][k][l];
						}

						gamma_mom_final[k][4] = X_min.GetMatrixArray()[k * 5];
						gamma_mom_final[k][5] = X_min.GetMatrixArray()[k * 5 + 1];
						gamma_mom_final[k][6] = X_min.GetMatrixArray()[k * 5 + 2];
						gamma_mom_final[k][7] = X_min.GetMatrixArray()[k * 5 + 3];
					}

					fourKnetri_kinfit[0] = fourKnetri_tmp[0][0];
					fourKnetri_kinfit[1] = fourKnetri_tmp[0][1];
					fourKnetri_kinfit[2] = fourKnetri_tmp[0][2];
					fourKnetri_kinfit[3] = fourKnetri_tmp[0][3];
					fourKnetri_kinfit[4] = fourKnetri_tmp[0][4];
					fourKnetri_kinfit[5] = fourKnetri_tmp[0][5];
					fourKnetri_kinfit[6] = neu_vtx_min[0];
					fourKnetri_kinfit[7] = neu_vtx_min[1];
					fourKnetri_kinfit[8] = neu_vtx_min[2];
					fourKnetri_kinfit[9] = neu_vtx_min[3];

					iptri_kinfit[0] = ip_tmp[0][0];
					iptri_kinfit[1] = ip_tmp[0][1];
					iptri_kinfit[2] = ip_tmp[0][2];
				}
				else
				{
					neu_vtx_min[0] = neu_vtx[1][0];
					neu_vtx_min[1] = neu_vtx[1][1];
					neu_vtx_min[2] = neu_vtx[1][2];
					neu_vtx_min[3] = neu_vtx[1][3];

					for (Int_t k = 0; k < 4; k++)
					{
						for (Int_t l = 0; l < 4; l++)
						{
							gamma_mom_final[k][l] = gamma_mom_tmp[1][k][l];
						}

						gamma_mom_final[k][4] = X_min.GetMatrixArray()[k * 5];
						gamma_mom_final[k][5] = X_min.GetMatrixArray()[k * 5 + 1];
						gamma_mom_final[k][6] = X_min.GetMatrixArray()[k * 5 + 2];
						gamma_mom_final[k][7] = X_min.GetMatrixArray()[k * 5 + 3];
					}

					fourKnetri_kinfit[0] = fourKnetri_tmp[1][0];
					fourKnetri_kinfit[1] = fourKnetri_tmp[1][1];
					fourKnetri_kinfit[2] = fourKnetri_tmp[1][2];
					fourKnetri_kinfit[3] = fourKnetri_tmp[1][3];
					fourKnetri_kinfit[4] = fourKnetri_tmp[1][4];
					fourKnetri_kinfit[5] = fourKnetri_tmp[1][5];
					fourKnetri_kinfit[6] = neu_vtx_min[0];
					fourKnetri_kinfit[7] = neu_vtx_min[1];
					fourKnetri_kinfit[8] = neu_vtx_min[2];
					fourKnetri_kinfit[9] = neu_vtx_min[3];

					iptri_kinfit[0] = ip_tmp[1][0];
					iptri_kinfit[1] = ip_tmp[1][1];
					iptri_kinfit[2] = ip_tmp[1][2];
				}
				//chi2->Fill(TMath::Prob(CHISQRMIN, M));

			}
		}

		tree->Fill();
	}

	/*gStyle->SetOptStat("iMr");
	TF1 *func = new TF1("chi2dist", chi2dist, 0, 100, 2);

	func->SetParameters(3000.0, (Double_t)M);
	func->SetParNames("Norm", "Degrees of Freedom");

	func->SetParLimits(0, 0.001, 5000.0);
	func->SetParLimits(1, 1.0, 10.0);

	TCanvas *c1 = new TCanvas("c1", "", 750, 750);
	//chi2->Fit(func);
	c1->SetLogy(1);

	chi2->GetYaxis()->SetRangeUser(1, 1.2 * chi2->GetMaximum());
	chi2->Draw();
	c1->Print("test_chi2.png");*/

	tree->Print();

	file->Write();
	file->Close();
	delete file;
	return 0;
}