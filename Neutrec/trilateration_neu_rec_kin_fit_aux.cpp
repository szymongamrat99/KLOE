#include <string.h>

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TF1.h"
#include "TH2.h"
#include "TCanvas.h"
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

const Int_t N_free = 24, N_const = 7, M = 2;
const Float_t Trf = 2.715; // ns - time of a bunch (correction)

const Int_t loopcount = 5;

TF1 *constraints[M];

Float_t det, CHISQR[2], CHISQRTMP;
Int_t fail[2];

Reconstructor R;
Solution S;

Int_t selected[4] = {1, 2, 3, 4};

Int_t main(int argc, char *argv[])
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
	Double_t reinitialize[(N_free + N_const)*(N_free + N_const)] = {0.};

	TMatrixD V(N_free + N_const, N_free + N_const), D(M, N_free + N_const), D_T(N_free + N_const, M), V_final[2], Aux1(M, N_free + N_const), Aux(M,M), Aux2(N_free + N_const,M);
	TVectorD X(N_free + N_const), C(M), X_final[2], L(M), CORR(N_free+N_const), X_init(N_free + N_const), X_min(N_free + N_const);

	constraints[0] = new TF1("Ene consv", &ene_consv, 0, 1, N_free + N_const);
	constraints[1] = new TF1("Minv consv", &minv_consv, 0, 1, N_free + N_const);
	constraints[2] = new TF1("x consv", &x_consv, 0, 1, N_free + N_const);
	constraints[3] = new TF1("y consv", &y_consv, 0, 1, N_free + N_const);
	constraints[4] = new TF1("z consv", &z_consv, 0, 1, N_free + N_const);

	X_final[0].ResizeTo(N_free + N_const);
	X_final[1].ResizeTo(N_free + N_const);

	V_final[0].ResizeTo(N_free + N_const, N_free + N_const);
	V_final[1].ResizeTo(N_free + N_const, N_free + N_const);

	TH1 *chi2 = new TH1F("chi2", "", 50, 0, 50);

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
							V.SetMatrixArray(reinitialize);

							CHISQR[0] = 0.;
							CHISQR[1] = 0.;

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

									V(k * 5, k * 5) = pow(1.2,2);
									V(k * 5 + 1, k * 5 + 1) = pow(1.2,2);
									V(k * 5 + 2, k * 5 + 2) = pow(1.2 / sqrt(X(k * 5 + 4) / 1000.),2); // cm
									V(k * 5 + 3, k * 5 + 3) = pow(clu_time_error(X(k * 5 + 4)),2);			// ns
									V(k * 5 + 4, k * 5 + 4) = pow(clu_ene_error(X(k * 5 + 4)),2);			// MeV
								}

								X(20) = bhabha_mom[0];
								X(21) = bhabha_mom[1];
								X(22) = bhabha_mom[2];
								X(23) = bhabha_mom[3];

								V(20, 20) = pow(bhabha_mom_err[0],2);
								V(21, 21) = pow(bhabha_mom_err[1],2);
								V(22, 22) = pow(bhabha_mom_err[2],2);
								V(23, 23) = pow(bhabha_mom_err[3],2);

								V(24, 24) = 0.;
								V(25, 25) = 0.;
								V(26, 26) = 0.;
								V(27, 27) = 0.;

								X(28) = bhabha_vtx[0];
								X(29) = bhabha_vtx[1];
								X(30) = bhabha_vtx[2];

								V(28, 28) = 0.;
								V(29, 29) = 0.;
								V(30, 30) = 0.;

								X_init = X;

								for (Int_t j = 0; j < loopcount; j++)
								{
									// Setting clusters for a solution
									for (Int_t k = 0; k < 4; k++)
										R.SetClu(k, X(k * 5),
														 X(k * 5 + 1),
														 X(k * 5 + 2),
														 X(k * 5 + 3),
														 X(k * 5 + 4));

									R.SetClu(4, 0., 0., 0., 0., 0.);
									R.SetClu(5, 0., 0., 0., 0., 0.);

									S = R.MySolve(selected);

									for (Int_t k = 0; k < 2; k++)
									{
										fail[k] = 0;

										X(24) = S.sol[k][0];
										X(25) = S.sol[k][1];
										X(26) = S.sol[k][2];
										X(27) = S.sol[k][3];

										if (S.error[k] == false && abs(X(24)) < 200 && abs(X(25)) < 200 && abs(X(26)) < 160)
										{

											for (Int_t l = 0; l < M; l++)
											{
												constraints[l]->SetParameters(X.GetMatrixArray());
												C(l) = constraints[l]->Eval(0);
												for (Int_t m = 0; m < N_free + N_const; m++)
												{
													D(l, m) = constraints[l]->GradientPar(m, 0, 0.001);
												}
											}

											D_T.Transpose(D);

											Aux1 = D * V;
											Aux = Aux1 * D_T;

											//Aux.Print();

											det = Aux.Determinant();

											if (!TMath::IsNaN(det) && det != 0)
											{
												Aux.Invert();

												L = Aux * C;

												Aux2 = V * D_T;

												CORR = Aux2 * L;

												X_final[k] = X - CORR;
												V_final[k] = V - V * D_T * Aux * D * V;

												if(X_final[k](4) < 20) X_final[k](4) = 20;
												if(X_final[k](9) < 20) X_final[k](9) = 20;
												if(X_final[k](14) < 20) X_final[k](14) = 20;
												if(X_final[k](19) < 20) X_final[k](19) = 20;

												CHISQR[k] = Dot((X_final[k] - X_init).GetSub(0,N_free-1), V.GetSub(0,N_free - 1, 0, N_free - 1).Invert() * (X_final[k] - X_init).GetSub(0,N_free-1));

											}
											else
											{
												fail[k] = 1;
											}
										}
										else
											fail[k] = 1;

									}

									if(fail[0] == 0 && fail[1] == 1)
									{
										X = X_final[0];
										// V = V_final[0];
										CHISQRTMP = CHISQR[0];
										isConverged = true;
									}
									else if(fail[0] == 1 && fail[1] == 0)
									{
										X = X_final[1];
										// V = V_final[1];
										CHISQRTMP = CHISQR[1];
										isConverged = true;
									}
									else if(fail[0] == 0 && fail[1] == 0)
									{
										if(CHISQR[0] < CHISQR[1])
										{
											X = X_final[0];
											// V = V_final[0];
											CHISQRTMP = CHISQR[0];
											isConverged = true;
										}
										else
										{
											X = X_final[1];
											// V = V_final[1];
											CHISQRTMP = CHISQR[1];
											isConverged = true;
										}
									}
									else
									{
										CHISQRTMP = 999999.;
										isConverged = false;
										break;
									}

								}

								if(CHISQRTMP < CHISQRMIN)
								{
									isConverged = true;
									CHISQRMIN = CHISQRTMP;
									X_min = X;
								}
							}
						}

						if(isConverged)
						{
							cout << constraints[1]->EvalPar(0, X_min.GetMatrixArray()) << endl;
							chi2->Fill(CHISQRMIN);
						}
		}
	}

	TF1 *func = new TF1("chi2dist", chi2dist, 0, 50, 2);

	func->SetParameters(3000.0, (Double_t)M);
	func->SetParNames("Norm", "Degrees of Freedom");

	func->SetParLimits(0, 0.001, 5000.0);
	func->SetParLimits(1, 1.0, 10.0);

	TCanvas *c1 = new TCanvas("c1","",750,750);
	chi2->Fit(func);

	chi2->GetYaxis()->SetRangeUser(0,1.2*chi2->GetMaximum());
	chi2->GetXaxis()->SetRangeUser(0,50);
	chi2->Draw();
	c1->Print("test_chi2.png");

	return 0;
}