#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"

#include "../../../Include/Codes/reconstructor.h"
#include "../../../Include/const.h"
#include "../../../Include/Codes/plane_intersection.h"
#include "../inc/trilateration.hpp"

void tri_neurec(int data_type, int first_file, int last_file, int good_clus) //	arguments are: 1. Number of points
																															//                2. First analyzed file
																															//                3. Last analyzed file

{
	TChain *chain = new TChain("INTERF/h1");
	chain_init(chain, first_file, last_file);

	TFile *file_gen = new TFile("../Generated_vars/gen_vars_1_56.root");
	TTree *tree_gen = (TTree *)file_gen->Get("h_gen_vars");

	TString name = "neuvtx_tri_rec_" + std::to_string(first_file) + "_" + std::to_string(last_file) + ".root";

	TFile *file = new TFile(name, "recreate");
	TTree *tree = new TTree("h_tri", "Neu vtx rec with trilateration");

	// Branches' addresses
	// Bhabha vars
	Float_t bhabha_mom[4], bhabha_vtx[3];

	chain->SetBranchAddress("Bpx", &bhabha_mom[0]);
	chain->SetBranchAddress("Bpy", &bhabha_mom[1]);
	chain->SetBranchAddress("Bpz", &bhabha_mom[2]);
	chain->SetBranchAddress("Broots", &bhabha_mom[3]);

	chain->SetBranchAddress("Bx", &bhabha_vtx[0]);
	chain->SetBranchAddress("By", &bhabha_vtx[1]);
	chain->SetBranchAddress("Bz", &bhabha_vtx[2]);

	// Cluster vars
	Int_t nclu;
	UChar_t mctruth;
	Float_t cluster[5][500], Kchboost[9], Knereclor[9], Knemc[9];

	chain->SetBranchAddress("nclu", &nclu);
	chain->SetBranchAddress("Xacl", cluster[0]);
	chain->SetBranchAddress("Yacl", cluster[1]);
	chain->SetBranchAddress("Zacl", cluster[2]);
	chain->SetBranchAddress("Tcl", cluster[3]);
	chain->SetBranchAddress("Enecl", cluster[4]);

	chain->SetBranchAddress("mctruth", &mctruth);

	Int_t good_clus_ind[4];
	Float_t pgammc[4][8];

	tree_gen->SetBranchAddress("pgammc1", pgammc[0]);
	tree_gen->SetBranchAddress("pgammc2", pgammc[1]);
	tree_gen->SetBranchAddress("pgammc3", pgammc[2]);
	tree_gen->SetBranchAddress("pgammc4", pgammc[3]);

	tree_gen->SetBranchAddress("clusindgood", good_clus_ind);

	chain->AddFriend(tree_gen);

	Int_t nentries = (Int_t)chain->GetEntries();

	// Variables
	Bool_t cond_ene, cond_clus[4], cond_sol[2][2], onesol;
	Int_t ind_gam[4], error = 0, index = 0;
	Float_t ene_sum, gamma_mom[2][6][4], gamma_len[2][6],
			kaon_mom[2][4], kaon_vel[2], kaon_len[2], path_diff[2],
			sol_tmp[4], total_err = 0., total_err_def = 0.,
									gamma_len_tmp[4], kaon_mom_tmp[4], kaon_inv_mass, ene_sum_tmp;
	Double_t solution[4];
	Float_t gammatri_tmp[4][8], Knetri_tmp[10];

	Float_t y_axis[3], neu_vtx[3];

	Int_t selected[4] = {1, 2, 3, 4};

	Reconstructor R;
	Solution S;

	Int_t done, fourg4taken[4], chosen;
	Float_t gammatri[4][8], Knetri[10], totalerr, ip_tri[3], sol1[4], sol2[4], sol1err, sol2err;

	TBranch *b_gamma1tri = tree->Branch("fourgamma1tri", gammatri[0], "fourgamma1tri[8]/F");
	TBranch *b_gamma2tri = tree->Branch("fourgamma2tri", gammatri[1], "fourgamma2tri[8]/F");
	TBranch *b_gamma3tri = tree->Branch("fourgamma3tri", gammatri[2], "fourgamma3tri[8]/F");
	TBranch *b_gamma4tri = tree->Branch("fourgamma4tri", gammatri[3], "fourgamma4tri[8]/F");

	TBranch *b_Knetri = tree->Branch("fourKnetri", Knetri, "fourKnetri[10]/F");
	TBranch *b_iptri = tree->Branch("iptri", ip_tri, "iptri[3]/F");

	TBranch *b_done = tree->Branch("done4", &done, "done4/I");

	TBranch *b_sol1 = tree->Branch("sol1", sol1, "sol1[4]/F");
	TBranch *b_sol2 = tree->Branch("sol2", sol2, "sol2[4]/F");

	TBranch *b_sol1err = tree->Branch("sol1err", &sol1err, "sol1err/F");
	TBranch *b_sol2err = tree->Branch("sol2err", &sol2err, "sol2err/F");

	TBranch *b_chosen = tree->Branch("chosen", &chosen, "chosen/I");

	TBranch *b_fourg4taken = tree->Branch("fourg4taken", fourg4taken, "fourg4taken[4]/I");

	Float_t distance = 0., distance_mc = 0., distance_diff = 0.;
	Float_t mc_dist[90], sigmas[90];

	Float_t eqn_check[4], eqn_check_tmp[2][4], eqn_check_tmp_tot[2], kaon_inv_mass_tmp[2];

	for (Int_t i = 0; i < nentries; i++)
	{
		chain->GetEntry(i);

		for (Int_t j1 = 0; j1 < 4; j1++)
			for (Int_t j2 = 0; j2 < 8; j2++)
				gammatri[j1][j2] = -999.;

		for (Int_t j = 0; j < 10; j++)
			Knetri[j] = -999.;

		for (Int_t j = 0; j < 4; j++)
			fourg4taken[j] = -999;

		totalerr = -999.;

		done = 0;

		total_err_def = 999999.;

		if (nclu >= 4 && (mctruth == 1 || mctruth == 2))
		{
			if (good_clus == 0)
			{
				for (Int_t j1 = 0; j1 < nclu - 3; j1++)
					for (Int_t j2 = j1 + 1; j2 < nclu - 2; j2++)
						for (Int_t j3 = j2 + 1; j3 < nclu - 1; j3++)
							for (Int_t j4 = j3 + 1; j4 < nclu; j4++)
							{
								ene_sum = 0.;
								error = 2;
								total_err = 999999.;

								ind_gam[0] = j1;
								ind_gam[1] = j2;
								ind_gam[2] = j3;
								ind_gam[3] = j4;

								cond_ene = cluster[4][ind_gam[0]] > MIN_CLU_ENE && cluster[4][ind_gam[1]] > MIN_CLU_ENE &&
													 cluster[4][ind_gam[2]] > MIN_CLU_ENE && cluster[4][ind_gam[3]] > MIN_CLU_ENE;

								cond_clus[0] = cluster[3][ind_gam[0]] > 0 && cluster[0][ind_gam[0]] != 0 && cluster[1][ind_gam[0]] != 0 && cluster[2][ind_gam[0]] != 0;
								cond_clus[1] = cluster[3][ind_gam[1]] > 0 && cluster[0][ind_gam[1]] != 0 && cluster[1][ind_gam[1]] != 0 && cluster[2][ind_gam[1]] != 0;
								cond_clus[2] = cluster[3][ind_gam[2]] > 0 && cluster[0][ind_gam[2]] != 0 && cluster[1][ind_gam[2]] != 0 && cluster[2][ind_gam[2]] != 0;
								cond_clus[3] = cluster[3][ind_gam[3]] > 0 && cluster[0][ind_gam[3]] != 0 && cluster[1][ind_gam[3]] != 0 && cluster[2][ind_gam[3]] != 0;

								if (cond_ene == true && cond_clus[0] && cond_clus[1] && cond_clus[2] && cond_clus[3])
								{
									for (Int_t l = 0; l < 4; l++)
									{
										R.SetClu(l, cluster[0][ind_gam[l]],
														 cluster[1][ind_gam[l]],
														 cluster[2][ind_gam[l]],
														 cluster[3][ind_gam[l]],
														 cluster[4][ind_gam[l]]);
									}

									R.SetClu(4, 0, 0, 0, 0, 0);
									R.SetClu(5, 0, 0, 0, 0, 0);

									S = R.MySolve(selected);

									if (S.error[0] == false || S.error[1] == false) // && cond_ene_sum == true)
									{

										kaon_mom[0][0] = 0.;
										kaon_mom[0][1] = 0.;
										kaon_mom[0][2] = 0.;
										kaon_mom[0][3] = 0.;

										kaon_mom[1][0] = 0.;
										kaon_mom[1][1] = 0.;
										kaon_mom[1][2] = 0.;
										kaon_mom[1][3] = 0.;

										error = 1;

										for (Int_t l1 = 0; l1 < 2; l1++)
										{
											for (Int_t l2 = 0; l2 < 4; l2++)
											{
												gamma_len[l1][l2] = sqrt(pow(cluster[0][ind_gam[l2]] - S.sol[l1][0], 2) +
																								 pow(cluster[1][ind_gam[l2]] - S.sol[l1][1], 2) +
																								 pow(cluster[2][ind_gam[l2]] - S.sol[l1][2], 2));

												gamma_mom[l1][l2][0] = cluster[4][ind_gam[l2]] * (cluster[0][ind_gam[l2]] - S.sol[l1][0]) / gamma_len[l1][l2];
												gamma_mom[l1][l2][1] = cluster[4][ind_gam[l2]] * (cluster[1][ind_gam[l2]] - S.sol[l1][1]) / gamma_len[l1][l2];
												gamma_mom[l1][l2][2] = cluster[4][ind_gam[l2]] * (cluster[2][ind_gam[l2]] - S.sol[l1][2]) / gamma_len[l1][l2];
												gamma_mom[l1][l2][3] = cluster[4][ind_gam[l2]];

												kaon_mom[l1][0] += gamma_mom[l1][l2][0];
												kaon_mom[l1][1] += gamma_mom[l1][l2][1];
												kaon_mom[l1][2] += gamma_mom[l1][l2][2];
												kaon_mom[l1][3] += gamma_mom[l1][l2][3];
											}

											kaon_inv_mass_tmp[l1] = abs(sqrt(pow(kaon_mom[l1][3],2) - pow(kaon_mom[l1][0],2) - pow(kaon_mom[l1][1],2) - pow(kaon_mom[l1][2],2)) - mK0);

											y_axis[0] = 0.;
											y_axis[1] = bhabha_mom[1];
											y_axis[2] = 0.;

											plane_intersection(bhabha_vtx, y_axis, S.sol[l1], kaon_mom[l1], ip_tri); //! Plane rec

											if (abs(ip_tri[2] - bhabha_vtx[2]) > 2)
												ip_tri[2] = bhabha_vtx[2];

											ip_tri[0] = bhabha_vtx[0];
											ip_tri[1] = bhabha_vtx[1];

											kaon_vel[l1] = cVel * sqrt(pow(kaon_mom[l1][0], 2) + pow(kaon_mom[l1][1], 2) + pow(kaon_mom[l1][2], 2)) / (kaon_mom[l1][3]);
											kaon_len[l1] = sqrt(pow(S.sol[l1][0] - ip_tri[0], 2) + pow(S.sol[l1][1] - ip_tri[1], 2) + pow(S.sol[l1][2] - ip_tri[2], 2));
										}

										cond_sol[0][0] = (S.sol[0][3] < cluster[3][ind_gam[0]]) && (S.sol[0][3] < cluster[3][ind_gam[1]]) && (S.sol[0][3] < cluster[3][ind_gam[2]]) && (S.sol[0][3] < cluster[3][ind_gam[3]]);
										cond_sol[0][1] = (sqrt(pow(S.sol[0][0], 2) + pow(S.sol[0][1], 2)) < 200) && (abs(S.sol[0][2]) < 165);
										path_diff[0] = kaon_vel[0] * S.sol[0][3] - kaon_len[0];


										cond_sol[1][0] = (S.sol[1][3] < cluster[3][ind_gam[0]]) && (S.sol[1][3] < cluster[3][ind_gam[1]]) && (S.sol[1][3] < cluster[3][ind_gam[2]]) && (S.sol[1][3] < cluster[3][ind_gam[3]]);
										cond_sol[1][1] = (sqrt(pow(S.sol[1][0], 2) + pow(S.sol[1][1], 2)) < 200) && (abs(S.sol[1][2]) < 165);
										path_diff[1] = kaon_vel[1] * S.sol[1][3] - kaon_len[1];

										sol1[0] = S.sol[0][0];
										sol1[1] = S.sol[0][1];
										sol1[2] = S.sol[0][2];
										sol1[3] = S.sol[0][3];

										sol2[0] = S.sol[1][0];
										sol2[1] = S.sol[1][1];
										sol2[2] = S.sol[1][2];
										sol2[3] = S.sol[1][3];

										sol1err = path_diff[0];
										sol2err = path_diff[1];

										// We have to choose solutions

										if ((cond_sol[1][0] == true && cond_sol[1][1] == true) &&
															abs(path_diff[1]) + kaon_inv_mass_tmp[1] < abs(path_diff[0]) + kaon_inv_mass_tmp[0])
										{
											solution[0] = S.sol[1][0];
											solution[1] = S.sol[1][1];
											solution[2] = S.sol[1][2];
											solution[3] = S.sol[1][3];

											total_err = abs(path_diff[1]);

											chosen = 2;
											error = 0;
										}
										else if ((cond_sol[0][0] == true && cond_sol[0][1] == true) && abs(path_diff[0]) + kaon_inv_mass_tmp[0] < abs(path_diff[1]) + kaon_inv_mass_tmp[1])
										{
											solution[0] = S.sol[0][0];
											solution[1] = S.sol[0][1];
											solution[2] = S.sol[0][2];
											solution[3] = S.sol[0][3];

											total_err = abs(path_diff[0]);

											chosen = 1;
											error = 0;
										}
										else
										{
											error = 1;
										}

										if (error == 0)
										{
											Knetri_tmp[0] = 0.;
											Knetri_tmp[1] = 0.;
											Knetri_tmp[2] = 0.;
											Knetri_tmp[3] = 0.;

											for (Int_t l = 0; l < 4; l++)
											{
												gamma_len_tmp[l] = sqrt(pow(cluster[0][ind_gam[l]] - solution[0], 2) +
																								pow(cluster[1][ind_gam[l]] - solution[1], 2) +
																								pow(cluster[2][ind_gam[l]] - solution[2], 2));

												gammatri_tmp[l][0] = cluster[4][ind_gam[l]] * (cluster[0][ind_gam[l]] - solution[0]) / gamma_len_tmp[l];
												gammatri_tmp[l][1] = cluster[4][ind_gam[l]] * (cluster[1][ind_gam[l]] - solution[1]) / gamma_len_tmp[l];
												gammatri_tmp[l][2] = cluster[4][ind_gam[l]] * (cluster[2][ind_gam[l]] - solution[2]) / gamma_len_tmp[l];
												gammatri_tmp[l][3] = cluster[4][ind_gam[l]];
												gammatri_tmp[l][4] = cluster[0][ind_gam[l]];
												gammatri_tmp[l][5] = cluster[1][ind_gam[l]];
												gammatri_tmp[l][6] = cluster[2][ind_gam[l]];
												gammatri_tmp[l][7] = cluster[3][ind_gam[l]];

												Knetri_tmp[0] += gammatri_tmp[l][0];
												Knetri_tmp[1] += gammatri_tmp[l][1];
												Knetri_tmp[2] += gammatri_tmp[l][2];
												Knetri_tmp[3] += gammatri_tmp[l][3];

												eqn_check[l] = abs(cVel * (cluster[3][ind_gam[l]] - solution[3]) - gamma_len_tmp[l]);
											}

											Knetri_tmp[4] = sqrt(pow(Knetri_tmp[0], 2) + pow(Knetri_tmp[1], 2) + pow(Knetri_tmp[2], 2));
											Knetri_tmp[5] = sqrt(pow(Knetri_tmp[3], 2) - pow(Knetri_tmp[0], 2) - pow(Knetri_tmp[1], 2) - pow(Knetri_tmp[2], 2));
											Knetri_tmp[6] = solution[0];
											Knetri_tmp[7] = solution[1];
											Knetri_tmp[8] = solution[2];
											Knetri_tmp[9] = solution[3];

											total_err += abs(Knetri_tmp[5] - mK0);
										}
										else
										{
											error = 3;
										}
									}

									if (total_err < total_err_def)
									{
										total_err_def = total_err;

										done = 1;

										fourg4taken[0] = ind_gam[0];
										fourg4taken[1] = ind_gam[1];
										fourg4taken[2] = ind_gam[2];
										fourg4taken[3] = ind_gam[3];

										Knetri[0] = 0.;
										Knetri[1] = 0.;
										Knetri[2] = 0.;
										Knetri[3] = 0.;

										for (Int_t l = 0; l < 4; l++)
										{
											gammatri[l][0] = gammatri_tmp[l][0];
											gammatri[l][1] = gammatri_tmp[l][1];
											gammatri[l][2] = gammatri_tmp[l][2];
											gammatri[l][3] = gammatri_tmp[l][3];
											gammatri[l][4] = gammatri_tmp[l][4];
											gammatri[l][5] = gammatri_tmp[l][5];
											gammatri[l][6] = gammatri_tmp[l][6];
											gammatri[l][7] = gammatri_tmp[l][7];

											Knetri[0] += gammatri[l][0];
											Knetri[1] += gammatri[l][1];
											Knetri[2] += gammatri[l][2];
											Knetri[3] += gammatri[l][3];
										}

										Knetri[4] = Knetri_tmp[4];
										Knetri[5] = Knetri_tmp[5];
										Knetri[6] = Knetri_tmp[6];
										Knetri[7] = Knetri_tmp[7];
										Knetri[8] = Knetri_tmp[8];
										Knetri[9] = Knetri_tmp[9];

										y_axis[0] = 0.;
										y_axis[1] = bhabha_mom[1];
										y_axis[2] = 0.;

										neu_vtx[0] = Knetri[6];
										neu_vtx[1] = Knetri[7];
										neu_vtx[2] = Knetri[8];

										plane_intersection(bhabha_vtx, y_axis, neu_vtx, Knetri, ip_tri); //! Plane rec

										if (abs(ip_tri[2] - bhabha_vtx[2]) > 2)
											ip_tri[2] = bhabha_vtx[2];

										ip_tri[0] = bhabha_vtx[0];
										ip_tri[1] = bhabha_vtx[1];
									}

									R.SetClu(0, 0, 0, 0, 0, 0);
									R.SetClu(1, 0, 0, 0, 0, 0);
									R.SetClu(2, 0, 0, 0, 0, 0);
									R.SetClu(3, 0, 0, 0, 0, 0);
									R.SetClu(4, 0, 0, 0, 0, 0);
									R.SetClu(5, 0, 0, 0, 0, 0);
								}
							}
			}
			else if (good_clus == 1)
			{
				ene_sum = 0.;
				error = 2;
				total_err = 999999.;

				ind_gam[0] = good_clus_ind[0];
				ind_gam[1] = good_clus_ind[1];
				ind_gam[2] = good_clus_ind[2];
				ind_gam[3] = good_clus_ind[3];

				cond_ene = cluster[4][ind_gam[0]] > MIN_CLU_ENE && cluster[4][ind_gam[1]] > MIN_CLU_ENE &&
									 cluster[4][ind_gam[2]] > MIN_CLU_ENE && cluster[4][ind_gam[3]] > MIN_CLU_ENE;

				cond_clus[0] = cluster[3][ind_gam[0]] > 0 && cluster[0][ind_gam[0]] != 0 && cluster[1][ind_gam[0]] != 0 && cluster[2][ind_gam[0]] != 0;
				cond_clus[1] = cluster[3][ind_gam[1]] > 0 && cluster[0][ind_gam[1]] != 0 && cluster[1][ind_gam[1]] != 0 && cluster[2][ind_gam[1]] != 0;
				cond_clus[2] = cluster[3][ind_gam[2]] > 0 && cluster[0][ind_gam[2]] != 0 && cluster[1][ind_gam[2]] != 0 && cluster[2][ind_gam[2]] != 0;
				cond_clus[3] = cluster[3][ind_gam[3]] > 0 && cluster[0][ind_gam[3]] != 0 && cluster[1][ind_gam[3]] != 0 && cluster[2][ind_gam[3]] != 0;

				if (cond_ene == true && cond_clus[0] && cond_clus[1] && cond_clus[2] && cond_clus[3])
				{
					for (Int_t l = 0; l < 4; l++)
					{
						R.SetClu(l, cluster[0][ind_gam[l]],
										 cluster[1][ind_gam[l]],
										 cluster[2][ind_gam[l]],
										 cluster[3][ind_gam[l]],
										 cluster[4][ind_gam[l]]);
					}

					R.SetClu(4, 0, 0, 0, 0, 0);
					R.SetClu(5, 0, 0, 0, 0, 0);

					S = R.MySolve(selected);

					if (S.error[0] == false || S.error[1] == false) // && cond_ene_sum == true)
					{

						kaon_mom[0][0] = 0.;
						kaon_mom[0][1] = 0.;
						kaon_mom[0][2] = 0.;
						kaon_mom[0][3] = 0.;

						kaon_mom[1][0] = 0.;
						kaon_mom[1][1] = 0.;
						kaon_mom[1][2] = 0.;
						kaon_mom[1][3] = 0.;

						error = 1;

						for (Int_t l1 = 0; l1 < 2; l1++)
						{
							for (Int_t l2 = 0; l2 < 4; l2++)
							{
								gamma_len[l1][l2] = sqrt(pow(cluster[0][ind_gam[l2]] - S.sol[l1][0], 2) +
																				 pow(cluster[1][ind_gam[l2]] - S.sol[l1][1], 2) +
																				 pow(cluster[2][ind_gam[l2]] - S.sol[l1][2], 2));

								gamma_mom[l1][l2][0] = cluster[4][ind_gam[l2]] * (cluster[0][ind_gam[l2]] - S.sol[l1][0]) / gamma_len[l1][l2];
								gamma_mom[l1][l2][1] = cluster[4][ind_gam[l2]] * (cluster[1][ind_gam[l2]] - S.sol[l1][1]) / gamma_len[l1][l2];
								gamma_mom[l1][l2][2] = cluster[4][ind_gam[l2]] * (cluster[2][ind_gam[l2]] - S.sol[l1][2]) / gamma_len[l1][l2];
								gamma_mom[l1][l2][3] = cluster[4][ind_gam[l2]];

								kaon_mom[l1][0] += gamma_mom[l1][l2][0];
								kaon_mom[l1][1] += gamma_mom[l1][l2][1];
								kaon_mom[l1][2] += gamma_mom[l1][l2][2];
								kaon_mom[l1][3] += gamma_mom[l1][l2][3];
							}

							y_axis[0] = 0.;
							y_axis[1] = bhabha_mom[1];
							y_axis[2] = 0.;

							plane_intersection(bhabha_vtx, y_axis, S.sol[l1], kaon_mom[l1], ip_tri); //! Plane rec

							if (abs(ip_tri[2] - bhabha_vtx[2]) > 2)
								ip_tri[2] = bhabha_vtx[2];

							ip_tri[0] = bhabha_vtx[0];
							ip_tri[1] = bhabha_vtx[1];

							kaon_vel[l1] = cVel * sqrt(pow(kaon_mom[l1][0], 2) + pow(kaon_mom[l1][1], 2) + pow(kaon_mom[l1][2], 2)) / (kaon_mom[l1][3]);
							kaon_len[l1] = sqrt(pow(S.sol[l1][0] - ip_tri[0], 2) + pow(S.sol[l1][1] - ip_tri[1], 2) + pow(S.sol[l1][2] - ip_tri[2], 2));
						}

						cond_sol[0][0] = (S.sol[0][3] < cluster[3][ind_gam[0]]) && (S.sol[0][3] < cluster[3][ind_gam[1]]) && (S.sol[0][3] < cluster[3][ind_gam[2]]) && (S.sol[0][3] < cluster[3][ind_gam[3]]);

						cond_sol[0][1] = (sqrt(pow(S.sol[0][0] - bhabha_vtx[0], 2) + pow(S.sol[0][1] - bhabha_vtx[1], 2)) < 200) && (abs(S.sol[0][2] - bhabha_vtx[2]) < 165);
						path_diff[0] = kaon_vel[0] * S.sol[0][3] - kaon_len[0];
						eqn_check_tmp[0][0] = cVel * (cluster[3][ind_gam[0]] - S.sol[0][3]) - gamma_len[0][0];
						eqn_check_tmp[0][1] = cVel * (cluster[3][ind_gam[1]] - S.sol[0][3]) - gamma_len[0][1];
						eqn_check_tmp[0][2] = cVel * (cluster[3][ind_gam[2]] - S.sol[0][3]) - gamma_len[0][2];
						eqn_check_tmp[0][3] = cVel * (cluster[3][ind_gam[3]] - S.sol[0][3]) - gamma_len[0][3];

						eqn_check_tmp_tot[0] = sqrt(pow(eqn_check_tmp[0][0], 2) +
																				pow(eqn_check_tmp[0][1], 2) +
																				pow(eqn_check_tmp[0][2], 2) +
																				pow(eqn_check_tmp[0][3], 2));

						cond_sol[1][0] = (S.sol[1][3] < cluster[3][ind_gam[0]]) && (S.sol[1][3] < cluster[3][ind_gam[1]]) && (S.sol[1][3] < cluster[3][ind_gam[2]]) && (S.sol[1][3] < cluster[3][ind_gam[3]]);
						cond_sol[1][1] = (sqrt(pow(S.sol[1][0] - bhabha_vtx[0], 2) + pow(S.sol[1][1] - bhabha_vtx[1], 2)) < 200) && (abs(S.sol[1][2] - bhabha_vtx[2]) < 165);
						path_diff[1] = kaon_vel[1] * S.sol[1][3] - kaon_len[1];
						eqn_check_tmp[1][0] = cVel * (cluster[3][ind_gam[0]] - S.sol[1][3]) - gamma_len[1][0];
						eqn_check_tmp[1][1] = cVel * (cluster[3][ind_gam[1]] - S.sol[1][3]) - gamma_len[1][1];
						eqn_check_tmp[1][2] = cVel * (cluster[3][ind_gam[2]] - S.sol[1][3]) - gamma_len[1][2];
						eqn_check_tmp[1][3] = cVel * (cluster[3][ind_gam[3]] - S.sol[1][3]) - gamma_len[1][3];

						eqn_check_tmp_tot[1] = sqrt(pow(eqn_check_tmp[1][0], 2) +
																				pow(eqn_check_tmp[1][1], 2) +
																				pow(eqn_check_tmp[1][2], 2) +
																				pow(eqn_check_tmp[1][3], 2));

						sol1[0] = S.sol[0][0];
						sol1[1] = S.sol[0][1];
						sol1[2] = S.sol[0][2];
						sol1[3] = S.sol[0][3];

						sol2[0] = S.sol[1][0];
						sol2[1] = S.sol[1][1];
						sol2[2] = S.sol[1][2];
						sol2[3] = S.sol[1][3];

						sol1err = path_diff[0];
						sol2err = path_diff[1];

						/*std::cout << "2. Script: " << std::endl;

						std::cout << S.sol[0][0] << " " << S.sol[0][1] << " " << S.sol[0][2] << " " << S.sol[0][3] << std::endl;
						std::cout << S.sol[1][0] << " " << S.sol[1][1] << " " << S.sol[1][2] << " " << S.sol[1][3] << std::endl;

						std::cout << eqn_check_tmp[0][0] << " " << eqn_check_tmp[1][0] << std::endl << std::endl;*/

						// We have to choose solutions

						if ((cond_sol[1][0] == true && cond_sol[1][1] == true) &&
										 abs(path_diff[1]) < abs(path_diff[0]))
						{
							solution[0] = S.sol[1][0];
							solution[1] = S.sol[1][1];
							solution[2] = S.sol[1][2];
							solution[3] = S.sol[1][3];

							total_err = abs(path_diff[1]);

							chosen = 2;
							error = 0;
						}
						else if ((cond_sol[0][0] == true && cond_sol[0][1] == true) &&
										 abs(path_diff[0]) < abs(path_diff[1]))
						{
							solution[0] = S.sol[0][0];
							solution[1] = S.sol[0][1];
							solution[2] = S.sol[0][2];
							solution[3] = S.sol[0][3];

							total_err = abs(path_diff[0]);

							chosen = 1;
							error = 0;
						}
						else
						{
							error = 1;
						}

						if (error == 0)
						{
							Knetri_tmp[0] = 0.;
							Knetri_tmp[1] = 0.;
							Knetri_tmp[2] = 0.;
							Knetri_tmp[3] = 0.;

							for (Int_t l = 0; l < 4; l++)
							{
								gamma_len_tmp[l] = sqrt(pow(cluster[0][ind_gam[l]] - solution[0], 2) +
																				pow(cluster[1][ind_gam[l]] - solution[1], 2) +
																				pow(cluster[2][ind_gam[l]] - solution[2], 2));

								gammatri_tmp[l][0] = cluster[4][ind_gam[l]] * (cluster[0][ind_gam[l]] - solution[0]) / gamma_len_tmp[l];
								gammatri_tmp[l][1] = cluster[4][ind_gam[l]] * (cluster[1][ind_gam[l]] - solution[1]) / gamma_len_tmp[l];
								gammatri_tmp[l][2] = cluster[4][ind_gam[l]] * (cluster[2][ind_gam[l]] - solution[2]) / gamma_len_tmp[l];
								gammatri_tmp[l][3] = cluster[4][ind_gam[l]];
								gammatri_tmp[l][4] = cluster[0][ind_gam[l]];
								gammatri_tmp[l][5] = cluster[1][ind_gam[l]];
								gammatri_tmp[l][6] = cluster[2][ind_gam[l]];
								gammatri_tmp[l][7] = cluster[3][ind_gam[l]];

								Knetri_tmp[0] += gammatri_tmp[l][0];
								Knetri_tmp[1] += gammatri_tmp[l][1];
								Knetri_tmp[2] += gammatri_tmp[l][2];
								Knetri_tmp[3] += gammatri_tmp[l][3];

								eqn_check[l] = cVel * (cluster[3][ind_gam[l]] - solution[3]) - gamma_len_tmp[l];
							}

							Knetri_tmp[4] = sqrt(pow(Knetri_tmp[0], 2) + pow(Knetri_tmp[1], 2) + pow(Knetri_tmp[2], 2));
							Knetri_tmp[5] = sqrt(pow(Knetri_tmp[3], 2) - pow(Knetri_tmp[0], 2) - pow(Knetri_tmp[1], 2) - pow(Knetri_tmp[2], 2));
							Knetri_tmp[6] = solution[0];
							Knetri_tmp[7] = solution[1];
							Knetri_tmp[8] = solution[2];
							Knetri_tmp[9] = solution[3];

							total_err += abs(Knetri_tmp[5] - mK0);
															// + pow(eqn_check[0], 2) +
															//   pow(eqn_check[1], 2) +
															//   pow(eqn_check[2], 2) +
															//   pow(eqn_check[3], 2));
						}
						else
						{
							error = 3;
						}
					}

					if (total_err < total_err_def)
					{
						total_err_def = total_err;

						done = 1;

						fourg4taken[0] = ind_gam[0];
						fourg4taken[1] = ind_gam[1];
						fourg4taken[2] = ind_gam[2];
						fourg4taken[3] = ind_gam[3];

						Knetri[0] = 0.;
						Knetri[1] = 0.;
						Knetri[2] = 0.;
						Knetri[3] = 0.;

						for (Int_t l = 0; l < 4; l++)
						{
							gammatri[l][0] = gammatri_tmp[l][0];
							gammatri[l][1] = gammatri_tmp[l][1];
							gammatri[l][2] = gammatri_tmp[l][2];
							gammatri[l][3] = gammatri_tmp[l][3];
							gammatri[l][4] = gammatri_tmp[l][4];
							gammatri[l][5] = gammatri_tmp[l][5];
							gammatri[l][6] = gammatri_tmp[l][6];
							gammatri[l][7] = gammatri_tmp[l][7];

							Knetri[0] += gammatri[l][0];
							Knetri[1] += gammatri[l][1];
							Knetri[2] += gammatri[l][2];
							Knetri[3] += gammatri[l][3];
						}

						Knetri[4] = Knetri_tmp[4];
						Knetri[5] = Knetri_tmp[5];
						Knetri[6] = Knetri_tmp[6];
						Knetri[7] = Knetri_tmp[7];
						Knetri[8] = Knetri_tmp[8];
						Knetri[9] = Knetri_tmp[9];

						y_axis[0] = 0.;
						y_axis[1] = bhabha_mom[1];
						y_axis[2] = 0.;

						neu_vtx[0] = Knetri[6];
						neu_vtx[1] = Knetri[7];
						neu_vtx[2] = Knetri[8];

						plane_intersection(bhabha_vtx, y_axis, neu_vtx, Knetri, ip_tri); //! Plane rec

						if (abs(ip_tri[2] - bhabha_vtx[2]) > 2)
							ip_tri[2] = bhabha_vtx[2];

						ip_tri[0] = bhabha_vtx[0];
						ip_tri[1] = bhabha_vtx[1];
					}

					R.SetClu(0, 0, 0, 0, 0, 0);
					R.SetClu(1, 0, 0, 0, 0, 0);
					R.SetClu(2, 0, 0, 0, 0, 0);
					R.SetClu(3, 0, 0, 0, 0, 0);
					R.SetClu(4, 0, 0, 0, 0, 0);
					R.SetClu(5, 0, 0, 0, 0, 0);
				}
			}
			else
			{
				ene_sum = 0.;
				error = 2;
				total_err = 999999.;

				cond_ene = pgammc[0][3] > MIN_CLU_ENE && pgammc[1][3] > MIN_CLU_ENE &&
									 pgammc[2][3] > MIN_CLU_ENE && pgammc[3][3] > MIN_CLU_ENE;

				cond_clus[0] = pgammc[0][7] > 0 && pgammc[0][4] != 0 && pgammc[0][5] != 0 && pgammc[0][6] != 0;
				cond_clus[1] = pgammc[1][7] > 0 && pgammc[1][4] != 0 && pgammc[1][5] != 0 && pgammc[1][6] != 0;
				cond_clus[2] = pgammc[2][7] > 0 && pgammc[2][4] != 0 && pgammc[2][5] != 0 && pgammc[2][6] != 0;
				cond_clus[3] = pgammc[3][7] > 0 && pgammc[3][4] != 0 && pgammc[3][5] != 0 && pgammc[3][6] != 0;

				if (cond_ene == true && cond_clus[0] && cond_clus[1] && cond_clus[2] && cond_clus[3])
				{
					for (Int_t l = 0; l < 4; l++)
					{
						R.SetClu(l, pgammc[l][4],
										 pgammc[l][5],
										 pgammc[l][6],
										 pgammc[l][7],
										 pgammc[l][3]);
					}

					R.SetClu(4, 0, 0, 0, 0, 0);
					R.SetClu(5, 0, 0, 0, 0, 0);

					S = R.MySolve(selected);

					if (S.error[0] == false || S.error[1] == false) // && cond_ene_sum == true)
					{

						kaon_mom[0][0] = 0.;
						kaon_mom[0][1] = 0.;
						kaon_mom[0][2] = 0.;
						kaon_mom[0][3] = 0.;

						kaon_mom[1][0] = 0.;
						kaon_mom[1][1] = 0.;
						kaon_mom[1][2] = 0.;
						kaon_mom[1][3] = 0.;

						error = 1;

						for (Int_t l1 = 0; l1 < 2; l1++)
						{
							for (Int_t l2 = 0; l2 < 4; l2++)
							{
								gamma_len[l1][l2] = sqrt(pow(pgammc[l2][4] - S.sol[l1][0], 2) +
																				 pow(pgammc[l2][5] - S.sol[l1][1], 2) +
																				 pow(pgammc[l2][6] - S.sol[l1][2], 2));

								gamma_mom[l1][l2][0] = pgammc[l2][3] * (pgammc[l2][4] - S.sol[l1][0]) / gamma_len[l1][l2];
								gamma_mom[l1][l2][1] = pgammc[l2][3] * (pgammc[l2][5] - S.sol[l1][1]) / gamma_len[l1][l2];
								gamma_mom[l1][l2][2] = pgammc[l2][3] * (pgammc[l2][6] - S.sol[l1][2]) / gamma_len[l1][l2];
								gamma_mom[l1][l2][3] = pgammc[l2][3];

								kaon_mom[l1][0] += gamma_mom[l1][l2][0];
								kaon_mom[l1][1] += gamma_mom[l1][l2][1];
								kaon_mom[l1][2] += gamma_mom[l1][l2][2];
								kaon_mom[l1][3] += gamma_mom[l1][l2][3];
							}

							y_axis[0] = 0.;
							y_axis[1] = bhabha_mom[1];
							y_axis[2] = 0.;

							plane_intersection(bhabha_vtx, y_axis, S.sol[l1], kaon_mom[l1], ip_tri); //! Plane rec

							if (abs(ip_tri[2] - bhabha_vtx[2]) > 2)
								ip_tri[2] = bhabha_vtx[2];

							ip_tri[0] = bhabha_vtx[0];
							ip_tri[1] = bhabha_vtx[1];

							kaon_vel[l1] = cVel * sqrt(pow(kaon_mom[l1][0], 2) + pow(kaon_mom[l1][1], 2) + pow(kaon_mom[l1][2], 2)) / (kaon_mom[l1][3]);
							kaon_len[l1] = sqrt(pow(S.sol[l1][0] - ip_tri[0], 2) + pow(S.sol[l1][1] - ip_tri[1], 2) + pow(S.sol[l1][2] - ip_tri[2], 2));
						}

						std::cout << pgammc[0][7] << " " << pgammc[1][7] << " " << pgammc[2][7] << " " << pgammc[3][7] << std::endl; 

						cond_sol[0][0] = 1;//(S.sol[0][3] < pgammc[0][7]) && (S.sol[0][3] < pgammc[1][7]) && (S.sol[0][3] < pgammc[2][7]) && (S.sol[0][3] < pgammc[3][7]);

						cond_sol[0][1] = (sqrt(pow(S.sol[0][0] - bhabha_vtx[0], 2) + pow(S.sol[0][1] - bhabha_vtx[1], 2)) < 200) && (abs(S.sol[0][2] - bhabha_vtx[2]) < 165);
						path_diff[0] = kaon_vel[0] * S.sol[0][3] - kaon_len[0];
						eqn_check_tmp[0][0] = cVel * (pgammc[0][7] - S.sol[0][3]) - gamma_len[0][0];
						eqn_check_tmp[0][1] = cVel * (pgammc[1][7] - S.sol[0][3]) - gamma_len[0][1];
						eqn_check_tmp[0][2] = cVel * (pgammc[2][7] - S.sol[0][3]) - gamma_len[0][2];
						eqn_check_tmp[0][3] = cVel * (pgammc[3][7] - S.sol[0][3]) - gamma_len[0][3];

						eqn_check_tmp_tot[0] = sqrt(pow(eqn_check_tmp[0][0], 2) +
																				pow(eqn_check_tmp[0][1], 2) +
																				pow(eqn_check_tmp[0][2], 2) +
																				pow(eqn_check_tmp[0][3], 2));

						cond_sol[1][0] = 1;//(S.sol[1][3] < pgammc[0][7]) && (S.sol[1][3] < pgammc[1][7]) && (S.sol[1][3] < pgammc[2][7]) && (S.sol[1][3] < pgammc[3][7]);
						cond_sol[1][1] = (sqrt(pow(S.sol[1][0] - bhabha_vtx[0], 2) + pow(S.sol[1][1] - bhabha_vtx[1], 2)) < 200) && (abs(S.sol[1][2] - bhabha_vtx[2]) < 165);
						path_diff[1] = kaon_vel[1] * S.sol[1][3] - kaon_len[1];
						eqn_check_tmp[1][0] = cVel * (pgammc[0][7] - S.sol[1][3]) - gamma_len[1][0];
						eqn_check_tmp[1][1] = cVel * (pgammc[1][7] - S.sol[1][3]) - gamma_len[1][1];
						eqn_check_tmp[1][2] = cVel * (pgammc[2][7] - S.sol[1][3]) - gamma_len[1][2];
						eqn_check_tmp[1][3] = cVel * (pgammc[3][7] - S.sol[1][3]) - gamma_len[1][3];

						eqn_check_tmp_tot[1] = sqrt(pow(eqn_check_tmp[1][0], 2) +
																				pow(eqn_check_tmp[1][1], 2) +
																				pow(eqn_check_tmp[1][2], 2) +
																				pow(eqn_check_tmp[1][3], 2));

						sol1[0] = S.sol[0][0];
						sol1[1] = S.sol[0][1];
						sol1[2] = S.sol[0][2];
						sol1[3] = S.sol[0][3];

						sol2[0] = S.sol[1][0];
						sol2[1] = S.sol[1][1];
						sol2[2] = S.sol[1][2];
						sol2[3] = S.sol[1][3];

						sol1err = path_diff[0];
						sol2err = path_diff[1];

						/*std::cout << "2. Script: " << std::endl;

						std::cout << S.sol[0][0] << " " << S.sol[0][1] << " " << S.sol[0][2] << " " << S.sol[0][3] << std::endl;
						std::cout << S.sol[1][0] << " " << S.sol[1][1] << " " << S.sol[1][2] << " " << S.sol[1][3] << std::endl;

						std::cout << eqn_check_tmp[0][0] << " " << eqn_check_tmp[1][0] << std::endl << std::endl;*/

						// We have to choose solutions

						if ((cond_sol[1][0] == true && cond_sol[1][1] == true) &&
										 abs(path_diff[1]) < abs(path_diff[0]))
						{
							solution[0] = S.sol[1][0];
							solution[1] = S.sol[1][1];
							solution[2] = S.sol[1][2];
							solution[3] = S.sol[1][3];

							total_err = abs(path_diff[1]);

							chosen = 2;
							error = 0;
						}
						else if ((cond_sol[0][0] == true && cond_sol[0][1] == true) &&
										 abs(path_diff[0]) < abs(path_diff[1]))
						{
							solution[0] = S.sol[0][0];
							solution[1] = S.sol[0][1];
							solution[2] = S.sol[0][2];
							solution[3] = S.sol[0][3];

							total_err = abs(path_diff[0]);

							chosen = 1;
							error = 0;
						}
						else
						{
							error = 1;
						}

						if (error == 0)
						{
							Knetri_tmp[0] = 0.;
							Knetri_tmp[1] = 0.;
							Knetri_tmp[2] = 0.;
							Knetri_tmp[3] = 0.;

							for (Int_t l = 0; l < 4; l++)
							{
								gamma_len_tmp[l] = sqrt(pow(pgammc[l][4] - solution[0], 2) +
																				pow(pgammc[l][5] - solution[1], 2) +
																				pow(pgammc[l][6] - solution[2], 2));

								gammatri_tmp[l][0] = pgammc[l][3] * (pgammc[l][4] - solution[0]) / gamma_len_tmp[l];
								gammatri_tmp[l][1] = pgammc[l][3] * (pgammc[l][5] - solution[1]) / gamma_len_tmp[l];
								gammatri_tmp[l][2] = pgammc[l][3] * (pgammc[l][6] - solution[2]) / gamma_len_tmp[l];
								gammatri_tmp[l][3] = pgammc[l][3];
								gammatri_tmp[l][4] = pgammc[l][4];
								gammatri_tmp[l][5] = pgammc[l][5];
								gammatri_tmp[l][6] = pgammc[l][6];
								gammatri_tmp[l][7] = pgammc[l][7];

								Knetri_tmp[0] += gammatri_tmp[l][0];
								Knetri_tmp[1] += gammatri_tmp[l][1];
								Knetri_tmp[2] += gammatri_tmp[l][2];
								Knetri_tmp[3] += gammatri_tmp[l][3];

								eqn_check[l] = cVel * (pgammc[l][7] - solution[3]) - gamma_len_tmp[l];
							}

							Knetri_tmp[4] = sqrt(pow(Knetri_tmp[0], 2) + pow(Knetri_tmp[1], 2) + pow(Knetri_tmp[2], 2));
							Knetri_tmp[5] = sqrt(pow(Knetri_tmp[3], 2) - pow(Knetri_tmp[0], 2) - pow(Knetri_tmp[1], 2) - pow(Knetri_tmp[2], 2));
							Knetri_tmp[6] = solution[0];
							Knetri_tmp[7] = solution[1];
							Knetri_tmp[8] = solution[2];
							Knetri_tmp[9] = solution[3];

							total_err += abs(Knetri_tmp[5] - mK0);
															// + pow(eqn_check[0], 2) +
															//   pow(eqn_check[1], 2) +
															//   pow(eqn_check[2], 2) +
															//   pow(eqn_check[3], 2));
						}
						else
						{
							error = 3;
						}
					}

					if (total_err < total_err_def)
					{
						total_err_def = total_err;

						done = 1;

						fourg4taken[0] = ind_gam[0];
						fourg4taken[1] = ind_gam[1];
						fourg4taken[2] = ind_gam[2];
						fourg4taken[3] = ind_gam[3];

						Knetri[0] = 0.;
						Knetri[1] = 0.;
						Knetri[2] = 0.;
						Knetri[3] = 0.;

						for (Int_t l = 0; l < 4; l++)
						{
							gammatri[l][0] = gammatri_tmp[l][0];
							gammatri[l][1] = gammatri_tmp[l][1];
							gammatri[l][2] = gammatri_tmp[l][2];
							gammatri[l][3] = gammatri_tmp[l][3];
							gammatri[l][4] = gammatri_tmp[l][4];
							gammatri[l][5] = gammatri_tmp[l][5];
							gammatri[l][6] = gammatri_tmp[l][6];
							gammatri[l][7] = gammatri_tmp[l][7];

							Knetri[0] += gammatri[l][0];
							Knetri[1] += gammatri[l][1];
							Knetri[2] += gammatri[l][2];
							Knetri[3] += gammatri[l][3];
						}

						Knetri[4] = Knetri_tmp[4];
						Knetri[5] = Knetri_tmp[5];
						Knetri[6] = Knetri_tmp[6];
						Knetri[7] = Knetri_tmp[7];
						Knetri[8] = Knetri_tmp[8];
						Knetri[9] = Knetri_tmp[9];

						y_axis[0] = 0.;
						y_axis[1] = bhabha_mom[1];
						y_axis[2] = 0.;

						neu_vtx[0] = Knetri[6];
						neu_vtx[1] = Knetri[7];
						neu_vtx[2] = Knetri[8];

						plane_intersection(bhabha_vtx, y_axis, neu_vtx, Knetri, ip_tri); //! Plane rec

						if (abs(ip_tri[2] - bhabha_vtx[2]) > 2)
							ip_tri[2] = bhabha_vtx[2];

						ip_tri[0] = bhabha_vtx[0];
						ip_tri[1] = bhabha_vtx[1];
					}

					R.SetClu(0, 0, 0, 0, 0, 0);
					R.SetClu(1, 0, 0, 0, 0, 0);
					R.SetClu(2, 0, 0, 0, 0, 0);
					R.SetClu(3, 0, 0, 0, 0, 0);
					R.SetClu(4, 0, 0, 0, 0, 0);
					R.SetClu(5, 0, 0, 0, 0, 0);
				}
			}
		}

		tree->Fill();
	}

	tree->Print();

	tree->Write();
	file->Close();
	delete file;
}
