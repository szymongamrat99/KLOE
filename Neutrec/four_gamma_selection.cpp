#include "TFile.h"
#include "TTree.h"
#include "TH1.h"

#include "../../Include/Codes/reconstructor.h"
#include "../../Include/const.h"
#include "chain_init.C"

void four_gamma_selection(UInt_t filenumber = 1, TString directory = "230531_data", TString rootname = "data_stream42_")
{
    TString treedir = "INTERF", treename = "h1", fulltree = "",
            extension = ".root", fullname = "";

    fullname = "../../ROOT_files/" + directory + "/backup/" + rootname + filenumber + extension;
    fulltree = treedir + "/" + treename;

    TFile *file = new TFile(fullname, "UPDATE");
    TTree *tree = (TTree*)file->Get(fulltree);

    //Branches' addresses
    //Bhabha vars
    Float_t bhabha_mom[4], bhabha_vtx[3];

    tree->SetBranchAddress("Bpx", &bhabha_mom[0]);
    tree->SetBranchAddress("Bpy", &bhabha_mom[1]);
    tree->SetBranchAddress("Bpz", &bhabha_mom[2]);
    tree->SetBranchAddress("Broots", &bhabha_mom[3]);

    tree->SetBranchAddress("Bx", &bhabha_vtx[0]);
    tree->SetBranchAddress("By", &bhabha_vtx[1]);
    tree->SetBranchAddress("Bz", &bhabha_vtx[2]);

    //Cluster vars
    Int_t nclu;
    UChar_t mctruth;
    Float_t cluster[5][500], Kchboost[9], Knereclor[9], Knemc[9];

    tree->SetBranchAddress("nclu", &nclu);
    tree->SetBranchAddress("Xcl", cluster[0]);
    tree->SetBranchAddress("Ycl", cluster[1]);
    tree->SetBranchAddress("Zcl", cluster[2]);
    tree->SetBranchAddress("Tcl", cluster[3]);
    tree->SetBranchAddress("Enecl", cluster[4]);

    Int_t nentries = (Int_t)tree->GetEntries();

    //Variables
    Bool_t cond_ene, cond_ene_sum, cond_sol[2][2], onesol;
    Int_t ind_gam[4], error = 0, index = 0;
    Float_t ene_sum, gamma_mom[2][6][4], gamma_len[2][6],
            kaon_mom[2][4], kaon_vel[2], kaon_len[2], path_diff[2], 
            sol_tmp[4], total_err = 0., total_err_def = 0.,
            gamma_len_tmp[4], kaon_mom_tmp[4], kaon_inv_mass, ene_sum_tmp;
    Double_t solution[4];
    Float_t gammatri_tmp[4][8], Knetri_tmp[10];

    int selected[4] = {1, 2, 3, 4};

    Reconstructor R;
    Solution S;

    Int_t done, fourg4taken[4];
    Float_t gammatri[4][8], Knetri[10], totalerr;

    TBranch *b_gamma1tri = tree->Branch("fourgamma1tri", gammatri[0], "fourgamma1tri[8]/F");
    TBranch *b_gamma2tri = tree->Branch("fourgamma2tri", gammatri[1], "fourgamma2tri[8]/F");
    TBranch *b_gamma3tri = tree->Branch("fourgamma3tri", gammatri[2], "fourgamma3tri[8]/F");
    TBranch *b_gamma4tri = tree->Branch("fourgamma4tri", gammatri[3], "fourgamma4tri[8]/F");

    TBranch *b_Knetri = tree->Branch("fourKnetri", Knetri, "fourKnetri[10]/F");

    TBranch *b_done = tree->Branch("done4", &done, "done4/I");

    TBranch *b_fourg4taken = tree->Branch("fourg4taken", fourg4taken, "fourg4taken[4]/I");


    Float_t distance = 0., distance_mc = 0., distance_diff = 0.;
    Float_t mc_dist[90],sigmas[90];

    for(Int_t i = 0; i < nentries; i++)
    {
        tree->GetEntry(i);

        for(Int_t j1 = 0; j1 < 4; j1++)
            for(Int_t j2 = 0; j2 < 8; j2++) gammatri[j1][j2] = -999.;

        for(Int_t j = 0; j < 10; j++) Knetri[j] = -999.;

        for(Int_t j = 0; j < 4; j++) fourg4taken[j] = -999;

        totalerr = -999.;

        done = 0;

        total_err_def = 999999.;

        if(nclu >= 4)
        {
            for(Int_t j1 = 0; j1 < nclu - 3; j1++)
                for(Int_t j2 = j1 + 1; j2 < nclu - 2; j2++)
                    for(Int_t j3 = j2 + 1; j3 < nclu - 1; j3++)
                        for(Int_t j4 = j3 + 1; j4 < nclu; j4++)
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

                            for(Int_t l = 0; l < 4; l++)
                            {
                                R.SetClu(l, cluster[0][ind_gam[l]],
                                            cluster[1][ind_gam[l]], 
                                            cluster[2][ind_gam[l]],
                                            cluster[3][ind_gam[l]], 
                                            cluster[4][ind_gam[l]]);
                            }

                            R.SetClu(4, 0, 0, 0, 0, 0);
                            R.SetClu(5, 0, 0, 0, 0, 0);

                            cond_ene_sum = R.CombEnergy(selected) > 350 && R.CombEnergy(selected) < 700;                               

                            S = R.MySolve(selected);

                            if((S.error[0] == false || S.error[1] == false) && cond_ene == true && cond_ene_sum == true)
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

                                for(Int_t l1 = 0; l1 < 2; l1++)
                                {
                                    for(Int_t l2 = 0; l2 < 4; l2++)
                                    {
                                        gamma_len[l1][l2] = sqrt(pow(cluster[0][ind_gam[l2]] - S.sol[l1][0],2) +
                                                                 pow(cluster[1][ind_gam[l2]] - S.sol[l1][1],2) +
                                                                 pow(cluster[2][ind_gam[l2]] - S.sol[l1][2],2) );

                                        gamma_mom[l1][l2][0] = cluster[4][ind_gam[l2]]*(cluster[0][ind_gam[l2]] - S.sol[l1][0])/gamma_len[l1][l2];
                                        gamma_mom[l1][l2][1] = cluster[4][ind_gam[l2]]*(cluster[1][ind_gam[l2]] - S.sol[l1][1])/gamma_len[l1][l2];
                                        gamma_mom[l1][l2][2] = cluster[4][ind_gam[l2]]*(cluster[2][ind_gam[l2]] - S.sol[l1][2])/gamma_len[l1][l2];
                                        gamma_mom[l1][l2][3] = cluster[4][ind_gam[l2]];

                                        kaon_mom[l1][0] += gamma_mom[l1][l2][0];
                                        kaon_mom[l1][1] += gamma_mom[l1][l2][1];
                                        kaon_mom[l1][2] += gamma_mom[l1][l2][2];
                                        kaon_mom[l1][3] += gamma_mom[l1][l2][3];
                                    }
                                                            
                                    kaon_vel[l1] = c_vel*sqrt(pow(kaon_mom[l1][0],2) + pow(kaon_mom[l1][1],2) + pow(kaon_mom[l1][2],2))/(kaon_mom[l1][3]); 
                                    kaon_len[l1] = sqrt(pow(S.sol[l1][0] - bhabha_vtx[0],2) + pow(S.sol[l1][1] - bhabha_vtx[1],2) + pow(S.sol[l1][2] - bhabha_vtx[2],2));
                                }

                                cond_sol[0][0] = (S.sol[0][3] >= 0) && (S.sol[0][3] < 60);
                                cond_sol[0][1] = (sqrt(pow(S.sol[0][0],2) + pow(S.sol[0][1],2)) < 200) && (abs(S.sol[0][2]) < 169);
                                path_diff[0] = kaon_vel[0]*S.sol[0][3] - kaon_len[0];

                                cond_sol[1][0] = (S.sol[1][3] >= 0) && (S.sol[1][3] < 60);
                                cond_sol[1][1] = (sqrt(pow(S.sol[1][0],2) + pow(S.sol[1][1],2)) < 200) && (abs(S.sol[1][2]) < 169);
                                path_diff[1] = kaon_vel[1]*S.sol[1][3] - kaon_len[1];

                                onesol = (abs(S.sol[0][3] - S.sol[1][3]) < 1);

                                //We have to choose solutions

                                if(onesol == 1 && cond_sol[0][0] == 1 && cond_sol[0][1] == 1
                                               && cond_sol[1][0] == 1 && cond_sol[1][1] == 1 && 
                                               abs(path_diff[0]) < 10. && abs(path_diff[1]) < 10.)
                                {
                                    solution[0] = (S.sol[0][0] + S.sol[1][0])/2.;
                                    solution[1] = (S.sol[0][1] + S.sol[1][1])/2.;
                                    solution[2] = (S.sol[0][2] + S.sol[1][2])/2.;
                                    solution[3] = (S.sol[0][3] + S.sol[1][3])/2.;


                                    error = 0;
                                }
                                else if(onesol != 1 && (cond_sol[0][0] == 1 && cond_sol[0][1] == 1) &&
                                         abs(path_diff[0]) < 10. && abs(path_diff[1]) > 10.)
                                {
                                    solution[0] = S.sol[0][0];
                                    solution[1] = S.sol[0][1];
                                    solution[2] = S.sol[0][2];
                                    solution[3] = S.sol[0][3];


                                    error = 0;
                                }
                                else if(onesol != 1 && (cond_sol[1][0] == 1 && cond_sol[1][1] == 1) &&
                                        abs(path_diff[1]) < 10. && abs(path_diff[0]) > 10.)
                                {
                                    solution[0] = S.sol[1][0];
                                    solution[1] = S.sol[1][1];
                                    solution[2] = S.sol[1][2];
                                    solution[3] = S.sol[1][3];

                                    error = 0;
                                }
                                else
                                {
                                    error = 1;
                                }

                                if( error == 0 )
                                {
                                    Knetri_tmp[0] = 0.;
                                    Knetri_tmp[1] = 0.;
                                    Knetri_tmp[2] = 0.;
                                    Knetri_tmp[3] = 0.;

                                    for(Int_t l = 0; l < 4; l++)
                                    {
                                        gamma_len_tmp[l] = sqrt(pow(cluster[0][ind_gam[l]] - solution[0],2) +
                                                                pow(cluster[1][ind_gam[l]] - solution[1],2) +
                                                                pow(cluster[2][ind_gam[l]] - solution[2],2) );

                                        gammatri_tmp[l][0] = cluster[4][ind_gam[l]]*(cluster[0][ind_gam[l]] - solution[0])/gamma_len_tmp[l];
                                        gammatri_tmp[l][1] = cluster[4][ind_gam[l]]*(cluster[1][ind_gam[l]] - solution[1])/gamma_len_tmp[l];
                                        gammatri_tmp[l][2] = cluster[4][ind_gam[l]]*(cluster[2][ind_gam[l]] - solution[2])/gamma_len_tmp[l];
                                        gammatri_tmp[l][3] = cluster[4][ind_gam[l]];
                                        gammatri_tmp[l][4] = cluster[0][ind_gam[l]];
                                        gammatri_tmp[l][5] = cluster[1][ind_gam[l]];
                                        gammatri_tmp[l][6] = cluster[2][ind_gam[l]];
                                        gammatri_tmp[l][7] = cluster[3][ind_gam[l]];

                                        Knetri_tmp[0] += gammatri_tmp[l][0];
                                        Knetri_tmp[1] += gammatri_tmp[l][1];
                                        Knetri_tmp[2] += gammatri_tmp[l][2];
                                        Knetri_tmp[3] += gammatri_tmp[l][3];
                                    }

                                        Knetri_tmp[4] = sqrt(pow(Knetri_tmp[0],2) + pow(Knetri_tmp[1],2) + pow(Knetri_tmp[2],2));
                                        Knetri_tmp[5] = sqrt(pow(Knetri_tmp[3],2) - pow(Knetri_tmp[0],2) - pow(Knetri_tmp[1],2) - pow(Knetri_tmp[2],2));
                                        Knetri_tmp[6] = solution[0];
                                        Knetri_tmp[7] = solution[1];
                                        Knetri_tmp[8] = solution[2];
                                        Knetri_tmp[9] = solution[3];

                                        total_err = abs(Knetri_tmp[5] - m_k0);

                                }
                                else
                                {
                                    error = 3;
                                }

                            }

                            if(total_err < total_err_def)
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

                                for(Int_t l = 0; l < 4; l++)
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


                            } 
                                                    

                        }         
        }

        b_gamma1tri->Fill();
        b_gamma2tri->Fill();
        b_gamma3tri->Fill();
        b_gamma4tri->Fill();

        b_Knetri->Fill();

        b_done->Fill();

        b_fourg4taken->Fill();

    }


    file->cd(treedir);
    tree->Write();
    file->Close();
    delete file;


}
