#include "../../Include/Codes/reconstructor.h"
#include "../../Include/const.h"

void six_gamma_selection(UInt_t filenumber = 1, TString directory = "230623_mc", TString rootname = "mc_stream62_mccard2_")
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
    Float_t cluster[5][200], Knemc[9];

    tree->SetBranchAddress("nclu", &nclu);
    tree->SetBranchAddress("Xcl", cluster[0]);
    tree->SetBranchAddress("Ycl", cluster[1]);
    tree->SetBranchAddress("Zcl", cluster[2]);
    tree->SetBranchAddress("Tcl", cluster[3]);
    tree->SetBranchAddress("Enecl", cluster[4]);

    Int_t nentries = (Int_t)tree->GetEntries();

    //Variables
    Bool_t cond_ene, cond_ene_sum, cond_sol[2][2], onesol;
    Int_t indgam[6], selected[4], error = 0, index = 0;
    Float_t ene_sum, gamma_mom[2][6][4], gamma_len[2][6],
            kaon_mom[2][4], kaon_vel[2], kaon_len[2], path_diff[2], 
            res_err[15], sol_tmp[4], total_err = 0., total_err_def = 0.,
            gamma_mom_tmp[6][4], gamma_len_tmp[6], kaon_mom_tmp[4], kaon_inv_mass,
            ene_sum_tmp;
    Double_t solution[15][4];

    Reconstructor R;
    Solution S;

    Int_t done, g6taken[6];
    Float_t gammatri[6][8], Knetri[10], totalerr;

    TBranch *b_gamma1tri = tree->Branch("gamma1tri", gammatri[0], "gamma1tri[8]/F");
    TBranch *b_gamma2tri = tree->Branch("gamma2tri", gammatri[1], "gamma2tri[8]/F");
    TBranch *b_gamma3tri = tree->Branch("gamma3tri", gammatri[2], "gamma3tri[8]/F");
    TBranch *b_gamma4tri = tree->Branch("gamma4tri", gammatri[3], "gamma4tri[8]/F");
    TBranch *b_gamma5tri = tree->Branch("gamma5tri", gammatri[4], "gamma5tri[8]/F");
    TBranch *b_gamma6tri = tree->Branch("gamma6tri", gammatri[5], "gamma6tri[8]/F");

    TBranch *b_Knetri = tree->Branch("Knetri", Knetri, "Knetri[10]/F");

    TBranch *b_done = tree->Branch("done", &done, "done/I");
    TBranch *b_totalerr = tree->Branch("totalerr", &totalerr, "totalerr/F");

    TBranch *b_g6taken = tree->Branch("g6taken", g6taken, "g6taken[6]/I");

    Float_t distance = 0., distance_mc = 0., distance_diff = 0.;
    Float_t mc_dist[90],sigmas[90];

    /*TH2 *hist2d = new TH2F("hist2d", "hist2d", 11, 0, 50., 5, -50, 50.);
    TH1 *hist1d = new TH1F("hist1d", "hist1d", 11, 0, 50.);
    TH1 *hist = new TH1F("hist", "hist", 200, 0, 7000.);*/

    for(Int_t i = 0; i < nentries; i++)
    {
        tree->GetEntry(i);

        for(Int_t j1 = 0; j1 < 6; j1++)
            for(Int_t j2 = 0; j2 < 8; j2++) gammatri[j1][j2] = -999.;

        for(Int_t j = 0; j < 10; j++) Knetri[j] = -999.;

        for(Int_t j = 0; j < 6; j++) g6taken[j] = -999;

        totalerr = -999.;
        done = 0;

        if(nclu >= 6)
        {
            for(Int_t j1 = 0; j1 < nclu - 5; j1++)
                for(Int_t j2 = j1 + 1; j2 < nclu - 4; j2++)
                    for(Int_t j3 = j2 + 1; j3 < nclu - 3; j3++)
                        for(Int_t j4 = j3 + 1; j4 < nclu - 2; j4++)
                            for(Int_t j5 = j4 + 1; j5 < nclu - 1; j5++)
                                for(Int_t j6 = j5 + 1; j6 < nclu; j6++)
                                {

                                    ene_sum = 0.;
                                    index = 0;
                                    total_err = 999999.;
                                    error = 2;

                                    total_err_def = 999999.;

                                    indgam[0] = j1;
                                    indgam[1] = j2;
                                    indgam[2] = j3;
                                    indgam[3] = j4;
                                    indgam[4] = j5;
                                    indgam[5] = j6;

                                    cond_ene = cluster[4][indgam[0]] > MIN_CLU_ENE && cluster[4][indgam[1]] > MIN_CLU_ENE &&
                                               cluster[4][indgam[2]] > MIN_CLU_ENE && cluster[4][indgam[3]] > MIN_CLU_ENE &&
                                               cluster[4][indgam[4]] > MIN_CLU_ENE && cluster[4][indgam[5]] > MIN_CLU_ENE;

                                    for(Int_t k = 0; k < 6; k++) ene_sum += cluster[4][indgam[k]];
                                    cond_ene_sum = ene_sum > 350 && ene_sum < 700;                               

                                    if(cond_ene == true && cond_ene_sum == true)
                                    {
                                        for(Int_t l = 0; l < 6; l++) R.SetClu(l, cluster[0][indgam[l]],
                                                                            cluster[1][indgam[l]], cluster[2][indgam[l]],
                                                                            cluster[3][indgam[l]], cluster[4][indgam[l]]);

                                        for(Int_t k1 = 0; k1 < 3; k1++)
                                            for(Int_t k2 = k1 + 1; k2 < 4; k2++)
                                                for(Int_t k3 = k2 + 1; k3 < 5; k3++)
                                                    for(Int_t k4 = k3 + 1; k4 < 6; k4++)
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

                                                        selected[0] = k1 + 1;
                                                        selected[1] = k2 + 1;
                                                        selected[2] = k3 + 1;
                                                        selected[3] = k4 + 1;

                                                        S = R.MySolve(selected);

                                                        for(Int_t l1 = 0; l1 < 2; l1++){
                                                            for(Int_t l2 = 0; l2 < 6; l2++)
                                                            {
                                                                gamma_len[l1][l2] = sqrt(pow(cluster[0][indgam[l2]] - S.sol[l1][0],2) +
                                                                                         pow(cluster[1][indgam[l2]] - S.sol[l1][1],2) +
                                                                                         pow(cluster[2][indgam[l2]] - S.sol[l1][2],2) );

                                                                gamma_mom[l1][l2][0] = cluster[4][indgam[l2]]*(cluster[0][indgam[l2]] - S.sol[l1][0])/gamma_len[l1][l2];
                                                                gamma_mom[l1][l2][1] = cluster[4][indgam[l2]]*(cluster[1][indgam[l2]] - S.sol[l1][1])/gamma_len[l1][l2];
                                                                gamma_mom[l1][l2][2] = cluster[4][indgam[l2]]*(cluster[2][indgam[l2]] - S.sol[l1][2])/gamma_len[l1][l2];
                                                                gamma_mom[l1][l2][3] = cluster[4][indgam[l2]];

                                                                kaon_mom[l1][0] += gamma_mom[l1][l2][0];
                                                                kaon_mom[l1][1] += gamma_mom[l1][l2][1];
                                                                kaon_mom[l1][2] += gamma_mom[l1][l2][2];
                                                                kaon_mom[l1][3] += gamma_mom[l1][l2][3];
                                                            }
                                                            
                                                            kaon_vel[l1] = cVel*sqrt(pow(kaon_mom[l1][0],2) + pow(kaon_mom[l1][1],2) + pow(kaon_mom[l1][2],2))/(kaon_mom[l1][3]); 
                                                            kaon_len[l1] = sqrt(pow(S.sol[l1][0] - bhabha_vtx[0],2) + pow(S.sol[l1][1] - bhabha_vtx[1],2) + pow(S.sol[l1][2] - bhabha_vtx[2],2));

                                                            cond_sol[l1][0] = (S.sol[l1][3] >= 0) && (S.sol[l1][3] < 60);
                                                            cond_sol[l1][1] = (sqrt(pow(S.sol[l1][0],2) + pow(S.sol[l1][1],2)) < 200) && (abs(S.sol[l1][2]) < 169);
                                                            path_diff[l1] = kaon_vel[l1]*S.sol[l1][3] - kaon_len[l1];
                                                        }

                                                        onesol = (abs(S.sol[0][3] - S.sol[1][3]) < 1);

                                                        //We have to choose solutions

                                                        if(onesol == 1 && cond_sol[0][0] == 1 && cond_sol[0][1] == 1
                                                                       && cond_sol[1][0] == 1 && cond_sol[1][1] == 1 &&
                                                                       abs(path_diff[0]) < 0.1 && abs(path_diff[1]) < 0.1)
                                                        {
                                                            solution[index][0] = (S.sol[0][0] + S.sol[1][0])/2.;
                                                            solution[index][1] = (S.sol[0][1] + S.sol[1][1])/2.;
                                                            solution[index][2] = (S.sol[0][2] + S.sol[1][2])/2.;
                                                            solution[index][3] = (S.sol[0][3] + S.sol[1][3])/2.;

                                                            res_err[index] = R.ResidualErrTot(solution[index]);

                                                            index++;

                                                            error = 0;
                                                        }
                                                        else if((cond_sol[0][0] == 1 && cond_sol[0][1] == 1) &&
                                                             abs(path_diff[0]) < 0.1 && abs(path_diff[1]) > 0.1)
                                                        {
                                                            solution[index][0] = S.sol[0][0];
                                                            solution[index][1] = S.sol[0][1];
                                                            solution[index][2] = S.sol[0][2];
                                                            solution[index][3] = S.sol[0][3];

                                                            res_err[index] = R.ResidualErrTot(solution[index]);

                                                            index++;

                                                            error = 0;
                                                        }
                                                        else if((cond_sol[1][0] == 1 && cond_sol[1][1] == 1) &&
                                                            abs(path_diff[1]) < 0.1 && abs(path_diff[0]) > 0.1)
                                                        {
                                                            solution[index][0] = S.sol[1][0];
                                                            solution[index][1] = S.sol[1][1];
                                                            solution[index][2] = S.sol[1][2];
                                                            solution[index][3] = S.sol[1][3];

                                                            res_err[index] = R.ResidualErrTot(solution[index]);

                                                            index++;

                                                            error = 0;
                                                        }
                                                        else
                                                        {
                                                            error = 1;
                                                        }

                                                    }

                                                    if(index > 0)
                                                    {
                                                        total_err = 0.;
                                                        for(Int_t l = 0; l < index; l++) total_err += res_err[index];

                                                        if( total_err < total_err_def )
                                                        {
                                                            total_err_def = total_err;

                                                            totalerr = total_err_def;

                                                            g6taken[0] = indgam[0];
                                                            g6taken[1] = indgam[1];
                                                            g6taken[2] = indgam[2];
                                                            g6taken[3] = indgam[3];
                                                            g6taken[4] = indgam[4];
                                                            g6taken[5] = indgam[5];

                                                            sol_tmp[0] = 0.;
                                                            sol_tmp[1] = 0.;
                                                            sol_tmp[2] = 0.;
                                                            sol_tmp[3] = 0.;

                                                            Knetri[0] = 0.;
                                                            Knetri[1] = 0.;
                                                            Knetri[2] = 0.;
                                                            Knetri[3] = 0.;

                                                            done = 1;

                                                            for(Int_t l = 0; l < index; l++)
                                                            {
                                                                sol_tmp[0] += solution[l][0]/index;
                                                                sol_tmp[1] += solution[l][1]/index;
                                                                sol_tmp[2] += solution[l][2]/index;
                                                                sol_tmp[3] += solution[l][3]/index;
                                                            }



                                                            for(Int_t l = 0; l < 6; l++)
                                                            {
                                                                gamma_len_tmp[l] = sqrt(pow(cluster[0][indgam[l]] - sol_tmp[0],2) +
                                                                                            pow(cluster[1][indgam[l]] - sol_tmp[1],2) +
                                                                                            pow(cluster[2][indgam[l]] - sol_tmp[2],2) );

                                                                gammatri[l][0] = cluster[4][indgam[l]]*(cluster[0][indgam[l]] - sol_tmp[0])/gamma_len_tmp[l];
                                                                gammatri[l][1] = cluster[4][indgam[l]]*(cluster[1][indgam[l]] - sol_tmp[1])/gamma_len_tmp[l];
                                                                gammatri[l][2] = cluster[4][indgam[l]]*(cluster[2][indgam[l]] - sol_tmp[2])/gamma_len_tmp[l];
                                                                gammatri[l][3] = cluster[4][indgam[l]];
                                                                gammatri[l][4] = cluster[0][indgam[l]];
                                                                gammatri[l][5] = cluster[1][indgam[l]];
                                                                gammatri[l][6] = cluster[2][indgam[l]];
                                                                gammatri[l][7] = cluster[3][indgam[l]];

                                                                Knetri[0] += gammatri[l][0];
                                                                Knetri[1] += gammatri[l][1];
                                                                Knetri[2] += gammatri[l][2];
                                                                Knetri[3] += gammatri[l][3];
                                                            }

                                                                Knetri[4] = sqrt(pow(Knetri[0],2) + pow(Knetri[1],2) + pow(Knetri[2],2));
                                                                Knetri[5] = sqrt(pow(Knetri[3],2) - pow(Knetri[0],2) - pow(Knetri[1],2) - pow(Knetri[2],2));
                                                                Knetri[6] = sol_tmp[0];
                                                                Knetri[7] = sol_tmp[1];
                                                                Knetri[8] = sol_tmp[2];
                                                                Knetri[9] = sol_tmp[3];

                                                        }
                                                    }

                                    }
                                    else
                                    {
                                        error = 2;
                                    }

                                }
            
        }

        /*if(total_err_def < 999999.)
        {   
            distance = sqrt(pow(Knetri[6] - bhabha_vtx[0],2) + pow(Knetri[7] - bhabha_vtx[1],2));
            distance_mc = sqrt(pow(Knemc[6] - bhabha_vtx[0],2) + pow(Knemc[7] - bhabha_vtx[1],2));
            distance_diff = sqrt(pow(Knetri[6] - Knemc[6],2) + pow(Knetri[7] - Knemc[7],2));

            hist->Fill(Knetri[5]);
            hist2d->Fill(distance_mc, distance_diff);    
        } */

        b_gamma1tri->Fill();
        b_gamma2tri->Fill();
        b_gamma3tri->Fill();
        b_gamma4tri->Fill();
        b_gamma5tri->Fill();
        b_gamma6tri->Fill();

        b_Knetri->Fill();

        b_done->Fill();
        b_totalerr->Fill();

        b_g6taken->Fill();

    }

    /*TFitResultPtr r;

    for(Int_t i = 1; i <= 11; i++)
    {
        r = hist2d->ProjectionY("_py", i, i)->Fit("gaus","S");
        hist1d->SetBinContent(i, r->Parameter(2));
        hist1d->SetBinError(i, r->ParError(2));
    }*/

    //hist1d->Draw();

    file->cd(treedir);
    tree->Write();
    file->Close();
    delete file;


}
