#include "../../Include/const.h"
#include "../../Include/Codes/charged_mom.h"

void doublepipi_selection(UInt_t filenumber = 1, TString directory = "230531_data", TString rootname = "data_stream42_")
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

    //Charged vars
    Float_t charged_vtx[3][MaxNumVtx], charged_vars[3][MaxNumtrkv];
    Int_t ntv, nv, qualv[MaxNumVtx];
    UInt_t iv[MaxNumtrkv], vtx[MaxNumVtx];

    tree->SetBranchAddress("ntv", &ntv);
    tree->SetBranchAddress("nv", &nv);

    tree->SetBranchAddress("xv", charged_vtx[0]);
    tree->SetBranchAddress("yv", charged_vtx[1]);
    tree->SetBranchAddress("zv", charged_vtx[2]);

    tree->SetBranchAddress("Curv", charged_vars[0]);
    tree->SetBranchAddress("Phiv", charged_vars[1]);
    tree->SetBranchAddress("Cotv", charged_vars[2]);

    tree->SetBranchAddress("qualv", qualv);
    tree->SetBranchAddress("iv", iv);
    tree->SetBranchAddress("vtx", vtx);

    UInt_t nentries = tree->GetEntries();

    Int_t vtakenpi1[3], vtakenpi2[3], done[2];

    TBranch *b_donepipi = tree->Branch("donepipi", done, "donepipi[2]/I");

    TBranch *b_vtakenpi1 = tree->Branch("vtakenpi1", vtakenpi1, "vtakenpi1[3]/I");
    TBranch *b_vtakenpi2 = tree->Branch("vtakenpi2", vtakenpi2, "vtakenpi2[3]/I");

    Float_t trkpi11[4], trkpi21[4], trkpi12[4], trkpi22[4], Kchrec1[9], Kchrec2[9];

    TBranch *b_trkpi11 = tree->Branch("trkpi11", trkpi11, "trkpi11[4]/F");
    TBranch *b_trkpi12 = tree->Branch("trkpi12", trkpi12, "trkpi12[4]/F");
    TBranch *b_trkpi21 = tree->Branch("trkpi21", trkpi21, "trkpi21[4]/F");
    TBranch *b_trkpi22 = tree->Branch("trkpi22", trkpi22, "trkpi22[4]/F");

    TBranch *b_kchrec1 = tree->Branch("Kchrec1", Kchrec1, "Kchrec1[9]/F");
    TBranch *b_kchrec2 = tree->Branch("Kchrec2", Kchrec2, "Kchrec2[9]/F");

    Int_t sign[2];
    Float_t mom_vec[2][4], radius, height, kaon_mom[2][9], inv_mass, inv_mass_diff, inv_mass_diff_def1, inv_mass_diff_def2;

    for(Int_t i = 0; i < nentries; i++)
    {
        tree->GetEntry(i);

        inv_mass_diff_def1 = 9999999.;
        inv_mass_diff_def2 = 9999999.;

        done[0] = 0;
        done[1] = 0;

        for(Int_t j = 0; j < 4; j++)
        {
            trkpi11[j] = -999.;
            trkpi12[j] = -999.;
            trkpi21[j] = -999.;
            trkpi22[j] = -999.;
        }

        for(Int_t j = 0; j < 9; j++)
        {
            Kchrec1[j] = -999.;
            Kchrec2[j] = -999.;
        }

        for(Int_t j = 0; j < 3; j++)
        {
            vtakenpi1[j] = -999;
            vtakenpi2[j] = -999;
        }

        if(nv >= 2)
        {
            for(Int_t j = 0; j < nv; j++)
            {
                radius = sqrt(pow(charged_vtx[0][j] - bhabha_vtx[0],2) + pow(charged_vtx[1][j] - bhabha_vtx[1],2));
                height = abs(charged_vtx[2][j] - bhabha_vtx[2]);

                if( qualv[j] == 1 && vtx[j] == 2 && radius < 10 && height < 20 )
                {

                    for(Int_t k1 = 0; k1 < ntv - 1; k1++)
                    {
                        
                        sign[0] = charged_vars[0][k1]/abs(charged_vars[0][k1]);
                        if(iv[k1]-1 == j) 
                            charged_mom(charged_vars[0][k1], charged_vars[1][k1], charged_vars[2][k1], mom_vec[0], 1);

                        for(Int_t k2 = k1 + 1; k2 < ntv; k2++)
                        {
                            sign[1] = charged_vars[0][k2]/abs(charged_vars[0][k2]);
                            if(iv[k2]-1 == j && sign[0] != sign[1])
                            { 
                                charged_mom(charged_vars[0][k2], charged_vars[1][k2], charged_vars[2][k2], mom_vec[1], 1);
                                inv_mass = sqrt(pow(mom_vec[0][3] + mom_vec[1][3],2) - 
                                                pow(mom_vec[0][0] + mom_vec[1][0],2) - 
                                                pow(mom_vec[0][1] + mom_vec[1][1],2) -
                                                pow(mom_vec[0][2] + mom_vec[1][2],2) );
                                inv_mass_diff = abs(inv_mass - mK0);
                            }
                            else
                            {
                                inv_mass_diff = 9999999.;
                            }

                            if(inv_mass_diff < inv_mass_diff_def1 || vtakenpi1[0] == -999)
                            {
                                inv_mass_diff_def1 = inv_mass_diff;

                                done[0] = 1;

                                vtakenpi1[0] = j;
                                vtakenpi1[1] = k1;
                                vtakenpi1[2] = k2;

                                trkpi11[0] = mom_vec[0][0];
                                trkpi11[1] = mom_vec[0][1];
                                trkpi11[2] = mom_vec[0][2];
                                trkpi11[3] = mom_vec[0][3];

                                trkpi12[0] = mom_vec[1][0];
                                trkpi12[1] = mom_vec[1][1];
                                trkpi12[2] = mom_vec[1][2];
                                trkpi12[3] = mom_vec[1][3];

                                Kchrec1[0] = trkpi11[0] + trkpi12[0];
                                Kchrec1[1] = trkpi11[1] + trkpi12[1];
                                Kchrec1[2] = trkpi11[2] + trkpi12[2];
                                Kchrec1[3] = trkpi11[3] + trkpi12[3];
                                Kchrec1[4] = sqrt(pow(trkpi11[0] + trkpi12[0],2) + 
                                                pow(trkpi11[1] + trkpi12[1],2) +
                                                pow(trkpi11[2] + trkpi12[2],2) );
                                Kchrec1[5] = sqrt(pow(Kchrec1[3] ,2) -
                                                pow(Kchrec1[4], 2) );
                                Kchrec1[6] = charged_vtx[0][vtakenpi1[0]];
                                Kchrec1[7] = charged_vtx[1][vtakenpi1[0]];
                                Kchrec1[8] = charged_vtx[2][vtakenpi1[0]];
                            }
                        }
                    }
                    
                }
            }

            for(Int_t j = 0; j < nv; j++)
            {

                if( qualv[j] == 1 && vtx[j] == 2 && vtakenpi1[0] != j )
                {

                    for(Int_t k1 = 0; k1 < ntv - 1; k1++)
                    {
                        
                        sign[0] = charged_vars[0][k1]/abs(charged_vars[0][k1]);
                        if(iv[k1]-1 == j) 
                            charged_mom(charged_vars[0][k1], charged_vars[1][k1], charged_vars[2][k1], mom_vec[0], 1);

                        for(Int_t k2 = k1 + 1; k2 < ntv; k2++)
                        {
                            sign[1] = charged_vars[0][k2]/abs(charged_vars[0][k2]);
                            if(iv[k2]-1 == j && sign[0] != sign[1])
                            { 
                                charged_mom(charged_vars[0][k2], charged_vars[1][k2], charged_vars[2][k2], mom_vec[1], 1);
                                inv_mass = sqrt(pow(mom_vec[0][3] + mom_vec[1][3],2) - 
                                                pow(mom_vec[0][0] + mom_vec[1][0],2) - 
                                                pow(mom_vec[0][1] + mom_vec[1][1],2) -
                                                pow(mom_vec[0][2] + mom_vec[1][2],2) );
                                inv_mass_diff = abs(inv_mass - mK0);
                            }
                            else
                            {
                                inv_mass_diff = 9999999.;
                            }

                            if(inv_mass_diff < inv_mass_diff_def2 || vtakenpi2[0] == -999)
                            {
                                inv_mass_diff_def2 = inv_mass_diff;

                                done[1] = 1;

                                vtakenpi2[0] = j;
                                vtakenpi2[1] = k1;
                                vtakenpi2[2] = k2;

                                trkpi21[0] = mom_vec[0][0];
                                trkpi21[1] = mom_vec[0][1];
                                trkpi21[2] = mom_vec[0][2];
                                trkpi21[3] = mom_vec[0][3];

                                trkpi22[0] = mom_vec[1][0];
                                trkpi22[1] = mom_vec[1][1];
                                trkpi22[2] = mom_vec[1][2];
                                trkpi22[3] = mom_vec[1][3];

                                Kchrec2[0] = trkpi21[0] + trkpi22[0];
                                Kchrec2[1] = trkpi21[1] + trkpi22[1];
                                Kchrec2[2] = trkpi21[2] + trkpi22[2];
                                Kchrec2[3] = trkpi21[3] + trkpi22[3];
                                Kchrec2[4] = sqrt(pow(trkpi21[0] + trkpi22[0],2) + 
                                                pow(trkpi21[1] + trkpi22[1],2) +
                                                pow(trkpi21[2] + trkpi22[2],2) );
                                Kchrec2[5] = sqrt(pow(Kchrec2[3] ,2) -
                                                pow(Kchrec2[4], 2) );
                                Kchrec2[6] = charged_vtx[0][vtakenpi2[0]];
                                Kchrec2[7] = charged_vtx[1][vtakenpi2[0]];
                                Kchrec2[8] = charged_vtx[2][vtakenpi2[0]];
                            }
                        }
                    }
                    
                }
            }

        }

        b_donepipi->Fill();
        b_kchrec1->Fill();
        b_kchrec2->Fill();
        b_trkpi11->Fill();
        b_trkpi12->Fill();
        b_trkpi21->Fill();
        b_trkpi22->Fill();
        b_vtakenpi1->Fill();
        b_vtakenpi2->Fill();

    }

    file->cd(treedir);
    tree->Write();
    file->Close();
    delete file;

}
