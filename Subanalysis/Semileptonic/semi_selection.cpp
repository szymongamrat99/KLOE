#include "../../Include/const.h"

void semi_selection(UInt_t filenumber = 1, TString directory = "230531_data", TString rootname = "data_stream42_")
{
    TString treedir = "INTERF", treename = "h1", fulltree = "",
            extension = ".root", fullname = "";

    fullname = "../../ROOT_files/" + directory + "/backup/" + rootname + filenumber + extension;
    fulltree = treedir + "/" + treename;

    TFile *file = new TFile(fullname, "UPDATE");
    TTree *tree = (TTree*)file->Get(fulltree);

    Float_t bhabha_mom[4], ip[3];
    tree->SetBranchAddress("Bpx", &bhabha_mom[0]);
    tree->SetBranchAddress("Bpy", &bhabha_mom[1]);
    tree->SetBranchAddress("Bpz", &bhabha_mom[2]);
    tree->SetBranchAddress("Broots", &bhabha_mom[3]);

    tree->SetBranchAddress("ip", ip);

    Float_t Kchboost[9], Kchrec[9], trk1[4], trk2[4];
    tree->SetBranchAddress("trk1", trk1);
    tree->SetBranchAddress("trk2", trk2);
    tree->SetBranchAddress("Kchboost", Kchboost);
    tree->SetBranchAddress("Kchrec", Kchrec);



    Int_t nentries = tree->GetEntries();

    Float_t Kchanother[4], kch_length_LAB, pich1_momlength, pich2_momlength, nominator, denominator;

    Float_t Qmiss_inv = 0, anglepipi_CM_kch = 0;

    TBranch *b_qmiss = tree->Branch("Qmiss_inv", &Qmiss_inv, "Qmiss_inv/F");
    TBranch *b_angle = tree->Branch("anglepipi_CM_kch", &anglepipi_CM_kch, "anglepipi_CM_kch/F");

    TLorentzVector phi_mom, pich1_mom, pich2_mom, kchanother_mom;

    TVector3 boost_phi, boost_kaon;
    
    for(Int_t i = 0; i < nentries; i++)
    {
        tree->GetEntry(i);

        Qmiss_inv = -999.;
        anglepipi_CM_kch = -999.;

        if(1)
        {
            //Lorentz invariant quantity
            Qmiss_inv = sqrt( pow(Kchboost[3] - Kchrec[3], 2) -
                        pow(Kchboost[0] - Kchrec[0], 2) -
                        pow(Kchboost[1] - Kchrec[1], 2) -
                        pow(Kchboost[2] - Kchrec[2], 2) );

            //Another method to calculate charged kaon's momentum
            //Calculation using boosted variables, so semileptonic should be worse than others

            kch_length_LAB = sqrt( pow(Kchboost[6] - ip[0],2) + pow(Kchboost[7] - ip[1],2) + pow(Kchboost[8] - ip[2],2) );

            Kchanother[0] = sqrt(pow(Kchboost[3],2) - pow(mK0,2))*(Kchboost[6] - ip[0])/kch_length_LAB;
            Kchanother[1] = sqrt(pow(Kchboost[3],2) - pow(mK0,2))*(Kchboost[7] - ip[1])/kch_length_LAB;
            Kchanother[2] = sqrt(pow(Kchboost[3],2) - pow(mK0,2))*(Kchboost[8] - ip[2])/kch_length_LAB;
            Kchanother[3] = Kchboost[3];

            //Initialization of Lorentz vectors

            phi_mom(0) = bhabha_mom[0];
            phi_mom(1) = bhabha_mom[1];
            phi_mom(2) = bhabha_mom[2];
            phi_mom(3) = bhabha_mom[3];

            pich1_mom(0) = trk1[0];
            pich1_mom(1) = trk1[1];
            pich1_mom(2) = trk1[2];
            pich1_mom(3) = trk1[3];

            pich2_mom(0) = trk2[0];
            pich2_mom(1) = trk2[1];
            pich2_mom(2) = trk2[2];
            pich2_mom(3) = trk2[3];

            kchanother_mom(0) = Kchanother[0];
            kchanother_mom(1) = Kchanother[1];
            kchanother_mom(2) = Kchanother[2];
            kchanother_mom(3) = Kchanother[3];

            boost_phi(0) = -phi_mom(0)/phi_mom(3);
            boost_phi(1) = -phi_mom(1)/phi_mom(3);
            boost_phi(2) = -phi_mom(2)/phi_mom(3);

            //Doing boost to CM phi
            phi_mom.Boost(boost_phi);
            pich1_mom.Boost(boost_phi);
            pich2_mom.Boost(boost_phi);
            kchanother_mom.Boost(boost_phi);

            //Doing boost to CM kch
            boost_kaon(0) = -kchanother_mom(0)/kchanother_mom(3);
            boost_kaon(1) = -kchanother_mom(1)/kchanother_mom(3);
            boost_kaon(2) = -kchanother_mom(2)/kchanother_mom(3);

            phi_mom.Boost(boost_kaon);
            pich1_mom.Boost(boost_kaon);
            pich2_mom.Boost(boost_kaon);
            kchanother_mom.Boost(boost_kaon);

            pich1_momlength = sqrt(pow(pich1_mom(0),2)+pow(pich1_mom(1),2)+pow(pich1_mom(2),2));
            pich2_momlength = sqrt(pow(pich2_mom(0),2)+pow(pich2_mom(1),2)+pow(pich2_mom(2),2));

            denominator = pich1_momlength*pich2_momlength;
            nominator = pich1_mom(0)*pich2_mom(0) + pich1_mom(1)*pich2_mom(1) + pich1_mom(2)*pich2_mom(2);

            anglepipi_CM_kch = 180.*acos( nominator/denominator )/M_PI;
        }

        b_qmiss->Fill();
        b_angle->Fill();

    }

    file->cd(treedir);
    tree->Write();
    file->Close();
    delete file;

}
