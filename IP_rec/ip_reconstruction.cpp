#include "../../Include/Codes/plane_intersection.h"
#include "../../Include/Codes/closest_approach.h"

void ip_reconstruction(UInt_t filenumber = 1, TString directory = "data_root_files_vtx_part", 
                                                TString rootname = "data_stream42_" )
{

    TString treedir = "INTERF", treename = "h1", fulltree = "",
    extension = ".root", fullname = "";

    fullname = "../../ROOT_files/" + directory + "/" + rootname + filenumber + extension;
    fulltree = treedir + "/" + treename;

    TFile *file = new TFile(fullname, "UPDATE");
    TTree *tree = (TTree*)file->Get(fulltree);

    //Branches' addresses
    //Bhabha vars
    UInt_t mctruth;
    Float_t bhabha_mom[3], bhabha_vtx[3], Kchboost[9], ipmc[3];

    tree->SetBranchAddress("Bpx", &bhabha_mom[0]);
    tree->SetBranchAddress("Bpy", &bhabha_mom[1]);
    tree->SetBranchAddress("Bpz", &bhabha_mom[2]);

    tree->SetBranchAddress("Bx", &bhabha_vtx[0]);
    tree->SetBranchAddress("By", &bhabha_vtx[1]);
    tree->SetBranchAddress("Bz", &bhabha_vtx[2]);

    tree->SetBranchAddress("Kchboost", Kchboost);

    Int_t nentries = (Int_t)tree->GetEntries();

    Float_t beamline[3] = {0., 0., 1.}, kaon_vtx[3], kaon_mom[3], ip_plane[3], ip_closest[3];

    TBranch *b_ipplane = tree->Branch("ip_plane", ip_plane, "ip_plane[3]/F");
    TBranch *b_ipclosest = tree->Branch("ip_closest", ip_closest, "ip_closest[3]/F");

    for(Int_t i = 0; i < nentries; i++)
    {
        tree->GetEntry(i);

        for(Int_t j = 0; j < 3; j++)
        {
            kaon_mom[j] = Kchboost[j];
            kaon_vtx[j] = Kchboost[6 + j];
        }

        plane_intersection(bhabha_vtx, bhabha_mom, beamline,
                           kaon_vtx, kaon_mom,
                           ip_plane);

        closest_approach(bhabha_vtx, beamline,
                         kaon_vtx, kaon_mom,
                         ip_closest);

        b_ipplane->Fill();
        b_ipclosest->Fill();
        
    }

    file->cd(treedir);
    tree->Write();
    file->Close();
    delete file;

}