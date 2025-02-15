// Author: Szymon Gamrat
// Date of last update: 03.02.2025

#include <boost/optional.hpp>

#include "../inc/kchrec.hpp"

int kchrec_Kmass(TChain &chain, Controls::DataType &dataType, ErrorHandling::ErrorLogs &logger, KLOE::pm00 &Obj)
{

  // Structure from const.h to ease navigation
  BaseKinematics baseKin;

  // Bhabha variables - avg per run
  chain.SetBranchAddress("Bx", &baseKin.bhabha_vtx[0]);
  chain.SetBranchAddress("By", &baseKin.bhabha_vtx[1]);
  chain.SetBranchAddress("Bz", &baseKin.bhabha_vtx[2]);

  chain.SetBranchAddress("Bpx", &baseKin.phi_mom[0]);
  chain.SetBranchAddress("Bpy", &baseKin.phi_mom[1]);
  chain.SetBranchAddress("Bpz", &baseKin.phi_mom[2]);
  chain.SetBranchAddress("Broots", &baseKin.phi_mom[3]);
  // -----------------------------------------------------------

  // Track - vtx info
  chain.SetBranchAddress("nv", &baseKin.nv);
  chain.SetBranchAddress("ntv", &baseKin.ntv);
  chain.SetBranchAddress("iv", baseKin.iv);
  // -----------------------------------------------------------

  // Track properties
  chain.SetBranchAddress("Curv", baseKin.Curv);
  chain.SetBranchAddress("Phiv", baseKin.Phiv);
  chain.SetBranchAddress("Cotv", baseKin.Cotv);
  // -----------------------------------------------------------

  // Vertex position
  chain.SetBranchAddress("xv", baseKin.xv);
  chain.SetBranchAddress("yv", baseKin.yv);
  chain.SetBranchAddress("zv", baseKin.zv);
  // -----------------------------------------------------------

  // mcflag and mctruth position
  chain.SetBranchAddress("mcflag", &baseKin.mcflag);
  chain.SetBranchAddress("mctruth", &baseKin.mctruth);
  // -----------------------------------------------------------

  // Creation of a ttree with new variables
  TFile *file = new TFile("invmass_test_pdg.root", "RECREATE");
  TTree *tree = new TTree("h1", "Test of pdg const file");

  Float_t
      invMass = 0.;

  TBranch *b = tree->Branch("invMass", &invMass, "invMass/F");

  Int_t nentries = chain.GetEntries();

  Int_t mode = 1; // Model for pi+pi-

  // Initialization of Charged part of decay reconstruction class
  // Constructor is below, in the loop
  boost::optional<KLOE::ChargedVtxRec<>> eventAnalysis;
  // -------------------------------------------------------------

  for (Int_t i = 0; i < nentries; i++)
  {
    chain.GetEntry(i);

    // Finding Ks from KSKL->pi+pi-pi+pi-
    int
        findKl = 0,
        findKs = 1,
        findClose = 0;

    int
        last_vtx = 0;

    // Set the proper data type to be analyzed
    Bool_t data_flag = false;
    Int_t mctruth_int = int(baseKin.mctruth), mcflag_int = int(baseKin.mcflag);
    Obj.dataFlagSetter(dataType, data_flag, mcflag_int, mctruth_int);
    // -------------------------------------------------------------------

    if(data_flag)
    {
			// Construction of the charged rec class object
			eventAnalysis.emplace(baseKin.nv, baseKin.ntv, baseKin.iv, baseKin.bhabha_vtx, baseKin.Curv, baseKin.Phiv, baseKin.Cotv, baseKin.xv, baseKin.yv, baseKin.zv, mode);

      // KMASS HYPOTHESIS
      eventAnalysis->findKchRec(baseKin.Kchrec, baseKin.trk[0], baseKin.trk[1], baseKin.vtaken, baseKin.errFlag);
      // ------------------------------------------------------------------

      // KSL HYPOTHESIS
      eventAnalysis->findKchRec(baseKin.Kchrec, baseKin.trk[0], baseKin.trk[1], baseKin.vtaken, baseKin.errFlag);

      eventAnalysis->findKchRec(baseKin.Kchrec, baseKin.trk[0], baseKin.trk[1], baseKin.vtaken, baseKin.errFlag);
      // ------------------------------------------------------------------

      // CLOSEST TO IP
			eventAnalysis->findKClosestRec(baseKin.Kchrec, baseKin.trk[0], baseKin.trk[1], baseKin.vtaken, baseKin.errFlag);
			// ------------------------------------------------------------------

    }

    tree->Fill();
  }

  tree->Print();

  file->Write();
  file->Close();

  return 0;
}