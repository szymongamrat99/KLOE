#define four_pi_cxx
// The class definition in four_pi.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.

// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// root> T->Process("four_pi.C")
// root> T->Process("four_pi.C","some options")
// root> T->Process("four_pi.C+")
//

#include "four_pi.h"
#include <TH2.h>
#include <TStyle.h>
#include <HistManager.h>
#include <interf_function.h>
#include <StatisticalCutter.h>
#include <kloe_class.h>

HistManager *histMgr;
StatisticalCutter *cutter;
HistManager::HistConfig invMassKneConfig;
HistManager::HistConfig invMassKneMCConfig;
HistManager::HistConfig invMassKchConfig;
HistManager::HistConfig invMassKchMCConfig;
HistManager::HistConfig timeDiffConfig;
HistManager::HistConfig timeDiffMCConfig;

Float_t EmissKS = 0.0;
Float_t EmissKL = 0.0;
Float_t PmissKS = 0.0;
Float_t PmissKL = 0.0;
Float_t KchrecKSMom = 0.0;
Float_t KchrecKLMom = 0.0;

KLOE::pm00 *Obj;

void four_pi::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

   histMgr = new HistManager(channNum, channColor, channNames, kFullCircle, kBlack, kOrange);

   std::string cutFileName = "/data/ssd/gamrat/KLOE/Subanalysis/Properties/cut-limits-final.json";

   cutter = new StatisticalCutter(cutFileName, 7, KLOE::HypothesisCode::FOUR_PI);

   Obj = new KLOE::pm00();

   Float_t pKTwoBody = Obj->TwoBodyDecayMass(mPhi, mK0, mK0);

   ///////////////////////////////////////////////////////////////////
   cutter->RegisterVariableGetter("InvMassKch", [&]()
                                  { return KchrecKS[5]; });
   cutter->RegisterCentralValueGetter("InvMassKch", [&]()
                                      { return mK0; });
   ///////////////////////////////////////////////////////////////////
   cutter->RegisterVariableGetter("InvMassKne", [&]()
                                  { return KchrecKL[5]; });
   cutter->RegisterCentralValueGetter("InvMassKne", [&]()
                                      { return mK0; });
   ///////////////////////////////////////////////////////////////////
   cutter->RegisterVariableGetter("TwoBodyMomKS", [&]()
                                  { return KchrecKSMom; });
   cutter->RegisterCentralValueGetter("TwoBodyMomKS", [&]()
                                      { return pKTwoBody; });
   ///////////////////////////////////////////////////////////////////
   cutter->RegisterVariableGetter("TwoBodyMomKL", [&]()
                                  { return KchrecKLMom; });
   cutter->RegisterCentralValueGetter("TwoBodyMomKL", [&]()
                                      { return pKTwoBody; });
   ///////////////////////////////////////////////////////////////////
   cutter->RegisterVariableGetter("MissTotKS", [&]()
                                  { return sqrt(pow(PmissKS, 2) + pow(EmissKS, 2)); });
   ///////////////////////////////////////////////////////////////////
   cutter->RegisterVariableGetter("MissHigherKS", [&]()
                                  { return (pow(EmissKS, 2) - pow(PmissKS, 2)); });
   cutter->RegisterVariableGetter("MissLowerKS", [&]()
                                  { return (pow(EmissKS, 2) - pow(PmissKS, 2)); });
   ///////////////////////////////////////////////////////////////////
   cutter->RegisterVariableGetter("MissTotKL", [&]()
                                  { return sqrt(pow(PmissKL, 2) + pow(EmissKL, 2)); });
   ///////////////////////////////////////////////////////////////////
   cutter->RegisterVariableGetter("MissHigherKL", [&]()
                                  { return (pow(EmissKL, 2) - pow(PmissKL, 2)); });
   cutter->RegisterVariableGetter("MissLowerKL", [&]()
                                  { return (pow(EmissKL, 2) - pow(PmissKL, 2)); });

   invMassKchConfig.name = "invMassKch";
   invMassKchConfig.xtitle = "m^{inv}_{K#rightarrow#pi^{+}#pi^{-}} [MeV/c^{2}]";
   invMassKchConfig.ytitle = "Counts/2";
   invMassKchConfig.bins = 100;
   invMassKchConfig.xmin = 490;
   invMassKchConfig.xmax = 500;
   invMassKchConfig.logy = true;
   invMassKchConfig.showStats = false;

   invMassKchMCConfig.name = "invMassKchMC";
   invMassKchMCConfig.xtitle = "m^{inv,MC}_{K#rightarrow#pi^{+}#pi^{-}} [MeV/c^{2}]";
   invMassKchMCConfig.ytitle = "Counts/2";
   invMassKchMCConfig.bins = 100;
   invMassKchMCConfig.xmin = -5;
   invMassKchMCConfig.xmax = 5;
   invMassKchMCConfig.logy = true;
   invMassKchMCConfig.showStats = false;

   invMassKneConfig.name = "invMassKne";
   invMassKneConfig.xtitle = "m^{inv}_{K#rightarrow#pi^{0}#pi^{0}} [MeV/c^{2}]";
   invMassKneConfig.ytitle = "Counts/10";
   invMassKneConfig.bins = 100;
   invMassKneConfig.xmin = 490;
   invMassKneConfig.xmax = 500;
   invMassKneConfig.logy = true;
   invMassKneConfig.showStats = false;

   invMassKneMCConfig.name = "invMassKneMC";
   invMassKneMCConfig.xtitle = "m^{inv,MC}_{K#rightarrow#pi^{0}#pi^{0}} [MeV/c^{2}]";
   invMassKneMCConfig.ytitle = "Counts/10";
   invMassKneMCConfig.bins = 100;
   invMassKneMCConfig.xmin = -5;
   invMassKneMCConfig.xmax = 5;
   invMassKneMCConfig.logy = true;
   invMassKneMCConfig.showStats = false;

   timeDiffConfig.name = "timeDiff";
   timeDiffConfig.xtitle = "#DeltaT [#tau_{S}]";
   timeDiffConfig.ytitle = "Counts/2";
   timeDiffConfig.bins = 20;
   timeDiffConfig.xmin = -20;
   timeDiffConfig.xmax = 20;
   timeDiffConfig.logy = true;
   timeDiffConfig.showStats = false;

   timeDiffMCConfig.name = "timeDiffMC";
   timeDiffMCConfig.xtitle = "#DeltaT [#tau_{S}]";
   timeDiffMCConfig.ytitle = "Counts/2";
   timeDiffMCConfig.bins = 301;
   timeDiffMCConfig.xmin = -10;
   timeDiffMCConfig.xmax = 10;
   timeDiffMCConfig.logy = true;
   timeDiffMCConfig.showStats = false;

   histMgr->CreateHistSet1D("invMassKch", invMassKchConfig);
   histMgr->CreateHistSet1D("invMassKne", invMassKneConfig);
   histMgr->CreateHistSet1D("timeDiff", timeDiffConfig);
   histMgr->CreateHistSet1D("invMassKchMC", invMassKchMCConfig);
   histMgr->CreateHistSet1D("invMassKneMC", invMassKneMCConfig);
   histMgr->CreateHistSet1D("timeDiffMC", timeDiffMCConfig);
}

void four_pi::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();
}

Bool_t four_pi::Process(Long64_t entry)
{
   // The Process() function is called for each entry in the tree (or possibly
   // keyed object in the case of PROOF) to be processed. The entry argument
   // specifies which entry in the currently loaded tree is to be processed.
   // When processing keyed objects with PROOF, the object is already loaded
   // and is available via the fObject pointer.
   //
   // This function should contain the \"body\" of the analysis. It can contain
   // simple or elaborate selection criteria, run algorithms on the data
   // of the event and typically fill histograms.
   //
   // The processing can be stopped by calling Abort().
   //
   // Use fStatus to set the return value of TTree::Process().
   //
   // The return value is currently not used.

   fReader.SetLocalEntry(entry);

   Float_t
       PhiMom[3] = {*Bpx, *Bpy, *Bpz},
       MissMomKS[3] = {},
       MissMomKL[3] = {};

   for (Int_t comp = 0; comp < 3; comp++)
   {
      MissMomKS[comp] = PhiMom[comp] - KchboostKS[comp] - KchrecKL[comp];
      MissMomKL[comp] = PhiMom[comp] - KchboostKL[comp] - KchrecKS[comp];
   }

   PmissKS = sqrt(pow(MissMomKS[0], 2) + pow(MissMomKS[1], 2) + pow(MissMomKS[2], 2));
   PmissKL = sqrt(pow(MissMomKL[0], 2) + pow(MissMomKL[1], 2) + pow(MissMomKL[2], 2));

   EmissKS = KchboostKS[3] - KchrecKS[3];
   EmissKL = KchboostKL[3] - KchrecKL[3];

   Float_t
       boostPhi[3] = {
           -*Bpx / *Broots,
           -*Bpy / *Broots,
           -*Bpz / *Broots},
       phiMom[4] = {*Bpx, *Bpy, *Bpz, *Broots}, trkKS_PhiCM[2][4] = {}, KchrecKS_PhiCM[4] = {}, trkKL_PhiCM[2][4], KchrecKL_PhiCM[4] = {};

   Obj->lorentz_transf(boostPhi, &trk1KS[0], trkKS_PhiCM[0]);
   Obj->lorentz_transf(boostPhi, &trk2KS[0], trkKS_PhiCM[1]);
   Obj->lorentz_transf(boostPhi, &trk1KL[0], trkKL_PhiCM[0]);
   Obj->lorentz_transf(boostPhi, &trk2KL[0], trkKL_PhiCM[1]);

   for (Int_t part = 0; part < 2; part++)
      for (Int_t comp = 0; comp < 4; comp++)
      {
         KchrecKS_PhiCM[comp] += trkKS_PhiCM[part][comp];
         KchrecKL_PhiCM[comp] += trkKL_PhiCM[part][comp];
      }

   KchrecKSMom = sqrt(pow(KchrecKS_PhiCM[0], 2) + pow(KchrecKS_PhiCM[1], 2) + pow(KchrecKS_PhiCM[2], 2));
   KchrecKLMom = sqrt(pow(KchrecKL_PhiCM[0], 2) + pow(KchrecKL_PhiCM[1], 2) + pow(KchrecKL_PhiCM[2], 2));

   Float_t weight = 1.0;

   if (1)
   {
      if (*mctruth == 1)
         weight = interf_function(*KaonChTimeCMMC - *KaonNeTimeCMMC);

      histMgr->Fill1D("invMassKch", *mctruth, KchrecKS[5], weight);
      histMgr->Fill1D("invMassKne", *mctruth, KchrecKL[5], weight);
      histMgr->Fill1D("timeDiff", *mctruth, *KaonChTimeCM - *KaonNeTimeCM, weight);

      histMgr->Fill1D("invMassKchMC", *mctruth, Kchmc[5] - KchrecKS[5], weight);
      histMgr->Fill1D("invMassKneMC", *mctruth, Knemc[5] - KchrecKL[5], weight);
      histMgr->Fill1D("timeDiffMC", *mctruth, (*KaonChTimeCMMC - *KaonNeTimeCMMC) - (*KaonChTimeCM - *KaonNeTimeCM), weight);
   }

   cutter->UpdateStats(*mctruth);

   return kTRUE;
}

void four_pi::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.
}

void four_pi::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.

   // 1D histogramy z danymi
   histMgr->DrawSet1D("invMassKch", "HIST", false);
   histMgr->DrawSet1D("invMassKne", "HIST", false);
   histMgr->DrawSet1D("timeDiff", "HIST", false);

   histMgr->DrawSet1D("invMassKchMC", "HIST", false);
   histMgr->DrawSet1D("invMassKneMC", "HIST", false);
   histMgr->DrawSet1D("timeDiffMC", "HIST", false);

   // histMgr->SaveToRoot("analysis_results.root");

   histMgr->SaveSet("invMassKch", "invMassKch");
   histMgr->SaveSet("invMassKne", "invMassKne");
   histMgr->SaveSet("timeDiff", "timeDiff");

   histMgr->SaveSet("invMassKchMC", "invMassKchMC");
   histMgr->SaveSet("invMassKneMC", "invMassKneMC");
   histMgr->SaveSet("timeDiffMC", "timeDiffMC");

   // Wyniki
   for (size_t i = 0; i < cutter->GetCuts().size(); ++i)
   {
      std::cout << "Cut " << i << ": Eff=" << cutter->GetEfficiency(i) << " +- " << cutter->GetEfficiencyError(i)
                << " Purity=" << cutter->GetPurity(i) << " +- " << cutter->GetPurityError(i)
                << " S/B=" << cutter->GetSignalToBackground(i) << " +- " << cutter->GetSignalToBackgroundError(i) << "\n";
   }

   delete histMgr;
}