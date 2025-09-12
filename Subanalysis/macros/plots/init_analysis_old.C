#define init_analysis_cxx
// The class definition in init_analysis.h has been generated automatically
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
// root> T->Process("init_analysis.C")
// root> T->Process("init_analysis.C","some options")
// root> T->Process("init_analysis.C+")
//

#include "init_analysis.h"
#include <TH2.h>
#include <TStyle.h>
#include <HistManager.h>
#include <interf_function.h>
#include <StatisticalCutter.h>

HistManager *histMgr;
StatisticalCutter *cutter;
HistManager::HistConfig invMassKneConfig;
HistManager::HistConfig invMassKneMCConfig;
HistManager::HistConfig invMassKchConfig;
HistManager::HistConfig invMassKchMCConfig;
HistManager::HistConfig timeDiffConfig;
HistManager::HistConfig timeDiffMCConfig;
HistManager::HistConfig chi2TriKinFitConfig;
HistManager::Hist2DConfig TimeNeutral2dConfig;

void init_analysis::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

   histMgr = new HistManager(channNum, channColor, channNames, kFullCircle, kBlack, kOrange);

   std::string cutFileName = "/data/ssd/gamrat/KLOE/Subanalysis/Properties/cut-limits-final.json";

   cutter = new StatisticalCutter(cutFileName, 1, KLOE::HypothesisCode::SIGNAL);
   HistManager::FitConstraints constraints(7);
   constraints.lowerBounds = {0.2, 0.0, 0.2, 0.2, 0.2, 0.0, 0.2};
   constraints.upperBounds = {1.0, 0.1, 1.0, 1.0, 1.0, 0.1, 1.0};

   cutter->RegisterVariableGetter("InvMassKch", [&]()
                                  { return Kchrec[5]; });
   cutter->RegisterCentralValueGetter("InvMassKch", [&]()
                                      { return mK0; });

   cutter->RegisterVariableGetter("Qmiss", [&]()
                                  { return *Qmiss; });

   cutter->RegisterVariableGetter("InvMassKne", [&]()
                                  { return *minv4gam; });
   cutter->RegisterCentralValueGetter("InvMassKne", [&]()
                                      { return mK0; });

   invMassKchConfig.name = "invMassKch";
   invMassKchConfig.xtitle = "m^{inv}_{K#rightarrow#pi^{+}#pi^{-}} [MeV/c^{2}]";
   invMassKchConfig.ytitle = "Counts/2";
   invMassKchConfig.bins = 100;
   invMassKchConfig.xmin = 400;
   invMassKchConfig.xmax = 600;
   invMassKchConfig.logy = true;
   invMassKchConfig.showStats = false;

   invMassKchMCConfig.name = "invMassKchMC";
   invMassKchMCConfig.xtitle = "m^{inv,MC}_{K#rightarrow#pi^{+}#pi^{-}} [MeV/c^{2}]";
   invMassKchMCConfig.ytitle = "Counts/2";
   invMassKchMCConfig.bins = 100;
   invMassKchMCConfig.xmin = -100;
   invMassKchMCConfig.xmax = 100;
   invMassKchMCConfig.logy = true;
   invMassKchMCConfig.showStats = false;

   invMassKneConfig.name = "invMassKne";
   invMassKneConfig.xtitle = "m^{inv}_{K#rightarrow#pi^{0}#pi^{0}} [MeV/c^{2}]";
   invMassKneConfig.ytitle = "Counts/10";
   invMassKneConfig.bins = 100;
   invMassKneConfig.xmin = 0;
   invMassKneConfig.xmax = 1000;
   invMassKneConfig.logy = true;
   invMassKneConfig.showStats = false;

   invMassKneMCConfig.name = "invMassKneMC";
   invMassKneMCConfig.xtitle = "m^{inv,MC}_{K#rightarrow#pi^{0}#pi^{0}} [MeV/c^{2}]";
   invMassKneMCConfig.ytitle = "Counts/10";
   invMassKneMCConfig.bins = 100;
   invMassKneMCConfig.xmin = -100;
   invMassKneMCConfig.xmax = 100;
   invMassKneMCConfig.logy = true;
   invMassKneMCConfig.showStats = false;

   timeDiffConfig.name = "timeDiff";
   timeDiffConfig.xtitle = "#DeltaT [#tau_{S}]";
   timeDiffConfig.ytitle = "Counts/2";
   timeDiffConfig.bins = 200;
   timeDiffConfig.xmin = -200;
   timeDiffConfig.xmax = 200;
   timeDiffConfig.logy = true;
   timeDiffConfig.showStats = false;

   timeDiffMCConfig.name = "timeDiffMC";
   timeDiffMCConfig.xtitle = "#DeltaT [#tau_{S}]";
   timeDiffMCConfig.ytitle = "Counts/2";
   timeDiffMCConfig.bins = 20;
   timeDiffMCConfig.xmin = -20;
   timeDiffMCConfig.xmax = 20;
   timeDiffMCConfig.logy = true;
   timeDiffMCConfig.showStats = false;

   chi2TriKinFitConfig.name = "chi2TriKinFit";
   chi2TriKinFitConfig.xtitle = "#DeltaT [#tau_{S}]";
   chi2TriKinFitConfig.ytitle = "Counts/2";
   chi2TriKinFitConfig.bins = 50;
   chi2TriKinFitConfig.xmin = -10;
   chi2TriKinFitConfig.xmax = 100;
   chi2TriKinFitConfig.logy = true;
   chi2TriKinFitConfig.showStats = false;

   TimeNeutral2dConfig.name = "TimeNeutral2d";
   TimeNeutral2dConfig.xtitle = "t^{MC}_{neu} [#tau_{S}]";
   TimeNeutral2dConfig.ytitle = "t^{rec}_{neu} [#tau_{S}]";
   TimeNeutral2dConfig.bins = 50;
   TimeNeutral2dConfig.xmin = -10;
   TimeNeutral2dConfig.xmax = 100;
   TimeNeutral2dConfig.binsy = 50;
   TimeNeutral2dConfig.ymin = -10;
   TimeNeutral2dConfig.ymax = 100;
   TimeNeutral2dConfig.logy = false;
   TimeNeutral2dConfig.showStats = false;

   histMgr->CreateHistSet1D("invMassKch", invMassKchConfig);
   histMgr->CreateHistSet1D("invMassKne", invMassKneConfig);
   histMgr->CreateHistSet1D("timeDiff", timeDiffConfig);
   histMgr->CreateHistSet1D("invMassKchMC", invMassKchMCConfig);
   histMgr->CreateHistSet1D("invMassKneMC", invMassKneMCConfig);
   histMgr->CreateHistSet1D("timeDiffMC", timeDiffMCConfig);
   histMgr->CreateHistSet1D("chi2TriKinFit", chi2TriKinFitConfig);
   histMgr->CreateHistSet2D("TimeNeutral2d", TimeNeutral2dConfig);
}

void init_analysis::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();
}

Bool_t init_analysis::Process(Long64_t entry)
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

   Float_t weight = 1.0;

   std::cout << "Processing entry: " << entry << "\r" << std::flush;

   if (1)
   {
      if (*mcflag == 1)
      {
         if (*mctruth == 1)
            weight = interf_function(*KaonChTimeCMMC - *KaonNeTimeCMMC);

         histMgr->Fill1D("invMassKch", *mctruth, Kchrec[5], weight);
         histMgr->Fill1D("invMassKne", *mctruth, *minv4gam, weight);
         histMgr->Fill1D("timeDiff", *mctruth, *KaonChTimeCM - *KaonNeTimeCM, weight);
         histMgr->Fill1D("chi2TriKinFit", *mctruth, *Chi2TriKinFit, weight);

         histMgr->Fill1D("invMassKchMC", *mctruth, Kchmc[5] - Kchrec[5], weight);
         histMgr->Fill1D("invMassKneMC", *mctruth, Knemc[5] - *minv4gam, weight);
         histMgr->Fill1D("timeDiffMC", *mctruth, (*KaonChTimeCMMC - *KaonNeTimeCMMC) - (*KaonChTimeCM - *KaonNeTimeCM), weight);
      }

      if (*mcflag == 0)
      {
         histMgr->FillData1D("invMassKch", Kchrec[5], weight);
         histMgr->FillData1D("invMassKne", *minv4gam, weight);
         histMgr->FillData1D("timeDiff", *KaonChTimeCM - *KaonNeTimeCM, weight);
         histMgr->FillData1D("chi2TriKinFit", *Chi2TriKinFit, weight);

         histMgr->FillData1D("invMassKchMC", Kchmc[5] - Kchrec[5], weight);
         histMgr->FillData1D("invMassKneMC", Knemc[5] - *minv4gam, weight);
         histMgr->FillData1D("timeDiffMC", (*KaonChTimeCMMC - *KaonNeTimeCMMC) - (*KaonChTimeCM - *KaonNeTimeCM), weight);
      }
   }

   if (*mcflag == 1)
      cutter->UpdateStats(*mctruth);

   return kTRUE;
}

void init_analysis::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.
}

void init_analysis::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.

   histMgr->SetNormalizationType(HistManager::NormalizationType::SIMPLE_SCALE); // Ustaw proste skalowanie --- IGNORE ---

   // 1D histogramy z danymi
   histMgr->DrawSet1D("invMassKch", "HIST", true);
   histMgr->DrawSet1D("invMassKne", "HIST", true);
   histMgr->DrawSet1D("timeDiff", "HIST", true);
   histMgr->DrawSet1D("chi2TriKinFit", "HIST", true);

   histMgr->DrawSet2D("TimeNeutral2d", "COLZ", true);

   histMgr->DrawSet1D("invMassKchMC", "HIST", true);
   histMgr->DrawSet1D("invMassKneMC", "HIST", true);
   histMgr->DrawSet1D("timeDiffMC", "HIST", true);

   // histMgr->SaveToRoot("analysis_results.root");

   histMgr->SaveSet("invMassKch", "invMassKch");
   histMgr->SaveSet("invMassKne", "invMassKne");
   histMgr->SaveSet("timeDiff", "timeDiff");
   histMgr->SaveSet("chi2TriKinFit", "chi2TriKinFit");

   histMgr->SaveSet("TimeNeutral2d", "TimeNeutral2d");

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