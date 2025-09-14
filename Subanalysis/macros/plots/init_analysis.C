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
HistManager::HistConfig pullsTriKinFitConfig;
HistManager::Hist2DConfig TimeNeutral2dConfig;

// Konfiguracje array histogramów
HistManager::ArrayConfig momentumArrayConfig;
HistManager::ArrayConfig momentumNeutralArrayConfig;
HistManager::ArrayConfig pullsArrayConfig;
HistManager::ArrayConfig neuVtxArrayConfig;

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
   constraints.lowerBounds = {0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.0};
   constraints.upperBounds = {0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 1.0};
   constraints.initialValues = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
   constraints.fitRangeMin = 480.0;
   constraints.fitRangeMax = 520.0;

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
   invMassKchConfig.xmin = 480;
   invMassKchConfig.xmax = 520;
   invMassKchConfig.logy = false;
   invMassKchConfig.showStats = false;

   invMassKchMCConfig.name = "invMassKchMC";
   invMassKchMCConfig.xtitle = "m^{inv,MC}_{K#rightarrow#pi^{+}#pi^{-}} [MeV/c^{2}]";
   invMassKchMCConfig.ytitle = "Counts/2";
   invMassKchMCConfig.bins = 100;
   invMassKchMCConfig.xmin = -100;
   invMassKchMCConfig.xmax = 100;
   invMassKchMCConfig.logy = false;
   invMassKchMCConfig.showStats = false;

   invMassKneConfig.name = "invMassKne";
   invMassKneConfig.xtitle = "m^{inv}_{K#rightarrow#pi^{0}#pi^{0}} [MeV/c^{2}]";
   invMassKneConfig.ytitle = "Counts/10";
   invMassKneConfig.bins = 100;
   invMassKneConfig.xmin = 0;
   invMassKneConfig.xmax = 1000;
   invMassKneConfig.logy = false;
   invMassKneConfig.showStats = false;

   invMassKneMCConfig.name = "invMassKneMC";
   invMassKneMCConfig.xtitle = "m^{inv,MC}_{K#rightarrow#pi^{0}#pi^{0}} [MeV/c^{2}]";
   invMassKneMCConfig.ytitle = "Counts/10";
   invMassKneMCConfig.bins = 100;
   invMassKneMCConfig.xmin = -100;
   invMassKneMCConfig.xmax = 100;
   invMassKneMCConfig.logy = false;
   invMassKneMCConfig.showStats = false;

   timeDiffConfig.name = "timeDiff";
   timeDiffConfig.xtitle = "#DeltaT [#tau_{S}]";
   timeDiffConfig.ytitle = "Counts/2";
   timeDiffConfig.bins = 200;
   timeDiffConfig.xmin = -200;
   timeDiffConfig.xmax = 200;
   timeDiffConfig.logy = false;
   timeDiffConfig.showStats = false;

   timeDiffMCConfig.name = "timeDiffMC";
   timeDiffMCConfig.xtitle = "#DeltaT [#tau_{S}]";
   timeDiffMCConfig.ytitle = "Counts/2";
   timeDiffMCConfig.bins = 20;
   timeDiffMCConfig.xmin = -20;
   timeDiffMCConfig.xmax = 20;
   timeDiffMCConfig.logy = false;
   timeDiffMCConfig.showStats = false;

   chi2TriKinFitConfig.name = "chi2TriKinFit";
   chi2TriKinFitConfig.xtitle = "#DeltaT [#tau_{S}]";
   chi2TriKinFitConfig.ytitle = "Counts/2";
   chi2TriKinFitConfig.bins = 50;
   chi2TriKinFitConfig.xmin = -0.2;
   chi2TriKinFitConfig.xmax = 1.2;
   chi2TriKinFitConfig.logy = true;
   chi2TriKinFitConfig.showStats = false;

   pullsTriKinFitConfig.name = "pullsTriKinFit";
   pullsTriKinFitConfig.xtitle = "Pull value [-]";
   pullsTriKinFitConfig.ytitle = "Counts";
   pullsTriKinFitConfig.bins = 50;
   pullsTriKinFitConfig.xmin = -10;
   pullsTriKinFitConfig.xmax = 10;
   pullsTriKinFitConfig.logy = false;
   pullsTriKinFitConfig.showStats = false;

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

   // Konfiguracja array histogramów dla składowych pędu Kchrec
   HistManager::HistConfig momentumConfig;
   momentumConfig.bins = 50;
   momentumConfig.xmin = -10;
   momentumConfig.xmax = 10;
   momentumConfig.xtitle = "p [MeV/c]";
   momentumConfig.ytitle = "Counts";
   momentumConfig.logy = false;
   momentumConfig.showStats = false;

   momentumArrayConfig.baseName = "KchrecMomentum";
   momentumArrayConfig.baseTitle = "K_{ch} Reconstructed Momentum Components";
   momentumArrayConfig.arraySize = 4; // px, py, pz, E
   momentumArrayConfig.varNames = {"px", "py", "pz", "E"};
   momentumArrayConfig.varTitles = {"p_{x}", "p_{y}", "p_{z}", "E"};
   momentumArrayConfig.commonConfig = momentumConfig;
   momentumArrayConfig.commonConfig.xtitle = "p [MeV/c]"; // Ogólny tytuł dla px,py,pz

   // Konfiguracja array histogramów dla składowych pędu Kchrec
   HistManager::HistConfig momentumNeutralConfig;
   momentumNeutralConfig.bins = 50;
   momentumNeutralConfig.xmin = -100;
   momentumNeutralConfig.xmax = 100;
   momentumNeutralConfig.xtitle = "p [MeV/c]";
   momentumNeutralConfig.ytitle = "Counts";
   momentumNeutralConfig.logy = false;
   momentumNeutralConfig.showStats = false;

   momentumNeutralArrayConfig.baseName = "KnerecMomentum";
   momentumNeutralArrayConfig.baseTitle = "K_{ne} Reconstructed Momentum Components";
   momentumNeutralArrayConfig.arraySize = 4; // px, py, pz, E
   momentumNeutralArrayConfig.varNames = {"px", "py", "pz", "E"};
   momentumNeutralArrayConfig.varTitles = {"p_{x}", "p_{y}", "p_{z}", "E"};
   momentumNeutralArrayConfig.commonConfig = momentumNeutralConfig;
   momentumNeutralArrayConfig.commonConfig.xtitle = "p [MeV/c]"; // Ogólny tytuł dla px,py,pz

   // Konfiguracja array histogramów dla składowych pędu Kchrec
   HistManager::HistConfig neuVtxConfig;
   neuVtxConfig.bins = 50;
   neuVtxConfig.xmin = -1;
   neuVtxConfig.xmax = 1;
   neuVtxConfig.xtitle = "x [MeV/c]";
   neuVtxConfig.ytitle = "Counts";
   neuVtxConfig.logy = false;
   neuVtxConfig.showStats = false;

   neuVtxArrayConfig.baseName = "NeuVtx";
   neuVtxArrayConfig.baseTitle = "Neutral vertex resolution";
   neuVtxArrayConfig.arraySize = 4; // px, py, pz, E
   neuVtxArrayConfig.varNames = {"x", "y", "z", "t"};
   neuVtxArrayConfig.varTitles = {"x", "y", "z", "t"};
   neuVtxArrayConfig.commonConfig = neuVtxConfig;
   neuVtxArrayConfig.commonConfig.xtitle = "x [cm]"; // Ogólny tytuł dla px,py,pz
   
   // Specjalna konfiguracja dla energii
   // momentumArrayConfig.commonConfig.xmax = 1000; // E ma większy zakres

   // Konfiguracja array histogramów dla pullsTriKinFit
   HistManager::HistConfig pullArrayConfig;
   pullArrayConfig.bins = 50;
   pullArrayConfig.xmin = -10;
   pullArrayConfig.xmax = 10;
   pullArrayConfig.xtitle = "Pull value [-]";
   pullArrayConfig.ytitle = "Counts";
   pullArrayConfig.logy = false;
   pullArrayConfig.showStats = false;

   pullsArrayConfig.baseName = "PullsTriKinFit";
   pullsArrayConfig.baseTitle = "Trilateration Kinematic Fit Pulls";
   pullsArrayConfig.arraySize = 24; // Assuming 5 pull values
   pullsArrayConfig.varNames = {"pull0", "pull1", "pull2", "pull3", "pull4"};
   pullsArrayConfig.varTitles = {"Pull_{0}", "Pull_{1}", "Pull_{2}", "Pull_{3}", "Pull_{4}"};
   pullsArrayConfig.commonConfig = pullArrayConfig;

   histMgr->CreateHistSet1D("invMassKch", invMassKchConfig);
   histMgr->CreateHistSet1D("invMassKne", invMassKneConfig);
   histMgr->CreateHistSet1D("timeDiff", timeDiffConfig);
   histMgr->CreateHistSet1D("invMassKchMC", invMassKchMCConfig);
   histMgr->CreateHistSet1D("invMassKneMC", invMassKneMCConfig);
   histMgr->CreateHistSet1D("timeDiffMC", timeDiffMCConfig);
   histMgr->CreateHistSet1D("chi2TriKinFit", chi2TriKinFitConfig);
   histMgr->CreateHistSet1D("pullsTriKinFit", pullsTriKinFitConfig);
   histMgr->CreateHistSet2D("TimeNeutral2d", TimeNeutral2dConfig);

   // Stwórz array histogramy
   histMgr->CreateHistArray1D(momentumArrayConfig);
   histMgr->CreateHistArray1D(momentumNeutralArrayConfig);
   histMgr->CreateHistArray1D(neuVtxArrayConfig);
   histMgr->CreateHistArray1D(pullsArrayConfig);

   histMgr->SetFitConstraints("invMassKch", constraints);
   histMgr->SetFitConstraints("invMassKne", constraints);
   histMgr->SetFitConstraints("timeDiff", constraints);
   histMgr->SetFitConstraints("invMassKchMC", constraints);
   histMgr->SetFitConstraints("invMassKneMC", constraints);
   histMgr->SetFitConstraints("timeDiffMC", constraints);
   histMgr->SetFitConstraints("chi2TriKinFit", constraints);
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

   if (1)
   {
      if (*mcflag == 1)
      {
         // if (*mctruth == 1)
         // weight = interf_function(*KaonChTimeCMMC - *KaonNeTimeCMMC);

         Float_t minv4gam_tri = sqrt(pow(KnetriKinFit[3],2) - pow(KnetriKinFit[0],2) - pow(KnetriKinFit[1],2) - pow(KnetriKinFit[2],2));

         histMgr->Fill1D("invMassKch", *mctruth, Kchrec[5], weight);
         histMgr->Fill1D("invMassKne", *mctruth, minv4gam_tri, weight);
         histMgr->Fill1D("timeDiff", *mctruth, *KaonChTimeCM - *KaonNeTimeCM, weight);
         histMgr->Fill1D("chi2TriKinFit", *mctruth, TMath::Prob(*Chi2TriKinFit, 5), weight);
         histMgr->Fill1D("pullsTriKinFit", *mctruth, pullsTriKinFit[4], weight);

         // Wypełnij array histogramy dla składowych pędu Kchrec
         std::vector<Double_t> momentum = {Kchboost[0] - Kchmc[0], Kchboost[1] - Kchmc[1], Kchboost[2] - Kchmc[2], Kchboost[3] - Kchmc[3]};
         histMgr->FillArrayAll1D("KchrecMomentum", *mctruth, momentum, weight);

         // Wypełnij array histogramy dla składowych pędu Kchrec
         // momentum = {KneTriangle[0] - Knemc[0], KneTriangle[1] - Knemc[1], KneTriangle[2] - Knemc[2], KneTriangle[3] - Knemc[3]};
         momentum = {KnetriKinFit[0] - Knemc[0], KnetriKinFit[1] - Knemc[1], KnetriKinFit[2] - Knemc[2], KnetriKinFit[3] - Knemc[3]};
         histMgr->FillArrayAll1D("KnerecMomentum", *mctruth, momentum, weight);

         // Wypełnij array histogramy dla składowych neu vtx
         // momentum = {KneTriangle[0] - Knemc[0], KneTriangle[1] - Knemc[1], KneTriangle[2] - Knemc[2], KneTriangle[3] - Knemc[3]};
         momentum = {KnetriKinFit[6] - Knemc[6], KnetriKinFit[7] - Knemc[7], KnetriKinFit[8] - Knemc[8], KnetriKinFit[9] - *KaonNeTimeLAB * tau_S};
         histMgr->FillArrayAll1D("NeuVtx", *mctruth, momentum, weight);

         // Wypełnij array histogramy dla pulls (wszystkie 5 składowych)
         std::vector<Double_t> pulls = {pullsTriKinFit[0], pullsTriKinFit[1], pullsTriKinFit[2], 
                                       pullsTriKinFit[3], pullsTriKinFit[4], pullsTriKinFit[5], pullsTriKinFit[6], pullsTriKinFit[7], 
                                       pullsTriKinFit[8], pullsTriKinFit[9], pullsTriKinFit[10], pullsTriKinFit[11], pullsTriKinFit[12], 
                                       pullsTriKinFit[13], pullsTriKinFit[14], pullsTriKinFit[15], pullsTriKinFit[16], pullsTriKinFit[17], 
                                       pullsTriKinFit[18], pullsTriKinFit[19], pullsTriKinFit[20], pullsTriKinFit[21], pullsTriKinFit[22],
                                    pullsTriKinFit[23]};
         histMgr->FillArrayAll1D("PullsTriKinFit", *mctruth, pulls, weight);

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

         // Wypełnij array histogramy dla danych - składowe pędu Kchrec
         std::vector<Double_t> momentumData = {Kchboost[0] - Kchmc[0], Kchboost[1] - Kchmc[1], Kchboost[2] - Kchmc[2], Kchboost[3] - Kchmc[3]};
         histMgr->FillArrayAllData1D("KchrecMomentum", momentumData, weight);

         // Wypełnij array histogramy dla danych - pulls (jeśli dostępne dla danych)
         std::vector<Double_t> pullsData = {pullsTriKinFit[0], pullsTriKinFit[1], pullsTriKinFit[2], 
                                           pullsTriKinFit[3], pullsTriKinFit[4]};
         histMgr->FillArrayAllData1D("PullsTriKinFit", pullsData, weight);
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

   histMgr->SetNormalizationType(HistManager::NormalizationType::FRACTION_FIT); // Ustaw proste skalowanie --- IGNORE ---

   // 1D histogramy z danymi
   histMgr->DrawSet1D("invMassKch", "HIST", true);
   histMgr->DrawSet1D("invMassKne", "HIST", true);
   histMgr->DrawSet1D("timeDiff", "HIST", true);
   histMgr->DrawSet1D("chi2TriKinFit", "HIST", true);
   histMgr->DrawSet1D("pullsTriKinFit", "HIST", true);

   histMgr->DrawSet2D("TimeNeutral2d", "COLZ", true);

   histMgr->DrawSet1D("invMassKchMC", "HIST", true);
   histMgr->DrawSet1D("invMassKneMC", "HIST", true);
   histMgr->DrawSet1D("timeDiffMC", "HIST", true);

   // Rysuj array histogramy
   histMgr->DrawArray1D("KchrecMomentum", true); // z danymi
   histMgr->DrawArray1D("KnerecMomentum", true); // z danymi
   histMgr->DrawArray1D("NeuVtx", true); // z danymi
   histMgr->DrawArray1D("PullsTriKinFit", true); // z danymi

   // histMgr->SaveToRoot("analysis_results.root");

   histMgr->SaveSet("invMassKch", "invMassKch");
   histMgr->SaveSet("invMassKne", "invMassKne");
   histMgr->SaveSet("timeDiff", "timeDiff");
   histMgr->SaveSet("chi2TriKinFit", "chi2TriKinFit");
   histMgr->SaveSet("pullsTriKinFit", "pullsTriKinFit");

   histMgr->SaveSet("TimeNeutral2d", "TimeNeutral2d");

   histMgr->SaveSet("invMassKchMC", "invMassKchMC");
   histMgr->SaveSet("invMassKneMC", "invMassKneMC");
   histMgr->SaveSet("timeDiffMC", "timeDiffMC");

   // Zapisz array histogramy
   histMgr->ExportSet("KchrecMomentum", "KchrecMomentum", HistManager::ImageFormat::SVG);
   histMgr->ExportSet("KnerecMomentum", "KnerecMomentum", HistManager::ImageFormat::SVG);
   histMgr->ExportSet("NeuVtx", "NeuVtx", HistManager::ImageFormat::SVG);
   histMgr->ExportSet("PullsTriKinFit", "PullsTriKinFit", HistManager::ImageFormat::SVG);

   // Wyniki
   for (size_t i = 0; i < cutter->GetCuts().size(); ++i)
   {
      std::cout << "Cut " << i << ": Eff=" << cutter->GetEfficiency(i) << " +- " << cutter->GetEfficiencyError(i)
                << " Purity=" << cutter->GetPurity(i) << " +- " << cutter->GetPurityError(i)
                << " S/B=" << cutter->GetSignalToBackground(i) << " +- " << cutter->GetSignalToBackgroundError(i) << "\n";
   }

   delete histMgr;
}