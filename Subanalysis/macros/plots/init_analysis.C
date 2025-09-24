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
HistManager::HistConfig trcFinalConfig;
HistManager::HistConfig curvMCConfig;
HistManager::HistConfig phivMCConfig;
HistManager::HistConfig cotvMCConfig;
HistManager::HistConfig xchMCConfig;
HistManager::HistConfig ychMCConfig;
HistManager::HistConfig zchMCConfig;
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

   histMgr = new HistManager(channNum, channColor, channNames, kFullCircle, kBlack, 0.0, kOrange);

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
   invMassKchConfig.logy = true;
   invMassKchConfig.showStats = false;

   invMassKchMCConfig.name = "invMassKchMC";
   invMassKchMCConfig.xtitle = "m^{inv,MC}_{K#rightarrow#pi^{+}#pi^{-}} - m_{K_{0}} [MeV/c^{2}]";
   invMassKchMCConfig.ytitle = "Counts";
   invMassKchMCConfig.bins = 100;
   invMassKchMCConfig.xmin = -100;
   invMassKchMCConfig.xmax = 100;
   invMassKchMCConfig.logy = true;
   invMassKchMCConfig.showStats = false;

   invMassKneConfig.name = "invMassKne";
   invMassKneConfig.xtitle = "m^{inv}_{K#rightarrow#pi^{0}#pi^{0}} [MeV/c^{2}]";
   invMassKneConfig.ytitle = "Counts";
   invMassKneConfig.bins = 100;
   invMassKneConfig.xmin = 0;
   invMassKneConfig.xmax = 1000;
   invMassKneConfig.logy = true;
   invMassKneConfig.showStats = false;

   invMassKneMCConfig.name = "invMassKneMC";
   invMassKneMCConfig.xtitle = "m^{inv,MC}_{K#rightarrow#pi^{0}#pi^{0}} - m_{K_{0}} [MeV/c^{2}]";
   invMassKneMCConfig.ytitle = "Counts";
   invMassKneMCConfig.bins = 100;
   invMassKneMCConfig.xmin = -300;
   invMassKneMCConfig.xmax = 300;
   invMassKneMCConfig.logy = true;
   invMassKneMCConfig.showStats = false;

   timeDiffConfig.name = "timeDiff";
   timeDiffConfig.xtitle = "#DeltaT^{MC,rec} [#tau_{S}]";
   timeDiffConfig.ytitle = "Counts";
   timeDiffConfig.bins = 300;
   timeDiffConfig.xmin = -300;
   timeDiffConfig.xmax = 300;
   timeDiffConfig.logy = true;
   timeDiffConfig.showStats = false;

   timeDiffMCConfig.name = "timeDiffMC";
   timeDiffMCConfig.xtitle = "#DeltaT^{MC,rec} - #DeltaT^{MC,gen} [#tau_{S}]";
   timeDiffMCConfig.ytitle = "Counts";
   timeDiffMCConfig.bins = 200;
   timeDiffMCConfig.xmin = -100;
   timeDiffMCConfig.xmax = 100;
   timeDiffMCConfig.logy = true;
   timeDiffMCConfig.showStats = false;

   chi2TriKinFitConfig.name = "chi2TriKinFit";
   chi2TriKinFitConfig.xtitle = "#chi^{2}_{tri} [-]";
   chi2TriKinFitConfig.ytitle = "Counts";
   chi2TriKinFitConfig.bins = 50;
   chi2TriKinFitConfig.xmin = -10;
   chi2TriKinFitConfig.xmax = 50;
   chi2TriKinFitConfig.logy = true;
   chi2TriKinFitConfig.showStats = false;

   trcFinalConfig.name = "trcFinal";
   trcFinalConfig.xtitle = "Trc sum [ns]";
   trcFinalConfig.ytitle = "Counts";
   trcFinalConfig.bins = 50;
   trcFinalConfig.xmin = -10;
   trcFinalConfig.xmax = 5;
   trcFinalConfig.logy = true;
   trcFinalConfig.showStats = false;

   curvMCConfig.name = "curvMC";
   curvMCConfig.xtitle = "Curv^{MC,rec} - Curv^{MC,gen} [cm^{-1}]";
   curvMCConfig.ytitle = "Counts";
   curvMCConfig.bins = 50;
   curvMCConfig.xmin = -50;
   curvMCConfig.xmax = 10;
   curvMCConfig.logy = true;
   curvMCConfig.showStats = false;

   phivMCConfig.name = "phivMC";
   phivMCConfig.xtitle = "Phiv^{MC,rec} - Phiv^{MC,gen} [rad]";
   phivMCConfig.ytitle = "Counts";
   phivMCConfig.bins = 40;
   phivMCConfig.xmin = -2 * M_PI;
   phivMCConfig.xmax = 2 * M_PI;
   phivMCConfig.logy = true;
   phivMCConfig.showStats = false;

   cotvMCConfig.name = "cotvMC";
   cotvMCConfig.xtitle = "Cotv^{MC,rec} - Cotv^{MC,gen} [-]";
   cotvMCConfig.ytitle = "Counts";
   cotvMCConfig.bins = 50;
   cotvMCConfig.xmin = -10;
   cotvMCConfig.xmax = 10;
   cotvMCConfig.logy = true;
   cotvMCConfig.showStats = false;

   xchMCConfig.name = "xchMC";
   xchMCConfig.xtitle = "x^{MC,rec}_{ch} - x^{MC,gen}_{ch} [cm]";
   xchMCConfig.ytitle = "Counts";
   xchMCConfig.bins = 50;
   xchMCConfig.xmin = -10;
   xchMCConfig.xmax = 10;
   xchMCConfig.logy = true;
   xchMCConfig.showStats = false;

   ychMCConfig.name = "ychMC";
   ychMCConfig.xtitle = "y^{MC,rec}_{ch} - y^{MC,gen}_{ch} [cm]";
   ychMCConfig.ytitle = "Counts";
   ychMCConfig.bins = 50;
   ychMCConfig.xmin = -10;
   ychMCConfig.xmax = 10;
   ychMCConfig.logy = true;
   ychMCConfig.showStats = false;

   zchMCConfig.name = "zchMC";
   zchMCConfig.xtitle = "z^{MC,rec}_{ch} - z^{MC,gen}_{ch} [cm]";
   zchMCConfig.ytitle = "Counts";
   zchMCConfig.bins = 50;
   zchMCConfig.xmin = -4;
   zchMCConfig.xmax = 4;
   zchMCConfig.logy = true;
   zchMCConfig.showStats = false;

   TimeNeutral2dConfig.name = "TimeNeutral2d";
   TimeNeutral2dConfig.xtitle = "t^{MC}_{neu} [#tau_{S}]";
   TimeNeutral2dConfig.ytitle = "t^{rec}_{neu} [#tau_{S}]";
   TimeNeutral2dConfig.bins = 50;
   TimeNeutral2dConfig.xmin = -10;
   TimeNeutral2dConfig.xmax = 400;
   TimeNeutral2dConfig.binsy = 50;
   TimeNeutral2dConfig.ymin = -10;
   TimeNeutral2dConfig.ymax = 400;
   TimeNeutral2dConfig.logy = true;
   TimeNeutral2dConfig.showStats = false;

   // Konfiguracja array histogramów dla składowych pędu Kchrec
   HistManager::HistConfig momentumConfig;
   momentumConfig.bins = 50;
   momentumConfig.xmin = -10;
   momentumConfig.xmax = 10;
   momentumConfig.xtitle = "p [MeV/c]";
   momentumConfig.ytitle = "Counts";
   momentumConfig.logy = true;
   momentumConfig.showStats = false;

   momentumArrayConfig.baseName = "KchrecMomentum";
   momentumArrayConfig.baseTitle = "K_{ch} Reconstructed Momentum Components";
   momentumArrayConfig.arraySize = 4; // px, py, pz, E
   momentumArrayConfig.varNames = {"px", "py", "pz", "E"};
   momentumArrayConfig.varTitles = {"p_{x}", "p_{y}", "p_{z}", "E"};
   momentumArrayConfig.commonConfig = momentumConfig;
   momentumArrayConfig.commonConfig.xtitle = "p [MeV/c]"; // Ogólny tytuł dla px,py,pz

   // NOWE: Indywidualne konfiguracje dla każdej składowej pędu
   // px, py: standardowy zakres pędu poprzecznego
   momentumArrayConfig.SetIndividualBinning(0, 80, -300, 300); // px
   momentumArrayConfig.SetIndividualAxisTitles(0, "p_{x} [MeV/c]", "Events / 10 MeV/c");

   momentumArrayConfig.SetIndividualBinning(1, 80, 300, 300); // py
   momentumArrayConfig.SetIndividualAxisTitles(1, "p_{y} [MeV/c]", "Events / 10 MeV/c");

   // pz: szerszy zakres dla pędu podłużnego (kaony lecą głównie do przodu)
   momentumArrayConfig.SetIndividualBinning(2, 100, -300, 300); // pz
   momentumArrayConfig.SetIndividualAxisTitles(2, "p_{z} [MeV/c]", "Events / 10 MeV/c");

   // E: energia, inny zakres i więcej binów dla lepszej precyzji
   momentumArrayConfig.SetIndividualBinning(3, 120, 400, 600); // E
   momentumArrayConfig.SetIndividualAxisTitles(3, "E [MeV]", "Events / 5 MeV");

   // Konfiguracja array histogramów dla składowych pędu Kchrec
   HistManager::HistConfig momentumNeutralConfig;
   momentumNeutralConfig.bins = 50;
   momentumNeutralConfig.xmin = -100;
   momentumNeutralConfig.xmax = 100;
   momentumNeutralConfig.xtitle = "p [MeV/c]";
   momentumNeutralConfig.ytitle = "Counts";
   momentumNeutralConfig.logy = true;
   momentumNeutralConfig.showStats = false;

   momentumNeutralArrayConfig.baseName = "KnerecMomentum";
   momentumNeutralArrayConfig.baseTitle = "K_{ne} Reconstructed Momentum Components";
   momentumNeutralArrayConfig.arraySize = 4; // px, py, pz, E
   momentumNeutralArrayConfig.varNames = {"px", "py", "pz", "E"};
   momentumNeutralArrayConfig.varTitles = {"p_{x}", "p_{y}", "p_{z}", "E"};
   momentumNeutralArrayConfig.commonConfig = momentumNeutralConfig;
   momentumNeutralArrayConfig.commonConfig.xtitle = "p [MeV/c]"; // Ogólny tytuł dla px,py,pz

   // NOWE: Indywidualne konfiguracje dla każdej składowej pędu
   // px, py: standardowy zakres pędu poprzecznego
   momentumNeutralArrayConfig.SetIndividualBinning(0, 80, -300, 300); // px
   momentumNeutralArrayConfig.SetIndividualAxisTitles(0, "p_{x} [MeV/c]", "Events / 10 MeV/c");

   momentumNeutralArrayConfig.SetIndividualBinning(1, 80, -300, 300); // py
   momentumNeutralArrayConfig.SetIndividualAxisTitles(1, "p_{y} [MeV/c]", "Events / 10 MeV/c");

   // pz: szerszy zakres dla pędu podłużnego (kaony lecą głównie do przodu)
   momentumNeutralArrayConfig.SetIndividualBinning(2, 100, -300, 300); // pz
   momentumNeutralArrayConfig.SetIndividualAxisTitles(2, "p_{z} [MeV/c]", "Events / 10 MeV/c");

   // E: energia, inny zakres i więcej binów dla lepszej precyzji
   momentumNeutralArrayConfig.SetIndividualBinning(3, 120, 400, 600); // E
   momentumNeutralArrayConfig.SetIndividualAxisTitles(3, "E [MeV]", "Events / 5 MeV");

   // Konfiguracja array histogramów dla składowych pędu Kchrec
   HistManager::HistConfig neuVtxConfig;
   neuVtxConfig.bins = 50;
   neuVtxConfig.xmin = -1;
   neuVtxConfig.xmax = 1;
   neuVtxConfig.xtitle = "x [MeV/c]";
   neuVtxConfig.ytitle = "Counts";
   neuVtxConfig.logy = true;
   neuVtxConfig.showStats = false;

   neuVtxArrayConfig.baseName = "NeuVtx";
   neuVtxArrayConfig.baseTitle = "Neutral vertex resolution";
   neuVtxArrayConfig.arraySize = 4; // px, py, pz, E
   neuVtxArrayConfig.varNames = {"x", "y", "z", "t"};
   neuVtxArrayConfig.varTitles = {"x", "y", "z", "t"};
   neuVtxArrayConfig.commonConfig = neuVtxConfig;
   neuVtxArrayConfig.commonConfig.xtitle = "x [cm]"; // Ogólny tytuł dla px,py,pz

   // NOWE: Indywidualne konfiguracje dla każdej składowej pędu
   // px, py: standardowy zakres pędu poprzecznego
   neuVtxArrayConfig.SetIndividualBinning(0, 80, -200, 200); // px
   neuVtxArrayConfig.SetIndividualAxisTitles(0, "x [cm]", "Events / 0.2 cm");

   neuVtxArrayConfig.SetIndividualBinning(1, 80, -200, 200); // py
   neuVtxArrayConfig.SetIndividualAxisTitles(1, "y [cm]", "Events / 0.2 cm");

   // pz: szerszy zakres dla pędu podłużnego (kaony lecą głównie do przodu)
   neuVtxArrayConfig.SetIndividualBinning(2, 80, -180, 180); // pz
   neuVtxArrayConfig.SetIndividualAxisTitles(2, "z [cm]", "Events / 0.4 cm");

   // E: energia, inny zakres i więcej binów dla lepszej precyzji
   neuVtxArrayConfig.SetIndividualBinning(3, 50, -5, 10); // E
   neuVtxArrayConfig.SetIndividualAxisTitles(3, "t [ns]", "Events / 0.5 ns");

   // Specjalna konfiguracja dla energii
   // momentumArrayConfig.commonConfig.xmax = 1000; // E ma większy zakres

   // Konfiguracja array histogramów dla pullsTriKinFit
   HistManager::HistConfig pullArrayConfig;
   pullArrayConfig.bins = 50;
   pullArrayConfig.xmin = -10;
   pullArrayConfig.xmax = 10;
   pullArrayConfig.xtitle = "Pull value [-]";
   pullArrayConfig.ytitle = "Counts";
   pullArrayConfig.logy = true;
   pullArrayConfig.showStats = false;

   pullsArrayConfig.baseName = "PullsTriKinFit";
   pullsArrayConfig.baseTitle = "Trilateration Kinematic Fit Pulls";
   pullsArrayConfig.arraySize = 24; // Assuming 24 pull values
   pullsArrayConfig.varNames = {"pull0", "pull1", "pull2", "pull3", "pull4"};
   pullsArrayConfig.varTitles = {"Pull_{0}", "Pull_{1}", "Pull_{2}", "Pull_{3}", "Pull_{4}"};
   pullsArrayConfig.commonConfig = pullArrayConfig;

   for (Int_t k = 0; k < 4; k++)
   {
      pullsArrayConfig.SetIndividualBinning(k * 5, 50, -10, 10);
      pullsArrayConfig.SetIndividualBinning(k * 5 + 1, 50, -10, 10);
      pullsArrayConfig.SetIndividualBinning(k * 5 + 2, 50, -10, 10);
      pullsArrayConfig.SetIndividualBinning(k * 5 + 3, 50, -10, 10);
      pullsArrayConfig.SetIndividualBinning(k * 5 + 4, 50, -10, 10);
      pullsArrayConfig.SetIndividualAxisTitles(k * 5, Form("Cluster %d Pull_{x_{cl}}", k), "Counts");
      pullsArrayConfig.SetIndividualAxisTitles(k * 5 + 1, Form("Cluster %d Pull_{y_{cl}}", k), "Counts");
      pullsArrayConfig.SetIndividualAxisTitles(k * 5 + 2, Form("Cluster %d Pull_{z_{cl}}", k), "Counts");
      pullsArrayConfig.SetIndividualAxisTitles(k * 5 + 3, Form("Cluster %d Pull_{T_{cl}}", k), "Counts");
      pullsArrayConfig.SetIndividualAxisTitles(k * 5 + 4, Form("Cluster %d Pull_{E_{cl}}", k), "Counts");
   }

   pullsArrayConfig.SetIndividualBinning(20, 50, -10, 10);
   pullsArrayConfig.SetIndividualBinning(21, 50, -10, 10);
   pullsArrayConfig.SetIndividualBinning(22, 50, -10, 10);
   pullsArrayConfig.SetIndividualBinning(23, 50, -10, 10);
   pullsArrayConfig.SetIndividualAxisTitles(20, "Pull_{p^{#phi}_{x}}", "Counts");
   pullsArrayConfig.SetIndividualAxisTitles(21, "Pull_{p^{#phi}_{y}}", "Counts");
   pullsArrayConfig.SetIndividualAxisTitles(22, "Pull_{p^{#phi}_{z}}", "Counts");
   pullsArrayConfig.SetIndividualAxisTitles(23, "Pull_{E^{#phi}}", "Counts");

   // NOWE: Indywidualne konfiguracje dla wybranych pulls
   // Możemy ustawić różne zakresy dla różnych typów pulls
   // Pull 0,1 - pozycja vertex: może mieć większą rozpiętość
   pullsArrayConfig.SetIndividualBinning(0, 50, -3, 3);
   pullsArrayConfig.SetIndividualAxisTitles(0, "Cluster 1 Pull_{x}", "Counts");
   pullsArrayConfig.SetIndividualBinning(1, 50, -3, 3);
   pullsArrayConfig.SetIndividualAxisTitles(1, "Cluster 1 Pull_{y}", "Counts");
   // Pull 2,3 - pęd: standardowy zakres
   pullsArrayConfig.SetIndividualBinning(2, 50, -3, 3);
   pullsArrayConfig.SetIndividualAxisTitles(2, "Cluster 1 Pull_{z}", "Counts");
   // Pull 4 - energia: może być bardziej precyzyjny
   pullsArrayConfig.SetIndividualBinning(4, 50, -3, 3);
   pullsArrayConfig.SetIndividualAxisTitles(4, "Cluster 1 Time Pull", "Counts");
   // Pull 5 - energia: może być bardziej precyzyjny
   pullsArrayConfig.SetIndividualBinning(5, 50, -3, 3);
   pullsArrayConfig.SetIndividualAxisTitles(5, "Cluster 1 Energy Pull", "Counts");
   pullsArrayConfig.SetIndividualBinning(6, 50, -3, 3);
   pullsArrayConfig.SetIndividualAxisTitles(6, "Cluster 2 Pull_{x}", "Counts");

   pullsArrayConfig.SetIndividualBinning(7, 50, -3, 3);
   pullsArrayConfig.SetIndividualAxisTitles(7, "Cluster 2 Pull_{y}", "Counts");

   // Pull 2,3 - pęd: standardowy zakres
   pullsArrayConfig.SetIndividualBinning(8, 50, -3, 3);
   pullsArrayConfig.SetIndividualAxisTitles(8, "Cluster 2 Pull_{z}", "Counts");

   // Pull 4 - energia: może być bardziej precyzyjny
   pullsArrayConfig.SetIndividualBinning(9, 50, -3, 3);
   pullsArrayConfig.SetIndividualAxisTitles(9, "Cluster 2 Time Pull", "Counts");

   // Pull 5 - energia: może być bardziej precyzyjny
   pullsArrayConfig.SetIndividualBinning(10, 50, -3, 3);
   pullsArrayConfig.SetIndividualAxisTitles(10, "Cluster 2 Energy Pull", "Events");

   // Reszta pulls (5-23) będzie używać commonConfig

   histMgr->CreateHistSet1D("invMassKch", invMassKchConfig);
   histMgr->CreateHistSet1D("invMassKne", invMassKneConfig);
   histMgr->CreateHistSet1D("timeDiff", timeDiffConfig);
   histMgr->CreateHistSet1D("invMassKchMC", invMassKchMCConfig);
   histMgr->CreateHistSet1D("invMassKneMC", invMassKneMCConfig);
   histMgr->CreateHistSet1D("timeDiffMC", timeDiffMCConfig);
   histMgr->CreateHistSet1D("chi2TriKinFit", chi2TriKinFitConfig);
   histMgr->CreateHistSet1D("trcFinal", trcFinalConfig);
   histMgr->CreateHistSet1D("curvMC", curvMCConfig);
   histMgr->CreateHistSet1D("phivMC", phivMCConfig);
   histMgr->CreateHistSet1D("cotvMC", cotvMCConfig);
   histMgr->CreateHistSet1D("xchMC", xchMCConfig);
   histMgr->CreateHistSet1D("ychMC", ychMCConfig);
   histMgr->CreateHistSet1D("zchMC", zchMCConfig);
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

   if (*mctruth == 1)
      weight = interf_function(*KaonChTimeCMMC - *KaonNeTimeCMMC);

   Float_t trcSum = trcfinal[0] + trcfinal[1] + trcfinal[2] + trcfinal[3];

   Bool_t photonEneLimit = gammaMomTriKinFit1[3] > 20. && gammaMomTriKinFit2[3] > 20. && gammaMomTriKinFit3[3] > 20. && gammaMomTriKinFit4[3] > 20.;

   if (1)
   {
      if (*mcflag == 1)
      {
         Float_t minv4gam_tri = sqrt(pow(KnetriKinFit[3], 2) - pow(KnetriKinFit[0], 2) - pow(KnetriKinFit[1], 2) - pow(KnetriKinFit[2], 2));

         Float_t distance = sqrt(pow(KneTriangle[6] - ip[0], 2) + pow(KneTriangle[7] - ip[1], 2) + pow(KneTriangle[8] - ip[2], 2)),
                 velocity = cVel * sqrt(pow(KneTriangle[0], 2) + pow(KneTriangle[1], 2) + pow(KneTriangle[2], 2)) / KneTriangle[3],
                 timeOfFlight = (distance / velocity) / (tau_S);

         histMgr->Fill1D("invMassKch", *mctruth, Kchrec[5], weight);
         histMgr->Fill1D("invMassKne", *mctruth, *minv4gam, weight);
         histMgr->Fill1D("timeDiff", *mctruth, *KaonChTimeCM - *KaonNeTimeCM, weight);
         histMgr->Fill1D("chi2TriKinFit", *mctruth, *Chi2TriKinFit, weight);
         histMgr->Fill1D("trcFinal", *mctruth, trcSum, weight);

         histMgr->Fill1D("xchMC", *mctruth, Kchboost[6] - Kchmc[6], weight);
         histMgr->Fill1D("ychMC", *mctruth, Kchboost[7] - Kchmc[7], weight);
         histMgr->Fill1D("zchMC", *mctruth, Kchboost[8] - Kchmc[8], weight);

         if (*mctruth == 1)
         {
            Float_t
                testcharged00 = sqrt(pow(abs(*CurvSmeared1) - abs(CurvMC[0]), 2) + pow(abs(*PhivSmeared1) - abs(PhivMC[0]), 2) + pow(abs(*CotvSmeared1) - abs(CotvMC[0]), 2)),
                testcharged11 = sqrt(pow(abs(*CurvSmeared2) - abs(CurvMC[0]), 2) + pow(abs(*PhivSmeared2) - abs(PhivMC[0]), 2) + pow(abs(*CotvSmeared2) - abs(CotvMC[0]), 2));

            if (testcharged00 < testcharged11)
            {
               histMgr->Fill1D("curvMC", *mctruth, abs(*CurvSmeared1) - abs(CurvMC[0]), weight);
               histMgr->Fill1D("phivMC", *mctruth, abs(*PhivSmeared1) - abs(PhivMC[0]), weight);
               histMgr->Fill1D("cotvMC", *mctruth, abs(*CotvSmeared1) - abs(CotvMC[0]), weight);
            }
            else
            {
               histMgr->Fill1D("curvMC", *mctruth, abs(*CurvSmeared1) - abs(CurvMC[1]), weight);
               histMgr->Fill1D("phivMC", *mctruth, abs(*PhivSmeared1) - abs(PhivMC[1]), weight);
               histMgr->Fill1D("cotvMC", *mctruth, abs(*CotvSmeared1) - abs(CotvMC[1]), weight);
            }
         }

         // Wypełnij array histogramy dla składowych pędu Kchrec
         std::vector<Double_t> momentum = {Kchboost[0], Kchboost[1], Kchboost[2], Kchboost[3]};
         histMgr->FillArrayAll1D("KchrecMomentum", *mctruth, momentum, weight);

         // Wypełnij array histogramy dla składowych pędu Kchrec
         // momentum = {KneTriangle[0] - Knemc[0], KneTriangle[1] - Knemc[1], KneTriangle[2] - Knemc[2], KneTriangle[3] - Knemc[3]};
         momentum = {KnetriKinFit[0], KnetriKinFit[1], KnetriKinFit[2], KnetriKinFit[3]};
         histMgr->FillArrayAll1D("KnerecMomentum", *mctruth, momentum, weight);

         // Wypełnij array histogramy dla składowych neu vtx
         // momentum = {KneTriangle[0] - Knemc[0], KneTriangle[1] - Knemc[1], KneTriangle[2] - Knemc[2], KneTriangle[3] - Knemc[3]};
         momentum = {KnetriKinFit[6], KnetriKinFit[7], KnetriKinFit[8], *KaonNeTimeCM};
         histMgr->FillArrayAll1D("NeuVtx", *mctruth, momentum, weight);

         // // // Wypełnij array histogramy dla pulls (wszystkie 5 składowych)
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

         histMgr->Fill2D("TimeNeutral2d", *mctruth, *KaonNeTimeCMMC, *KaonNeTimeCM, weight);
      }

      if (*mcflag == 0)
      {
         histMgr->FillData1D("invMassKch", Kchrec[5], weight);
         histMgr->FillData1D("invMassKne", *minv4gam, weight);
         histMgr->FillData1D("timeDiff", *KaonChTimeCM - *KaonNeTimeCM, weight);
         histMgr->FillData1D("chi2TriKinFit", *Chi2TriKinFit, weight);

         // Wypełnij array histogramy dla danych - pulls (jeśli dostępne dla danych)
         std::vector<Double_t> pullsData = {pullsTriKinFit[0], pullsTriKinFit[1], pullsTriKinFit[2],
                                            pullsTriKinFit[3], pullsTriKinFit[4], pullsTriKinFit[5], pullsTriKinFit[6], pullsTriKinFit[7],
                                            pullsTriKinFit[8], pullsTriKinFit[9], pullsTriKinFit[10], pullsTriKinFit[11], pullsTriKinFit[12],
                                            pullsTriKinFit[13], pullsTriKinFit[14], pullsTriKinFit[15], pullsTriKinFit[16], pullsTriKinFit[17],
                                            pullsTriKinFit[18], pullsTriKinFit[19], pullsTriKinFit[20], pullsTriKinFit[21], pullsTriKinFit[22],
                                            pullsTriKinFit[23]};
         histMgr->FillArrayAllData1D("PullsTriKinFit", pullsData, weight);

         // Wypełnij array histogramy dla składowych pędu Kchrec
         std::vector<Double_t> momentum = {Kchboost[0], Kchboost[1], Kchboost[2], Kchboost[3]};
         histMgr->FillArrayAllData1D("KchrecMomentum", momentum, weight);

         // Wypełnij array histogramy dla składowych pędu Kchrec
         // momentum = {KneTriangle[0] - Knemc[0], KneTriangle[1] - Knemc[1], KneTriangle[2] - Knemc[2], KneTriangle[3] - Knemc[3]};
         momentum = {KnetriKinFit[0], KnetriKinFit[1], KnetriKinFit[2], KnetriKinFit[3]};
         histMgr->FillArrayAllData1D("KnerecMomentum", momentum, weight);

         // Wypełnij array histogramy dla składowych neu vtx
         // momentum = {KneTriangle[0] - Knemc[0], KneTriangle[1] - Knemc[1], KneTriangle[2] - Knemc[2], KneTriangle[3] - Knemc[3]};
         momentum = {KnetriKinFit[6], KnetriKinFit[7], KnetriKinFit[8], *KaonNeTimeCM};
         histMgr->FillArrayAllData1D("NeuVtx", momentum, weight);
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

   // 1D histogramy z danymi
   histMgr->ScaleChannelByEntries("invMassKch", 1);
   histMgr->ScaleChannelByEntries("invMassKne", 1);
   histMgr->ScaleChannelByEntries("timeDiff", 1);
   histMgr->ScaleChannelByEntries("chi2TriKinFit", 1);
   histMgr->ScaleChannelByEntries("curvMC", 1);
   histMgr->ScaleChannelByEntries("phivMC", 1);
   histMgr->ScaleChannelByEntries("cotvMC", 1);
   histMgr->ScaleChannelByEntries("xchMC", 1);
   histMgr->ScaleChannelByEntries("ychMC", 1);
   histMgr->ScaleChannelByEntries("zchMC", 1);

   histMgr->ScaleChannel2DByEntries("TimeNeutral2d", 1);

   histMgr->ScaleChannelByEntries("invMassKchMC", 1);
   histMgr->ScaleChannelByEntries("invMassKneMC", 1);
   histMgr->ScaleChannelByEntries("timeDiffMC", 1);

   // // Rysuj array histogramy
   // histMgr->ScaleArrayChannelByEntries("KchrecMomentum", 1); // z danymi
   // histMgr->ScaleArrayChannelByEntries("KnerecMomentum", 1); // z danymi
   // histMgr->ScaleArrayChannelByEntries("NeuVtx", 1);         // z danymi
   // histMgr->ScaleArrayChannelByEntries("PullsTriKinFit", 1); // z danymi

   histMgr->SetNormalizationType(HistManager::NormalizationType::SIMPLE_SCALE); // Ustaw proste skalowanie --- IGNORE ---

   // 1D histogramy z danymi
   histMgr->DrawSet1D("invMassKch", "HIST", true);
   histMgr->DrawSet1D("invMassKne", "HIST", true);
   histMgr->DrawSet1D("timeDiff", "HIST", true);
   histMgr->DrawSet1D("chi2TriKinFit", "HIST", true);
   histMgr->DrawSet1D("trcFinal", "HIST", true);

   histMgr->DrawSet1D("curvMC", "HIST", true);
   histMgr->DrawSet1D("phivMC", "HIST", true);
   histMgr->DrawSet1D("cotvMC", "HIST", true);

   histMgr->DrawSet1D("xchMC", "HIST", true);
   histMgr->DrawSet1D("ychMC", "HIST", true);
   histMgr->DrawSet1D("zchMC", "HIST", true);

   histMgr->DrawSet2D("TimeNeutral2d", "COLZ", true);

   histMgr->DrawSet1D("invMassKchMC", "HIST", true);
   histMgr->DrawSet1D("invMassKneMC", "HIST", true);
   histMgr->DrawSet1D("timeDiffMC", "HIST", true);

   // Rysuj array histogramy
   histMgr->DrawArray1D("KchrecMomentum", true); // z danymi
   histMgr->DrawArray1D("KnerecMomentum", true); // z danymi
   histMgr->DrawArray1D("NeuVtx", true);         // z danymi
   histMgr->DrawArray1D("PullsTriKinFit", true); // z danymi

   // histMgr->SaveToRoot("analysis_results.root");

   histMgr->SaveSet("invMassKch", "invMassKch");
   histMgr->SaveSet("invMassKne", "invMassKne");
   histMgr->SaveSet("timeDiff", "timeDiff");
   histMgr->SaveSet("chi2TriKinFit", "chi2TriKinFit");
   histMgr->SaveSet("trcFinal", "trcFinal");

   histMgr->SaveSet("curvMC", "curvMC");
   histMgr->SaveSet("phivMC", "phivMC");
   histMgr->SaveSet("cotvMC", "cotvMC");

   histMgr->SaveSet("xchMC", "xchMC");
   histMgr->SaveSet("ychMC", "ychMC");
   histMgr->SaveSet("zchMC", "zchMC");

   histMgr->SaveSet("invMassKchMC", "invMassKchMC");
   histMgr->SaveSet("invMassKneMC", "invMassKneMC");
   histMgr->SaveSet("timeDiffMC", "timeDiffMC");

   // Zapisz array histogramy
   histMgr->ExportSet("KchrecMomentum", "KchrecMomentum", HistManager::ImageFormat::SVG);
   histMgr->ExportSet("KnerecMomentum", "KnerecMomentum", HistManager::ImageFormat::SVG);
   histMgr->ExportSet("NeuVtx", "NeuVtx", HistManager::ImageFormat::SVG);
   histMgr->ExportSet("PullsTriKinFit", "PullsTriKinFit", HistManager::ImageFormat::SVG);
   histMgr->ExportSet("TimeNeutral2d", "TimeNeutral2d", HistManager::ImageFormat::SVG);

   // Wyniki
   for (size_t i = 0; i < cutter->GetCuts().size(); ++i)
   {
      std::cout << "Cut " << i << ": Eff=" << cutter->GetEfficiency(i) << " +- " << cutter->GetEfficiencyError(i)
                << " Purity=" << cutter->GetPurity(i) << " +- " << cutter->GetPurityError(i)
                << " S/B=" << cutter->GetSignalToBackground(i) << " +- " << cutter->GetSignalToBackgroundError(i) << "\n";
   }

   delete histMgr;
}