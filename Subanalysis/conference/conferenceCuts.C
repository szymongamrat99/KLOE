#define conferenceCuts_cxx
// The class definition in conferenceCuts.h has been generated automatically
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
// root> T->Process("conferenceCuts.C")
// root> T->Process("conferenceCuts.C","some options")
// root> T->Process("conferenceCuts.C+")
//

#include "conferenceCuts.h"
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
HistManager::HistConfig QmissConfig;
HistManager::HistConfig chi2KinFitConfig;
HistManager::HistConfig fittedPi0MassesConfig;
HistManager::HistConfig trcvSumConfig;
HistManager::HistConfig pi0Mass1Config;
HistManager::HistConfig pi0Mass2Config;

TTree *tree;
TFile *file;

Double_t
	CombinedPi0Masses = 0,
	trcvsum = 0;

void conferenceCuts::Begin(TTree * /*tree*/)
{
	// The Begin() function is called at the start of the query.
	// When running with PROOF Begin() is only called on the client.
	// The tree argument is deprecated (on PROOF 0 is passed).

	TString option = GetOption();

	histMgr = new HistManager(channNum, channColor, channNames, kFullCircle, kBlack, kOrange);

	std::string cutFileName = "/data/ssd/gamrat/KLOE/Subanalysis/Properties/cut-limits-final.json";

	cutter = new StatisticalCutter(cutFileName, 1, KLOE::HypothesisCode::SIGNAL);

	cutter->RegisterVariableGetter("InvMassKch", [&]()
								   { return Kchrec[5]; });
	cutter->RegisterCentralValueGetter("InvMassKch", [&]()
									   { return mK0; });

	cutter->RegisterVariableGetter("Qmiss", [&]()
								   { return *Qmiss; });

	cutter->RegisterVariableGetter("Chi2KinFit", [&]()
								   { return *Chi2; });

	cutter->RegisterVariableGetter("FittedPi0Masses", [&]()
								   { return CombinedPi0Masses; });

	cutter->RegisterVariableGetter("TrcvSum", [&]()
								   { return trcvsum; });

	cutter->RegisterVariableGetter("InvMassKne", [&]()
								   { return *minv4gam; });
	cutter->RegisterCentralValueGetter("InvMassKne", [&]()
									   { return mK0; });

	invMassKchConfig.name = "invMassKch";
	invMassKchConfig.xtitle = "m^{inv}_{K#rightarrow#pi^{+}#pi^{-}} [MeV/c^{2}]";
	invMassKchConfig.ytitle = "Counts/2";
	invMassKchConfig.bins = 50;
	invMassKchConfig.xmin = 490;
	invMassKchConfig.xmax = 505;
	invMassKchConfig.logy = false;
	invMassKchConfig.showStats = false;

	QmissConfig.name = "Qmiss";
	QmissConfig.xtitle = "Q_{miss} [MeV/c^{2}]";
	QmissConfig.ytitle = "Counts/2";
	QmissConfig.bins = 50;
	QmissConfig.xmin = 0;
	QmissConfig.xmax = 6;
	QmissConfig.logy = false;
	QmissConfig.showStats = false;

	chi2KinFitConfig.name = "chi2KinFit";
	chi2KinFitConfig.xtitle = "#chi^{2} Kinematic Fit";
	chi2KinFitConfig.ytitle = "Counts";
	chi2KinFitConfig.bins = 50;
	chi2KinFitConfig.xmin = 0;
	chi2KinFitConfig.xmax = 60;
	chi2KinFitConfig.logy = false;
	chi2KinFitConfig.showStats = false;

	fittedPi0MassesConfig.name = "fittedPi0Masses";
	fittedPi0MassesConfig.xtitle = "Combined #pi^{0} Mass Deviation [MeV/c^{2}]";
	fittedPi0MassesConfig.ytitle = "Counts";
	fittedPi0MassesConfig.bins = 50;
	fittedPi0MassesConfig.xmin = 0;
	fittedPi0MassesConfig.xmax = 50;
	fittedPi0MassesConfig.logy = false;
	fittedPi0MassesConfig.showStats = false;

	trcvSumConfig.name = "trcvSum";
	trcvSumConfig.xtitle = "Time deviation of neutral vertex [ns]";
	trcvSumConfig.ytitle = "Counts";
	trcvSumConfig.bins = 50;
	trcvSumConfig.xmin = -0.5;
	trcvSumConfig.xmax = 0.5;
	trcvSumConfig.logy = false;
	trcvSumConfig.showStats = false;

	invMassKneConfig.name = "invMassKne";
	invMassKneConfig.xtitle = "m^{inv}_{K#rightarrow#pi^{0}#pi^{0}} [MeV/c^{2}]";
	invMassKneConfig.ytitle = "Counts";
	invMassKneConfig.bins = 50;
	invMassKneConfig.xmin = 400;
	invMassKneConfig.xmax = 600;
	invMassKneConfig.logy = false;
	invMassKneConfig.showStats = false;

	pi0Mass1Config.name = "pi0Mass1";
	pi0Mass1Config.xtitle = "m^{inv}_{#pi^{0}} first [MeV/c^{2}]";
	pi0Mass1Config.ytitle = "Counts";
	pi0Mass1Config.bins = 50;
	pi0Mass1Config.xmin = 120;
	pi0Mass1Config.xmax = 150;
	pi0Mass1Config.logy = false;
	pi0Mass1Config.showStats = false;

	pi0Mass2Config.name = "pi0Mass2";
	pi0Mass2Config.xtitle = "m^{inv}_{#pi^{0}} second [MeV/c^{2}]";
	pi0Mass2Config.ytitle = "Counts";
	pi0Mass2Config.bins = 50;
	pi0Mass2Config.xmin = 120;
	pi0Mass2Config.xmax = 150;
	pi0Mass2Config.logy = false;
	pi0Mass2Config.showStats = false;

	histMgr->CreateHistSet1D("invMassKch", invMassKchConfig);
	histMgr->CreateHistSet1D("Qmiss", QmissConfig);
	histMgr->CreateHistSet1D("chi2KinFit", chi2KinFitConfig);
	histMgr->CreateHistSet1D("fittedPi0Masses", fittedPi0MassesConfig);
	histMgr->CreateHistSet1D("trcvSum", trcvSumConfig);
	histMgr->CreateHistSet1D("invMassKne", invMassKneConfig);
	histMgr->CreateHistSet1D("pi0Mass1", pi0Mass1Config);
	histMgr->CreateHistSet1D("pi0Mass2", pi0Mass2Config);

}

void conferenceCuts::SlaveBegin(TTree * /*tree*/)
{
	// The SlaveBegin() function is called after the Begin() function.
	// When running with PROOF SlaveBegin() is called on each slave server.
	// The tree argument is deprecated (on PROOF 0 is passed).

	TString option = GetOption();
}

Bool_t conferenceCuts::Process(Long64_t entry)
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

	Int_t mctruth_corr;

	fReader.SetLocalEntry(entry);

	CombinedPi0Masses = sqrt(pow(pi0fit[0] - mPi0, 2) + pow(pi0fit[1] - mPi0, 2));

	trcvsum = trcv[g4taken[0]] + trcv[g4taken[1]] + trcv[g4taken[2]] + trcv[g4taken[3]];

	Float_t weight = 1.0;

	if (*mctruth_int > 1)
		mctruth_corr = *mctruth_int - 1;
	else
		mctruth_corr = *mctruth_int;

	if (cutter->PassAllCuts())
	{
		if (*mcflag == 1 && *mctruth_int != 0)
		{
			histMgr->Fill1D("invMassKch", mctruth_corr, Kchrec[5], weight);
			histMgr->Fill1D("Qmiss", mctruth_corr, *Qmiss, weight);
			histMgr->Fill1D("chi2KinFit", mctruth_corr, *Chi2, weight);
			histMgr->Fill1D("fittedPi0Masses", mctruth_corr, CombinedPi0Masses, weight);
			histMgr->Fill1D("trcvSum", mctruth_corr, trcvsum, weight);
			histMgr->Fill1D("invMassKne", mctruth_corr, *minv4gam, weight);
			histMgr->Fill1D("pi0Mass1", mctruth_corr, pi0fit[0], weight);
			histMgr->Fill1D("pi0Mass2", mctruth_corr, pi0fit[1], weight);
		}
		else if (*mcflag == 0)
		{
			histMgr->FillData1D("invMassKch", Kchrec[5]);
			histMgr->FillData1D("Qmiss", *Qmiss);
			histMgr->FillData1D("chi2KinFit", *Chi2);
			histMgr->FillData1D("fittedPi0Masses", CombinedPi0Masses);
			histMgr->FillData1D("trcvSum", trcvsum);
			histMgr->FillData1D("invMassKne", *minv4gam);
			histMgr->FillData1D("pi0Mass1", pi0fit[0]);
			histMgr->FillData1D("pi0Mass2", pi0fit[1]);
		}
	}

	if (*mcflag == 1 && *mctruth_int != 0)
		cutter->UpdateStats(mctruth_corr);

	return kTRUE;
}

void conferenceCuts::SlaveTerminate()
{
	// The SlaveTerminate() function is called after all entries or objects
	// have been processed. When running with PROOF SlaveTerminate() is called
	// on each slave server.
}

void conferenceCuts::Terminate()
{
	// The Terminate() function is the last function to be called during
	// a query. It always runs on the client, it can be used to present
	// the results graphically or save the results to file.

	histMgr->SetNormalizationType(HistManager::NormalizationType::SIMPLE_SCALE);

	// 1D histogramy z danymi
	histMgr->DrawSet1D("invMassKch", "HIST", true);
	histMgr->DrawSet1D("Qmiss", "HIST", true);
	histMgr->DrawSet1D("chi2KinFit", "HIST", true);
	histMgr->DrawSet1D("fittedPi0Masses", "HIST", true);
	histMgr->DrawSet1D("trcvSum", "HIST", true);
	histMgr->DrawSet1D("invMassKne", "HIST", true);
	histMgr->DrawSet1D("pi0Mass1", "HIST", true);
	histMgr->DrawSet1D("pi0Mass2", "HIST", true);

	// histMgr->SaveToRoot("analysis_results.root");

	histMgr->SaveSet("invMassKch", "invMassKch");
	histMgr->SaveSet("Qmiss", "Qmiss");
	histMgr->SaveSet("chi2KinFit", "chi2KinFit");
	histMgr->SaveSet("fittedPi0Masses", "fittedPi0Masses");
	histMgr->SaveSet("trcvSum", "trcvSum");
	histMgr->SaveSet("invMassKne", "invMassKne");
	histMgr->SaveSet("pi0Mass1", "pi0Mass1");
	histMgr->SaveSet("pi0Mass2", "pi0Mass2");

	// Wyniki
	for (size_t i = 0; i < cutter->GetCuts().size(); ++i)
	{
		std::cout << "Cut " << i << ": Eff=" << cutter->GetEfficiency(i) << " +- " << cutter->GetEfficiencyError(i)
				  << " Purity=" << cutter->GetPurity(i) << " +- " << cutter->GetPurityError(i)
				  << " S/B=" << cutter->GetSignalToBackground(i) << " +- " << cutter->GetSignalToBackgroundError(i) << "\n";
	}
}