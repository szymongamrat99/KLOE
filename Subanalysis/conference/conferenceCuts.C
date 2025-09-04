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
	QmissConfig.xtitle = "m^{inv}_{K#rightarrow#pi^{+}#pi^{-}} [MeV/c^{2}]";
	QmissConfig.ytitle = "Counts/2";
	QmissConfig.bins = 50;
	QmissConfig.xmin = 0;
	QmissConfig.xmax = 20;
	QmissConfig.logy = false;
	QmissConfig.showStats = false;

	histMgr->CreateHistSet1D("invMassKch", invMassKchConfig);
	histMgr->CreateHistSet1D("Qmiss", QmissConfig);
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
		if (*mcflag == 1 && *mctruth_int != 0 && *mctruth_int != 2)
		{

			// if (*mctruth == 1 || *mctruth == 2)
			// 	weight = interf_function(*Dtmc);

			histMgr->Fill1D("invMassKch", mctruth_corr, Kchrec[5], weight);
			histMgr->Fill1D("Qmiss", mctruth_corr, *Qmiss, weight);
			histMgr->Fill1D("Chi2", mctruth_corr, *Qmiss, weight);
		}
		else if (*mcflag == 0)
		{
			histMgr->FillData1D("invMassKch", Kchrec[5]);
			histMgr->FillData1D("Qmiss", *Qmiss);
			
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

	// 1D histogramy z danymi
	histMgr->DrawSet1D("invMassKch", "HIST", true);
	histMgr->DrawSet1D("Qmiss", "HIST", true);

	// histMgr->SaveToRoot("analysis_results.root");

	histMgr->SaveSet("invMassKch", "invMassKch");
	histMgr->SaveSet("Qmiss", "Qmiss");

	// Wyniki
	for (size_t i = 0; i < cutter->GetCuts().size(); ++i)
	{
		std::cout << "Cut " << i << ": Eff=" << cutter->GetEfficiency(i) << " +- " << cutter->GetEfficiencyError(i)
				  << " Purity=" << cutter->GetPurity(i) << " +- " << cutter->GetPurityError(i)
				  << " S/B=" << cutter->GetSignalToBackground(i) << " +- " << cutter->GetSignalToBackgroundError(i) << "\n";
	}
}