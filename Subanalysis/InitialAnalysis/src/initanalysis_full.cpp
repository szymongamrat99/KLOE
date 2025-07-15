#include <iostream>
#include <fstream>

#include <TROOT.h>
#include <TTree.h>
#include <TFile.h>
#include <TChain.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>

#include <event_data.h>
#include <GeneratedVariables.h>
#include <boost/optional.hpp>
#include <SplitFileWriter.h>
#include <charged_mom.h>

#include "../inc/initialanalysis.hpp"
#include "initialanalysis.hpp"

int InitialAnalysis_full(TChain &chain, ErrorHandling::ErrorLogs &logger, KLOE::pm00 &Obj)
{
	TTreeReader reader(&chain);
	GeneralEventProperties generalProps(reader);
	ClusterProperties clusterProps(reader);
	ChargedVertexProperties chVtxProps(reader);
	BhabhaIP bhabhaProps(reader);

	BaseKinematics baseKin;

	GeneratedVariables genVarClassifier;

	// Set flag for initial analysis
	Bool_t MonteCarloInitAnalysis = properties["flags"]["initialAnalysisExec"]["MC"];

	// Which analysis to follow
	KLOE::HypothesisCode hypoCode = Obj.StringToHypothesisCode((std::string)properties["flags"]["analysisCode"]);

	if (hypoCode == KLOE::HypothesisCode::INVALID_VALUE)
		return 1;

	std::string
			base_filename = "mk0_all_phys3",
			dirname = (std::string)initialanalysis_dir + (std::string)root_files_dir,
			dated_folder = Obj.CreateDatedFolder(dirname);

	SplitFileWriter writer(base_filename, 1.5 * 1024 * 1024 * 1024, false, dated_folder);

	Int_t mcflag = 0, mctruth = 0, NCLMIN = 4; // Assuming NCLMIN is 4, adjust as needed;
	std::vector<Int_t> neuclulist;

	UInt_t mctruth_num[8] = {0, 0, 0, 0, 0, 0, 0, 0}; // Array to hold mctruth values

	// Progress bar
	boost::progress_display show_progress(reader.GetEntries());
	// ---------------------------------------------------

	ErrorHandling::ErrorCodes errorCode;

	// Initialization of Charged part of decay reconstruction class
	// Constructor is below, in the loop
	boost::optional<KLOE::ChargedVtxRec<>> eventAnalysis;
	// -------------------------------------------------------------

	Int_t mode = 1; // Model for pi+pi-

	GeneralEventPropertiesMC *eventProps;

	if (MonteCarloInitAnalysis)
		eventProps = new GeneralEventPropertiesMC(reader);

	while (reader.Next())
	{
		// Here you would process each entry in the tree.
		// For example, you can read values from the tree and perform calculations.
		// This is a placeholder for your actual analysis logic.

		// Initial values of mcflag and mctruth
		mcflag = 0;
		mctruth = 0;

		baseKin.vtaken.clear();
		baseKin.vtakenClosest.clear();

		baseKin.Kchrecnew.clear();
		baseKin.KchrecClosest.clear();

		for (Int_t i = 0; i < 2; i++)
		{
			baseKin.trknew[i].clear();
			baseKin.trkClosest[i].clear();
		}

		baseKin.vtaken.resize(3);
		baseKin.vtakenKS.resize(3);
		baseKin.vtakenKL.resize(3);
		baseKin.vtakenClosest.resize(3);

		baseKin.Kchrecnew.resize(9);
		baseKin.KchrecKS.resize(9);
		baseKin.KchrecKL.resize(9);
		baseKin.KchrecClosest.resize(9);
		baseKin.KchrecKLTwoBody.resize(9);

		for (Int_t i = 0; i < 2; i++)
		{
			baseKin.trknew[i].resize(4);
			baseKin.trkKS[i].resize(4);
			baseKin.trkKL[i].resize(4);
			baseKin.trkClosest[i].resize(4);
		}

		if (MonteCarloInitAnalysis)
		{
			mcflag = 1;

			genVarClassifier.classifyChannel(
					*eventProps->ntmc,
					*eventProps->nvtxmc,
					&eventProps->pidmc[0],
					&eventProps->vtxmc[0],
					&eventProps->mother[0],
					mcflag, // Assuming mcflag is 1 for MC events
					mctruth);

			MctruthCounter(mctruth, mctruth_num);
			// -------------------------------------------------------------------
		}

		if (hypoCode == KLOE::HypothesisCode::FOUR_PI) // If we look for pipipipi - clusters do not matter
			errorCode = ErrorHandling::ErrorCodes::NO_ERROR;
		else
		{
			errorCode = genVarClassifier.FindNeutralCluster(*clusterProps.nclu,
																										*clusterProps.ntcl,
																										&clusterProps.asscl[0],
																										NCLMIN,
																										logger,
																										neuclulist);
		}

		if (errorCode != ErrorHandling::ErrorCodes::NO_ERROR)
			logger.getErrLog(errorCode, "", mctruth);
		else
		{
			// Construction of the charged rec class object
			Float_t bhabha_vtx[3] = {*bhabhaProps.x, *bhabhaProps.y, *bhabhaProps.z};

			eventAnalysis.emplace(*chVtxProps.nv, *chVtxProps.ntv, &chVtxProps.iv[0], bhabha_vtx, &chVtxProps.Curv[0], &chVtxProps.Phiv[0], &chVtxProps.Cotv[0], &chVtxProps.xv[0], &chVtxProps.yv[0], &chVtxProps.zv[0], mode);

			// --------------------------------------------------------------------------------
			// Error codes for different hypotheses
			std::map<KLOE::HypothesisCode, ErrorHandling::ErrorCodes> hypoMap;

			// KMASS HYPOTHESIS - FOR SIGNAL
			hypoMap[KLOE::HypothesisCode::SIGNAL] = eventAnalysis->findKchRec(baseKin.Kchrecnew, baseKin.trknew[0], baseKin.trknew[1], baseKin.vtaken, logger);

			// VTX CLOSEST TO BHABHA IP - FOR OMEGAPI
			hypoMap[KLOE::HypothesisCode::OMEGAPI] = eventAnalysis->findKClosestRec(baseKin.KchrecClosest, baseKin.trkClosest[0], baseKin.trkClosest[1], baseKin.vtakenClosest, logger);


			ErrorHandling::ErrorCodes errTmp[2];

			// VTX OF KS - FOR PIPIPIPI
			errTmp[0] = eventAnalysis->findKSLRec(16, -1, baseKin.KchrecKS, baseKin.trkKS[0], baseKin.trkKS[1], baseKin.vtakenKS, logger);

			// VTX OF KL - FOR PIPIPIPI
			errTmp[1] = eventAnalysis->findKSLRec(10, baseKin.vtakenKS[0], baseKin.KchrecKL, baseKin.trkKL[0], baseKin.trkKL[1], baseKin.vtakenKL, logger);
			// --------------------------------------------------------------------------------

			if	(errTmp[0] != ErrorHandling::ErrorCodes::NO_ERROR)
				hypoMap[KLOE::HypothesisCode::FOUR_PI] = errTmp[0];
			else if (errTmp[1] != ErrorHandling::ErrorCodes::NO_ERROR)
				hypoMap[KLOE::HypothesisCode::FOUR_PI] = errTmp[1];
			else
				hypoMap[KLOE::HypothesisCode::FOUR_PI] = ErrorHandling::ErrorCodes::NO_ERROR;

			errorCode = hypoMap[hypoCode]; // error code based on the hypothesis

			if (errorCode != ErrorHandling::ErrorCodes::NO_ERROR)
				logger.getErrLog(errorCode, "", mctruth);
			else
			{
				if (abs(baseKin.KchrecKS[5] - mK0) < 5 && abs(baseKin.KchrecKL[5] - mK0) < 5)
				{
					errorCode = ErrorHandling::ErrorCodes::NO_ERROR;

					// Clone of the branches of the old tree
					// General properties of the event
					baseKin.nrun = *generalProps.nrun;
					baseKin.nev = *generalProps.nev;

					baseKin.necls = *generalProps.necls;
					baseKin.eclfilfo = *generalProps.eclfilfo;
					baseKin.eclfilfoword = *generalProps.eclfilfoword;

					baseKin.eclstream.assign(generalProps.eclstream.begin(), generalProps.eclstream.end());
					// -------------------------------------------------------------------------------------
					// Bhabha interaction point and momentum
					baseKin.Bx = *bhabhaProps.x;
					baseKin.By = *bhabhaProps.y;
					baseKin.Bz = *bhabhaProps.z;
					baseKin.Bpx = *bhabhaProps.px;
					baseKin.Bpy = *bhabhaProps.py;
					baseKin.Bpz = *bhabhaProps.pz;
					baseKin.Broots = *bhabhaProps.energy;
					// -------------------------------------------------------------------------------------
					// Cluster data
					baseKin.nclu = *clusterProps.nclu;
					baseKin.ntcl = *clusterProps.ntcl;
					baseKin.T0step1 = *clusterProps.t0step1;
					baseKin.Asscl.assign(clusterProps.asscl.begin(), clusterProps.asscl.end());
					baseKin.Xcl.assign(clusterProps.xcl.begin(), clusterProps.xcl.end());
					baseKin.Ycl.assign(clusterProps.ycl.begin(), clusterProps.ycl.end());
					baseKin.Zcl.assign(clusterProps.zcl.begin(), clusterProps.zcl.end());
					baseKin.Tcl.assign(clusterProps.tcl.begin(), clusterProps.tcl.end());
					baseKin.Enecl.assign(clusterProps.enecl.begin(), clusterProps.enecl.end());
					// -------------------------------------------------------------------------------------
					// Charged decay data
					baseKin.nv = *chVtxProps.nv;
					baseKin.ntv = *chVtxProps.ntv;
					baseKin.iv.assign(chVtxProps.iv.begin(), chVtxProps.iv.end());
					baseKin.Curv.assign(chVtxProps.Curv.begin(), chVtxProps.Curv.end());
					baseKin.Phiv.assign(chVtxProps.Phiv.begin(), chVtxProps.Phiv.end());
					baseKin.Cotv.assign(chVtxProps.Cotv.begin(), chVtxProps.Cotv.end());
					baseKin.xv.assign(chVtxProps.xv.begin(), chVtxProps.xv.end());
					baseKin.yv.assign(chVtxProps.yv.begin(), chVtxProps.yv.end());
					baseKin.zv.assign(chVtxProps.zv.begin(), chVtxProps.zv.end());
					// -------------------------------------------------------------------------------------
					// Monte carlo data
					if (MonteCarloInitAnalysis)
					{
						baseKin.ntmc = *eventProps->ntmc;
						baseKin.nvtxmc = *eventProps->nvtxmc;
						baseKin.vtxmc.assign(eventProps->vtxmc.begin(), eventProps->vtxmc.end());
						baseKin.pidmc.assign(eventProps->pidmc.begin(), eventProps->pidmc.end());
						baseKin.mother.assign(eventProps->mother.begin(), eventProps->mother.end());
						baseKin.xvmc.assign(eventProps->xvmc.begin(), eventProps->xvmc.end());
						baseKin.yvmc.assign(eventProps->yvmc.begin(), eventProps->yvmc.end());
						baseKin.zvmc.assign(eventProps->zvmc.begin(), eventProps->zvmc.end());
						baseKin.pxmc.assign(eventProps->pxmc.begin(), eventProps->pxmc.end());
						baseKin.pymc.assign(eventProps->pymc.begin(), eventProps->pymc.end());
						baseKin.pzmc.assign(eventProps->zvmc.begin(), eventProps->zvmc.end());
					}
					else
					{
						baseKin.ntmc = 0;
						baseKin.nvtxmc = 0;
						baseKin.vtxmc = {};
						baseKin.pidmc = {};
						baseKin.mother = {};
						baseKin.xvmc = {};
						baseKin.yvmc = {};
						baseKin.zvmc = {};
						baseKin.pxmc = {};
						baseKin.pymc = {};
						baseKin.pzmc = {};
					}
					// -------------------------------------------------------------------------------------

					// Int_t zmienne
					std::map<std::string, Int_t> intVars = {
							{"nrun", baseKin.nrun},									// Number of run
							{"nev", baseKin.nev},										// Number of event
							{"necls", baseKin.necls},								// Number of ECL words
							{"Eclfilfo", baseKin.eclfilfo},					// Which filfo was used
							{"Eclfilfoword", baseKin.eclfilfoword}, // Filfo word
							{"mcflag", mcflag},											// If event from MC of Data
							{"mctruth", mctruth},										// What event type
							{"nclu", baseKin.nclu},
							{"ntcl", baseKin.ntcl},
							{"nv", baseKin.nv},
							{"ntv", baseKin.ntv},
							{"ntmc", baseKin.ntmc},
							{"nvtxmc", baseKin.nvtxmc}};

					// Float_t zmienne
					std::map<std::string, Float_t> floatVars = {
							{"T0step1", baseKin.T0step1},
							{"Bx", baseKin.Bx},
							{"By", baseKin.By},
							{"Bz", baseKin.Bz},
							{"Bpx", baseKin.Bpx},
							{"Bpy", baseKin.Bpy},
							{"Bpz", baseKin.Bpz},
							{"Broots", baseKin.Broots}};

					// Tablice
					std::map<std::string, std::vector<Int_t>> intArrays = {
							{"eclstream", baseKin.eclstream},
							{"Asscl", baseKin.Asscl},
							{"iv", baseKin.iv},
							{"vtxmc", baseKin.vtxmc},
							{"pidmc", baseKin.pidmc},
							{"mother", baseKin.mother},
							{"vtakenClosest", baseKin.vtakenClosest},
							{"vtaken", baseKin.vtaken}};

					std::map<std::string, std::vector<Float_t>> floatArrays = {
							{"Xcl", baseKin.Xcl},
							{"Ycl", baseKin.Ycl},
							{"Zcl", baseKin.Zcl},
							{"Tcl", baseKin.Tcl},
							{"Enecl", baseKin.Enecl},
							{"Curv", baseKin.Curv},
							{"Phiv", baseKin.Phiv},
							{"Cotv", baseKin.Cotv},
							{"xv", baseKin.xv},
							{"yv", baseKin.yv},
							{"zv", baseKin.zv},
							{"xvmc", baseKin.xvmc},
							{"yvmc", baseKin.yvmc},
							{"zvmc", baseKin.zvmc},
							{"pxmc", baseKin.pxmc},
							{"pymc", baseKin.pymc},
							{"pzmc", baseKin.pzmc},
							{"KchrecClosest", baseKin.KchrecClosest},
							{"trk1Closest", baseKin.trkClosest[0]},
							{"trk2Closest", baseKin.trkClosest[1]},
							{"Kchrec", baseKin.Kchrecnew},
							{"trk1", baseKin.trknew[0]},
							{"trk2", baseKin.trknew[1]},
							{"KchrecKS", baseKin.KchrecKS},
							{"trk1KS", baseKin.trkKS[0]},
							{"trk2KS", baseKin.trkKS[1]},
							{"KchrecKL", baseKin.KchrecKL},
							{"trk1KL", baseKin.trkKL[0]},
							{"trk2KL", baseKin.trkKL[1]}};

					writer.Fill(intVars, floatVars, intArrays, floatArrays);
				}
				else
				{
					errorCode = ErrorHandling::ErrorCodes::CHARGED_KAON_MASS_PRE;
				}
			}

			// ------------------------------------------------------------------
		}

		++show_progress; // Progress of the loading bar
	}

	std::map<ErrorHandling::ErrorCodes, int> physicsErrorCountsPerMctruth[8];

	for (int i = 0; i < 8; ++i)
	{
		physicsErrorCountsPerMctruth[i] = logger.getPhysicsErrorCountsForMctruth(i);
	};

	logger.printPhysicsErrorStatsPerMctruth(false);

	writer.Close();

	return 0;
}

void MctruthCounter(Int_t mctruth, UInt_t mctruth_num[8])
{
	switch (mctruth)
	{
	case 0:
		mctruth_num[0]++;
		break;
	case 1:
		mctruth_num[1]++;
		break;
	case 2:
		mctruth_num[2]++;
		break;
	case 3:
		mctruth_num[3]++;
		break;
	case 4:
		mctruth_num[4]++;
		break;
	case 5:
		mctruth_num[5]++;
		break;
	case 6:
		mctruth_num[6]++;
		break;
	case 7:
		mctruth_num[7]++;
		break;
	default:
		std::cerr << "Unknown mctruth value: " << mctruth << std::endl;
	}
}