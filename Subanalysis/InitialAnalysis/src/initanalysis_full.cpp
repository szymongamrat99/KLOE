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
#include <StatisticalCutter.h>

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

	std::ifstream file(cutlimitsName);
	json j = json::parse(file);

	StatisticalCutter cutter(cutlimitsName, 7);

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

	Float_t
			KchrecKSMom = 0, KchrecKLMom = 0;

	cutter.RegisterVariableGetter("InvMassKch", [&]()
																{ return baseKin.KchrecKS[5]; });
	cutter.RegisterVariableGetter("InvMassKne", [&]()
																{ return baseKin.KchrecKL[5]; });
	cutter.RegisterVariableGetter("TwoBodyMomKS", [&]()
																{ return KchrecKSMom; });
	cutter.RegisterVariableGetter("TwoBodyMomKL", [&]()
																{ return KchrecKLMom; });

	while (reader.Next())
	{
		// Here you would process each entry in the tree.
		// For example, you can read values from the tree and perform calculations.
		// This is a placeholder for your actual analysis logic.

		// Initial values of mcflag and mctruth
		mcflag = 0;
		mctruth = 0;

		baseKin.vtaken.clear();
		baseKin.vtakenKS.clear();
		baseKin.vtakenKL.clear();
		baseKin.vtakenClosest.clear();

		baseKin.Kchrecnew.clear();
		baseKin.KchrecKS.clear();
		baseKin.KchrecKL.clear();
		baseKin.KchrecClosest.clear();

		baseKin.KchboostKS.clear();
		baseKin.KchboostKL.clear();

		baseKin.ipKS.clear();
		baseKin.ipKL.clear();

		for (Int_t i = 0; i < 2; i++)
		{
			baseKin.trknew[i].clear();
			baseKin.trkKS[i].clear();
			baseKin.trkKL[i].clear();
			baseKin.trkClosest[i].clear();
			baseKin.trkKLmc[i].clear();
			baseKin.trkKSmc[i].clear();
		}

		baseKin.vtaken.resize(3);
		baseKin.vtakenKS.resize(3);
		baseKin.vtakenKL.resize(3);
		baseKin.vtakenClosest.resize(3);

		baseKin.Kchrecnew.resize(9);
		baseKin.KchrecKS.resize(9);
		baseKin.KchrecKL.resize(9);
		baseKin.KchrecClosest.resize(9);

		baseKin.KchboostKS.resize(9);
		baseKin.KchboostKL.resize(9);

		baseKin.ipKS.resize(3);
		baseKin.ipKL.resize(3);

		for (Int_t i = 0; i < 2; i++)
		{
			baseKin.trknew[i].resize(4);
			baseKin.trkKS[i].resize(4);
			baseKin.trkKL[i].resize(4);
			baseKin.trkClosest[i].resize(4);
		}

		std::vector<std::vector<Float_t>>
				trkMC,
				pgammaMC,
				clusterMC;

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

			genVarClassifier.genVars(*eventProps->ntmc,
															 *eventProps->nvtxmc,
															 *clusterProps.nclu,
															 &eventProps->pidmc[0],
															 &eventProps->vtxmc[0],
															 &eventProps->mother[0],
															 &eventProps->xvmc[0],
															 &eventProps->yvmc[0],
															 &eventProps->zvmc[0],
															 &eventProps->pxmc[0],
															 &eventProps->pymc[0],
															 &eventProps->pzmc[0],
															 mcflag,
															 mctruth,
															 baseKin.ipmc,
															 baseKin.Knemc,
															 baseKin.Kchmc,
															 trkMC,
															 4,
															 pgammaMC,
															 baseKin.goodClusIndex,
															 clusterMC);

			if (mctruth == 7)
			{
				for (Int_t iter = 0; iter < 4; iter++)
				{
					if (trkMC[iter][4] == 10)
					{
						if (baseKin.trkKLmc[0].size() == 0)
							baseKin.trkKLmc[0].assign(trkMC[iter].begin(), trkMC[iter].end() - 1);
						else
							baseKin.trkKLmc[1].assign(trkMC[iter].begin(), trkMC[iter].end() - 1);
					}
					else if (trkMC[iter][4] == 16)
					{
						if (baseKin.trkKSmc[0].size() == 0)
							baseKin.trkKSmc[0].assign(trkMC[iter].begin(), trkMC[iter].end() - 1);
						else
							baseKin.trkKSmc[1].assign(trkMC[iter].begin(), trkMC[iter].end() - 1);
					}
				}
			}
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

			if (errTmp[0] != ErrorHandling::ErrorCodes::NO_ERROR)
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
				Bool_t cutCombined = false;
				std::vector<Bool_t> cutOrdered;

				if (hypoCode == KLOE::HypothesisCode::FOUR_PI)
				{
					cutOrdered.push_back(abs(baseKin.KchrecKS[5] - mK0) < 5);
					cutOrdered.push_back(abs(baseKin.KchrecKL[5] - mK0) < 5);

					if (cutOrdered[0] && cutOrdered[1])
					{

						Float_t
								pKTwoBody = Obj.TwoBodyDecayMass(mPhi, mK0, mK0),
								boostPhi[3] = {
										-*bhabhaProps.px / *bhabhaProps.energy,
										-*bhabhaProps.py / *bhabhaProps.energy,
										-*bhabhaProps.pz / *bhabhaProps.energy},
								phiMom[4] = {*bhabhaProps.px, *bhabhaProps.py, *bhabhaProps.pz, *bhabhaProps.energy}, trkKS_PhiCM[2][4] = {}, KchrecKS_PhiCM[4] = {}, trkKL_PhiCM[2][4], KchrecKL_PhiCM[4] = {};

						Obj.lorentz_transf(boostPhi, baseKin.trkKS[0].data(), trkKS_PhiCM[0]);
						Obj.lorentz_transf(boostPhi, baseKin.trkKS[1].data(), trkKS_PhiCM[1]);
						Obj.lorentz_transf(boostPhi, baseKin.trkKL[0].data(), trkKL_PhiCM[0]);
						Obj.lorentz_transf(boostPhi, baseKin.trkKL[1].data(), trkKL_PhiCM[1]);

						for (Int_t part = 0; part < 2; part++)
							for (Int_t comp = 0; comp < 4; comp++)
							{
								KchrecKS_PhiCM[comp] += trkKS_PhiCM[part][comp];
								KchrecKL_PhiCM[comp] += trkKL_PhiCM[part][comp];
							}

						KchrecKSMom = sqrt(pow(KchrecKS_PhiCM[0], 2) + pow(KchrecKS_PhiCM[1], 2) + pow(KchrecKS_PhiCM[2], 2));
						KchrecKLMom = sqrt(pow(KchrecKL_PhiCM[0], 2) + pow(KchrecKL_PhiCM[1], 2) + pow(KchrecKL_PhiCM[2], 2));

						cutOrdered.push_back(abs(KchrecKSMom - pKTwoBody) < 10);
						cutOrdered.push_back(abs(KchrecKLMom - pKTwoBody) < 10);

						eventAnalysis->KaonMomFromBoost(baseKin.KchrecKS, phiMom, baseKin.KchboostKS);
						eventAnalysis->KaonMomFromBoost(baseKin.KchrecKL, phiMom, baseKin.KchboostKL);

						Float_t X_lineKS[3] = {baseKin.KchboostKS[6],
																	 baseKin.KchboostKS[7],
																	 baseKin.KchboostKS[8]}, // Vertex laying on the line
								X_lineKL[3] = {baseKin.KchboostKL[6],
															 baseKin.KchboostKL[7],
															 baseKin.KchboostKL[8]}, // Vertex laying on the line
								pKS[3] = {baseKin.KchboostKS[0],
													baseKin.KchboostKS[1],
													baseKin.KchboostKS[2]}, // Direction of the line
								pKL[3] = {baseKin.KchboostKL[0],
													baseKin.KchboostKL[1],
													baseKin.KchboostKL[2]}, // Direction of the line
								xB[3] = {baseKin.bhabha_vtx[0],
												 baseKin.bhabha_vtx[1],
												 baseKin.bhabha_vtx[2]}, // Bhabha vertex - laying on the plane
								plane_perp[3] = {0.,
																 baseKin.phi_mom[1],
																 0.}; // Vector perpendicular to the plane from Bhabha momentum

						// Corrected IP event by event
						eventAnalysis->IPBoostCorr(X_lineKS, pKL, xB, plane_perp, baseKin.ipKS);
						eventAnalysis->IPBoostCorr(X_lineKL, pKL, xB, plane_perp, baseKin.ipKL);

						baseKin.ipKS[0] = baseKin.bhabha_vtx[0];
						baseKin.ipKS[1] = baseKin.bhabha_vtx[1];
						// z coordinate of the IP is set to the Bhabha vertex z coordinate if it differs by more than 2 cm
						if (abs(baseKin.ipKS[2] - baseKin.bhabha_vtx[2]) > 2.0)
							baseKin.ipKS[2] = baseKin.bhabha_vtx[2];

						baseKin.ipKL[0] = baseKin.bhabha_vtx[0];
						baseKin.ipKL[1] = baseKin.bhabha_vtx[1];
						// z coordinate of the IP is set to the Bhabha vertex z coordinate if it differs by more than 2 cm
						if (abs(baseKin.ipKL[2] - baseKin.bhabha_vtx[2]) > 2.0)
							baseKin.ipKL[2] = baseKin.bhabha_vtx[2];

						Float_t
								PhiMom[3] = {*bhabhaProps.px, *bhabhaProps.py, *bhabhaProps.pz},
								MissMomKS[3] = {},
								MissMomKL[3] = {},
								PmissKS = 0,
								PmissKL = 0,
								EmissKS = 0,
								EmissKL = 0;

						for (Int_t comp = 0; comp < 3; comp++)
						{
							MissMomKS[comp] = PhiMom[comp] - baseKin.KchboostKS[comp] - baseKin.KchrecKL[comp];
							MissMomKL[comp] = PhiMom[comp] - baseKin.KchboostKL[comp] - baseKin.KchrecKS[comp];
						}

						PmissKS = sqrt(pow(MissMomKS[0], 2) + pow(MissMomKS[1], 2) + pow(MissMomKS[2], 2));
						PmissKL = sqrt(pow(MissMomKL[0], 2) + pow(MissMomKL[1], 2) + pow(MissMomKL[2], 2));

						EmissKS = baseKin.KchboostKS[3] - baseKin.KchrecKS[3];
						EmissKL = baseKin.KchboostKL[3] - baseKin.KchrecKL[3];

						cutOrdered.push_back(sqrt(pow(EmissKS, 2) + pow(PmissKS, 2)) < 10);
						cutOrdered.push_back(sqrt(pow(EmissKL, 2) + pow(PmissKL, 2)) < 10);

						cutOrdered.push_back((pow(EmissKL, 2) - pow(PmissKL, 2) < 10) && (pow(EmissKL, 2) - pow(PmissKL, 2) > -50));

						cutOrdered.push_back((pow(EmissKS, 2) - pow(PmissKS, 2) < 10) && (pow(EmissKS, 2) - pow(PmissKS, 2) > -50));

						cutCombined = cutOrdered[0] && cutOrdered[1] && cutOrdered[2] && cutOrdered[3] && cutOrdered[4] && cutOrdered[5] && cutOrdered[6] && cutOrdered[7];

						cutter.UpdateStats(mctruth);

						baseKin.cuts.clear();
						baseKin.cuts.resize(cutOrdered.size());

						for (Int_t iter = 0; iter < cutOrdered.size(); iter++)
							if (cutOrdered[iter])
								baseKin.cuts[iter] = 1;
							else if (!cutOrdered[iter])
							{
								baseKin.cuts[iter] = 0;

								if (mctruth == 7)
									mctruth = 0;
							}
					}
				}

				if (cutCombined || mctruth == 7 || mctruth == 0)
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
							{"vtaken", baseKin.vtaken},
							{"cutsApplied", baseKin.cuts}};

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
							{"trk2KL", baseKin.trkKL[1]},
							{"KchboostKS", baseKin.KchboostKS},
							{"KchboostKL", baseKin.KchboostKL},
							{"ipKS", baseKin.ipKS},
							{"ipKL", baseKin.ipKL},
							{"ipmc", baseKin.ipmc},
							{"Kchmc", baseKin.Kchmc},
							{"Knemc", baseKin.Knemc},
							{"trk1KSmc", baseKin.trkKSmc[0]},
							{"trk2KSmc", baseKin.trkKSmc[1]},
							{"trk1KLmc", baseKin.trkKLmc[0]},
							{"trk2KLmc", baseKin.trkKLmc[1]}};

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

	// Wyniki
	for (size_t i = 0; i < 4; ++i)
	{
		std::cout << "Cut " << i << ": Eff=" << cutter.GetEfficiency(i)
							<< " Purity=" << cutter.GetPurity(i)
							<< " S/B=" << cutter.GetSignalToBackground(i) << "\n";
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