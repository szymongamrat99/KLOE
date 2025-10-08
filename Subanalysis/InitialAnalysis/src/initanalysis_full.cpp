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
#include <ConfigManager.h>
#include <NeutralReconstruction.h>

#include <DataAccessWrapper.h>

#include <trilaterationKinFit.h>
#include <signalKinFit.h>

#include "../../Neutrec/inc/trilateration.hpp"

#include "../inc/initialanalysis.hpp"
#include "initialanalysis.hpp"

int InitialAnalysis_full(TChain &chain, Controls::FileType &fileTypeOpt, ErrorHandling::ErrorLogs &logger, KLOE::pm00 &Obj)
{
  ConfigManager &config = ConfigManager::getInstance();

  // TTreeReader reader(&chain);
  // GeneralEventProperties generalProps(reader);
  // ClusterProperties clusterProps(reader);
  // ChargedVertexProperties chVtxProps(reader);
  // BhabhaIP bhabhaProps(reader);

  // --------------- DataAccessWrapper initialization ----------------

  KLOE::DataAccessWrapper dataAccess(chain);

  // Inicjalizacja wrapper'a
  if (!dataAccess.Initialize())
  {
    std::cerr << "ERROR: Failed to initialize DataAccessWrapper" << std::endl;
    return 1;
  }

  std::map<Int_t, Int_t> totEventsPerMctruth = {
      {1, 0}, {2, 0}, {3, 0}, {4, 0}, {5, 0}, {6, 0}, {7, 0}};

  Int_t totEvents = 0;

  BaseKinematics baseKin;

  GeneratedVariables genVarClassifier;
  // Set flag for initial analysis
  Bool_t MonteCarloInitAnalysis = config.getProperty<Bool_t>("flags.initialAnalysisExec.MC");

  // Set flag for covariance matrix type
  std::string covMatrixType = config.getProperty<std::string>("flags.covMatrixType");
  std::string covMatrixName = "momSmearing.covarianceMatrix" + covMatrixType;
  std::string covMatrixNameMC = "momSmearing.covarianceMatrixMC" + covMatrixType;

  std::vector<double> elems = config.getProperty<std::vector<double>>(covMatrixName + ".fElements");
  std::vector<double> elemsMC = config.getProperty<std::vector<double>>(covMatrixNameMC + ".fElements");

  Int_t nRows = config.getProperty<Int_t>(covMatrixName + ".fNrows"),
        nCols = config.getProperty<Int_t>(covMatrixName + ".fNcols");

  Int_t nRowsMC = config.getProperty<Int_t>(covMatrixNameMC + ".fNrows"),
        nColsMC = config.getProperty<Int_t>(covMatrixNameMC + ".fNcols");

  TMatrixT<Double_t>
      covMatrix(nRows, nCols, elems.data()),
      covMatrixMC(nRowsMC, nColsMC, elemsMC.data()),
      covMatrixTot(nRows, nCols);

  // Only difference between data and MC covariance matrices is used

  covMatrixTot = covMatrix - covMatrixMC;

  // --------------------------------------------------------------------------------

  // Which analysis to follow
  std::string hypoCodeStr = config.getProperty<std::string>("flags.analysisCode");
  KLOE::HypothesisCode hypoCode = Obj.StringToHypothesisCode(hypoCodeStr);

  if (hypoCode == KLOE::HypothesisCode::INVALID_VALUE)
    return 1;

  Int_t nPhotons = 0, nPions = 0;

  if (hypoCode == KLOE::HypothesisCode::SIGNAL ||
      hypoCode == KLOE::HypothesisCode::OMEGAPI)
  {
    nPhotons = 4;
    nPions = 2;
  }

  std::ifstream file(cutlimitsName);
  json j = json::parse(file);

  StatisticalCutter cutter(cutlimitsName, 7, hypoCode);

  std::ifstream rootFiles(rootfilesName);
  json filePaths = json::parse(rootFiles);

  std::vector<std::string> baseFilenames = {filePaths["Data"]["filenameBase"],
                                            filePaths["MC"]["filenameBase"][0],
                                            filePaths["MC"]["filenameBase"][1],
                                            filePaths["MC"]["filenameBase"][2]};

  for (Int_t i = 0; i < baseFilenames.size(); i++)
  {
    baseFilenames[i] = baseFilenames[i] + "_" + hypoCodeStr + "_" + covMatrixType;
  }

  std::string
      dirname = (std::string)initialanalysis_dir + (std::string)root_files_dir,
      dated_folder = Obj.CreateDatedFolder(dirname);

  SplitFileWriter writer(baseFilenames[int(fileTypeOpt)], 1.5 * 1024 * 1024 * 1024 * 0.01, false, dated_folder);

  Int_t mcflag = 0, mctruth = 0, NCLMIN = 4; // Assuming NCLMIN is 4, adjust as needed;
  std::vector<Int_t> neuclulist;

  UInt_t mctruth_num[8] = {0, 0, 0, 0, 0, 0, 0, 0}; // Array to hold mctruth values

  // Progress bar
  boost::progress_display show_progress(dataAccess.GetEntries());
  // ---------------------------------------------------

  ErrorHandling::ErrorCodes errorCode = ErrorHandling::ErrorCodes::NO_ERROR;

  // Initialization of Charged part of decay reconstruction class
  // Constructor is below, in the loop
  // boost::optional<KLOE::ChargedVtxRec<>> eventAnalysis;
  KLOE::ChargedVtxRec<> *eventAnalysis = nullptr;
  // -------------------------------------------------------------

  // Initialization of Neutral part of decay reconstruction class
  KLOE::NeutralReconstruction neutRec;

  if (hypoCode == KLOE::HypothesisCode::SIGNAL ||
      hypoCode == KLOE::HypothesisCode::OMEGAPI)
  {
    neutRec.SetNumberOfPhotons(4);
  }
  // -------------------------------------------------------------

  Int_t mode = 1; // Model for pi+pi-

  // GeneralEventPropertiesMC *eventProps;

  // if (MonteCarloInitAnalysis)
  //   eventProps = new GeneralEventPropertiesMC(reader);

  Float_t
      KchrecKSMom = 0,
      KchrecKLMom = 0,
      PmissKS = 0,
      PmissKL = 0,
      EmissKS = 0,
      EmissKL = 0,
      pKTwoBody = 0,
      TrcSum = 0;

  Float_t Emiss = 0., Pmiss = 0., MissMom[3] = {};

  std::vector<KLOE::neutralParticle> photons(nPhotons), pions(nPions);
  KLOE::neutralParticle omega;
  std::vector<KLOE::chargedParticle> chargedPions(2);

  // Cuts application

  if (hypoCode == KLOE::HypothesisCode::FOUR_PI)
  {
    ///////////////////////////////////////////////////////////////////
    cutter.RegisterVariableGetter("InvMassKch", [&]()
                                  { return baseKin.KchrecKS[5]; });
    cutter.RegisterCentralValueGetter("InvMassKch", [&]()
                                      { return mK0; });
    ///////////////////////////////////////////////////////////////////
    cutter.RegisterVariableGetter("InvMassKne", [&]()
                                  { return baseKin.KchrecKL[5]; });
    cutter.RegisterCentralValueGetter("InvMassKne", [&]()
                                      { return mK0; });
    ///////////////////////////////////////////////////////////////////
    cutter.RegisterVariableGetter("TwoBodyMomKS", [&]()
                                  { return KchrecKSMom; });
    cutter.RegisterCentralValueGetter("TwoBodyMomKS", [&]()
                                      { return pKTwoBody; });
    ///////////////////////////////////////////////////////////////////
    cutter.RegisterVariableGetter("TwoBodyMomKL", [&]()
                                  { return KchrecKLMom; });
    cutter.RegisterCentralValueGetter("TwoBodyMomKL", [&]()
                                      { return pKTwoBody; });
    ///////////////////////////////////////////////////////////////////
    cutter.RegisterVariableGetter("MissTotKS", [&]()
                                  { return sqrt(pow(PmissKS, 2) + pow(EmissKS, 2)); });
    ///////////////////////////////////////////////////////////////////
    cutter.RegisterVariableGetter("MissHigherKS", [&]()
                                  { return (pow(EmissKS, 2) - pow(PmissKS, 2)); });
    cutter.RegisterVariableGetter("MissLowerKS", [&]()
                                  { return (pow(EmissKS, 2) - pow(PmissKS, 2)); });
    ///////////////////////////////////////////////////////////////////
    cutter.RegisterVariableGetter("MissTotKL", [&]()
                                  { return sqrt(pow(PmissKL, 2) + pow(EmissKL, 2)); });
    ///////////////////////////////////////////////////////////////////
    cutter.RegisterVariableGetter("MissHigherKL", [&]()
                                  { return (pow(EmissKL, 2) - pow(PmissKL, 2)); });
    cutter.RegisterVariableGetter("MissLowerKL", [&]()
                                  { return (pow(EmissKL, 2) - pow(PmissKL, 2)); });
  }
  else if (hypoCode == KLOE::HypothesisCode::SIGNAL)
  {
    ///////////////////////////////////////////////////////////////////
    cutter.RegisterVariableGetter("InvMassKch", [&]()
                                  { return baseKin.Kchrecnew[5]; });
    cutter.RegisterCentralValueGetter("InvMassKch", [&]()
                                      { return mK0; });

    cutter.RegisterVariableGetter("Chi2Signal", [&]()
                                  { return baseKin.Chi2SignalKinFit; });

    cutter.RegisterVariableGetter("TrcSum", [&]()
                                  { return TrcSum; });

    cutter.RegisterVariableGetter("Qmiss", [&]()
                                  { return sqrt(pow(Pmiss, 2) + pow(Emiss, 2)); });

    cutter.RegisterVariableGetter("InvMassKne", [&]()
                                  { return baseKin.KneTriangle[5]; });
    cutter.RegisterCentralValueGetter("InvMassKne", [&]()
                                      { return mK0; });
  }

  // Initialization of momentum smearing
  // -------------------------------------------------------------
  KLOE::ChargedVtxRec<Float_t, UChar_t> BoostMethodObj;
  // -------------------------------------------------------------

  Bool_t
      good_clus = (Bool_t)properties["variables"]["KinFit"]["Trilateration"]["goodClus"];

  const Short_t
      loopcount = (Short_t)properties["variables"]["KinFit"]["Trilateration"]["loopCount"],
      jmin = (Short_t)properties["variables"]["KinFit"]["Trilateration"]["bunchMin"],
      jmax = (Short_t)properties["variables"]["KinFit"]["Trilateration"]["bunchMax"],
      M = (Short_t)properties["variables"]["KinFit"]["Trilateration"]["numOfConstraints"],
      N_const = (Short_t)properties["variables"]["KinFit"]["Trilateration"]["fixedVars"],
      N_free = (Short_t)properties["variables"]["KinFit"]["Trilateration"]["freeVars"],
      range = Int_t(jmax - jmin) + 1;

  const Short_t
      loopcountSignal = (Short_t)properties["variables"]["KinFit"]["Signal"]["loopCount"],
      MSignal = (Short_t)properties["variables"]["KinFit"]["Signal"]["numOfConstraints"],
      N_constSignal = (Short_t)properties["variables"]["KinFit"]["Signal"]["fixedVars"],
      N_freeSignal = (Short_t)properties["variables"]["KinFit"]["Signal"]["freeVars"];

  const Double_t
      chiSqrStep = (Double_t)properties["variables"]["KinFit"]["Trilateration"]["chiSqrStep"],
      chiSqrStepSignal = (Double_t)properties["variables"]["KinFit"]["Signal"]["chiSqrStep"];

  KLOE::TrilaterationReconstructionKinFit trilatKinFitObj(N_free, N_const, M, loopcount, chiSqrStep, jmin, jmax, logger);
  KLOE::SignalKinFit signalKinFitObj(N_freeSignal, N_constSignal, MSignal, loopcountSignal, chiSqrStepSignal, logger);

  while (dataAccess.Next())
  {
    // Here you would process each entry in the tree.
    // For example, you can read values from the tree and perform calculations.
    // This is a placeholder for your actual analysis logic.
    Bool_t noError = true;
    Bool_t cutCombined = false, passed = false;

    baseKin.Chi2SignalKinFit = 999999.;

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
    baseKin.Kchrecsmeared.clear();

    baseKin.KchboostKS.clear();
    baseKin.KchboostKL.clear();
    baseKin.Kchboostsmeared.clear();
    baseKin.Kchboostnew.clear();

    baseKin.ipKS.clear();
    baseKin.ipKL.clear();
    baseKin.ipnew.clear();

    baseKin.pullsTriKinFit.clear();
    baseKin.pullsSignalFit.clear();

    for (Int_t i = 0; i < 2; i++)
    {
      baseKin.trknew[i].clear();
      baseKin.trkKS[i].clear();
      baseKin.trkKL[i].clear();
      baseKin.trkClosest[i].clear();
      baseKin.trkKLmc[i].clear();
      baseKin.trkKSmc[i].clear();
      baseKin.trksmeared[i].clear();
    }

    baseKin.vtaken.resize(3);
    baseKin.vtakenKS.resize(3);
    baseKin.vtakenKL.resize(3);
    baseKin.vtakenClosest.resize(3);

    baseKin.Kchrecnew.resize(9);
    baseKin.KchrecKS.resize(9);
    baseKin.KchrecKL.resize(9);
    baseKin.KchrecClosest.resize(9);
    baseKin.Kchrecsmeared.resize(9);

    baseKin.KchboostKS.resize(9);
    baseKin.KchboostKL.resize(9);
    baseKin.Kchboostsmeared.resize(9);
    baseKin.Kchboostnew.resize(9);

    baseKin.ipKS.resize(3);
    baseKin.ipKL.resize(3);
    baseKin.ipnew.resize(3);

    baseKin.pullsTriKinFit.resize(0);
    baseKin.pullsSignalFit.resize(0);

    for (Int_t i = 0; i < 2; i++)
    {
      baseKin.trknew[i].resize(4);
      baseKin.trkKS[i].resize(4);
      baseKin.trkKL[i].resize(4);
      baseKin.trkClosest[i].resize(4);
      baseKin.trksmeared[i].resize(4);
    }

    baseKin.ParamSignalFit.resize(N_freeSignal + N_constSignal);
    baseKin.ErrorsSignalFit.resize(N_freeSignal + N_constSignal);

    for (Int_t i = 0; i < 4; i++)
      baseKin.photonFit[i].resize(8);

    baseKin.ipFit.resize(3);
    baseKin.KchrecFit.resize(10);
    baseKin.KchboostFit.resize(10);
    baseKin.KnerecFit.resize(10);
    baseKin.KnereclorFit.resize(10);
    for (Int_t i = 0; i < 2; i++)
      baseKin.trkFit[i].resize(4);

    baseKin.CurvMC.clear();
    baseKin.PhivMC.clear();
    baseKin.CotvMC.clear();

    // Construction of the charged rec class object
    Float_t bhabha_vtx[3] = {dataAccess.GetBx(),
                             dataAccess.GetBy(),
                             dataAccess.GetBz()};

    // ZMIANA: Stwórz zmienne lokalne dla referencji
    Int_t nv_local = dataAccess.GetNV();
    Int_t ntv_local = dataAccess.GetNTV();
    Int_t mode_local = mode;

    // Skopiuj dane iv do lokalnej tablicy (jeśli potrzeba)
    std::vector<Int_t> iv_data = dataAccess.GetIv();

    // Skopiuj dane do lokalnych tablic dla wskaźników
    std::vector<Float_t> curv_data = dataAccess.GetCurv();
    std::vector<Float_t> phiv_data = dataAccess.GetPhiv();
    std::vector<Float_t> cotv_data = dataAccess.GetCotv();
    std::vector<Float_t> xv_data = dataAccess.GetXv();
    std::vector<Float_t> yv_data = dataAccess.GetYv();
    std::vector<Float_t> zv_data = dataAccess.GetZv();

    baseKin.bhabha_vtx[0] = dataAccess.GetBx();
    baseKin.bhabha_vtx[1] = dataAccess.GetBy();
    baseKin.bhabha_vtx[2] = dataAccess.GetBz();

    baseKin.phi_mom[0] = dataAccess.GetBpx();
    baseKin.phi_mom[1] = dataAccess.GetBpy();
    baseKin.phi_mom[2] = dataAccess.GetBpz();
    baseKin.phi_mom[3] = dataAccess.GetBRoots();

    std::vector<std::vector<Float_t>>
        trkMC,
        pgammaMC,
        clusterMC;

    if (MonteCarloInitAnalysis)
    {
      mcflag = 1;

      genVarClassifier.classifyChannel(
          dataAccess.GetNTMC(),
          dataAccess.GetNVtxMC(),
          dataAccess.GetPidMC().data(),
          dataAccess.GetVtxMC().data(),
          dataAccess.GetMother().data(),
          mcflag, // Assuming mcflag is 1 for MC events
          mctruth);

      MctruthCounter(mctruth, mctruth_num);
      // -------------------------------------------------------------------

      genVarClassifier.genVars(dataAccess.GetNTMC(),
                               dataAccess.GetNVtxMC(),
                               dataAccess.GetNClu(),
                               dataAccess.GetPidMC().data(),
                               dataAccess.GetVtxMC().data(),
                               dataAccess.GetMother().data(),
                               dataAccess.GetXvMC().data(),
                               dataAccess.GetYvMC().data(),
                               dataAccess.GetZvMC().data(),
                               dataAccess.GetPxMC().data(),
                               dataAccess.GetPyMC().data(),
                               dataAccess.GetPzMC().data(),
                               mcflag,
                               mctruth,
                               baseKin.ipmc,
                               baseKin.Knemc,
                               baseKin.Kchmc,
                               trkMC,
                               4,
                               pgammaMC,
                               baseKin.CurvMC,
                               baseKin.PhivMC,
                               baseKin.CotvMC,
                               baseKin.goodClusIndex,
                               clusterMC);

      Float_t
          knemcVel = sqrt(pow(baseKin.Knemc[0], 2) + pow(baseKin.Knemc[1], 2) + pow(baseKin.Knemc[2], 2)) / baseKin.Knemc[3],
          knemcPath = sqrt(pow(baseKin.Knemc[6] - baseKin.ipmc[0], 2) + pow(baseKin.Knemc[7] - baseKin.ipmc[1], 2) + pow(baseKin.Knemc[8] - baseKin.ipmc[2], 2)),
          kchmcVel = sqrt(pow(baseKin.Kchmc[0], 2) + pow(baseKin.Kchmc[1], 2) + pow(baseKin.Kchmc[2], 2)) / baseKin.Kchmc[3],
          kchmcPath = sqrt(pow(baseKin.Kchmc[6] - baseKin.ipmc[0], 2) + pow(baseKin.Kchmc[7] - baseKin.ipmc[1], 2) + pow(baseKin.Kchmc[8] - baseKin.ipmc[2], 2));

      baseKin.kaonNeTimeLABMC = knemcPath / (knemcVel);
      baseKin.kaonChTimeLABMC = kchmcPath / (kchmcVel);

      // Go to Kaon CM frame to get the proper time
      TVector3
          kaonNeMomLAB = {-baseKin.Knemc[0] / baseKin.Knemc[3],
                          -baseKin.Knemc[1] / baseKin.Knemc[3],
                          -baseKin.Knemc[2] / baseKin.Knemc[3]},
          kaonChMomLAB = {-baseKin.Kchmc[0] / baseKin.Kchmc[3],
                          -baseKin.Kchmc[1] / baseKin.Kchmc[3],
                          -baseKin.Kchmc[2] / baseKin.Kchmc[3]};

      TLorentzVector
          KaonNe4VecLAB = {baseKin.Knemc[6] - baseKin.ipmc[0], // cm
                           baseKin.Knemc[7] - baseKin.ipmc[1], // cm
                           baseKin.Knemc[8] - baseKin.ipmc[2], // cm
                           baseKin.kaonNeTimeLABMC},           // cm
          KaonNe4VecKaonCM = {0., 0., 0., 0.},
          KaonCh4VecLAB = {baseKin.Kchmc[6] - baseKin.ipmc[0], // cm
                           baseKin.Kchmc[7] - baseKin.ipmc[1], // cm
                           baseKin.Kchmc[8] - baseKin.ipmc[2], // cm
                           baseKin.kaonChTimeLABMC},           // cm
          KaonCh4VecKaonCM = {0., 0., 0., 0.};

      Obj.lorentz_transf(kaonNeMomLAB, KaonNe4VecLAB, KaonNe4VecKaonCM);
      Obj.lorentz_transf(kaonChMomLAB, KaonCh4VecLAB, KaonCh4VecKaonCM);

      baseKin.kaonNeTimeCMMC = KaonNe4VecKaonCM.T() / (cVel * tau_S_nonCPT);
      baseKin.kaonNeTimeLABMC = baseKin.kaonNeTimeLABMC / (cVel * tau_S_nonCPT);

      baseKin.kaonChTimeCMMC = KaonCh4VecKaonCM.T() / (cVel * tau_S_nonCPT);
      baseKin.kaonChTimeLABMC = baseKin.kaonChTimeLABMC / (cVel * tau_S_nonCPT);

      totEventsPerMctruth[mctruth]++;
      totEvents++;

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

    if (mctruth == 1)
    {
      if (hypoCode == KLOE::HypothesisCode::FOUR_PI) // If we look for pipipipi - clusters do not matter
        errorCode = ErrorHandling::ErrorCodes::NO_ERROR;
      else
      {
        errorCode = genVarClassifier.FindNeutralCluster(dataAccess.GetNClu(),
                                                        dataAccess.GetNTCl(),
                                                        dataAccess.GetAssCl().data(),
                                                        NCLMIN,
                                                        logger,
                                                        neuclulist);
      }

      if (errorCode != ErrorHandling::ErrorCodes::NO_ERROR)
      {
        logger.getErrLog(errorCode, "", mctruth);
        noError = false;

        if (mctruth == 1)
        {
          passed = true;
          mctruth = -1;
        }
      }
      else
      {
        // Correction of cluster times based on T0
        baseKin.TclCorr.assign(dataAccess.GetTCl().begin(), dataAccess.GetTCl().end());
        // Obj.CorrectClusterTime(*clusterProps.t0step1, baseKin.TclCorr);

        if (dataAccess.GetNV() >= 1 && dataAccess.GetNTV() >= 1)
        {
          if (eventAnalysis != nullptr)
          {
            delete eventAnalysis; // Usuń poprzedni obiekt jeśli istnieje
          }

          eventAnalysis = new KLOE::ChargedVtxRec<>(nv_local, ntv_local, iv_data.data(), bhabha_vtx, curv_data.data(), phiv_data.data(), cotv_data.data(), xv_data.data(), yv_data.data(), zv_data.data(), mode_local);

          // --------------------------------------------------------------------------------
          // Error codes for different hypotheses
          std::map<KLOE::HypothesisCode, ErrorHandling::ErrorCodes> hypoMap;

          // KMASS HYPOTHESIS - FOR SIGNAL
          hypoMap[KLOE::HypothesisCode::SIGNAL] = eventAnalysis->findKchRec(mcflag, 1, covMatrixTot, baseKin.Kchrecnew, baseKin.trknew[0], baseKin.trknew[1], baseKin.vtaken, logger);

          baseKin.CurvSmeared1 = 1000. / sqrt(pow(baseKin.trknew[0][0], 2) + pow(baseKin.trknew[0][1], 2));
          baseKin.PhivSmeared1 = acos(baseKin.trknew[0][0] / sqrt(pow(baseKin.trknew[0][0], 2) + pow(baseKin.trknew[0][1], 2)));
          baseKin.CotvSmeared1 = baseKin.trknew[0][2] / sqrt(pow(baseKin.trknew[0][0], 2) + pow(baseKin.trknew[0][1], 2));

          baseKin.CurvSmeared2 = 1000. / sqrt(pow(baseKin.trknew[1][0], 2) + pow(baseKin.trknew[1][1], 2));
          baseKin.PhivSmeared2 = acos(baseKin.trknew[1][0] / sqrt(pow(baseKin.trknew[1][0], 2) + pow(baseKin.trknew[1][1], 2)));
          baseKin.CotvSmeared2 = baseKin.trknew[1][2] / sqrt(pow(baseKin.trknew[1][0], 2) + pow(baseKin.trknew[1][1], 2));

          if (Obj.signum(dataAccess.GetCurv()[baseKin.vtaken[1]]) != Obj.signum(baseKin.CurvSmeared1))
          {
            baseKin.CurvSmeared1 = -baseKin.CurvSmeared1;
          }

          if (Obj.signum(dataAccess.GetCurv()[baseKin.vtaken[2]]) != Obj.signum(baseKin.CurvSmeared2))
          {
            baseKin.CurvSmeared2 = -baseKin.CurvSmeared2;
          }

          if (Obj.signum(dataAccess.GetPhiv()[baseKin.vtaken[1]]) != Obj.signum(baseKin.PhivSmeared1))
          {
            baseKin.PhivSmeared1 = -baseKin.PhivSmeared1;
          }

          if (Obj.signum(dataAccess.GetPhiv()[baseKin.vtaken[2]]) != Obj.signum(baseKin.PhivSmeared2))
          {
            baseKin.PhivSmeared2 = -baseKin.PhivSmeared2;
          }

          // Unsmeared versions of vtx variables

          baseKin.Curv1 = dataAccess.GetCurv()[baseKin.vtaken[1]];
          baseKin.Phiv1 = dataAccess.GetPhiv()[baseKin.vtaken[1]];
          baseKin.Cotv1 = dataAccess.GetCotv()[baseKin.vtaken[1]];

          baseKin.Curv2 = dataAccess.GetCurv()[baseKin.vtaken[2]];
          baseKin.Phiv2 = dataAccess.GetPhiv()[baseKin.vtaken[2]];
          baseKin.Cotv2 = dataAccess.GetCotv()[baseKin.vtaken[2]];

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
          {
            logger.getErrLog(errorCode, "", mctruth);
            noError = false;

            if (mctruth == 1)
            {
              passed = true;
              mctruth = -1;
            }
          }
          else
          {
            if (hypoCode == KLOE::HypothesisCode::FOUR_PI)
            {

              if (cutter.PassCut(0) && cutter.PassCut(1))
              {

                Float_t
                    boostPhi[3] = {
                        -dataAccess.GetBpx() / dataAccess.GetBRoots(),
                        -dataAccess.GetBpy() / dataAccess.GetBRoots(),
                        -dataAccess.GetBpz() / dataAccess.GetBRoots()},
                    trkKS_PhiCM[2][4] = {}, KchrecKS_PhiCM[4] = {}, trkKL_PhiCM[2][4], KchrecKL_PhiCM[4] = {};

                pKTwoBody = Obj.TwoBodyDecayMass(mPhi, mK0, mK0);

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

                eventAnalysis->KaonMomFromBoost(baseKin.KchrecKS, baseKin.phi_mom, baseKin.KchboostKS);
                eventAnalysis->KaonMomFromBoost(baseKin.KchrecKL, baseKin.phi_mom, baseKin.KchboostKL);

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
                    MissMomKS[3] = {},
                    MissMomKL[3] = {};

                for (Int_t comp = 0; comp < 3; comp++)
                {
                  MissMomKS[comp] = baseKin.phi_mom[comp] - baseKin.KchboostKS[comp] - baseKin.KchrecKL[comp];
                  MissMomKL[comp] = baseKin.phi_mom[comp] - baseKin.KchboostKL[comp] - baseKin.KchrecKS[comp];
                }

                PmissKS = sqrt(pow(MissMomKS[0], 2) + pow(MissMomKS[1], 2) + pow(MissMomKS[2], 2));
                PmissKL = sqrt(pow(MissMomKL[0], 2) + pow(MissMomKL[1], 2) + pow(MissMomKL[2], 2));

                EmissKS = baseKin.KchboostKS[3] - baseKin.KchrecKS[3];
                EmissKL = baseKin.KchboostKL[3] - baseKin.KchrecKL[3];

                cutter.UpdateStats(mctruth);

                baseKin.cuts.clear();
                baseKin.cuts.resize(cutter.GetCuts().size());

                for (Int_t iter = 0; iter < cutter.GetCuts().size(); iter++)
                  if (cutter.PassCut(iter))
                    baseKin.cuts[iter] = 1;
                  else if (!cutter.PassCut(iter))
                  {
                    baseKin.cuts[iter] = 0;

                    if (mctruth == 7)
                    {
                      passed = true;
                      mctruth = 0;
                    }
                  }
              }
            }
            else if (hypoCode == KLOE::HypothesisCode::SIGNAL)
            {
              // -----------------------------------------------------------------------
              // Boost of the charged part of the decay

              BoostMethodObj.KaonMomFromBoost(baseKin.Kchrecnew, baseKin.phi_mom, baseKin.Kchboostnew);

              Float_t X_line[3] = {baseKin.Kchboostnew[6],
                                   baseKin.Kchboostnew[7],
                                   baseKin.Kchboostnew[8]}, // Vertex laying on the line
                  p[3] = {baseKin.Kchboostnew[0],
                          baseKin.Kchboostnew[1],
                          baseKin.Kchboostnew[2]}, // Direction of the line
                  xB[3] = {baseKin.bhabha_vtx[0],
                           baseKin.bhabha_vtx[1],
                           baseKin.bhabha_vtx[2]}, // Bhabha vertex - laying on the plane
                  plane_perp[3] = {0.,
                                   baseKin.phi_mom[1],
                                   0.}; // Vector perpendicular to the plane from Bhabha momentum

              eventAnalysis->IPBoostCorr(X_line, p, xB, plane_perp, baseKin.ipnew);

              baseKin.ipnew[0] = baseKin.bhabha_vtx[0];
              baseKin.ipnew[1] = baseKin.bhabha_vtx[1];
              // z coordinate of the IP is set to the Bhabha vertex z coordinate if it differs by more than 2 cm
              if (abs(baseKin.ipnew[2] - baseKin.bhabha_vtx[2]) > 2.0)
                baseKin.ipnew[2] = baseKin.bhabha_vtx[2];

              // Go to Kaon CM frame to get the proper time
              TVector3
                  kaonMomLAB = {-baseKin.Kchboostnew[0] / baseKin.Kchboostnew[3],
                                -baseKin.Kchboostnew[1] / baseKin.Kchboostnew[3],
                                -baseKin.Kchboostnew[2] / baseKin.Kchboostnew[3]};

              Double_t
                  KaonPathLAB = sqrt(pow(baseKin.Kchboostnew[6] - baseKin.ipnew[0], 2) +
                                     pow(baseKin.Kchboostnew[7] - baseKin.ipnew[1], 2) +
                                     pow(baseKin.Kchboostnew[8] - baseKin.ipnew[2], 2)),
                  KaonVelocityLAB = kaonMomLAB.Mag();

              baseKin.kaonChTimeLAB = KaonPathLAB / KaonVelocityLAB;

              TLorentzVector
                  Kaon4VecLAB = {baseKin.Kchboostnew[6] - baseKin.ipnew[0], // cm
                                 baseKin.Kchboostnew[7] - baseKin.ipnew[1], // cm
                                 baseKin.Kchboostnew[8] - baseKin.ipnew[2], // cm
                                 baseKin.kaonChTimeLAB},                    // cm
                  Kaon4VecKaonCM = {0., 0., 0., 0.};

              Obj.lorentz_transf(kaonMomLAB, Kaon4VecLAB, Kaon4VecKaonCM);

              baseKin.kaonChTimeCM = Kaon4VecKaonCM.T() / (cVel * tau_S_nonCPT);
              baseKin.kaonChTimeLAB = baseKin.kaonChTimeLAB / (cVel * tau_S_nonCPT);
              // ----------------------------------------------------------------------

              // Application of the cuts
              // 1. Invariant mass of the charged kaon
              // 2. Missing energy Qmiss

              for (Int_t comp = 0; comp < 3; comp++)
              {
                MissMom[comp] = baseKin.Kchboostnew[comp] - baseKin.Kchrecnew[comp];
              }

              Pmiss = sqrt(pow(MissMom[0], 2) + pow(MissMom[1], 2) + pow(MissMom[2], 2));
              Emiss = baseKin.Kchboostnew[3] - baseKin.Kchrecnew[3];

              baseKin.Qmiss = sqrt(pow(Emiss, 2) + pow(Pmiss, 2));

              ErrorHandling::ErrorCodes codeTri = ErrorHandling::ErrorCodes::CHARGED_KAON_MASS_PRE;

              if (cutter.PassCut(2) && cutter.PassCut(3))
              {

                std::vector<Float_t> cluster[5];

                cluster[0] = dataAccess.GetXCl();
                cluster[1] = dataAccess.GetYCl();
                cluster[2] = dataAccess.GetZCl();
                cluster[3] = dataAccess.GetTCl();
                cluster[4] = dataAccess.GetEneCl();

                std::vector<Float_t>
                    bhabha_mom_err = {dataAccess.GetBpxErr(),
                                      dataAccess.GetBpyErr(),
                                      dataAccess.GetBpzErr(),
                                      dataAccess.GetBRootsErr()},
                    bhabha_mom = {dataAccess.GetBpx(),
                                  dataAccess.GetBpy(),
                                  dataAccess.GetBpz(),
                                  dataAccess.GetBRoots()},
                    bhabha_vtx = {dataAccess.GetBx(),
                                  dataAccess.GetBy(),
                                  dataAccess.GetBz()},
                    gamma_mom_final[4];

                gamma_mom_final[0].resize(8);
                gamma_mom_final[1].resize(8);
                gamma_mom_final[2].resize(8);
                gamma_mom_final[3].resize(8);

                // Trilateration Kin Fit + Results

                trilatKinFitObj.SetParameters(cluster, neuclulist, bhabha_mom, bhabha_mom_err, bhabha_vtx);
                errorCode = trilatKinFitObj.Reconstruct();
                trilatKinFitObj.GetResults(baseKin.bunchnum, baseKin.ipTriKinFit, baseKin.g4takenTriKinFit, gamma_mom_final, baseKin.KnetriKinFit, baseKin.neuVtxTriKinFit, baseKin.Chi2TriKinFit, baseKin.pullsTriKinFit);

                baseKin.gammaMomTriKinFit1.assign(gamma_mom_final[0].begin(), gamma_mom_final[0].end());
                baseKin.gammaMomTriKinFit2.assign(gamma_mom_final[1].begin(), gamma_mom_final[1].end());
                baseKin.gammaMomTriKinFit3.assign(gamma_mom_final[2].begin(), gamma_mom_final[2].end());
                baseKin.gammaMomTriKinFit4.assign(gamma_mom_final[3].begin(), gamma_mom_final[3].end());

                if (errorCode != ErrorHandling::ErrorCodes::NO_ERROR)
                {
                  logger.getErrLog(errorCode, "", mctruth);
                  noError = false;

                  if (mctruth == 1)
                  {
                    passed = true;
                    mctruth = -1;
                  }
                }
                else
                {
                  errorCode = TriangleRec(baseKin.g4takenTriKinFit, cluster, neuclulist, bhabha_mom, baseKin.Kchboostnew, baseKin.ipnew, baseKin.KneTriangle, gamma_mom_final, baseKin.minv4gam, baseKin.trcfinal, logger);

                  if (errorCode != ErrorHandling::ErrorCodes::NO_ERROR)
                  {
                    logger.getErrLog(errorCode, "", mctruth);
                    noError = false;

                    TrcSum = -999.;

                    if (mctruth == 1)
                    {
                      passed = true;
                      mctruth = -1;
                    }
                  }
                  else
                  {
                    baseKin.gammaMomTriangle1.assign(gamma_mom_final[0].begin(), gamma_mom_final[0].end());
                    baseKin.gammaMomTriangle2.assign(gamma_mom_final[1].begin(), gamma_mom_final[1].end());
                    baseKin.gammaMomTriangle3.assign(gamma_mom_final[2].begin(), gamma_mom_final[2].end());
                    baseKin.gammaMomTriangle4.assign(gamma_mom_final[3].begin(), gamma_mom_final[3].end());

                    TrcSum = baseKin.trcfinal[0] + baseKin.trcfinal[1] + baseKin.trcfinal[2] + baseKin.trcfinal[3];

                    // Pairing of photons to pions and pion reconstruction

                    for (Int_t i = 0; i < nPhotons; i++)
                    {
                      photons[i].FillFourMom(gamma_mom_final[i][0],
                                             gamma_mom_final[i][1],
                                             gamma_mom_final[i][2],
                                             gamma_mom_final[i][3]);
                    }

                    std::vector<Int_t> bestPairingIndex1;

                    neutRec.PhotonPairingToPi0(photons, bestPairingIndex1);
                    neutRec.Pi0Reconstruction(pions);

                    for (Int_t i = 0; i < 4; i++)
                    {
                      baseKin.pi01[i] = pions[0].fourMom[i];
                      baseKin.pi02[i] = pions[1].fourMom[i];
                    }

                    baseKin.pi01[4] = pions[0].totalMomentum;
                    baseKin.pi01[5] = pions[0].mass;

                    baseKin.pi02[4] = pions[1].totalMomentum;
                    baseKin.pi02[5] = pions[1].mass;
                    ///////////////////////////////////////////////////////////////////

                    // Go to Kaon CM frame to get the proper time

                    std::vector<Float_t> Knereclor = {bhabha_mom[0] - baseKin.Kchboostnew[0],
                                                      bhabha_mom[1] - baseKin.Kchboostnew[1],
                                                      bhabha_mom[2] - baseKin.Kchboostnew[2],
                                                      bhabha_mom[3] - baseKin.Kchboostnew[3]};

                    TVector3
                        kaonMomTriangleLAB = {-Knereclor[0] / Knereclor[3],
                                              -Knereclor[1] / Knereclor[3],
                                              -Knereclor[2] / Knereclor[3]};

                    Double_t
                        KaonPathTriangleLAB = sqrt(pow(baseKin.KneTriangle[6] - baseKin.ipnew[0], 2) +
                                                   pow(baseKin.KneTriangle[7] - baseKin.ipnew[1], 2) +
                                                   pow(baseKin.KneTriangle[8] - baseKin.ipnew[2], 2)),
                        KaonVelocityTriangleLAB = kaonMomTriangleLAB.Mag();

                    baseKin.kaonNeTimeLAB = KaonPathTriangleLAB / KaonVelocityTriangleLAB;

                    TLorentzVector
                        Kaon4VecTriangleLAB = {baseKin.KneTriangle[6] - baseKin.ipnew[0], // cm
                                               baseKin.KneTriangle[7] - baseKin.ipnew[1], // cm
                                               baseKin.KneTriangle[8] - baseKin.ipnew[2], // cm
                                               baseKin.kaonNeTimeLAB},                    // cm						 // cm
                        Kaon4VecKaonTriangleCM = {0., 0., 0., 0.};

                    Obj.lorentz_transf(kaonMomTriangleLAB, Kaon4VecTriangleLAB, Kaon4VecKaonTriangleCM);

                    baseKin.kaonNeTimeCM = Kaon4VecKaonTriangleCM.T() / (cVel * tau_S_nonCPT);
                    baseKin.kaonNeTimeLAB = baseKin.kaonNeTimeLAB / (cVel * tau_S_nonCPT);

                    // Signal Global Kinematic Fit

                    std::vector<Float_t>
                        trackParameters[2],
                        trackParametersErr[2],
                        clusterChosen[4],
                        chargedVtx,
                        chargedVtxErr;

                    trackParameters[0].push_back(baseKin.CurvSmeared1);
                    trackParameters[0].push_back(baseKin.PhivSmeared1);
                    trackParameters[0].push_back(baseKin.CotvSmeared1);
                    trackParameters[1].push_back(baseKin.CurvSmeared2);
                    trackParameters[1].push_back(baseKin.PhivSmeared2);
                    trackParameters[1].push_back(baseKin.CotvSmeared2);

                    trackParametersErr[0].push_back(0.025);
                    trackParametersErr[0].push_back(0.0078);
                    trackParametersErr[0].push_back(0.0076);
                    trackParametersErr[1].push_back(0.020);
                    trackParametersErr[1].push_back(0.022);
                    trackParametersErr[1].push_back(0.0045);

                    for (Int_t k = 0; k < 4; k++)
                    {
                      clusterChosen[k].push_back(dataAccess.GetXCl()[neuclulist[baseKin.g4takenTriKinFit[k]] - 1]);
                      clusterChosen[k].push_back(dataAccess.GetYCl()[neuclulist[baseKin.g4takenTriKinFit[k]] - 1]);
                      clusterChosen[k].push_back(dataAccess.GetZCl()[neuclulist[baseKin.g4takenTriKinFit[k]] - 1]);
                      clusterChosen[k].push_back(dataAccess.GetTCl()[neuclulist[baseKin.g4takenTriKinFit[k]] - 1]);
                      clusterChosen[k].push_back(dataAccess.GetEneCl()[neuclulist[baseKin.g4takenTriKinFit[k]] - 1]);
                    }

                    for (Int_t k = 6; k < 9; k++)
                    {
                      chargedVtx.push_back(baseKin.Kchboostnew[k]);
                    }

                    chargedVtxErr.push_back(0.271);
                    chargedVtxErr.push_back(0.301);
                    chargedVtxErr.push_back(0.212);

                    signalKinFitObj.SetParameters(trackParameters, trackParametersErr, clusterChosen, chargedVtx, chargedVtxErr, bhabha_mom, bhabha_mom_err, bhabha_vtx);
                    errorCode = signalKinFitObj.Reconstruct();
                    signalKinFitObj.GetResults(baseKin.ParamSignal,
                                               baseKin.ErrorsSignal,
                                               baseKin.ParamSignalFit,
                                               baseKin.ErrorsSignalFit,
                                               baseKin.trkFit,
                                               baseKin.KchrecFit,
                                               baseKin.KchboostFit,
                                               baseKin.ipFit,
                                               baseKin.photonFit,
                                               baseKin.KnerecFit,
                                               baseKin.KnereclorFit,
                                               baseKin.Chi2SignalKinFit,
                                               baseKin.pullsSignalFit);

                    // Pairing of photons to pions and pion reconstruction

                    for (Int_t i = 0; i < nPhotons; i++)
                    {
                      photons[i].FillFourMom(baseKin.photonFit[i][0],
                                             baseKin.photonFit[i][1],
                                             baseKin.photonFit[i][2],
                                             baseKin.photonFit[i][3]);
                    }

                    std::vector<Int_t> bestPairingIndex;

                    neutRec.PhotonPairingToPi0(photons, bestPairingIndex);
                    neutRec.Pi0Reconstruction(pions);

                    for (Int_t i = 0; i < 4; i++)
                    {
                      baseKin.pi01Fit[i] = pions[0].fourMom[i];
                      baseKin.pi02Fit[i] = pions[1].fourMom[i];
                    }

                    baseKin.pi01Fit[4] = pions[0].totalMomentum;
                    baseKin.pi01Fit[5] = pions[0].mass;

                    baseKin.pi02Fit[4] = pions[1].totalMomentum;
                    baseKin.pi02Fit[5] = pions[1].mass;
                    ///////////////////////////////////////////////////////////////////

                    if (errorCode != ErrorHandling::ErrorCodes::NO_ERROR)
                    {
                      logger.getErrLog(errorCode, "", mctruth);
                      noError = false;

                      TrcSum = -999.;

                      if (mctruth == 1)
                      {
                        passed = true;
                        mctruth = -1;
                      }
                    }
                    else
                    {

                      baseKin.cuts.clear();
                      baseKin.cuts.resize(cutter.GetCuts().size());

                      for (Int_t iter = 0; iter < cutter.GetCuts().size(); iter++)
                        if (cutter.PassCut(iter))
                          baseKin.cuts[iter] = 1;
                        else if (!cutter.PassCut(iter))
                        {
                          baseKin.cuts[iter] = 0;

                          if (mctruth == 1)
                          {
                            passed = true;
                            mctruth = 0;
                          }
                        }
                    }
                  }
                }
              }

              cutter.UpdateStats(mctruth);
            }

            if ((cutter.PassAllCuts() && noError) || passed)
            {
              errorCode = ErrorHandling::ErrorCodes::NO_ERROR;

              // Clone of the branches of the old tree
              // General properties of the event
              baseKin.nrun = dataAccess.GetNRun();
              baseKin.nev = dataAccess.GetNEv();

              baseKin.necls = dataAccess.GetNECls();
              baseKin.eclfilfo = dataAccess.GetEclFilfo();
              // baseKin.eclfilfoword = *generalProps.eclfilfoword;

              baseKin.eclstream.assign(dataAccess.GetEclStream().begin(), dataAccess.GetEclStream().end());
              // -------------------------------------------------------------------------------------
              // Bhabha interaction point and momentum
              baseKin.Bx = dataAccess.GetBx();
              baseKin.By = dataAccess.GetBy();
              baseKin.Bz = dataAccess.GetBz();
              baseKin.Bsx = dataAccess.GetBxErr();
              baseKin.Bsy = dataAccess.GetByErr();
              baseKin.Bsz = dataAccess.GetBzErr();
              baseKin.Bpx = dataAccess.GetBpx();
              baseKin.Bpy = dataAccess.GetBpy();
              baseKin.Bpz = dataAccess.GetBpz();
              baseKin.Bpxerr = dataAccess.GetBpxErr();
              baseKin.Bpyerr = dataAccess.GetBpyErr();
              baseKin.Bpzerr = dataAccess.GetBpzErr();
              baseKin.Broots = dataAccess.GetBRoots();
              baseKin.BrootsErr = dataAccess.GetBRootsErr();
              // -------------------------------------------------------------------------------------
              // Cluster data
              baseKin.nclu = dataAccess.GetNClu();
              baseKin.ntcl = dataAccess.GetNTCl();
              baseKin.T0step1 = dataAccess.GetT0Step1();
              baseKin.Asscl.assign(dataAccess.GetAssCl().begin(), dataAccess.GetAssCl().end());
              baseKin.Xcl.assign(dataAccess.GetXCl().begin(), dataAccess.GetXCl().end());
              baseKin.Ycl.assign(dataAccess.GetYCl().begin(), dataAccess.GetYCl().end());
              baseKin.Zcl.assign(dataAccess.GetZCl().begin(), dataAccess.GetZCl().end());
              baseKin.Tcl.assign(baseKin.TclCorr.begin(), baseKin.TclCorr.end());
              baseKin.Enecl.assign(dataAccess.GetEneCl().begin(), dataAccess.GetEneCl().end());
              // -------------------------------------------------------------------------------------
              // Charged decay data
              baseKin.nv = dataAccess.GetNV();
              baseKin.ntv = dataAccess.GetNTV();
              baseKin.iv.assign(dataAccess.GetIv().begin(), dataAccess.GetIv().end());
              baseKin.Curv.assign(dataAccess.GetCurv().begin(), dataAccess.GetCurv().end());
              baseKin.Phiv.assign(dataAccess.GetPhiv().begin(), dataAccess.GetPhiv().end());
              baseKin.Cotv.assign(dataAccess.GetCotv().begin(), dataAccess.GetCotv().end());
              baseKin.xv.assign(dataAccess.GetXv().begin(), dataAccess.GetXv().end());
              baseKin.yv.assign(dataAccess.GetYv().begin(), dataAccess.GetYv().end());
              baseKin.zv.assign(dataAccess.GetZv().begin(), dataAccess.GetZv().end());
              // -------------------------------------------------------------------------------------
              // Monte carlo data
              if (MonteCarloInitAnalysis)
              {
                baseKin.ntmc = dataAccess.GetNTMC();
                baseKin.nvtxmc = dataAccess.GetNVtxMC();
                baseKin.vtxmc.assign(dataAccess.GetVtxMC().begin(), dataAccess.GetVtxMC().end());
                baseKin.pidmc.assign(dataAccess.GetPidMC().begin(), dataAccess.GetPidMC().end());
                baseKin.mother.assign(dataAccess.GetMother().begin(), dataAccess.GetMother().end());
                baseKin.xvmc.assign(dataAccess.GetXvMC().begin(), dataAccess.GetXvMC().end());
                baseKin.yvmc.assign(dataAccess.GetYvMC().begin(), dataAccess.GetYvMC().end());
                baseKin.zvmc.assign(dataAccess.GetZvMC().begin(), dataAccess.GetZvMC().end());
                baseKin.pxmc.assign(dataAccess.GetPxMC().begin(), dataAccess.GetPxMC().end());
                baseKin.pymc.assign(dataAccess.GetPyMC().begin(), dataAccess.GetPyMC().end());
                baseKin.pzmc.assign(dataAccess.GetPzMC().begin(), dataAccess.GetPzMC().end());
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
                baseKin.ipmc = {};
                baseKin.Kchmc = {};
                baseKin.Knemc = {};
                baseKin.trkKSmc[0] = {};
                baseKin.trkKSmc[1] = {};
                baseKin.trkKLmc[0] = {};
                baseKin.trkKLmc[1] = {};
              }
              // -------------------------------------------------------------------------------------

              // Int_t zmienne
              std::map<std::string, Int_t> intVars = {
                  {"nrun", baseKin.nrun},                 // Number of run
                  {"nev", baseKin.nev},                   // Number of event
                  {"necls", baseKin.necls},               // Number of ECL words
                  {"Eclfilfo", baseKin.eclfilfo},         // Which filfo was used
                  {"Eclfilfoword", baseKin.eclfilfoword}, // Filfo word
                  {"mcflag", mcflag},                     // If event from MC of Data
                  {"mctruth", mctruth},                   // What event type
                  {"nclu", baseKin.nclu},
                  {"ntcl", baseKin.ntcl},
                  {"nv", baseKin.nv},
                  {"ntv", baseKin.ntv},
                  {"ntmc", baseKin.ntmc},
                  {"nvtxmc", baseKin.nvtxmc},
                  {"bunchnum", baseKin.bunchnum}};

              // Float_t zmienne
              std::map<std::string, Float_t> floatVars = {
                  {"T0step1", baseKin.T0step1},
                  {"Bx", baseKin.Bx},
                  {"By", baseKin.By},
                  {"Bz", baseKin.Bz},
                  {"Bpx", baseKin.Bpx},
                  {"Bpy", baseKin.Bpy},
                  {"Bpz", baseKin.Bpz},
                  {"Broots", baseKin.Broots},
                  {"KaonChTimeLAB", baseKin.kaonChTimeLAB},
                  {"KaonChTimeCM", baseKin.kaonChTimeCM},
                  {"KaonNeTimeLAB", baseKin.kaonNeTimeLAB},
                  {"KaonNeTimeCM", baseKin.kaonNeTimeCM},
                  {"KaonNeTimeLABMC", baseKin.kaonNeTimeLABMC},
                  {"KaonNeTimeCMMC", baseKin.kaonNeTimeCMMC},
                  {"KaonChTimeLABMC", baseKin.kaonChTimeLABMC},
                  {"KaonChTimeCMMC", baseKin.kaonChTimeCMMC},
                  {"Qmiss", baseKin.Qmiss},
                  {"minv4gam", baseKin.minv4gam},
                  {"Chi2TriKinFit", baseKin.Chi2TriKinFit},
                  {"CurvSmeared1", baseKin.CurvSmeared1},
                  {"PhivSmeared1", baseKin.PhivSmeared1},
                  {"CotvSmeared1", baseKin.CotvSmeared1},
                  {"CurvSmeared2", baseKin.CurvSmeared2},
                  {"PhivSmeared2", baseKin.PhivSmeared2},
                  {"CotvSmeared2", baseKin.CotvSmeared2},
                  {"Chi2SignalKinFit", baseKin.Chi2SignalKinFit},
                  {"TrcSum", TrcSum},
                  {"Curv1", baseKin.Curv1},
                  {"Phiv1", baseKin.Phiv1},
                  {"Cotv1", baseKin.Cotv1},
                  {"Curv2", baseKin.Curv2},
                  {"Phiv2", baseKin.Phiv2},
                  {"Cotv2", baseKin.Cotv2}};

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
                  {"cutsApplied", baseKin.cuts},
                  {"g4takenTriKinFit", baseKin.g4takenTriKinFit}};

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
                  {"Kchboost", baseKin.Kchboostnew},
                  {"ip", baseKin.ipnew},
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
                  {"trk2KLmc", baseKin.trkKLmc[1]},
                  {"KnetriKinFit", baseKin.KnetriKinFit},
                  {"ipTriKinFit", baseKin.ipTriKinFit},
                  {"neuVtxTriKinFit", baseKin.neuVtxTriKinFit},
                  {"gammaMomTriKinFit1", baseKin.gammaMomTriKinFit1},
                  {"gammaMomTriKinFit2", baseKin.gammaMomTriKinFit2},
                  {"gammaMomTriKinFit3", baseKin.gammaMomTriKinFit3},
                  {"gammaMomTriKinFit4", baseKin.gammaMomTriKinFit4},
                  {"KneTriangle", baseKin.KneTriangle},
                  {"gammaMomTriangle1", baseKin.gammaMomTriangle1},
                  {"gammaMomTriangle2", baseKin.gammaMomTriangle2},
                  {"gammaMomTriangle3", baseKin.gammaMomTriangle3},
                  {"gammaMomTriangle4", baseKin.gammaMomTriangle4},
                  {"trcfinal", baseKin.trcfinal},
                  {"PhivMC", baseKin.PhivMC},
                  {"CurvMC", baseKin.CurvMC},
                  {"CotvMC", baseKin.CotvMC},
                  {"pullsTriKinFit", baseKin.pullsTriKinFit},
                  {"trk1Fit", baseKin.trkFit[0]},
                  {"trk2Fit", baseKin.trkFit[1]},
                  {"KchrecFit", baseKin.KchrecFit},
                  {"KchboostFit", baseKin.KchboostFit},
                  {"ipFit", baseKin.ipFit},
                  {"photonFit1", baseKin.photonFit[0]},
                  {"photonFit2", baseKin.photonFit[1]},
                  {"photonFit3", baseKin.photonFit[2]},
                  {"photonFit4", baseKin.photonFit[3]},
                  {"KnerecFit", baseKin.KnerecFit},
                  {"KnereclorFit", baseKin.KnereclorFit},
                  {"pullsSignalFit", baseKin.pullsSignalFit},
                  {"ParamSignal", baseKin.ParamSignal},
                  {"ErrorsSignal", baseKin.ErrorsSignal},
                  {"ParamSignalFit", baseKin.ParamSignalFit},
                  {"ErrorsSignalFit", baseKin.ErrorsSignalFit},
                  {"pi01", baseKin.pi01},
                  {"pi02", baseKin.pi02},
                  {"pi01Fit", baseKin.pi01Fit},
                  {"pi02Fit", baseKin.pi02Fit}};

              writer.Fill(intVars, floatVars, intArrays, floatArrays);

              neuclulist.clear(); // Clear the list of neutral clusters for the next event
            }
            else
            {
              errorCode = ErrorHandling::ErrorCodes::CHARGED_KAON_MASS_PRE;
            }
          }
        }
        else
        {
          errorCode = ErrorHandling::ErrorCodes::NO_CHARGED_VTX_OR_TRACKS;
          logger.getErrLog(errorCode, "", mctruth);
          noError = false;
        }

        // ------------------------------------------------------------------
      }
    }

    else
    {
      noError = false;
      passed = false;
    }

    ++show_progress; // Progress of the loading bar
  }

  // Wyniki
  for (size_t i = 0; i < cutter.GetCuts().size(); ++i)
  {
    std::cout << "Cut " << i << ": Eff=" << cutter.GetEfficiency(i) << " +- " << cutter.GetEfficiencyError(i)
              << " Purity=" << cutter.GetPurity(i) << " +- " << cutter.GetPurityError(i)
              << " S/B=" << cutter.GetSignalToBackground(i) << " +- " << cutter.GetSignalToBackgroundError(i) << "\n";
  }

  std::map<ErrorHandling::ErrorCodes, int> physicsErrorCountsPerMctruth[8];

  for (int i = 0; i < 8; ++i)
  {
    physicsErrorCountsPerMctruth[i] = logger.getPhysicsErrorCountsForMctruth(i);
  };

  for (Int_t i = 1; i <= 7; i++)
  {
    Double_t efficiency = 100 * totEventsPerMctruth[i] / (Float_t)totEvents,
             efficiencyError = efficiency / sqrt((Float_t)totEventsPerMctruth[i]);
    std::cout << "Mctruth " << i << " count: " << totEventsPerMctruth[i] << " (" << efficiency << " +- " << efficiencyError << ") %\n";
  }

  logger.printPhysicsErrorStatsPerMctruth(false);

  writer.Close();

  config.setProperty<std::string>("lastUpdate", Obj.getCurrentDate());
  config.setProperty<std::string>("lastScript", "Initial analysis");

  config.saveProperties();

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