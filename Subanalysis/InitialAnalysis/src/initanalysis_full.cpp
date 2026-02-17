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
#include <boost/progress.hpp>
#include <chrono>

#include <event_data.h>
#include <GeneratedVariables.h>
#include <boost/optional.hpp>
#include <SplitFileWriter.h>
#include <charged_mom.h>
#include <StatisticalCutter.h>
#include <ConfigManager.h>
#include <NeutralReconstruction.h>
#include <FileManager.h>

#include <DataAccessWrapper.h>

#include <trilaterationKinFit.h>
#include <signalKinFit.h>
#include <omegaKinFit.h>

#include <AnalysisManager.h>

#include "../../Neutrec/inc/trilateration.hpp"

#include "../inc/initialanalysis.hpp"
#include "initialanalysis.hpp"

int InitialAnalysis_full(TChain &chain, Controls::FileType &fileTypeOpt, ErrorHandling::ErrorLogs &logger, KLOE::pm00 &Obj, bool singleFile)
{
  ConfigManager &config = ConfigManager::getInstance();
  config.loadConfig();

  KLOE::AnalysisConfig &analysisConfig = KLOE::AnalysisConfig::getInstance();
  analysisConfig.Print();
  // --------------- DataAccessWrapper initialization ----------------

  KLOE::DataAccessWrapper dataAccess(chain, logger);

  // Inicjalizacja wrapper'a
  if (!dataAccess.Initialize())
  {
    std::cerr << "ERROR: Failed to initialize DataAccessWrapper" << std::endl;
    return 1;
  }

  //////////////////////////////////////////////////////////////////////

  Int_t totEvents = 0;

  KLOE::BaseKinematics baseKin;

  GeneratedVariables genVarClassifier;
  // Set flag if MC variables should be classified - false for Data
  Bool_t MonteCarloInitAnalysis = analysisConfig.GetActiveHypothesisConfig().modules.classifyMCVariables;

  Int_t mctruthSignal = analysisConfig.GetActiveHypothesisConfig().signal;
  Bool_t SignalOnly = analysisConfig.GetActiveHypothesisConfig().modules.signalOnly;
  Bool_t smearing = analysisConfig.GetActiveHypothesisConfig().modules.momentumSmearing;
  Bool_t trilaterationKinFit = analysisConfig.GetActiveHypothesisConfig().modules.trilaterationKinFit;
  Bool_t signalKinFit = analysisConfig.GetActiveHypothesisConfig().modules.signalKinFit;
  Bool_t omegaKinFit = analysisConfig.GetActiveHypothesisConfig().modules.omegaKinFit;

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
  KLOE::HypothesisCode hypoCode = analysisConfig.activeHypothesis;
  std::string hypoCodeStr = Obj.HypothesisCodeToString(hypoCode);

  if (hypoCode == KLOE::HypothesisCode::INVALID_VALUE)
    return 1;

  Int_t nPhotons = 0, nPions = 0;

  if (hypoCode == KLOE::HypothesisCode::SIGNAL ||
      hypoCode == KLOE::HypothesisCode::OMEGAPI || hypoCode == KLOE::HypothesisCode::FOUR_PI)
  {
    nPhotons = 4;
    nPions = 2;
  }

  std::ifstream file(Paths::cutlimitsName);
  json j = json::parse(file);

  StatisticalCutter cutter(Paths::cutlimitsName, mctruthSignal, hypoCode, logger);

  std::ifstream rootFiles(Paths::rootfilesName);
  json filePaths = json::parse(rootFiles);

  std::vector<std::string> baseFilenames = {filePaths["Data"]["filenameBase"],
                                            filePaths["MC"]["filenameBase"][0],
                                            filePaths["MC"]["filenameBase"][1],
                                            filePaths["MC"]["filenameBase"][2]},
                           baseFilenamesTot = {"",
                                               "",
                                               "",
                                               ""};

  std::string smearingName = "NoSmearing";

  if (analysisConfig.GetActiveHypothesisConfig().modules.momentumSmearing)
  {
    smearingName = covMatrixType;
  }

  if (SignalOnly)
  {
    for (Int_t i = 0; i < baseFilenames.size(); i++)
    {
      baseFilenamesTot[i] = baseFilenames[i] + "_" + hypoCodeStr + "_" + smearingName + "_" + KLOE::channName.at(int(mctruthSignal));
    }
  }
  else
  {
    for (Int_t i = 0; i < baseFilenames.size(); i++)
    {
      baseFilenamesTot[i] = baseFilenames[i] + "_" + hypoCodeStr + "_" + smearingName;
    }
  }

  // Helper function to convert FileType to string
  auto fileTypeToString = [](Controls::FileType fileType) -> std::string
  {
    switch (fileType)
    {
    case Controls::FileType::DATA:
      return "DATA";
    case Controls::FileType::ALL_PHYS:
      return "ALL_PHYS";
    case Controls::FileType::ALL_PHYS2:
      return "ALL_PHYS2";
    case Controls::FileType::ALL_PHYS3:
      return "ALL_PHYS3";
    default:
      return "UNKNOWN";
    }
  };

  std::string fileTypeStr = fileTypeToString(fileTypeOpt);

  std::string
      dirname = (std::string)Paths::initialanalysis_dir + (std::string)Paths::root_files_dir,
      dated_folder = Obj.CreateDatedFolder(dirname),
      log_file_writer_lumi = "";

  if (SignalOnly)
  {
    log_file_writer_lumi = "file_lumi_" + fileTypeStr + "_" + hypoCodeStr + "_" + smearingName + "_" + KLOE::channName.at(int(mctruthSignal)) + ".log";
  }
  else
  {
    log_file_writer_lumi = "file_lumi_" + fileTypeStr + "_" + hypoCodeStr + "_" + smearingName + ".log";
  }

  SplitFileWriter writer(baseFilenamesTot[int(fileTypeOpt)], 1.5 * 1024 * 1024 * 1024 * 0.01, false, dated_folder, log_file_writer_lumi, fileTypeOpt, singleFile);

  KLOE::FileManager fileManager(logger);
  std::string inputLumiLog = "";

  if (SignalOnly)
  {
    inputLumiLog = dated_folder + "/input_luminosity_" + fileTypeStr + "_" + hypoCodeStr + "_" + smearingName + "_" + KLOE::channName.at(int(mctruthSignal)) + ".log";
  }
  else
  {
    inputLumiLog = dated_folder + "/input_luminosity_" + fileTypeStr + "_" + hypoCodeStr + "_" + smearingName + ".log";
  }

  fileManager.LogChainLuminosity(chain, logger, inputLumiLog);

  // Oblicz całkowitą luminozność
  double totalInputLuminosity = 0.0;
  TObjArray *fileElements = chain.GetListOfFiles();
  TIter next(fileElements);
  TChainElement *chEl = nullptr;

  while ((chEl = (TChainElement *)next()))
  {
    Long64_t entries = chEl->GetEntries();
    totalInputLuminosity += KLOE::FileManager::EventsToLuminosity(entries);
  }

  // Współczynnik konwersji [nb^-1 / event]
  const double luminosityPerEvent = 0.000908;

  // W pętli - śledź zmiany plików
  std::string currentInputFile = "";
  Long64_t currentFileEvents = 0;

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

  std::vector<KLOE::neutralParticle> photons(nPhotons), pions(nPions), pionsOmega(nPions);
  KLOE::kaonNeutral Knerec, Knereclor;
  KLOE::neutralParticle omega, omegaFit, omegaFitSignal;
  std::vector<KLOE::chargedParticle> chargedPions(2);

  // Cuts application

  if (hypoCode == KLOE::HypothesisCode::FOUR_PI)
  {
    ///////////////////////////////////////////////////////////////////
    cutter.RegisterVariableGetter("InvMassKch", [&]()
                                  { return baseKin.KchrecKS[5]; });
    cutter.RegisterCentralValueGetter("InvMassKch", [&]()
                                      { return PhysicsConstants::mK0; });
    ///////////////////////////////////////////////////////////////////
    cutter.RegisterVariableGetter("InvMassKne", [&]()
                                  { return baseKin.KchrecKL[5]; });
    cutter.RegisterCentralValueGetter("InvMassKne", [&]()
                                      { return PhysicsConstants::mK0; });
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
                                      { return PhysicsConstants::mK0; });

    cutter.RegisterVariableGetter("InvMassKne", [&]()
                                  { return baseKin.minv4gam; });
    cutter.RegisterCentralValueGetter("InvMassKne", [&]()
                                      { return PhysicsConstants::mK0; });
  }

  // Initialization of momentum smearing
  // -------------------------------------------------------------
  KLOE::ChargedVtxRec<Float_t, UChar_t> BoostMethodObj(logger);
  // -------------------------------------------------------------

  Bool_t
      good_clus = (Bool_t)Utils::properties["variables"]["KinFit"]["Trilateration"]["goodClus"];

  std::map<std::string, Short_t>
      N_free,
      N_const,
      M,
      loopcount;

  std::map<std::string, Double_t>
      chiSqrStep;

  std::vector<std::string> kinFitMethods = {"Trilateration", "Signal", "Omega"};

  for (const auto &method : kinFitMethods)
  {
    loopcount[method] = (Short_t)Utils::properties["variables"]["KinFit"][method]["loopCount"];
    N_free[method] = (Short_t)Utils::properties["variables"]["KinFit"][method]["freeVars"];
    N_const[method] = (Short_t)Utils::properties["variables"]["KinFit"][method]["fixedVars"];
    M[method] = (Short_t)Utils::properties["variables"]["KinFit"][method]["numOfConstraints"];
    chiSqrStep[method] = (Float_t)Utils::properties["variables"]["KinFit"][method]["chiSqrStep"];
  }

  const Short_t
      jmin = (Short_t)Utils::properties["variables"]["KinFit"]["Trilateration"]["bunchMin"],
      jmax = (Short_t)Utils::properties["variables"]["KinFit"]["Trilateration"]["bunchMax"],
      range = Int_t(jmax - jmin) + 1;

  KLOE::TrilaterationReconstructionKinFit trilatKinFitObj(N_free["Trilateration"], N_const["Trilateration"], M["Trilateration"], loopcount["Trilateration"], chiSqrStep["Trilateration"], jmin, jmax, logger);
  KLOE::SignalKinFit signalKinFitObj(N_free["Signal"], N_const["Signal"], M["Signal"], loopcount["Signal"], chiSqrStep["Signal"], logger);
  KLOE::OmegaKinFit omegaKinFitObj(N_free["Omega"], N_const["Omega"], M["Omega"], loopcount["Omega"], chiSqrStep["Omega"], logger);

  // Skopiuj dane iv do lokalnej tablicy (jeśli potrzeba)
  std::vector<Int_t> iv_data;

  // Skopiuj dane do lokalnych tablic dla wskaźników
  std::vector<Float_t> curv_data;
  std::vector<Float_t> phiv_data;
  std::vector<Float_t> cotv_data;
  std::vector<Float_t> xv_data;
  std::vector<Float_t> yv_data;
  std::vector<Float_t> zv_data;

  std::vector<std::vector<Float_t>>
      trkMC,
      pgammaMC,
      clusterMC;

  std::vector<Int_t> ivTmp;
  std::map<Int_t, Int_t> mapTmp;

  // Error codes for different hypotheses
  std::map<KLOE::HypothesisCode, ErrorHandling::ErrorCodes> hypoMap;
  // Kaon times to calculate

  KLOE::KaonProperTimes
      kaonTimesMC,
      kaonTimesRecTriKinFit,
      kaonTimesBoostTriKinFit,
      kaonTimesTriangleRecRec,
      kaonTimesTriangleBoostRec,
      kaonTimesTriangleRecLor,
      kaonTimesTriangleBoostLor,
      kaonTimesSignalKinFit;

  // Auxiliary reconstruction of Kneutral with 6 gammas (for studies)
  KLOE::kaonNeutral KnerecSix;
  std::vector<KLOE::neutralParticle> photonFourMomSix(6);
  Float_t bestError;
  std::vector<Int_t> bestIndicesSix;

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

    baseKin.resize();

    baseKin.ParamSignalFit.resize(N_free["Signal"] + N_const["Signal"]);
    baseKin.ErrorsSignalFit.resize(N_free["Signal"] + N_const["Signal"]);

    baseKin.ParamOmegaFit.resize(N_free["Omega"] + N_const["Omega"]);
    baseKin.ErrorsOmegaFit.resize(N_free["Omega"] + N_const["Omega"]);

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

    TFile *currentFile = chain.GetFile();
    if (currentFile)
    {
      std::string newInputFile = currentFile->GetName();
      if (newInputFile != currentInputFile)
      {
        currentInputFile = newInputFile;

        // Znajdź liczbę zdarzeń w tym pliku
        TObjArray *files = chain.GetListOfFiles();
        TIter it(files);
        TChainElement *el = nullptr;
        while ((el = (TChainElement *)it()))
        {
          if (std::string(el->GetTitle()) == currentInputFile)
          {
            currentFileEvents = el->GetEntries();
            writer.SetCurrentInputFile(currentInputFile, currentFileEvents, luminosityPerEvent);
            break;
          }
        }
      }
    }

    // Construction of the charged rec class object
    Float_t bhabha_vtx[3] = {dataAccess.GetBx(),
                             dataAccess.GetBy(),
                             dataAccess.GetBz()};

    // ZMIANA: Stwórz zmienne lokalne dla referencji
    Int_t nv_local = dataAccess.GetNV();
    Int_t ntv_local = dataAccess.GetNTV();
    Int_t mode_local = mode;

    // Skopiuj dane iv do lokalnej tablicy (jeśli potrzeba)
    iv_data = dataAccess.GetIv();

    // Skopiuj dane do lokalnych tablic dla wskaźników
    curv_data = dataAccess.GetCurv();
    phiv_data = dataAccess.GetPhiv();
    cotv_data = dataAccess.GetCotv();
    xv_data = dataAccess.GetXv();
    yv_data = dataAccess.GetYv();
    zv_data = dataAccess.GetZv();

    baseKin.pxtv = dataAccess.GetPxtv();
    baseKin.pytv = dataAccess.GetPytv();
    baseKin.pztv = dataAccess.GetPztv();
    baseKin.vtxcov[0] = dataAccess.GetVtxCov1();
    baseKin.vtxcov[1] = dataAccess.GetVtxCov2();
    baseKin.vtxcov[2] = dataAccess.GetVtxCov3();
    baseKin.vtxcov[3] = dataAccess.GetVtxCov4();
    baseKin.vtxcov[4] = dataAccess.GetVtxCov5();
    baseKin.vtxcov[5] = dataAccess.GetVtxCov6();

    baseKin.bhabha_vtx[0] = dataAccess.GetBx();
    baseKin.bhabha_vtx[1] = dataAccess.GetBy();
    baseKin.bhabha_vtx[2] = dataAccess.GetBz();

    baseKin.phi_mom[0] = dataAccess.GetBpx();
    baseKin.phi_mom[1] = dataAccess.GetBpy();
    baseKin.phi_mom[2] = dataAccess.GetBpz();
    baseKin.phi_mom[3] = dataAccess.GetBRoots();

    // Transverse momenta of the two charged pions
    Double_t pT1 = 0, pT2 = 0;

    ivTmp = std::vector<Int_t>(dataAccess.GetIv().begin(), dataAccess.GetIv().end());
    mapTmp = Obj.CountRepeatingElements(ivTmp);

    // Sprawdź czy jest przynajmniej jeden wierzchołek z dwoma dołączonymi śladami
    bool hasOne = false;
    for (const auto &pair : mapTmp)
    {
      if (pair.second == 2)
      {
        hasOne = true;
        break;
      }
    }

    // Sprawdź czy są przynajmniej dwa wierzchołki z dwoma dołączonymi śladami
    bool hasTwo = false;
    Int_t countTmp = 0;
    for (const auto &pair : mapTmp)
    {
      if (pair.second == 2)
      {
        if (countTmp == 1)
        {
          hasTwo = true;
          break;
        }
        countTmp++;
      }
    }

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
                               clusterMC,
                               baseKin.muonAlertPlus,
                               baseKin.muonAlertMinus);

      std::vector<Float_t> kaonChMom = {baseKin.Kchmc[0], baseKin.Kchmc[1], baseKin.Kchmc[2], baseKin.Kchmc[3]},
                           kaonChPos = {baseKin.Kchmc[6], baseKin.Kchmc[7], baseKin.Kchmc[8]},
                           kaonNeMom = {baseKin.Knemc[0], baseKin.Knemc[1], baseKin.Knemc[2], baseKin.Knemc[3]},
                           kaonNePos = {baseKin.Knemc[6], baseKin.Knemc[7], baseKin.Knemc[8]},
                           ipPos = {baseKin.ipmc[0], baseKin.ipmc[1], baseKin.ipmc[2]};

      kaonTimesMC = Obj.CalculateKaonProperTimes(kaonChMom, kaonChPos, kaonNeMom, kaonNePos, ipPos);

      KLOE::channEventCount[mctruth]++;

      if (mctruth == 1)
      {
        baseKin.trk1MC[0] = trkMC[0][0];
        baseKin.trk1MC[1] = trkMC[0][1];
        baseKin.trk1MC[2] = trkMC[0][2];
        baseKin.trk1MC[3] = trkMC[0][3];

        baseKin.trk2MC[0] = trkMC[1][0];
        baseKin.trk2MC[1] = trkMC[1][1];
        baseKin.trk2MC[2] = trkMC[1][2];
        baseKin.trk2MC[3] = trkMC[1][3];
      }

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
    else
    {
      // Default values for data events
      mcflag = 0;
      mctruth = 0;
      baseKin.ipmc = {0.0f, 0.0f, 0.0f};
      baseKin.Kchmc = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
      baseKin.Knemc = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
      baseKin.trkMC[0] = {0.0f, 0.0f, 0.0f, 0.0f};
      baseKin.trkMC[1] = {0.0f, 0.0f, 0.0f, 0.0f};
      // // baseKin.goodClusIndex = {0, 0, 0, 0};
      baseKin.trkKSmc[0] = {0.0f, 0.0f, 0.0f, 0.0f};
      baseKin.trkKSmc[1] = {0.0f, 0.0f, 0.0f, 0.0f};
      baseKin.trkKLmc[0] = {0.0f, 0.0f, 0.0f, 0.0f};
      baseKin.trkKLmc[1] = {0.0f, 0.0f, 0.0f, 0.0f};
      baseKin.CurvMC = {0.0f, 0.0f};
      baseKin.PhivMC = {0.0f, 0.0f};
      baseKin.CotvMC = {0.0f, 0.0f};

      kaonTimesMC = KLOE::KaonProperTimes();
    }

    Bool_t badMcTruth = (mctruth != mctruthSignal);

    // If only signal MC is to be analyzed, skip bad mctruth events
    if (SignalOnly && badMcTruth)
    {
      noError = false;
      passed = false;
      ++show_progress;
      continue;
    }

    // --------------------------------------------------------------------------------

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
      LOG_PHYSICS_ERROR(logger, errorCode, mctruth, ErrorHandling::LogFiles::LogType::ERROR);

      noError = false;

      if (mctruth == mctruthSignal)
      {
        passed = true;
        mctruth = -1;

        goto skipEvent;
      }
      else
      {
        neuclulist.clear();
        ++show_progress;
        continue;
      }
    }

    if (!hasOne)
    {
      errorCode = ErrorHandling::ErrorCodes::NO_VTX_WITH_TWO_TRACKS;
      LOG_PHYSICS_ERROR(logger, errorCode, mctruth, ErrorHandling::LogFiles::LogType::ERROR);
      noError = false;

      if (mctruth == mctruthSignal)
      {
        passed = true;
        mctruth = -1;

        goto skipEvent;
      }
      else
      {
        neuclulist.clear();
        ++show_progress;
        continue;
      }
    }

    if (eventAnalysis != nullptr)
    {
      delete eventAnalysis; // Usuń poprzedni obiekt jeśli istnieje
    }

    // Another constructor to be used for KLOE-2 (including helix parameters)
    eventAnalysis = new KLOE::ChargedVtxRec<>(nv_local, ntv_local, iv_data.data(), bhabha_vtx, curv_data.data(), baseKin.pxtv.data(), baseKin.pytv.data(), baseKin.pztv.data(), xv_data.data(), yv_data.data(), zv_data.data(), mode_local, logger);

    // --------------------------------------------------------------------------------
    if (smearing)
    {
      // KMASS HYPOTHESIS WITH MOM SMEARING - FOR SIGNAL
      hypoMap[KLOE::HypothesisCode::SIGNAL] = eventAnalysis->findKchRec(mcflag, 1, covMatrixTot, baseKin.Kchrecnew, baseKin.trknew[0], baseKin.trknew[1], baseKin.vtaken, logger);
    }
    else
    {
      // NO MOMENTUM SMEARING - FOR SIGNAL
      hypoMap[KLOE::HypothesisCode::SIGNAL] = eventAnalysis->findKchRec(mcflag, 0, covMatrixTot, baseKin.Kchrecnew, baseKin.trknew[0], baseKin.trknew[1], baseKin.vtaken, logger);
    }

    pT1 = sqrt(pow(baseKin.trknew[0][0], 2) + pow(baseKin.trknew[0][1], 2)),
    pT2 = sqrt(pow(baseKin.trknew[1][0], 2) + pow(baseKin.trknew[1][1], 2));

    baseKin.CurvSmeared1 = 1000. / pT1;
    baseKin.PhivSmeared1 = atan2(baseKin.trknew[0][1], baseKin.trknew[0][0]);
    baseKin.CotvSmeared1 = baseKin.trknew[0][2] / pT1;

    baseKin.CurvSmeared2 = 1000. / pT2;
    baseKin.PhivSmeared2 = atan2(baseKin.trknew[1][1], baseKin.trknew[1][0]);
    baseKin.CotvSmeared2 = baseKin.trknew[1][2] / pT2;

    if (baseKin.vtaken[1] >= 0 && baseKin.vtaken[2] >= 0)
    {
      if (std::signbit(dataAccess.GetCurv()[baseKin.vtaken[1]]) != std::signbit(baseKin.CurvSmeared1))
      {
        baseKin.CurvSmeared1 = -baseKin.CurvSmeared1;
      }

      if (std::signbit(dataAccess.GetCurv()[baseKin.vtaken[2]]) != std::signbit(baseKin.CurvSmeared2))
      {
        baseKin.CurvSmeared2 = -baseKin.CurvSmeared2;
      }

      // Unsmeared versions of vtx variables
      baseKin.Curv1 = dataAccess.GetCurv()[baseKin.vtaken[1]];
      baseKin.Phiv1 = dataAccess.GetPhiv()[baseKin.vtaken[1]];
      baseKin.Cotv1 = dataAccess.GetCotv()[baseKin.vtaken[1]];

      baseKin.Curv2 = dataAccess.GetCurv()[baseKin.vtaken[2]];
      baseKin.Phiv2 = dataAccess.GetPhiv()[baseKin.vtaken[2]];
      baseKin.Cotv2 = dataAccess.GetCotv()[baseKin.vtaken[2]];
    }
    else
    {
      // If for some reason the track indices are not valid, set smeared variables to 0
      baseKin.CurvSmeared1 = 0;
      baseKin.PhivSmeared1 = 0;
      baseKin.CotvSmeared1 = 0;

      baseKin.CurvSmeared2 = 0;
      baseKin.PhivSmeared2 = 0;
      baseKin.CotvSmeared2 = 0;

      // Also set unsmeared versions to 0
      baseKin.Curv1 = 0;
      baseKin.Phiv1 = 0;
      baseKin.Cotv1 = 0;

      baseKin.Curv2 = 0;
      baseKin.Phiv2 = 0;
      baseKin.Cotv2 = 0;
    }

    // VTX CLOSEST TO BHABHA IP - FOR OMEGAPI
    hypoMap[KLOE::HypothesisCode::OMEGAPI] = eventAnalysis->findKClosestRec(baseKin.KchrecClosest, baseKin.trkClosest[0], baseKin.trkClosest[1], baseKin.vtakenClosest, logger);

    ErrorHandling::ErrorCodes errTmp[2];

    // VTX OF KS - FOR PIPIPIPI
    errTmp[0] = eventAnalysis->findKSLRec(16, -1, baseKin.KchrecKS, baseKin.trkKS[0], baseKin.trkKS[1], baseKin.vtakenKS, logger);

    // --------------------------------------------------------------------------------
    if (hasTwo)
    {
      // VTX OF KL - FOR PIPIPIPI
      errTmp[1] = eventAnalysis->findKSLRec(10, baseKin.vtakenKS[0], baseKin.KchrecKL, baseKin.trkKL[0], baseKin.trkKL[1], baseKin.vtakenKL, logger);
      // --------------------------------------------------------------------------------
    }
    else if (!hasTwo && hypoCode == KLOE::HypothesisCode::FOUR_PI)
      errTmp[1] = ErrorHandling::ErrorCodes::NO_TWO_VTX_WITH_TWO_TRACKS; // Special error if looking for 4pi but only one vertex found
    else
      errTmp[1] = ErrorHandling::ErrorCodes::NO_ERROR; // No error for other hypotheses

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
      LOG_PHYSICS_ERROR(logger, errorCode, mctruth, ErrorHandling::LogFiles::LogType::ERROR);
      noError = false;

      if (mctruth == mctruthSignal)
      {
        passed = true;
        mctruth = -1;

        goto skipEvent;
      }
      else
      {
        neuclulist.clear();
        ++show_progress;
        continue;
      }
    }

    if (!cutter.PassCut(0))
    {
      if (mctruth == mctruthSignal)
      {
        mctruth = 0;
        passed = true;
      }
      else
      {
        neuclulist.clear();
        ++show_progress;
        continue;
      }
    }

    // Calculation of variables specific to 4pi hypothesis
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

        pKTwoBody = Obj.TwoBodyDecayMass(PhysicsConstants::mPhi, PhysicsConstants::mK0, PhysicsConstants::mK0);

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
            plane_perp[3] = {baseKin.phi_mom[0],
                             baseKin.phi_mom[1],
                             0.}; // Vector perpendicular to the plane from Bhabha momentum

        // Corrected IP event by event
        eventAnalysis->IPBoostCorr(X_lineKS, pKL, xB, plane_perp, baseKin.ipKS);
        eventAnalysis->IPBoostCorr(X_lineKL, pKL, xB, plane_perp, baseKin.ipKL);
        // Setting x and y coordinates of the IP to the Bhabha vertex

        baseKin.ipKS[0] = baseKin.bhabha_vtx[0];
        baseKin.ipKS[1] = baseKin.bhabha_vtx[1];
        // z coordinate of the IP is set to the Bhabha vertex z coordinate if it differs by more than 2 cm
        if (abs(baseKin.ipKS[2] - baseKin.bhabha_vtx[2]) > 2.8)
          baseKin.ipKS[2] = baseKin.bhabha_vtx[2];

        baseKin.ipKL[0] = baseKin.bhabha_vtx[0];
        baseKin.ipKL[1] = baseKin.bhabha_vtx[1];
        // z coordinate of the IP is set to the Bhabha vertex z coordinate if it differs by more than 2 cm
        if (abs(baseKin.ipKL[2] - baseKin.bhabha_vtx[2]) > 2.8)
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

            if (mctruth == mctruthSignal)
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
          plane_perp[3] = {baseKin.phi_mom[0],
                           baseKin.phi_mom[1],
                           0.}; // Vector perpendicular to the plane from Bhabha momentum

      eventAnalysis->IPBoostCorr(X_line, p, xB, plane_perp, baseKin.ipnew);

      baseKin.ipnew[0] = baseKin.bhabha_vtx[0];
      baseKin.ipnew[1] = baseKin.bhabha_vtx[1];
      // z coordinate of the IP is set to the Bhabha vertex z coordinate if it differs by more than 2 cm
      if (abs(baseKin.ipnew[2] - baseKin.bhabha_vtx[2]) > 2.8)
        baseKin.ipnew[2] = baseKin.bhabha_vtx[2];

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
        LOG_PHYSICS_ERROR(logger, errorCode, mctruth, ErrorHandling::LogFiles::LogType::ERROR);
        noError = false;

        if (mctruth == mctruthSignal)
        {
          passed = true;
          mctruth = -1;

          goto skipEvent;
        }
        else
        {
          neuclulist.clear();
          ++show_progress;
          continue;
        }
      }
      else
      {
        genVarClassifier.MCvsReconstructedClustersComparator(neuclulist, baseKin.g4takenTriKinFit, dataAccess.GetPNum1(), dataAccess.GetNTMC(), dataAccess.GetMother(), dataAccess.GetVtxMC(), dataAccess.GetPidMC(), dataAccess.GetKine(), dataAccess.GetKinMom(), baseKin.goodClustersTriKinFit);

        errorCode = TriangleRec(baseKin.g4takenTriKinFit, cluster, neuclulist, bhabha_mom, baseKin.Kchboostnew, baseKin.ipnew, baseKin.Knerec, gamma_mom_final, baseKin.minv4gam, baseKin.trcfinal, logger);

        if (!cutter.PassCut(1))
        {
          if (mctruth == mctruthSignal)
          {
            mctruth = 0;
            passed = true;
          }
          else
          {
            neuclulist.clear();
            ++show_progress;
            continue;
          }
        }

        if (errorCode != ErrorHandling::ErrorCodes::NO_ERROR)
        {
          LOG_PHYSICS_ERROR(logger, errorCode, mctruth, ErrorHandling::LogFiles::LogType::ERROR);
          noError = false;

          TrcSum = -999.;

          if (mctruth == mctruthSignal)
          {
            passed = true;
            mctruth = -1;

            goto skipEvent;
          }
          else
          {
            neuclulist.clear();
            ++show_progress;
            continue;
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

          // Looking for pion pairing to omega

          std::vector<Int_t> bestPairingOmegaNeutral, bestPairingOmegaCharged;

          chargedPions[0].FillFourMom(baseKin.trkClosest[0][0],
                                      baseKin.trkClosest[0][1],
                                      baseKin.trkClosest[0][2],
                                      baseKin.trkClosest[0][3]);
          chargedPions[1].FillFourMom(baseKin.trkClosest[1][0],
                                      baseKin.trkClosest[1][1],
                                      baseKin.trkClosest[1][2],
                                      baseKin.trkClosest[1][3]);

          neutRec.PhotonPairingToPi0WithOmega(photons, chargedPions, bestPairingOmegaNeutral, bestPairingOmegaCharged, omega, pionsOmega);

          omega.total[6] = baseKin.KchrecClosest[6];
          omega.total[7] = baseKin.KchrecClosest[7];
          omega.total[8] = baseKin.KchrecClosest[8]; // Use kaon decay z vertex as omega z vertex - better resolution

          if (neuclulist.size() >= 6)
          {
            // 4 clusters chosen with the trilateration + 2 remaining clusters
            // To keep the statistical independence of the samples
            neutRec.ReconstructSixGammaVertexWithFourTaken(cluster, neuclulist, baseKin.g4takenTriKinFit, bestIndicesSix, bestError, KnerecSix, photonFourMomSix);
          }
          else
          {
            bestError = 999999.;
            KnerecSix.total = {0., 0., 0., 0., 0., 0., 0., 0., 0.};
            photonFourMomSix = std::vector<KLOE::neutralParticle>(6, KLOE::neutralParticle());
          }
          ///////////////////////////////////////////////////////////////////

          // Signal Global Kinematic Fit

          std::vector<Float_t>
              trackParameters[2],
              trackParametersErr[2],
              clusterChosen[4],
              chargedVtx,
              chargedVtxErr,
              neuVtx,
              neuVtxErr,
              bhabhaVtxErr;

          trackParameters[0].push_back(baseKin.trknew[0][0]);
          trackParameters[0].push_back(baseKin.trknew[0][1]);
          trackParameters[0].push_back(baseKin.trknew[0][2]);
          trackParameters[1].push_back(baseKin.trknew[1][0]);
          trackParameters[1].push_back(baseKin.trknew[1][1]);
          trackParameters[1].push_back(baseKin.trknew[1][2]);

          trackParametersErr[0].push_back(pow(1.5, 2) / 2.0);
          trackParametersErr[0].push_back(pow(1.5, 2) / 2.0);
          trackParametersErr[0].push_back(pow(1.8, 2) / 2.0);
          trackParametersErr[1].push_back(pow(1.5, 2) / 2.0);
          trackParametersErr[1].push_back(pow(1.5, 2) / 2.0);
          trackParametersErr[1].push_back(pow(1.8, 2) / 2.0);

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

          chargedVtxErr.push_back(sqrt(baseKin.vtxcov[0][baseKin.vtaken[0]]));
          chargedVtxErr.push_back(sqrt(baseKin.vtxcov[3][baseKin.vtaken[0]]));
          chargedVtxErr.push_back(sqrt(baseKin.vtxcov[5][baseKin.vtaken[0]]));

          for (Int_t k = 6; k < 9; k++)
          {
            neuVtx.push_back(baseKin.Knerec[k]);
          }

          neuVtxErr.push_back(0.523);
          neuVtxErr.push_back(0.520);
          neuVtxErr.push_back(1.334);

          // bhabha_vtx errors from the data access
          bhabhaVtxErr.push_back(sqrt(pow(dataAccess.GetBxErr(), 2) + pow(dataAccess.GetBlumx(), 2)));
          bhabhaVtxErr.push_back(dataAccess.GetByErr());
          bhabhaVtxErr.push_back(sqrt(pow(dataAccess.GetBzErr(), 2) + pow(dataAccess.GetBlumz(), 2)));

          if (signalKinFit)
          {
            signalKinFitObj.SetParameters(trackParameters, trackParametersErr, clusterChosen, chargedVtx, chargedVtxErr, bhabha_mom, bhabha_mom_err, neuVtx, neuVtxErr, bhabha_vtx, bhabhaVtxErr);
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

            for (Int_t i = 0; i < nPhotons; i++)
            {
              photons[i].FillFourMom(baseKin.photonFit[i][0],
                                     baseKin.photonFit[i][1],
                                     baseKin.photonFit[i][2],
                                     baseKin.photonFit[i][3]);

              photons[i].fourPos[0] = baseKin.photonFit[i][4];
              photons[i].fourPos[1] = baseKin.photonFit[i][5];
              photons[i].fourPos[2] = baseKin.photonFit[i][6];
              photons[i].fourPos[3] = baseKin.photonFit[i][7];
            }

            std::vector<Int_t> bestPairingIndex, bestPairingIndexNeutral, bestPairingIndexCharged;

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
          }

          // Go to Kaon CM frame to get the proper time
          baseKin.Knereclor[0] = bhabha_mom[0] - baseKin.Kchboostnew[0];
          baseKin.Knereclor[1] = bhabha_mom[1] - baseKin.Kchboostnew[1];
          baseKin.Knereclor[2] = bhabha_mom[2] - baseKin.Kchboostnew[2];
          baseKin.Knereclor[3] = bhabha_mom[3] - baseKin.Kchboostnew[3];
          baseKin.Knereclor[4] = sqrt(pow(baseKin.Knereclor[0], 2) + pow(baseKin.Knereclor[1], 2) + pow(baseKin.Knereclor[2], 2));
          baseKin.Knereclor[5] = PhysicsConstants::mK0;
          baseKin.Knereclor[6] = baseKin.Knerec[6];
          baseKin.Knereclor[7] = baseKin.Knerec[7];
          baseKin.Knereclor[8] = baseKin.Knerec[8];

          std::vector<Float_t>
              kaonChMomRec = {baseKin.Kchrecnew[0], baseKin.Kchrecnew[1], baseKin.Kchrecnew[2], baseKin.Kchrecnew[3]},
              kaonChMomBoost = {baseKin.Kchboostnew[0], baseKin.Kchboostnew[1], baseKin.Kchboostnew[2], baseKin.Kchboostnew[3]},
              kaonChMomSignalKinFit = {baseKin.KchboostFit[0], baseKin.KchboostFit[1], baseKin.KchboostFit[2], baseKin.KchboostFit[3]},
              kaonChPos = {baseKin.Kchboostnew[6], baseKin.Kchboostnew[7], baseKin.Kchboostnew[8]},
              kaonChPosSignalKinFit = {baseKin.KchboostFit[6], baseKin.KchboostFit[7], baseKin.KchboostFit[8]},
              kaonNeMomRec = {baseKin.Knerec[0], baseKin.Knerec[1], baseKin.Knerec[2], baseKin.Knerec[3]},
              kaonNeMomLor = {baseKin.Knereclor[0], baseKin.Knereclor[1], baseKin.Knereclor[2], baseKin.Knereclor[3]},
              kaonNeMomSignalKinFit = {baseKin.KnerecFit[0], baseKin.KnerecFit[1], baseKin.KnerecFit[2], baseKin.KnerecFit[3]},
              kaonNeMomTriKinFit = {baseKin.KnetriKinFit[0], baseKin.KnetriKinFit[1], baseKin.KnetriKinFit[2], baseKin.KnetriKinFit[3]},
              kaonNePos = {baseKin.Knerec[6], baseKin.Knerec[7], baseKin.Knerec[8]},
              kaonNePosSignalKinFit = {baseKin.KnerecFit[6], baseKin.KnerecFit[7], baseKin.KnerecFit[8]},
              kaonNePosTriKinFit = {baseKin.KnetriKinFit[6], baseKin.KnetriKinFit[7], baseKin.KnetriKinFit[8]},
              ipPos = {baseKin.ipnew[0], baseKin.ipnew[1], baseKin.ipnew[2]},
              ipPosSignalKinFit = {baseKin.ipFit[0], baseKin.ipFit[1], baseKin.ipFit[2]},
              ipPosTriKinFit = {baseKin.ipTriKinFit[0], baseKin.ipTriKinFit[1], baseKin.ipTriKinFit[2]};

          kaonTimesTriangleRecRec = Obj.CalculateKaonProperTimes(kaonChMomRec, kaonChPos, kaonNeMomRec, kaonNePos, ipPos);

          kaonTimesTriangleBoostRec = Obj.CalculateKaonProperTimes(kaonChMomBoost, kaonChPos, kaonNeMomRec, kaonNePos, ipPos);

          kaonTimesTriangleRecLor = Obj.CalculateKaonProperTimes(kaonChMomRec, kaonChPos, kaonNeMomLor, kaonNePos, ipPos);

          kaonTimesTriangleBoostLor = Obj.CalculateKaonProperTimes(kaonChMomBoost, kaonChPos, kaonNeMomLor, kaonNePos, ipPos);

          kaonTimesRecTriKinFit = Obj.CalculateKaonProperTimes(kaonChMomRec, kaonChPos, kaonNeMomTriKinFit, kaonNePosTriKinFit, ipPosTriKinFit);

          kaonTimesBoostTriKinFit = Obj.CalculateKaonProperTimes(kaonChMomBoost, kaonChPos, kaonNeMomTriKinFit, kaonNePosTriKinFit, ipPosTriKinFit);

          kaonTimesSignalKinFit = Obj.CalculateKaonProperTimes(kaonChMomSignalKinFit, kaonChPosSignalKinFit, kaonNeMomSignalKinFit, kaonNePosSignalKinFit, ipPosSignalKinFit);

          // Additional Omega-Pi0 fit for better bkg rejection
          if (omegaKinFit)
          {
            trackParameters[0].clear();
            trackParameters[1].clear();
            trackParametersErr[0].clear();
            trackParametersErr[1].clear();

            trackParameters[0].resize(0);
            trackParameters[1].resize(0);
            trackParametersErr[0].resize(0);
            trackParametersErr[1].resize(0);

            trackParameters[0].push_back(baseKin.trkClosest[0][0]);
            trackParameters[0].push_back(baseKin.trkClosest[0][1]);
            trackParameters[0].push_back(baseKin.trkClosest[0][2]);
            trackParameters[1].push_back(baseKin.trkClosest[1][0]);
            trackParameters[1].push_back(baseKin.trkClosest[1][1]);
            trackParameters[1].push_back(baseKin.trkClosest[1][2]);

            trackParametersErr[0].push_back(pow(1.5, 2) / 2.0);
            trackParametersErr[0].push_back(pow(1.5, 2) / 2.0);
            trackParametersErr[0].push_back(pow(1.8, 2) / 2.0);
            trackParametersErr[1].push_back(pow(1.5, 2) / 2.0);
            trackParametersErr[1].push_back(pow(1.5, 2) / 2.0);
            trackParametersErr[1].push_back(pow(1.8, 2) / 2.0);

            std::vector<Float_t> omegaVtx = {baseKin.KchrecClosest[6],
                                             baseKin.KchrecClosest[7],
                                             baseKin.KchrecClosest[8]},
                                 omegaVtxErr = {sqrt(baseKin.vtxcov[0][baseKin.vtakenClosest[0]]),
                                                sqrt(baseKin.vtxcov[3][baseKin.vtakenClosest[0]]),
                                                sqrt(baseKin.vtxcov[5][baseKin.vtakenClosest[0]])};

            omegaKinFitObj.SetParameters(trackParameters, trackParametersErr, clusterChosen, bhabha_mom, bhabha_mom_err, chargedVtx, chargedVtxErr, omegaVtx, omegaVtxErr);
            errorCode = omegaKinFitObj.Reconstruct();
            omegaKinFitObj.GetResults(baseKin.ParamOmega,
                                      baseKin.ErrorsOmega,
                                      baseKin.ParamOmegaFit,
                                      baseKin.ErrorsOmegaFit,
                                      baseKin.trkOmegaFit,
                                      baseKin.ipOmegaFit,
                                      baseKin.photonOmegaFit,
                                      baseKin.omegaFit,
                                      baseKin.pi0OmegaFit,
                                      baseKin.phiOmegaFit,
                                      baseKin.Chi2OmegaKinFit,
                                      baseKin.pullsOmegaFit);
          }

          if (errorCode != ErrorHandling::ErrorCodes::NO_ERROR)
          {
            LOG_PHYSICS_ERROR(logger, errorCode, mctruth, ErrorHandling::LogFiles::LogType::ERROR);
            noError = false;

            TrcSum = -999.;

            if (mctruth == mctruthSignal)
            {
              passed = true;
              mctruth = -1;
            }
            {
              neuclulist.clear();
              ++show_progress;
              continue;
            }
          }
          else
          {
            for (Int_t iter = 0; iter < cutter.GetCuts().size(); iter++)
              if (!cutter.PassCut(iter))
              {
                switch (iter)
                {
                case 0:
                {
                  baseKin.cut = Int_t(ErrorHandling::ErrorCodes::CUT_CHI2_SIGNAL);
                  break;
                }
                case 1:
                {
                  baseKin.cut = Int_t(ErrorHandling::ErrorCodes::CUT_TRCV);
                  break;
                }
                case 2:
                {
                  baseKin.cut = Int_t(ErrorHandling::ErrorCodes::CUT_INV_MASS_KCH);
                  break;
                }
                case 3:
                {
                  baseKin.cut = Int_t(ErrorHandling::ErrorCodes::CUT_INV_MASS_KNE);
                  break;
                }
                case 4:
                {
                  baseKin.cut = Int_t(ErrorHandling::ErrorCodes::CUT_QMISS);
                  break;
                }
                }

                if (mctruth == 1)
                {
                  passed = true;
                  mctruth = 0;
                }

                break;
              }
          }
        }
      }

      cutter.UpdateStats(mctruth);
    }

  // All mctruth = 1 with errors appear here
  skipEvent:

    if ((cutter.PassAllCuts() && noError) || passed)
    {

      // Assign Int_t value for errorCode
      baseKin.errorCode = static_cast<Int_t>(errorCode);

      // Clone of the branches of the old tree
      // General Utils::properties of the event
      baseKin.nrun = dataAccess.GetNRun();
      baseKin.nev = dataAccess.GetNEv();

      baseKin.necls = dataAccess.GetNECls();
      baseKin.eclfilfo = dataAccess.GetEclFilfo();

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
      baseKin.Tcl.assign(dataAccess.GetTCl().begin(), dataAccess.GetTCl().end());
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
        baseKin.vtxmc.assign(20, 0.0f);
        baseKin.pidmc.assign(20, 0);
        baseKin.mother.assign(20, 0);
        baseKin.xvmc.assign(20, 0.0f);
        baseKin.yvmc.assign(20, 0.0f);
        baseKin.zvmc.assign(20, 0.0f);
        baseKin.pxmc.assign(20, 0.0f);
        baseKin.pymc.assign(20, 0.0f);
        baseKin.pzmc.assign(20, 0.0f);
        baseKin.ipmc.assign(3, 0.0f);
        baseKin.Kchmc.assign(9, 0.0f);
        baseKin.Knemc.assign(9, 0.0f);
        baseKin.trkKSmc[0].assign(4, 0.0f);
        baseKin.trkKSmc[1].assign(4, 0.0f);
        baseKin.trkKLmc[0].assign(4, 0.0f);
        baseKin.trkKLmc[1].assign(4, 0.0f);
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
          {"bunchnum", baseKin.bunchnum},
          {"errorcode", baseKin.errorCode},
          {"goodClustersTriKinFitSize", baseKin.goodClustersTriKinFit.size()},
          {"cutApplied", baseKin.cut},
          {"muonAlertPlus", baseKin.muonAlertPlus},
          {"muonAlertMinus", baseKin.muonAlertMinus}};

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
          {"KaonChTimeLABBoostLor", kaonTimesTriangleBoostLor.kaon1TimeLAB},
          {"KaonChTimeCMBoostLor", kaonTimesTriangleBoostLor.kaon1TimeCM},
          {"KaonNeTimeLABBoostLor", kaonTimesTriangleBoostLor.kaon2TimeLAB},
          {"KaonNeTimeCMBoostLor", kaonTimesTriangleBoostLor.kaon2TimeCM},
          {"KaonChTimeLABRecLor", kaonTimesTriangleRecLor.kaon1TimeLAB},
          {"KaonChTimeCMRecLor", kaonTimesTriangleRecLor.kaon1TimeCM},
          {"KaonNeTimeLABRecLor", kaonTimesTriangleRecLor.kaon2TimeLAB},
          {"KaonNeTimeCMRecLor", kaonTimesTriangleRecLor.kaon2TimeCM},
          {"KaonChTimeLABBoostRec", kaonTimesTriangleBoostRec.kaon1TimeLAB},
          {"KaonChTimeCMBoostRec", kaonTimesTriangleBoostRec.kaon1TimeCM},
          {"KaonNeTimeLABBoostRec", kaonTimesTriangleBoostRec.kaon2TimeLAB},
          {"KaonNeTimeCMBoostRec", kaonTimesTriangleBoostRec.kaon2TimeCM},
          {"KaonChTimeLABRecRec", kaonTimesTriangleRecRec.kaon1TimeLAB},
          {"KaonChTimeCMRecRec", kaonTimesTriangleRecRec.kaon1TimeCM},
          {"KaonNeTimeLABRecRec", kaonTimesTriangleRecRec.kaon2TimeLAB},
          {"KaonNeTimeCMRecRec", kaonTimesTriangleRecRec.kaon2TimeCM},
          {"KaonChTimeLABBoostTriFit", kaonTimesBoostTriKinFit.kaon1TimeLAB},
          {"KaonChTimeCMBoostTriFit", kaonTimesBoostTriKinFit.kaon1TimeCM},
          {"KaonNeTimeLABBoostTriFit", kaonTimesBoostTriKinFit.kaon2TimeLAB},
          {"KaonNeTimeCMBoostTriFit", kaonTimesBoostTriKinFit.kaon2TimeCM},
          {"KaonChTimeLABRecTriFit", kaonTimesRecTriKinFit.kaon1TimeLAB},
          {"KaonChTimeCMRecTriFit", kaonTimesRecTriKinFit.kaon1TimeCM},
          {"KaonNeTimeLABRecTriFit", kaonTimesRecTriKinFit.kaon2TimeLAB},
          {"KaonNeTimeCMRecTriFit", kaonTimesRecTriKinFit.kaon2TimeCM},
          {"KaonNeTimeLABMC", kaonTimesMC.kaon2TimeLAB},
          {"KaonNeTimeCMMC", kaonTimesMC.kaon2TimeCM},
          {"KaonChTimeLABMC", kaonTimesMC.kaon1TimeLAB},
          {"KaonChTimeCMMC", kaonTimesMC.kaon1TimeCM},
          {"KaonNeTimeLABSignalFit", kaonTimesSignalKinFit.kaon2TimeLAB},
          {"KaonNeTimeCMSignalFit", kaonTimesSignalKinFit.kaon2TimeCM},
          {"KaonChTimeLABSignalFit", kaonTimesSignalKinFit.kaon1TimeLAB},
          {"KaonChTimeCMSignalFit", kaonTimesSignalKinFit.kaon1TimeCM},
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
          {"Cotv2", baseKin.Cotv2},
          {"Chi2OmegaKinFit", baseKin.Chi2OmegaKinFit},
          {"bestErrorSixGamma", bestError}};

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
          {"g4takenTriKinFit", baseKin.g4takenTriKinFit},
          {"goodClustersTriKinFit", baseKin.goodClustersTriKinFit}};

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
          {"Knerec", baseKin.Knerec},
          {"Knereclor", baseKin.Knereclor},
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
          {"pi02Fit", baseKin.pi02Fit},
          {"ParamOmega", baseKin.ParamOmega},
          {"ErrorsOmega", baseKin.ErrorsOmega},
          {"ParamOmegaFit", baseKin.ParamOmegaFit},
          {"ErrorsOmegaFit", baseKin.ErrorsOmegaFit},
          {"trkOmegaFit1", baseKin.trkOmegaFit[0]},
          {"trkOmegaFit2", baseKin.trkOmegaFit[1]},
          {"omegaFit", baseKin.omegaFit},
          {"pi0OmegaFit1", baseKin.pi0OmegaFit[0]},
          {"pi0OmegaFit2", baseKin.pi0OmegaFit[1]},
          {"phiOmegaFit", baseKin.phiOmegaFit},
          {"omega", omega.total},
          {"pi0Omega1", pionsOmega[0].total},
          {"pi0Omega2", pionsOmega[1].total},
          {"photonOmegaFit1", baseKin.photonOmegaFit[0]},
          {"photonOmegaFit2", baseKin.photonOmegaFit[1]},
          {"photonOmegaFit3", baseKin.photonOmegaFit[2]},
          {"photonOmegaFit4", baseKin.photonOmegaFit[3]},
          {"ipOmegaFit", baseKin.ipOmegaFit},
          {"trk1MC", baseKin.trk1MC},
          {"trk2MC", baseKin.trk2MC},
          {"KnerecSix", KnerecSix.total},
          {"photonSix1", photonFourMomSix[0].total},
          {"photonSix2", photonFourMomSix[1].total},
          {"photonSix3", photonFourMomSix[2].total},
          {"photonSix4", photonFourMomSix[3].total},
          {"photonSix5", photonFourMomSix[4].total},
          {"photonSix6", photonFourMomSix[5].total}};

      writer.Fill(intVars, floatVars, intArrays, floatArrays);
    }
    else
    {
      errorCode = ErrorHandling::ErrorCodes::CHARGED_KAON_MASS_PRE;
    }
    // ------------------------------------------------------------------

    neuclulist.clear(); // Clear the list of neutral clusters for the next event

    ++show_progress; // Progress of the loading bar
  }

  // End of event loop
  totEvents = KLOE::TotalCountMC();

  for (Int_t i = 1; i <= 7; i++)
    std::cout << "Mctruth " << i << " count: " << KLOE::channEventCount[i] << std::endl;

  // Wyniki
  cutter.CutSummary();

  std::map<ErrorHandling::ErrorCodes, int> physicsErrorCountsPerMctruth[8];

  for (int i = 0; i < 8; ++i)
  {
    physicsErrorCountsPerMctruth[i] = logger.getPhysicsErrorCountsForMctruth(i);
  };

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