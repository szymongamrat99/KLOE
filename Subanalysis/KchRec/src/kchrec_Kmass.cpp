// Author: Szymon Gamrat
// Date of last update: 03.02.2025

#include <TGraph.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TFitResultPtr.h>
#include <TFitResult.h>
#include <TLatex.h>
#include <triple_gaus.h>
#include <TLegend.h>
#include <THStack.h>
#include <TTreeReader.h>

#include <boost/optional.hpp>
#include <boost/progress.hpp>
#include <SplitFileWriter.h>
#include <event_data.h>
#include <ConfigManager.h>
#include <StatisticalCutter.h>

#include "../inc/kchrec.hpp"

/**
 * @brief Reconstructs the kinematics of charged kaon decays and performs analysis for the KchRec control sample.
 *
 * This function processes entries from a ROOT TChain, reconstructs the kinematics of charged kaon decays,
 * calculates momenta using two-body decay kinematics, applies boost corrections, and fills histograms for
 * further analysis. The results are written to output files for later inspection.
 *
 * @param chain     Reference to the ROOT TChain containing the event data.
 * @param dataType  Reference to the Controls::DataType object specifying the data type (MC/data).
 * @param logger    Reference to the ErrorHandling::ErrorLogs object for logging errors and messages.
 * @param Obj       Reference to the KLOE::pm00 object providing utility and transformation methods.
 * @return int      Returns 0 on successful completion.
 *
 * @details
 * - Sets up branch addresses for all relevant event variables.
 * - Initializes histograms and output file writers.
 * - Loops over all entries in the TChain:
 *   - Reconstructs charged kaon decay vertices and momenta.
 *   - Calculates pion momenta in the kaon rest frame and boosts them to the lab frame.
 *   - Rotates momenta to align with the kaon momentum direction.
 *   - Compares reconstructed and theoretical momenta, filling histograms for differences.
 *   - Writes results to output files.
 * - At the end, saves all histograms as images and closes output files.
 */
int kchrec_Kmass(TChain &chain, Controls::DataType &dataType, ErrorHandling::ErrorLogs &logger, KLOE::pm00 &Obj)
{
  ConfigManager &config = ConfigManager::getInstance();

  // Structure from const.h to ease navigation
  KLOE::BaseKinematics baseKin;

  TTreeReader reader(&chain);

  BhabhaIP bhabhaProps(reader);

  chain.SetBranchAddress("nev", &baseKin.nev);
  chain.SetBranchAddress("nrun", &baseKin.nrun);

  // Bhabha variables - avg per run
  chain.SetBranchAddress("Bx", &baseKin.bhabha_vtx[0]);
  chain.SetBranchAddress("By", &baseKin.bhabha_vtx[1]);
  chain.SetBranchAddress("Bz", &baseKin.bhabha_vtx[2]);

  chain.SetBranchAddress("Bpx", &baseKin.phi_mom[0]);
  chain.SetBranchAddress("Bpy", &baseKin.phi_mom[1]);
  chain.SetBranchAddress("Bpz", &baseKin.phi_mom[2]);
  chain.SetBranchAddress("Broots", &baseKin.phi_mom[3]);

  std::vector<Float_t>
      *ipKS(&baseKin.ipKS),
      *ipKL(&baseKin.ipKL),
      *ipmc(&baseKin.ipmc),
      *trkKS[2],
      *trkKL[2],
      *trkKLmc[2],
      *trkKSmc[2],
      *KchrecKS(&baseKin.KchrecKS),
      *KchrecKL(&baseKin.KchrecKL),
      *KchboostKS(&baseKin.KchboostKS),
      *KchboostKL(&baseKin.KchboostKL),
      *Knemc(&baseKin.Knemc),
      *Kchmc(&baseKin.Kchmc);

  trkKS[0] = &baseKin.trkKS[0];
  trkKS[1] = &baseKin.trkKS[1];
  trkKL[0] = &baseKin.trkKL[0];
  trkKL[1] = &baseKin.trkKL[1];
  trkKSmc[0] = &baseKin.trkKSmc[0];
  trkKSmc[1] = &baseKin.trkKSmc[1];
  trkKLmc[0] = &baseKin.trkKLmc[0];
  trkKLmc[1] = &baseKin.trkKLmc[1];

  baseKin.ipKL.resize(3);
  baseKin.ipKL.resize(3);
  baseKin.ipmc.resize(3);

  chain.SetBranchAddress("ipKS", &ipKS);
  chain.SetBranchAddress("ipKL", &ipKL);
  chain.SetBranchAddress("ipmc", &ipmc);
  // -----------------------------------------------------------

  // mcflag and mctruth position
  Int_t mcflag, mctruth;
  chain.SetBranchAddress("mcflag", &mcflag);
  chain.SetBranchAddress("mctruth", &mctruth);
  // -----------------------------------------------------------

  // Trk momentum
  baseKin.KchrecKS.resize(9);
  baseKin.KchrecKL.resize(9);
  baseKin.KchboostKS.resize(9);
  baseKin.KchboostKL.resize(9);
  baseKin.trkKS[0].resize(4);
  baseKin.trkKS[1].resize(4);
  baseKin.trkKL[0].resize(4);
  baseKin.trkKL[1].resize(4);

  chain.SetBranchAddress("KchrecKS", &KchrecKS);
  chain.SetBranchAddress("KchrecKL", &KchrecKL);
  chain.SetBranchAddress("KchboostKS", &KchboostKS);
  chain.SetBranchAddress("KchboostKL", &KchboostKL);
  chain.SetBranchAddress("trk1KS", &trkKS[0]);
  chain.SetBranchAddress("trk2KS", &trkKS[1]);
  chain.SetBranchAddress("trk1KL", &trkKL[0]);
  chain.SetBranchAddress("trk2KL", &trkKL[1]);

  baseKin.Kchmc.resize(9);
  baseKin.Knemc.resize(9);
  baseKin.trkKSmc[0].resize(4);
  baseKin.trkKSmc[1].resize(4);
  baseKin.trkKLmc[0].resize(4);
  baseKin.trkKLmc[1].resize(4);

  chain.SetBranchAddress("Kchmc", &Kchmc);
  chain.SetBranchAddress("Knemc", &Knemc);
  chain.SetBranchAddress("trk1KSmc", &trkKSmc[0]);
  chain.SetBranchAddress("trk2KSmc", &trkKSmc[1]);
  chain.SetBranchAddress("trk1KLmc", &trkKLmc[0]);
  chain.SetBranchAddress("trk2KLmc", &trkKLmc[1]);
  // -----------------------------------------------------------

  Float_t EmissKS = 0.0;
  Float_t EmissKL = 0.0;
  Float_t PmissKS = 0.0;
  Float_t PmissKL = 0.0;
  Float_t KchrecKSMom = 0.0;
  Float_t KchrecKLMom = 0.0;

  std::string cutFileName = "/data/ssd/gamrat/KLOE/Subanalysis/Properties/cut-limits-final.json";
  StatisticalCutter cutter(cutFileName, 7, KLOE::HypothesisCode::FOUR_PI);

  Float_t pKTwoBody = Obj.TwoBodyDecayMass(PhysicsConstants::mPhi, PhysicsConstants::mK0, PhysicsConstants::mK0);

  ///////////////////////////////////////////////////////////////////
  cutter.RegisterVariableGetter("InvMassKS", [&]()
                                { return KchrecKS->at(5); });
  cutter.RegisterCentralValueGetter("InvMassKS", [&]()
                                    { return PhysicsConstants::mK0; });
  ///////////////////////////////////////////////////////////////////
  cutter.RegisterVariableGetter("InvMassKL", [&]()
                                { return KchrecKL->at(5); });
  cutter.RegisterCentralValueGetter("InvMassKL", [&]()
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

  std::string
      dirname = (std::string)Paths::charged_dir + (std::string)Paths::root_files_dir,
      dated_folder = Obj.CreateDatedFolder(dirname);

  std::string base_filename = "KchRec_FourPiTwoBody";

  switch (dataType)
  {
  case Controls::DataType::DATA_ONLY:
    base_filename = base_filename + "_DATA";
    config.setProperty<std::string>("variables.tree.filepath.Data.fourPiTwoBody", dated_folder + "/" + base_filename + "*.root");
    break;
  case Controls::DataType::MC_ONLY:
    base_filename = base_filename + "_MC";
    config.setProperty<std::string>("variables.tree.filepath.MC.fourPiTwoBody", dated_folder + "/" + base_filename + "*.root");
    break;
  }

  SplitFileWriter writer(base_filename, 1.5 * 1024 * 1024 * 1024, false, dated_folder);

  Float_t
      invMass = 0.;

  Int_t nentries = chain.GetEntries();

  Int_t mode = 1; // Model for pi+pi-

  baseKin.KchrecKLTwoBody.resize(9);

  std::vector<std::vector<TH1 *>> histMomPiKLTwoBody(2);
  std::vector<std::vector<TH1 *>> histMomtrkKL(2);

  std::vector<TH1 *> histKLTwoBody;
  std::vector<TH1 *> histKLBoost;
  std::vector<TH1 *> histKL;

  for (Int_t i = 0; i < 2; i++)
    for (Int_t j = 0; j < 4; j++)
    {
      if (j < 3)
      {
        histMomPiKLTwoBody[i].push_back(new TH1F(Form("histMomPiKLTwoBody_%d_%d", i, j), Form("Histogram Pi KL Two Body %d %d", i, j), 50, -5, 5));
        histMomtrkKL[i].push_back(new TH1F(Form("histMomtrkKL_%d_%d", i, j), Form("Histogram trk KL %d %d", i, j), 50, -5, 5));
        histKLTwoBody.push_back(new TH1F(Form("histKLTwoBody_%d_%d", i, j), Form("Histogram KL Two Body %d %d", i, j), 50, -5, 5));
        histKLBoost.push_back(new TH1F(Form("histKLBoost_%d_%d", i, j), Form("Histogram KL Boost %d %d", i, j), 50, -5, 5));
        histKL.push_back(new TH1F(Form("histKL_%d_%d", i, j), Form("Histogram KL %d %d", i, j), 50, -5, 5));
      }
      else
      {
        histMomPiKLTwoBody[i].push_back(new TH1F(Form("histMomPiKLTwoBody_%d_%d", i, j), Form("Histogram Pi KL Two Body %d %d", i, j), 50, -5, 5));
        histMomtrkKL[i].push_back(new TH1F(Form("histMomtrkKL_%d_%d", i, j), Form("Histogram trk KL %d %d", i, j), 50, -5, 5));
        histKLTwoBody.push_back(new TH1F(Form("histKLTwoBody_%d_%d", i, j), Form("Histogram KL Two Body %d %d", i, j), 50, -5, 5));
        histKLBoost.push_back(new TH1F(Form("histKLBoost_%d_%d", i, j), Form("Histogram KL Boost %d %d", i, j), 50, -5, 5));
        histKL.push_back(new TH1F(Form("histKL_%d_%d", i, j), Form("Histogram KL %d %d", i, j), 50, -5, 5));
      }
    }

  std::vector<std::vector<TH1 *>> histMomPiKSTwoBody(2);
  std::vector<std::vector<TH1 *>> histMomtrkKS(2);

  std::vector<TH1 *> histKSTwoBody;
  std::vector<TH1 *> histKSBoost;
  std::vector<TH1 *> histKS;

  for (Int_t i = 0; i < 2; i++)
    for (Int_t j = 0; j < 4; j++)
    {
      if (j < 3)
      {
        histMomPiKSTwoBody[i].push_back(new TH1F(Form("histMomPiKSTwoBody_%d_%d", i, j), Form("Histogram Pi KS Two Body %d %d", i, j), 50, -5, 5));
        histMomtrkKS[i].push_back(new TH1F(Form("histMomtrkKS_%d_%d", i, j), Form("Histogram trk KS %d %d", i, j), 50, -5, 5));
        histKSTwoBody.push_back(new TH1F(Form("histKSTwoBody_%d_%d", i, j), Form("Histogram KS Two Body %d %d", i, j), 50, -5, 5));
        histKSBoost.push_back(new TH1F(Form("histKSBoost_%d_%d", i, j), Form("Histogram KS Boost %d %d", i, j), 50, -5, 5));
        histKS.push_back(new TH1F(Form("histKS_%d_%d", i, j), Form("Histogram KS %d %d", i, j), 50, -5, 5));
      }
      else
      {
        histMomPiKSTwoBody[i].push_back(new TH1F(Form("histMomPiKSTwoBody_%d_%d", i, j), Form("Histogram Pi KS Two Body %d %d", i, j), 50, -5, 5));
        histMomtrkKS[i].push_back(new TH1F(Form("histMomtrkKS_%d_%d", i, j), Form("Histogram trk KS %d %d", i, j), 50, -5, 5));
        histKSTwoBody.push_back(new TH1F(Form("histKSTwoBody_%d_%d", i, j), Form("Histogram KS Two Body %d %d", i, j), 50, -5, 5));
        histKSBoost.push_back(new TH1F(Form("histKSBoost_%d_%d", i, j), Form("Histogram KS Boost %d %d", i, j), 50, -5, 5));
        histKS.push_back(new TH1F(Form("histKS_%d_%d", i, j), Form("Histogram KS %d %d", i, j), 50, -5, 5));
      }
    }

  Int_t graph_flag = 0;

  Float_t minDiff = 999999.;

  // Progress bar
  boost::progress_display show_progress(nentries);
  // --------------------------------------------------------------------

  std::vector<Float_t>
      trkKLTwoBody1(4),
      trkKLTwoBody2(4),
      trkKSTwoBody1(4),
      trkKSTwoBody2(4);

  Float_t gammaKS = 0.0, gammaKL = 0.0;

  for (Int_t i = 0; i < nentries; i++)
  {
    chain.GetEntry(i);

    Float_t
        PhiMom[3] = {baseKin.phi_mom[0], baseKin.phi_mom[1], baseKin.phi_mom[2]},
        MissMomKS[3] = {},
        MissMomKL[3] = {};

    for (Int_t comp = 0; comp < 3; comp++)
    {
      MissMomKS[comp] = PhiMom[comp] - KchboostKS->at(comp) - KchrecKL->at(comp);
      MissMomKL[comp] = PhiMom[comp] - KchboostKL->at(comp) - KchrecKS->at(comp);
    }

    PmissKS = sqrt(pow(MissMomKS[0], 2) + pow(MissMomKS[1], 2) + pow(MissMomKS[2], 2));
    PmissKL = sqrt(pow(MissMomKL[0], 2) + pow(MissMomKL[1], 2) + pow(MissMomKL[2], 2));

    EmissKS = KchboostKS->at(3) - KchrecKS->at(3);
    EmissKL = KchboostKL->at(3) - KchrecKL->at(3);

    Float_t
        boostPhi[3] = {
            -baseKin.phi_mom[0] / baseKin.phi_mom[3],
            -baseKin.phi_mom[1] / baseKin.phi_mom[3],
            -baseKin.phi_mom[2] / baseKin.phi_mom[3]},
        phiMom[4] = {baseKin.phi_mom[0], baseKin.phi_mom[1], baseKin.phi_mom[2], baseKin.phi_mom[3]}, trkKS_PhiCM[2][4] = {}, KchrecKS_PhiCM[4] = {}, trkKL_PhiCM[2][4], KchrecKL_PhiCM[4] = {};

    Obj.lorentz_transf(boostPhi, trkKS[0]->data(), trkKS_PhiCM[0]);
    Obj.lorentz_transf(boostPhi, trkKS[1]->data(), trkKS_PhiCM[1]);
    Obj.lorentz_transf(boostPhi, trkKL[0]->data(), trkKL_PhiCM[0]);
    Obj.lorentz_transf(boostPhi, trkKL[1]->data(), trkKL_PhiCM[1]);

    for (Int_t part = 0; part < 2; part++)
      for (Int_t comp = 0; comp < 4; comp++)
      {
        KchrecKS_PhiCM[comp] += trkKS_PhiCM[part][comp];
        KchrecKL_PhiCM[comp] += trkKL_PhiCM[part][comp];
      }

    KchrecKSMom = sqrt(pow(KchrecKS_PhiCM[0], 2) + pow(KchrecKS_PhiCM[1], 2) + pow(KchrecKS_PhiCM[2], 2));
    KchrecKLMom = sqrt(pow(KchrecKL_PhiCM[0], 2) + pow(KchrecKL_PhiCM[1], 2) + pow(KchrecKL_PhiCM[2], 2));

    // Finding Ks from KSKL->pi+pi-pi+pi-
    int
        findKl = 0,
        findKs = 1,
        findClose = 0;

    int
        last_vtx = 0;

    minDiff = 999999.;

    // Set the proper data type to

    Bool_t data_flag = false;

    Obj.dataFlagSetter(dataType, data_flag, mcflag, mctruth);
    // -------------------------------------------------------------------

    baseKin.errFlag = 1;
    baseKin.errFlagKS = 1;
    baseKin.errFlagKL = 1;
    baseKin.errFlagClosest = 1;

    baseKin.KchrecKLTwoBody.clear();

    for (Int_t i = 0; i < 2; i++)
    {
      baseKin.trkKLTwoBody[i].clear();
    }

    TLorentzVector
        PiKL4VecLAB[2],
        trkKL4VecLAB[2],
        PiKS4VecLAB[2],
        trkKS4VecLAB[2];

    baseKin.KchrecKLTwoBody.resize(9);
    baseKin.KchrecKSTwoBody.resize(9);

    if (cutter.PassAllCuts())
    {

      TwoBodyReconstruction(KchboostKL, ipKL, trkKL, Obj, baseKin.KchrecKLTwoBody, gammaKL, PiKL4VecLAB, trkKL4VecLAB);

      TwoBodyReconstruction(KchboostKS, ipKS, trkKS, Obj, baseKin.KchrecKSTwoBody, gammaKS, PiKS4VecLAB, trkKS4VecLAB);

      for (Int_t k = 0; k < 4; k++)
      {
        trkKLTwoBody1[k] = PiKL4VecLAB[0][k];
        trkKLTwoBody2[k] = PiKL4VecLAB[1][k];

        trkKSTwoBody1[k] = PiKS4VecLAB[0][k];
        trkKSTwoBody2[k] = PiKS4VecLAB[1][k];

        baseKin.KchrecKLTwoBody[k] = trkKLTwoBody1[k] + trkKLTwoBody2[k];
        baseKin.KchrecKSTwoBody[k] = trkKSTwoBody1[k] + trkKSTwoBody2[k];
      }

      baseKin.KchrecKLTwoBody[4] = sqrt(pow(baseKin.KchrecKLTwoBody[0], 2) +
                                        pow(baseKin.KchrecKLTwoBody[1], 2) +
                                        pow(baseKin.KchrecKLTwoBody[2], 2));
      baseKin.KchrecKLTwoBody[5] = sqrt(pow(baseKin.KchrecKLTwoBody[3], 2) -
                                        pow(baseKin.KchrecKLTwoBody[4], 2));

      baseKin.KchrecKSTwoBody[4] = sqrt(pow(baseKin.KchrecKSTwoBody[0], 2) +
                                        pow(baseKin.KchrecKSTwoBody[1], 2) +
                                        pow(baseKin.KchrecKSTwoBody[2], 2));
      baseKin.KchrecKSTwoBody[5] = sqrt(pow(baseKin.KchrecKSTwoBody[3], 2) -
                                        pow(baseKin.KchrecKSTwoBody[4], 2));

      for (Int_t k = 0; k < 4; k++)
      {

        if (mcflag == 1 && mctruth == 7)
        {
          if (trkKLTwoBody1[k] - trkKLmc[0]->at(k) < trkKLTwoBody2[k] - trkKLmc[0]->at(k))
          {
            histMomPiKLTwoBody[0][k]->Fill(trkKLTwoBody1[k] - trkKLmc[0]->at(k));
            histMomPiKLTwoBody[1][k]->Fill(trkKLTwoBody2[k] - trkKLmc[1]->at(k));
          }
          else
          {
            histMomPiKLTwoBody[0][k]->Fill(trkKLTwoBody1[k] - trkKLmc[1]->at(k));
            histMomPiKLTwoBody[1][k]->Fill(trkKLTwoBody2[k] - trkKLmc[0]->at(k));
          }

          if (trkKL4VecLAB[0][k] - trkKLmc[0]->at(k) < trkKL4VecLAB[1][k] - trkKLmc[0]->at(k))
          {
            histMomtrkKL[0][k]->Fill(trkKL4VecLAB[0][k] - trkKLmc[0]->at(k));
            histMomtrkKL[1][k]->Fill(trkKL4VecLAB[1][k] - trkKLmc[1]->at(k));
          }
          else
          {
            histMomtrkKL[0][k]->Fill(trkKL4VecLAB[0][k] - trkKLmc[1]->at(k));
            histMomtrkKL[1][k]->Fill(trkKL4VecLAB[1][k] - trkKLmc[0]->at(k));
          }

          histKLTwoBody[k]->Fill(baseKin.KchrecKLTwoBody[k] - Knemc->at(k));
          histKLBoost[k]->Fill(KchboostKL->at(k) - Knemc->at(k));
          histKL[k]->Fill(KchrecKL->at(k) - Knemc->at(k));
        }

        if (mcflag == 1 && mctruth == 7)
        {
          if (trkKSTwoBody1[k] - trkKSmc[0]->at(k) < trkKSTwoBody2[k] - trkKSmc[0]->at(k))
          {
            histMomPiKSTwoBody[0][k]->Fill(trkKSTwoBody1[k] - trkKSmc[0]->at(k));
            histMomPiKSTwoBody[1][k]->Fill(trkKSTwoBody2[k] - trkKSmc[1]->at(k));
          }
          else
          {
            histMomPiKSTwoBody[0][k]->Fill(trkKSTwoBody1[k] - trkKSmc[1]->at(k));
            histMomPiKSTwoBody[1][k]->Fill(trkKSTwoBody2[k] - trkKSmc[0]->at(k));
          }

          if (trkKS4VecLAB[0][k] - trkKSmc[0]->at(k) < trkKS4VecLAB[1][k] - trkKSmc[0]->at(k))
          {
            histMomtrkKS[0][k]->Fill(trkKS4VecLAB[0][k] - trkKSmc[0]->at(k));
            histMomtrkKS[1][k]->Fill(trkKS4VecLAB[1][k] - trkKSmc[1]->at(k));
          }
          else
          {
            histMomtrkKS[0][k]->Fill(trkKS4VecLAB[0][k] - trkKSmc[1]->at(k));
            histMomtrkKS[1][k]->Fill(trkKS4VecLAB[1][k] - trkKSmc[0]->at(k));
          }

          histKSTwoBody[k]->Fill(baseKin.KchrecKSTwoBody[k] - Kchmc->at(k));
          histKSBoost[k]->Fill(KchboostKS->at(k) - Kchmc->at(k));
          histKS[k]->Fill(KchrecKS->at(k) - Kchmc->at(k));
        }
      }

      // Int_t zmienne
      std::map<std::string, Int_t> intVars = {
          {"nrun", baseKin.nrun},
          {"nev", baseKin.nev},
          {"mcflag", mcflag},
          {"mctruth", mctruth}};

      std::map<std::string, Float_t> floatVars = {
          {"gammaKS", gammaKS},
          {"gammaKL", gammaKL}};

      // Tablice
      std::map<std::string, std::vector<Int_t>> intArrays = {};

      std::map<std::string, std::vector<Float_t>> floatArrays = {
          {"Kchmc", *Kchmc},
          {"Knemc", *Knemc},
          {"KchrecKS", *KchrecKS},
          {"KchrecKL", *KchrecKL},
          {"KchboostKS", *KchboostKS},
          {"KchboostKL", *KchboostKL},
          {"KchrecKSTwoBody", baseKin.KchrecKSTwoBody},
          {"KchrecKLTwoBody", baseKin.KchrecKLTwoBody},
          {"trk1KS", *trkKS[0]},
          {"trk2KS", *trkKS[1]},
          {"trk1KL", *trkKL[0]},
          {"trk2KL", *trkKL[1]},
          {"trk1KSTwoBody", trkKSTwoBody1},
          {"trk2KSTwoBody", trkKSTwoBody2},
          {"trk1KLTwoBody", trkKLTwoBody1},
          {"trk2KLTwoBody", trkKLTwoBody2},
          {"trk1KSmc", *trkKSmc[0]},
          {"trk2KSmc", *trkKSmc[1]},
          {"trk1KLmc", *trkKLmc[0]},
          {"trk2KLmc", *trkKLmc[1]}};

      writer.Fill(intVars, floatVars, intArrays, floatArrays);
    }

    ++show_progress; // Progress of the loading bar
  }

  if (dataType == Controls::DataType::MC_ONLY)
  {
    std::string
        hist_dir = (std::string)Paths::charged_dir + (std::string)Paths::img_dir;

    std::string dated_hist_dir = Obj.CreateDatedFolder(hist_dir);

    TCanvas *canvas = new TCanvas("canvas", "Canvas", 790, 790);

    std::string xTitle[4] = {"p_{x} [MeV/c]", "p_{y} [MeV/c]", "p_{z} [MeV/c]", "E [MeV]"};

    gStyle->SetOptStat(0);

    for (Int_t j = 0; j < 4; j++)
    {
      for (Int_t i = 0; i < 2; i++)
      {
        TLegend *legendPi = new TLegend(0.15, 0.8, 0.6, 0.9);
        legendPi->AddEntry(histMomtrkKL[i][j], Form("Reconstructed charged pion %d", i), "l");
        legendPi->AddEntry(histMomPiKLTwoBody[i][j], Form("Charged pion %d from 2-body decay", i), "l");
        histMomPiKLTwoBody[i][j]->SetTitle(Form("Histogram Pi Two Body %d %d", i, j));
        histMomPiKLTwoBody[i][j]->GetXaxis()->SetTitle(xTitle[j].c_str());
        histMomPiKLTwoBody[i][j]->GetYaxis()->SetTitle("Counts");
        histMomPiKLTwoBody[i][j]->GetYaxis()->SetRangeUser(0, histMomPiKLTwoBody[i][j]->GetMaximum() * 2.0);
        histMomPiKLTwoBody[i][j]->GetYaxis()->SetMaxDigits(3);
        histMomPiKLTwoBody[i][j]->SetLineColor(kBlue);
        histMomPiKLTwoBody[i][j]->Draw();
        histMomtrkKL[i][j]->SetLineColor(kRed);
        histMomtrkKL[i][j]->Draw("SAME");
        // Triple Gaussian fit
        TF1 *fitFunc = new TF1("fitFunc", triple_gaus, histMomPiKLTwoBody[i][j]->GetXaxis()->GetXmin(), histMomPiKLTwoBody[i][j]->GetXaxis()->GetXmax(), 9);
        fitFunc->SetNpx(1000);

        fitFunc->SetParLimits(0, 0.0, 1000000.0);
        fitFunc->SetParLimits(1, -3.0, 3.0);
        fitFunc->SetParLimits(2, 0.0, 10.0);

        fitFunc->SetParLimits(3, 0.0, 1000000.0);
        fitFunc->SetParLimits(4, -3.0, 3.0);
        fitFunc->SetParLimits(5, 0.0, 10.0);

        fitFunc->SetParLimits(6, 0.0, 1000000.0);
        fitFunc->SetParLimits(7, -3.0, 3.0);
        fitFunc->SetParLimits(8, 0.0, 10.0);

        fitFunc->SetParameter(0, 1000.0);
        fitFunc->SetParameter(1, 0.0);
        fitFunc->SetParameter(2, 1.0);

        fitFunc->SetParameter(3, 1000.0);
        fitFunc->SetParameter(4, 2.0);
        fitFunc->SetParameter(5, 1.0);

        fitFunc->SetParameter(6, 1000.0);
        fitFunc->SetParameter(7, -2.0);
        fitFunc->SetParameter(8, 1.0);

        // Triple Gaussian fit
        TF1 *fitFunc1 = new TF1("fitFunc1", triple_gaus, histMomPiKLTwoBody[i][j]->GetXaxis()->GetXmin(), histMomPiKLTwoBody[i][j]->GetXaxis()->GetXmax(), 9);
        fitFunc1->SetNpx(1000);

        fitFunc1->SetParLimits(0, 0.0, 1000000.0);
        fitFunc1->SetParLimits(1, -3.0, 3.0);
        fitFunc1->SetParLimits(2, 0.0, 10.0);

        fitFunc1->SetParLimits(3, 0.0, 1000000.0);
        fitFunc1->SetParLimits(4, -3.0, 3.0);
        fitFunc1->SetParLimits(5, 0.0, 10.0);

        fitFunc1->SetParLimits(6, 0.0, 1000000.0);
        fitFunc1->SetParLimits(7, -3.0, 3.0);
        fitFunc1->SetParLimits(8, 0.0, 10.0);

        fitFunc1->SetParameter(0, 1000.0);
        fitFunc1->SetParameter(1, 0.0);
        fitFunc1->SetParameter(2, 1.0);

        fitFunc1->SetParameter(3, 1000.0);
        fitFunc1->SetParameter(4, 2.0);
        fitFunc1->SetParameter(5, 1.0);

        fitFunc1->SetParameter(6, 1000.0);
        fitFunc1->SetParameter(7, -2.0);
        fitFunc1->SetParameter(8, 1.0);

        TFitResultPtr fitRes = histMomPiKLTwoBody[i][j]->Fit(fitFunc, "S");
        TFitResultPtr fitRes1 = histMomtrkKL[i][j]->Fit(fitFunc1, "S");
        // Get fit parameters and errors
        Double_t params[9], errors[9];
        for (int k = 0; k < 9; ++k)
        {
          params[k] = fitFunc->GetParameter(k);
          errors[k] = fitFunc->GetParError(k);
        }

        // Annotate fit results
        TLatex latex;
        latex.SetNDC();
        latex.SetTextSize(0.03);
        latex.DrawLatex(0.15, 0.75 - 0.05, Form("Two body / boost: #mu=%.2f#pm%.2f, #sigma=%.2f#pm%.2f", comb_mean(params, errors), comb_mean_err(params, errors), comb_std_dev(params, errors), comb_std_dev_err(params, errors)));

        Double_t params1[9], errors1[9];
        for (int k = 0; k < 9; ++k)
        {
          params1[k] = fitFunc1->GetParameter(k);
          errors1[k] = fitFunc1->GetParError(k);
        }

        // Annotate fit results
        TLatex latex1;
        latex1.SetNDC();
        latex1.SetTextSize(0.03);
        latex1.DrawLatex(0.15, 0.75, Form("Tracks: #mu=%.2f#pm%.2f, #sigma=%.2f#pm%.2f", comb_mean(params1, errors1), comb_mean_err(params1, errors1), comb_std_dev(params1, errors1), comb_std_dev_err(params1, errors1)));

        legendPi->Draw();
        canvas->Print(Form((dated_hist_dir + "/hist_mom_piKL_two_body_%d_%d.png").c_str(), i, j));
      }

      TLegend *legendKaon = new TLegend(0.15, 0.8, 0.6, 0.9);
      legendKaon->AddEntry(histKL[j], "Reconstructed KL", "l");
      legendKaon->AddEntry(histKLBoost[j], "KL from boost method", "l");
      legendKaon->AddEntry(histKLTwoBody[j], "KL from 2-body decay", "l");

      histKLTwoBody[j]->SetTitle(Form("Histogram KL Two Body %d", j));
      histKLTwoBody[j]->GetXaxis()->SetTitle(xTitle[j].c_str());
      histKLTwoBody[j]->GetYaxis()->SetTitle("Counts");
      histKLTwoBody[j]->GetYaxis()->SetRangeUser(0, histKLTwoBody[j]->GetMaximum() * 2.0);
      histKLTwoBody[j]->GetYaxis()->SetMaxDigits(3);
      histKLTwoBody[j]->SetLineColor(kBlue);
      histKLTwoBody[j]->Draw();
      histKL[j]->SetLineColor(kRed);
      histKL[j]->Draw("SAME");
      // Triple Gaussian fit
      TF1 *fitFuncKL = new TF1("fitFuncKL", triple_gaus, histKLTwoBody[j]->GetXaxis()->GetXmin(), histKLTwoBody[j]->GetXaxis()->GetXmax(), 9);
      fitFuncKL->SetNpx(1000);

      fitFuncKL->SetParLimits(0, 0.0, 1000000.0);
      fitFuncKL->SetParLimits(1, -3.0, 3.0);
      fitFuncKL->SetParLimits(2, 0.0, 10.0);

      fitFuncKL->SetParLimits(3, 0.0, 1000000.0);
      fitFuncKL->SetParLimits(4, -3.0, 3.0);
      fitFuncKL->SetParLimits(5, 0.0, 10.0);

      fitFuncKL->SetParLimits(6, 0.0, 1000000.0);
      fitFuncKL->SetParLimits(7, -3.0, 3.0);
      fitFuncKL->SetParLimits(8, 0.0, 10.0);

      fitFuncKL->SetParameter(0, 1000.0);
      fitFuncKL->SetParameter(1, 0.0);
      fitFuncKL->SetParameter(2, 1.0);

      fitFuncKL->SetParameter(3, 1000.0);
      fitFuncKL->SetParameter(4, 2.0);
      fitFuncKL->SetParameter(5, 1.0);

      fitFuncKL->SetParameter(6, 1000.0);
      fitFuncKL->SetParameter(7, -2.0);
      fitFuncKL->SetParameter(8, 1.0);

      TFitResultPtr fitResKL = histKLTwoBody[j]->Fit(fitFuncKL, "S");
      Double_t paramsKL[9], errorsKL[9];
      for (int k = 0; k < 9; ++k)
      {
        paramsKL[k] = fitFuncKL->GetParameter(k);
        errorsKL[k] = fitFuncKL->GetParError(k);
      }
      TLatex latexKL;
      latexKL.SetNDC();
      latexKL.SetTextSize(0.03);
      latexKL.DrawLatex(0.15, 0.75 - 0.05, Form("Two body / boost: #mu=%.2f#pm%.2f, #sigma=%.2f#pm%.2f", comb_mean(paramsKL, errorsKL), comb_mean_err(paramsKL, errorsKL), comb_std_dev(paramsKL, errorsKL), comb_std_dev_err(paramsKL, errorsKL)));

      // Triple Gaussian fit
      TF1 *fitFuncKL1 = new TF1("fitFuncKL1", triple_gaus, histKL[j]->GetXaxis()->GetXmin(), histKL[j]->GetXaxis()->GetXmax(), 9);
      fitFuncKL1->SetNpx(1000);

      fitFuncKL1->SetParLimits(0, 0.0, 1000000.0);
      fitFuncKL1->SetParLimits(1, -3.0, 3.0);
      fitFuncKL1->SetParLimits(2, 0.0, 10.0);

      fitFuncKL1->SetParLimits(3, 0.0, 1000000.0);
      fitFuncKL1->SetParLimits(4, -3.0, 3.0);
      fitFuncKL1->SetParLimits(5, 0.0, 10.0);

      fitFuncKL1->SetParLimits(6, 0.0, 1000000.0);
      fitFuncKL1->SetParLimits(7, -3.0, 3.0);
      fitFuncKL1->SetParLimits(8, 0.0, 10.0);

      fitFuncKL1->SetParameter(0, 1000.0);
      fitFuncKL1->SetParameter(1, 0.0);
      fitFuncKL1->SetParameter(2, 1.0);

      fitFuncKL1->SetParameter(3, 1000.0);
      fitFuncKL1->SetParameter(4, 2.0);
      fitFuncKL1->SetParameter(5, 1.0);

      fitFuncKL1->SetParameter(6, 1000.0);
      fitFuncKL1->SetParameter(7, -2.0);
      fitFuncKL1->SetParameter(8, 1.0);

      TFitResultPtr fitRes1 = histKL[j]->Fit(fitFuncKL1, "S");

      Double_t params1[9], errors1[9];
      for (int k = 0; k < 9; ++k)
      {
        params1[k] = fitFuncKL1->GetParameter(k);
        errors1[k] = fitFuncKL1->GetParError(k);
      }

      // Annotate fit results
      TLatex latex1;
      latex1.SetNDC();
      latex1.SetTextSize(0.03);
      latex1.DrawLatex(0.15, 0.75, Form("Tracks: #mu=%.2f#pm%.2f, #sigma=%.2f#pm%.2f", comb_mean(params1, errors1), comb_mean_err(params1, errors1), comb_std_dev(params1, errors1), comb_std_dev_err(params1, errors1)));

      histKLBoost[j]->SetLineColor(kGreen);
      histKLBoost[j]->Draw("SAME");
      legendKaon->Draw();
      canvas->Print(Form((dated_hist_dir + "/hist_kl_two_body_%d.png").c_str(), j));
    }

    for (Int_t j = 0; j < 4; j++)
    {
      for (Int_t i = 0; i < 2; i++)
      {
        TLegend *legendPi = new TLegend(0.15, 0.8, 0.6, 0.9);
        legendPi->AddEntry(histMomtrkKS[i][j], Form("Reconstructed charged pion %d", i), "l");
        legendPi->AddEntry(histMomPiKSTwoBody[i][j], Form("Charged pion %d from 2-body decay", i), "l");
        histMomPiKSTwoBody[i][j]->SetTitle(Form("Histogram Pi Two Body %d %d", i, j));
        histMomPiKSTwoBody[i][j]->GetXaxis()->SetTitle(xTitle[j].c_str());
        histMomPiKSTwoBody[i][j]->GetYaxis()->SetTitle("Counts");
        histMomPiKSTwoBody[i][j]->GetYaxis()->SetRangeUser(0, histMomPiKSTwoBody[i][j]->GetMaximum() * 2.0);
        histMomPiKSTwoBody[i][j]->GetYaxis()->SetMaxDigits(3);
        histMomPiKSTwoBody[i][j]->SetLineColor(kBlue);
        histMomPiKSTwoBody[i][j]->Draw();
        histMomtrkKS[i][j]->SetLineColor(kRed);
        histMomtrkKS[i][j]->Draw("SAME");
        // Triple Gaussian fit
        TF1 *fitFuncKS = new TF1("fitFuncKS", triple_gaus, histMomPiKSTwoBody[i][j]->GetXaxis()->GetXmin(), histMomPiKSTwoBody[i][j]->GetXaxis()->GetXmax(), 9);
        fitFuncKS->SetNpx(1000);

        fitFuncKS->SetParLimits(0, 0.0, 1000000.0);
        fitFuncKS->SetParLimits(1, -3.0, 3.0);
        fitFuncKS->SetParLimits(2, 0.0, 10.0);

        fitFuncKS->SetParLimits(3, 0.0, 1000000.0);
        fitFuncKS->SetParLimits(4, -3.0, 3.0);
        fitFuncKS->SetParLimits(5, 0.0, 10.0);

        fitFuncKS->SetParLimits(6, 0.0, 1000000.0);
        fitFuncKS->SetParLimits(7, -3.0, 3.0);
        fitFuncKS->SetParLimits(8, 0.0, 10.0);

        fitFuncKS->SetParameter(0, 1000.0);
        fitFuncKS->SetParameter(1, 0.0);
        fitFuncKS->SetParameter(2, 1.0);

        fitFuncKS->SetParameter(3, 1000.0);
        fitFuncKS->SetParameter(4, 2.0);
        fitFuncKS->SetParameter(5, 1.0);

        fitFuncKS->SetParameter(6, 1000.0);
        fitFuncKS->SetParameter(7, -2.0);
        fitFuncKS->SetParameter(8, 1.0);

        // Triple Gaussian fit
        TF1 *fitFunc1 = new TF1("fitFunc1", triple_gaus, histMomPiKLTwoBody[i][j]->GetXaxis()->GetXmin(), histMomPiKLTwoBody[i][j]->GetXaxis()->GetXmax(), 9);
        fitFunc1->SetNpx(1000);

        fitFunc1->SetParLimits(0, 0.0, 1000000.0);
        fitFunc1->SetParLimits(1, -3.0, 3.0);
        fitFunc1->SetParLimits(2, 0.0, 10.0);

        fitFunc1->SetParLimits(3, 0.0, 1000000.0);
        fitFunc1->SetParLimits(4, -3.0, 3.0);
        fitFunc1->SetParLimits(5, 0.0, 10.0);

        fitFunc1->SetParLimits(6, 0.0, 1000000.0);
        fitFunc1->SetParLimits(7, -3.0, 3.0);
        fitFunc1->SetParLimits(8, 0.0, 10.0);

        fitFunc1->SetParameter(0, 1000.0);
        fitFunc1->SetParameter(1, 0.0);
        fitFunc1->SetParameter(2, 1.0);

        fitFunc1->SetParameter(3, 1000.0);
        fitFunc1->SetParameter(4, 2.0);
        fitFunc1->SetParameter(5, 1.0);

        fitFunc1->SetParameter(6, 1000.0);
        fitFunc1->SetParameter(7, -2.0);
        fitFunc1->SetParameter(8, 1.0);

        TFitResultPtr fitRes = histMomPiKSTwoBody[i][j]->Fit(fitFuncKS, "S");
        TFitResultPtr fitRes1 = histMomtrkKS[i][j]->Fit(fitFunc1, "S");

        // Get fit parameters and errors
        Double_t params[9], errors[9];
        for (int k = 0; k < 9; ++k)
        {
          params[k] = fitFuncKS->GetParameter(k);
          errors[k] = fitFuncKS->GetParError(k);
        }

        // Annotate fit results
        TLatex latex;
        latex.SetNDC();
        latex.SetTextSize(0.03);
        latex.DrawLatex(0.15, 0.75 - 0.05, Form("Two body / boost: #mu=%.2f#pm%.2f, #sigma=%.2f#pm%.2f", comb_mean(params, errors), comb_mean_err(params, errors), comb_std_dev(params, errors), comb_std_dev_err(params, errors)));

        Double_t params1[9], errors1[9];
        for (int k = 0; k < 9; ++k)
        {
          params1[k] = fitFunc1->GetParameter(k);
          errors1[k] = fitFunc1->GetParError(k);
        }

        // Annotate fit results
        TLatex latex1;
        latex1.SetNDC();
        latex1.SetTextSize(0.03);
        latex1.DrawLatex(0.15, 0.75, Form("Tracks: #mu=%.2f#pm%.2f, #sigma=%.2f#pm%.2f", comb_mean(params1, errors1), comb_mean_err(params1, errors1), comb_std_dev(params1, errors1), comb_std_dev_err(params1, errors1)));

        legendPi->Draw();
        canvas->Print(Form((dated_hist_dir + "/hist_mom_piKS_two_body_%d_%d.png").c_str(), i, j));
      }

      TLegend *legendKaon = new TLegend(0.15, 0.8, 0.6, 0.9);
      legendKaon->AddEntry(histKS[j], "Reconstructed KS", "l");
      legendKaon->AddEntry(histKSBoost[j], "KS from boost method", "l");
      legendKaon->AddEntry(histKSTwoBody[j], "KS from 2-body decay", "l");

      histKSTwoBody[j]->SetTitle(Form("Histogram KS Two Body %d", j));
      histKSTwoBody[j]->GetXaxis()->SetTitle(xTitle[j].c_str());
      histKSTwoBody[j]->GetYaxis()->SetTitle("Counts");
      histKSTwoBody[j]->GetYaxis()->SetRangeUser(0, histKSTwoBody[j]->GetMaximum() * 2.0);
      histKSTwoBody[j]->GetYaxis()->SetMaxDigits(3);
      histKSTwoBody[j]->SetLineColor(kBlue);
      histKSTwoBody[j]->Draw();
      histKS[j]->SetLineColor(kRed);
      histKS[j]->Draw("SAME");
      // Triple Gaussian fit
      TF1 *fitFuncKSTwo = new TF1("fitFuncKSTwo", triple_gaus, histKSTwoBody[j]->GetXaxis()->GetXmin(), histKSTwoBody[j]->GetXaxis()->GetXmax(), 9);
      fitFuncKSTwo->SetNpx(1000);

      fitFuncKSTwo->SetParLimits(0, 0.0, 1000000.0);
      fitFuncKSTwo->SetParLimits(1, -3.0, 3.0);
      fitFuncKSTwo->SetParLimits(2, 0.0, 10.0);

      fitFuncKSTwo->SetParLimits(3, 0.0, 1000000.0);
      fitFuncKSTwo->SetParLimits(4, -3.0, 3.0);
      fitFuncKSTwo->SetParLimits(5, 0.0, 10.0);

      fitFuncKSTwo->SetParLimits(6, 0.0, 1000000.0);
      fitFuncKSTwo->SetParLimits(7, -3.0, 3.0);
      fitFuncKSTwo->SetParLimits(8, 0.0, 10.0);

      fitFuncKSTwo->SetParameter(0, 1000.0);
      fitFuncKSTwo->SetParameter(1, 0.0);
      fitFuncKSTwo->SetParameter(2, 1.0);

      fitFuncKSTwo->SetParameter(3, 1000.0);
      fitFuncKSTwo->SetParameter(4, 2.0);
      fitFuncKSTwo->SetParameter(5, 1.0);

      fitFuncKSTwo->SetParameter(6, 1000.0);
      fitFuncKSTwo->SetParameter(7, -2.0);
      fitFuncKSTwo->SetParameter(8, 1.0);

      TFitResultPtr fitResKSTwo = histKSTwoBody[j]->Fit(fitFuncKSTwo, "S");
      Double_t paramsKSTwo[9], errorsKSTwo[9];
      for (int k = 0; k < 9; ++k)
      {
        paramsKSTwo[k] = fitFuncKSTwo->GetParameter(k);
        errorsKSTwo[k] = fitFuncKSTwo->GetParError(k);
      }
      TLatex latexKSTwo;
      latexKSTwo.SetNDC();
      latexKSTwo.SetTextSize(0.03);

      latexKSTwo.DrawLatex(0.15, 0.75 - 0.05, Form("Two body / boost: #mu=%.2f#pm%.2f, #sigma=%.2f#pm%.2f", comb_mean(paramsKSTwo, errorsKSTwo), comb_mean_err(paramsKSTwo, errorsKSTwo), comb_std_dev(paramsKSTwo, errorsKSTwo), comb_std_dev_err(paramsKSTwo, errorsKSTwo)));

      // Triple Gaussian fit
      TF1 *fitFuncKS1 = new TF1("fitFuncKS1", triple_gaus, histKS[j]->GetXaxis()->GetXmin(), histKS[j]->GetXaxis()->GetXmax(), 9);
      fitFuncKS1->SetNpx(1000);

      fitFuncKS1->SetParLimits(0, 0.0, 1000000.0);
      fitFuncKS1->SetParLimits(1, -3.0, 3.0);
      fitFuncKS1->SetParLimits(2, 0.0, 10.0);

      fitFuncKS1->SetParLimits(3, 0.0, 1000000.0);
      fitFuncKS1->SetParLimits(4, -3.0, 3.0);
      fitFuncKS1->SetParLimits(5, 0.0, 10.0);

      fitFuncKS1->SetParLimits(6, 0.0, 1000000.0);
      fitFuncKS1->SetParLimits(7, -3.0, 3.0);
      fitFuncKS1->SetParLimits(8, 0.0, 10.0);

      fitFuncKS1->SetParameter(0, 1000.0);
      fitFuncKS1->SetParameter(1, 0.0);
      fitFuncKS1->SetParameter(2, 1.0);

      fitFuncKS1->SetParameter(3, 1000.0);
      fitFuncKS1->SetParameter(4, 2.0);
      fitFuncKS1->SetParameter(5, 1.0);

      fitFuncKS1->SetParameter(6, 1000.0);
      fitFuncKS1->SetParameter(7, -2.0);
      fitFuncKS1->SetParameter(8, 1.0);

      TFitResultPtr fitRes1 = histKS[j]->Fit(fitFuncKS1, "S");

      Double_t params1[9], errors1[9];
      for (int k = 0; k < 9; ++k)
      {
        params1[k] = fitFuncKS1->GetParameter(k);
        errors1[k] = fitFuncKS1->GetParError(k);
      }

      // Annotate fit results
      TLatex latex1;
      latex1.SetNDC();
      latex1.SetTextSize(0.03);
      latex1.DrawLatex(0.15, 0.75, Form("Tracks: #mu=%.2f#pm%.2f, #sigma=%.2f#pm%.2f", comb_mean(params1, errors1), comb_mean_err(params1, errors1), comb_std_dev(params1, errors1), comb_std_dev_err(params1, errors1)));

      histKSBoost[j]->SetLineColor(kGreen);
      histKSBoost[j]->Draw("SAME");
      legendKaon->Draw();
      canvas->Print(Form((dated_hist_dir + "/hist_ks_two_body_%d.png").c_str(), j));
    }
  }

  writer.Close();

  config.setProperty<std::string>("lastUpdate", Obj.getCurrentTimestamp());
  config.setProperty<std::string>("lastScript", "Two body KS-KL reconstruction.");

  config.saveProperties();

  return 0;
}

ErrorHandling::ErrorCodes TwoBodyReconstruction(std::vector<Float_t> *Kchboost, std::vector<Float_t> *ip, std::vector<Float_t> *trk[2], KLOE::pm00 &Obj, std::vector<Float_t> &KchrecTwoBody, Float_t &gamma, TLorentzVector PiKaon4VecLAB[2], TLorentzVector trk4VecLAB[2])
{
  // 2. Calculation of KL flight direction
  // Double_t KLpath = sqrt(pow(Kchboost->at(6) - ip->at(0), 2) +
  //                        pow(Kchboost->at(7) - ip->at(1), 2) +
  //                        pow(Kchboost->at(8) - ip->at(2), 2));
  // TVector3 KLflightDirection = {(Kchboost->at(6) - ip->at(0)) / KLpath,
  //                               (Kchboost->at(7) - ip->at(1)) / KLpath,
  //                               (Kchboost->at(8) - ip->at(2)) / KLpath};

  // // // 2.1 Calculation of KL momentum magnitude
  // Double_t KLmomMag = sqrt(pow(Kchboost->at(0), 2) +
  //                          pow(Kchboost->at(1), 2) +
  //                          pow(Kchboost->at(2), 2));

  // 2.2 Calculation of KL momentum from 2 body decay
  for (Int_t j = 0; j < 3; j++)
  {
    // Components
    KchrecTwoBody[j] = Kchboost->at(j); // KLflightDirection[j] * KLmomMag;
  }
  // Energy
  KchrecTwoBody[3] = Kchboost->at(3);

  // Calculate the magnitude of the KL momentum
  KchrecTwoBody[4] = sqrt(pow(KchrecTwoBody[0], 2) +
                          pow(KchrecTwoBody[1], 2) +
                          pow(KchrecTwoBody[2], 2));
  // Calculate the invariant mass of the KL
  KchrecTwoBody[5] = sqrt(pow(KchrecTwoBody[3], 2) -
                          pow(KchrecTwoBody[4], 2));

  // Copy the rest of the KL vertex components from the original KL vertex
  KchrecTwoBody[6] = Kchboost->at(6);
  KchrecTwoBody[7] = Kchboost->at(7);
  KchrecTwoBody[8] = Kchboost->at(8);
  // ------------------------------------------------------------------

  // 2. Go to KL CM frame
  TVector3
      PiMomKaonCM[2],           // Pion momenta in the KL CM frame
      PiMomKaonLAB[2],          // Pion momenta in the LAB frame
      z_axis = {0.0, 0.0, 1.0}, // z-axis for rotation,
      x_axis = {1.0, 0.0, 0.0}, // x-axis for rotation
      kaonMomLAB = {-KchrecTwoBody[0] / KchrecTwoBody[3],
                    -KchrecTwoBody[1] / KchrecTwoBody[3],
                    -KchrecTwoBody[2] / KchrecTwoBody[3]},
      Pi1MomLAB = {trk[0]->at(0),
                   trk[0]->at(1),
                   trk[0]->at(2)},
      Pi2MomLAB = {trk[1]->at(0),
                   trk[1]->at(1),
                   trk[1]->at(2)},
      nVec = Pi1MomLAB.Cross(Pi2MomLAB), // Normal vector to the plane of the two pions
      cross = nVec.Cross(z_axis),
      trkMomVecLAB[2], // Track momenta in the LAB frame
      trkMomVecKaonCM[2];

  TLorentzVector
      kaon4VecLAB,
      kaon4VecKaonCM,
      trk4VecKaonCM[2],
      PiKaon4VecKaonCM[2];

  Double_t rotAngle = nVec.Angle(z_axis);

  kaonMomLAB.Rotate(rotAngle, cross);

  Double_t rotAngleKaonX = kaonMomLAB.Angle(x_axis);

  kaonMomLAB.Rotate(rotAngleKaonX, z_axis);

  for (Int_t j = 0; j < 3; j++)
  {
    // Components of the track momenta in the LAB frame
    trkMomVecLAB[0][j] = trk[0]->at(j);
    trkMomVecLAB[1][j] = trk[1]->at(j);
  }

  trkMomVecLAB[0].Rotate(rotAngle, cross);
  trkMomVecLAB[0].Rotate(rotAngleKaonX, z_axis);

  trkMomVecLAB[1].Rotate(rotAngle, cross);
  trkMomVecLAB[1].Rotate(rotAngleKaonX, z_axis);

  trk4VecLAB[0].SetPxPyPzE(trkMomVecLAB[0][0],
                           trkMomVecLAB[0][1],
                           trkMomVecLAB[0][2],
                           trk[0]->at(3));

  trk4VecLAB[1].SetPxPyPzE(trkMomVecLAB[1][0],
                           trkMomVecLAB[1][1],
                           trkMomVecLAB[1][2],
                           trk[1]->at(3));

  Obj.lorentz_transf(kaonMomLAB, trk4VecLAB[0], trk4VecKaonCM[0]);
  Obj.lorentz_transf(kaonMomLAB, trk4VecLAB[1], trk4VecKaonCM[1]);

  trkMomVecKaonCM[0].SetXYZ(trk4VecKaonCM[0][0],
                            trk4VecKaonCM[0][1],
                            trk4VecKaonCM[0][2]);
  trkMomVecKaonCM[1].SetXYZ(trk4VecKaonCM[1][0],
                            trk4VecKaonCM[1][1],
                            trk4VecKaonCM[1][2]);

  Double_t
      PiMomMagKaonCM1 = Obj.TwoBodyDecayMass(KchrecTwoBody[5], PhysicsConstants::mPiCh, PhysicsConstants::mPiCh);

  Double_t angle1 = trkMomVecKaonCM[0].Phi(),
           angle2 = trkMomVecKaonCM[1].Phi(),
           theta1 = trkMomVecKaonCM[0].Theta(),
           theta2 = trkMomVecKaonCM[1].Theta();

  gamma = (M_PI_2 - 0.5 * abs(angle1) - 0.5 * abs(angle2));

  if (angle1 < 0.0)
    PiMomKaonCM[0].SetXYZ(PiMomMagKaonCM1 * sin(theta1) * cos(angle1 - gamma),
                          PiMomMagKaonCM1 * sin(theta1) * sin(angle1 - gamma),
                          PiMomMagKaonCM1 * cos(theta1));
  else
    PiMomKaonCM[0].SetXYZ(PiMomMagKaonCM1 * sin(theta1) * cos(angle1 + gamma),
                          PiMomMagKaonCM1 * sin(theta1) * sin(angle1 + gamma),
                          PiMomMagKaonCM1 * cos(theta1));
  if (angle2 < 0.0)
    PiMomKaonCM[1].SetXYZ(PiMomMagKaonCM1 * sin(theta2) * cos(angle2 - gamma),
                          PiMomMagKaonCM1 * sin(theta2) * sin(angle2 - gamma),
                          PiMomMagKaonCM1 * cos(theta2));
  else
    PiMomKaonCM[1].SetXYZ(PiMomMagKaonCM1 * sin(theta2) * cos(angle2 + gamma),
                          PiMomMagKaonCM1 * sin(theta2) * sin(angle2 + gamma),
                          PiMomMagKaonCM1 * cos(theta2));

  PiKaon4VecKaonCM[0].SetPxPyPzE(PiMomKaonCM[0][0], PiMomKaonCM[0][1], PiMomKaonCM[0][2], sqrt(pow(PiMomKaonCM[0][0], 2) + pow(PiMomKaonCM[0][1], 2) + pow(PiMomKaonCM[0][2], 2) + pow(PhysicsConstants::mPiCh, 2)));
  PiKaon4VecKaonCM[1].SetPxPyPzE(PiMomKaonCM[1][0], PiMomKaonCM[1][1], PiMomKaonCM[1][2], sqrt(pow(PiMomKaonCM[1][0], 2) + pow(PiMomKaonCM[1][1], 2) + pow(PiMomKaonCM[1][2], 2) + pow(PhysicsConstants::mPiCh, 2)));

  kaonMomLAB = -kaonMomLAB; // Invert the boost direction for the pions
  Obj.lorentz_transf(kaonMomLAB, PiKaon4VecKaonCM[0], PiKaon4VecLAB[0]);
  Obj.lorentz_transf(kaonMomLAB, PiKaon4VecKaonCM[1], PiKaon4VecLAB[1]);

  PiMomKaonLAB[0].SetXYZ(PiKaon4VecLAB[0][0],
                         PiKaon4VecLAB[0][1],
                         PiKaon4VecLAB[0][2]);
  PiMomKaonLAB[1].SetXYZ(PiKaon4VecLAB[1][0],
                         PiKaon4VecLAB[1][1],
                         PiKaon4VecLAB[1][2]);

  // Rotate the pion momenta to align with the kaon momentum direction

  PiMomKaonLAB[0].Rotate(-rotAngleKaonX, z_axis);
  PiMomKaonLAB[0].Rotate(-rotAngle, cross);

  PiMomKaonLAB[1].Rotate(-rotAngleKaonX, z_axis);
  PiMomKaonLAB[1].Rotate(-rotAngle, cross);

  for (Int_t j = 0; j < 3; j++)
  {
    // Components of the track momenta in the LAB frame
    trkMomVecLAB[0][j] = trk[0]->at(j);
    trkMomVecLAB[1][j] = trk[1]->at(j);
  }

  PiKaon4VecLAB[0].SetPxPyPzE(PiMomKaonLAB[0][0],
                              PiMomKaonLAB[0][1],
                              PiMomKaonLAB[0][2],
                              sqrt(PiMomKaonLAB[0].Mag2() + pow(PhysicsConstants::mPiCh, 2)));
  PiKaon4VecLAB[1].SetPxPyPzE(PiMomKaonLAB[1][0],
                              PiMomKaonLAB[1][1],
                              PiMomKaonLAB[1][2],
                              sqrt(PiMomKaonLAB[1].Mag2() + pow(PhysicsConstants::mPiCh, 2)));

  trk4VecLAB[0].SetPxPyPzE(trkMomVecLAB[0][0],
                           trkMomVecLAB[0][1],
                           trkMomVecLAB[0][2],
                           trk[0]->at(3));
  trk4VecLAB[1].SetPxPyPzE(trkMomVecLAB[1][0],
                           trkMomVecLAB[1][1],
                           trkMomVecLAB[1][2],
                           trk[1]->at(3));

  return ErrorHandling::ErrorCodes::NO_ERROR;
};