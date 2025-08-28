// Author: Szymon Gamrat
// Date of last update: 03.02.2025

#include <TGraph.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TFitResultPtr.h>
#include <TFitResult.h>
#include <TLegend.h>
#include <THStack.h>
#include <TTreeReader.h>

#include <boost/optional.hpp>
#include <SplitFileWriter.h>
#include <event_data.h>
#include <ConfigManager.h>

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
  // Structure from const.h to ease navigation
  BaseKinematics baseKin;

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

  std::string
      base_filename = "KchRec_4pi",
      dirname = (std::string)charged_dir + (std::string)root_files_dir,
      dated_folder = Obj.CreateDatedFolder(dirname);

  SplitFileWriter writer(base_filename, 1.5 * 1024 * 1024 * 1024, false, dated_folder);

  Float_t
      invMass = 0.;

  Int_t nentries = chain.GetEntries();

  Int_t mode = 1; // Model for pi+pi-

  baseKin.KchrecKLTwoBody.resize(9);

  TH1 *histDifferences = new TH1F("histDifferences", "Histogram", 100, 0, M_PI);

  std::vector<std::vector<TH1 *>> histMomPiTwoBody(2);
  std::vector<std::vector<TH1 *>> histMomtrkKL(2);

  std::vector<TH1 *> histKLTwoBody;
  std::vector<TH1 *> histKL;

  for (Int_t i = 0; i < 2; i++)
    for (Int_t j = 0; j < 4; j++)
    {
      if (j < 3)
      {
        histMomPiTwoBody[i].push_back(new TH1F(Form("histMomPiTwoBody_%d_%d", i, j), Form("Histogram Pi Two Body %d %d", i, j), 50, -5, 5));
        histMomtrkKL[i].push_back(new TH1F(Form("histMomtrkKL_%d_%d", i, j), Form("Histogram trk KL %d %d", i, j), 50, -5, 5));
        histKLTwoBody.push_back(new TH1F(Form("histKLTwoBody_%d_%d", i, j), Form("Histogram KL Two Body %d %d", i, j), 50, -5, 5));
        histKL.push_back(new TH1F(Form("histKL_%d_%d", i, j), Form("Histogram KL %d %d", i, j), 50, -5, 5));
      }
      else
      {
        histMomPiTwoBody[i].push_back(new TH1F(Form("histMomPiTwoBody_%d_%d", i, j), Form("Histogram Pi Two Body %d %d", i, j), 50, -5, 5));
        histMomtrkKL[i].push_back(new TH1F(Form("histMomtrkKL_%d_%d", i, j), Form("Histogram trk KL %d %d", i, j), 50, -5, 5));
        histKLTwoBody.push_back(new TH1F(Form("histKLTwoBody_%d_%d", i, j), Form("Histogram KL Two Body %d %d", i, j), 50, -5, 5));
        histKL.push_back(new TH1F(Form("histKL_%d_%d", i, j), Form("Histogram KL %d %d", i, j), 50, -5, 5));
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

  Float_t gamma = 0.0, gamma_theta = 0.0;

  for (Int_t i = 0; i < nentries; i++)
  {
    chain.GetEntry(i);

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

    // 2. Calculation of KL flight direction
    Double_t KLpath = sqrt(pow(KchboostKL->at(6) - ipKL->at(0), 2) +
                           pow(KchboostKL->at(7) - ipKL->at(1), 2) +
                           pow(KchboostKL->at(8) - ipKL->at(2), 2));
    TVector3 KLflightDirection = {(KchboostKL->at(6) - ipKL->at(0)) / KLpath,
                                  (KchboostKL->at(7) - ipKL->at(1)) / KLpath,
                                  (KchboostKL->at(8) - ipKL->at(2)) / KLpath};

    // 2.1 Calculation of KL momentum magnitude
    Double_t KLmomMag = sqrt(pow(KchboostKL->at(0), 2) +
                             pow(KchboostKL->at(1), 2) +
                             pow(KchboostKL->at(2), 2));

    // 2.2 Calculation of KL momentum from 2 body decay
    for (Int_t j = 0; j < 3; j++)
    {
      // Components
      baseKin.KchrecKLTwoBody[j] = KLflightDirection[j] * KLmomMag;
    }
    // Energy
    baseKin.KchrecKLTwoBody[3] = KchboostKL->at(3);

    // Calculate the magnitude of the KL momentum
    baseKin.KchrecKLTwoBody[4] = sqrt(pow(baseKin.KchrecKLTwoBody[0], 2) +
                                      pow(baseKin.KchrecKLTwoBody[1], 2) +
                                      pow(baseKin.KchrecKLTwoBody[2], 2));
    // Calculate the invariant mass of the KL
    baseKin.KchrecKLTwoBody[5] = sqrt(pow(baseKin.KchrecKLTwoBody[3], 2) -
                                      pow(baseKin.KchrecKLTwoBody[4], 2));

    // Copy the rest of the KL vertex components from the original KL vertex
    baseKin.KchrecKLTwoBody[6] = KchboostKL->at(6);
    baseKin.KchrecKLTwoBody[7] = KchboostKL->at(7);
    baseKin.KchrecKLTwoBody[8] = KchboostKL->at(8);
    // ------------------------------------------------------------------

    // 2. Go to KL CM frame
    TVector3
        PiMomKaonCM[2],           // Pion momenta in the KL CM frame
        PiMomKaonLAB[2],          // Pion momenta in the LAB frame
        z_axis = {0.0, 0.0, 1.0}, // z-axis for rotation,
        x_axis = {1.0, 0.0, 0.0}, // x-axis for rotation
        kaonMomLAB = {-baseKin.KchrecKLTwoBody[0] / baseKin.KchrecKLTwoBody[3],
                      -baseKin.KchrecKLTwoBody[1] / baseKin.KchrecKLTwoBody[3],
                      -baseKin.KchrecKLTwoBody[2] / baseKin.KchrecKLTwoBody[3]},
        Pi1MomLAB = {trkKL[0]->at(0),
                     trkKL[0]->at(1),
                     trkKL[0]->at(2)},
        Pi2MomLAB = {trkKL[1]->at(0),
                     trkKL[1]->at(1),
                     trkKL[1]->at(2)},
        nVec = Pi1MomLAB.Cross(Pi2MomLAB), // Normal vector to the plane of the two pions
        cross = nVec.Cross(z_axis),
        trkKLMomVecLAB[2], // Track momenta in the LAB frame
        trkKLMomVecKaonCM[2];

    TLorentzVector
        kaon4VecLAB,
        kaon4VecKaonCM,
        trkKL4VecLAB[2],
        trkKL4VecKaonCM[2],
        PiKaon4VecKaonCM[2],
        PiKaon4VecLAB[2];

    Double_t rotAngle = nVec.Angle(z_axis);

    kaonMomLAB.Rotate(rotAngle, cross);

    Double_t rotAngleKaonX = kaonMomLAB.Angle(x_axis);

    kaonMomLAB.Rotate(rotAngleKaonX, z_axis);

    for (Int_t j = 0; j < 3; j++)
    {
      // Components of the track momenta in the LAB frame
      trkKLMomVecLAB[0][j] = trkKL[0]->at(j);
      trkKLMomVecLAB[1][j] = trkKL[1]->at(j);
    }

    trkKLMomVecLAB[0].Rotate(rotAngle, cross);
    trkKLMomVecLAB[0].Rotate(rotAngleKaonX, z_axis);

    trkKLMomVecLAB[1].Rotate(rotAngle, cross);
    trkKLMomVecLAB[1].Rotate(rotAngleKaonX, z_axis);

    trkKL4VecLAB[0].SetPxPyPzE(trkKLMomVecLAB[0][0],
                               trkKLMomVecLAB[0][1],
                               trkKLMomVecLAB[0][2],
                               trkKL[0]->at(3));

    trkKL4VecLAB[1].SetPxPyPzE(trkKLMomVecLAB[1][0],
                               trkKLMomVecLAB[1][1],
                               trkKLMomVecLAB[1][2],
                               trkKL[1]->at(3));

    Obj.lorentz_transf(kaonMomLAB, trkKL4VecLAB[0], trkKL4VecKaonCM[0]);
    Obj.lorentz_transf(kaonMomLAB, trkKL4VecLAB[1], trkKL4VecKaonCM[1]);

    trkKLMomVecKaonCM[0].SetXYZ(trkKL4VecKaonCM[0][0],
                                trkKL4VecKaonCM[0][1],
                                trkKL4VecKaonCM[0][2]);
    trkKLMomVecKaonCM[1].SetXYZ(trkKL4VecKaonCM[1][0],
                                trkKL4VecKaonCM[1][1],
                                trkKL4VecKaonCM[1][2]);

    Double_t
        PiMomMagKaonCM1 = Obj.TwoBodyDecayMass(baseKin.KchrecKLTwoBody[5], mPiCh, mPiCh),
        PiMomMagKaonCM2 = Obj.TwoBodyDecayMass(baseKin.KchrecKLTwoBody[5], mPiCh, mPiCh);

    Double_t angle1 = trkKLMomVecKaonCM[0].Phi(),
             angle2 = trkKLMomVecKaonCM[1].Phi(),
             theta1 = trkKLMomVecKaonCM[0].Theta(),
             theta2 = trkKLMomVecKaonCM[1].Theta();

    gamma = (M_PI_2 - 0.5 * abs(angle1) - 0.5 * abs(angle2));
    gamma_theta = 0.0; //(M_PI_2 - 0.5 * theta1 - 0.5 * theta2);

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

    PiKaon4VecKaonCM[0].SetPxPyPzE(PiMomKaonCM[0][0], PiMomKaonCM[0][1], PiMomKaonCM[0][2], sqrt(pow(PiMomKaonCM[0][0], 2) + pow(PiMomKaonCM[0][1], 2) + pow(PiMomKaonCM[0][2], 2) + pow(mPiCh, 2)));
    PiKaon4VecKaonCM[1].SetPxPyPzE(PiMomKaonCM[1][0], PiMomKaonCM[1][1], PiMomKaonCM[1][2], sqrt(pow(PiMomKaonCM[1][0], 2) + pow(PiMomKaonCM[1][1], 2) + pow(PiMomKaonCM[1][2], 2) + pow(mPiCh, 2)));

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
      trkKLMomVecLAB[0][j] = trkKL[0]->at(j);
      trkKLMomVecLAB[1][j] = trkKL[1]->at(j);
    }

    PiKaon4VecLAB[0].SetPxPyPzE(PiMomKaonLAB[0][0],
                                PiMomKaonLAB[0][1],
                                PiMomKaonLAB[0][2],
                                sqrt(PiMomKaonLAB[0].Mag2() + pow(mPiCh, 2)));
    PiKaon4VecLAB[1].SetPxPyPzE(PiMomKaonLAB[1][0],
                                PiMomKaonLAB[1][1],
                                PiMomKaonLAB[1][2],
                                sqrt(PiMomKaonLAB[1].Mag2() + pow(mPiCh, 2)));

    trkKL4VecLAB[0].SetPxPyPzE(trkKLMomVecLAB[0][0],
                               trkKLMomVecLAB[0][1],
                               trkKLMomVecLAB[0][2],
                               trkKL[0]->at(3));
    trkKL4VecLAB[1].SetPxPyPzE(trkKLMomVecLAB[1][0],
                               trkKLMomVecLAB[1][1],
                               trkKLMomVecLAB[1][2],
                               trkKL[1]->at(3));

    for (Int_t k = 0; k < 4; k++)
    {

      trkKLTwoBody1[k] = PiKaon4VecLAB[0][k];
      trkKLTwoBody2[k] = PiKaon4VecLAB[1][k];

      if (mcflag == 1 && mctruth == 7)
      {
        if (trkKLTwoBody1[k] - trkKLmc[0]->at(k) < trkKLTwoBody2[k] - trkKLmc[0]->at(k))
        {
          histMomPiTwoBody[0][k]->Fill(trkKLTwoBody1[k] - trkKLmc[0]->at(k));
          histMomPiTwoBody[1][k]->Fill(trkKLTwoBody2[k] - trkKLmc[1]->at(k));
        }
        else
        {
          histMomPiTwoBody[0][k]->Fill(trkKLTwoBody1[k] - trkKLmc[1]->at(k));
          histMomPiTwoBody[1][k]->Fill(trkKLTwoBody2[k] - trkKLmc[0]->at(k));
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

        histKLTwoBody[k]->Fill(baseKin.KchrecKLTwoBody[k] - Kchmc->at(k));
        histKL[k]->Fill(KchrecKL->at(k) - Kchmc->at(k));
      }
    }

    // Int_t zmienne
    std::map<std::string, Int_t> intVars = {
        {"nrun", baseKin.nrun},
        {"nev", baseKin.nev},
        {"mcflag", mcflag},
        {"mctruth", mctruth}};

    std::map<std::string, Float_t> floatVars = {
        {"gamma", gamma},
        {"gamma_theta", gamma_theta}};

    // Tablice
    std::map<std::string, std::vector<Int_t>> intArrays = {};

    std::map<std::string, std::vector<Float_t>> floatArrays = {
        {"Kchmc", *Kchmc},
        {"Knemc", *Knemc},
        {"KchrecKS", *KchrecKS},
        {"KchrecKL", *KchrecKL},
        {"KchboostKS", *KchboostKS},
        {"KchboostKL", *KchboostKL},
        {"KchrecKLTwoBody", baseKin.KchrecKLTwoBody},
        {"KchrecKSTwoBody", baseKin.KchrecKSTwoBody},
        {"trk1KS", *trkKS[0]},
        {"trk2KS", *trkKS[1]},
        {"trk1KL", *trkKL[0]},
        {"trk2KL", *trkKL[1]},
        {"trk1KLTwoBody", trkKLTwoBody1},
        {"trk2KLTwoBody", trkKLTwoBody2},
        {"trk1KSTwoBody", trkKSTwoBody1},
        {"trk2KSTwoBody", trkKSTwoBody2},
        {"trk1KSmc", *trkKSmc[0]},
        {"trk2KSmc", *trkKSmc[1]},
        {"trk1KLmc", *trkKLmc[0]},
        {"trk2KLmc", *trkKLmc[1]}};

    writer.Fill(intVars, floatVars, intArrays, floatArrays);

    ++show_progress; // Progress of the loading bar
  }

  TCanvas *canvas = new TCanvas("canvas", "Canvas", 800, 600);

  std::string xTitle[4] = {"p_{x} [MeV/c]", "p_{y} [MeV/c]", "p_{z} [MeV/c]", "E [MeV]"};

  gStyle->SetOptStat(0);

  for (Int_t j = 0; j < 4; j++)
  {
    for (Int_t i = 0; i < 2; i++)
    {
      TLegend *legendPi = new TLegend(0.1, 0.7, 0.3, 0.9);
      legendPi->AddEntry(histMomtrkKL[i][j], Form("Reconstructed charged pion %d", i), "l");
      legendPi->AddEntry(histMomPiTwoBody[i][j], Form("Charged pion %d from 2-body decay", i), "l");
      histMomPiTwoBody[i][j]->SetTitle(Form("Histogram Pi Two Body %d %d", i, j));
      histMomPiTwoBody[i][j]->GetXaxis()->SetTitle(xTitle[j].c_str());
      histMomPiTwoBody[i][j]->GetYaxis()->SetTitle("Counts");
      histMomPiTwoBody[i][j]->GetYaxis()->SetRangeUser(0, histMomPiTwoBody[i][j]->GetMaximum() * 2.0);
      histMomPiTwoBody[i][j]->SetLineColor(kBlue);
      histMomPiTwoBody[i][j]->Draw();
      histMomtrkKL[i][j]->SetLineColor(kRed);
      histMomtrkKL[i][j]->Draw("SAME");
      legendPi->Draw();
      canvas->Print(Form("hist_mom_pi_two_body_%d_%d.png", i, j));
    }

    TLegend *legendKaon = new TLegend(0.1, 0.7, 0.3, 0.9);
    legendKaon->AddEntry(histKL[j], "Reconstructed KL", "l");
    legendKaon->AddEntry(histKLTwoBody[j], "KL from 2-body decay", "l");

    histKLTwoBody[j]->SetTitle(Form("Histogram KL Two Body %d", j));
    histKLTwoBody[j]->GetXaxis()->SetTitle(xTitle[j].c_str());
    histKLTwoBody[j]->GetYaxis()->SetTitle("Counts");
    histKLTwoBody[j]->GetYaxis()->SetRangeUser(0, histKLTwoBody[j]->GetMaximum() * 2.0);
    histKLTwoBody[j]->SetLineColor(kBlue);
    histKLTwoBody[j]->Draw();
    histKL[j]->SetLineColor(kRed);
    histKL[j]->Draw("SAME");
    legendKaon->Draw();
    canvas->Print(Form("hist_kl_two_body_%d.png", j));
  }

  writer.Close();

  return 0;
}