// Author: Szymon Gamrat
// Date of last update: 03.02.2025

#include <TGraph.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TFitResultPtr.h>
#include <TFitResult.h>
#include <TLegend.h>
#include <THStack.h>

#include <boost/optional.hpp>
#include <SplitFileWriter.h>

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

  chain.SetBranchAddress("nev", &baseKin.nevent);

  // Bhabha variables - avg per run
  chain.SetBranchAddress("Bx", &baseKin.bhabha_vtx[0]);
  chain.SetBranchAddress("By", &baseKin.bhabha_vtx[1]);
  chain.SetBranchAddress("Bz", &baseKin.bhabha_vtx[2]);

  chain.SetBranchAddress("Bpx", &baseKin.phi_mom[0]);
  chain.SetBranchAddress("Bpy", &baseKin.phi_mom[1]);
  chain.SetBranchAddress("Bpz", &baseKin.phi_mom[2]);
  chain.SetBranchAddress("Broots", &baseKin.phi_mom[3]);
  // -----------------------------------------------------------

  // Track - vtx info
  chain.SetBranchAddress("nv", &baseKin.nv);
  chain.SetBranchAddress("ntv", &baseKin.ntv);
  chain.SetBranchAddress("ivOld", baseKin.ivOld);
  // -----------------------------------------------------------

  // Track properties
  chain.SetBranchAddress("CurvOld", baseKin.CurvOld);
  chain.SetBranchAddress("PhivOld", baseKin.PhivOld);
  chain.SetBranchAddress("CotvOld", baseKin.CotvOld);
  // -----------------------------------------------------------

  // Vertex position
  chain.SetBranchAddress("xvOld", baseKin.xvOld);
  chain.SetBranchAddress("yvOld", baseKin.yvOld);
  chain.SetBranchAddress("zvOld", baseKin.zvOld);
  // -----------------------------------------------------------

  // mcflag and mctruth position
  chain.SetBranchAddress("mcflag", &baseKin.mcflag);
  chain.SetBranchAddress("mctruth", &baseKin.mctruth);
  // -----------------------------------------------------------

  // Trk momentum
  baseKin.KchrecKS.resize(9);
  baseKin.KchrecKL.resize(9);
  baseKin.trkKS[0].resize(4);
  baseKin.trkKS[1].resize(4);
  baseKin.trkKL[0].resize(4);
  baseKin.trkKL[1].resize(4);
  chain.SetBranchAddress("Kchrecks", baseKin.KchrecKS.data());
  chain.SetBranchAddress("Kchreckl", baseKin.KchrecKL.data());
  chain.SetBranchAddress("trk1ks", baseKin.trkKS[0].data());
  chain.SetBranchAddress("trk2ks", baseKin.trkKS[1].data());
  chain.SetBranchAddress("trk1kl", baseKin.trkKL[0].data());
  chain.SetBranchAddress("trk2kl", baseKin.trkKL[1].data());
  // -----------------------------------------------------------

  std::string
      base_filename = "KchRec_Control_Sample",
      dirname = (std::string)charged_dir + (std::string)root_files_dir,
      dated_folder = Obj.CreateDatedFolder(dirname);

  SplitFileWriter writer(base_filename, 1.5 * 1024 * 1024 * 1024, false, dated_folder);

  Float_t
      invMass = 0.;

  Int_t nentries = chain.GetEntries();

  Int_t mode = 1; // Model for pi+pi-

  // Initialization of Charged part of decay reconstruction class
  // Constructor is below, in the loop
  boost::optional<KLOE::ChargedVtxRec<Float_t, UChar_t>> eventAnalysis;
  // -------------------------------------------------------------

  baseKin.vtaken.resize(3);
  // baseKin.vtakenKS.resize(3);
  // baseKin.vtakenKL.resize(3);
  baseKin.vtakenClosest.resize(3);

  baseKin.Kchrecnew.resize(9);
  // baseKin.KchrecKS.resize(9);
  // baseKin.KchrecKL.resize(9);
  baseKin.KchrecClosest.resize(9);
  baseKin.KchrecKLTwoBody.resize(9);

  for (Int_t i = 0; i < 2; i++)
  {
    baseKin.trknew[i].resize(4);
    // baseKin.trkKS[i].resize(4);
    // baseKin.trkKL[i].resize(4);
    baseKin.trkClosest[i].resize(4);
  }

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
        histMomPiTwoBody[i].push_back(new TH1F(Form("histMomPiTwoBody_%d_%d", i, j), Form("Histogram Pi Two Body %d %d", i, j), 100, -400, 400));
        histMomtrkKL[i].push_back(new TH1F(Form("histMomtrkKL_%d_%d", i, j), Form("Histogram trk KL %d %d", i, j), 100, -400, 400));
        histKLTwoBody.push_back(new TH1F(Form("histKLTwoBody_%d_%d", i, j), Form("Histogram KL Two Body %d %d", i, j), 100, -600, 600));
        histKL.push_back(new TH1F(Form("histKL_%d_%d", i, j), Form("Histogram KL %d %d", i, j), 100, -600, 600));
      }
      else
      {
        histMomPiTwoBody[i].push_back(new TH1F(Form("histMomPiTwoBody_%d_%d", i, j), Form("Histogram Pi Two Body %d %d", i, j), 100, 100, 400));
        histMomtrkKL[i].push_back(new TH1F(Form("histMomtrkKL_%d_%d", i, j), Form("Histogram trk KL %d %d", i, j), 100, 100, 400));
        histKLTwoBody.push_back(new TH1F(Form("histKLTwoBody_%d_%d", i, j), Form("Histogram KL Two Body %d %d", i, j), 100, 100, 1000));
        histKL.push_back(new TH1F(Form("histKL_%d_%d", i, j), Form("Histogram KL %d %d", i, j), 100, 200, 600));
      }
    }

  Int_t graph_flag = 0;

  Float_t minDiff = 999999.;

  // Progress bar
  boost::progress_display show_progress(nentries);
  // --------------------------------------------------------------------

  std::vector<Float_t>
      trkKLTwoBody1(4),
      trkKLTwoBody2(4);

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

    // Set the proper data type to be analyzed
    Bool_t data_flag = false;
    Int_t mctruth_int = int(baseKin.mctruth), mcflag_int = int(baseKin.mcflag);

    Obj.dataFlagSetter(dataType, data_flag, mcflag_int, mctruth_int);
    // -------------------------------------------------------------------

    baseKin.errFlag = 1;
    baseKin.errFlagKS = 1;
    baseKin.errFlagKL = 1;
    baseKin.errFlagClosest = 1;

    baseKin.vtaken.clear();
    baseKin.vtakenClosest.clear();

    baseKin.Kchrecnew.clear();
    baseKin.KchrecClosest.clear();
    baseKin.KchrecKLTwoBody.clear();

    for (Int_t i = 0; i < 2; i++)
    {
      baseKin.trknew[i].clear();
      baseKin.trkClosest[i].clear();
      baseKin.trkKLTwoBody[i].clear();
    }

    if (data_flag)
    {
      // Construction of the charged rec class object
      eventAnalysis.emplace(baseKin.nv, baseKin.ntv, baseKin.ivOld, baseKin.bhabha_vtx, baseKin.CurvOld, baseKin.PhivOld, baseKin.CotvOld, baseKin.xvOld, baseKin.yvOld, baseKin.zvOld, mode);

      // KMASS HYPOTHESIS
      eventAnalysis->findKchRec(baseKin.Kchrecnew.data(), baseKin.trknew[0].data(), baseKin.trknew[1].data(), baseKin.vtaken.data(), logger);
      // ------------------------------------------------------------------

      // // CLOSEST TO IP
      // eventAnalysis->findKClosestRec(baseKin.KchrecClosest.data(), baseKin.trkClosest[0].data(), baseKin.trkClosest[1].data(), baseKin.vtakenClosest.data(), baseKin.errFlagClosest);
      // // ------------------------------------------------------------------

      if (1)
      {
        // 1. Calculation of pions' momenta from two-body decay using the boost method
        // Using the K1 side of the event
        eventAnalysis->KaonMomFromBoost(baseKin.KchrecKS.data(), baseKin.phi_mom, baseKin.KchboostKS);

        Float_t X_line[3] = {baseKin.KchboostKS[6],
                             baseKin.KchboostKS[7],
                             baseKin.KchboostKS[8]}, // Vertex laying on the line
            p[3] = {baseKin.KchboostKS[0],
                    baseKin.KchboostKS[1],
                    baseKin.KchboostKS[2]}, // Direction of the line
            xB[3] = {baseKin.bhabha_vtx[0],
                     baseKin.bhabha_vtx[1],
                     baseKin.bhabha_vtx[2]}, // Bhabha vertex - laying on the plane
            plane_perp[3] = {0.,
                             baseKin.phi_mom[1],
                             0.}; // Vector perpendicular to the plane from Bhabha momentum

        // Corrected IP event by event
        eventAnalysis->IPBoostCorr(X_line, p, xB, plane_perp, baseKin.ip);

        baseKin.ip[0] = baseKin.bhabha_vtx[0];
        baseKin.ip[1] = baseKin.bhabha_vtx[1];
        // z coordinate of the IP is set to the Bhabha vertex z coordinate if it differs by more than 2 cm
        if (abs(baseKin.ip[2] - baseKin.bhabha_vtx[2]) > 2.0)
          baseKin.ip[2] = baseKin.bhabha_vtx[2];

        // 2. Calculation of KL flight direction
        Double_t KLpath = sqrt(pow(baseKin.KchrecKL[6] - baseKin.ip[0], 2) +
                               pow(baseKin.KchrecKL[7] - baseKin.ip[1], 2) +
                               pow(baseKin.KchrecKL[8] - baseKin.ip[2], 2));
        TVector3 KLflightDirection = {(baseKin.KchrecKL[6] - baseKin.ip[0]) / KLpath,
                                      (baseKin.KchrecKL[7] - baseKin.ip[1]) / KLpath,
                                      (baseKin.KchrecKL[8] - baseKin.ip[2]) / KLpath};

        eventAnalysis->KaonMomFromBoost(baseKin.KchrecKL.data(),
                                        baseKin.phi_mom,
                                        baseKin.KchboostKL);

        // 2.1 Calculation of KL momentum magnitude
        Double_t KLmomMag = sqrt(pow(baseKin.KchboostKL[0], 2) +
                                 pow(baseKin.KchboostKL[1], 2) +
                                 pow(baseKin.KchboostKL[2], 2));

        // 2.2 Calculation of KL momentum from 2 body decay
        for (Int_t j = 0; j < 3; j++)
        {
          // Components
          baseKin.KchrecKLTwoBody[j] = KLflightDirection[j] * KLmomMag;
        }
        // Energy
        baseKin.KchrecKLTwoBody[3] = baseKin.KchboostKL[3];

        // Calculate the magnitude of the KL momentum
        baseKin.KchrecKLTwoBody[4] = sqrt(pow(baseKin.KchrecKLTwoBody[0], 2) +
                                          pow(baseKin.KchrecKLTwoBody[1], 2) +
                                          pow(baseKin.KchrecKLTwoBody[2], 2));
        // Calculate the invariant mass of the KL
        baseKin.KchrecKLTwoBody[5] = sqrt(pow(baseKin.KchrecKLTwoBody[3], 2) -
                                          pow(baseKin.KchrecKLTwoBody[4], 2));

        // Copy the rest of the KL vertex components from the original KL vertex
        baseKin.KchrecKLTwoBody[6] = baseKin.KchrecKL[6];
        baseKin.KchrecKLTwoBody[7] = baseKin.KchrecKL[7];
        baseKin.KchrecKLTwoBody[8] = baseKin.KchrecKL[8];
        // ------------------------------------------------------------------

        // 2. Go to KL CM frame
        TVector3
            PiMomKaonCM[2],           // Pion momenta in the KL CM frame
            PiMomKaonLAB[2],          // Pion momenta in the LAB frame
            x_axis = {1.0, 0.0, 0.0}, // x-axis for rotation
            kaonMomLAB = {baseKin.KchrecKLTwoBody[0] / baseKin.KchrecKLTwoBody[4],
                          baseKin.KchrecKLTwoBody[1] / baseKin.KchrecKLTwoBody[4],
                          baseKin.KchrecKLTwoBody[2] / baseKin.KchrecKLTwoBody[4]},
            cross = kaonMomLAB.Cross(x_axis),
            trkKLMomVecLAB[2], // Track momenta in the LAB frame
            trkKLMomVecKaonCM[2];

        TLorentzVector
            kaon4VecLAB,
            kaon4VecKaonCM,
            trkKL4VecLAB[2],
            trkKL4VecKaonCM[2],
            PiKaon4VecKaonCM[2],
            PiKaon4VecLAB[2];

        Double_t rotAngle = kaonMomLAB.Angle(x_axis);

        kaonMomLAB.Rotate(rotAngle, cross);

        TVector3
          boost_KaonTwoBody = { -baseKin.KchrecKLTwoBody[0] / baseKin.KchrecKLTwoBody[3],
                                -baseKin.KchrecKLTwoBody[1] / baseKin.KchrecKLTwoBody[3],
                                -baseKin.KchrecKLTwoBody[2] / baseKin.KchrecKLTwoBody[3]};
            // boost_KaonTwoBody = {-baseKin.KchrecKLTwoBody[4] / baseKin.KchrecKLTwoBody[3],
            //                      0.0,
            //                      0.0};

        for (Int_t j = 0; j < 3; j++)
        {
          // Components of the track momenta in the LAB frame
          trkKLMomVecLAB[0][j] = baseKin.trkKL[0][j];
          trkKLMomVecLAB[1][j] = baseKin.trkKL[1][j];
        }

        // std::cout << "Momenta before rotation: " << std::endl;
        // trkKLMomVecLAB[0].Print();
        // trkKLMomVecLAB[1].Print();

        // trkKLMomVecLAB[0].Rotate(rotAngle, cross);
        // trkKLMomVecLAB[1].Rotate(rotAngle, cross);

        // std::cout << "Momenta after rotation: " << std::endl;
        // trkKLMomVecLAB[0].Print();
        // trkKLMomVecLAB[1].Print();

        trkKL4VecLAB[0].SetPxPyPzE(trkKLMomVecLAB[0][0],
                                   trkKLMomVecLAB[0][1],
                                   trkKLMomVecLAB[0][2],
                                   baseKin.trkKL[0][3]);

        trkKL4VecLAB[1].SetPxPyPzE(trkKLMomVecLAB[1][0],
                                   trkKLMomVecLAB[1][1],
                                   trkKLMomVecLAB[1][2],
                                   baseKin.trkKL[1][3]);

        Obj.lorentz_transf(boost_KaonTwoBody, trkKL4VecLAB[0], trkKL4VecKaonCM[0]);
        Obj.lorentz_transf(boost_KaonTwoBody, trkKL4VecLAB[1], trkKL4VecKaonCM[1]);

        trkKLMomVecKaonCM[0].SetXYZ(trkKL4VecKaonCM[0][0],
                                    trkKL4VecKaonCM[0][1],
                                    trkKL4VecKaonCM[0][2]);
        trkKLMomVecKaonCM[1].SetXYZ(trkKL4VecKaonCM[1][0],
                                    trkKL4VecKaonCM[1][1],
                                    trkKL4VecKaonCM[1][2]);

        Double_t  
                  PiMomMagKaonCM1 = Obj.TwoBodyDecayMass(baseKin.KchrecKLTwoBody[5],  mPiCh, mPiCh),
                  PiMomMagKaonCM2 = Obj.TwoBodyDecayMass(baseKin.KchrecKLTwoBody[5], mPiCh, mPiCh);

        Double_t angle1 = trkKLMomVecKaonCM[0].Phi(),
                 angle2 = trkKLMomVecKaonCM[1].Phi(),
                 theta1 = trkKLMomVecKaonCM[0].Theta(),
                 theta2 = trkKLMomVecKaonCM[1].Theta();

        gamma = (M_PI_2 - 0.5 * abs(angle1) - 0.5 * abs(angle2));
        gamma_theta = (M_PI_2 - 0.5 * theta1 - 0.5 * theta2);

        if (angle1 < 0.0)
          PiMomKaonCM[0].SetXYZ(  PiMomMagKaonCM1 * sin(theta1 + gamma_theta) * cos(angle1 - gamma),
                                  PiMomMagKaonCM1 * sin(theta1 + gamma_theta) * sin(angle1 - gamma),
                                  PiMomMagKaonCM1 * cos(theta1 + gamma_theta) );
        else
          PiMomKaonCM[0].SetXYZ(  PiMomMagKaonCM1 * sin(theta1 + gamma_theta) * cos(angle1 + gamma),
                                  PiMomMagKaonCM1 * sin(theta1 + gamma_theta) * sin(angle1 + gamma),
                                  PiMomMagKaonCM1 * cos(theta1 + gamma_theta) );
        if (angle2 < 0.0)
          PiMomKaonCM[1].SetXYZ(  PiMomMagKaonCM1 * sin(theta2 + gamma_theta) * cos(angle2 - gamma),
                                  PiMomMagKaonCM1 * sin(theta2 + gamma_theta) * sin(angle2 - gamma),
                                  PiMomMagKaonCM1 * cos(theta2 + gamma_theta) );
        else
          PiMomKaonCM[1].SetXYZ(  PiMomMagKaonCM1 * sin(theta2 + gamma_theta) * cos(angle2 + gamma),
                                  PiMomMagKaonCM1 * sin(theta2 + gamma_theta) * sin(angle2 + gamma),
                                  PiMomMagKaonCM1 * cos(theta2 + gamma_theta) );

        PiKaon4VecKaonCM[0].SetPxPyPzE(PiMomKaonCM[0][0], PiMomKaonCM[0][1], PiMomKaonCM[0][2], sqrt(pow(PiMomKaonCM[0][0], 2) + pow(PiMomKaonCM[0][1], 2) + pow(PiMomKaonCM[0][2], 2) + pow(mPiCh, 2)));
        PiKaon4VecKaonCM[1].SetPxPyPzE(PiMomKaonCM[1][0], PiMomKaonCM[1][1], PiMomKaonCM[1][2], sqrt(pow(PiMomKaonCM[1][0], 2) + pow(PiMomKaonCM[1][1], 2) + pow(PiMomKaonCM[1][2], 2) + pow(mPiCh, 2)));

        boost_KaonTwoBody = -boost_KaonTwoBody; // Invert the boost direction for the pions
        Obj.lorentz_transf(boost_KaonTwoBody, PiKaon4VecKaonCM[0], PiKaon4VecLAB[0]);
        Obj.lorentz_transf(boost_KaonTwoBody, PiKaon4VecKaonCM[1], PiKaon4VecLAB[1]);

        PiMomKaonLAB[0].SetXYZ(PiKaon4VecLAB[0][0],
                               PiKaon4VecLAB[0][1],
                               PiKaon4VecLAB[0][2]);
        PiMomKaonLAB[1].SetXYZ(PiKaon4VecLAB[1][0],
                               PiKaon4VecLAB[1][1],
                               PiKaon4VecLAB[1][2]);

        // PiMomKaonLAB[0].Rotate(-rotAngle, cross);
        // PiMomKaonLAB[1].Rotate(-rotAngle, cross);

        for (Int_t j = 0; j < 3; j++)
        {
          // Components of the track momenta in the LAB frame
          trkKLMomVecLAB[0][j] = baseKin.trkKL[0][j];
          trkKLMomVecLAB[1][j] = baseKin.trkKL[1][j];
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
                                   baseKin.trkKL[0][3]);
        trkKL4VecLAB[1].SetPxPyPzE(trkKLMomVecLAB[1][0],
                                   trkKLMomVecLAB[1][1],
                                   trkKLMomVecLAB[1][2],
                                   baseKin.trkKL[1][3]);

        for (Int_t k = 0; k < 4; k++)
        {

          trkKLTwoBody1[k] = PiKaon4VecLAB[0][k];
          trkKLTwoBody2[k] = PiKaon4VecLAB[1][k];

          if (mcflag_int == 0)
          {
            histMomPiTwoBody[0][k]->Fill(trkKLTwoBody1[k]);
            histMomPiTwoBody[1][k]->Fill(trkKLTwoBody2[k]);
            histMomtrkKL[0][k]->Fill(trkKL4VecLAB[0][k]);
            histMomtrkKL[1][k]->Fill(trkKL4VecLAB[1][k]);

            histKLTwoBody[k]->Fill(baseKin.KchrecKLTwoBody[k]);
            histKL[k]->Fill(baseKin.KchrecKL[k]);
          }
        }
      }

      // Int_t zmienne
      std::map<std::string, Int_t> intVars = {
          {"nrun", 0},
          {"nev", baseKin.nevent},
          {"mcflag", mcflag_int},
          {"mctruth", mctruth_int},
          {"errflag", baseKin.errFlag},
          {"errflagks", baseKin.errFlagKS},
          {"errflagkl", baseKin.errFlagKL},
          {"errflagclosest", baseKin.errFlagClosest}};

      std::map<std::string, Float_t> floatVars = {
          {"gamma", gamma}};

      // Tablice
      std::map<std::string, std::vector<Int_t>> intArrays = {
          {"vtaken", baseKin.vtaken},
          {"vtakenKS", baseKin.vtakenKS},
          {"vtakenKL", baseKin.vtakenKL},
          {"vtakenClosest", baseKin.vtakenClosest}};

      std::map<std::string, std::vector<Float_t>> floatArrays = {
          {"Kchrec", baseKin.Kchrecnew},
          {"KchrecKS", baseKin.KchrecKS},
          {"KchrecKL", baseKin.KchrecKL},
          {"KchrecClosest", baseKin.KchrecClosest},
          {"KchrecKLTwoBody", baseKin.KchrecKLTwoBody},
          {"trk1", baseKin.trknew[0]},
          {"trk2", baseKin.trknew[1]},
          {"trk1KS", baseKin.trkKS[0]},
          {"trk2KS", baseKin.trkKS[1]},
          {"trk1KL", baseKin.trkKL[0]},
          {"trk2KL", baseKin.trkKL[1]},
          {"trk1Closest", baseKin.trkClosest[0]},
          {"trk2Closest", baseKin.trkClosest[1]},
          {"trk1TwoBody", trkKLTwoBody1},
          {"trk2TwoBody", trkKLTwoBody2},
      };
      writer.Fill(intVars, floatVars, intArrays, floatArrays);
    }

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