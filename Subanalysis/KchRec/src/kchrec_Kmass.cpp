// Author: Szymon Gamrat
// Date of last update: 03.02.2025

#include <TGraph.h>
#include <TCanvas.h>
#include <TF1.h>

#include <boost/optional.hpp>
#include <SplitFileWriter.h>

#include "../inc/kchrec.hpp"

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
  chain.SetBranchAddress("iv", baseKin.iv);
  // -----------------------------------------------------------

  // Track properties
  chain.SetBranchAddress("Curv", baseKin.Curv);
  chain.SetBranchAddress("Phiv", baseKin.Phiv);
  chain.SetBranchAddress("Cotv", baseKin.Cotv);
  // -----------------------------------------------------------

  // Vertex position
  chain.SetBranchAddress("xv", baseKin.xv);
  chain.SetBranchAddress("yv", baseKin.yv);
  chain.SetBranchAddress("zv", baseKin.zv);
  // -----------------------------------------------------------

  // mcflag and mctruth position
  chain.SetBranchAddress("mcflag", &baseKin.mcflag);
  chain.SetBranchAddress("mctruth", &baseKin.mctruth);
  // -----------------------------------------------------------

  std::string
      base_filename = "KchRec_Control_Sample",
      dirname = (std::string)charged_dir + (std::string)root_files_dir;

  SplitFileWriter writer(base_filename, 1.5 * 1024 * 1024 * 1024, false, dirname);

  Float_t
      invMass = 0.;

  Int_t nentries = chain.GetEntries();

  TGraph *angleGraph;

  Int_t dupa = 0;

  Int_t mode = 1; // Model for pi+pi-

  // Initialization of Charged part of decay reconstruction class
  // Constructor is below, in the loop
  boost::optional<KLOE::ChargedVtxRec<>> eventAnalysis;
  // -------------------------------------------------------------

  baseKin.vtaken.resize(3);
  baseKin.vtakenKS.resize(3);
  baseKin.vtakenKL.resize(3);
  baseKin.vtakenClosest.resize(3);

  baseKin.Kchrecnew.resize(9);
  baseKin.KchrecKS.resize(9);
  baseKin.KchrecKL.resize(9);
  baseKin.KchrecClosest.resize(9);

  for (Int_t i = 0; i < 2; i++)
  {
    baseKin.trknew[i].resize(4);
    baseKin.trkKS[i].resize(4);
    baseKin.trkKL[i].resize(4);
    baseKin.trkClosest[i].resize(4);
  }

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
    baseKin.vtakenKS.clear();
    baseKin.vtakenKL.clear();
    baseKin.vtakenClosest.clear();

    baseKin.Kchrecnew.clear();
    baseKin.KchrecKS.clear();
    baseKin.KchrecKL.clear();
    baseKin.KchrecClosest.clear();
    baseKin.KchrecKLTwoBody.clear();

    for (Int_t i = 0; i < 2; i++)
    {
      baseKin.trknew[i].clear();
      baseKin.trkKS[i].clear();
      baseKin.trkKL[i].clear();
      baseKin.trkClosest[i].clear();
      baseKin.trkKLTwoBody[i].clear();
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
      baseKin.trkKLTwoBody[i].resize(4);
    }

    if (data_flag)
    {
      // Construction of the charged rec class object
      eventAnalysis.emplace(baseKin.nv, baseKin.ntv, baseKin.iv, baseKin.bhabha_vtx, baseKin.Curv, baseKin.Phiv, baseKin.Cotv, baseKin.xv, baseKin.yv, baseKin.zv, mode);

      // KMASS HYPOTHESIS
      eventAnalysis->findKchRec(baseKin.Kchrecnew.data(), baseKin.trknew[0].data(), baseKin.trknew[1].data(), baseKin.vtaken.data(), baseKin.errFlag);
      // ------------------------------------------------------------------

      // KSL HYPOTHESIS
      eventAnalysis->findKSLRec(16, -1, baseKin.KchrecKS.data(), baseKin.trkKS[0].data(), baseKin.trkKS[1].data(), baseKin.vtakenKS.data(), baseKin.errFlagKS);

      eventAnalysis->findKSLRec(10, baseKin.vtakenKS[0], baseKin.KchrecKL.data(), baseKin.trkKL[0].data(), baseKin.trkKL[1].data(), baseKin.vtakenKL.data(), baseKin.errFlagKL);
      // ------------------------------------------------------------------

      // CLOSEST TO IP
      eventAnalysis->findKClosestRec(baseKin.KchrecClosest.data(), baseKin.trkClosest[0].data(), baseKin.trkClosest[1].data(), baseKin.vtakenClosest.data(), baseKin.errFlagClosest);
      // ------------------------------------------------------------------

      // CALCULATION OF PIONS' MOMENTA FROM TWO BODY DECAY
      for (Int_t j = 0; j < 4; j++)
      {
        baseKin.KchrecKLTwoBody[j] = baseKin.phi_mom[j] - baseKin.KchrecKS[j];
      }
      baseKin.KchrecKLTwoBody[4] = sqrt(pow(baseKin.KchrecKLTwoBody[0], 2) +
                                        pow(baseKin.KchrecKLTwoBody[1], 2) +
                                        pow(baseKin.KchrecKLTwoBody[2], 2));
      baseKin.KchrecKLTwoBody[5] = sqrt(pow(baseKin.KchrecKLTwoBody[3], 2) -
                                        pow(baseKin.KchrecKLTwoBody[4], 2));
      baseKin.KchrecKLTwoBody[6] = baseKin.KchrecKL[6];
      baseKin.KchrecKLTwoBody[7] = baseKin.KchrecKL[7];
      baseKin.KchrecKLTwoBody[8] = baseKin.KchrecKL[8];

      // 1. Calculate the axis of kaon
      TVector3
          boost_Kaon = {-baseKin.KchrecKL[4] / baseKin.KchrecKL[3],
                        0.0,
                        0.0},
          boost_KaonTwoBody = {-baseKin.KchrecKLTwoBody[4] / baseKin.KchrecKLTwoBody[3],
                               0.0,
                               0.0},
          trkKLMomVecLAB[2],
          trkKLMomVecKaonCM[2];

      TLorentzVector
          kaon4VecLAB,
          kaon4VecKaonCM,
          trkKL4VecLAB[2],
          trkKL4VecKaonCM[2],
          PiKaon4VecKaonCM[2],
          PiKaon4VecLAB[2];

      kaon4VecLAB.SetPxPyPzE(baseKin.KchrecKL[0],
                             baseKin.KchrecKL[1],
                             baseKin.KchrecKL[2],
                             baseKin.KchrecKL[3]);

      Double_t angleKaonTrk,
          magTrkLAB;

      for (Int_t j = 0; j < 2; j++)
      {
        trkKLMomVecLAB[j].SetXYZ(baseKin.trkKL[j][0],
                                 baseKin.trkKL[j][1],
                                 baseKin.trkKL[j][2]);

        angleKaonTrk = kaon4VecLAB.Angle(trkKLMomVecLAB[j]);
        magTrkLAB = trkKLMomVecLAB[j].Mag();

        if (j == 0)
          trkKL4VecLAB[j].SetPxPyPzE(magTrkLAB * cos(angleKaonTrk),
                                     magTrkLAB * sin(angleKaonTrk),
                                     0.0,
                                     baseKin.trkKL[j][3]);
        else
          trkKL4VecLAB[j].SetPxPyPzE(magTrkLAB * cos(2 * M_PI - angleKaonTrk),
                                     magTrkLAB * sin(2 * M_PI - angleKaonTrk),
                                     0.0,
                                     baseKin.trkKL[j][3]);

        Obj.lorentz_transf(boost_Kaon, trkKL4VecLAB[j], trkKL4VecKaonCM[j]);

        trkKLMomVecKaonCM[j].SetXYZ(trkKL4VecKaonCM[j][0],
                                    trkKL4VecKaonCM[j][1],
                                    trkKL4VecKaonCM[j][2]);
      }

      // 2. Calculation of momentum magnitude using two body decay (based on theory)

      Double_t
          PiMomMagKaonCM = 0.5 * sqrt(pow(baseKin.KchrecKLTwoBody[5], 2) - 4. * pow(mPiCh, 2));

      // 3. Check for what angle accordance is the best

      TVector3
          PiMomKaonCM[2],
          PiMomKaonLAB[2];

      TF1 *func = new TF1("minFunc", KLOE::pm00::MomMinAngleFunction, 0., M_PI, 5, 1);

      func->SetParameter(0, trkKLMomVecKaonCM[0][0]);
      func->SetParameter(1, trkKLMomVecKaonCM[0][1]);
      func->SetParameter(2, trkKLMomVecKaonCM[1][0]);
      func->SetParameter(3, trkKLMomVecKaonCM[1][1]);
      func->SetParameter(4, PiMomMagKaonCM);

      TVector3
          x_axis = {1.0, 0.0, 0.0},
          kaonMomLAB = {baseKin.KchrecKLTwoBody[0] / baseKin.KchrecKLTwoBody[4],
                        baseKin.KchrecKLTwoBody[1] / baseKin.KchrecKLTwoBody[4],
                        baseKin.KchrecKLTwoBody[2] / baseKin.KchrecKLTwoBody[4]},
          cross = kaonMomLAB.Cross(x_axis);

      Double_t
          bestAngle = func->GetMinimumX(0.0, M_PI),
          rotAngle = kaonMomLAB.Angle(x_axis);

      // TMatrixD K(3, 3);

      // K(0, 0) = 0;
      // K(0, 1) = -cross[2];
      // K(0, 2) = cross[1];
      // K(1, 0) = cross[2];
      // K(1, 1) = 0;
      // K(1, 2) = -cross[0];
      // K(2, 0) = -cross[1];
      // K(2, 1) = cross[0];
      // K(2, 2) = 0;

      // TMatrixD
      //     K2 = K * K,
      //     I(TMatrixD::kUnit, K),
      //     R = I;

      // R += sin(rotAngle) * K;
      // R += (1 - cos(rotAngle)) * K2;

      // R = R.Invert();

      PiMomKaonCM[0].SetXYZ(PiMomMagKaonCM * cos(bestAngle),
                            PiMomMagKaonCM * sin(bestAngle),
                            0.0);
      PiMomKaonCM[1].SetXYZ(PiMomMagKaonCM * cos(M_PI + bestAngle),
                            PiMomMagKaonCM * sin(M_PI + bestAngle),
                            0.0);

      PiKaon4VecKaonCM[0].SetPxPyPzE(PiMomKaonCM[0][0], PiMomKaonCM[0][1], PiMomKaonCM[0][2], sqrt(pow(PiMomMagKaonCM, 2) + pow(mPiCh, 2)));
      PiKaon4VecKaonCM[1].SetPxPyPzE(PiMomKaonCM[1][0], PiMomKaonCM[1][1], PiMomKaonCM[1][2], sqrt(pow(PiMomMagKaonCM, 2) + pow(mPiCh, 2)));

      boost_KaonTwoBody = -boost_KaonTwoBody;

      Obj.lorentz_transf(boost_KaonTwoBody, PiKaon4VecKaonCM[0], PiKaon4VecLAB[0]);
      Obj.lorentz_transf(boost_KaonTwoBody, PiKaon4VecKaonCM[1], PiKaon4VecLAB[1]);

      PiMomKaonLAB[0].SetXYZ(PiKaon4VecLAB[0][0],
                             PiKaon4VecLAB[0][1],
                             PiKaon4VecLAB[0][2]);
      PiMomKaonLAB[1].SetXYZ(PiKaon4VecLAB[1][0],
                             PiKaon4VecLAB[1][1],
                             PiKaon4VecLAB[1][2]);

      PiMomKaonLAB[0].Rotate(-rotAngle, cross);
      PiMomKaonLAB[1].Rotate(-rotAngle, cross);

      PiKaon4VecLAB[0].SetPxPyPzE(PiMomKaonLAB[0][0],
                                  PiMomKaonLAB[0][1],
                                  PiMomKaonLAB[0][2],
                                  PiKaon4VecLAB[0][3]);
      PiKaon4VecLAB[1].SetPxPyPzE(PiMomKaonLAB[1][0],
                                  PiMomKaonLAB[1][1],
                                  PiMomKaonLAB[1][2],
                                  PiKaon4VecLAB[1][3]);

      trkKL4VecLAB[0].SetPxPyPzE(baseKin.trkKL[0][0],
                                 baseKin.trkKL[0][1],
                                 baseKin.trkKL[0][2],
                                 baseKin.trkKL[0][3]);

      trkKL4VecLAB[1].SetPxPyPzE(baseKin.trkKL[1][0],
                                 baseKin.trkKL[1][1],
                                 baseKin.trkKL[1][2],
                                 baseKin.trkKL[1][3]);

      std::vector<Float_t>
          trkKLTwoBody1(4),
          trkKLTwoBody2(4);

      Double_t
          metric1 = sqrt((PiKaon4VecLAB[0] - trkKL4VecLAB[0]).Mag2() + (PiKaon4VecLAB[1] - trkKL4VecLAB[1]).Mag2()),
          metric2 = sqrt((PiKaon4VecLAB[0] - trkKL4VecLAB[1]).Mag2() + (PiKaon4VecLAB[1] - trkKL4VecLAB[0]).Mag2());

      for (Int_t k = 0; k < 4; k++)
      {
        if (metric1 < metric2)
        {
          trkKLTwoBody1[k] = PiKaon4VecLAB[0][k];
          trkKLTwoBody2[k] = PiKaon4VecLAB[1][k];
        }
        else
        {
          trkKLTwoBody1[k] = PiKaon4VecLAB[1][k];
          trkKLTwoBody2[k] = PiKaon4VecLAB[0][k];
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

      std::map<std::string, Float_t> floatVars;

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
          {"trk2TwoBody", trkKLTwoBody2}};
          writer.Fill(intVars, floatVars, intArrays, floatArrays);
    }
  }

  writer.Close();

  return 0;
}