// Author: Szymon Gamrat
// Date of last update: 03.02.2025

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
                                        pow(baseKin.KchrecKLTwoBody[4], 2) );
      baseKin.KchrecKLTwoBody[6] = baseKin.KchrecKL[6];
      baseKin.KchrecKLTwoBody[7] = baseKin.KchrecKL[7];
      baseKin.KchrecKLTwoBody[8] = baseKin.KchrecKL[8];

      
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
        {"trk2Closest", baseKin.trkClosest[1]}};

    writer.Fill(intVars, floatVars, intArrays, floatArrays);
  }

  writer.Close();

  return 0;
}