//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Oct  1 19:16:57 2025 by ROOT version 6.24/07
// from TChain h1/
//////////////////////////////////////////////////////////

#ifndef MC_fit_comparison_h
#define MC_fit_comparison_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

// Headers needed by this particular selector
#include <vector>

class MC_fit_comparison : public TSelector
{
public:
  TTreeReader fReader; //! the tree reader
  TTree *fChain = 0;   //! pointer to the analyzed TTree or TChain

  // Readers to access the data (delete the ones you do not need).
  TTreeReaderValue<Int_t> Eclfilfo = {fReader, "Eclfilfo"};
  TTreeReaderValue<Int_t> Eclfilfoword = {fReader, "Eclfilfoword"};
  TTreeReaderValue<Int_t> bunchnum = {fReader, "bunchnum"};
  TTreeReaderValue<Int_t> mcflag = {fReader, "mcflag"};
  TTreeReaderValue<Int_t> mctruth = {fReader, "mctruth"};
  TTreeReaderValue<Int_t> nclu = {fReader, "nclu"};
  TTreeReaderValue<Int_t> necls = {fReader, "necls"};
  TTreeReaderValue<Int_t> nev = {fReader, "nev"};
  TTreeReaderValue<Int_t> nrun = {fReader, "nrun"};
  TTreeReaderValue<Int_t> ntcl = {fReader, "ntcl"};
  TTreeReaderValue<Int_t> ntmc = {fReader, "ntmc"};
  TTreeReaderValue<Int_t> ntv = {fReader, "ntv"};
  TTreeReaderValue<Int_t> nv = {fReader, "nv"};
  TTreeReaderValue<Int_t> nvtxmc = {fReader, "nvtxmc"};
  TTreeReaderValue<Int_t> errorcode = {fReader, "errorcode"};
  TTreeReaderValue<Int_t> goodClustersTriKinFitSize = {fReader, "goodClustersTriKinFitSize"};
  TTreeReaderValue<Float_t> Bpx = {fReader, "Bpx"};
  TTreeReaderValue<Float_t> Bpy = {fReader, "Bpy"};
  TTreeReaderValue<Float_t> Bpz = {fReader, "Bpz"};
  TTreeReaderValue<Float_t> Broots = {fReader, "Broots"};
  TTreeReaderValue<Float_t> Bx = {fReader, "Bx"};
  TTreeReaderValue<Float_t> By = {fReader, "By"};
  TTreeReaderValue<Float_t> Bz = {fReader, "Bz"};
  TTreeReaderValue<Float_t> Chi2SignalKinFit = {fReader, "Chi2SignalKinFit"};
  TTreeReaderValue<Float_t> Chi2TriKinFit = {fReader, "Chi2TriKinFit"};
  TTreeReaderValue<Float_t> Cotv1 = {fReader, "Cotv1"};
  TTreeReaderValue<Float_t> Cotv2 = {fReader, "Cotv2"};
  TTreeReaderValue<Float_t> CotvSmeared1 = {fReader, "CotvSmeared1"};
  TTreeReaderValue<Float_t> CotvSmeared2 = {fReader, "CotvSmeared2"};
  TTreeReaderValue<Float_t> Curv1 = {fReader, "Curv1"};
  TTreeReaderValue<Float_t> Curv2 = {fReader, "Curv2"};
  TTreeReaderValue<Float_t> CurvSmeared1 = {fReader, "CurvSmeared1"};
  TTreeReaderValue<Float_t> CurvSmeared2 = {fReader, "CurvSmeared2"};
  TTreeReaderValue<Float_t> KaonChTimeCMBoostLor = {fReader, "KaonChTimeCMBoostLor"};
  TTreeReaderValue<Float_t> KaonChTimeCMMC = {fReader, "KaonChTimeCMMC"};
  TTreeReaderValue<Float_t> KaonChTimeLABBoostLor = {fReader, "KaonChTimeLABBoostLor"};
  TTreeReaderValue<Float_t> KaonChTimeLABMC = {fReader, "KaonChTimeLABMC"};
  TTreeReaderValue<Float_t> KaonNeTimeCMBoostLor = {fReader, "KaonNeTimeCMBoostLor"};
  TTreeReaderValue<Float_t> KaonNeTimeCMMC = {fReader, "KaonNeTimeCMMC"};
  TTreeReaderValue<Float_t> KaonNeTimeLABBoostLor = {fReader, "KaonNeTimeLABBoostLor"};
  TTreeReaderValue<Float_t> KaonNeTimeLABMC = {fReader, "KaonNeTimeLABMC"};
  TTreeReaderValue<Float_t> Phiv1 = {fReader, "Phiv1"};
  TTreeReaderValue<Float_t> Phiv2 = {fReader, "Phiv2"};
  TTreeReaderValue<Float_t> PhivSmeared1 = {fReader, "PhivSmeared1"};
  TTreeReaderValue<Float_t> PhivSmeared2 = {fReader, "PhivSmeared2"};
  TTreeReaderValue<Float_t> Qmiss = {fReader, "Qmiss"};
  TTreeReaderValue<Float_t> T0step1 = {fReader, "T0step1"};
  TTreeReaderValue<Float_t> TrcSum = {fReader, "TrcSum"};
  TTreeReaderValue<Float_t> minv4gam = {fReader, "minv4gam"};
  TTreeReaderArray<int> Asscl = {fReader, "Asscl"};
  TTreeReaderArray<int> cutsApplied = {fReader, "cutsApplied"};
  TTreeReaderArray<int> eclstream = {fReader, "eclstream"};
  TTreeReaderArray<int> g4takenTriKinFit = {fReader, "g4takenTriKinFit"};
  TTreeReaderArray<int> iv = {fReader, "iv"};
  TTreeReaderArray<int> mother = {fReader, "mother"};
  TTreeReaderArray<int> pidmc = {fReader, "pidmc"};
  TTreeReaderArray<int> vtaken = {fReader, "vtaken"};
  TTreeReaderArray<int> vtakenClosest = {fReader, "vtakenClosest"};
  TTreeReaderArray<int> vtxmc = {fReader, "vtxmc"};
  TTreeReaderArray<int> goodClustersTriKinFit = {fReader, "goodClustersTriKinFit"};
  TTreeReaderArray<float> Cotv = {fReader, "Cotv"};
  TTreeReaderArray<float> CotvMC = {fReader, "CotvMC"};
  TTreeReaderArray<float> Curv = {fReader, "Curv"};
  TTreeReaderArray<float> CurvMC = {fReader, "CurvMC"};
  TTreeReaderArray<float> Enecl = {fReader, "Enecl"};
  TTreeReaderArray<float> Kchboost = {fReader, "Kchboost"};
  TTreeReaderArray<float> KchboostFit = {fReader, "KchboostFit"};
  TTreeReaderArray<float> KchboostKL = {fReader, "KchboostKL"};
  TTreeReaderArray<float> KchboostKS = {fReader, "KchboostKS"};
  TTreeReaderArray<float> Kchmc = {fReader, "Kchmc"};
  TTreeReaderArray<float> Kchrec = {fReader, "Kchrec"};
  TTreeReaderArray<float> KchrecClosest = {fReader, "KchrecClosest"};
  TTreeReaderArray<float> KchrecFit = {fReader, "KchrecFit"};
  TTreeReaderArray<float> KchrecKL = {fReader, "KchrecKL"};
  TTreeReaderArray<float> KchrecKS = {fReader, "KchrecKS"};
  TTreeReaderArray<float> Knerec = {fReader, "Knerec"};
  TTreeReaderArray<float> Knemc = {fReader, "Knemc"};
  TTreeReaderArray<float> KnerecFit = {fReader, "KnerecFit"};
  TTreeReaderArray<float> KnereclorFit = {fReader, "KnereclorFit"};
  TTreeReaderArray<float> KnetriKinFit = {fReader, "KnetriKinFit"};
  TTreeReaderArray<float> Phiv = {fReader, "Phiv"};
  TTreeReaderArray<float> PhivMC = {fReader, "PhivMC"};
  TTreeReaderArray<float> Tcl = {fReader, "Tcl"};
  TTreeReaderArray<float> Xcl = {fReader, "Xcl"};
  TTreeReaderArray<float> Ycl = {fReader, "Ycl"};
  TTreeReaderArray<float> Zcl = {fReader, "Zcl"};
  TTreeReaderArray<float> gammaMomTriKinFit1 = {fReader, "gammaMomTriKinFit1"};
  TTreeReaderArray<float> gammaMomTriKinFit2 = {fReader, "gammaMomTriKinFit2"};
  TTreeReaderArray<float> gammaMomTriKinFit3 = {fReader, "gammaMomTriKinFit3"};
  TTreeReaderArray<float> gammaMomTriKinFit4 = {fReader, "gammaMomTriKinFit4"};
  TTreeReaderArray<float> gammaMomTriangle1 = {fReader, "gammaMomTriangle1"};
  TTreeReaderArray<float> gammaMomTriangle2 = {fReader, "gammaMomTriangle2"};
  TTreeReaderArray<float> gammaMomTriangle3 = {fReader, "gammaMomTriangle3"};
  TTreeReaderArray<float> gammaMomTriangle4 = {fReader, "gammaMomTriangle4"};
  TTreeReaderArray<float> ip = {fReader, "ip"};
  TTreeReaderArray<float> ipFit = {fReader, "ipFit"};
  TTreeReaderArray<float> ipKL = {fReader, "ipKL"};
  TTreeReaderArray<float> ipKS = {fReader, "ipKS"};
  TTreeReaderArray<float> ipTriKinFit = {fReader, "ipTriKinFit"};
  TTreeReaderArray<float> ipmc = {fReader, "ipmc"};
  TTreeReaderArray<float> neuVtxTriKinFit = {fReader, "neuVtxTriKinFit"};
  TTreeReaderArray<float> photonFit1 = {fReader, "photonFit1"};
  TTreeReaderArray<float> photonFit2 = {fReader, "photonFit2"};
  TTreeReaderArray<float> photonFit3 = {fReader, "photonFit3"};
  TTreeReaderArray<float> photonFit4 = {fReader, "photonFit4"};
  TTreeReaderArray<float> pullsSignalFit = {fReader, "pullsSignalFit"};
  TTreeReaderArray<float> pullsTriKinFit = {fReader, "pullsTriKinFit"};
  TTreeReaderArray<float> pxmc = {fReader, "pxmc"};
  TTreeReaderArray<float> pymc = {fReader, "pymc"};
  TTreeReaderArray<float> pzmc = {fReader, "pzmc"};
  TTreeReaderArray<float> trcfinal = {fReader, "trcfinal"};
  TTreeReaderArray<float> trk1 = {fReader, "trk1"};
  TTreeReaderArray<float> trk1Closest = {fReader, "trk1Closest"};
  TTreeReaderArray<float> trk1Fit = {fReader, "trk1Fit"};
  TTreeReaderArray<float> trk1KL = {fReader, "trk1KL"};
  TTreeReaderArray<float> trk1KLmc = {fReader, "trk1KLmc"};
  TTreeReaderArray<float> trk1KS = {fReader, "trk1KS"};
  TTreeReaderArray<float> trk1KSmc = {fReader, "trk1KSmc"};
  TTreeReaderArray<float> trk2 = {fReader, "trk2"};
  TTreeReaderArray<float> trk2Closest = {fReader, "trk2Closest"};
  TTreeReaderArray<float> trk2Fit = {fReader, "trk2Fit"};
  TTreeReaderArray<float> trk2KL = {fReader, "trk2KL"};
  TTreeReaderArray<float> trk2KLmc = {fReader, "trk2KLmc"};
  TTreeReaderArray<float> trk2KS = {fReader, "trk2KS"};
  TTreeReaderArray<float> trk2KSmc = {fReader, "trk2KSmc"};
  TTreeReaderArray<float> xv = {fReader, "xv"};
  TTreeReaderArray<float> xvmc = {fReader, "xvmc"};
  TTreeReaderArray<float> yv = {fReader, "yv"};
  TTreeReaderArray<float> yvmc = {fReader, "yvmc"};
  TTreeReaderArray<float> zv = {fReader, "zv"};
  TTreeReaderArray<float> zvmc = {fReader, "zvmc"};
  TTreeReaderArray<float> ParamSignal = {fReader, "ParamSignal"};
  TTreeReaderArray<float> ErrorsSignal = {fReader, "ErrorsSignal"};
  TTreeReaderArray<float> ParamSignalFit = {fReader, "ParamSignalFit"};
  TTreeReaderArray<float> ErrorsSignalFit = {fReader, "ErrorsSignalFit"};
  TTreeReaderArray<float> pi01 = {fReader, "pi01"};
  TTreeReaderArray<float> pi02 = {fReader, "pi02"};
  TTreeReaderArray<float> pi01Fit = {fReader, "pi01Fit"};
  TTreeReaderArray<float> pi02Fit = {fReader, "pi02Fit"};

  MC_fit_comparison(TTree * /*tree*/ = 0) {}
  virtual ~MC_fit_comparison() {}
  virtual Int_t Version() const { return 2; }
  virtual void Begin(TTree *tree);
  virtual void SlaveBegin(TTree *tree);
  virtual void Init(TTree *tree);
  virtual Bool_t Notify();
  virtual Bool_t Process(Long64_t entry);
  virtual Int_t GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
  virtual void SetOption(const char *option) { fOption = option; }
  virtual void SetObject(TObject *obj) { fObject = obj; }
  virtual void SetInputList(TList *input) { fInput = input; }
  virtual TList *GetOutputList() const { return fOutput; }
  virtual void SlaveTerminate();
  virtual void Terminate();

  ClassDef(MC_fit_comparison, 0);
};

#endif

#ifdef MC_fit_comparison_cxx
void MC_fit_comparison::Init(TTree *tree)
{
  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the reader is initialized.
  // It is normally not necessary to make changes to the generated
  // code, but the routine can be extended by the user if needed.
  // Init() will be called many times when running on PROOF
  // (once per file to be processed).

  fReader.SetTree(tree);
}

Bool_t MC_fit_comparison::Notify()
{
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. It is normally not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed. The return value is currently not used.

  return kTRUE;
}

#endif // #ifdef MC_fit_comparison_cxx
