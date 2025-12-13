//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Dec  1 21:27:39 2025 by ROOT version 6.24/07
// from TChain h1/
//////////////////////////////////////////////////////////

#ifndef SignalCutOptimizer_h
#define SignalCutOptimizer_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

// Headers needed by this particular selector
#include <vector>
#include <TTree.h>
#include <TCutG.h>
#include <TMath.h>
#include <TString.h>

class SignalCutOptimizer : public TSelector
{
public:
  // Enumeracja do określenia trybu cięcia (niezależne cięcie)
  enum ECutType
  {
    kSquareOneSided = 0, // Cięcie jednostronne (np. VarX < P)
    kSquareSymmetric,    // Cięcie symetryczne (np. |VarX - Nominal| < P)
    kGeometric           // Cięcie TCutG na płaszczyźnie (VarX, VarY)
  };

  // Parametry wejściowe przekazane z makra sterującego
  ECutType fCutType;      // Wybrany typ cięcia
  TCutG *fCurrentCutG;    // Wskaźnik do obiektu TCutG (używany tylko dla kGeometric)
  TString fVarXName;      // Nazwa zmiennej X (np. "Minv")
  TString fVarYName;      // Nazwa zmiennej Y (np. "DeltaChi2")
  Double_t fParamP;       // Główny parametr cięcia (np. próg, sigma, itp.)
  Double_t fNominalValue; // Wartość nominalna dla cięcia symetrycznego (np. masa wzorcowa)

  std::map<TString, Int_t> categoryMap = {
      {"Signal", 1},
      {"Regeneration", 2},
      {"Omega", 3},
      {"3pi0", 4},
      {"Semileptonic", 5},
      {"Other", 6}};

  Int_t fmctruth = 0;

  // Wyniki
  Long64_t fNCut; // Licznik zdarzeń, które przeszły cięcie

  TTreeReader fReader; //! the tree reader
  TTree *fChain = 0;   //! pointer to the analyzed TTree or TChain

  // Readers to access the data (delete the ones you do not need).
  TTreeReaderValue<Int_t> Eclfilfo = {fReader, "Eclfilfo"};
  TTreeReaderValue<Int_t> Eclfilfoword = {fReader, "Eclfilfoword"};
  TTreeReaderValue<Int_t> bunchnum = {fReader, "bunchnum"};
  TTreeReaderValue<Int_t> cutApplied = {fReader, "cutApplied"};
  TTreeReaderValue<Int_t> errorcode = {fReader, "errorcode"};
  TTreeReaderValue<Int_t> goodClustersTriKinFitSize = {fReader, "goodClustersTriKinFitSize"};
  TTreeReaderValue<Int_t> mcflag = {fReader, "mcflag"};
  TTreeReaderValue<Int_t> mctruth = {fReader, "mctruth"};
  TTreeReaderValue<Int_t> muonAlertMinus = {fReader, "muonAlertMinus"};
  TTreeReaderValue<Int_t> muonAlertPlus = {fReader, "muonAlertPlus"};
  TTreeReaderValue<Int_t> nclu = {fReader, "nclu"};
  TTreeReaderValue<Int_t> necls = {fReader, "necls"};
  TTreeReaderValue<Int_t> nev = {fReader, "nev"};
  TTreeReaderValue<Int_t> nrun = {fReader, "nrun"};
  TTreeReaderValue<Int_t> ntcl = {fReader, "ntcl"};
  TTreeReaderValue<Int_t> ntmc = {fReader, "ntmc"};
  TTreeReaderValue<Int_t> ntv = {fReader, "ntv"};
  TTreeReaderValue<Int_t> nv = {fReader, "nv"};
  TTreeReaderValue<Int_t> nvtxmc = {fReader, "nvtxmc"};
  TTreeReaderValue<Float_t> Bpx = {fReader, "Bpx"};
  TTreeReaderValue<Float_t> Bpy = {fReader, "Bpy"};
  TTreeReaderValue<Float_t> Bpz = {fReader, "Bpz"};
  TTreeReaderValue<Float_t> Broots = {fReader, "Broots"};
  TTreeReaderValue<Float_t> Bx = {fReader, "Bx"};
  TTreeReaderValue<Float_t> By = {fReader, "By"};
  TTreeReaderValue<Float_t> Bz = {fReader, "Bz"};
  TTreeReaderValue<Float_t> Chi2OmegaKinFit = {fReader, "Chi2OmegaKinFit"};
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
  TTreeReaderValue<Float_t> KaonChTimeCMBoostRec = {fReader, "KaonChTimeCMBoostRec"};
  TTreeReaderValue<Float_t> KaonChTimeCMBoostTriFit = {fReader, "KaonChTimeCMBoostTriFit"};
  TTreeReaderValue<Float_t> KaonChTimeCMMC = {fReader, "KaonChTimeCMMC"};
  TTreeReaderValue<Float_t> KaonChTimeCMRecLor = {fReader, "KaonChTimeCMRecLor"};
  TTreeReaderValue<Float_t> KaonChTimeCMRecRec = {fReader, "KaonChTimeCMRecRec"};
  TTreeReaderValue<Float_t> KaonChTimeCMRecTriFit = {fReader, "KaonChTimeCMRecTriFit"};
  TTreeReaderValue<Float_t> KaonChTimeCMSignalFit = {fReader, "KaonChTimeCMSignalFit"};
  TTreeReaderValue<Float_t> KaonChTimeLABBoostLor = {fReader, "KaonChTimeLABBoostLor"};
  TTreeReaderValue<Float_t> KaonChTimeLABBoostRec = {fReader, "KaonChTimeLABBoostRec"};
  TTreeReaderValue<Float_t> KaonChTimeLABBoostTriFit = {fReader, "KaonChTimeLABBoostTriFit"};
  TTreeReaderValue<Float_t> KaonChTimeLABMC = {fReader, "KaonChTimeLABMC"};
  TTreeReaderValue<Float_t> KaonChTimeLABRecLor = {fReader, "KaonChTimeLABRecLor"};
  TTreeReaderValue<Float_t> KaonChTimeLABRecRec = {fReader, "KaonChTimeLABRecRec"};
  TTreeReaderValue<Float_t> KaonChTimeLABRecTriFit = {fReader, "KaonChTimeLABRecTriFit"};
  TTreeReaderValue<Float_t> KaonChTimeLABSignalFit = {fReader, "KaonChTimeLABSignalFit"};
  TTreeReaderValue<Float_t> KaonNeTimeCMBoostLor = {fReader, "KaonNeTimeCMBoostLor"};
  TTreeReaderValue<Float_t> KaonNeTimeCMBoostRec = {fReader, "KaonNeTimeCMBoostRec"};
  TTreeReaderValue<Float_t> KaonNeTimeCMBoostTriFit = {fReader, "KaonNeTimeCMBoostTriFit"};
  TTreeReaderValue<Float_t> KaonNeTimeCMMC = {fReader, "KaonNeTimeCMMC"};
  TTreeReaderValue<Float_t> KaonNeTimeCMRecLor = {fReader, "KaonNeTimeCMRecLor"};
  TTreeReaderValue<Float_t> KaonNeTimeCMRecRec = {fReader, "KaonNeTimeCMRecRec"};
  TTreeReaderValue<Float_t> KaonNeTimeCMRecTriFit = {fReader, "KaonNeTimeCMRecTriFit"};
  TTreeReaderValue<Float_t> KaonNeTimeCMSignalFit = {fReader, "KaonNeTimeCMSignalFit"};
  TTreeReaderValue<Float_t> KaonNeTimeLABBoostLor = {fReader, "KaonNeTimeLABBoostLor"};
  TTreeReaderValue<Float_t> KaonNeTimeLABBoostRec = {fReader, "KaonNeTimeLABBoostRec"};
  TTreeReaderValue<Float_t> KaonNeTimeLABBoostTriFit = {fReader, "KaonNeTimeLABBoostTriFit"};
  TTreeReaderValue<Float_t> KaonNeTimeLABMC = {fReader, "KaonNeTimeLABMC"};
  TTreeReaderValue<Float_t> KaonNeTimeLABRecLor = {fReader, "KaonNeTimeLABRecLor"};
  TTreeReaderValue<Float_t> KaonNeTimeLABRecRec = {fReader, "KaonNeTimeLABRecRec"};
  TTreeReaderValue<Float_t> KaonNeTimeLABRecTriFit = {fReader, "KaonNeTimeLABRecTriFit"};
  TTreeReaderValue<Float_t> KaonNeTimeLABSignalFit = {fReader, "KaonNeTimeLABSignalFit"};
  TTreeReaderValue<Float_t> Phiv1 = {fReader, "Phiv1"};
  TTreeReaderValue<Float_t> Phiv2 = {fReader, "Phiv2"};
  TTreeReaderValue<Float_t> PhivSmeared1 = {fReader, "PhivSmeared1"};
  TTreeReaderValue<Float_t> PhivSmeared2 = {fReader, "PhivSmeared2"};
  TTreeReaderValue<Float_t> Qmiss = {fReader, "Qmiss"};
  TTreeReaderValue<Float_t> T0step1 = {fReader, "T0step1"};
  TTreeReaderValue<Float_t> TrcSum = {fReader, "TrcSum"};
  TTreeReaderValue<Float_t> bestErrorSixGamma = {fReader, "bestErrorSixGamma"};
  TTreeReaderValue<Float_t> minv4gam = {fReader, "minv4gam"};
  TTreeReaderArray<int> Asscl = {fReader, "Asscl"};
  TTreeReaderArray<int> eclstream = {fReader, "eclstream"};
  TTreeReaderArray<int> g4takenTriKinFit = {fReader, "g4takenTriKinFit"};
  TTreeReaderArray<int> goodClustersTriKinFit = {fReader, "goodClustersTriKinFit"};
  TTreeReaderArray<int> iv = {fReader, "iv"};
  TTreeReaderArray<int> mother = {fReader, "mother"};
  TTreeReaderArray<int> pidmc = {fReader, "pidmc"};
  TTreeReaderArray<int> vtaken = {fReader, "vtaken"};
  TTreeReaderArray<int> vtakenClosest = {fReader, "vtakenClosest"};
  TTreeReaderArray<int> vtxmc = {fReader, "vtxmc"};
  TTreeReaderArray<float> Cotv = {fReader, "Cotv"};
  TTreeReaderArray<float> CotvMC = {fReader, "CotvMC"};
  TTreeReaderArray<float> Curv = {fReader, "Curv"};
  TTreeReaderArray<float> CurvMC = {fReader, "CurvMC"};
  TTreeReaderArray<float> Enecl = {fReader, "Enecl"};
  TTreeReaderArray<float> ErrorsOmega = {fReader, "ErrorsOmega"};
  TTreeReaderArray<float> ErrorsOmegaFit = {fReader, "ErrorsOmegaFit"};
  TTreeReaderArray<float> ErrorsSignal = {fReader, "ErrorsSignal"};
  TTreeReaderArray<float> ErrorsSignalFit = {fReader, "ErrorsSignalFit"};
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
  TTreeReaderArray<float> Knemc = {fReader, "Knemc"};
  TTreeReaderArray<float> Knerec = {fReader, "Knerec"};
  TTreeReaderArray<float> KnerecFit = {fReader, "KnerecFit"};
  TTreeReaderArray<float> KnerecSix = {fReader, "KnerecSix"};
  TTreeReaderArray<float> Knereclor = {fReader, "Knereclor"};
  TTreeReaderArray<float> KnereclorFit = {fReader, "KnereclorFit"};
  TTreeReaderArray<float> KnetriKinFit = {fReader, "KnetriKinFit"};
  TTreeReaderArray<float> ParamOmega = {fReader, "ParamOmega"};
  TTreeReaderArray<float> ParamOmegaFit = {fReader, "ParamOmegaFit"};
  TTreeReaderArray<float> ParamSignal = {fReader, "ParamSignal"};
  TTreeReaderArray<float> ParamSignalFit = {fReader, "ParamSignalFit"};
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
  TTreeReaderArray<float> ipOmegaFit = {fReader, "ipOmegaFit"};
  TTreeReaderArray<float> ipTriKinFit = {fReader, "ipTriKinFit"};
  TTreeReaderArray<float> ipmc = {fReader, "ipmc"};
  TTreeReaderArray<float> neuVtxTriKinFit = {fReader, "neuVtxTriKinFit"};
  TTreeReaderArray<float> omega = {fReader, "omega"};
  TTreeReaderArray<float> omegaFit = {fReader, "omegaFit"};
  TTreeReaderArray<float> phiOmegaFit = {fReader, "phiOmegaFit"};
  TTreeReaderArray<float> photonFit1 = {fReader, "photonFit1"};
  TTreeReaderArray<float> photonFit2 = {fReader, "photonFit2"};
  TTreeReaderArray<float> photonFit3 = {fReader, "photonFit3"};
  TTreeReaderArray<float> photonFit4 = {fReader, "photonFit4"};
  TTreeReaderArray<float> photonOmegaFit1 = {fReader, "photonOmegaFit1"};
  TTreeReaderArray<float> photonOmegaFit2 = {fReader, "photonOmegaFit2"};
  TTreeReaderArray<float> photonOmegaFit3 = {fReader, "photonOmegaFit3"};
  TTreeReaderArray<float> photonOmegaFit4 = {fReader, "photonOmegaFit4"};
  TTreeReaderArray<float> photonSix1 = {fReader, "photonSix1"};
  TTreeReaderArray<float> photonSix2 = {fReader, "photonSix2"};
  TTreeReaderArray<float> photonSix3 = {fReader, "photonSix3"};
  TTreeReaderArray<float> photonSix4 = {fReader, "photonSix4"};
  TTreeReaderArray<float> photonSix5 = {fReader, "photonSix5"};
  TTreeReaderArray<float> photonSix6 = {fReader, "photonSix6"};
  TTreeReaderArray<float> pi01 = {fReader, "pi01"};
  TTreeReaderArray<float> pi01Fit = {fReader, "pi01Fit"};
  TTreeReaderArray<float> pi02 = {fReader, "pi02"};
  TTreeReaderArray<float> pi02Fit = {fReader, "pi02Fit"};
  TTreeReaderArray<float> pi0Omega1 = {fReader, "pi0Omega1"};
  TTreeReaderArray<float> pi0Omega2 = {fReader, "pi0Omega2"};
  TTreeReaderArray<float> pi0OmegaFit1 = {fReader, "pi0OmegaFit1"};
  TTreeReaderArray<float> pi0OmegaFit2 = {fReader, "pi0OmegaFit2"};
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
  TTreeReaderArray<float> trk1MC = {fReader, "trk1MC"};
  TTreeReaderArray<float> trk2 = {fReader, "trk2"};
  TTreeReaderArray<float> trk2Closest = {fReader, "trk2Closest"};
  TTreeReaderArray<float> trk2Fit = {fReader, "trk2Fit"};
  TTreeReaderArray<float> trk2KL = {fReader, "trk2KL"};
  TTreeReaderArray<float> trk2KLmc = {fReader, "trk2KLmc"};
  TTreeReaderArray<float> trk2KS = {fReader, "trk2KS"};
  TTreeReaderArray<float> trk2KSmc = {fReader, "trk2KSmc"};
  TTreeReaderArray<float> trk2MC = {fReader, "trk2MC"};
  TTreeReaderArray<float> trkOmegaFit1 = {fReader, "trkOmegaFit1"};
  TTreeReaderArray<float> trkOmegaFit2 = {fReader, "trkOmegaFit2"};
  TTreeReaderArray<float> xv = {fReader, "xv"};
  TTreeReaderArray<float> xvmc = {fReader, "xvmc"};
  TTreeReaderArray<float> yv = {fReader, "yv"};
  TTreeReaderArray<float> yvmc = {fReader, "yvmc"};
  TTreeReaderArray<float> zv = {fReader, "zv"};
  TTreeReaderArray<float> zvmc = {fReader, "zvmc"};

  SignalCutOptimizer(TChain *chain = 0, Double_t limitCut = 0, Double_t centralValue = 0) : fParamP(limitCut), fNominalValue(centralValue), fNCut(0)
  {
  }

  virtual ~SignalCutOptimizer() {}
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

  // Metody do odczytu wyników z makra sterującego
  Long64_t GetNCut() const { return fNCut; }

  std::map<TString, Float_t> channLumi = {
      {"Data", 8876},
      {"Signal", 1549784},
      {"Regeneration", 1180120},
      {"Omega", 735724},
      {"3pi0", 28701},
      {"Semileptonic", 104630},
      {"Other", 42747}};

  std::map<TString, Float_t> channFactor;
  std::map<TString, Int_t> channEventsTotal, channEventsCut;

  ClassDef(SignalCutOptimizer, 0);

  private:
    std::function<bool(Double_t cutValue)> defaultFunc = [](Double_t) { return true; };

    std::map<TString, std::function<bool(Double_t cutValue)>> 
      fCutsOneSided = {
        {"Chi2SignalKinFit", defaultFunc},
        {"CombinedPi0MassErr", defaultFunc}
      }, 
      fCutsSymmetric;
};
#endif
