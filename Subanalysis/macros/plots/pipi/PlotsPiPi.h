//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Jul 17 01:12:51 2025 by ROOT version 6.24/07
// from TTree h1/Split data
// found on file: ../../InitialAnalysis/root_files/2025-07-16/mk0_all_phys3_1.root
//////////////////////////////////////////////////////////

#ifndef PlotsPiPi_h
#define PlotsPiPi_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

// Headers needed by this particular selector
#include <vector>



class PlotsPiPi : public TSelector {
public :
   TTreeReader     fReader;  //!the tree reader
   TTree          *fChain = 0;   //!pointer to the analyzed TTree or TChain

   // Readers to access the data (delete the ones you do not need).
   TTreeReaderValue<Int_t> Eclfilfo = {fReader, "Eclfilfo"};
   TTreeReaderValue<Int_t> Eclfilfoword = {fReader, "Eclfilfoword"};
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
   TTreeReaderValue<Float_t> Bpx = {fReader, "Bpx"};
   TTreeReaderValue<Float_t> Bpy = {fReader, "Bpy"};
   TTreeReaderValue<Float_t> Bpz = {fReader, "Bpz"};
   TTreeReaderValue<Float_t> Broots = {fReader, "Broots"};
   TTreeReaderValue<Float_t> Bx = {fReader, "Bx"};
   TTreeReaderValue<Float_t> By = {fReader, "By"};
   TTreeReaderValue<Float_t> Bz = {fReader, "Bz"};
   TTreeReaderValue<Float_t> T0step1 = {fReader, "T0step1"};
   TTreeReaderArray<int> Asscl = {fReader, "Asscl"};
   TTreeReaderArray<int> eclstream = {fReader, "eclstream"};
   TTreeReaderArray<int> iv = {fReader, "iv"};
   TTreeReaderArray<int> mother = {fReader, "mother"};
   TTreeReaderArray<int> pidmc = {fReader, "pidmc"};
   TTreeReaderArray<int> vtaken = {fReader, "vtaken"};
   TTreeReaderArray<int> vtakenClosest = {fReader, "vtakenClosest"};
   TTreeReaderArray<int> vtxmc = {fReader, "vtxmc"};
   TTreeReaderArray<float> Cotv = {fReader, "Cotv"};
   TTreeReaderArray<float> Curv = {fReader, "Curv"};
   TTreeReaderArray<float> Enecl = {fReader, "Enecl"};
   TTreeReaderArray<float> Kchrec = {fReader, "Kchrec"};
   TTreeReaderArray<float> KchrecClosest = {fReader, "KchrecClosest"};
   TTreeReaderArray<float> KchrecKL = {fReader, "KchrecKL"};
   TTreeReaderArray<float> KchrecKS = {fReader, "KchrecKS"};
   TTreeReaderArray<float> Phiv = {fReader, "Phiv"};
   TTreeReaderArray<float> Tcl = {fReader, "Tcl"};
   TTreeReaderArray<float> Xcl = {fReader, "Xcl"};
   TTreeReaderArray<float> Ycl = {fReader, "Ycl"};
   TTreeReaderArray<float> Zcl = {fReader, "Zcl"};
   TTreeReaderArray<float> pxmc = {fReader, "pxmc"};
   TTreeReaderArray<float> pymc = {fReader, "pymc"};
   TTreeReaderArray<float> pzmc = {fReader, "pzmc"};
   TTreeReaderArray<float> trk1 = {fReader, "trk1"};
   TTreeReaderArray<float> trk1Closest = {fReader, "trk1Closest"};
   TTreeReaderArray<float> trk1KL = {fReader, "trk1KL"};
   TTreeReaderArray<float> trk1KS = {fReader, "trk1KS"};
   TTreeReaderArray<float> trk2 = {fReader, "trk2"};
   TTreeReaderArray<float> trk2Closest = {fReader, "trk2Closest"};
   TTreeReaderArray<float> trk2KL = {fReader, "trk2KL"};
   TTreeReaderArray<float> trk2KS = {fReader, "trk2KS"};
   TTreeReaderArray<float> xv = {fReader, "xv"};
   TTreeReaderArray<float> xvmc = {fReader, "xvmc"};
   TTreeReaderArray<float> yv = {fReader, "yv"};
   TTreeReaderArray<float> yvmc = {fReader, "yvmc"};
   TTreeReaderArray<float> zv = {fReader, "zv"};
   TTreeReaderArray<float> zvmc = {fReader, "zvmc"};


   PlotsPiPi(TTree * /*tree*/ =0) { }
   virtual ~PlotsPiPi() { }
   virtual Int_t   Version() const { return 2; }
   virtual void    Begin(TTree *tree);
   virtual void    SlaveBegin(TTree *tree);
   virtual void    Init(TTree *tree);
   virtual Bool_t  Notify();
   virtual Bool_t  Process(Long64_t entry);
   virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }
   virtual TList  *GetOutputList() const { return fOutput; }
   virtual void    SlaveTerminate();
   virtual void    Terminate();

   ClassDef(PlotsPiPi,0);

};

#endif

#ifdef PlotsPiPi_cxx
void PlotsPiPi::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the reader is initialized.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   fReader.SetTree(tree);
}

Bool_t PlotsPiPi::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}


#endif // #ifdef PlotsPiPi_cxx
