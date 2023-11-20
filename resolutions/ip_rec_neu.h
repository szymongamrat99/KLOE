//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sun Oct 29 01:24:13 2023 by ROOT version 6.24/07
// from TChain INTERF/h1/
//////////////////////////////////////////////////////////

#ifndef ip_rec_neu_h
#define ip_rec_neu_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

// Headers needed by this particular selector


class ip_rec_neu : public TSelector {
public :
   TTreeReader     fReader;  //!the tree reader
   TTree          *fChain = 0;   //!pointer to the analyzed TTree or TChain

   // Readers to access the data (delete the ones you do not need).
   TTreeReaderValue<Int_t> nev = {fReader, "nev"};
   TTreeReaderValue<UChar_t> pileup = {fReader, "pileup"};
   TTreeReaderValue<UChar_t> gcod = {fReader, "gcod"};
   TTreeReaderValue<UChar_t> phid = {fReader, "phid"};
   TTreeReaderValue<UChar_t> a1typ = {fReader, "a1typ"};
   TTreeReaderValue<UChar_t> a2typ = {fReader, "a2typ"};
   TTreeReaderValue<UChar_t> a3typ = {fReader, "a3typ"};
   TTreeReaderValue<UChar_t> b1typ = {fReader, "b1typ"};
   TTreeReaderValue<UChar_t> b2typ = {fReader, "b2typ"};
   TTreeReaderValue<UChar_t> b3typ = {fReader, "b3typ"};
   TTreeReaderValue<Float_t> T0step1 = {fReader, "T0step1"};
   TTreeReaderValue<UChar_t> mcflag = {fReader, "mcflag"};
   TTreeReaderValue<Float_t> Bx = {fReader, "Bx"};
   TTreeReaderValue<Float_t> By = {fReader, "By"};
   TTreeReaderValue<Float_t> Bz = {fReader, "Bz"};
   TTreeReaderValue<Float_t> Bsx = {fReader, "Bsx"};
   TTreeReaderValue<Float_t> Bsy = {fReader, "Bsy"};
   TTreeReaderValue<Float_t> Bsz = {fReader, "Bsz"};
   TTreeReaderValue<Float_t> Bpx = {fReader, "Bpx"};
   TTreeReaderValue<Float_t> Bpy = {fReader, "Bpy"};
   TTreeReaderValue<Float_t> Bpz = {fReader, "Bpz"};
   TTreeReaderValue<Float_t> Bwidpx = {fReader, "Bwidpx"};
   TTreeReaderValue<Float_t> Bwidpy = {fReader, "Bwidpy"};
   TTreeReaderValue<Float_t> Bwidpz = {fReader, "Bwidpz"};
   TTreeReaderValue<Float_t> Blumx = {fReader, "Blumx"};
   TTreeReaderValue<Float_t> Blumz = {fReader, "Blumz"};
   TTreeReaderValue<Float_t> Broots = {fReader, "Broots"};
   TTreeReaderValue<Float_t> Brootserr = {fReader, "Brootserr"};
   TTreeReaderValue<Int_t> necls2 = {fReader, "necls2"};
   TTreeReaderValue<Int_t> Ecltrgw2 = {fReader, "Ecltrgw2"};
   TTreeReaderValue<Int_t> Eclfilfo2 = {fReader, "Eclfilfo2"};
   TTreeReaderArray<Int_t> Eclword2 = {fReader, "Eclword2"};
   TTreeReaderArray<Int_t> Eclstream2 = {fReader, "Eclstream2"};
   TTreeReaderArray<Int_t> Ecltagnum2 = {fReader, "Ecltagnum2"};
   TTreeReaderArray<Int_t> Eclevtype2 = {fReader, "Eclevtype2"};
   TTreeReaderValue<Int_t> nclu = {fReader, "nclu"};
   TTreeReaderArray<Float_t> Enecl = {fReader, "Enecl"};
   TTreeReaderArray<Float_t> Tcl = {fReader, "Tcl"};
   TTreeReaderArray<Float_t> Xcl = {fReader, "Xcl"};
   TTreeReaderArray<Float_t> Ycl = {fReader, "Ycl"};
   TTreeReaderArray<Float_t> Zcl = {fReader, "Zcl"};
   TTreeReaderArray<Float_t> Xacl = {fReader, "Xacl"};
   TTreeReaderArray<Float_t> Yacl = {fReader, "Yacl"};
   TTreeReaderArray<Float_t> Zacl = {fReader, "Zacl"};
   TTreeReaderArray<Float_t> Xrmcl = {fReader, "Xrmcl"};
   TTreeReaderArray<Float_t> Yrmscl = {fReader, "Yrmscl"};
   TTreeReaderArray<Float_t> Zrmscl = {fReader, "Zrmscl"};
   TTreeReaderArray<Float_t> Trmscl = {fReader, "Trmscl"};
   TTreeReaderArray<Int_t> Flagcl = {fReader, "Flagcl"};
   TTreeReaderValue<Int_t> nv = {fReader, "nv"};
   TTreeReaderArray<Float_t> chivtx = {fReader, "chivtx"};
   TTreeReaderArray<Int_t> qualv = {fReader, "qualv"};
   TTreeReaderArray<Int_t> fitidv = {fReader, "fitidv"};
   TTreeReaderValue<Int_t> nt = {fReader, "nt"};
   TTreeReaderArray<Float_t> chi2fit = {fReader, "chi2fit"};
   TTreeReaderArray<Float_t> chi2ms = {fReader, "chi2ms"};
   TTreeReaderValue<Int_t> ntmc = {fReader, "ntmc"};
   TTreeReaderArray<UChar_t> kine = {fReader, "kine"};
   TTreeReaderArray<UChar_t> pidmc = {fReader, "pidmc"};
   TTreeReaderArray<UChar_t> virmom = {fReader, "virmom"};
   TTreeReaderArray<Float_t> pxmc = {fReader, "pxmc"};
   TTreeReaderArray<Float_t> pymc = {fReader, "pymc"};
   TTreeReaderArray<Float_t> pzmc = {fReader, "pzmc"};
   TTreeReaderArray<Float_t> themc = {fReader, "themc"};
   TTreeReaderArray<Float_t> phimc = {fReader, "phimc"};
   TTreeReaderArray<UChar_t> vtxmc = {fReader, "vtxmc"};
   TTreeReaderValue<Int_t> nvtxmc = {fReader, "nvtxmc"};
   TTreeReaderArray<UChar_t> kinmom = {fReader, "kinmom"};
   TTreeReaderArray<UChar_t> mother = {fReader, "mother"};
   TTreeReaderArray<Float_t> xvmc = {fReader, "xvmc"};
   TTreeReaderArray<Float_t> yvmc = {fReader, "yvmc"};
   TTreeReaderArray<Float_t> zvmc = {fReader, "zvmc"};
   TTreeReaderArray<Float_t> ntvtx = {fReader, "ntvtx"};
   TTreeReaderValue<UChar_t> mctruth = {fReader, "mctruth"};
   TTreeReaderValue<UChar_t> Errid = {fReader, "Errid"};
   TTreeReaderValue<UChar_t> Cutid = {fReader, "Cutid"};
   TTreeReaderArray<Float_t> Kchmc = {fReader, "Kchmc"};
   TTreeReaderArray<Float_t> Knemc = {fReader, "Knemc"};
   TTreeReaderArray<Float_t> ipmc = {fReader, "ipmc"};
   TTreeReaderArray<Float_t> Kchrec = {fReader, "Kchrec"};
   TTreeReaderArray<Float_t> Kchboost = {fReader, "Kchboost"};
   TTreeReaderArray<Float_t> Knerec = {fReader, "Knerec"};
   TTreeReaderArray<Float_t> Knereclor = {fReader, "Knereclor"};
   TTreeReaderArray<Float_t> ip = {fReader, "ip"};
   TTreeReaderValue<Float_t> chdist = {fReader, "chdist"};
   TTreeReaderValue<Float_t> cldist = {fReader, "cldist"};
   TTreeReaderValue<UChar_t> mcisr = {fReader, "mcisr"};
   TTreeReaderArray<UChar_t> vtaken = {fReader, "vtaken"};
   TTreeReaderValue<Int_t> ncl = {fReader, "ncl"};
   TTreeReaderArray<UChar_t> ncll = {fReader, "ncll"};
   TTreeReaderValue<Int_t> nclwrong = {fReader, "nclwrong"};
   TTreeReaderArray<UChar_t> ncllwrong = {fReader, "ncllwrong"};
   TTreeReaderValue<Float_t> Dlboostrec = {fReader, "Dlboostrec"};
   TTreeReaderValue<Float_t> Dtboostrec = {fReader, "Dtboostrec"};
   TTreeReaderValue<Float_t> Dlboostlor = {fReader, "Dlboostlor"};
   TTreeReaderValue<Float_t> Dtboostlor = {fReader, "Dtboostlor"};
   TTreeReaderValue<Float_t> Dlrec = {fReader, "Dlrec"};
   TTreeReaderValue<Float_t> Dtrec = {fReader, "Dtrec"};
   TTreeReaderValue<Float_t> Dlmc = {fReader, "Dlmc"};
   TTreeReaderValue<Float_t> Dtmc = {fReader, "Dtmc"};
   TTreeReaderValue<Float_t> Rc = {fReader, "Rc"};
   TTreeReaderValue<Float_t> Rtc = {fReader, "Rtc"};
   TTreeReaderValue<Float_t> Rn = {fReader, "Rn"};
   TTreeReaderValue<Float_t> Rtn = {fReader, "Rtn"};
   TTreeReaderValue<Float_t> Rcmc = {fReader, "Rcmc"};
   TTreeReaderValue<Float_t> Rtcmc = {fReader, "Rtcmc"};
   TTreeReaderValue<Float_t> Rnmc = {fReader, "Rnmc"};
   TTreeReaderValue<Float_t> Rtnmc = {fReader, "Rtnmc"};
   TTreeReaderArray<Float_t> trc = {fReader, "trc"};
   TTreeReaderArray<Float_t> trcv = {fReader, "trcv"};
   TTreeReaderArray<Float_t> pi0 = {fReader, "pi0"};
   TTreeReaderValue<Float_t> minv4gam = {fReader, "minv4gam"};
   TTreeReaderArray<Float_t> Pgamrec1 = {fReader, "Pgamrec1"};
   TTreeReaderArray<Float_t> Pgamrec2 = {fReader, "Pgamrec2"};
   TTreeReaderArray<Float_t> Pgamrec3 = {fReader, "Pgamrec3"};
   TTreeReaderArray<Float_t> Pgamrec4 = {fReader, "Pgamrec4"};
   TTreeReaderArray<Int_t> gpairtaken = {fReader, "gpairtaken"};
   TTreeReaderArray<Float_t> g4vtxerr = {fReader, "g4vtxerr"};
   TTreeReaderArray<Float_t> ominv = {fReader, "ominv"};
   TTreeReaderArray<Float_t> Ppioomega = {fReader, "Ppioomega"};
   TTreeReaderArray<Float_t> P4prirest = {fReader, "P4prirest"};
   TTreeReaderArray<Float_t> trk1 = {fReader, "trk1"};
   TTreeReaderArray<Float_t> trk2 = {fReader, "trk2"};
   TTreeReaderValue<Float_t> costrk = {fReader, "costrk"};
   TTreeReaderValue<Float_t> Qmiss = {fReader, "Qmiss"};
   TTreeReaderArray<UChar_t> g4taken = {fReader, "g4taken"};
   TTreeReaderValue<Int_t> nparfit = {fReader, "nparfit"};
   TTreeReaderArray<Float_t> Kchfit = {fReader, "Kchfit"};
   TTreeReaderArray<Float_t> Knefit = {fReader, "Knefit"};
   TTreeReaderArray<Float_t> Knereclorfit = {fReader, "Knereclorfit"};
   TTreeReaderValue<Float_t> cldistfit = {fReader, "cldistfit"};
   TTreeReaderArray<Float_t> ipfit = {fReader, "ipfit"};
   TTreeReaderValue<Int_t> nclwrongfit = {fReader, "nclwrongfit"};
   TTreeReaderArray<UChar_t> ncllwrongfit = {fReader, "ncllwrongfit"};
   TTreeReaderValue<Float_t> Dlfit = {fReader, "Dlfit"};
   TTreeReaderValue<Float_t> Dtfit = {fReader, "Dtfit"};
   TTreeReaderValue<Float_t> Rcfit = {fReader, "Rcfit"};
   TTreeReaderValue<Float_t> Rtcfit = {fReader, "Rtcfit"};
   TTreeReaderValue<Float_t> Rnfit = {fReader, "Rnfit"};
   TTreeReaderValue<Float_t> Rtnfit = {fReader, "Rtnfit"};
   TTreeReaderValue<Int_t> nclfit = {fReader, "nclfit"};
   TTreeReaderArray<Float_t> trcfit = {fReader, "trcfit"};
   TTreeReaderArray<Float_t> trcvfit = {fReader, "trcvfit"};
   TTreeReaderArray<Float_t> pi0fit = {fReader, "pi0fit"};
   TTreeReaderArray<Float_t> Pgamrec1fit = {fReader, "Pgamrec1fit"};
   TTreeReaderArray<Float_t> Pgamrec2fit = {fReader, "Pgamrec2fit"};
   TTreeReaderArray<Float_t> Pgamrec3fit = {fReader, "Pgamrec3fit"};
   TTreeReaderArray<Float_t> Pgamrec4fit = {fReader, "Pgamrec4fit"};
   TTreeReaderArray<Int_t> gpairtakenfit = {fReader, "gpairtakenfit"};
   TTreeReaderArray<Float_t> g4vtxerrfit = {fReader, "g4vtxerrfit"};
   TTreeReaderValue<Float_t> minv4gamfit = {fReader, "minv4gamfit"};
   TTreeReaderArray<Float_t> ominvfit = {fReader, "ominvfit"};
   TTreeReaderArray<Float_t> Ppioomegafit = {fReader, "Ppioomegafit"};
   TTreeReaderArray<Float_t> P4prirestfit = {fReader, "P4prirestfit"};
   TTreeReaderArray<UChar_t> g4takenfit = {fReader, "g4takenfit"};
   TTreeReaderValue<Int_t> Niter = {fReader, "Niter"};
   TTreeReaderValue<Int_t> Niterw = {fReader, "Niterw"};
   TTreeReaderValue<Float_t> Chi2 = {fReader, "Chi2"};
   TTreeReaderValue<Float_t> Chi2w = {fReader, "Chi2w"};
   TTreeReaderArray<Float_t> Fitpar = {fReader, "Fitpar"};
   TTreeReaderArray<Float_t> Errpar = {fReader, "Errpar"};
   TTreeReaderArray<Float_t> Bkgfitpar = {fReader, "Bkgfitpar"};
   TTreeReaderArray<Float_t> Bkgerrpar = {fReader, "Bkgerrpar"};
   TTreeReaderArray<Float_t> Fitparstart = {fReader, "Fitparstart"};
   TTreeReaderArray<Float_t> Errparstart = {fReader, "Errparstart"};
   TTreeReaderValue<Int_t> simok = {fReader, "simok"};
   TTreeReaderValue<UChar_t> Chvtxid = {fReader, "Chvtxid"};
   TTreeReaderArray<UChar_t> Trkidx = {fReader, "Trkidx"};
   TTreeReaderArray<UChar_t> Cluidx = {fReader, "Cluidx"};
   TTreeReaderArray<Float_t> Trkk1 = {fReader, "Trkk1"};
   TTreeReaderArray<Float_t> Trkk2 = {fReader, "Trkk2"};
   TTreeReaderArray<Float_t> Chavtx = {fReader, "Chavtx"};
   TTreeReaderArray<Float_t> Neuvtx = {fReader, "Neuvtx"};
   TTreeReaderArray<Float_t> Phivtx = {fReader, "Phivtx"};
   TTreeReaderArray<Float_t> test = {fReader, "test"};


   ip_rec_neu(TTree * /*tree*/ =0) { }
   virtual ~ip_rec_neu() { }
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

   ClassDef(ip_rec_neu,0);

};

#endif

#ifdef ip_rec_neu_cxx
void ip_rec_neu::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the reader is initialized.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   fReader.SetTree(tree);
}

Bool_t ip_rec_neu::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}


#endif // #ifdef ip_rec_neu_cxx
