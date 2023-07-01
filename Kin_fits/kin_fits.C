#define kin_fits_cxx
// The class definition in kin_fits.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.


// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// root> T->Process("kin_fits.C")
// root> T->Process("kin_fits.C","some options")
// root> T->Process("kin_fits.C+")
//


#include "kin_fits.h"
#include <TH2.h>
#include <TStyle.h>

#include "../../Include/Codes/uncertainties.h"

const Int_t N = 36, M = 10; 
Double_t P[N], DP[N], CHISQR, C[M], *D[N], *MTX[N+M], Z[N+M];
Int_t loopcount;

extern "C"
{
    void fit_interf_signal_(Int_t N, Int_t M, Double_t P[36], Double_t DP[36], Double_t CHISQR, Double_t C[10],
                            Double_t *D[36], Double_t *MTX[36+10], Double_t Z[36+10], Int_t loopcount);
}

void kin_fits::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();
}

void kin_fits::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

}

Bool_t kin_fits::Process(Long64_t entry)
{
   // The Process() function is called for each entry in the tree (or possibly
   // keyed object in the case of PROOF) to be processed. The entry argument
   // specifies which entry in the currently loaded tree is to be processed.
   // When processing keyed objects with PROOF, the object is already loaded
   // and is available via the fObject pointer.
   //
   // This function should contain the \"body\" of the analysis. It can contain
   // simple or elaborate selection criteria, run algorithms on the data
   // of the event and typically fill histograms.
   //
   // The processing can be stopped by calling Abort().
   //
   // Use fStatus to set the return value of TTree::Process().
   //
   // The return value is currently not used.

   CHISQR = 0.;
   loopcount = 0;

   fReader.SetLocalEntry(entry);

   for(Int_t i = 0; i < 2; i++)
   {
      P[i*3] = Curv[vtaken[i] - 1];
      P[i*3 + 1] = Phiv[vtaken[i] - 1];
      P[i*3 + 2] = Cotv[vtaken[i] - 1];

      DP[i*3] = 0.025;
      DP[i*3 + 1] = 0.013;
      DP[i*3 + 2] = 0.012;
   }

   for(Int_t i = 0; i < 4; i++)
   {
      P[6 + i*4] = Enecl[ncll[g4taken[i] - 1] - 1];
      P[6 + i*4 + 1] = Xcl[ncll[g4taken[i] - 1] - 1];
      P[6 + i*4 + 2] = Ycl[ncll[g4taken[i] - 1] - 1];
      P[6 + i*4 + 3] = Zcl[ncll[g4taken[i] - 1] - 1];
      P[6 + i*4 + 4] = Tcl[ncll[g4taken[i] - 1] - 1];

      DP[6 + i*4] = clu_ene_error(Enecl[ncll[g4taken[i] - 1] - 1]);
      DP[6 + i*4 + 1] = 1.2;
      DP[6 + i*4 + 2] = 1.2;
      DP[6 + i*4 + 3] = 1.2;
      DP[6 + i*4 + 4] = clu_time_error(Enecl[ncll[g4taken[i] - 1] - 1]);
   }

   P[26] = Knereclor[6];
   P[27] = Knereclor[7];
   P[28] = Knereclor[8];

   DP[26] = 0.8;
   DP[27] = 0.8;
   DP[28] = 1.1;

   P[29] = ip[0];
   P[30] = ip[1];
   P[31] = ip[2];

   DP[29] = sqrt( pow(*Bwidpx,2) + pow(*Blumx,2) );
   DP[30] = *Bsy;
   DP[31] = 0.6;

   P[32] = *Broots;
   P[33] = *Bpx;
   P[34] = *Bpy;
   P[35] = *Bpz;

   DP[32] = *Brootserr;
   DP[33] = *Bwidpx;
   DP[34] = *Bwidpy;
   DP[35] = *Bwidpz;

   return kTRUE;
}

void kin_fits::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.

}

void kin_fits::Terminate()
{

   fit_interf_signal_(N, M, P, DP, CHISQR, C, D, MTX, Z, loopcount);
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.

}