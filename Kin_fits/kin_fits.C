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
#include <TVectorD.h>
#include <TMatrixD.h>
#include <TCanvas.h>

#include "../../Include/Codes/uncertainties.h"
#include "../../Include/Codes/kinematic_fits.h"
#include "../../Include/Codes/charged_mom.h"
#include "../../Include/Codes/neutral_mom.h"

TH1 *chi2 = new TH1D("chi2", "", 50, 0, 1);

TH1 *pull = new TH1D("pull", "", 50, -10, 10);
TH1 *pull_init = new TH1D("pull_init", "", 50, -10, 10);

const Int_t N = 36, M = 10;
Int_t loopcount = 10;

Double_t P[N], DP[N];

TMatrixD V_final(N,N);
TVectorD P_final(N);

Double_t CHISQR = 0;

KLOE::kin_fits event(N, M, loopcount);

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

   fReader.SetLocalEntry(entry);

   if(*mctruth == 1 && *mcflag == 1)
   {

      P[0] = *Bpx;
      P[1] = *Bpy;
      P[2] = *Bpz;
      P[3] = *Broots;

      DP[0] = *Bwidpx;
      DP[1] = *Bwidpy;
      DP[2] = *Bwidpz;
      DP[3] = *Brootserr;

      for(Int_t i = 0; i < 2; i++)
      {
         P[i*3 + 4] = Curv[vtaken[i + 1] - 1];
         P[i*3 + 5] = Phiv[vtaken[i + 1] - 1];
         P[i*3 + 6] = Cotv[vtaken[i + 1] - 1];

         DP[i*3 + 4] = 0.025;
         DP[i*3 + 5] = 0.013;
         DP[i*3 + 6] = 0.012;
      }

      for(Int_t i = 0; i < 4; i++)
      {
         P[10 + i*5] = Xcl[ncll[g4taken[i] - 1] - 1];
         P[10 + i*5 + 1] = Ycl[ncll[g4taken[i] - 1] - 1];
         P[10 + i*5 + 2] = Zcl[ncll[g4taken[i] - 1] - 1];
         P[10 + i*5 + 3] = Tcl[ncll[g4taken[i] - 1] - 1];
         P[10 + i*5 + 4] = Enecl[ncll[g4taken[i] - 1] - 1];

         DP[10 + i*5] = 1.2;
         DP[10 + i*5 + 1] = 1.2;
         DP[10 + i*5 + 2] = 1.2;
         DP[10 + i*5 + 3] = clu_time_error(Enecl[ncll[g4taken[i] - 1] - 1]);
         DP[10 + i*5 + 4] = clu_ene_error(Enecl[ncll[g4taken[i] - 1] - 1]);
      }

      P[30] = Knerec[6];
      P[31] = Knerec[7];
      P[32] = Knerec[8];

      DP[30] = 0.8;
      DP[31] = 0.8;
      DP[32] = 1.1;

      P[33] = ip[0];
      P[34] = ip[1];
      P[35] = ip[2];

      DP[33] = sqrt( pow(*Bwidpx,2) + pow(*Blumx,2) );
      DP[34] = *Bsy;
      DP[35] = 0.6;

      event.FillConstructors(P, DP);

      event.solution();

      CHISQR = event.Chi2Final();

      P_final = event.FinalPars();
      V_final = event.FinalCovMtx();

      if(1){
         std::cout << TMath::Prob(CHISQR, M) << std::endl;

         chi2->Fill(TMath::Prob(CHISQR, M));
         pull->Fill((P[10] - P_final(10))/sqrt(pow(DP[10],2) - V_final(10,10)));
         pull_init->Fill((P[32] - Knemc[8])/DP[32]);
      }

   }

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
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.

   TCanvas *c1 = new TCanvas("c1", "", 790, 790);
   chi2->Draw();

   c1->Print("chi2.png");

   TCanvas *c2 = new TCanvas("c2", "", 790, 790);

   pull_init->SetLineColor(kRed);
   pull->Fit("gaus");
   pull->Draw();
   //pull_init->Draw("SAME");

   c2->Print("pull.png");

}