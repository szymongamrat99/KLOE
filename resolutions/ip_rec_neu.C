#define ip_rec_neu_cxx
// The class definition in ip_rec_neu.h has been generated automatically
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
// root> T->Process("ip_rec_neu.C")
// root> T->Process("ip_rec_neu.C","some options")
// root> T->Process("ip_rec_neu.C+")
//

#include "ip_rec_neu.h"
#include <TH2.h>
#include <TH1.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TMath.h>

#include "../../Include/const.h"
#include "../../Include/Codes/uncertainties.h"
#include "../../Include/Codes/charged_mom.h"
#include "../../Include/Codes/neutral_mom.h"
#include "../../Include/Codes/plane_intersection.h"
#include "../../Include/Codes/closest_approach.h"
#include "../../Include/Codes/lorentz_transf.h"
#include "chain_init.C"

const UInt_t N = 3;
const Float_t bin_width = 0.01;
const Float_t x_min[N] = {-1.0, -0.1, -5.0}, x_max[N] = {1.0, 0.1, 5.0};

UInt_t binn[N];
TString name;

TH1 *ip_bhabha_res[N], *ip_plane_res[N], *ip_closest_res[N], *mom_res_bhabha[N];
TCanvas *canva[N];
void ip_rec_neu::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   for (UInt_t i = 0; i < N; i++)
   {
      binn[i] = Int_t(abs(x_max[i] - x_min[i]) / bin_width);

      name = "Canva" + std::to_string(i);
      canva[i] = new TCanvas(name, "", 750, 750);

      name = "IP Bhabha res" + std::to_string(i);
      ip_bhabha_res[i] = new TH1F(name, "", binn[i], x_min[i], x_max[i]);

      name = "IP Plane res" + std::to_string(i);
      ip_plane_res[i] = new TH1F(name, "", binn[i], x_min[i], x_max[i]);

      name = "IP Closest res" + std::to_string(i);
      ip_closest_res[i] = new TH1F(name, "", binn[i], x_min[i], x_max[i]);

      name = "Mom res" + std::to_string(i);
      mom_res_bhabha[i] = new TH1F(name, "", 100, -10.0, 10.0);
   }

   TString option = GetOption();
}

void ip_rec_neu::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();
}

Bool_t ip_rec_neu::Process(Long64_t entry)
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

   Float_t ip_closest[3], ip_plane[3];
   Float_t z_axis[3] = {0., 0., 1.}, vec_kaon[3], vtx_kaon[3], vtx_bhabha[3], vec_bhabha[3], beta_bhabha[3], phi_mom[4], phi_mom_lor[4], perp_vec[3];

   fReader.SetLocalEntry(entry);

   if (*mctruth != 0 && *mctruth != 3)
   {
      for (UInt_t i = 0; i < 3; i++)
      {
         vec_kaon[i] = Knerec[i];
         vtx_kaon[i] = Knerec[i + 6];
      }

      vtx_bhabha[0] = *Bx;
      vtx_bhabha[1] = *By;
      vtx_bhabha[2] = *Bz;

      vec_bhabha[0] = *Bpx;
      vec_bhabha[1] = *Bpy;
      vec_bhabha[2] = *Bpz;

      beta_bhabha[0] = *Bpx / *Broots;
      beta_bhabha[1] = *Bpy / *Broots;
      beta_bhabha[2] = *Bpz / *Broots;

      phi_mom[0] = *Bpx;
      phi_mom[1] = *Bpy;
      phi_mom[2] = *Bpz;
      phi_mom[3] = *Broots;

      perp_vec[0] = 0.;
      perp_vec[1] = 1.;
      perp_vec[2] = 0.;

      plane_intersection(vtx_bhabha, perp_vec, vtx_kaon, vec_kaon, ip_plane);
      closest_approach(vtx_bhabha, z_axis, vtx_kaon, vec_kaon, ip_closest);

      for (UInt_t i = 0; i < N; i++)
      {
         ip_bhabha_res[i]->Fill(vtx_bhabha[i] - ipmc[i]);
         ip_plane_res[i]->Fill(ip_plane[i] - ipmc[i]);
         ip_closest_res[i]->Fill(ip_closest[i] - ipmc[i]);
      }

      
   }
   return kTRUE;
}

void ip_rec_neu::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.
}

void ip_rec_neu::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.

   TLegend *legend[N];
   Color_t color[N] = {kBlack, kRed, kBlue};

   for (UInt_t i = 0; i < N; i++)
   {
      legend[i] = new TLegend(0.2, 0.6,0.4,0.9);

      ip_bhabha_res[i]->SetLineColor(color[0]);
      ip_plane_res[i]->SetLineColor(color[1]);
      ip_closest_res[i]->SetLineColor(color[2]);

      legend[i]->AddEntry(ip_bhabha_res[i],"Bhabha res","l");
      legend[i]->AddEntry(ip_plane_res[i],"Plane res","l");
      legend[i]->AddEntry(ip_closest_res[i],"Closest res","l");

      canva[i]->cd();

      ip_bhabha_res[i]->GetXaxis()->SetTitle("X^{method}_{IP,i} - X^{gen}_{IP,i} [cm]");
      ip_bhabha_res[i]->GetYaxis()->SetTitle("Counts [cm]");
      ip_bhabha_res[i]->GetYaxis()->SetRangeUser(0.0,1.2*ip_bhabha_res[i]->GetMaximum());
      //ip_bhabha_res[i]->Fit("gaus");
      ip_bhabha_res[i]->Draw();
      ip_plane_res[i]->Draw("SAMES");
      ip_closest_res[i]->Draw("SAMES");
      legend[i]->Draw();

      name = "Coor" + std::to_string(i) + "res.png";
      canva[i]->Print(name);
   }
}