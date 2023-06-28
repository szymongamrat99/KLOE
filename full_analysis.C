#define full_analysis_cxx
// The class definition in full_analysis.h has been generated automatically
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
// root> T->Process("full_analysis.C")
// root> T->Process("full_analysis.C","some options")
// root> T->Process("full_analysis.C+")
//


#include "full_analysis.h"
#include <TH2.h>
#include <TStyle.h>

#include "../Include/const.h"

TH1 *hist_signal, *hist_signal_nocuts;

void full_analysis::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   hist_signal = new TH1F("hist_signal", "", 201, -100.0, 100.0);
   hist_signal_nocuts = new TH1F("hist_signal_nocuts", "", 201, -100.0, 100.0);

   TString option = GetOption();
}

void full_analysis::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

}

Bool_t full_analysis::Process(Long64_t entry)
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

   //Lorentz vectors
   TLorentzVector mom_kaon_pm, mom_kaon_00_std, mom_kaon_00_tri, mom_kaon_s, mom_kaon_l, 
                  mom_kaon_pm_alt, mom_kaon_00_std_alt, mom_kaon_00_tri_alt, mom_kaon_s_alt, mom_kaon_l_alt, mom_phi;
   TLorentzVector pos_kaon_pm, pos_kaon_00_std, pos_kaon_00_tri, pos_kaon_s, pos_kaon_l, pos_phi;

   TVector3 phi_boost, kaon_00_std_boost, kaon_00_tri_boost, kaon_s_boost, kaon_l_boost, kaon_pm_boost;

   //Kaon vars
   Float_t k_pathpm, k_betapm;
   Float_t k_path00_std, k_beta00_std;
   Float_t k_path00_tri, k_beta00_tri;
   Float_t k_pathks, k_betaks;
   Float_t k_betakl, k_pathkl;

   //trcv
   Float_t TRCV[4], trcv_sum_tri, trcv_sum_std;

   //Boolean cuts
   Bool_t cuts_signal_charged_cs[3], cuts_signal_neutral_cs[2], cuts_semi[2], cuts_three, cuts_pipi;
   Bool_t tot_cuts_semi, tot_cuts_neutral, tot_cuts_charged;

   //Auxilliary vars
   Float_t phi_vel, kaon_pm_vel, kaon_s_vel, kaon_l_vel, kaon_00_std_vel, kaon_00_tri_vel;

   Float_t tpm, t00_std, t00_tri, ts, tl;

   Float_t DeltaT_signal, DeltaT_control, DeltaT_pipi;

   Float_t trcv_sum, trcv_sum_signal;

   fReader.SetLocalEntry(entry);

   ////////////////////////////////////////////////////////////////////////////////////////////////
   //Momentum and position of Phi in LAB
   mom_phi(0) = *Bpx;
   mom_phi(1) = *Bpy;
   mom_phi(2) = *Bpz;
   mom_phi(3) = *Broots;

   phi_vel = c_vel*(sqrt(pow(mom_phi(0),2) + pow(mom_phi(1),2) + pow(mom_phi(2),2))/mom_phi(3));

   pos_phi(0) = ip[0];
   pos_phi(1) = ip[1];
   pos_phi(2) = ip[2];
   pos_phi(3) = sqrt(pow(pos_phi(0),2) + pow(pos_phi(1),2) + pow(pos_phi(2),2))/phi_vel;

   //Momentum and position of charged kaon
   mom_kaon_pm(0) = Kchboost[0];
   mom_kaon_pm(1) = Kchboost[1];
   mom_kaon_pm(2) = Kchboost[2];
   mom_kaon_pm(3) = Kchboost[3];

   kaon_pm_vel = c_vel*(sqrt(pow(mom_kaon_pm(0),2) + pow(mom_kaon_pm(1),2) + pow(mom_kaon_pm(2),2))/mom_kaon_pm(3));

   pos_kaon_pm(0) = Kchboost[6];
   pos_kaon_pm(1) = Kchboost[7];
   pos_kaon_pm(2) = Kchboost[8];
   pos_kaon_pm(3) = pos_phi(3) + sqrt(pow(pos_kaon_pm(0),2) + pow(pos_kaon_pm(1),2) + pow(pos_kaon_pm(2),2))/kaon_pm_vel;

   //Momentum and position of charged kaon
   mom_kaon_pm_alt(0) = sqrt(pow(mom_kaon_pm(3),2) - pow(m_k0,2))*(pos_kaon_pm(0) - pos_phi(0))/sqrt(pow(pos_kaon_pm(0) - pos_phi(0),2) + pow(pos_kaon_pm(1) - pos_phi(1),2) + pow(pos_kaon_pm(2) - pos_phi(2),2));
   mom_kaon_pm_alt(1) = sqrt(pow(mom_kaon_pm(3),2) - pow(m_k0,2))*(pos_kaon_pm(1) - pos_phi(1))/sqrt(pow(pos_kaon_pm(0) - pos_phi(0),2) + pow(pos_kaon_pm(1) - pos_phi(1),2) + pow(pos_kaon_pm(2) - pos_phi(2),2));;
   mom_kaon_pm_alt(2) = sqrt(pow(mom_kaon_pm(3),2) - pow(m_k0,2))*(pos_kaon_pm(2) - pos_phi(2))/sqrt(pow(pos_kaon_pm(0) - pos_phi(0),2) + pow(pos_kaon_pm(1) - pos_phi(1),2) + pow(pos_kaon_pm(2) - pos_phi(2),2));;
   mom_kaon_pm_alt(3) = mom_kaon_pm(3);

   if(sqrt(pow(Kchrec1[6] - *Bx,2) + pow(Kchrec1[7] - *By,2) + pow(Kchrec1[8] - *Bz,2)) < sqrt(pow(Kchrec2[6] - *Bx,2) + pow(Kchrec2[7] - *By,2) + pow(Kchrec2[8] - *Bz,2)))
   {
      //Momentum and position of charged kaon_s
      mom_kaon_s(0) = Kchrec1[0];
      mom_kaon_s(1) = Kchrec1[1];
      mom_kaon_s(2) = Kchrec1[2];
      mom_kaon_s(3) = Kchrec1[3];

      kaon_s_vel = c_vel*(sqrt(pow(mom_kaon_s(0),2) + pow(mom_kaon_s(1),2) + pow(mom_kaon_s(2),2))/mom_kaon_s(3));

      pos_kaon_s(0) = Kchrec1[6];
      pos_kaon_s(1) = Kchrec1[7];
      pos_kaon_s(2) = Kchrec1[8];
      pos_kaon_s(3) = pos_phi(3) + sqrt(pow(pos_kaon_s(0),2) + pow(pos_kaon_s(1),2) + pow(pos_kaon_s(2),2))/kaon_s_vel;

      //Momentum and position of charged kaon_s alternative
      mom_kaon_s_alt(0) = sqrt(pow(mom_kaon_s(3),2) - pow(m_k0,2))*(pos_kaon_s(0) - pos_phi(0))/sqrt(pow(pos_kaon_s(0) - pos_phi(0),2) + pow(pos_kaon_s(1) - pos_phi(1),2) + pow(pos_kaon_s(2) - pos_phi(2),2));
      mom_kaon_s_alt(1) = sqrt(pow(mom_kaon_s(3),2) - pow(m_k0,2))*(pos_kaon_s(1) - pos_phi(1))/sqrt(pow(pos_kaon_s(0) - pos_phi(0),2) + pow(pos_kaon_s(1) - pos_phi(1),2) + pow(pos_kaon_s(2) - pos_phi(2),2));;
      mom_kaon_s_alt(2) = sqrt(pow(mom_kaon_s(3),2) - pow(m_k0,2))*(pos_kaon_s(2) - pos_phi(2))/sqrt(pow(pos_kaon_s(0) - pos_phi(0),2) + pow(pos_kaon_s(1) - pos_phi(1),2) + pow(pos_kaon_s(2) - pos_phi(2),2));;
      mom_kaon_s_alt(3) = mom_kaon_s(3);

      //Momentum and position of charged kaon_l
      mom_kaon_l(0) = Kchrec2[0];
      mom_kaon_l(1) = Kchrec2[1];
      mom_kaon_l(2) = Kchrec2[2];
      mom_kaon_l(3) = Kchrec2[3];

      kaon_l_vel = c_vel*(sqrt(pow(mom_kaon_l(0),2) + pow(mom_kaon_l(1),2) + pow(mom_kaon_l(2),2))/mom_kaon_l(3));

      pos_kaon_l(0) = Kchrec2[6];
      pos_kaon_l(1) = Kchrec2[7];
      pos_kaon_l(2) = Kchrec2[8];
      pos_kaon_l(3) = pos_phi(3) + sqrt(pow(pos_kaon_l(0),2) + pow(pos_kaon_l(1),2) + pow(pos_kaon_l(2),2))/kaon_l_vel;

      //Momentum and position of charged kaon_l alternative
      mom_kaon_l_alt(0) = sqrt(pow(mom_kaon_l(3),2) - pow(m_k0,2))*(pos_kaon_l(0) - pos_phi(0))/sqrt(pow(pos_kaon_l(0) - pos_phi(0),2) + pow(pos_kaon_l(1) - pos_phi(1),2) + pow(pos_kaon_l(2) - pos_phi(2),2));
      mom_kaon_l_alt(1) = sqrt(pow(mom_kaon_l(3),2) - pow(m_k0,2))*(pos_kaon_l(1) - pos_phi(1))/sqrt(pow(pos_kaon_l(0) - pos_phi(0),2) + pow(pos_kaon_l(1) - pos_phi(1),2) + pow(pos_kaon_l(2) - pos_phi(2),2));;
      mom_kaon_l_alt(2) = sqrt(pow(mom_kaon_l(3),2) - pow(m_k0,2))*(pos_kaon_l(2) - pos_phi(2))/sqrt(pow(pos_kaon_l(0) - pos_phi(0),2) + pow(pos_kaon_l(1) - pos_phi(1),2) + pow(pos_kaon_l(2) - pos_phi(2),2));;
      mom_kaon_l_alt(3) = mom_kaon_l(3);
   }
   else
   {
      //Momentum and position of charged kaon_s
      mom_kaon_s(0) = Kchrec2[0];
      mom_kaon_s(1) = Kchrec2[1];
      mom_kaon_s(2) = Kchrec2[2];
      mom_kaon_s(3) = Kchrec2[3];

      kaon_s_vel = c_vel*(sqrt(pow(mom_kaon_s(0),2) + pow(mom_kaon_s(1),2) + pow(mom_kaon_s(2),2))/mom_kaon_s(3));

      pos_kaon_s(0) = Kchrec2[6];
      pos_kaon_s(1) = Kchrec2[7];
      pos_kaon_s(2) = Kchrec2[8];
      pos_kaon_s(3) = pos_phi(3) + sqrt(pow(pos_kaon_s(0),2) + pow(pos_kaon_s(1),2) + pow(pos_kaon_s(2),2))/kaon_s_vel;

      //Momentum and position of charged kaon_s alternative
      mom_kaon_s_alt(0) = sqrt(pow(mom_kaon_s(3),2) - pow(m_k0,2))*(pos_kaon_s(0) - pos_phi(0))/sqrt(pow(pos_kaon_s(0) - pos_phi(0),2) + pow(pos_kaon_s(1) - pos_phi(1),2) + pow(pos_kaon_s(2) - pos_phi(2),2));
      mom_kaon_s_alt(1) = sqrt(pow(mom_kaon_s(3),2) - pow(m_k0,2))*(pos_kaon_s(1) - pos_phi(1))/sqrt(pow(pos_kaon_s(0) - pos_phi(0),2) + pow(pos_kaon_s(1) - pos_phi(1),2) + pow(pos_kaon_s(2) - pos_phi(2),2));;
      mom_kaon_s_alt(2) = sqrt(pow(mom_kaon_s(3),2) - pow(m_k0,2))*(pos_kaon_s(2) - pos_phi(2))/sqrt(pow(pos_kaon_s(0) - pos_phi(0),2) + pow(pos_kaon_s(1) - pos_phi(1),2) + pow(pos_kaon_s(2) - pos_phi(2),2));;
      mom_kaon_s_alt(3) = mom_kaon_s(3);

      //Momentum and position of charged kaon_l
      mom_kaon_l(0) = Kchrec1[0];
      mom_kaon_l(1) = Kchrec1[1];
      mom_kaon_l(2) = Kchrec1[2];
      mom_kaon_l(3) = Kchrec1[3];

      kaon_l_vel = c_vel*(sqrt(pow(mom_kaon_l(0),2) + pow(mom_kaon_l(1),2) + pow(mom_kaon_l(2),2))/mom_kaon_l(3));

      pos_kaon_l(0) = Kchrec1[6];
      pos_kaon_l(1) = Kchrec1[7];
      pos_kaon_l(2) = Kchrec1[8];
      pos_kaon_l(3) = pos_phi(3) + sqrt(pow(pos_kaon_l(0),2) + pow(pos_kaon_l(1),2) + pow(pos_kaon_l(2),2))/kaon_l_vel;

      //Momentum and position of charged kaon_l alternative
      mom_kaon_l_alt(0) = sqrt(pow(mom_kaon_l(3),2) - pow(m_k0,2))*(pos_kaon_l(0) - pos_phi(0))/sqrt(pow(pos_kaon_l(0) - pos_phi(0),2) + pow(pos_kaon_l(1) - pos_phi(1),2) + pow(pos_kaon_l(2) - pos_phi(2),2));
      mom_kaon_l_alt(1) = sqrt(pow(mom_kaon_l(3),2) - pow(m_k0,2))*(pos_kaon_l(1) - pos_phi(1))/sqrt(pow(pos_kaon_l(0) - pos_phi(0),2) + pow(pos_kaon_l(1) - pos_phi(1),2) + pow(pos_kaon_l(2) - pos_phi(2),2));;
      mom_kaon_l_alt(2) = sqrt(pow(mom_kaon_l(3),2) - pow(m_k0,2))*(pos_kaon_l(2) - pos_phi(2))/sqrt(pow(pos_kaon_l(0) - pos_phi(0),2) + pow(pos_kaon_l(1) - pos_phi(1),2) + pow(pos_kaon_l(2) - pos_phi(2),2));;
      mom_kaon_l_alt(3) = mom_kaon_l(3);
   }

   //Momentum and position of neutral standard kaon
   mom_kaon_00_std(0) = Knereclor[0];
   mom_kaon_00_std(1) = Knereclor[1];
   mom_kaon_00_std(2) = Knereclor[2];
   mom_kaon_00_std(3) = Knereclor[3];

   kaon_00_std_vel = c_vel*(sqrt(pow(mom_kaon_00_std(0),2) + pow(mom_kaon_00_std(1),2) + pow(mom_kaon_00_std(2),2))/mom_kaon_00_std(3));

   pos_kaon_00_std(0) = Knereclor[6];
   pos_kaon_00_std(1) = Knereclor[7];
   pos_kaon_00_std(2) = Knereclor[8];
   pos_kaon_00_std(3) = pos_phi(3) + sqrt(pow(pos_kaon_00_std(0),2) + pow(pos_kaon_00_std(1),2) + pow(pos_kaon_00_std(2),2))/kaon_00_std_vel;

   //Momentum and position of neutral standard kaon
   mom_kaon_00_std_alt(0) = sqrt(pow(mom_kaon_00_std(3),2) - pow(m_k0,2))*(pos_kaon_00_std(0) - pos_phi(0))/sqrt(pow(pos_kaon_00_std(0) - pos_phi(0),2) + pow(pos_kaon_00_std(1) - pos_phi(1),2) + pow(pos_kaon_00_std(2) - pos_phi(2),2));
   mom_kaon_00_std_alt(1) = sqrt(pow(mom_kaon_00_std(3),2) - pow(m_k0,2))*(pos_kaon_00_std(1) - pos_phi(1))/sqrt(pow(pos_kaon_00_std(0) - pos_phi(0),2) + pow(pos_kaon_00_std(1) - pos_phi(1),2) + pow(pos_kaon_00_std(2) - pos_phi(2),2));;
   mom_kaon_00_std_alt(2) = sqrt(pow(mom_kaon_00_std(3),2) - pow(m_k0,2))*(pos_kaon_00_std(2) - pos_phi(2))/sqrt(pow(pos_kaon_00_std(0) - pos_phi(0),2) + pow(pos_kaon_00_std(1) - pos_phi(1),2) + pow(pos_kaon_00_std(2) - pos_phi(2),2));;
   mom_kaon_00_std_alt(3) = mom_kaon_00_std(3);

   //Momentum and position of neutral trilateration kaon
   mom_kaon_00_tri(0) = fourKnetri[0];
   mom_kaon_00_tri(1) = fourKnetri[1];
   mom_kaon_00_tri(2) = fourKnetri[2];
   mom_kaon_00_tri(3) = fourKnetri[3];

   kaon_00_tri_vel = c_vel*(sqrt(pow(mom_kaon_00_tri(0),2) + pow(mom_kaon_00_tri(1),2) + pow(mom_kaon_00_tri(2),2))/mom_kaon_00_tri(3));

   pos_kaon_00_tri(0) = fourKnetri[6];
   pos_kaon_00_tri(1) = fourKnetri[7];
   pos_kaon_00_tri(2) = fourKnetri[8];
   pos_kaon_00_tri(3) = pos_phi(3) + sqrt(pow(pos_kaon_00_tri(0),2) + pow(pos_kaon_00_tri(1),2) + pow(pos_kaon_00_tri(2),2))/kaon_00_tri_vel;

   //Momentum and position of neutral standard kaon
   mom_kaon_00_tri_alt(0) = sqrt(pow(mom_kaon_00_tri(3),2) - pow(m_k0,2))*(pos_kaon_00_tri(0) - pos_phi(0))/sqrt(pow(pos_kaon_00_tri(0) - pos_phi(0),2) + pow(pos_kaon_00_tri(1) - pos_phi(1),2) + pow(pos_kaon_00_tri(2) - pos_phi(2),2));
   mom_kaon_00_tri_alt(1) = sqrt(pow(mom_kaon_00_tri(3),2) - pow(m_k0,2))*(pos_kaon_00_tri(1) - pos_phi(1))/sqrt(pow(pos_kaon_00_tri(0) - pos_phi(0),2) + pow(pos_kaon_00_tri(1) - pos_phi(1),2) + pow(pos_kaon_00_tri(2) - pos_phi(2),2));;
   mom_kaon_00_tri_alt(2) = sqrt(pow(mom_kaon_00_tri(3),2) - pow(m_k0,2))*(pos_kaon_00_tri(2) - pos_phi(2))/sqrt(pow(pos_kaon_00_tri(0) - pos_phi(0),2) + pow(pos_kaon_00_tri(1) - pos_phi(1),2) + pow(pos_kaon_00_tri(2) - pos_phi(2),2));;
   mom_kaon_00_tri_alt(3) = mom_kaon_00_tri(3);
   ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   // Lorentz transformation to CM of Phi

   phi_boost(0) = -mom_phi(0)/mom_phi(3);
   phi_boost(1) = -mom_phi(1)/mom_phi(3);
   phi_boost(2) = -mom_phi(2)/mom_phi(3);

   pos_phi.Boost(phi_boost);
   pos_kaon_pm.Boost(phi_boost);
   pos_kaon_00_std.Boost(phi_boost);
   pos_kaon_00_tri.Boost(phi_boost);
   pos_kaon_s.Boost(phi_boost);
   pos_kaon_l.Boost(phi_boost);

   mom_phi.Boost(phi_boost);
   mom_kaon_pm.Boost(phi_boost);
   mom_kaon_00_std.Boost(phi_boost);
   mom_kaon_00_tri.Boost(phi_boost);
   mom_kaon_s.Boost(phi_boost);
   mom_kaon_l.Boost(phi_boost);

   mom_phi.Boost(phi_boost);
   mom_kaon_pm_alt.Boost(phi_boost);
   mom_kaon_00_std_alt.Boost(phi_boost);
   mom_kaon_00_tri_alt.Boost(phi_boost);
   mom_kaon_s_alt.Boost(phi_boost);
   mom_kaon_l_alt.Boost(phi_boost);

   kaon_pm_boost(0) = -mom_kaon_pm_alt(0)/mom_kaon_pm_alt(3);
   kaon_pm_boost(1) = -mom_kaon_pm_alt(1)/mom_kaon_pm_alt(3);
   kaon_pm_boost(2) = -mom_kaon_pm_alt(2)/mom_kaon_pm_alt(3);

   kaon_00_std_boost(0) = -mom_kaon_00_std_alt(0)/mom_kaon_00_std_alt(3);
   kaon_00_std_boost(1) = -mom_kaon_00_std_alt(1)/mom_kaon_00_std_alt(3);
   kaon_00_std_boost(2) = -mom_kaon_00_std_alt(2)/mom_kaon_00_std_alt(3);

   kaon_00_tri_boost(0) = -mom_kaon_00_tri_alt(0)/mom_kaon_00_tri_alt(3);
   kaon_00_tri_boost(1) = -mom_kaon_00_tri_alt(1)/mom_kaon_00_tri_alt(3);
   kaon_00_tri_boost(2) = -mom_kaon_00_tri_alt(2)/mom_kaon_00_tri_alt(3);

   kaon_s_boost(0) = -mom_kaon_s_alt(0)/mom_kaon_s_alt(3);
   kaon_s_boost(1) = -mom_kaon_s_alt(1)/mom_kaon_s_alt(3);
   kaon_s_boost(2) = -mom_kaon_s_alt(2)/mom_kaon_s_alt(3);

   kaon_l_boost(0) = -mom_kaon_l_alt(0)/mom_kaon_l_alt(3);
   kaon_l_boost(1) = -mom_kaon_l_alt(1)/mom_kaon_l_alt(3);
   kaon_l_boost(2) = -mom_kaon_l_alt(2)/mom_kaon_l_alt(3);

   pos_kaon_pm.Boost(kaon_pm_boost);
   pos_kaon_00_std.Boost(kaon_00_std_boost);
   pos_kaon_00_tri.Boost(kaon_00_tri_boost);
   pos_kaon_s.Boost(kaon_s_boost);
   pos_kaon_l.Boost(kaon_l_boost);

   mom_kaon_pm.Boost(kaon_pm_boost);
   mom_kaon_00_std.Boost(kaon_00_std_boost);
   mom_kaon_00_tri.Boost(kaon_00_tri_boost);
   mom_kaon_s.Boost(kaon_s_boost);
   mom_kaon_l.Boost(kaon_l_boost);

   // Calculation of time difference

   tpm = pos_kaon_pm(3)/tau_S_nonCPT;
   t00_std = pos_kaon_00_std(3)/tau_S_nonCPT;
   t00_tri = pos_kaon_00_tri(3)/tau_S_nonCPT;
   ts = pos_kaon_s(3)/tau_S_nonCPT;
   tl = pos_kaon_l(3)/tau_S_nonCPT; 

   DeltaT_signal = tpm - t00_std;
   DeltaT_control = tpm - t00_tri;
   DeltaT_pipi = ts - tl;

   ///////////////////////////////////////////////////////////////////////////////////////////////

   //for(Int_t i = 0; i < 4; i++) TRCV[i] = Tcl[fourg4taken[i]] - (sqrt(pow(Xcl[fourg4taken[i]] - fourKnetri[6],2) + pow(Ycl[fourg4taken[i]] - fourKnetri[7],2) + pow(Zcl[fourg4taken[i]] - fourKnetri[8],2))/c_vel) - (k_path00/(k_beta00*c_vel));

   //trcv_sum = TRCV[0] + TRCV[1] + TRCV[2] + TRCV[3];

   cuts_semi[0] = abs(*Qmiss_inv - 115.31) < 13;
   cuts_semi[1] = abs(*anglepipi_CM_kch - 141) < 30;

   cuts_signal_neutral_cs[0] = 1;//(abs(fourKnetri[5] - m_k0) < 76);
   cuts_signal_neutral_cs[1] = 1;//(trcv_sum > -5);

   cuts_signal_charged_cs[0] = 1;//(abs(Kchrec[5] - m_k0) < 1.2);
   cuts_signal_charged_cs[1] = 1;//(*Qmiss_inv < 3.75);
   cuts_signal_charged_cs[2] = 1;//= (cos(M_PI**anglepipi_CM_kch/180.) < -0.8);

   tot_cuts_semi = cuts_semi[0] && cuts_semi[1];
   tot_cuts_neutral = cuts_signal_neutral_cs[0] && cuts_signal_neutral_cs[1];
   tot_cuts_charged = cuts_signal_charged_cs[0] && cuts_signal_charged_cs[1] && cuts_signal_charged_cs[2];


   ///////////////////////////////////////////////////////////////////////////////////////////////////
   //trcv_sum_signal = trcv[ncll[g4taken[0]-1]-1] + trcv[ncll[g4taken[1]-1]-1] + trcv[ncll[g4taken[2]-1]-1] + trcv[ncll[g4taken[3]-1]-1];
   if(*mctruth == 1)
	{
		hist_signal_nocuts->Fill(DeltaT_signal);
   }

   if(/*trcv_sum_signal > -1 &&*/ abs(*minv4gam - m_k0) < 76 && abs(Kchrec[5] - m_k0) < 1.2 /*&& *Qmiss_inv < 3.75 && cos(M_PI*(*anglepipi_CM_kch/180.)) < -0.8*/)
   {
   	if(*mctruth == 1)
	   {
		   hist_signal->Fill(DeltaT_signal);
	   }
   }

   return kTRUE;
}

void full_analysis::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.

}

void full_analysis::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.

   TCanvas *canva = new TCanvas("canva","", 750, 750);

   //hist_signal->Divide(hist_signal_nocuts);

   hist_signal->Draw("HIST");

   canva->Print("efficiency.png");

}
