#include <iostream>

#include "TChain.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TTree.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TStyle.h"

#include "../../Include/const.h"
#include "../../Include/Codes/lorentz_transf.h"

#include "chain_init.C"

using namespace std;

int bunch_crossing()
{
  TChain *chain = new TChain("INTERF/h1");

  chain_init(chain,1,56);

  Float_t Knemc[9], Kchmc[9], Knereclor[9], Tcl[50], Xcl[50], Ycl[50], Zcl[50], ip[3], ipmc[3];
  UChar_t mctruth, g4taken[4], ncll[50];

  chain->SetBranchAddress("Kchmc",Kchmc);
  chain->SetBranchAddress("Knemc",Knemc);
  chain->SetBranchAddress("Xcl",Xcl);
  chain->SetBranchAddress("Ycl",Ycl);
  chain->SetBranchAddress("Zcl",Zcl);
  chain->SetBranchAddress("Tcl",Tcl);
  chain->SetBranchAddress("g4taken",g4taken);
  chain->SetBranchAddress("ncll",ncll);
  chain->SetBranchAddress("ip",ip);
  chain->SetBranchAddress("ipmc",ipmc);
  chain->SetBranchAddress("Knerec",Knereclor);
  chain->SetBranchAddress("mctruth",&mctruth);

  UInt_t nentries = (UInt_t)chain->GetEntries();

  //! Variables
  Float_t kaon_vel, kaon_mom, kaon_time, kaon_path, Knemc_init[4], Knemc_lor[4], gamma_path[4];
  Float_t phi_vel[3], t0, t_bunch = 2.715;
  //!

  //! Histograms
  TH2 *vel_histo = new TH2F("vel_histo","",100,-20.0,5.0, 100.0, 0.0, 20.0);
  //!

  for(UInt_t i = 0; i < nentries; i++)
  {
    chain->GetEntry(i);

    if(mctruth == 1 || mctruth == 2)
    {
      phi_vel[0] = (Kchmc[0] + Knemc[0])/(Kchmc[3] + Knemc[3]);
      phi_vel[1] = (Kchmc[1] + Knemc[1])/(Kchmc[3] + Knemc[3]);
      phi_vel[2] = (Kchmc[2] + Knemc[2])/(Kchmc[3] + Knemc[3]);

      Knemc_init[0] = Knemc[0];
      Knemc_init[1] = Knemc[1];
      Knemc_init[2] = Knemc[2];
      Knemc_init[3] = Knemc[3];

      lorentz_transf(phi_vel,Knemc_init,Knemc_lor);

      kaon_mom = sqrt(pow(Knemc_init[0],2) + pow(Knemc_init[1],2) + pow(Knemc_init[2],2));
      kaon_vel = c_vel*kaon_mom/Knemc_init[3];
      kaon_path = sqrt(pow(Knemc[6] - ipmc[0],2) + pow(Knemc[7] - ipmc[1],2) + pow(Knemc[8] - ipmc[2],2));

      for(Int_t j = 0; j < 4; j++)
      {
        gamma_path[j] = sqrt(pow(Xcl[ncll[g4taken[j] - 1] - 1] - Knemc[6],2) + 
                            pow(Ycl[ncll[g4taken[j] - 1] - 1] - Knemc[7],2) +
                            pow(Zcl[ncll[g4taken[j] - 1] - 1] - Knemc[8],2));

        // t0 = -1.*t_bunch;
        kaon_time = Tcl[ncll[g4taken[j] - 1] - 1] - (gamma_path[j]/c_vel) - kaon_path/kaon_vel;

        // if(kaon_time < 2.715/2. && kaon_time > -2.715/2.) vel_histo->Fill(kaon_time);
        // else if(kaon_time < 2.715 + 2.715/2. && kaon_time > 2.715/2.) vel_histo->Fill(kaon_time - 2.715);
        // else if(kaon_time < 2*2.715 + 2.715/2. && kaon_time > 2.715 + 2.715/2.) vel_histo->Fill(kaon_time - 2*2.715);
        // else if(kaon_time < 3*2.715 + 2.715/2. && kaon_time > 2*2.715 + 2.715/2.) vel_histo->Fill(kaon_time - 3*2.715);
        // else if(kaon_time < 4*2.715 + 2.715/2. && kaon_time > 3*2.715 + 2.715/2.) vel_histo->Fill(kaon_time - 4*2.715);
        // else if(kaon_time < -2.715/2. && kaon_time > -2.715 - 2.715/2.) vel_histo->Fill(kaon_time + 2.715);
        // else if(kaon_time < -2.715 - 2.715/2. && kaon_time > -2*2.715 - 2.715/2.) vel_histo->Fill(kaon_time + 2*2.715);
        // else vel_histo->Fill(kaon_time);
        vel_histo->Fill(kaon_time/t_bunch,(gamma_path[j]/c_vel));
      }
    }
  }

  gStyle->SetOptStat("MiR");

  TFitResultPtr r;
  Float_t parameter[9];

  TF1 *func = new TF1("f1", "[0]*exp(-0.5*pow((x-[1])/[2],2)) + [3]*exp(-0.5*pow((x-[4])/[5],2)) + [6]*exp(-0.5*pow((x-[7])/[8],2))", -20.0, 20.0);
  func->SetParNames("Constant_l", "Mean_l", "Sigma_l", "Constant_r", "Mean_r", "Sigma_r");

  parameter[0] = 8E4;
  parameter[1] = 0.0;
  parameter[2] = 0.1;

  parameter[3] = 8E4;
  parameter[4] = 2.715;
  parameter[5] = 0.1;

  parameter[6] = 8E4;
  parameter[7] = 2*2.715;
  parameter[8] = 0.1;

  func->SetParameters(parameter[0], parameter[1], parameter[2], parameter[3], parameter[4], parameter[5], parameter[6], parameter[7], parameter[8]);

  func->SetParLimits(0,1E2,1E6);
  func->SetParLimits(1,-0.2,0.2);
  func->SetParLimits(2,0.0001,100.);
  func->SetParLimits(3,1E2,1E6);
  func->SetParLimits(4,2.515,2.915);
  func->SetParLimits(5,0.0001,100.);
  func->SetParLimits(6,1E2,1E6);
  func->SetParLimits(7,2*2.715 - 0.1,2*2.715 + 0.1);
  func->SetParLimits(8,0.0001,100.);

  TCanvas *c1 = new TCanvas("c1","",750,750);
  c1->SetLogz(1);

  //vel_histo->GetYaxis()->SetRangeUser(1., 1.2*vel_histo->GetMaximum());
  // vel_histo->Fit(func);
  vel_histo->Draw("COLZ");

  c1->Print("kaon_vel.png");

  return 0;
}