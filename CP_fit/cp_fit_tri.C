#include <iostream>
#include <fstream>

#include <TROOT.h>
#include <TTree.h>
#include <TFile.h>
#include <TChain.h>
#include <Math/Minimizer.h>
#include <Math/Factory.h>
#include <Math/Functor.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TRatioPlot.h>
#include <TStyle.h>
#include <TGraphAsymmErrors.h>

#include "../../Include/Codes/interference.h"
#include "../../Include/Codes/kloe_class.h"
#include "../../Include/Codes/lorentz_transf.h"
#include "../../Include/const.h"
#include "chain_init.C"

using namespace KLOE;

void cp_fit(Bool_t check_corr = false, TString mode = "")
{
    TFile file_corr("../Efficiency_analysis/correction_factor.root");

    Double_t *eff_vals;

    TGraphAsymmErrors *eff_signal = (TGraphAsymmErrors *)file_corr.Get("correction_factor");
    eff_vals = eff_signal->GetY();

    TFile file("../../../old_root_files/mctruth.root");
    TTree *tree = (TTree*)file.Get("h1");

    Int_t mctruth = 0;

    tree->SetBranchAddress("mctruth", &mctruth);

    TFile file_tri("../Neutrec/neuvtx_tri_kin_fit_1_56_10_9_1_2_bunch_corr_automatic.root");
    TTree *tree_tri = (TTree*)file_tri.Get("h_tri_kin_fit");

    Int_t done4 = 0, g4takentri_kinfit[4];
    Float_t fourKnetri[10] = {0.};

    tree_tri->SetBranchAddress("done4_kinfit", &done4);
    tree_tri->SetBranchAddress("fourKnetri_kinfit", fourKnetri);
    tree_tri->SetBranchAddress("g4takentri_kinfit", g4takentri_kinfit);

    TChain *chain = new TChain("INTERF/h1");

    chain_init(chain);

    UChar_t mcflag;
    Float_t Dtboostlor, Dtmc, Chi2, minv4gam, Kchrec[9], Qmiss, ip[3], Kchboost[9], phi_mom[4], Xcl[50], Ycl[50], Zcl[50], Tcl[50];;

    chain->SetBranchAddress("mcflag", &mcflag);

    chain->SetBranchAddress("Dtmc", &Dtmc);

    chain->SetBranchAddress("Chi2", &Chi2);
    chain->SetBranchAddress("minv4gam", &minv4gam);
    chain->SetBranchAddress("Kchboost", Kchboost);
    chain->SetBranchAddress("Qmiss", &Qmiss);

    chain->SetBranchAddress("Xcl", &Xcl);
    chain->SetBranchAddress("Ycl", &Ycl);
    chain->SetBranchAddress("Zcl", &Zcl);
    chain->SetBranchAddress("Tcl", &Tcl);

    chain->SetBranchAddress("ip", ip);

    chain->SetBranchAddress("Bpx", &phi_mom[0]);
    chain->SetBranchAddress("Bpy", &phi_mom[1]);
    chain->SetBranchAddress("Bpz", &phi_mom[2]);
    chain->SetBranchAddress("Broots", &phi_mom[3]);

    chain->AddFriend(tree);
    chain->AddFriend(tree_tri);

    UInt_t nentries = chain->GetEntries();

    Double_t x_min = -90.0, x_max = 90.0;
    UInt_t nbins = 1 + ((x_max - x_min)/2.);

    Double_t split[3] = {-30.0, 0.0, 30.0};

    interference event(mode, check_corr, nbins, x_min, x_max, split);
    
    Double_t velocity_kch, velocity_kne, tch_LAB, tch_CM, tne_LAB, tne_CM, trcv_sum, TRCV[4];
    Float_t Kch_LAB[4], Kne_LAB[4], Kch_CM[4], Kne_CM[4], Kch_CMCM[4], Kne_CMCM[4], Kne_boost[3], Kch_boost[3], Phi_boost[3], Kchmom_LAB[4], Knemom_LAB[4], Kchmom_CM[4], Knemom_CM[4];

    for(UInt_t i = 0; i < nentries; i++)
    {
        chain->GetEntry(i);

        if(done4 == 1)
        {
            velocity_kch = c_vel*sqrt(pow(Kchboost[0],2) + pow(Kchboost[1],2) + pow(Kchboost[2],2))/Kchboost[3];
            
            tch_LAB = sqrt(pow(Kchboost[6] - ip[0],2) + pow(Kchboost[7] - ip[1],2) + pow(Kchboost[8] - ip[2],2))/velocity_kch;
            tne_LAB = fourKnetri[9];

            Kch_LAB[0] = Kchboost[6] - ip[0];
            Kch_LAB[1] = Kchboost[7] - ip[1];
            Kch_LAB[2] = Kchboost[8] - ip[2];
            Kch_LAB[3] = tch_LAB * c_vel;

            Kchmom_LAB[0] = Kchboost[0];
            Kchmom_LAB[1] = Kchboost[1];
            Kchmom_LAB[2] = Kchboost[2];
            Kchmom_LAB[3] = Kchboost[3];

            Kne_LAB[0] = fourKnetri[6] - ip[0];
            Kne_LAB[1] = fourKnetri[7] - ip[1];
            Kne_LAB[2] = fourKnetri[8] - ip[2];
            Kne_LAB[3] = tne_LAB * c_vel;

            Knemom_LAB[0] = fourKnetri[0];
            Knemom_LAB[1] = fourKnetri[1];
            Knemom_LAB[2] = fourKnetri[2];
            Knemom_LAB[3] = fourKnetri[3];

            Phi_boost[0] = -phi_mom[0] / phi_mom[3];
            Phi_boost[1] = -phi_mom[1] / phi_mom[3];
            Phi_boost[2] = -phi_mom[2] / phi_mom[3]; 

            lorentz_transf(Phi_boost, Kch_LAB, Kch_CM);
            lorentz_transf(Phi_boost, Kne_LAB, Kne_CM);
            lorentz_transf(Phi_boost, Kchmom_LAB, Kchmom_CM);
            lorentz_transf(Phi_boost, Knemom_LAB, Knemom_CM);

            Kch_boost[0] = -Kchmom_CM[0] / Kchmom_CM[3];
            Kch_boost[1] = -Kchmom_CM[1] / Kchmom_CM[3];
            Kch_boost[2] = -Kchmom_CM[2] / Kchmom_CM[3];

            Kne_boost[0] = Kchmom_CM[0] / Kchmom_CM[3];
            Kne_boost[1] = Kchmom_CM[1] / Kchmom_CM[3];
            Kne_boost[2] = Kchmom_CM[2] / Kchmom_CM[3];

            // std::cout << Kch_boost[0] << " " << Kne_boost[0] << std::endl; 
            // std::cout << Kch_boost[1] << " " << Kne_boost[1] << std::endl; 
            // std::cout << Kch_boost[2] << " " << Kne_boost[2] << std::endl << std::endl; 

            lorentz_transf(Kch_boost, Kch_CM, Kch_CMCM);
            lorentz_transf(Kch_boost, Kne_CM, Kne_CMCM);

            Dtboostlor = (Kch_CMCM[3] - Kne_CMCM[3])/(c_vel * tau_S_nonCPT);

            for(Int_t i = 0; i < 4; i++) TRCV[i] = Tcl[g4takentri_kinfit[i]-1] - (sqrt(pow(Xcl[g4takentri_kinfit[i]-1] - fourKnetri[6],2) + pow(Ycl[g4takentri_kinfit[i]-1] - fourKnetri[7],2) + pow(Zcl[g4takentri_kinfit[i]-1] - fourKnetri[8],2))/c_vel) - fourKnetri[9];

            trcv_sum = (TRCV[0] + TRCV[1] + TRCV[2] + TRCV[3]);

            if(mcflag == 1 && trcv_sum > -5. && (minv4gam - m_k0) < 76 && abs(Kchboost[5] - m_k0) < 1.2 && Qmiss < 3.75)
            {
                if(mctruth == 1)
                {
                    event.time_diff_gen.push_back(Dtmc);
                    event.time_diff[0].push_back(Dtboostlor);
                }

                if(mctruth == 3)
                {
                    event.time_diff[1].push_back(Dtboostlor);
                }

                if(mctruth == 4)
                {
                    event.time_diff[2].push_back(Dtboostlor);
                }

                if(mctruth == 5)
                {
                    event.time_diff[3].push_back(Dtboostlor);
                }

                if(mctruth == 6)
                {
                    event.time_diff[4].push_back(Dtboostlor);
                }

                if(mctruth == 7)
                {
                    event.time_diff[5].push_back(Dtboostlor);
                }
            }

            if(mcflag == 0 && (minv4gam - m_k0) < 76 && abs(Kchboost[5] - m_k0) < 1.2 && Qmiss < 3.75)
            {
                event.time_diff_data.push_back(Dtboostlor);
            }
        }

    }

    ROOT::Math::Minimizer* minimum = 
        ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");

    // set tolerance , etc...
    minimum->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2
    minimum->SetTolerance(0.1);
    minimum->SetPrintLevel(1);
    minimum->SetStrategy(2);

    const UInt_t num_of_vars = 11;

    ROOT::Math::Functor minimized_function(&event, &interference::interf_chi2, num_of_vars);
 
    minimum->SetFunction(minimized_function);

    const Double_t init_vars[num_of_vars] = {Re, M_PI*Im_nonCPT/180., 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0}, 
                   step[num_of_vars] = {1E-5, 1E-5, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05};

    Double_t limit_span_signal = 0.3, limit_span = 0.3 , limit_pars = 1.0;

    minimum->SetVariable(0, "Real part",      init_vars[0], step[0]);//, init_vars[0] - limit_pars*init_vars[0], init_vars[0] + limit_pars*init_vars[0]);
    minimum->SetVariable(1, "Imaginary part", init_vars[1], step[1]);//, init_vars[1] - limit_pars*init_vars[1], init_vars[1] + limit_pars*init_vars[1]);
    minimum->SetLimitedVariable(2, "Norm signal",    init_vars[2], step[2], init_vars[2] - limit_span_signal*init_vars[2], init_vars[2] + limit_span_signal*init_vars[2]);

    if(x_min > -30.0 && x_max < 30.0)
    {
        minimum->SetFixedVariable(3, "Norm left DC wall",     init_vars[3]);
        minimum->SetLowerLimitedVariable(4, "Norm left beam pipe",     init_vars[4], step[4], 0.0);
        minimum->SetLowerLimitedVariable(5, "Norm right beam pipe",     init_vars[5], step[5], 0.0);
        minimum->SetFixedVariable(6, "Norm right DC wall",     init_vars[6]);
    }
    else
    {
        minimum->SetLowerLimitedVariable(3, "Norm left DC wall",     init_vars[3], step[3], 0.0);
        minimum->SetLowerLimitedVariable(4, "Norm left beam pipe",     init_vars[4], step[4], 0.0);
        minimum->SetLowerLimitedVariable(5, "Norm right beam pipe",     init_vars[5], step[5], 0.0);
        minimum->SetLowerLimitedVariable(6, "Norm right DC wall",     init_vars[6], step[6], 0.0);
    }


    minimum->SetLimitedVariable(7, "Norm omega",     init_vars[7], step[7], init_vars[7] - limit_span*init_vars[7], init_vars[7] + limit_span*init_vars[7]);
    minimum->SetLimitedVariable(8, "Norm three",     init_vars[8], step[8], init_vars[8] - limit_span*init_vars[8], init_vars[8] + limit_span*init_vars[8]);
    minimum->SetLimitedVariable(9, "Norm semi",      init_vars[9], step[9], init_vars[9] - limit_span*init_vars[9], init_vars[9] + limit_span*init_vars[9]);
    minimum->SetLimitedVariable(10, "Norm other bcg", init_vars[10], step[10], init_vars[10] - limit_span*init_vars[10], init_vars[10] + limit_span*init_vars[10]);

    minimum->Minimize();

        /*std::cout << "Real part" << " " << minimum->X()[0] << std::endl;
    std::cout << "Imaginary part" << " " << minimum->X()[1] << std::endl;
    std::cout << "Norm signal" << " " << minimum->X()[2] << std::endl;
    std::cout << "Norm left DC wall" << " " << minimum->X()[3] << std::endl;
    std::cout << "Norm left beam pipe" << " " << minimum->X()[4] << std::endl;
    std::cout << "Norm right beam pipe" << " " << minimum->X()[5] << std::endl;
    std::cout << "Norm right DC wall" << " " << minimum->X()[6] << std::endl;
    std::cout << "Norm omega" << " " << minimum->X()[7] << std::endl;
    std::cout << "Norm three" << " " << minimum->X()[8] << std::endl;
    std::cout << "Norm semi" << " " << minimum->X()[9] << std::endl;
    std::cout << "Norm other bcg" << " " << minimum->X()[10] << std::endl;*/

    Double_t sum_of_events = 0., fractions[6] = {0.};

    for(Int_t i = 0; i < chann_num; i++)
      sum_of_events += event.time_diff[i].size();

    for(Int_t i = 0; i < chann_num; i++)
     fractions[i] = 100*event.time_diff[i].size()/sum_of_events;

    std::ofstream myfile_num;
    myfile_num.open ("num_of_events.csv");
    myfile_num << "Channel,Number of events,Fraction\n";
    myfile_num << "Signal," << event.time_diff[0].size() << "," << fractions[0] << "%,\n";
    myfile_num << "Regeneration," << event.time_diff[1].size() << "," << fractions[1] << "%,\n";
    myfile_num << "Omega," << event.time_diff[2].size() << "," << fractions[2] << "%,\n";
    myfile_num << "Three," << event.time_diff[3].size() << "," << fractions[3] << "%,\n";
    myfile_num << "Semi," << event.time_diff[4].size() << "," << fractions[4] << "%,\n";
    myfile_num << "Other bcg," << event.time_diff[5].size() << "," << fractions[5] << "%,\n";
    myfile_num.close();


    Double_t par[2] = {minimum->X()[0], minimum->X()[1]};


    for(UInt_t i = 0; i < chann_num; i++)
    {
        for(UInt_t j = 0; j < event.time_diff[i].size(); j++)
        {
            if(i == 0){
                event.frac[i]->Fill(event.time_diff[i][j], event.interf_function(event.time_diff_gen[j], 0, par));
            }
            else if(i == 1)
            {  
                if(event.time_diff[i][j] < event.left_x_split)
                    event.frac[i]->Fill(event.time_diff[i][j], minimum->X()[3]);
                else if(event.time_diff[i][j] > event.left_x_split && event.time_diff[i][j] < event.center_x_split)
                    event.frac[i]->Fill(event.time_diff[i][j], minimum->X()[4]);
                else if(event.time_diff[i][j] > event.center_x_split && event.time_diff[i][j] < event.right_x_split)
                    event.frac[i]->Fill(event.time_diff[i][j], minimum->X()[5]);
                else if(event.time_diff[i][j] > event.right_x_split)
                    event.frac[i]->Fill(event.time_diff[i][j], minimum->X()[6]);
            }
            else
            {
                event.frac[i]->Fill(event.time_diff[i][j]);
            }
        }
    }


    for(UInt_t j = 0; j < event.time_diff_data.size(); j++)
    {
        event.data->Fill(event.time_diff_data[j]);
    }
    
    
    event.frac[0]->Scale(minimum->X()[2]*event.frac[0]->GetEntries() / event.frac[0]->Integral(0, nbins + 1) );

    for(Int_t i = 0; i < nbins; i++)
    {
        event.frac[0]->SetBinContent(i + 1, event.frac[0]->GetBinContent(i + 1)*event.corr_vals[i] );
    }

    event.frac[2]->Scale(minimum->X()[7]*event.frac[2]->GetEntries() / event.frac[2]->Integral(0, nbins + 1) );
    event.frac[3]->Scale(minimum->X()[8]*event.frac[3]->GetEntries() / event.frac[3]->Integral(0, nbins + 1) );
    event.frac[4]->Scale(minimum->X()[9]*event.frac[4]->GetEntries() / event.frac[4]->Integral(0, nbins + 1) );
    event.frac[5]->Scale(minimum->X()[10]*event.frac[5]->GetEntries() / event.frac[5]->Integral(0, nbins + 1) );

    for(UInt_t i = 0; i < chann_num; i++)
    {
        event.mc_sum->Add(event.frac[i]);

        event.frac[i]->SetLineWidth(3);
        event.frac[i]->SetLineColor(chann_color[i]);
    }


    std::ofstream myfile;
    myfile.open ("results_split_fit.csv");
    myfile << "Parameter,Value,Error\n";
    myfile << "Real part," << minimum->X()[0] << "," << minimum->Errors()[0] << ",\n";
    myfile << "Imaginary part," << minimum->X()[1] << "," << minimum->Errors()[1] << ",\n";
    myfile << "Norm signal," << minimum->X()[2] << "," << minimum->Errors()[2] << ",\n";
    myfile << "Norm left DC wall," << minimum->X()[3] << "," << minimum->Errors()[3] << ",\n";
    myfile << "Norm left beam pipe," << minimum->X()[4] << "," << minimum->Errors()[4] << ",\n";
    myfile << "Norm right beam pipe," << minimum->X()[5] << "," << minimum->Errors()[5] << ",\n";
    myfile << "Norm right DC wall," << minimum->X()[6] << "," << minimum->Errors()[6] << ",\n";
    myfile << "Norm omega," << minimum->X()[7] << "," << minimum->Errors()[7] << ",\n";
    myfile << "Norm three," << minimum->X()[8] << "," << minimum->Errors()[8] << ",\n";
    myfile << "Norm semi," << minimum->X()[9] << "," << minimum->Errors()[9] << ",\n";
    myfile << "Norm other bcg," << minimum->X()[10] << "," << minimum->Errors()[10] << ",\n";
    myfile << "\u03C7\u00B2," << event.data->Chi2Test(event.mc_sum,"UW CHI2") << ",-,\n";
    myfile << "\u03C7\u00B2/" << (UInt_t)nbins << "," << event.data->Chi2Test(event.mc_sum,"UW CHI2/NDF") << ",-,\n";
    myfile.close();

    event.mc_sum->SetLineWidth(3);
    event.mc_sum->SetLineColor(mcsum_color);

    event.data->SetLineWidth(3);
    event.data->SetLineColor(data_color);

    TCanvas *c1 = new TCanvas("c1", "", 790, 790);

    gStyle->SetOptStat(0);

    TRatioPlot *rp = new TRatioPlot(event.mc_sum, event.data, "diffsig");

    rp->Draw();

    rp->GetLowerRefGraph()->SetMinimum(-5);
    rp->GetLowerRefGraph()->SetMaximum(5);
    rp->GetLowerRefGraph()->SetLineWidth(3);

    rp->GetUpperRefXaxis()->SetTitle("#Deltat [#tau_{S}]");

    rp->SetLowBottomMargin(0.5);
    rp->SetLeftMargin(0.15);

    rp->GetLowerRefYaxis()->SetLabelSize(0.02);

    Double_t max_height = event.mc_sum->GetMaximum();

    rp->GetUpperRefYaxis()->SetRangeUser(0.0,1.2*max_height);
    rp->GetUpperRefYaxis()->SetTitle("Counts/2#tau_{S}");

    rp->GetLowerRefYaxis()->SetTitleSize(0.03);
    rp->GetLowerRefYaxis()->SetTitle("Residuals");

    rp->GetUpperPad()->cd();

    TLegend *legend_chann = new TLegend(0.6,0.5,0.9,0.9);
    legend_chann->SetFillColor(kWhite);
    for(UInt_t i = 0; i < chann_num; i++)
    {
      legend_chann->AddEntry(event.frac[i],chann_name[i],"l");
      event.frac[i]->Draw("HISTSAME");
    }

    legend_chann->AddEntry(event.mc_sum,mcsum_name,"l");
    legend_chann->AddEntry(event.data,data_name,"le");
    legend_chann->Draw();

    c1->Print("split_fit_with_corr.png");

    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Residuals graph

    TCanvas *c2 = new TCanvas("c2", "", 790, 790);
    TH1 *residuals_hist = new TH1D("Residuals hist", "", 11, -5., 5.);

    event.resi_vals = rp->GetLowerRefGraph()->GetY();

    for(Int_t i = 0; i < event.bin_number; i++)
    {
      residuals_hist->Fill( event.resi_vals[i] );
    }

    residuals_hist->GetYaxis()->SetRangeUser(0,25);
    residuals_hist->Fit("gaus");

    c2->cd();

    residuals_hist->SetXTitle("Residuals");
    residuals_hist->SetYTitle("Counts");
    residuals_hist->SetLineWidth(5);
    residuals_hist->Draw();
    residuals_hist->SetStats(1);
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(1);

    gStyle->SetFitFormat("6.2g");
    gStyle->SetStatFormat("6.2g");

    residuals_hist->Draw();

    c2->Print("residuals_hist.png");

    delete residuals_hist;
    delete c1;
    delete c2;

    delete rp;

//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////////////////////////

    interference event_final("final", true, nbins, x_min, x_max, split);

    for(Int_t i = 0; i < 8; i++)
        event_final.tmp_norm[i] = minimum->X()[i+3];

    minimum->Clear();

    for(UInt_t i = 0; i < nentries; i++)
    {
        chain->GetEntry(i);

        if(done4 == 1)
        {
            velocity_kch = c_vel*sqrt(pow(Kchboost[0],2) + pow(Kchboost[1],2) + pow(Kchboost[2],2))/Kchboost[3];
            
            tch_LAB = sqrt(pow(Kchboost[6] - ip[0],2) + pow(Kchboost[7] - ip[1],2) + pow(Kchboost[8] - ip[2],2))/velocity_kch;
            tne_LAB = fourKnetri[9];

            Kch_LAB[0] = Kchboost[6] - ip[0];
            Kch_LAB[1] = Kchboost[7] - ip[1];
            Kch_LAB[2] = Kchboost[8] - ip[2];
            Kch_LAB[3] = tch_LAB * c_vel;

            Kchmom_LAB[0] = Kchboost[0];
            Kchmom_LAB[1] = Kchboost[1];
            Kchmom_LAB[2] = Kchboost[2];
            Kchmom_LAB[3] = Kchboost[3];

            Kne_LAB[0] = fourKnetri[6] - ip[0];
            Kne_LAB[1] = fourKnetri[7] - ip[1];
            Kne_LAB[2] = fourKnetri[8] - ip[2];
            Kne_LAB[3] = tne_LAB * c_vel;

            Knemom_LAB[0] = fourKnetri[0];
            Knemom_LAB[1] = fourKnetri[1];
            Knemom_LAB[2] = fourKnetri[2];
            Knemom_LAB[3] = fourKnetri[3];

            Phi_boost[0] = -phi_mom[0] / phi_mom[3];
            Phi_boost[1] = -phi_mom[1] / phi_mom[3];
            Phi_boost[2] = -phi_mom[2] / phi_mom[3]; 

            lorentz_transf(Phi_boost, Kch_LAB, Kch_CM);
            lorentz_transf(Phi_boost, Kne_LAB, Kne_CM);
            lorentz_transf(Phi_boost, Kchmom_LAB, Kchmom_CM);
            lorentz_transf(Phi_boost, Knemom_LAB, Knemom_CM);

            Kch_boost[0] = -Kchmom_CM[0] / Kchmom_CM[3];
            Kch_boost[1] = -Kchmom_CM[1] / Kchmom_CM[3];
            Kch_boost[2] = -Kchmom_CM[2] / Kchmom_CM[3];

            Kne_boost[0] = Kchmom_CM[0] / Kchmom_CM[3];
            Kne_boost[1] = Kchmom_CM[1] / Kchmom_CM[3];
            Kne_boost[2] = Kchmom_CM[2] / Kchmom_CM[3];

            // std::cout << Kch_boost[0] << " " << Kne_boost[0] << std::endl; 
            // std::cout << Kch_boost[1] << " " << Kne_boost[1] << std::endl; 
            // std::cout << Kch_boost[2] << " " << Kne_boost[2] << std::endl << std::endl; 

            lorentz_transf(Kch_boost, Kch_CM, Kch_CMCM);
            lorentz_transf(Kch_boost, Kne_CM, Kne_CMCM);

            Dtboostlor = (Kch_CMCM[3] - Kne_CMCM[3])/(c_vel * tau_S_nonCPT);

            for(Int_t i = 0; i < 4; i++) TRCV[i] = Tcl[g4takentri_kinfit[i]-1] - (sqrt(pow(Xcl[g4takentri_kinfit[i]-1] - fourKnetri[6],2) + pow(Ycl[g4takentri_kinfit[i]-1] - fourKnetri[7],2) + pow(Zcl[g4takentri_kinfit[i]-1] - fourKnetri[8],2))/c_vel) - fourKnetri[9];

            trcv_sum = (TRCV[0] + TRCV[1] + TRCV[2] + TRCV[3]);

            if(mcflag == 1 && (minv4gam - m_k0) < 76 && abs(Kchboost[5] - m_k0) < 1.2 && Qmiss < 3.75)
            {
                if(mctruth == 1)
                {
                    event_final.time_diff_gen.push_back(Dtmc);
                    event_final.time_diff[0].push_back(Dtboostlor);
                }

                if(mctruth == 3)
                {
                    event_final.time_diff[1].push_back(Dtboostlor);
                }

                if(mctruth == 4)
                {
                    event_final.time_diff[2].push_back(Dtboostlor);
                }

                if(mctruth == 5)
                {
                    event_final.time_diff[3].push_back(Dtboostlor);
                }

                if(mctruth == 6)
                {
                    event_final.time_diff[4].push_back(Dtboostlor);
                }

                if(mctruth == 7)
                {
                    event_final.time_diff[5].push_back(Dtboostlor);
                }
            }

            if(mcflag == 0 && (minv4gam - m_k0) < 76 && abs(Kchboost[5] - m_k0) < 1.2 && Qmiss < 3.75)
            {
                event_final.time_diff_data.push_back(Dtboostlor);
            }
        }

    }

    const UInt_t num_of_vars_final = 3;

    ROOT::Math::Functor minimized_function_final(&event_final, &interference::interf_chi2, num_of_vars_final);
 
    minimum->SetFunction(minimized_function_final);

    const Double_t init_vars_final[num_of_vars_final] = {Re, M_PI*Im_nonCPT/180., 1.0}, 
                   step_final[num_of_vars_final] = {1E-5, 1E-5, 0.05};

    limit_span_signal = 0.2;

    minimum->SetVariable(0, "Real part",      init_vars_final[0], step_final[0]);
    minimum->SetVariable(1, "Imaginary part", init_vars_final[1], step_final[1]);
    minimum->SetLimitedVariable(2, "Norm signal",    init_vars_final[2], step_final[2], init_vars_final[2] - limit_span_signal*init_vars_final[2], init_vars_final[2] + limit_span_signal*init_vars_final[2]);

    minimum->Minimize();


    par[0] = minimum->X()[0];
    par[1] = minimum->X()[1];

    for(UInt_t j = 0; j < event_final.time_diff_data.size(); j++)
    {
        event_final.data->Fill(event_final.time_diff_data[j]);
    }


    for (Int_t i = 0; i < chann_num; i++)
    {
        for(UInt_t j = 0; j < event_final.time_diff[i].size(); j++)
        {   
            if(i == 0)
            { 
                event_final.frac[i]->Fill(event_final.time_diff[i][j], event_final.interf_function(event_final.time_diff_gen[j], 0, par));
            }
            else if(i == 1)
            { 
                if(event_final.time_diff[i][j] < event_final.left_x_split)
                    event_final.frac[i]->Fill(event_final.time_diff[i][j], event_final.tmp_norm[0]);
                else if(event_final.time_diff[i][j] > event_final.left_x_split && event_final.time_diff[i][j] < event_final.center_x_split)
                    event_final.frac[i]->Fill(event_final.time_diff[i][j], event_final.tmp_norm[1]);
                else if(event_final.time_diff[i][j] > event_final.center_x_split && event_final.time_diff[i][j] < event_final.right_x_split)
                    event_final.frac[i]->Fill(event_final.time_diff[i][j], event_final.tmp_norm[2]);
                else if(event_final.time_diff[i][j] > event_final.right_x_split)
                    event_final.frac[i]->Fill(event_final.time_diff[i][j], event_final.tmp_norm[3]);
            }
            else
            {
                event_final.frac[i]->Fill(event_final.time_diff[i][j], event_final.tmp_norm[i + 2]);
            }
        }

        if(i != 0)
            event_final.data->Add(event_final.frac[i],-1.);
    }


    event_final.frac[0]->Scale(minimum->X()[2] * event_final.frac[0]->GetEntries() / event_final.frac[0]->Integral());
    event_final.frac[0]->SetLineWidth(3);
    event_final.frac[0]->SetLineColor(chann_color[0]);

    for(Int_t k = 1; k <= event_final.bin_number; k++)
					event_final.frac[0]->SetBinContent(k, event_final.frac[0]->GetBinContent(k) / event_final.corr_vals[k-1]);

    for(Int_t k = 1; k <= event_final.bin_number; k++)
					event_final.data->SetBinContent(k, event_final.data->GetBinContent(k) / event_final.corr_vals[k-1]);


    std::ofstream myfile1;
    myfile.open ("results_split_fit_final.csv");
    myfile << "Parameter,Value,Error\n";
    myfile << "Real part," << minimum->X()[0] << "," << minimum->Errors()[0] << ",\n";
    myfile << "Imaginary part," << minimum->X()[1] << "," << minimum->Errors()[1] << ",\n";
    myfile << "Norm signal," << minimum->X()[2] << "," << minimum->Errors()[2] << ",\n";
    myfile << "\u03C7\u00B2," << event_final.data->Chi2Test(event_final.frac[0],"UW CHI2") << ",-,\n";
    myfile << "\u03C7\u00B2/" << (UInt_t)nbins << "," << event_final.data->Chi2Test(event_final.frac[0],"UW CHI2/NDF") << ",-,\n";
    myfile.close();

    event_final.data->SetLineWidth(3);
    event_final.data->SetLineColor(data_color);

    TCanvas *c1_final = new TCanvas("c1_final", "", 790, 790);

    gStyle->SetOptStat(0);

    TRatioPlot *rp_final = new TRatioPlot(event_final.frac[0], event_final.data, "diffsig");

    rp_final->Draw();

    rp_final->GetLowerRefGraph()->SetMinimum(-5);
    rp_final->GetLowerRefGraph()->SetMaximum(5);
    rp_final->GetLowerRefGraph()->SetLineWidth(3);

    rp_final->GetUpperRefXaxis()->SetTitle("#Deltat [#tau_{S}]");

    rp_final->SetLowBottomMargin(0.5);
    rp_final->SetLeftMargin(0.15);

    rp_final->GetLowerRefYaxis()->SetLabelSize(0.02);

    max_height = event_final.data->GetMaximum();

    rp_final->GetUpperRefYaxis()->SetRangeUser(0.0,1.5*max_height);
    rp_final->GetUpperRefYaxis()->SetTitle("Counts/2#tau_{S}");

    rp_final->GetLowerRefYaxis()->SetTitleSize(0.03);
    rp_final->GetLowerRefYaxis()->SetTitle("Residuals");

    rp_final->GetUpperPad()->cd();

    TLegend *legend_chann_final = new TLegend(0.6,0.5,0.9,0.9);
    legend_chann_final->SetFillColor(kWhite);
      
    legend_chann_final->AddEntry(event_final.frac[0],chann_name[0],"l");
    event_final.frac[0]->Draw("HISTSAME");

    legend_chann_final->AddEntry(event_final.data,data_name,"le");
    legend_chann_final->Draw();

    c1_final->Print("split_fit_with_corr_final.png");

    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Residuals graph

    TCanvas *c2_final = new TCanvas("c2_final", "", 790, 790);
    TH1 *residuals_hist_final = new TH1D("Residuals hist", "", 11, -5., 5.);

    event_final.resi_vals = rp_final->GetLowerRefGraph()->GetY();

    for(Int_t i = 0; i < event_final.bin_number; i++)
    {
      residuals_hist_final->Fill( event_final.resi_vals[i] );
    }

    residuals_hist_final->GetYaxis()->SetRangeUser(0,25);
    residuals_hist_final->Fit("gaus");

    c2_final->cd();

    residuals_hist_final->SetXTitle("Residuals");
    residuals_hist_final->SetYTitle("Counts");
    residuals_hist_final->SetLineWidth(5);
    residuals_hist_final->Draw();
    residuals_hist_final->SetStats(1);
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(1);

    gStyle->SetFitFormat("6.2g");
    gStyle->SetStatFormat("6.2g");

    residuals_hist_final->Draw();

    c2_final->Print("residuals_hist_final.png");

    delete rp_final;

}