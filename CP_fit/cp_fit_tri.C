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

#include "../../Include/Codes/interference.h"
#include "../../Include/Codes/kloe_class.h"
#include "../../Include/const.h"
#include "chain_init.C"

using namespace KLOE;

void cp_fit(Bool_t check_corr = false, TString mode = "")
{
    TFile file("../../../old_root_files/mctruth.root");
    TTree *tree = (TTree*)file.Get("h1");

    Int_t mctruth = 0;

    tree->SetBranchAddress("mctruth", &mctruth);

    TFile file_tri("../../../old_root_files/neuvtx_tri_rec.root");
    TTree *tree_tri = (TTree*)file_tri.Get("h_tri");

    Int_t done4 = 0;
    Float_t fourKnetri[10] = {0.};

    tree_tri->SetBranchAddress("done4", &done4);
    tree_tri->SetBranchAddress("fourKnetri", fourKnetri);

    TChain *chain = new TChain("INTERF/h1");

    chain_init(chain);

    UChar_t mcflag;
    Float_t Dtboostlor, Dtmc, Chi2, minv4gam, Kchrec[9], Qmiss, ip[3], Kchboost[9];

    chain->SetBranchAddress("mcflag", &mcflag);

    chain->SetBranchAddress("Dtmc", &Dtmc);

    chain->SetBranchAddress("Chi2", &Chi2);
    chain->SetBranchAddress("minv4gam", &minv4gam);
    chain->SetBranchAddress("Kchrec", Kchrec);
    chain->SetBranchAddress("Qmiss", &Qmiss);

    chain->SetBranchAddress("ip", ip);

    chain->AddFriend(tree);
    chain->AddFriend(tree_tri);

    UInt_t nentries = chain->GetEntries();

    Double_t x_min = -90.0, x_max = 90.0;
    UInt_t nbins = 1 + ((x_max - x_min)/2.);

    Double_t split[3] = {-30.0, 0.0, 30.0};

    interference event(mode, check_corr, nbins, x_min, x_max, split);
    
    Double_t velocity_kch, velocity_kne;
    
    for(UInt_t i = 0; i < nentries; i++)
    {
        chain->GetEntry(i);

        if(done4 == 1)//Chi2 < 40 && abs(minv4gam - m_k0) < 76 && abs(Kchrec[5] - m_k0) < 1.2 && Qmiss < 3.75)
        {
            velocity_kch = c_vel*sqrt(pow(Kchboost[0],2) + pow(Kchboost[1],2) + pow(Kchboost[2],2))/Kchboost[3];

            velocity_kne = c_vel*sqrt(pow(fourKnetri[0],2) + pow(fourKnetri[1],2) + pow(fourKnetri[2],2))/fourKnetri[3];

            Dtboostlor = ((sqrt(pow(Kchboost[6] - ip[0],2) + pow(Kchboost[7] - ip[1],2) + pow(Kchboost[8] - ip[2],2))/velocity_kch) - (sqrt(pow(fourKnetri[6] - ip[0],2) + pow(fourKnetri[7] - ip[1],2) + pow(fourKnetri[8] - ip[2],2))/velocity_kne))/tau_S_nonCPT;

            if(mcflag == 1)
            {
                if(mctruth == 1)
                {
                    event.time_diff_gen.push_back(Dtmc);
                    event.time_diff[0].push_back(Dtboostlor);
                }

                if(mctruth == 3)
                {
                    event.time_diff[1].push_back( Dtboostlor);
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

            if(mcflag == 0)
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

    minimum->SetLimitedVariable(0, "Real part",      init_vars[0], step[0], init_vars[0] - limit_pars*init_vars[0], init_vars[0] + limit_pars*init_vars[0]);
    minimum->SetLimitedVariable(1, "Imaginary part", init_vars[1], step[1], init_vars[1] - limit_pars*init_vars[1], init_vars[1] + limit_pars*init_vars[1]);
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

    minimum->FixVariable(7);
    minimum->FixVariable(8);
    minimum->FixVariable(9);
    minimum->FixVariable(10);

    minimum->Minimize();

    std::cout << "Real part" << " " << minimum->X()[0] << std::endl;
    std::cout << "Imaginary part" << " " << minimum->X()[1] << std::endl;
    std::cout << "Norm signal" << " " << minimum->X()[2] << std::endl;
    std::cout << "Norm left DC wall" << " " << minimum->X()[3] << std::endl;
    std::cout << "Norm left beam pipe" << " " << minimum->X()[4] << std::endl;
    std::cout << "Norm right beam pipe" << " " << minimum->X()[5] << std::endl;
    std::cout << "Norm right DC wall" << " " << minimum->X()[6] << std::endl;
    std::cout << "Norm omega" << " " << minimum->X()[7] << std::endl;
    std::cout << "Norm three" << " " << minimum->X()[8] << std::endl;
    std::cout << "Norm semi" << " " << minimum->X()[9] << std::endl;
    std::cout << "Norm other bcg" << " " << minimum->X()[10] << std::endl;

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

    delete rp;

}