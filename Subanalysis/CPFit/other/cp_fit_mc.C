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

void cp_fit(Bool_t check = false, TString mode = "")
{
    TFile file("../../../old_root_files/mctruth.root");
    TTree *tree = (TTree*)file.Get("h1");

    Int_t mctruth = 0;

    tree->SetBranchAddress("mctruth", &mctruth);

    TChain *chain = new TChain("INTERF/h1");

    chain_init(chain);

    UChar_t mcflag;
    Float_t Dtboostlor, Dtmc, Chi2, minv4gam, Kchrec[9], Qmiss;

    chain->SetBranchAddress("mcflag", &mcflag);

    chain->SetBranchAddress("Dtboostlor", &Dtboostlor);
    chain->SetBranchAddress("Dtmc", &Dtmc);

    chain->SetBranchAddress("Dtmc", &Dtmc);

    chain->SetBranchAddress("Chi2", &Chi2);
    chain->SetBranchAddress("minv4gam", &minv4gam);
    chain->SetBranchAddress("Kchrec", Kchrec);
    chain->SetBranchAddress("Qmiss", &Qmiss);

    chain->AddFriend("h = h1", "../../../old_root_files/mctruth.root");

    UInt_t nentries = chain->GetEntries();

    Double_t x_min = -90.0, x_max = 90.0;
    UInt_t nbins = 1 + (x_max - x_min)/2.;

    Double_t split[3] = {-30.0, 0.0, 30.0};

    interference event(mode, check, nbins, x_min, x_max, split);
    
    
    for(UInt_t i = 0; i < nentries; i++)
    {
        chain->GetEntry(i);
        tree->GetEntry(i);

        if(Chi2 < 40 && abs(minv4gam - mK0) < 76 && abs(Kchrec[5] - mK0) < 1.2 && Qmiss < 3.75)
        {
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

    event.mc_randomize();

    ROOT::Math::Minimizer* minimum = 
        ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");

    // set tolerance , etc...
    minimum->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2
    minimum->SetTolerance(0.1);
    minimum->SetPrintLevel(1);
    minimum->SetStrategy(2);

    const UInt_t num_of_vars = 8;

    ROOT::Math::Functor minimized_function(&event, &interference::interf_chi2, num_of_vars);
 
    minimum->SetFunction(minimized_function);

    const Double_t init_vars[num_of_vars] = {Re, M_PI*Im_nonCPT/180., 1.0, 1.0, 1.0, 1.0, 1.0, 1.0}, 
                   step[num_of_vars] = {1E-5, 1E-5, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05};

    Double_t limit_span_signal = 0.3, limit_span = 0.3 , limit_pars = 1.0;

    minimum->SetLimitedVariable(0, "Real part",      init_vars[0], step[0], init_vars[0] - limit_pars*init_vars[0], init_vars[0] + limit_pars*init_vars[0]);
    minimum->SetLimitedVariable(1, "Imaginary part", init_vars[1], step[1], init_vars[1] - limit_pars*init_vars[1], init_vars[1] + limit_pars*init_vars[1]);
    minimum->SetLimitedVariable(2, "Norm signal",    init_vars[2], step[2], init_vars[2] - limit_span_signal*init_vars[2], init_vars[2] + limit_span_signal*init_vars[2]);
    minimum->SetLimitedVariable(3, "Norm regen",    init_vars[3], step[3], init_vars[3] - limit_span*init_vars[3], init_vars[3] + limit_span*init_vars[3]);
    minimum->SetLimitedVariable(4, "Norm omega",     init_vars[4], step[4], init_vars[4] - limit_span*init_vars[4], init_vars[4] + limit_span*init_vars[4]);
    minimum->SetLimitedVariable(5, "Norm three",     init_vars[5], step[5], init_vars[5] - limit_span*init_vars[5], init_vars[5] + limit_span*init_vars[5]);
    minimum->SetLimitedVariable(6, "Norm semi",      init_vars[6], step[6], init_vars[6] - limit_span*init_vars[6], init_vars[6] + limit_span*init_vars[6]);
    minimum->SetLimitedVariable(7, "Norm other bcg", init_vars[7], step[7], init_vars[7] - limit_span*init_vars[7], init_vars[7] + limit_span*init_vars[7]);

    minimum->FixVariable(4);
    minimum->FixVariable(5);
    minimum->FixVariable(6);
    minimum->FixVariable(7);

    minimum->Minimize();

    std::cout << "Real part" << " " << minimum->X()[0] << std::endl;
    std::cout << "Imaginary part" << " " << minimum->X()[1] << std::endl;
    std::cout << "Norm signal" << " " << minimum->X()[2] << std::endl;
    std::cout << "Norm regen" << " " << minimum->X()[3] << std::endl;
    std::cout << "Norm omega" << " " << minimum->X()[4] << std::endl;
    std::cout << "Norm three" << " " << minimum->X()[5] << std::endl;
    std::cout << "Norm semi" << " " << minimum->X()[6] << std::endl;
    std::cout << "Norm other bcg" << " " << minimum->X()[7] << std::endl;

    Double_t par[2] = {minimum->X()[0], minimum->X()[1]};

    for(UInt_t i = 0; i < channNum; i++)
    {

        std::cout <<  event.time_diff_rand_mc[i].size() << std::endl;
        std::cout <<  event.time_diff_rand_data[i].size() << std::endl;

        for(UInt_t j = 0; j < event.time_diff_rand_mc[i].size(); j++)
        {
            if(i == 0){
                event.frac[i]->Fill(event.time_diff_rand_mc[i][j], event.interf_function(event.time_diff_gen_rand_mc[j], 0, par));
            }
            else
            {
                event.frac[i]->Fill(event.time_diff_rand_mc[i][j]);
            }
        }

        for(UInt_t j = 0; j < event.time_diff_rand_data[i].size(); j++)
        {
            if(i == 0){
                event.frac_data[i]->Fill(event.time_diff_rand_data[i][j], event.interf_function(event.time_diff_gen_rand_data[j]));
            }
            else
            {
                event.frac_data[i]->Fill(event.time_diff_rand_data[i][j]);
            }
        }

        event.frac[i]->Scale(minimum->X()[i + 2]*event.frac[i]->GetEntries() / event.frac[i]->Integral(0, nbins + 1) );
        event.frac_data[i]->Scale(event.frac_data[i]->GetEntries() / event.frac_data[i]->Integral(0, nbins + 1) );

        event.frac[i]->SetLineColor(channColor[i]);

        event.mc_sum->Add(event.frac[i]);
        event.data->Add(event.frac_data[i]);
    }

    std::ofstream myfile;
    myfile.open ("results_mc_fit.csv");
    myfile << "Parameter,Value,Error\n";
    myfile << "Real part," << minimum->X()[0] << "," << minimum->Errors()[0] << ",\n";
    myfile << "Imaginary part," << minimum->X()[1] << "," << minimum->Errors()[1] << ",\n";
    myfile << "Norm signal," << minimum->X()[2] << "," << minimum->Errors()[2] << ",\n";
    myfile << "Norm regen," << minimum->X()[3] << "," << minimum->Errors()[3] << ",\n";
    myfile << "Norm omega," << minimum->X()[4] << "," << minimum->Errors()[4] << ",\n";
    myfile << "Norm three," << minimum->X()[5] << "," << minimum->Errors()[5] << ",\n";
    myfile << "Norm semi," << minimum->X()[6] << "," << minimum->Errors()[6] << ",\n";
    myfile << "Norm other bcg," << minimum->X()[7] << "," << minimum->Errors()[7] << ",\n";
    myfile << "\u03C7\u00B2," << event.data->Chi2Test(event.mc_sum,"WW CHI2") << ",-,\n";
    myfile << "\u03C7\u00B2/" << (UInt_t)nbins << "," << event.data->Chi2Test(event.mc_sum,"WW CHI2/NDF") << ",-,\n";
    myfile.close();

    event.mc_sum->SetLineColor(mcSumColor);
    event.data->SetLineColor(dataColor);

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

    for(UInt_t i = 0; i < channNum; i++)
    {
        event.frac[i]->Draw("HISTSAME");
    }

    TLegend *legend_chann = new TLegend(0.6,0.5,0.9,0.9);
    legend_chann->SetFillColor(kWhite);
    for(UInt_t i = 0; i < channNum; i++)
    {
      legend_chann->AddEntry(event.frac[i],channName[i],"l");
      event.frac[i]->Draw("HISTSAME");
    }

    legend_chann->AddEntry(event.mc_sum,mcSumName,"l");
    legend_chann->AddEntry(event.data,dataName,"le");
    legend_chann->Draw();

    //gStyle->SetOptStat(0);

    c1->Print("mc_bcg_fit_with_corr.png");

    // Residuals graph

    TCanvas *c2 = new TCanvas("c2", "", 790, 790);
    TH1 *residuals_hist = new TH1D("Residuals hist", "", 11, -5., 5.);

    event.resi_vals = rp->GetLowerRefGraph()->GetY();

    for(Int_t i = 0; i < event.bin_number; i++)
    {
      residuals_hist->Fill( event.resi_vals[i] );
    }

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


    c2->Print("residuals_hist_mc.png");

    delete rp;

}