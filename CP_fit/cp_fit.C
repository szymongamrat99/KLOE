#include "../../Include/Codes/interference.h"
#include "../../Include/Codes/kloe_class.h"
#include "../../Include/const.h"
#include <TROOT.h>
#include <TTree.h>
#include <TFile.h>
#include <TChain.h>
#include <THStack.h>
#include <Math/Minimizer.h>
#include <Math/Factory.h>
#include <Math/Functor.h>
#include <TCanvas.h>
#include <TMath.h>

using namespace KLOE;

inline void chain_init(TChain *chain_init)
{
    TString fullname = "",
    dirnamemc = "MONTE_CARLO", dirnamedata = "DATA",
    filenamemc = "mc_stream62_mccard2_", filenamedata = "data_stream42_",
    extension = ".root";

    for(Int_t i = 1; i <= 56; i++)
    {
    	fullname = "/internal/big_one/4/users/gamrat/old_root_files/" + dirnamedata + "/" + filenamedata + i + extension;
    	chain_init->Add(fullname);
    }
    for(Int_t i = 1; i <= 56; i++)
    {
        if( i != 7)
        {
            fullname = "/internal/big_one/4/users/gamrat/old_root_files/" + dirnamemc + "/" + filenamemc + i + extension;
    	    chain_init->Add(fullname);
        }
    }
};

void cp_fit(TString mode = "")
{
    TFile file("../../../old_root_files/mctruth.root");
    TTree *tree = (TTree*)file.Get("h1");

    Int_t mctruth = 0;

    tree->SetBranchAddress("mctruth", &mctruth);

    TChain *chain = new TChain("INTERF/h1");

    chain_init(chain);

    UChar_t mcflag;
    Float_t Dtboostlor, Dtmc;

    chain->SetBranchAddress("mcflag", &mcflag);

    chain->SetBranchAddress("Dtboostlor", &Dtboostlor);
    chain->SetBranchAddress("Dtmc", &Dtmc);

    chain->AddFriend("h = h1", "../../../old_root_files/mctruth.root");

    UInt_t nentries = chain->GetEntries();
    UInt_t *entries = new UInt_t[chann_num];

    entries[0] = chain->GetEntries("h.mctruth == 1");
    entries[1] = chain->GetEntries("h.mctruth == 3");
    entries[2] = chain->GetEntries("h.mctruth == 4");
    entries[3] = chain->GetEntries("h.mctruth == 5");
    entries[4] = chain->GetEntries("h.mctruth == 6");
    entries[5] = chain->GetEntries("h.mctruth == 7");
    entries[6] = chain->GetEntries("h1.mcflag == 0");

    UInt_t nbins = 181;
    Double_t x_min = -90.0, x_max = 90.0;

    interference event(mode, entries, nbins, x_min, x_max);

    for(UInt_t i = 0; i < chann_num; i++) entries[i] = 0;

    for(UInt_t i = 0; i < nentries; i++)
    {
        chain->GetEntry(i);
        tree->GetEntry(i);

        if(mcflag == 1)
        {
            if(mctruth == 1)
            {
                event.time_diff_gen[entries[0]] = Dtmc;
                event.time_diff[0][entries[0]] = Dtboostlor;

                entries[0]++;
            }

            if(mctruth == 3)
            {
                event.time_diff[1][entries[1]] = Dtboostlor;

                entries[1]++;
            }

            if(mctruth == 4)
            {
                event.time_diff[2][entries[2]] = Dtboostlor;

                entries[2]++;
            }

            if(mctruth == 5)
            {
                event.time_diff[3][entries[3]] = Dtboostlor;

                entries[3]++;
            }

            if(mctruth == 6)
            {
                event.time_diff[4][entries[4]] = Dtboostlor;

                entries[4]++;
            }

            if(mctruth == 7)
            {
                event.time_diff[5][entries[5]] = Dtboostlor;

                entries[5]++;
            }
        }

        if(mcflag == 0)
        {
            event.time_diff[6][entries[6]] = Dtboostlor;

            entries[6]++;
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

    Double_t limit_span = 1.0;

    minimum->SetLimitedVariable(0, "Real part",      init_vars[0], step[0], init_vars[0] - limit_span*init_vars[0], init_vars[0] + limit_span*init_vars[0]);
    minimum->SetLimitedVariable(1, "Imaginary part", init_vars[1], step[1], init_vars[1] - limit_span*init_vars[1], init_vars[1] + limit_span*init_vars[1]);
    minimum->SetLimitedVariable(2, "Norm signal",    init_vars[2], step[2], init_vars[2] - limit_span*init_vars[2], init_vars[2] + limit_span*init_vars[2]);
    minimum->SetLimitedVariable(3, "Norm left DC wall",     init_vars[3], step[3], init_vars[3] - limit_span*init_vars[3], init_vars[3] + limit_span*init_vars[3]);
    minimum->SetLimitedVariable(4, "Norm left beam pipe",     init_vars[4], step[4], init_vars[4] - limit_span*init_vars[4], init_vars[4] + limit_span*init_vars[4]);
    minimum->SetLimitedVariable(5, "Norm right beam pipe",     init_vars[5], step[5], init_vars[5] - limit_span*init_vars[5], init_vars[5] + limit_span*init_vars[5]);
    minimum->SetLimitedVariable(6, "Norm right DC wall",     init_vars[6], step[6], init_vars[6] - limit_span*init_vars[6], init_vars[6] + limit_span*init_vars[6]);
    minimum->SetLimitedVariable(7, "Norm omega",     init_vars[7], step[7], init_vars[7] - limit_span*init_vars[7], init_vars[7] + limit_span*init_vars[7]);
    minimum->SetLimitedVariable(8, "Norm three",     init_vars[8], step[8], init_vars[8] - limit_span*init_vars[8], init_vars[8] + limit_span*init_vars[8]);
    minimum->SetLimitedVariable(9, "Norm semi",      init_vars[9], step[9], init_vars[9] - limit_span*init_vars[9], init_vars[9] + limit_span*init_vars[9]);
    minimum->SetLimitedVariable(10, "Norm other bcg", init_vars[10], step[10], init_vars[10] - limit_span*init_vars[10], init_vars[10] + limit_span*init_vars[10]);

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

    TH1 *histo[chann_num - 1];

    Double_t par[2] = {minimum->X()[0], minimum->X()[1]};

    for(UInt_t i = 0; i < chann_num - 1; i++)
    {
        histo[i] = new TH1F(("histogram" + std::to_string(i)).c_str(), "", nbins, x_min, x_max);
    }

    for(UInt_t i = 0; i < chann_num - 2; i++)
    {
        for(UInt_t j = 0; j < entries[i]; j++)
        {
            if(i == 0){
                histo[i]->Fill(event.time_diff[i][j], event.interf_function(event.time_diff_gen[j]));
            }
            else if(i == 1)
            {  
                if(event.time_diff[i][j] < -30.0)
                    histo[i]->Fill(event.time_diff[i][j], minimum->X()[3]);
                else if(event.time_diff[i][j] > -30.0 && event.time_diff[i][j] < 0.0)
                    histo[i]->Fill(event.time_diff[i][j], minimum->X()[4]);
                else if(event.time_diff[i][j] > 0.0 && event.time_diff[i][j] < 30.0)
                    histo[i]->Fill(event.time_diff[i][j], minimum->X()[5]);
                else if(event.time_diff[i][j] > 30.0)
                    histo[i]->Fill(event.time_diff[i][j], minimum->X()[6]);
            }
            else
            {
                histo[i]->Fill(event.time_diff[i][j]);
            }
        }
    }

    histo[0]->Scale(minimum->X()[2]*histo[0]->GetEntries() / histo[0]->Integral(0, nbins + 1) );

    histo[2]->Scale(minimum->X()[7]*histo[2]->GetEntries() / histo[2]->Integral(0, nbins + 1) );
    histo[3]->Scale(minimum->X()[8]*histo[3]->GetEntries() / histo[3]->Integral(0, nbins + 1) );
    histo[4]->Scale(minimum->X()[9]*histo[4]->GetEntries() / histo[4]->Integral(0, nbins + 1) );
    histo[5]->Scale(minimum->X()[10]*histo[5]->GetEntries() / histo[5]->Integral(0, nbins + 1) );

    histo[chann_num - 2]->Add(histo[0]);
    histo[chann_num - 2]->Add(histo[1]);
    histo[chann_num - 2]->Add(histo[2]);
    histo[chann_num - 2]->Add(histo[3]);
    histo[chann_num - 2]->Add(histo[4]);
    histo[chann_num - 2]->Add(histo[5]);

    histo[0]->SetLineColor(kRed);
    histo[1]->SetLineColor(kGreen);
    histo[2]->SetLineColor(kViolet);
    histo[3]->SetLineColor(kCyan);
    histo[4]->SetLineColor(kBlue);
    histo[5]->SetLineColor(kGreen - 1);
    histo[6]->SetLineColor(kBlack);
    histo[7]->SetLineColor(kGray - 1);

    TCanvas *c1 = new TCanvas("c1", "", 790, 790);

    histo[7]->Draw("HIST");
    histo[6]->Draw("PE1SAME");
    histo[0]->Draw("HISTSAME");
    histo[1]->Draw("HISTSAME");
    histo[2]->Draw("HISTSAME");
    histo[3]->Draw("HISTSAME");
    histo[4]->Draw("HISTSAME");
    histo[5]->Draw("HISTSAME");


    c1->Print("dupa.png");


}