#include "interference.h"
#include "kloe_class.h"
#include "../../Include/const.h"
#include <TROOT.h>
#include <Math/Minimizer.h>
#include <Math/Factory.h>
#include <Math/Functor.h>

using namespace KLOE;

void cp_fit(TString mode = "")
{
    gROOT->ProcessLine(".x ../load_sel.C");

    interference event(mode);

    UInt_t mcflag, mctruth;
    Float_t Dtboostlor, Dtmc;

    chain.SetBranchAddress("mcflag", &mcflag);
    chain.SetBranchAddress("mctruth", &mctruth);

    chain.SetBranchAddress("Dtboostlor", &Dtboostlor);
    chain.SetBranchAddress("Dtmc", &Dtmc);

    UInt_t nentries = chain.GetEntries(), entries[chann_num] = {0};

    TH1 *histos[chann_num];

    Double_t *time_diff[chann_num], *time_diff_gen;

    for(Int_t i = 0; i < nentries; i++)
    {
        chain.GetEntry(i);

        if(mcflag == 1)
        {
            if(mctruth == 1)
            {
                time_diff_gen[entries[0]];
                time_diff[0][entries[0]];
                entries[0]++;
            }

            if(mctruth == 3)
            {
                time_diff[1][entries[1]];
                entries[1]++;
            }

            if(mctruth == 4)
            {
                time_diff[2][entries[2]];
                entries[2]++;
            }

            if(mctruth == 5)
            {
                time_diff[3][entries[3]];
                entries[3]++;
            }

            if(mctruth == 6)
            {
                time_diff[4][entries[4]];
                entries[4]++;
            }

            if(mctruth == 7)
            {
                time_diff[5][entries[5]];
                entries[5]++;
            }
        }

        if(mcflag == 0)
        {
            time_diff[6][entries[6]];
            entries[6]++;
        }

    }

    ROOT::Math::Minimizer* minimum = 
        ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");

    // set tolerance , etc...
    minimum->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2
    minimum->SetTolerance(0.001);
    minimum->SetPrintLevel(1);
 
    ROOT::Math::Functor minimized_function(interference::interf_chi2, 2);
 
    minimum->SetFunction(minimized_function);

    for(Int_t i = 0; i < 7; i++)
        for(Int_t j = 0; j < entries[i]; j++)
        {
            minimum->SetVariable(i*entries[i] + entries[i],"x",variable[0], step[0]);
        }


}