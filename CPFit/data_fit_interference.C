#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_vector.h>
#include <algorithm>
#include <initializer_list>
#include <iostream>
#include <string>
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TTree.h"
#include "TStyle.h"
#include "TGraphAsymmErrors.h"
#include "TFile.h"
#include "chisquare.h"
#include "TLine.h"
#include "TLegend.h"
#include "TF1.h"
#include "TPaveText.h"
#include "TGaxis.h"
#include "TMultiGraph.h"
#include "fits.h"
#include "my_methods.cpp"

using namespace std;

Double_t funkcja(Float_t dt)
{
    Double_t Value = 0;
    Double_t Epsilon = 0, RePart = 0, Dphi = 0, TauKs = 0, TauKl = 0, MassDiff = 0;
    Double_t ImPart = 0, GammaKl = 0, GammaKs = 0, Gamma = 0, DMass = 0;

    // Parameters from PDG2020

    Epsilon = 2.228E-3;
    RePart = 1.66E-3;
    Dphi = -0.34;         // phi(+-)-phi(00) (degrees)
    TauKs = 0.8954E-10;   // PDG fit assuming CPT (s)
    TauKl = 5.116E-8;     // Kl mean life (s)
    MassDiff = 0.5293E10; // M(Kl)-M(Ks) ( (h/2pi)s-1 ):
                          // PDG fit assuming CPT

    // All parameters are calculated taking Int_to account that DT is in TauKs units

    ImPart = (-0.34 * M_PI / 180.) / 3.; // Im(epsilon'/epsilon) = Dphi/3
    GammaKs = 1.;
    GammaKl = TauKs / TauKl;
    Gamma = GammaKs + GammaKl;
    DMass = MassDiff * TauKs;

    if (dt >= 0.)
    {
        Value = (1. + 2. * RePart) * exp(-GammaKl * dt) +
                (1. - 4. * RePart) * exp(-GammaKs * dt) -
                2. * exp(-0.5 * Gamma * dt) *
                    ((1. - RePart) * cos(DMass * dt) +
                     3. * ImPart * sin(DMass * dt));
    }
    else
    {
        Value = (1. + 2. * RePart) * exp(-GammaKs * abs(dt)) +
                (1. - 4. * RePart) * exp(-GammaKl * abs(dt)) -
                2. * exp(-Gamma / 2. * abs(dt)) *
                    ((1. - RePart) * cos(DMass * abs(dt)) -
                     3. * ImPart * sin(DMass * abs(dt)));
    }

    return (Epsilon * Epsilon) / (2. * Gamma) * Value * 100000;
}

Double_t fitf(const Float_t x, const Double_t *par)
{
    Double_t Value = 0;
    Double_t Epsilon = 0, RePart = 0, Dphi = 0, TauKs = 0, TauKl = 0, MassDiff = 0;
    Double_t ImPart = 0, GammaKl = 0, GammaKs = 0, Gamma = 0, DMass = 0;

    Float_t dt = x;
    Int_t howmuchregen = par[0];
    Int_t howmuchomega = par[1];
    Int_t howmuchthree = par[2];
    Int_t howmuchsemi = par[3];
    Int_t howmuchelsee = par[4];
    Int_t binn = par[5];
    Int_t howmuch = par[6];
    // Parameters from PDG2020

    Epsilon = 2.228E-3;
    RePart = par[2 * howmuch + 2 * binn + howmuchregen + howmuchomega + howmuchthree + howmuchsemi + howmuchelsee + 7];
    Dphi = -0.34;         // phi(+-)-phi(00) (degrees)
    TauKs = 0.8954E-10;   // PDG fit assuming CPT (s)
    TauKl = 5.116E-8;     // Kl mean life (s)
    MassDiff = 0.5293E10; // M(Kl)-M(Ks) ( (h/2pi)s-1 ):
                          // PDG fit assuming CPT

    // All parameters are calculated taking Int_to account that DT is in TauKs units

    ImPart = par[2 * howmuch + 8 + 2 * binn + howmuchregen + howmuchomega + howmuchthree + howmuchsemi + howmuchelsee]; // Im(epsilon'/epsilon) = Dphi/3
    GammaKs = 1.;
    GammaKl = TauKs / TauKl;
    Gamma = GammaKs + GammaKl;
    DMass = MassDiff * TauKs;

    if (dt >= 0.)
    {
        Value = (1. + 2. * RePart) * exp(-GammaKl * dt) +
                (1. - 4. * RePart) * exp(-GammaKs * dt) -
                2. * exp(-0.5 * Gamma * dt) *
                    ((1. - RePart) * cos(DMass * dt) +
                     3. * ImPart * sin(DMass * dt));
    }
    else
    {
        Value = (1. + 2. * RePart) * exp(-GammaKs * abs(dt)) +
                (1. - 4. * RePart) * exp(-GammaKl * abs(dt)) -
                2. * exp(-Gamma / 2. * abs(dt)) *
                    ((1. - RePart) * cos(DMass * abs(dt)) -
                     3. * ImPart * sin(DMass * abs(dt)));
    }

    return ((((Epsilon * Epsilon) / (2. * Gamma)) * Value * 100000));
}

Double_t fitowana(const Double_t *xx)
{

    TTree *t1 = new TTree("t1", "a treet1");
    t1->ReadFile("whole_efficiency.txt", "wydajnosc/D");
    Double_t wydaj = 0;

    t1->SetBranchAddress("wydajnosc", &wydaj);

    Int_t nentriesgatwydaj = (Int_t)t1->GetEntries();

    Double_t *wydajtab = new Double_t[nentriesgatwydaj];

    for (Int_t i = 0; i < nentriesgatwydaj; i++)
    {
        t1->GetEntry(i);
        wydajtab[i] = wydaj;
    }

    Int_t howmuchregen = xx[0];
    Int_t howmuchomega = xx[1];
    Int_t howmuchthree = xx[2];
    Int_t howmuchsemi = xx[3];
    Int_t howmuchelsee = xx[4];
    Int_t binn = xx[5];
    Int_t howmuch = xx[6];
    Double_t N_signal = xx[2 * howmuch + 2 * binn + howmuchregen + howmuchsemi + howmuchomega + howmuchthree + howmuchelsee + 9];
    Double_t N_regen = xx[2 * howmuch + 2 * binn + howmuchregen + howmuchsemi + howmuchomega + howmuchthree + howmuchelsee + 11];
    Double_t N_omega = xx[2 * howmuch + 2 * binn + howmuchregen + howmuchsemi + howmuchomega + howmuchthree + howmuchelsee + 12];
    Double_t N_three = xx[2 * howmuch + 2 * binn + howmuchregen + howmuchsemi + howmuchomega + howmuchthree + howmuchelsee + 12];
    Double_t N_semi = xx[2 * howmuch + 2 * binn + howmuchregen + howmuchsemi + howmuchomega + howmuchthree + howmuchelsee + 12];
    Double_t N_elsee = xx[2 * howmuch + 2 * binn + howmuchregen + howmuchsemi + howmuchomega + howmuchthree + howmuchelsee + 12];
    Double_t N_regen_left_far = xx[2 * howmuch + 2 * binn + howmuchregen + howmuchsemi + howmuchomega + howmuchthree + howmuchelsee + 15];
    Double_t N_regen_left_close = xx[2 * howmuch + 2 * binn + howmuchregen + howmuchsemi + howmuchomega + howmuchthree + howmuchelsee + 16];
    Double_t N_regen_right_far = xx[2 * howmuch + 2 * binn + howmuchregen + howmuchsemi + howmuchomega + howmuchthree + howmuchelsee + 15];
    Double_t N_regen_right_close = xx[2 * howmuch + 2 * binn + howmuchregen + howmuchsemi + howmuchomega + howmuchthree + howmuchelsee + 16];

    Double_t *b = new Double_t[binn];
    Double_t *e = new Double_t[binn];
    Double_t *b0 = new Double_t[binn];
    Double_t *e0 = new Double_t[binn];

    Double_t *bregen = new Double_t[binn];
    Double_t *eregen = new Double_t[binn];

    Double_t *bomega = new Double_t[binn];
    Double_t *eomega = new Double_t[binn];

    Double_t *bthree = new Double_t[binn];
    Double_t *ethree = new Double_t[binn];

    Double_t *bsemi = new Double_t[binn];
    Double_t *esemi = new Double_t[binn];

    Double_t *belsee = new Double_t[binn];
    Double_t *eelsee = new Double_t[binn];

    TH1 *h1 = new TH1F("Sig", "", binn, -90.0, 90.0);

    for (Int_t i = 7; i < howmuch + 7; i++)
        h1->Fill(xx[i], fitf(xx[i + howmuch], xx));

    h1->Scale(71532 / h1->Integral(1, binn));

    TH1 *h1regen = new TH1F("Regen", "", binn, -90.0, 90.0);

    for (Int_t i = 0; i < howmuchregen; i++)
        h1regen->Fill(xx[i + 2 * howmuch + 2 * binn + 7], N_regen);
    for (Int_t i = 0; i < binn; i++)
        bregen[i] = (h1regen->GetBinContent(i + 1));
    for (Int_t i = 0; i < binn; i++)
        eregen[i] = (h1regen->GetBinError(i + 1));

    TH1 *h1omega = new TH1F("Omega", "", binn, -90.0, 90.0);

    for (Int_t i = 0; i < howmuchomega; i++)
        h1omega->Fill(xx[i + 2 * howmuch + 2 * binn + 7 + howmuchregen], N_omega);
    for (Int_t i = 0; i < binn; i++)
        bomega[i] = (h1omega->GetBinContent(i + 1));
    for (Int_t i = 0; i < binn; i++)
        eomega[i] = (h1omega->GetBinError(i + 1));

    TH1 *h1three = new TH1F("Thre", "", binn, -90.0, 90.0);

    for (Int_t i = 0; i < howmuchthree; i++)
        h1three->Fill(xx[i + 2 * howmuch + 2 * binn + 7 + howmuchregen + howmuchomega], N_three);
    for (Int_t i = 0; i < binn; i++)
        bthree[i] = (h1three->GetBinContent(i + 1));
    for (Int_t i = 0; i < binn; i++)
        ethree[i] = (h1three->GetBinError(i + 1));
    TH1 *h1semi = new TH1F("Semi", "", binn, -90.0, 90.0);

    for (Int_t i = 0; i < howmuchsemi; i++)
        h1semi->Fill(xx[i + 2 * howmuch + 2 * binn + 7 + howmuchregen + howmuchomega + howmuchthree], N_semi);
    for (Int_t i = 0; i < binn; i++)
        bsemi[i] = (h1semi->GetBinContent(i + 1));
    for (Int_t i = 0; i < binn; i++)
        esemi[i] = (h1semi->GetBinError(i + 1));
    TH1 *h1elsee = new TH1F("elsee", "", binn, -90.0, 90.0);

    for (Int_t i = 0; i < howmuchelsee; i++)
        h1elsee->Fill(xx[i + 2 * howmuch + 2 * binn + 7 + howmuchregen + howmuchomega + howmuchthree + howmuchsemi], N_elsee);
    for (Int_t i = 0; i < binn; i++)
        belsee[i] = (h1elsee->GetBinContent(i + 1));
    for (Int_t i = 0; i < binn; i++)
        eelsee[i] = (h1elsee->GetBinError(i + 1));
    Double_t value = 0;

    for (Int_t i = 0; i < binn; i++)
        b0[i] = xx[2 * howmuch + 7 + i];
    for (Int_t i = 0; i < binn; i++)
        e0[i] = xx[2 * howmuch + 7 + binn + i];

    for (Int_t i = 0; i < binn; i++)
        b[i] = (h1->GetBinContent(i + 1));
    for (Int_t i = 0; i < binn; i++)
        e[i] = (h1->GetBinError(i + 1));

    Int_t far_left_regen_bin_end = h1->FindBin(-30.0);

    Int_t close_left_regen_bin_start = h1->FindBin(-30.0);
    Int_t close_left_regen_bin_end = h1->FindBin(0.0);

    Int_t close_right_regen_bin_start = h1->FindBin(0.0);
    Int_t close_right_regen_bin_end = h1->FindBin(30.0);

    Int_t far_right_regen_bin_start = h1->FindBin(30.0);

    for (Int_t i = 0; i < far_left_regen_bin_end; i++)
        value += pow(b0[i] - N_signal * b[i] - N_regen_left_far * bregen[i] - N_omega * bomega[i] -
                         N_three * bthree[i] - N_semi * bsemi[i] - N_elsee * belsee[i],
                     2) /
                 (pow(e0[i], 2) +
                  pow(N_signal * e[i], 2) + pow(N_regen_left_far, 2) * bregen[i] + pow(N_omega, 2) * bomega[i] +
                  pow(N_three, 2) * bthree[i] + pow(N_semi, 2) * bsemi[i] + pow(N_elsee, 2) * belsee[i]);

    for (Int_t i = close_left_regen_bin_start; i < close_left_regen_bin_end; i++)
        value += pow(b0[i] - N_signal * b[i] - N_regen_left_close * bregen[i] - N_omega * bomega[i] -
                         N_three * bthree[i] - N_semi * bsemi[i] - N_elsee * belsee[i],
                     2) /
                 (pow(e0[i], 2) +
                  pow(N_signal * e[i], 2) + pow(N_regen_left_close, 2) * bregen[i] + pow(N_omega, 2) * bomega[i] +
                  pow(N_three, 2) * bthree[i] + pow(N_semi, 2) * bsemi[i] + pow(N_elsee, 2) * belsee[i]);

    for (Int_t i = close_right_regen_bin_start; i < close_right_regen_bin_end; i++)
        value += pow(b0[i] - N_signal * b[i] - N_regen_right_close * bregen[i] - N_omega * bomega[i] -
                         N_three * bthree[i] - N_semi * bsemi[i] - N_elsee * belsee[i],
                     2) /
                 (pow(e0[i], 2) +
                  pow(N_signal * e[i], 2) + pow(N_regen_right_close, 2) * bregen[i] + pow(N_omega, 2) * bomega[i] +
                  pow(N_three, 2) * bthree[i] + pow(N_semi, 2) * bsemi[i] + pow(N_elsee, 2) * belsee[i]);

    for (Int_t i = far_right_regen_bin_start; i < binn; i++)
        value += pow(b0[i] - N_signal * b[i] - N_regen_right_far * bregen[i] - N_omega * bomega[i] -
                         N_three * bthree[i] - N_semi * bsemi[i] - N_elsee * belsee[i],
                     2) /
                 (pow(e0[i], 2) +
                  pow(N_signal * e[i], 2) + pow(N_regen_right_far, 2) * bregen[i] + pow(N_omega, 2) * bomega[i] +
                  pow(N_three, 2) * bthree[i] + pow(N_semi, 2) * bsemi[i] + pow(N_elsee, 2) * belsee[i]);

    h1->Delete();
    h1omega->Delete();
    h1semi->Delete();
    h1three->Delete();
    h1regen->Delete();
    h1elsee->Delete();

    return value;
}

Int_t recon_mc_fit_eff(Int_t minpoints = 1, Int_t maxpoints = 1, TString name = "signal_pas.root")
{

    TTree *t = new TTree("t", "a tree");
    t->ReadFile("whole_efficiency.txt", "wydajnosc/D");

    Double_t wydaj = 0;

    t->SetBranchAddress("wydajnosc", &wydaj);

    Int_t nentriesgatwydaj = (Int_t)t->GetEntries();

    Double_t *wydajtab = new Double_t[nentriesgatwydaj];

    for (Int_t i = 0; i < nentriesgatwydaj; i++)
    {
        t->GetEntry(i);
        wydajtab[i] = wydaj;
    }

    // cout << wydajtab[0]

    TCanvas *c1 = new TCanvas("c1", "c1", 1000, 760);
    Double_t a = c1->GetWw();

    Double_t leftmargin = 100 / a;
    Double_t bottommargin = 150 / a;

    c1->SetLeftMargin(leftmargin);
    c1->SetBottomMargin(bottommargin);
    c1->SetRightMargin(leftmargin);
    c1->SetTopMargin(bottommargin);

    TCanvas *c2 = new TCanvas("c2", "c2", 1000, 760);

    c2->SetLeftMargin(leftmargin);
    c2->SetBottomMargin(bottommargin);
    c2->SetRightMargin(leftmargin);
    c2->SetTopMargin(bottommargin);

    // Double_t xRe[100] = {0}; Double_t xIm[100] = {0};
    // Double_t yRe[100] = {0}; Double_t yIm[100] = {0};
    Int_t n = 36;
    Int_t k = 0;

    Double_t down = -90.0;
    Double_t up = 90.0;

    Int_t howmuch = 0;

    Int_t binn = floor(up - down);

    TH1 *signal_after_fit = new TH1F("Signal_after_fit", "", binn, down, up);
    TH1 *mimicized_data = new TH1F("Mimicized_data", "", binn, down, up);
    TH1 *signal_before_fit = new TH1F("Signal_before_fit", "", binn, down, up);

    TH1 *omega_start = new TH1F("Omega_start", "", binn, down, up);
    TH1 *regen_start = new TH1F("Regen_start", "", binn, down, up);
    TH1 *three_start = new TH1F("Three_start", "", binn, down, up);
    TH1 *semi_start = new TH1F("Semi_start", "", binn, down, up);
    TH1 *elsee_start = new TH1F("elsee_start", "", binn, down, up);

    TH1 *omega_after_fit = new TH1F("Omega_after_fit", "", binn, down, up);
    TH1 *regen_after_fit = new TH1F("Regen_after_fit", "", binn, down, up);
    TH1 *three_after_fit = new TH1F("Three_after_fit", "", binn, down, up);
    TH1 *semi_after_fit = new TH1F("Semi_after_fit", "", binn, down, up);
    TH1 *elsee_after_fit = new TH1F("elsee_after_fit", "", binn, down, up);

    TFile *mc_filesemi = new TFile("semileptonic.root");
    TTree *semi = (TTree *)mc_filesemi->Get("h1");

    TFile *mc_filethree = new TFile("threepions.root");
    TTree *three = (TTree *)mc_filethree->Get("h1");

    TFile *mc_fileregen = new TFile("regen.root");
    TTree *regen = (TTree *)mc_fileregen->Get("h1");

    TFile *mc_fileomega = new TFile("omega.root");
    TTree *omega = (TTree *)mc_fileomega->Get("h1");

    TFile *mc_file2 = new TFile(name);
    TTree *gathered = (TTree *)mc_file2->Get("h1;1");

    TFile *mc_filedata = new TFile("whole_data.root");
    TTree *data = (TTree *)mc_filedata->Get("h1");

    TFile *mc_fileelsee = new TFile("else.root");
    TTree *elsee = (TTree *)mc_fileelsee->Get("h1");

    Float_t Dtboostlor = 0;
    Float_t Dtmc = 0;
    gathered->SetBranchAddress("Dtboostlor", &Dtboostlor);
    gathered->SetBranchAddress("Dtmc", &Dtmc);

    Float_t Dtboostlorregen = 0;
    regen->SetBranchAddress("Dtboostlor", &Dtboostlorregen);

    Float_t Dtboostloromega = 0;
    omega->SetBranchAddress("Dtboostlor", &Dtboostloromega);

    Float_t Dtboostlorthree = 0;
    three->SetBranchAddress("Dtboostlor", &Dtboostlorthree);

    Float_t Dtboostlorsemi = 0;
    semi->SetBranchAddress("Dtboostlor", &Dtboostlorsemi);

    Float_t Dtboostlordata = 0;
    data->SetBranchAddress("Dtboostlor", &Dtboostlordata);

    Float_t Dtboostlorelsee = 0;
    elsee->SetBranchAddress("Dtboostlor", &Dtboostlorelsee);

    Int_t nentriesgat = (Int_t)gathered->GetEntries();

    Int_t nentriesregen = (Int_t)regen->GetEntries();

    Int_t nentriesomega = (Int_t)omega->GetEntries();

    Int_t nentriesthree = (Int_t)three->GetEntries();

    Int_t nentriessemi = (Int_t)semi->GetEntries();

    Int_t nentriesdata = (Int_t)data->GetEntries();

    Int_t nentrieselsee = (Int_t)elsee->GetEntries();

    Int_t howmuchregen = floor(nentriesregen);
    Int_t howmuchsemi = floor(nentriessemi);
    Int_t howmuchthree = floor(nentriesthree);
    Int_t howmuchomega = floor(nentriesomega);
    Int_t howmuchelsee = floor(nentrieselsee);

    Double_t errup = 0;
    Double_t errdown = 0;
    Int_t nu = 0;

    ROOT::Math::Minimizer *minimum =
        ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");

    minimum->SetMaxFunctionCalls(10000000); // for Minuit/Minuit2
    minimum->SetTolerance(0.5);
    minimum->SetPrintLevel(3);
    // minimum->SetPrecision(1.0E-12);
    minimum->SetStrategy(2);

    Double_t xRe[20] = {0};
    Double_t xIm[20] = {0};
    Double_t yRe[20] = {0};
    Double_t yIm[20] = {0};

    Double_t xReeff[20] = {0};
    Double_t xImeff[20] = {0};
    Double_t yReeff[20] = {0};
    Double_t yImeff[20] = {0};

    //////////////////////////////////////////////////////////////////////////////////////////
    Double_t *dtmc = new Double_t[nentriesgat];
    Double_t *dtrec = new Double_t[nentriesgat];

    Double_t *dtregen = new Double_t[nentriesregen];

    Double_t *dtomega = new Double_t[nentriesomega];

    Double_t *dtthree = new Double_t[nentriesthree];

    Double_t *dtsemi = new Double_t[nentriessemi];

    Double_t *dtelsee = new Double_t[nentrieselsee];
    //////////////////////////////////////////////////////////////////////////////////////////
    Double_t *b0 = new Double_t[binn];
    Double_t *e0 = new Double_t[binn];
    Double_t *b00 = new Double_t[binn];
    Double_t *e00 = new Double_t[binn];

    Double_t *bregen = new Double_t[binn];
    Double_t *eregen = new Double_t[binn];
    Double_t *bomega = new Double_t[binn];
    Double_t *eomega = new Double_t[binn];
    Double_t *bthree = new Double_t[binn];
    Double_t *ethree = new Double_t[binn];
    Double_t *bsemi = new Double_t[binn];
    Double_t *esemi = new Double_t[binn];
    Double_t *belsee = new Double_t[binn];
    Double_t *eelsee = new Double_t[binn];
    ///////////////////////////////////////////////////////////////////////////////////////////
    for (Int_t i = 0; i < nentriesgat; i++)
    {
        gathered->GetEntry(i);

        dtmc[i] = Dtmc;
        dtrec[i] = Dtboostlor;
    }

    for (Int_t i = 0; i < nentriesregen; i++)
    {
        regen->GetEntry(i);

        dtregen[i] = Dtboostlorregen;
    }

    for (Int_t i = 0; i < nentriesomega; i++)
    {
        omega->GetEntry(i);

        dtomega[i] = Dtboostloromega;
    }

    for (Int_t i = 0; i < nentriesthree; i++)
    {
        three->GetEntry(i);

        dtthree[i] = Dtboostlorthree;
    }

    for (Int_t i = 0; i < nentriessemi; i++)
    {
        semi->GetEntry(i);

        dtsemi[i] = Dtboostlorsemi;
    }

    for (Int_t i = 0; i < nentrieselsee; i++)
    {
        elsee->GetEntry(i);

        dtelsee[i] = Dtboostlorelsee;
    }
    //////////////////////////////////////////////////////////////////////////////////////////i

    for (Int_t i = 0; i < nentriesdata; i++)
    {
        data->GetEntry(i);

        mimicized_data->Fill(Dtboostlordata);
    }

    // Double_t factor = mimicized_data->Integral(1,binn);

    for (Int_t i = 0; i < floor(howmuchregen); i++)
    {
        regen_start->Fill(dtregen[i]);
    }

    for (Int_t i = 0; i < floor(howmuchomega); i++)
    {
        omega_start->Fill(dtomega[i]);
    }

    for (Int_t i = 0; i < floor(howmuchthree); i++)
    {
        three_start->Fill(dtthree[i]);
    }

    for (Int_t i = 0; i < floor(howmuchsemi); i++)
    {
        semi_start->Fill(dtsemi[i]);
    }

    for (Int_t i = 0; i < floor(howmuchelsee); i++)
    {
        elsee_start->Fill(dtelsee[i]);
    }
    /*mimicized_data->Add(regen_start);
    mimicized_data->Add(omega_start);
    mimicized_data->Add(three_start);
    mimicized_data->Add(semi_start);*/

    for (Int_t i = 0; i < binn; i++)
        b0[i] = (mimicized_data->GetBinContent(i + 1));
    for (Int_t i = 0; i < binn; i++)
        e0[i] = (mimicized_data->GetBinError(i + 1));

    for (Int_t i = 0; i < binn; i++)
        b00[i] = (mimicized_data->GetBinContent(i + 1));
    for (Int_t i = 0; i < binn; i++)
        e00[i] = (mimicized_data->GetBinError(i + 1));

    for (Int_t i = 0; i < binn; i++)
        bregen[i] = (regen_start->GetBinContent(i + 1));
    for (Int_t i = 0; i < binn; i++)
        eregen[i] = (regen_start->GetBinError(i + 1));

    for (Int_t i = 0; i < binn; i++)
        bomega[i] = (omega_start->GetBinContent(i + 1));
    for (Int_t i = 0; i < binn; i++)
        eomega[i] = (omega_start->GetBinError(i + 1));

    for (Int_t i = 0; i < binn; i++)
        bthree[i] = (three_start->GetBinContent(i + 1));
    for (Int_t i = 0; i < binn; i++)
        ethree[i] = (three_start->GetBinError(i + 1));

    for (Int_t i = 0; i < binn; i++)
        bsemi[i] = (semi_start->GetBinContent(i + 1));
    for (Int_t i = 0; i < binn; i++)
        esemi[i] = (semi_start->GetBinError(i + 1));

    for (Int_t i = 0; i < binn; i++)
        belsee[i] = (elsee_start->GetBinContent(i + 1));
    for (Int_t i = 0; i < binn; i++)
        eelsee[i] = (elsee_start->GetBinError(i + 1));

    for (Int_t i = minpoints; i <= maxpoints; i++)
    {
        howmuch = nentriesgat;

        for (Int_t i = 0; i < howmuch; i++)
        {
            signal_before_fit->Fill(dtrec[i]);
        }

        ROOT::Math::Functor f3(&fitowana, 2 * howmuch + 19 + 2 * binn + howmuchregen + howmuchsemi + howmuchthree + howmuchomega + howmuchelsee);
        minimum->SetFunction(f3);

        minimum->SetFixedVariable(0, "Nentriesregen", howmuchregen);
        minimum->SetFixedVariable(1, "Nentriesomega", howmuchomega);
        minimum->SetFixedVariable(2, "Nentriesthree", howmuchthree);
        minimum->SetFixedVariable(3, "Nentriessemi", howmuchsemi);
        minimum->SetFixedVariable(4, "Nentrieselsee", howmuchelsee);
        minimum->SetFixedVariable(5, "Binn", binn);
        minimum->SetFixedVariable(6, "Howmuch", howmuch);

        for (Int_t i = 7; i < howmuch + 7; i++)
        {
            std::string(name);

            std::string(t) = std::to_string(i);
            name = "Dtrec" + t;

            minimum->SetFixedVariable(i, name.c_str(), dtrec[i - 7]);
        }

        for (Int_t i = 7; i < howmuch + 7; i++)
        {
            std::string(name1);

            std::string(t) = std::to_string(i);
            name1 = "Dtmc" + t;

            minimum->SetFixedVariable(i + howmuch, name1.c_str(), dtmc[i - 7]);
        }

        for (Int_t i = 0; i < binn; i++)
        {
            std::string name2;
            std::string(t) = std::to_string(i);

            name2 = ("Bin") + t;
            minimum->SetFixedVariable(2 * howmuch + 7 + i, name2.c_str(), b00[i]);
        }

        for (Int_t i = 0; i < binn; i++)
        {
            std::string name3;
            std::string(t) = std::to_string(i);

            name3 = ("Error") + t;
            minimum->SetFixedVariable(2 * howmuch + 7 + binn + i, name3.c_str(), e00[i]);
        }

        for (Int_t i = 0; i < howmuchregen; i++)
        {
            std::string name4;
            std::string(t) = std::to_string(i);

            name4 = ("Dtregen") + t;
            minimum->SetFixedVariable(2 * howmuch + 7 + 2 * binn + i, name4.c_str(), dtregen[i]);
        }

        for (Int_t i = 0; i < howmuchomega; i++)
        {
            std::string name5;
            std::string(t) = std::to_string(i);

            name5 = ("Dtomega") + t;
            minimum->SetFixedVariable(2 * howmuch + 7 + 2 * binn + howmuchregen + i, name5.c_str(), dtomega[i]);
        }

        for (Int_t i = 0; i < howmuchthree; i++)
        {
            std::string name6;
            std::string(t) = std::to_string(i);

            name6 = ("Dtthree") + t;
            minimum->SetFixedVariable(2 * howmuch + 7 + 2 * binn + howmuchregen + howmuchomega + i, name6.c_str(), dtthree[i]);
        }

        for (Int_t i = 0; i < howmuchsemi; i++)
        {
            std::string name7;
            std::string(t) = std::to_string(i);

            name7 = ("Dtsemi") + t;
            minimum->SetFixedVariable(2 * howmuch + 7 + 2 * binn + howmuchregen + howmuchomega + howmuchthree + i, name7.c_str(), dtsemi[i]);
        }

        for (Int_t i = 0; i < howmuchelsee; i++)
        {
            std::string name8;
            std::string(t) = std::to_string(i);

            name8 = ("Dtelsee") + t;
            minimum->SetFixedVariable(2 * howmuch + 7 + 2 * binn + howmuchregen + howmuchomega + howmuchthree + howmuchsemi + i, name8.c_str(), dtelsee[i]);
        }
        Double_t step[12] = {1.0E-2, 1.0E-2, 0.1, 0.0, 0.0, 0.1, 0.0, 0.0, 0.1, 0.1, 0.0, 0.0};

        Double_t variable[12] = {1.66E-3, -1.98E-3, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
        minimum->SetLimitedVariable(2 * howmuch + 2 * binn + howmuchregen + howmuchomega + howmuchthree + howmuchsemi + howmuchelsee + 7, "Re", variable[0], step[0], variable[0] - 5. * abs(variable[0]), variable[0] + 5. * abs(variable[0]));
        minimum->SetLimitedVariable(2 * howmuch + 2 * binn + howmuchregen + howmuchomega + howmuchthree + howmuchsemi + howmuchelsee + 8, "Im", variable[1], step[1], variable[1] - 5. * abs(variable[1]), variable[1] + 5. * abs(variable[1]));
        minimum->SetLimitedVariable(2 * howmuch + 2 * binn + howmuchregen + howmuchomega + howmuchthree + howmuchsemi + howmuchelsee + 9, "A", variable[2], step[2], 0.8, 1.2);

        minimum->SetLimitedVariable(2 * howmuch + 2 * binn + howmuchregen + howmuchomega + howmuchthree + howmuchsemi + howmuchelsee + 10, "N_regen", variable[3], step[3], 0.8, 1.2);
        minimum->SetLimitedVariable(2 * howmuch + 2 * binn + howmuchregen + howmuchomega + howmuchthree + howmuchsemi + howmuchelsee + 11, "N_omega", variable[4], step[4], 0.8, 1.2);
        minimum->SetLimitedVariable(2 * howmuch + 2 * binn + howmuchregen + howmuchomega + howmuchthree + howmuchsemi + howmuchelsee + 12, "N_three", variable[5], step[5], 0.8, 1.2);
        minimum->SetLimitedVariable(2 * howmuch + 2 * binn + howmuchregen + howmuchomega + howmuchthree + howmuchsemi + howmuchelsee + 13, "N_semi", variable[6], step[6], 0.8, 1.2);
        minimum->SetLimitedVariable(2 * howmuch + 2 * binn + howmuchregen + howmuchomega + howmuchthree + howmuchsemi + howmuchelsee + 14, "N_elsee", variable[7], step[7], 0.8, 1.2);
        minimum->SetLimitedVariable(2 * howmuch + 2 * binn + howmuchregen + howmuchomega + howmuchthree + howmuchsemi + howmuchelsee + 15, "N_regen_left_far", variable[8], step[8], 0.0, 50.0);
        minimum->SetLimitedVariable(2 * howmuch + 2 * binn + howmuchregen + howmuchomega + howmuchthree + howmuchsemi + howmuchelsee + 16, "N_regen_left_close", variable[9], step[9], 0.0, 50.0);
        minimum->SetLimitedVariable(2 * howmuch + 2 * binn + howmuchregen + howmuchomega + howmuchthree + howmuchsemi + howmuchelsee + 17, "N_regen_right_far", variable[10], step[10], 0.8, 1.2);
        minimum->SetLimitedVariable(2 * howmuch + 2 * binn + howmuchregen + howmuchomega + howmuchthree + howmuchsemi + howmuchelsee + 18, "N_regen_right_close", variable[11], step[11], 0.8, 1.2);

        unsigned int nstep = 500;
        Double_t x[500] = {};
        Double_t y[500] = {};

        minimum->Hesse();

        // minimum->Scan(2*howmuch + 2*binn + howmuchregen + howmuchomega + howmuchthree + howmuchsemi + howmuchelsee + 9, nstep, x, y, 0.0, 50.0);
        // minimum->Scan(2*howmuch + 2*binn + howmuchregen + howmuchomega + howmuchthree + howmuchsemi + howmuchelsee + 10, nstep, x, y, 0.0, 50.0);
        // minimum->Scan(2*howmuch + 2*binn + howmuchregen + howmuchomega + howmuchthree + howmuchsemi + howmuchelsee + 12, nstep, x, y, 0.0, 50.0);
        minimum->Scan(2 * howmuch + 2 * binn + howmuchregen + howmuchomega + howmuchthree + howmuchsemi + howmuchelsee + 15, nstep, x, y, 0.0, 50.0);
        minimum->Scan(2 * howmuch + 2 * binn + howmuchregen + howmuchomega + howmuchthree + howmuchsemi + howmuchelsee + 16, nstep, x, y, 0.0, 50.0);
        // minimum->Scan(2*howmuch + 2*binn + howmuchregen + howmuchomega + howmuchthree + howmuchsemi + howmuchelsee + 7, nstep, x, y, -0.2, 0.2);
        // minimum->Scan(2*howmuch + 2*binn + howmuchregen + howmuchomega + howmuchthree + howmuchsemi + howmuchelsee + 8, nstep, x, y, -0.2, 0.2);

        minimum->Minimize();
        minimum->Hesse();

        cout << minimum->X()[2 * howmuch + 2 * binn + howmuchregen + howmuchsemi + howmuchomega + howmuchthree + howmuchelsee + 7] << "+-" << minimum->Errors()[2 * howmuch + 2 * binn + howmuchregen + howmuchsemi + howmuchomega + howmuchthree + howmuchelsee + 7] << endl;
        cout << minimum->X()[2 * howmuch + 2 * binn + howmuchregen + howmuchsemi + howmuchomega + howmuchthree + howmuchelsee + 8] << "+-" << minimum->Errors()[2 * howmuch + 2 * binn + howmuchregen + howmuchsemi + howmuchomega + howmuchthree + howmuchelsee + 8] << endl;
        cout << minimum->X()[2 * howmuch + 2 * binn + howmuchregen + howmuchsemi + howmuchomega + howmuchthree + howmuchelsee + 9] << "+-" << minimum->Errors()[2 * howmuch + 2 * binn + howmuchregen + howmuchsemi + howmuchomega + howmuchthree + howmuchelsee + 9] << endl;
        cout << minimum->X()[2 * howmuch + 2 * binn + howmuchregen + howmuchsemi + howmuchomega + howmuchthree + howmuchelsee + 10] << "+-" << minimum->Errors()[2 * howmuch + 2 * binn + howmuchregen + howmuchsemi + howmuchomega + howmuchthree + howmuchelsee + 10] << endl;
        cout << minimum->X()[2 * howmuch + 2 * binn + howmuchregen + howmuchsemi + howmuchomega + howmuchthree + howmuchelsee + 11] << "+-" << minimum->Errors()[2 * howmuch + 2 * binn + howmuchregen + howmuchsemi + howmuchomega + howmuchthree + howmuchelsee + 11] << endl;
        cout << minimum->X()[2 * howmuch + 2 * binn + howmuchregen + howmuchsemi + howmuchomega + howmuchthree + howmuchelsee + 12] << "+-" << minimum->Errors()[2 * howmuch + 2 * binn + howmuchregen + howmuchsemi + howmuchomega + howmuchthree + howmuchelsee + 12] << endl;
        cout << minimum->X()[2 * howmuch + 2 * binn + howmuchregen + howmuchsemi + howmuchomega + howmuchthree + howmuchelsee + 13] << "+-" << minimum->Errors()[2 * howmuch + 2 * binn + howmuchregen + howmuchsemi + howmuchomega + howmuchthree + howmuchelsee + 13] << endl;
        cout << minimum->X()[2 * howmuch + 2 * binn + howmuchregen + howmuchsemi + howmuchomega + howmuchthree + howmuchelsee + 14] << "+-" << minimum->Errors()[2 * howmuch + 2 * binn + howmuchregen + howmuchsemi + howmuchomega + howmuchthree + howmuchelsee + 14] << endl;
        cout << minimum->X()[2 * howmuch + 2 * binn + howmuchregen + howmuchsemi + howmuchomega + howmuchthree + howmuchelsee + 15] << "+-" << minimum->Errors()[2 * howmuch + 2 * binn + howmuchregen + howmuchsemi + howmuchomega + howmuchthree + howmuchelsee + 15] << endl;
        cout << minimum->X()[2 * howmuch + 2 * binn + howmuchregen + howmuchsemi + howmuchomega + howmuchthree + howmuchelsee + 16] << "+-" << minimum->Errors()[2 * howmuch + 2 * binn + howmuchregen + howmuchsemi + howmuchomega + howmuchthree + howmuchelsee + 16] << endl;

        cout << "Reduced chi2: " << minimum->MinValue() / (binn - 4.) << endl;

        signal_before_fit->Reset();
    }

    return 0;
}
