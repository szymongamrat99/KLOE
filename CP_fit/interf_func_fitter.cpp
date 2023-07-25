#include "TMath.h"
#include "TString.h"
#include "TH1.h"
#include <string>
#include <iostream>    

class interf_func_fitter
{
    public:
        UInt_t bin_number = 201;
        Double_t x_min = -100., x_max = 100., y_min = 0., y_max = 1000.;

        void bin_extraction(UInt_t channel, TH1* histogram);

        Double_t interf_chi2_split(const Double_t *xx);
        Double_t interf_chi2_window(const Double_t *xx);
        Double_t interf_chi2_excluded(const Double_t *xx);
        Double_t interf_chi2_all(const Double_t *xx);

        Double_t interf_chi2_mc(const Double_t *xx);
        Double_t interf_chi2_bcg(const Double_t *xx);

    private:
        UInt_t bin_number = 201;
        TString mode = "all";  // "split", "window", "excluded, "mc", "bcg", "all"

        TH1 *frac[7];
        Double_t *b[7], *e[7]; 

        void histos_init()
        {
            for(Int_t i = 0; i < 7; i++)
                frac[i] = new TH1F(("histo" + to_string(i)).c_str(), "", bin_number, x_min, x_max);
        };

};

void interf_func_fitter::bin_extraction(UInt_t channel, TH1* histogram)
{
    for (Int_t i = 0; i < bin_number; i++)
        b[channel][i] = (histogram->GetBinContent(i + 1));

    for (Int_t i = 0; i < bin_number; i++)
        e[channel][i] = (histogram->GetBinError(i + 1));
};

Double_t interf_func_fitter::interf_chi2_split(const Double_t *xx)
{
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

    TH1 *h1 = new TH1F("Sig", "", binn, -90.0, 90.0);

    for (Int_t i = 7; i < howmuch + 7; i++)
        h1->Fill(xx[i], fitf(xx[i + howmuch], xx));

    h1->Scale(71532 / h1->Integral(1, binn));

    TH1 *h1regen = new TH1F("Regen", "", binn, -90.0, 90.0);

    for (Int_t i = 0; i < howmuchregen; i++)
        h1regen->Fill(xx[i + 2 * howmuch + 2 * binn + 7], N_regen);
    interf_func_fitter::bin_extraction(1, frac[1]);

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
};

Double_t interf_func_fitter::interf_chi2_window(const Double_t *xx)
{

};

Double_t interf_func_fitter::interf_chi2_excluded(const Double_t *xx)
{

};

Double_t interf_func_fitter::interf_chi2_all(const Double_t *xx)
{

};

Double_t interf_func_fitter::interf_chi2_mc(const Double_t *xx)
{

};

Double_t interf_func_fitter::interf_chi2_bcg(const Double_t *xx)
{

};