#include "kloe_class.h"
#include "../../Include/const.h"
#include "interference.h"
#include <TMath.h>
#include <TLorentzVector.h>

namespace KLOE
{
    Double_t interference::interf_function(const Float_t x, Int_t check = 1, const Double_t *par = 0)
    {
            Double_t Value = 0;
            Double_t Epsilon = 0, RePart = 0, Dphi = 0, TauKs = 0, TauKl = 0, MassDiff = 0;
            Double_t ImPart = 0, GammaKl = 0, GammaKs = 0, Gamma = 0, DMass = 0;

            Float_t dt = x;

            //Parameters from PDG2023

            Epsilon   = mod_epsilon;
            Dphi      = phi_pm_nonCPT - phi_00_nonCPT;            // phi(+-)-phi(00) (degrees)
            TauKs     = tau_S_nonCPT; // PDG fit not assuming CPT (s)
            TauKl     = tau_L;        // Kl mean life (s)
            MassDiff  = delta_mass_nonCPT;     // M(Kl)-M(Ks) ( (h/2pi)s-1 ):
                                                    // PDG fit not assuming CPT

            if(check == 0){
                Int_t howmuchregen = par[0];
                Int_t howmuchomega = par[1];
                Int_t howmuchthree = par[2];
                Int_t howmuchsemi = par[3];
                Int_t howmuchelsee = par[4];
                Int_t binn = par[5];
                Int_t howmuch = par[6];

                RePart    = par[2*howmuch + 2*binn + howmuchregen + howmuchomega + howmuchthree + howmuchsemi + howmuchelsee + 7];   
                ImPart  = par[2*howmuch + 8 + 2*binn + howmuchregen + howmuchomega + howmuchthree + howmuchsemi + howmuchelsee];
            }
            else{
                RePart = Re;
                ImPart = M_PI*( Dphi/3. )/180.; //Im(epsilon'/epsilon) = Dphi/3;
            }

            //All parameters are calculated taking into account that DT is in TauKs units
            GammaKs = 1.;
            GammaKl = TauKs / TauKl;
            Gamma   = GammaKs + GammaKl;
            DMass   = MassDiff * TauKs;


            if( dt >= 0. )
            {
                Value = ( 1. + 2.*RePart ) * exp(-GammaKl*dt) +
                    ( 1. - 4.*RePart ) * exp(-GammaKs*dt) -
                    2. * exp(-0.5*Gamma*dt) * 
                    ( ( 1.-RePart ) * cos(DMass*dt) + 
                    3. * ImPart * sin(DMass*dt) );
            }
            else
            {
                Value = ( 1. + 2.*RePart ) * exp(-GammaKs*abs(dt)) +
                    ( 1. - 4.*RePart ) * exp(-GammaKl*abs(dt)) -
                    2. * exp(-0.5*Gamma*abs(dt)) * 
                    ( ( 1. - RePart ) * cos(DMass*abs(dt)) - 
                    3. * ImPart * sin(DMass*abs(dt)) );
            }
            
            

            return (pow(Epsilon,2)/(2.*Gamma)) * Value * 100000 ;
    }

    void interference::bin_extraction(UInt_t channel, TH1* histogram)
    {
        for (Int_t i = 0; i < bin_number; i++)
            b[channel][i] = (histogram->GetBinContent(i + 1));

        for (Int_t i = 0; i < bin_number; i++)
            e[channel][i] = (histogram->GetBinError(i + 1));
    };

    Double_t interference::interf_chi2_split(const Double_t *xx)
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

        TH1 *h1 = new TH1F("Sig", "", binn, x_min, x_max);

        for (Int_t i = 7; i < howmuch + 7; i++)
            h1->Fill(xx[i], interf_function(xx[i + howmuch], 0, xx));

        h1->Scale(71532 / h1->Integral(1, binn));

        TH1 *h1regen = new TH1F("Regen", "", binn, x_min, x_max);

        for (Int_t i = 0; i < howmuchregen; i++)
            h1regen->Fill(xx[i + 2 * howmuch + 2 * binn + 7], N_regen);
        interference::bin_extraction(1, frac[1]);

        TH1 *h1omega = new TH1F("Omega", "", binn, x_min, x_max);

        for (Int_t i = 0; i < howmuchomega; i++)
            h1omega->Fill(xx[i + 2 * howmuch + 2 * binn + 7 + howmuchregen], N_omega);
        for (Int_t i = 0; i < binn; i++)
            bomega[i] = (h1omega->GetBinContent(i + 1));
        for (Int_t i = 0; i < binn; i++)
            eomega[i] = (h1omega->GetBinError(i + 1));

        TH1 *h1three = new TH1F("Thre", "", binn, x_min, x_max);

        for (Int_t i = 0; i < howmuchthree; i++)
            h1three->Fill(xx[i + 2 * howmuch + 2 * binn + 7 + howmuchregen + howmuchomega], N_three);
        for (Int_t i = 0; i < binn; i++)
            bthree[i] = (h1three->GetBinContent(i + 1));
        for (Int_t i = 0; i < binn; i++)
            ethree[i] = (h1three->GetBinError(i + 1));
        TH1 *h1semi = new TH1F("Semi", "", binn, x_min, x_max);

        for (Int_t i = 0; i < howmuchsemi; i++)
            h1semi->Fill(xx[i + 2 * howmuch + 2 * binn + 7 + howmuchregen + howmuchomega + howmuchthree], N_semi);
        for (Int_t i = 0; i < binn; i++)
            bsemi[i] = (h1semi->GetBinContent(i + 1));
        for (Int_t i = 0; i < binn; i++)
            esemi[i] = (h1semi->GetBinError(i + 1));
        TH1 *h1elsee = new TH1F("elsee", "", binn, x_min, x_max);

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

    Double_t interference::interf_chi2_window(const Double_t *xx)
    {

    };

    Double_t interference::interf_chi2_excluded(const Double_t *xx)
    {

    };

    Double_t interference::interf_chi2_all(const Double_t *xx)
    {

    };

    Double_t interference::interf_chi2_mc(const Double_t *xx)
    {

    };

    Double_t interference::interf_chi2_bcg(const Double_t *xx)
    {

    };

    Double_t interference::interf_chi2(const Double_t *xx)
    {
        if(mode == "split") return interf_chi2_split(xx);
        else if(mode == "window") return interf_chi2_window(xx);
        else if(mode == "excluded") return interf_chi2_excluded(xx);
        else if(mode == "mc") return interf_chi2_mc(xx);
        else if(mode == "bcg") return interf_chi2_bcg(xx);
        else return -1.0;
    };
}