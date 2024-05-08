#include <TCanvas.h>
#include <TF1.h>
#include <TStyle.h>
#include <TAxis.h>
#include <TMath.h>

#include "../../Include/const.h"

Double_t interf_function(Double_t *x, Double_t *par)
{
    Double_t Value = 0;
    Double_t Epsilon = 0, RePart = 0, Dphi = 0, TauKs = 0, TauKl = 0, MassDiff = 0;
    Double_t ImPart = 0, GammaKl = 0, GammaKs = 0, Gamma = 0, DMass = 0;

    Double_t dt = x[0];

    // Parameters from PDG2023

    Epsilon = mod_epsilon;
    Dphi = phi_pm_nonCPT - phi_00_nonCPT; // phi(+-)-phi(00) (degrees)
    TauKs = tau_S_nonCPT * pow(10, -9);   // PDG fit not assuming CPT (s)
    TauKl = tau_L * pow(10, -9);          // Kl mean life (s)
    MassDiff = delta_mass_nonCPT;         // M(Kl)-M(Ks) ( (h/2pi)s-1 ):
                                          // PDG fit not assuming CPT
    RePart = par[0];
    ImPart = par[1]; // Im(epsilon'/epsilon) = Dphi/3;

    // All parameters are calculated taking into account that DT is in TauKs units
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
                2. * exp(-0.5 * Gamma * abs(dt)) *
                    ((1. - RePart) * cos(DMass * abs(dt)) -
                     3. * ImPart * sin(DMass * abs(dt)));
    }

    return (pow(Epsilon, 2) / (2. * Gamma)) * Value * 100000;
}

void interf_function_draw()
{
    TCanvas *canva = new TCanvas();
    TF1 *func_with_im = new TF1("Interference function2", &interf_function, -20.0, 20.0, 2);
    func_with_im->SetParameter(0, Re);
    func_with_im->SetParameter(1, M_PI * ((phi_pm_nonCPT - phi_00_nonCPT) / 3.) / 180.);

    TF1 *func_wo_im = new TF1("Interference function1", &interf_function, -20.0, 20.0, 2);
    func_wo_im->SetParameter(0, 0);
    func_wo_im->SetParameter(1, 0);

    canva->SetGrid(1, 1);

    gStyle->SetGridStyle(3);
    gStyle->SetGridWidth(2);

    func_with_im->SetNpx(1E6);
    func_with_im->SetLineColor(kRed);
    func_with_im->SetLineWidth(3);

    func_wo_im->SetNpx(1E6);
    func_wo_im->SetLineColor(kBlack);
    func_wo_im->SetLineWidth(3);

    func_with_im->SetTitle("");

    func_with_im->GetYaxis()->SetRangeUser(0, 0.4);
    func_with_im->GetYaxis()->SetTitle("I [-]");
    func_with_im->GetYaxis()->CenterTitle(1);

    func_with_im->GetXaxis()->SetTitle("#Deltat [#tau_{S}]");
    func_with_im->GetXaxis()->CenterTitle(1);

    func_with_im->Draw();
    func_wo_im->Draw("SAME");

    TLegend *legend = new TLegend(0.6, 0.73, 0.9, 0.9,"");

    legend->AddEntry(func_wo_im, "#varepsilon'/#varepsilon = 0", "l");
    legend->AddEntry(func_with_im, "#splitline{Re(#varepsilon'/#varepsilon) = 1.66#times10^{-3}}{Im(#varepsilon'/#varepsilon) = -1.73#times10^{-3}}", "l");

    legend->SetTextSize(0.025);
    legend->SetTextFont(42);

    legend->Draw();

    canva->Print("interf_func.png");
}