#include "../../Include/const.h"

Double_t interf_function(const Float_t x)
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

        RePart = Re;
        ImPart = M_PI*( Dphi/3. )/180.; //Im(epsilon'/epsilon) = Dphi/3;

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

void interf_function_draw()
{
    TCanvas *canva = new TCanvas();
    TF1 *func = new TF1("Interference function","interf_function(x)",-20,20);

    canva->SetGrid(1,1);

    gStyle->SetGridStyle(3);
    gStyle->SetGridWidth(2);

    func->SetNpx(1E6);
    func->SetLineColor(kBlack);
    func->SetLineWidth(3);

    func->SetTitle("");

    func->GetYaxis()->SetRangeUser(0,0.3); 
    func->GetYaxis()->SetTitle("I [-]");
    func->GetYaxis()->CenterTitle(1);

    func->GetXaxis()->SetTitle("#Deltat [#tau_{S}]");
    func->GetXaxis()->CenterTitle(1);

    func->Draw();
}