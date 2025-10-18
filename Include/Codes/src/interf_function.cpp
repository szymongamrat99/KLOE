#include <interf_function.h>
#include <const.h>

Double_t interf_function(const Float_t x)
{
    Double_t Value = 0;
    Double_t Epsilon = 0, RePart = 0, Dphi = 0, TauKs = 0, TauKl = 0, MassDiff = 0;
    Double_t ImPart = 0, GammaKl = 0, GammaKs = 0, Gamma = 0, DMass = 0;

    Double_t dt = x;

    // Parameters from PDG2023

    Epsilon = PhysicsConstants::mod_epsilon;
    Dphi = PhysicsConstants::phi_pm_nonCPT - PhysicsConstants::phi_00_nonCPT; // phi(+-)-phi(00) (degrees)
    TauKs = PhysicsConstants::tau_S_nonCPT * pow(10, -9);   // PDG fit not assuming CPT (s)
    TauKl = PhysicsConstants::tau_L * pow(10, -9);          // Kl mean life (s)
    MassDiff = PhysicsConstants::delta_mass_nonCPT;         // M(Kl)-M(Ks) ( (h/2pi)s-1 ):
                                          // PDG fit not assuming CPT
    RePart = PhysicsConstants::Re;
    ImPart = PhysicsConstants::Im_nonCPT; // Im(epsilon'/epsilon) = Dphi/3;

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
