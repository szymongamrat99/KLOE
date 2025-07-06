#include <interf_function.h>
#include <PhysicsConstants.h>
#include <SystemPaths.h>

Double_t interf_function(const Float_t x)
{
    PhysicsConstants constants(SystemPath::pdgConstFilePath);

    Double_t Value = 0;
    Double_t Epsilon = 0, RePart = 0, Dphi = 0, TauKs = 0, TauKl = 0, MassDiff = 0;
    Double_t ImPart = 0, GammaKl = 0, GammaKs = 0, Gamma = 0, DMass = 0;

    Float_t dt = x;

    // Parameters from PDG2023

    Epsilon = constants.getModulusEpsilon();
    Dphi = constants.getPhiPMNonCPT() - constants.getPhi00NonCPT(); // phi(+-)-phi(00) (degrees)
    TauKs = constants.getTauSNonCPT() * pow(10, -9);   // PDG fit not assuming CPT (s)
    TauKl = constants.getTauL() * pow(10, -9);          // Kl mean life (s)
    MassDiff = constants.getDeltaMassNonCPT();         // M(Kl)-M(Ks) ( (h/2pi)s-1 ):
                                          // PDG fit not assuming CPT
    RePart = constants.getRe(); // Re(epsilon'/epsilon) = RePart
    ImPart = M_PI * (Dphi / 3.) / 180.; // Im(epsilon'/epsilon) = Dphi/3;

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