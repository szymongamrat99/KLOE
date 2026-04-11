#include <interf_function.h>
#include <const.h>

#include <TF1.h>
#include <TF2.h>

Double_t interf_function(const Double_t *x, const Double_t *par)
{
  Double_t Value = 0;
  Double_t Epsilon = 0, RePart = 0, Dphi = 0, TauKs = 0, TauKl = 0, MassDiff = 0;
  Double_t ImPart = 0, GammaKl = 0, GammaKs = 0, Gamma = 0, DMass = 0;

  Double_t dt = x[0];

  // Parameters from PDG2023

  Epsilon = PhysicsConstants::mod_epsilon;
  Dphi = PhysicsConstants::phi_pm_nonCPT - PhysicsConstants::phi_00_nonCPT; // phi(+-)-phi(00) (degrees)
  TauKs = PhysicsConstants::tau_S_nonCPT;                            // PDG fit not assuming CPT (s)
  TauKl = PhysicsConstants::tau_L;                                   // Kl mean life (s)
  MassDiff = PhysicsConstants::delta_mass_nonCPT;                           // M(Kl)-M(Ks) ( (h/2pi)s-1 ):

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
    Value = (1. + 2. * RePart) * exp(-GammaKs * std::abs(dt)) +
            (1. - 4. * RePart) * exp(-GammaKl * std::abs(dt)) -
            2. * exp(-0.5 * Gamma * std::abs(dt)) *
                ((1. - RePart) * cos(DMass * std::abs(dt)) -
                 3. * ImPart * sin(DMass * std::abs(dt)));
  }

  return (std::pow(Epsilon, 2) / (2. * Gamma)) * Value * 1000000;
}

Double_t interf_function_00pm(const Double_t *x, const Double_t *par)
{
  Double_t Value = 0;
  Double_t Epsilon = 0, RePart = 0, Dphi = 0, TauKs = 0, TauKl = 0, MassDiff = 0;
  Double_t ImPart = 0, GammaKl = 0, GammaKs = 0, Gamma = 0, DMass = 0;

  Double_t t1 = x[0], t2 = x[1], dt;

  dt = t2 - t1;

  // Parameters from PDG2023

  Epsilon = PhysicsConstants::mod_epsilon;
  Dphi = PhysicsConstants::phi_pm_nonCPT - PhysicsConstants::phi_00_nonCPT; // phi(+-)-phi(00) (degrees)
  TauKs = PhysicsConstants::tau_S_nonCPT;                            // PDG fit not assuming CPT (s)
  TauKl = PhysicsConstants::tau_L;                                   // Kl mean life (s)
  MassDiff = PhysicsConstants::delta_mass_nonCPT;                           // M(Kl)-M(Ks) ( (h/2pi)s-1 ):
                                                                            // PDG fit not assuming CPT
  RePart = par[0];                                                          // PhysicsConstants::Re;
  ImPart = par[1];                                                          // PhysicsConstants::Im_nonCPT; // Im(epsilon'/epsilon) = Dphi/3;

  // All parameters are calculated taking into account that DT is in TauKs units
  GammaKs = 1.;
  GammaKl = TauKs / TauKl;
  Gamma = GammaKs + GammaKl;
  DMass = MassDiff * TauKs;

  Value = (1. - 4. * RePart) * exp(-GammaKl * t1 - GammaKs * t2) +
          (1. + 2. * RePart) * exp(-GammaKs * t1 - GammaKl * t2) -
          2. * exp(-0.5 * Gamma * (t1 + t2)) *
              ((1. - RePart) * cos(DMass * dt) +
               3. * ImPart * sin(DMass * dt));

  return 0.5 * std::pow(Epsilon, 2) * Value;// * PhysicsConstants::br_ks_pippim * PhysicsConstants::br_ks_pi0pi0;
}

Double_t interf_function_pm00(const Double_t *x, const Double_t *par)
{
  Double_t Value = 0;
  Double_t Epsilon = 0, RePart = 0, Dphi = 0, TauKs = 0, TauKl = 0, MassDiff = 0;
  Double_t ImPart = 0, GammaKl = 0, GammaKs = 0, Gamma = 0, DMass = 0;

  Double_t t1 = x[0], t2 = x[1], dt;

  dt = t2 - t1;

  // Parameters from PDG2023

  Epsilon = PhysicsConstants::mod_epsilon;
  Dphi = PhysicsConstants::phi_pm_nonCPT - PhysicsConstants::phi_00_nonCPT; // phi(+-)-phi(00) (degrees)
  TauKs = PhysicsConstants::tau_S_nonCPT;                            // PDG fit not assuming CPT (s)
  TauKl = PhysicsConstants::tau_L;                                   // Kl mean life (s)
  MassDiff = PhysicsConstants::delta_mass_nonCPT;                           // M(Kl)-M(Ks) ( (h/2pi)s-1 ):
                                                                            // PDG fit not assuming CPT
  RePart = par[0];                                                          // PhysicsConstants::Re;
  ImPart = par[1];                                                          // PhysicsConstants::Im_nonCPT; // Im(epsilon'/epsilon) = Dphi/3;

  // All parameters are calculated taking into account that DT is in TauKs units
  GammaKs = 1.;
  GammaKl = TauKs / TauKl;
  Gamma = GammaKs + GammaKl;
  DMass = MassDiff * TauKs;

  Value = (1. + 2. * RePart) * exp(-GammaKl * t1 - GammaKs * t2) +
          (1. - 4. * RePart) * exp(-GammaKs * t1 - GammaKl * t2) -
          2. * exp(-0.5 * Gamma * (t1 + t2)) *
              ((1. - RePart) * cos(DMass * dt) -
               3. * ImPart * sin(DMass * dt));

  return 0.5 * std::pow(Epsilon, 2) * Value;// * PhysicsConstants::br_ks_pippim * PhysicsConstants::br_ks_pi0pi0;
}

Double_t interf_function_pmpm(const Double_t *x, const Double_t *par)
{
  Double_t Value = 0;
  Double_t Epsilon = 0, RePart = 0, Dphi = 0, TauKs = 0, TauKl = 0, MassDiff = 0;
  Double_t ImPart = 0, GammaKl = 0, GammaKs = 0, Gamma = 0, DMass = 0;

  Double_t t1 = x[0], t2 = x[1], dt;

  dt = t2 - t1;

  // Parameters from PDG2023

  Epsilon = PhysicsConstants::mod_epsilon;
  Dphi = PhysicsConstants::phi_pm_nonCPT - PhysicsConstants::phi_00_nonCPT; // phi(+-)-phi(00) (degrees)
  TauKs = PhysicsConstants::tau_S_nonCPT;                            // PDG fit not assuming CPT (s)
  TauKl = PhysicsConstants::tau_L;                                   // Kl mean life (s)
  MassDiff = PhysicsConstants::delta_mass_nonCPT;                           // M(Kl)-M(Ks) ( (h/2pi)s-1 ):
                                                                            // PDG fit not assuming CPT
  RePart = par[0];                                                          // PhysicsConstants::Re;                                                       // PhysicsConstants::Im_nonCPT; // Im(epsilon'/epsilon) = Dphi/3;

  // All parameters are calculated taking into account that DT is in TauKs units
  GammaKs = 1.;
  GammaKl = TauKs / TauKl;
  Gamma = GammaKs + GammaKl;
  DMass = MassDiff * TauKs;

  Value = (1. + 2. * RePart) * (exp(-GammaKl * t1 - GammaKs * t2) +
                                exp(-GammaKs * t1 - GammaKl * t2) -
                                2. * exp(-0.5 * Gamma * (t1 + t2)) * cos(DMass * dt));

  return 0.5 * std::pow(Epsilon, 2) * Value;// * std::pow(PhysicsConstants::br_ks_pippim, 2);
}

Double_t double_exponential(const Double_t *x, const Double_t *par)
{
  Double_t gammaS = 1.0; // Wartość gamma_S do ustawienia zakresów
  Double_t gammaL = PhysicsConstants::tau_S_nonCPT / PhysicsConstants::tau_L;

  Double_t dt = x[0];

  Double_t value = 0.;

  if (dt >= 0)
  {
    value = (1. + 2. * PhysicsConstants::Re) * exp(-gammaL * dt) +
            (1. - 4. * PhysicsConstants::Re) * exp(-gammaS * dt);
  }
  else
  {
    value = (1. + 2. * PhysicsConstants::Re) * exp(-gammaS * std::abs(dt)) +
            (1. - 4. * PhysicsConstants::Re) * exp(-gammaL * std::abs(dt));
  }

  return value;
}