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
  TauKs = PhysicsConstants::tau_S_nonCPT * 1E-9;                            // PDG fit not assuming CPT (s)
  TauKl = PhysicsConstants::tau_L * 1E-9;                                   // Kl mean life (s)
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
    Value = (1. + 2. * RePart) * exp(-GammaKs * abs(dt)) +
            (1. - 4. * RePart) * exp(-GammaKl * abs(dt)) -
            2. * exp(-0.5 * Gamma * abs(dt)) *
                ((1. - RePart) * cos(DMass * abs(dt)) -
                 3. * ImPart * sin(DMass * abs(dt)));
  }

  return (pow(Epsilon, 2) / (2. * Gamma)) * Value * 10000000;
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
  TauKs = PhysicsConstants::tau_S_nonCPT * 1E-9;                            // PDG fit not assuming CPT (s)
  TauKl = PhysicsConstants::tau_L * 1E-9;                                   // Kl mean life (s)
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

  return 0.5 * pow(Epsilon, 2) * Value;// * PhysicsConstants::br_ks_pippim * PhysicsConstants::br_ks_pi0pi0;
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
  TauKs = PhysicsConstants::tau_S_nonCPT * 1E-9;                            // PDG fit not assuming CPT (s)
  TauKl = PhysicsConstants::tau_L * 1E-9;                                   // Kl mean life (s)
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

  return 0.5 * pow(Epsilon, 2) * Value;// * PhysicsConstants::br_ks_pippim * PhysicsConstants::br_ks_pi0pi0;
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
  TauKs = PhysicsConstants::tau_S_nonCPT * 1E-9;                            // PDG fit not assuming CPT (s)
  TauKl = PhysicsConstants::tau_L * 1E-9;                                   // Kl mean life (s)
  MassDiff = PhysicsConstants::delta_mass_nonCPT;                           // M(Kl)-M(Ks) ( (h/2pi)s-1 ):
                                                                            // PDG fit not assuming CPT
  RePart = par[0];                                                          // PhysicsConstants::Re;
  ImPart = par[1];                                                          // PhysicsConstants::Im_nonCPT; // Im(epsilon'/epsilon) = Dphi/3;

  // All parameters are calculated taking into account that DT is in TauKs units
  GammaKs = 1.;
  GammaKl = TauKs / TauKl;
  Gamma = GammaKs + GammaKl;
  DMass = MassDiff * TauKs;

  Value = (1. + 2. * RePart) * (exp(-GammaKl * t1 - GammaKs * t2) +
                                exp(-GammaKs * t1 - GammaKl * t2) -
                                2. * exp(-0.5 * Gamma * (t1 + t2)) * cos(DMass * dt));

  return 0.5 * pow(Epsilon, 2) * Value;// * pow(PhysicsConstants::br_ks_pippim, 2);
}

////////////////////////////////////////////////////////////////////////////////////
// 
////////////////////////////////////////////////////////////////////////////////////

Double_t interf_function_00pm_to_fit_mock(const Double_t *x, const Double_t *par)
{
  Double_t Value = 0;
  Double_t Epsilon = 0, RePart = 0, Dphi = 0, TauKs = 0, TauKl = 0, MassDiff = 0;
  Double_t ImPart = 0, GammaKl = 0, GammaKs = 0, Gamma = 0, DMass = 0;

  Double_t t1 = x[0], t2 = x[1], dt;

  dt = t2 - t1;

  // Parameters from PDG2023

  Epsilon = PhysicsConstants::mod_epsilon;
  Dphi = PhysicsConstants::phi_pm_nonCPT - PhysicsConstants::phi_00_nonCPT; // phi(+-)-phi(00) (degrees)
  TauKs = PhysicsConstants::tau_S_nonCPT * 1E-9;                            // PDG fit not assuming CPT (s)
  TauKl = PhysicsConstants::tau_L * 1E-9;                                   // Kl mean life (s)
  MassDiff = PhysicsConstants::delta_mass_nonCPT;                           // M(Kl)-M(Ks) ( (h/2pi)s-1 ):
                                                                            // PDG fit not assuming CPT

  // All parameters are calculated taking into account that DT is in TauKs units
  GammaKs = 1.;
  GammaKl = TauKs / TauKl;
  Gamma = GammaKs + GammaKl;
  DMass = MassDiff * TauKs;

  Value = par[0] * exp(-GammaKl * t1 - GammaKs * t2) +
          par[1] * exp(-GammaKs * t1 - GammaKl * t2) -
          2. * exp(-0.5 * Gamma * (t1 + t2)) *
              (par[2] * cos(DMass * dt) +
               par[3] * sin(DMass * dt));

  return par[4] * 0.5 * pow(Epsilon, 2) * Value;// * PhysicsConstants::br_ks_pippim * PhysicsConstants::br_ks_pi0pi0;
}

Double_t interf_function_pm00_to_fit_mock(const Double_t *x, const Double_t *par)
{
  Double_t Value = 0;
  Double_t Epsilon = 0, RePart = 0, Dphi = 0, TauKs = 0, TauKl = 0, MassDiff = 0;
  Double_t ImPart = 0, GammaKl = 0, GammaKs = 0, Gamma = 0, DMass = 0;

  Double_t t1 = x[0], t2 = x[1], dt;

  dt = t2 - t1;

  // Parameters from PDG2023

  Epsilon = PhysicsConstants::mod_epsilon;
  Dphi = PhysicsConstants::phi_pm_nonCPT - PhysicsConstants::phi_00_nonCPT; // phi(+-)-phi(00) (degrees)
  TauKs = PhysicsConstants::tau_S_nonCPT * 1E-9;                            // PDG fit not assuming CPT (s)
  TauKl = PhysicsConstants::tau_L * 1E-9;                                   // Kl mean life (s)
  MassDiff = PhysicsConstants::delta_mass_nonCPT;                           // M(Kl)-M(Ks) ( (h/2pi)s-1 ):
                                                                            // PDG fit not assuming CPT

  // All parameters are calculated taking into account that DT is in TauKs units
  GammaKs = 1.;
  GammaKl = TauKs / TauKl;
  Gamma = GammaKs + GammaKl;
  DMass = MassDiff * TauKs;

  Value = par[0] * exp(-GammaKl * t1 - GammaKs * t2) +
          par[1] * exp(-GammaKs * t1 - GammaKl * t2) -
          2. * exp(-0.5 * Gamma * (t1 + t2)) *
              (par[2] * cos(DMass * dt) -
               par[3] * sin(DMass * dt));

  return par[4] * 0.5 * pow(Epsilon, 2) * Value;// * PhysicsConstants::br_ks_pippim * PhysicsConstants::br_ks_pi0pi0;
}

Double_t interf_function_pmpm_to_fit_mock(const Double_t *x, const Double_t *par)
{
  Double_t Value = 0;
  Double_t Epsilon = 0, RePart = 0, Dphi = 0, TauKs = 0, TauKl = 0, MassDiff = 0;
  Double_t ImPart = 0, GammaKl = 0, GammaKs = 0, Gamma = 0, DMass = 0;

  Double_t t1 = x[0], t2 = x[1], dt;

  dt = t2 - t1;

  // Parameters from PDG2023

  Epsilon = PhysicsConstants::mod_epsilon;
  Dphi = PhysicsConstants::phi_pm_nonCPT - PhysicsConstants::phi_00_nonCPT; // phi(+-)-phi(00) (degrees)
  TauKs = PhysicsConstants::tau_S_nonCPT * 1E-9;                            // PDG fit not assuming CPT (s)
  TauKl = PhysicsConstants::tau_L * 1E-9;                                   // Kl mean life (s)
  MassDiff = PhysicsConstants::delta_mass_nonCPT;                           // M(Kl)-M(Ks) ( (h/2pi)s-1 ):
                                                                            // PDG fit not assuming CPT

  // All parameters are calculated taking into account that DT is in TauKs units
  GammaKs = 1.;
  GammaKl = TauKs / TauKl;
  Gamma = GammaKs + GammaKl;
  DMass = MassDiff * TauKs;

  Value =(par[0] * exp(-GammaKl * t1 - GammaKs * t2) +
          par[1] * exp(-GammaKs * t1 - GammaKl * t2) -
          2. * par[2] * exp(-0.5 * Gamma * (t1 + t2)) * cos(DMass * dt));

  return par[3] * 0.5 * pow(Epsilon, 2) * Value;// * pow(PhysicsConstants::br_ks_pippim, 2);
}
