#pragma once

#include <const.h>

using namespace PhysicsConstants;

namespace KLOE
{
  struct BRCorrectionFactors
  {
    Double_t br_phi_kskl_pm00 = br_phi_kskl * (br_ks_pi0pi0 * br_kl_pippim + br_ks_pippim * br_kl_pi0pi0);
    Double_t br_phi_kskl_pm00_old = br_phi_kskl_old * (br_ks_pi0pi0_old * br_kl_pippim_old + br_ks_pippim_old * br_kl_pi0pi0_old);

    Double_t br_phi_omegapi0_pm00 = br_phi_omegapi0 * br_omega_pippimpi0;
    Double_t br_phi_omegapi0_pm00_old = br_phi_omegapi0_old * br_omega_pippimpi0_old;

    Double_t br_phi_kskl_semileptonic = br_phi_kskl * (br_ks_piele * br_kl_pippim + br_ks_pippim * br_kl_piele + br_ks_pimu * br_kl_pippim + br_ks_pippim * br_kl_pimu);
    Double_t br_phi_kskl_semileptonic_old = br_phi_kskl_old * (br_ks_piele_old * br_kl_pippim_old + br_ks_pippim_old * br_kl_piele_old + br_ks_pimu_old * br_kl_pippim_old + br_ks_pippim_old * br_kl_pimu_old);

    Double_t br_phi_kskl_3pi0 = br_phi_kskl * (br_ks_pippim * br_kl_3pi0);
    Double_t br_phi_kskl_3pi0_old = br_phi_kskl_old * (br_ks_pippim_old * br_kl_3pi0_old);

    // Correction factors
    std::map<TString, Double_t> BRcorrectionFactors;

    BRCorrectionFactors()
    {
      BRcorrectionFactors["Signal"] = br_phi_kskl_pm00 / br_phi_kskl_pm00_old;
      BRcorrectionFactors["Omega"] = br_phi_omegapi0_pm00 / br_phi_omegapi0_pm00_old;
      BRcorrectionFactors["3pi0"] = br_phi_kskl_3pi0 / br_phi_kskl_3pi0_old;
      BRcorrectionFactors["Semileptonic"] = (br_phi_kskl_semileptonic / br_phi_kskl_semileptonic_old);

      BRcorrectionFactors["Regeneration"] = 1.0; // No correction for regeneration, as it is data-driven
      BRcorrectionFactors["Other"] = 1.0;        // No correction for other background, as it is data-driven
      BRcorrectionFactors["Data"] = 1.0;         // No correction for data
      BRcorrectionFactors["MC sum"] = 1.0;       // No correction for MC sum, as it is a combination of other channels
      BRcorrectionFactors["pi+pi-pi+pi-"] = 1.0; // No correction for this channel, as it is not used in the fit
    }
  };

};