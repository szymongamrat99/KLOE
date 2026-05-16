#pragma once

#include <const.h>

using namespace PhysicsConstants;

namespace KLOE
{
  struct BRCorrectionFactors
  {

    Double_t br_phi_kskl_err = 0.4E-2;
    Double_t br_ks_pi0pi0_err = 0.05E-2;
    Double_t br_kl_pippim_err = 0.01E-3;
    Double_t br_ks_pippim_err = 0.05E-2;
    Double_t br_kl_pi0pi0_err = 0.06E-4;
    Double_t br_phi_omegapi0_err = 0.5E-5;
    Double_t br_omega_pippimpi0_err = 0.7E-2;
    Double_t br_ks_piele_err = 0.06E-4;
    Double_t br_kl_piele_err = 0.11E-2;
    Double_t br_ks_pimu_err = 0.2E-4;
    Double_t br_kl_pimu_err = 0.07E-2;
    Double_t br_kl_3pi0_err = 0.12E-2;

    Double_t br_phi_kskl_pm00 = br_phi_kskl * (br_ks_pi0pi0 * br_kl_pippim + br_ks_pippim * br_kl_pi0pi0);
    Double_t br_phi_kskl_pm00_err = sqrt(pow(br_phi_kskl_err * (br_ks_pi0pi0 * br_kl_pippim + br_ks_pippim * br_kl_pi0pi0), 2) +
                                           pow(br_ks_pi0pi0_err * br_phi_kskl * br_kl_pippim, 2) +
                                           pow(br_kl_pippim_err * br_phi_kskl * br_ks_pi0pi0, 2) +
                                           pow(br_ks_pippim_err * br_phi_kskl * br_kl_pi0pi0, 2) +
                                           pow(br_kl_pi0pi0_err * br_phi_kskl * br_ks_pippim, 2));

    Double_t br_phi_kskl_pm00_old = br_phi_kskl_old * (br_ks_pi0pi0_old * br_kl_pippim_old + br_ks_pippim_old * br_kl_pi0pi0_old);

    Double_t br_phi_omegapi0_pm00 = br_phi_omegapi0 * br_omega_pippimpi0;
    Double_t br_phi_omegapi0_pm00_err = sqrt(pow(br_phi_omegapi0_err * br_omega_pippimpi0, 2) +
                                              pow(br_omega_pippimpi0_err * br_phi_omegapi0, 2));

    Double_t br_phi_omegapi0_pm00_old = br_phi_omegapi0_old * br_omega_pippimpi0_old;

    Double_t br_phi_kskl_semileptonic = br_phi_kskl * (br_ks_piele * br_kl_pippim + br_ks_pippim * br_kl_piele + br_ks_pimu * br_kl_pippim + br_ks_pippim * br_kl_pimu);
    Double_t br_phi_kskl_semileptonic_err = sqrt(pow(br_phi_kskl_err * (br_ks_piele * br_kl_pippim + br_ks_pippim * br_kl_piele + br_ks_pimu * br_kl_pippim + br_ks_pippim * br_kl_pimu), 2) +
                                                  pow(br_ks_piele_err * br_phi_kskl * br_kl_pippim, 2) +
                                                  pow(br_kl_pippim_err * br_phi_kskl * (br_ks_piele + br_ks_pimu), 2) +
                                                  pow(br_ks_pippim_err * br_phi_kskl * (br_kl_piele + br_kl_pimu), 2) +
                                                  pow(br_kl_piele_err * br_phi_kskl * br_ks_pippim, 2) +
                                                  pow(br_ks_pimu_err * br_phi_kskl * br_kl_pippim, 2) +
                                                  pow(br_kl_pimu_err * br_phi_kskl * br_ks_pippim, 2));

    Double_t br_phi_kskl_semileptonic_old = br_phi_kskl_old * (br_ks_piele_old * br_kl_pippim_old + br_ks_pippim_old * br_kl_piele_old + br_ks_pimu_old * br_kl_pippim_old + br_ks_pippim_old * br_kl_pimu_old);

    Double_t br_phi_kskl_3pi0 = br_phi_kskl * (br_ks_pippim * br_kl_3pi0);
    Double_t br_phi_kskl_3pi0_err = sqrt(pow(br_phi_kskl_err * (br_ks_pippim * br_kl_3pi0), 2) +
                                         pow(br_ks_pippim_err * br_phi_kskl * br_kl_3pi0, 2) +
                                         pow(br_kl_3pi0_err * br_phi_kskl * br_ks_pippim, 2));

    Double_t br_phi_kskl_3pi0_old = br_phi_kskl_old * (br_ks_pippim_old * br_kl_3pi0_old);

    // Correction factors
    std::map<TString, Double_t> BRcorrectionFactors;
    std::map<TString, Double_t> BRcorrectionFactors_err;
    std::map<TString, Double_t> BranchingRatios;
    std::map<TString, Double_t> BranchingRatios_err;
    std::map<TString, Double_t> BranchingRatios_old;

    BRCorrectionFactors()
    {
      BRcorrectionFactors["Signal"] = br_phi_kskl_pm00 / br_phi_kskl_pm00_old;
      BRcorrectionFactors["Omega"] = br_phi_omegapi0_pm00 / br_phi_omegapi0_pm00_old;
      BRcorrectionFactors["3pi0"] = br_phi_kskl_3pi0 / br_phi_kskl_3pi0_old;
      BRcorrectionFactors["Semileptonic"] = (br_phi_kskl_semileptonic / br_phi_kskl_semileptonic_old);

      BRcorrectionFactors_err["Signal"] = abs(br_phi_kskl_pm00_err / br_phi_kskl_pm00_old);
      BRcorrectionFactors_err["Omega"] = abs(br_phi_omegapi0_pm00_err / br_phi_omegapi0_pm00_old);
      BRcorrectionFactors_err["3pi0"] = abs(br_phi_kskl_3pi0_err / br_phi_kskl_3pi0_old);
      BRcorrectionFactors_err["Semileptonic"] = abs(br_phi_kskl_semileptonic_err / br_phi_kskl_semileptonic_old);

      BRcorrectionFactors["Regeneration"] = 1.0; // No correction for regeneration, as it is data-driven
      BRcorrectionFactors["Other"] = 1.0;        // No correction for other background, as it is data-driven
      BRcorrectionFactors["Data"] = 1.0;         // No correction for data
      BRcorrectionFactors["MC sum"] = 1.0;       // No correction for MC sum, as it is a combination of other channels
      BRcorrectionFactors["pi+pi-pi+pi-"] = 1.0; // No correction for this channel, as it is not used in the fit

      BRcorrectionFactors_err["Regeneration"] = 0.0;
      BRcorrectionFactors_err["Other"] = 0.0;
      BRcorrectionFactors_err["Data"] = 0.0;
      BRcorrectionFactors_err["MC sum"] = 0.0;
      BRcorrectionFactors_err["pi+pi-pi+pi-"] = 0.0;

      // ---
      BranchingRatios["Signal"] = br_phi_kskl_pm00;
      BranchingRatios["Omega"] = br_phi_omegapi0_pm00;
      BranchingRatios["3pi0"] = br_phi_kskl_3pi0;
      BranchingRatios["Semileptonic"] = br_phi_kskl_semileptonic;

      BranchingRatios_err["Signal"] = br_phi_kskl_pm00_err;
      BranchingRatios_err["Omega"] = br_phi_omegapi0_pm00_err;
      BranchingRatios_err["3pi0"] = br_phi_kskl_3pi0_err;
      BranchingRatios_err["Semileptonic"] = br_phi_kskl_semileptonic_err;

      BranchingRatios_old["Signal"] = br_phi_kskl_pm00_old;
      BranchingRatios_old["Omega"] = br_phi_omegapi0_pm00_old;
      BranchingRatios_old["3pi0"] = br_phi_kskl_3pi0_old;
      BranchingRatios_old["Semileptonic"] = br_phi_kskl_semileptonic_old;
    }
  };

};