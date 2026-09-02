#include <GeneratedVariables.h>
#include <const.h>

#include <algorithm>
#include <ChannelProperties.h>

bool GeneratedVariables::IsPrimaryFromPhi(Int_t kineId) const
{
  auto it = nodes.find(kineId);
  while (it != nodes.end())
  {
    if (it->second.motherID == 50)
      return true; // dotarliśmy do phi
    if (it->second.motherID == 10 && it->second.particleID == 16)
      return false; // regeneracja Kl->Ks
    it = nodes.find(it->second.parentKine);
  }
  return false;
}

bool GeneratedVariables::IsFromKs(Int_t kineId) const
{
  auto it = nodes.find(kineId);
  while (it != nodes.end())
  {
    if (it->second.particleID == 16)
      return true; // ten węzeł to Ks
    it = nodes.find(it->second.parentKine);
  }
  return false;
}

bool GeneratedVariables::IsFromKl(Int_t kineId) const
{
  auto it = nodes.find(kineId);
  while (it != nodes.end())
  {
    if (it->second.particleID == 10)
      return true; // ten węzeł to Kl
    it = nodes.find(it->second.parentKine);
  }
  return false;
}

void GeneratedVariables::classifyChannel(
    ErrorHandling::ErrorLogs &logger,
    Int_t ntmc,
    Int_t nvtxmc,
    const Int_t *pidmcOld,
    const Int_t *vtxmcOld,
    const Int_t *motherOld,
    const Int_t *kine,
    const Int_t *kinmom,
    UInt_t mcflag,
    Int_t &mctruth_int,
    Int_t &semileptonic_flag,
    Int_t &other_flag,
    Int_t &isr_flag,
    Int_t &has_pm,
    Int_t &has_00,
    Int_t &has_000,
    Int_t &has_semileptonic)
{
  UInt_t Ks = 0, Kl = 0, Ksregen = 0, piplusks = 0, pipluskl = 0, piminusks = 0, piminuskl = 0, muonplusks = 0, muonpluskl = 0, muonminusks = 0, muonminuskl = 0, electronks = 0, electronkl = 0, positronks = 0, positronkl = 0, pi0ks = 0, pi0kl = 0, pi0phi = 0, piplusphi = 0, piminusphi = 0, otherphi = 0, otherkl = 0, otherks = 0, gammaphi = 0, neutrinoks = 0, neutrinokl = 0, neutrinokplus = 0, neutrinokminus = 0, Kllambda = 0, Kslambda = 0, Ksgamma = 0, Klgamma = 0, Kplusgamma = 0, Kminusgamma = 0, pi0kplus = 0, pi0kminus = 0, pipluskplus = 0, pipluskminus = 0, piminuskplus = 0, piminuskminus = 0, muonpluskplus = 0, muonpluskminus = 0, muonminuskplus = 0, muonminuskminus = 0, electronkplus = 0, electronkminus = 0, positronkplus = 0, positronkminus = 0, otherkplus = 0, otherkminus = 0, Kplus = 0, Kminus = 0;

  isr_flag = 0;
  semileptonic_flag = 0;
  other_flag = 0;
  has_pm = 0;
  has_00 = 0;
  has_000 = 0;
  has_semileptonic = 0;

  nodes.clear();

  if (mcflag == 1)
  {
    for (Int_t j = 0; j < ntmc; ++j)
    {
      MCNode node;
      node.particleID = pidmcOld[j];
      node.motherID = motherOld[vtxmcOld[j] - 1];
      node.vertexID = vtxmcOld[j];
      node.parentKine = kinmom[vtxmcOld[j] - 1];
      nodes[kine[j]] = node;

      if (motherOld[vtxmcOld[j] - 1] == 50)
      {
        switch (pidmcOld[j])
        {
        case 10:
          Kl++;
          break;
        case 16:
          Ks++;
          break;
        case 11:
          Kplus++;
          break;
        case 12:
          Kminus++;
          break;
        case 7:
          pi0phi++;
          break;
        case 8:
          piplusphi++;
          break;
        case 9:
          piminusphi++;
          break;
        case 1:
          gammaphi++;
          break;
        default:
          otherphi++;
          break;
        }
      }
      else if (motherOld[vtxmcOld[j] - 1] == 10)
      {
        switch (pidmcOld[j])
        {
        case 16:
          Ksregen++;
          break;
        case 7:
          pi0kl++;
          break;
        case 8:
          pipluskl++;
          break;
        case 9:
          piminuskl++;
          break;
        case 5:
          muonpluskl++;
          break;
        case 6:
          muonminuskl++;
          break;
        case 2:
          positronkl++;
          break;
        case 3:
          electronkl++;
          break;
        case 4:
          neutrinokl++;
          break;
        case 1:
          Klgamma++;
          break;
        case 18:
          Kllambda++;
          break;
        default:
          otherkl++;
          break;
        }
      }
      else if (motherOld[vtxmcOld[j] - 1] == 16)
      {
        switch (pidmcOld[j])
        {
        case 7:
          pi0ks++;
          break;
        case 8:
          piplusks++;
          break;
        case 9:
          piminusks++;
          break;
        case 5:
          muonplusks++;
          break;
        case 6:
          muonminusks++;
          break;
        case 2:
          positronks++;
          break;
        case 3:
          electronks++;
          break;
        case 4:
          neutrinoks++;
          break;
        case 1:
          Ksgamma++;
          break;
        case 18:
          Kslambda++;
          break;
        default:
          otherks++;
          break;
        }
      }
      else if (motherOld[vtxmcOld[j] - 1] == 11)
      {
        switch (pidmcOld[j])
        {
        case 7:
          pi0kplus++;
          break;
        case 8:
          pipluskplus++;
          break;
        case 9:
          piminuskplus++;
          break;
        case 5:
          muonpluskplus++;
          break;
        case 6:
          muonminuskplus++;
          break;
        case 2:
          positronkplus++;
          break;
        case 3:
          electronkplus++;
          break;
        case 1:
          Kplusgamma++;
          break;
        case 4:
          neutrinokplus++;
          break;
        default:
        {
          otherkplus++;
          break;
        }
        }
      }
      else if (motherOld[vtxmcOld[j] - 1] == 12)
      {
        switch (pidmcOld[j])
        {
        case 7:
          pi0kminus++;
          break;
        case 8:
          pipluskminus++;
          break;
        case 9:
          piminuskminus++;
          break;
        case 5:
          muonpluskminus++;
          break;
        case 6:
          muonminuskminus++;
          break;
        case 2:
          positronkminus++;
          break;
        case 3:
          electronkminus++;
          break;
        case 1:
          Kminusgamma++;
          break;
        case 4:
          neutrinokminus++;
          break;
        default:
        {
          otherkminus++;
          break;
        }
        }
      }
    }

    Bool_t base_ksl = (pi0phi == 0 && piplusphi == 0 && piminusphi == 0 && otherphi == 0 && otherks == 0 && otherkl == 0 && Ksregen == 0 && Kl == 1 && Ks == 1);

    Bool_t signal_cond = (base_ksl && Ksgamma == 0 && Klgamma == 0 &&
                          positronkl + positronks == 0 && electronkl + electronks == 0 && muonminuskl + muonminusks == 0 &&
                          muonpluskl + muonplusks == 0 && neutrinokl + neutrinoks == 0 && neutrinokplus + neutrinokminus == 0 &&
                          ((pi0ks == 2 && pipluskl == 1 && piminuskl == 1 && pi0kl == 0 && piplusks == 0 && piminusks == 0) ||
                           (pi0kl == 2 && piplusks == 1 && piminusks == 1 && pi0ks == 0 && pipluskl == 0 && piminuskl == 0)));

    Bool_t pipi_cond = (base_ksl && Ksgamma == 0 && Klgamma == 0 &&
                        positronkl + positronks == 0 && electronkl + electronks == 0 && muonminuskl + muonminusks == 0 &&
                        muonpluskl + muonplusks == 0 &&
                        ((piminusks == 1 && piplusks == 1 && pipluskl == 1 && piminuskl == 1 && pi0kl == 0 && pi0ks == 0)));

    Bool_t regen_cond = (pi0phi == 0 && piplusphi == 0 && piminusphi == 0 && otherphi == 0 && otherks == 0 && otherkl == 0 && 
                         Ksregen == 1 && Ks == 1 && Kl == 1 && Ksgamma == 0 && Klgamma == 0 && 
                         positronkl + positronks == 0 && electronkl + electronks == 0 && muonminuskl + muonminusks == 0 &&
                         muonpluskl + muonplusks == 0 && neutrinokl + neutrinoks == 0 && neutrinokplus + neutrinokminus == 0 &&
                        (pi0ks == 2 && piplusks == 1 && piminusks == 1 && pi0kl == 0 && pipluskl == 0 && piminuskl == 0));

    Bool_t omega_cond = (pi0phi == 2 && piplusphi == 1 && piminusphi == 1 && otherphi == 0 && otherks == 0 && otherkl == 0 &&
                         positronkl + positronks == 0 && electronkl + electronks == 0 && muonminuskl + muonminusks == 0 &&
                         muonpluskl + muonplusks == 0 && Ksregen == 0 &&
                         Ks == 0 && Kl == 0 && pi0ks == 0 && pi0kl == 0 && pipluskl + piplusks == 0 && piminuskl + piminusks == 0);

    Bool_t three_cond = (base_ksl && Ksgamma == 0 && Klgamma == 0 &&
                         positronkl + positronks == 0 && electronkl + electronks == 0 && muonminuskl + muonminusks == 0 &&
                         muonpluskl + muonplusks == 0 && (pi0kl == 3 && piplusks == 1 && piminusks == 1 && pi0ks == 0 && pipluskl == 0 && piminuskl == 0 && Klgamma == 0 && Ksgamma == 0));

    Bool_t semi_cond = (base_ksl && Ksgamma == 0 && Klgamma == 0 &&
                        ((pi0ks == 2 && positronkl == 1 && piminuskl == 1 && pi0kl == 0 && neutrinokl >= 0) ||
                         (pi0ks == 2 && pipluskl == 1 && electronkl == 1 && pi0kl == 0 && neutrinokl >= 0) ||
                         (pi0ks == 2 && pipluskl == 1 && muonminuskl == 1 && pi0kl == 0 && neutrinokl >= 0) ||
                         (pi0ks == 2 && piminuskl == 1 && muonpluskl == 1 && pi0kl == 0 && neutrinokl >= 0) ||
                         (pi0kl == 2 && positronks == 1 && piminusks == 1 && pi0ks == 0 && neutrinoks >= 0) ||
                         (pi0kl == 2 && piplusks == 1 && electronks == 1 && pi0ks == 0 && neutrinoks >= 0) ||
                         (pi0kl == 2 && piplusks == 1 && muonminusks == 1 && pi0ks == 0 && neutrinoks >= 0) ||
                         (pi0kl == 2 && piminusks == 1 && muonplusks == 1 && pi0ks == 0 && neutrinoks >= 0)));

    Bool_t semi_ele_pos_cond = (base_ksl && Ksgamma == 0 && Klgamma == 0 &&
                                ((pi0ks == 2 && positronkl == 1 && piminuskl == 1 && pi0kl == 0 && neutrinokl >= 0) ||
                                 (pi0ks == 2 && pipluskl == 1 && electronkl == 1 && pi0kl == 0 && neutrinokl >= 0) ||
                                 (pi0kl == 2 && positronks == 1 && piminusks == 1 && pi0ks == 0 && neutrinoks >= 0) ||
                                 (pi0kl == 2 && piplusks == 1 && electronks == 1 && pi0ks == 0 && neutrinoks >= 0)));

    Bool_t semi_muon_cond = (base_ksl && Ksgamma == 0 && Klgamma == 0 &&
                             ((pi0ks == 2 && pipluskl == 1 && muonminuskl == 1 && pi0kl == 0 && neutrinokl >= 0) ||
                              (pi0ks == 2 && piminuskl == 1 && muonpluskl == 1 && pi0kl == 0 && neutrinokl >= 0) ||
                              (pi0kl == 2 && piplusks == 1 && muonminusks == 1 && pi0ks == 0 && neutrinoks >= 0) ||
                              (pi0kl == 2 && piminusks == 1 && muonplusks == 1 && pi0ks == 0 && neutrinoks >= 0)));

    // Conditions for other_flag and isr_flag

    Bool_t base_kplusminus = (pi0phi == 0 && piplusphi == 0 && piminusphi == 0 && otherphi == 0 && otherks == 0 && otherkl == 0 && otherkplus == 0 && otherkminus == 0 && Ksregen == 0 && Kplus == 1 && Kminus == 1);

    Bool_t kplusminus_pipluspi0_piminuspi0 = base_kplusminus && (pi0kplus == 1 && pipluskplus == 1 && pi0kminus == 1 && piminuskminus == 1);

    Bool_t kplusminus_muplus_muminus = base_kplusminus && (muonpluskplus == 1 && muonminuskminus == 1);

    Bool_t kplusminus_pipluspi0_muminus = base_kplusminus && (pi0kplus == 1 && pipluskplus == 1 && muonminuskminus == 1);

    Bool_t kplusminus_pipluspi0pi0_muminus = base_kplusminus && (pi0kplus == 2 && pipluskplus == 1 && muonminuskminus == 1);

    Bool_t kplusminus_pipluspi0pi0_piminuspi0pi0 = base_kplusminus && (pi0kplus == 2 && pipluskplus == 1 && pi0kminus == 2 && piminuskminus == 1);

    Bool_t kplusminus_mupluspi0_piminuspi0pi0 = base_kplusminus && (muonpluskplus == 1 && pi0kplus == 1 && pi0kminus == 2 && piminuskminus == 1);

    Bool_t kplusminus_pipluspi0_piminuspi0pi0 = base_kplusminus && (pi0kplus == 1 && pipluskplus == 1 && pi0kminus == 2 && piminuskminus == 1);

    Bool_t kplusminus_pipluspi0pi0_piminuspi0 = base_kplusminus && (pi0kplus == 2 && pipluskplus == 1 && pi0kminus == 1 && piminuskminus == 1);

    Bool_t kplusminus_muplus_piminuspi0 = base_kplusminus && (muonpluskplus == 1 && pi0kminus == 1 && piminuskminus == 1);

    Bool_t kplusminus_muplus_piminuspi0pi0 = base_kplusminus && (muonpluskplus == 1 && pi0kminus == 2 && piminuskminus == 1);

    Bool_t kplusminus_positron_piminuspi0 = base_kplusminus && (positronkplus == 1 && pi0kminus == 1 && piminuskminus == 1);

    Bool_t kplusminus_pi0pi0piplus_pi0electronpositron = base_kplusminus && ((pipluskplus == 1 && pi0kplus == 2 && pi0kminus == 1 && electronkminus == 1) || (piminuskminus == 1 && pi0kminus == 2 && pi0kplus == 1 && positronkplus == 1) || (pipluskplus == 1 && pi0kplus == 2 && pi0kminus == 1 && muonminuskminus == 1) || (piminuskminus == 1 && pi0kminus == 2 && pi0kplus == 1 && muonpluskplus == 1));

    Bool_t kplusminus_pi0piplus_pi0electronpositron = base_kplusminus && ((pipluskplus == 1 && pi0kplus == 1 && pi0kminus == 1 && electronkminus == 1) || (piminuskminus == 1 && pi0kminus == 1 && pi0kplus == 1 && positronkplus == 1) || (pipluskplus == 1 && pi0kplus == 1 && pi0kminus == 1 && muonminuskminus == 1) || (piminuskminus == 1 && pi0kminus == 1 && pi0kplus == 1 && muonpluskplus == 1));

    Bool_t ksl_3pi0_pipluspiminusgamma = base_ksl && ((pi0kl == 3 && piplusks == 1 && piminusks == 1 && Ksgamma == 1) || (pi0ks == 3 && pipluskl == 1 && piminuskl == 1 && Klgamma == 1));

    Bool_t ksl_3pi0_semileptonic = base_ksl && ((pi0kl == 3 && positronks == 1 && piminusks == 1) || (pi0kl == 3 && electronks == 1 && piplusks == 1) || (pi0ks == 3 && positronkl == 1 && piminuskl == 1) || (pi0ks == 3 && electronkl == 1 && pipluskl == 1) || (pi0kl == 3 && muonplusks == 1 && piminusks == 1) || (pi0kl == 3 && muonminusks == 1 && piplusks == 1) || (pi0ks == 3 && muonpluskl == 1 && piminuskl == 1) || (pi0ks == 3 && muonminuskl == 1 && pipluskl == 1));

    Bool_t ksl_3pi0_2pi0 = base_ksl && ((pi0kl == 3 && pi0ks == 2) || (pi0ks == 3 && pi0kl == 2));

    Bool_t ksl_pipluspiminuspi0 = base_ksl && ((pipluskl == 1 && piminuskl == 1 && pi0kl == 1) || (piplusks == 1 && piminusks == 1 && pi0ks == 1));

    Bool_t ksl_lambda = base_ksl && (Kllambda == 1 || Kslambda == 1);

    Bool_t omegapi0_pipluspiminuspi0 = (pi0phi == 1 && piplusphi == 1 && piminusphi == 1 && otherphi == 0);

    if (gammaphi > 0)
    {
      isr_flag = 1;
    }

    if (signal_cond)
      mctruth_int = 1; // Signal channel: KSKL -> pi+pi-pi0pi0
    else if (regen_cond)
      mctruth_int = 2; // Regeneration
    else if (omega_cond)
      mctruth_int = 3; // omega
    else if (three_cond)
      mctruth_int = 4; // 3pi0
    else if (semi_cond)
      mctruth_int = 5; // semi
    else if (pipi_cond)
      mctruth_int = 7; // pipi
    else
      mctruth_int = 6; // Other background

    if (semi_ele_pos_cond)
      semileptonic_flag = 23; // Semi-leptonic with electron-positron
    else if (semi_muon_cond)
      semileptonic_flag = 56; // Semi-leptonic with muon
    else
      semileptonic_flag = 0; // Not semi-leptonic

    if (kplusminus_pipluspi0_piminuspi0)
      other_flag = 1;
    else if (kplusminus_muplus_muminus)
      other_flag = 2;
    else if (kplusminus_pipluspi0_muminus)
      other_flag = 3;
    else if (kplusminus_muplus_piminuspi0)
      other_flag = 4;
    else if (ksl_3pi0_pipluspiminusgamma)
      other_flag = 5;
    else if (ksl_pipluspiminuspi0)
      other_flag = 6;
    else if (ksl_lambda)
      other_flag = 7;
    else if (kplusminus_pipluspi0pi0_piminuspi0pi0)
      other_flag = 8;
    else if (kplusminus_pipluspi0_piminuspi0pi0)
      other_flag = 9;
    else if (kplusminus_pipluspi0pi0_piminuspi0)
      other_flag = 10;
    else if (kplusminus_positron_piminuspi0)
      other_flag = 11;
    else if (kplusminus_mupluspi0_piminuspi0pi0)
      other_flag = 12;
    else if (omegapi0_pipluspiminuspi0)
      other_flag = 13;
    else if (kplusminus_muplus_piminuspi0pi0)
      other_flag = 14;
    else if (kplusminus_pipluspi0pi0_muminus)
      other_flag = 15;
    else if (ksl_3pi0_semileptonic)
      other_flag = 16;
    else if (ksl_3pi0_2pi0)
      other_flag = 17;
    else if (kplusminus_pi0pi0piplus_pi0electronpositron)
      other_flag = 18;
    else if (kplusminus_pi0piplus_pi0electronpositron)
      other_flag = 19;
    else
      other_flag = 20; // Not recognized

    // ---------------------------------------------------
    // Looking for primary KSKL decays
    // ---------------------------------------------------

    int piplusks_primary = 0, piminusks_primary = 0,
        pipluskl_primary = 0, piminuskl_primary = 0,
        pi0ks_primary = 0, pi0kl_primary = 0,
        muonplusks_primary = 0, muonminusks_primary = 0,
        positronks_primary = 0, electronks_primary = 0,
        neutrinoks_primary = 0,
        muonpluskl_primary = 0, muonminuskl_primary = 0,
        positronkl_primary = 0, electronkl_primary = 0,
        neutrinokl_primary = 0,
        otherks_primary = 0,
        otherkl_primary = 0;

    for (int j = 0; j < ntmc; j++)
    {
      if (!IsPrimaryFromPhi(kine[j]))
        continue;

      if (motherOld[vtxmcOld[j] - 1] == 16)
      {
        switch (pidmcOld[j])
        {
        case 7:
          pi0ks_primary++;
          break;
        case 8:
          piplusks_primary++;
          break;
        case 9:
          piminusks_primary++;
          break;
        case 5:
          muonplusks_primary++;
          break;
        case 6:
          muonminusks_primary++;
          break;
        case 2:
          positronks_primary++;
          break;
        case 3:
          electronks_primary++;
          break;
        case 4:
          neutrinoks_primary++;
          break;
        default:
          otherks_primary++;
          break;
        }
      }
      else if (motherOld[vtxmcOld[j] - 1] == 10)
      {
        switch (pidmcOld[j])
        {
        case 7:
          pi0kl_primary++;
          break;
        case 8:
          pipluskl_primary++;
          break;
        case 9:
          piminuskl_primary++;
          break;
        case 5:
          muonpluskl_primary++;
          break;
        case 6:
          muonminuskl_primary++;
          break;
        case 2:
          positronkl_primary++;
          break;
        case 3:
          electronkl_primary++;
          break;
        case 4:
          neutrinokl_primary++;
          break;
        default:
          otherkl_primary++;
          break;
        }
      }
    }

    bool primary_ks_pm_decay = (piplusks_primary == 1 && piminusks_primary == 1 && otherks_primary == 0 && pi0ks_primary == 0 && muonplusks_primary == 0 && muonminusks_primary == 0 && positronks_primary == 0 && electronks_primary == 0 && neutrinoks_primary == 0);

    bool primary_kl_pm_decay = (pipluskl_primary == 1 && piminuskl_primary == 1 && otherkl_primary == 0 && pi0kl_primary == 0 && muonpluskl_primary == 0 && muonminuskl_primary == 0 && positronkl_primary == 0 && electronkl_primary == 0 && neutrinokl_primary == 0);

    bool primary_ks_00_decay = (pi0ks_primary == 2 && otherks_primary == 0 && piplusks_primary == 0 && piminusks_primary == 0 && muonplusks_primary == 0 && muonminusks_primary == 0 && positronks_primary == 0 && electronks_primary == 0 && neutrinoks_primary == 0);

    bool primary_kl_00_decay = (pi0kl_primary == 2 && otherkl_primary == 0 && pipluskl_primary == 0 && piminuskl_primary == 0 && muonpluskl_primary == 0 && muonminuskl_primary == 0 && positronkl_primary == 0 && electronkl_primary == 0 && neutrinokl_primary == 0);

    bool primary_ks_000_decay = (pi0ks_primary == 3 && otherks_primary == 0 && piplusks_primary == 0 && piminusks_primary == 0 && muonplusks_primary == 0 && muonminusks_primary == 0 && positronks_primary == 0 && electronks_primary == 0 && neutrinoks_primary == 0);

    bool primary_kl_000_decay = (pi0kl_primary == 3 && otherkl_primary == 0 && pipluskl_primary == 0 && piminuskl_primary == 0 && muonpluskl_primary == 0 && muonminuskl_primary == 0 && positronkl_primary == 0 && electronkl_primary == 0 && neutrinokl_primary == 0);

    bool primary_ks_semileptonic_decay = (((piplusks_primary == 1 && (muonminusks_primary == 1 || electronks_primary == 1)) || (piminusks_primary == 1 && (muonplusks_primary == 1 || positronks_primary == 1))) && neutrinoks_primary >= 0 && pi0ks_primary == 0 && otherks_primary == 0);

    bool primary_kl_semileptonic_decay = (((pipluskl_primary == 1 && (muonminuskl_primary == 1 || electronkl_primary == 1)) || (piminuskl_primary == 1 && (muonpluskl_primary == 1 || positronkl_primary == 1))) && neutrinokl_primary >= 0 && pi0kl_primary == 0 && otherkl_primary == 0);

    UInt_t channelFlags = KLOE::ChannelFlags::kNone;
    if (primary_ks_pm_decay || primary_kl_pm_decay)
      channelFlags |= KLOE::ChannelFlags::kHasPM;
    if (primary_ks_00_decay || primary_kl_00_decay)
      channelFlags |= KLOE::ChannelFlags::kHas00;
    if (primary_ks_000_decay || primary_kl_000_decay)
      channelFlags |= KLOE::ChannelFlags::kHas000;
    if (primary_ks_semileptonic_decay || primary_kl_semileptonic_decay)
      channelFlags |= KLOE::ChannelFlags::kHasSemiLeptonic;

    has_pm = (channelFlags & KLOE::ChannelFlags::kHasPM) != 0;
    has_00 = (channelFlags & KLOE::ChannelFlags::kHas00) != 0;
    has_000 = (channelFlags & KLOE::ChannelFlags::kHas000) != 0;
    has_semileptonic = (channelFlags & KLOE::ChannelFlags::kHasSemiLeptonic) != 0;

    // -------------------------------------------------------------------------
    bool has_pm00_primary = (channelFlags & (KLOE::ChannelFlags::kHasPM | KLOE::ChannelFlags::kHas00)) == (KLOE::ChannelFlags::kHasPM | KLOE::ChannelFlags::kHas00);
    bool has_pm000_primary = (channelFlags & (KLOE::ChannelFlags::kHasPM | KLOE::ChannelFlags::kHas000)) == (KLOE::ChannelFlags::kHasPM | KLOE::ChannelFlags::kHas000);
    bool has_semileptonic00_primary = (channelFlags & (KLOE::ChannelFlags::kHasSemiLeptonic | KLOE::ChannelFlags::kHas00)) == (KLOE::ChannelFlags::kHasSemiLeptonic | KLOE::ChannelFlags::kHas00);

    // -------------------------------------------------------------------------

    ErrorHandling::ErrorCodes err = ErrorHandling::ErrorCodes::IMPROPER_MCTRUTH;
    if (signal_cond && !has_pm00_primary)
      LOG_CRITICAL(logger, err, "signal_cond bez oczekiwanej topologii pm+00", ErrorHandling::LogFiles::LogType::ERROR);

    if (three_cond && !has_pm000_primary)
      LOG_CRITICAL(logger, err, "three_cond bez has_000", ErrorHandling::LogFiles::LogType::ERROR);

    if (semi_cond && !has_semileptonic00_primary)
      LOG_CRITICAL(logger, err, "semi_cond bez oczekiwanej topologii semileptonic00", ErrorHandling::LogFiles::LogType::ERROR);
  }
  else if (mcflag == 0)
  {
    mctruth_int = 0; // Data event

    has_pm = false;
    has_00 = false;
    has_000 = false;
    has_semileptonic = false;
  }
}

ErrorHandling::ErrorCodes GeneratedVariables::FindNeutralCluster(
    Int_t nclu,
    Int_t ntcl,
    const Int_t *asscl,
    Int_t NCLMIN,
    ErrorHandling::ErrorLogs &logger,
    std::vector<Int_t> &neuclulist)
{
  neuclulist.clear();
  Int_t neucluind = 0;
  for (Int_t i = 1; i <= nclu; ++i)
  { // Fortran: 1-based
    Bool_t neuclu = true;
    for (Int_t j = 1; j <= ntcl; ++j)
    {
      if (asscl[j - 1] == i)
      { // C++: 0-based
        neuclu = false;
        break;
      }
    }
    if (neuclu)
    {
      ++neucluind;
      neuclulist.push_back(i);
    }
  }
  if (neucluind < NCLMIN)
  {
    auto err = ErrorHandling::ErrorCodes::NOT_RECOGNIZED;

    if (NCLMIN == 4)
      err = ErrorHandling::ErrorCodes::LESS_THAN_FOUR_NEUTRAL_CLUSTERS;
    else if (NCLMIN == 6)
      err = ErrorHandling::ErrorCodes::LESS_THAN_SIX_NEUTRAL_CLUSTERS;

    return err;
  }
  return ErrorHandling::ErrorCodes::NO_ERROR;
}

ErrorHandling::ErrorCodes GeneratedVariables::genVars(
    Int_t ntmc,
    Int_t nvtxmc,
    Int_t nclu,
    const Int_t *pidmc,
    const Int_t *vtxmc,
    const Int_t *mother,
    const Double_t *xvmc,
    const Double_t *yvmc,
    const Double_t *zvmc,
    const Double_t *pxmc,
    const Double_t *pymc,
    const Double_t *pzmc,
    Int_t mcflag,
    Int_t mctruth,
    std::vector<Double_t> &ipmc,
    std::vector<Double_t> &Knemc,
    std::vector<Double_t> &Kchmc,
    std::vector<std::vector<Double_t>> &trkMC,
    const Int_t numberOfClusters,
    std::vector<std::vector<Double_t>> &pgammaMC,
    std::vector<Double_t> &CurvMC,
    std::vector<Double_t> &PhivMC,
    std::vector<Double_t> &CotvMC,
    std::vector<Int_t> &good_clus_ind,
    std::vector<std::vector<Double_t>> cluster_rec,
    Int_t &muonAlertPlus,
    Int_t &muonAlertMinus)
{
  Double_t
      Ks[9] = {},
      Kl[9] = {},
      neu_vtx[3] = {},
      clus_diff_min;
  Int_t
      count = 0,
      ind_gam[4] = {},
      min_ind[4] = {};

  std::vector<Bool_t>
      clus_time;

  const Int_t
      max_count = TMath::Factorial(numberOfClusters);

  std::vector<Double_t>
      clus_diff,
      region,
      cluster;

  std::vector<Int_t>
      mc_ind(numberOfClusters);

  KLOE::CylinderIntersection CylIndObj;

  ipmc.clear();
  Knemc.clear();
  Kchmc.clear();
  trkMC.clear();

  cluster.clear();
  region.clear();
  clus_diff.clear();

  pgammaMC.clear();
  CurvMC.clear();
  PhivMC.clear();
  CotvMC.clear();

  ipmc.resize(3);
  Knemc.resize(9);
  Kchmc.resize(9);

  cluster.resize(3);

  IPGenerated(nvtxmc, mother, ipmc, xvmc, yvmc, zvmc);

  if (mctruth == 1 || mctruth == 2 || mctruth == 4 || mctruth == 5 || mctruth == 7)
  {
    KSLGenerated(nvtxmc, mother, Kl, xvmc, yvmc, zvmc, Ks, ntmc, pidmc, pxmc, pymc, pzmc);

    twoTracksFinder(ntmc, mother, vtxmc, pidmc, Knemc, Kl, Kchmc, Ks, trkMC, pxmc, pymc, pzmc, mctruth, CurvMC, PhivMC, CotvMC, muonAlertPlus, muonAlertMinus);

    if (mctruth == 1 || mctruth == 2 || mctruth == 4 || mctruth == 5)
    {
      ClusterVariableFinder(ntmc, mother, vtxmc, pidmc, pgammaMC, count, pxmc, pymc, pzmc, neu_vtx, Knemc, region, CylIndObj, cluster, ipmc);
    }
  };

  count = 0;

  return ErrorHandling::ErrorCodes::NO_ERROR;
}
void GeneratedVariables::ClusterVariableFinder(Int_t ntmc, const Int_t *mother, const Int_t *vtxmc, const Int_t *pidmc, std::vector<std::vector<Double_t>> &pgammaMC, Int_t &count, const Double_t *pxmc, const Double_t *pymc, const Double_t *pzmc, Double_t neu_vtx[3], std::vector<Double_t> &Knemc, std::vector<Double_t> &region, KLOE::CylinderIntersection &CylIndObj, std::vector<Double_t> &cluster, std::vector<Double_t> &ipmc)
{
  for (Int_t j = 0; j < ntmc; j++)
  {
    if ((mother[vtxmc[j] - 1] == 7) && pidmc[j] == 1)
    {
      Double_t auxEne = std::sqrt(std::pow(pxmc[j], 2) + std::pow(pymc[j], 2) + std::pow(pzmc[j], 2));
      std::vector<Double_t> pgammaAux = {pxmc[j], pymc[j], pzmc[j], auxEne};

      neu_vtx[0] = Knemc[6];
      neu_vtx[1] = Knemc[7];
      neu_vtx[2] = Knemc[8];

      region.push_back(CylIndObj.inter_point(pgammaAux.data(), neu_vtx, cluster.data()));

      Double_t
          beta_c = PhysicsConstants::cVel * Knemc[4] / Knemc[3],
          length = std::sqrt(std::pow(Knemc[6] - ipmc[0], 2) +
                             std::pow(Knemc[7] - ipmc[1], 2) +
                             std::pow(Knemc[8] - ipmc[2], 2)),
          time_K = length / beta_c,
          length_clus = std::sqrt(std::pow(cluster[0] - Knemc[6], 2) +
                                  std::pow(cluster[1] - Knemc[7], 2) +
                                  std::pow(cluster[2] - Knemc[8], 2));

      Double_t auxTim = time_K + (length_clus / PhysicsConstants::cVel);
      std::vector<Double_t> auxiliaryVec = {pxmc[j], pymc[j], pzmc[j], auxEne, cluster[0], cluster[1], cluster[2], auxTim};

      pgammaMC.push_back(auxiliaryVec);
    }
  }
}

void GeneratedVariables::GeneratedClusterFinder(Int_t nclu, Int_t ind_gam[4], const Int_t max_count, std::vector<Double_t> &clus_diff, std::vector<std::vector<Double_t>> &cluster_rec, std::vector<std::vector<Double_t>> &pgammaMC, Int_t mc_ind[4], std::vector<bool> &clus_time, Int_t min_ind[4], Double_t &clus_diff_min, std::vector<Int_t> &good_clus_ind)
{
  for (Int_t j1 = 0; j1 < nclu - 3; j1++)
    for (Int_t j2 = j1 + 1; j2 < nclu - 2; j2++)
      for (Int_t j3 = j2 + 1; j3 < nclu - 1; j3++)
        for (Int_t j4 = j3 + 1; j4 < nclu; j4++)
        {
          ind_gam[0] = j1;
          ind_gam[1] = j2;
          ind_gam[2] = j3;
          ind_gam[3] = j4;

          for (Int_t k = 0; k < max_count; k++)
          {

            clus_diff.push_back(std::sqrt(std::pow(cluster_rec[0][ind_gam[0]] - pgammaMC[mc_ind[0]][4], 2) +
                                          std::pow(cluster_rec[1][ind_gam[0]] - pgammaMC[mc_ind[0]][5], 2) +
                                          std::pow(cluster_rec[2][ind_gam[0]] - pgammaMC[mc_ind[0]][6], 2)) +
                                std::sqrt(std::pow(cluster_rec[0][ind_gam[1]] - pgammaMC[mc_ind[1]][4], 2) +
                                          std::pow(cluster_rec[1][ind_gam[1]] - pgammaMC[mc_ind[1]][5], 2) +
                                          std::pow(cluster_rec[2][ind_gam[1]] - pgammaMC[mc_ind[1]][6], 2)) +
                                std::sqrt(std::pow(cluster_rec[0][ind_gam[2]] - pgammaMC[mc_ind[2]][4], 2) +
                                          std::pow(cluster_rec[1][ind_gam[2]] - pgammaMC[mc_ind[2]][5], 2) +
                                          std::pow(cluster_rec[2][ind_gam[2]] - pgammaMC[mc_ind[2]][6], 2)) +
                                std::sqrt(std::pow(cluster_rec[0][ind_gam[3]] - pgammaMC[mc_ind[3]][4], 2) +
                                          std::pow(cluster_rec[1][ind_gam[3]] - pgammaMC[mc_ind[3]][5], 2) +
                                          std::pow(cluster_rec[2][ind_gam[3]] - pgammaMC[mc_ind[3]][6], 2)));

            clus_time.push_back(pgammaMC[mc_ind[0]][7] > 0. && pgammaMC[mc_ind[1]][7] > 0. && pgammaMC[mc_ind[2]][7] > 0. && pgammaMC[mc_ind[3]][7] > 0.);

            std::next_permutation(mc_ind, mc_ind + 4);
          }

          TMath::Sort(max_count, clus_diff.data(), min_ind, kFALSE);

          if (clus_diff_min > clus_diff[min_ind[0]])
          {
            clus_diff_min = clus_diff[min_ind[0]];

            good_clus_ind[0] = ind_gam[0];
            good_clus_ind[1] = ind_gam[1];
            good_clus_ind[2] = ind_gam[2];
            good_clus_ind[3] = ind_gam[3];
          }
        }
}

void GeneratedVariables::MCvsReconstructedClustersComparator(const std::vector<Int_t> neuclulist, const std::vector<Int_t> gtaken, const std::vector<Int_t> Pnum1, const Int_t ntmc, const std::vector<Int_t> mother, const std::vector<Int_t> vtxmc, const std::vector<Int_t> pidmc, const std::vector<Int_t> kine, const std::vector<Int_t> kinmom, std::vector<Int_t> &goodCluster)
{
  goodCluster.clear();

  for (Int_t i = 0; i < ntmc; ++i)
  {
    for (Int_t j = 0; j < gtaken.size(); ++j)
    {
      if (kine[i] == Pnum1[neuclulist[gtaken[j]] - 1])
      {
        if (mother[vtxmc[i] - 1] == 7 && pidmc[i] == 1)
        {
          Int_t kinPi0 = kinmom[vtxmc[i] - 1];
          for (Int_t k = 0; k < ntmc; ++k)
          {
            if (pidmc[k] == 7 && kine[k] == kinPi0 && (mother[vtxmc[k] - 1] == 10 || mother[vtxmc[k] - 1] == 16))
            {
              const Int_t clusterIndex = neuclulist[gtaken[j]];

              if (std::find(goodCluster.begin(), goodCluster.end(), clusterIndex) == goodCluster.end())
              {
                goodCluster.push_back(clusterIndex);
              }

              break;
            }
          }
        }
      }
    }
  }
}

void GeneratedVariables::twoTracksFinder(Int_t ntmc, const Int_t *mother, const Int_t *vtxmc, const Int_t *pidmc, std::vector<Double_t> &Knemc, Double_t Kl[9], std::vector<Double_t> &Kchmc, Double_t Ks[9], std::vector<std::vector<Double_t>> &trkMC, const Double_t *pxmc, const Double_t *pymc, const Double_t *pzmc, Int_t mctruth, std::vector<Double_t> &CurvMC, std::vector<Double_t> &PhivMC, std::vector<Double_t> &CotvMC, Int_t &muonAlertPlus, Int_t &muonAlertMinus)
{
  muonAlertPlus = false;
  muonAlertMinus = false;

  for (Int_t j = 0; j < ntmc; j++)
  {
    if (mother[vtxmc[j] - 1] == 10)
    {
      if (mctruth != 7)
      {
        if (pidmc[j] == 7)
        {
          std::copy(Ks, Ks + 9, Kchmc.begin());
          std::copy(Kl, Kl + 9, Knemc.begin());
        }
        else if (pidmc[j] == 8 || pidmc[j] == 9)
        {
          std::copy(Ks, Ks + 9, Knemc.begin());
          std::copy(Kl, Kl + 9, Kchmc.begin());
        }
      }
      else
      {
        std::copy(Ks, Ks + 9, Kchmc.begin());
        std::copy(Kl, Kl + 9, Knemc.begin());
      }
    }

    if ((mother[vtxmc[j] - 1] == 10) && (pidmc[j] == 8 || pidmc[j] == 9))
    {
      Double_t auxEne = std::sqrt(std::pow(pxmc[j], 2) +
                                  std::pow(pymc[j], 2) +
                                  std::pow(pzmc[j], 2) +
                                  std::pow(PhysicsConstants::mPiCh, 2));
      std::vector<Double_t> auxiliaryVec = {pxmc[j], pymc[j], pzmc[j], auxEne, 10};

      trkMC.push_back(auxiliaryVec);

      CurvMC.push_back(1000. / std::sqrt(std::pow(pxmc[j], 2) + std::pow(pymc[j], 2)));
      PhivMC.push_back(atan2(pymc[j], pxmc[j]));
      CotvMC.push_back(pzmc[j] / std::sqrt(std::pow(pxmc[j], 2) + std::pow(pymc[j], 2)));

      if (pidmc[j] == 9)
      {
        CurvMC.back() = -CurvMC.back();
      }

      for (Int_t k = 0; k < ntmc; k++)
      {
        if ((mother[vtxmc[k] - 1] == 8) && (pidmc[k] == 5 || pidmc[k] == 6))
        {
          muonAlertPlus = true;
        }

        if ((mother[vtxmc[k] - 1] == 9) && (pidmc[k] == 5 || pidmc[k] == 6))
        {
          muonAlertMinus = true;
        }
      }
    }

    if ((mother[vtxmc[j] - 1] == 16) && (pidmc[j] == 8 || pidmc[j] == 9))
    {
      Double_t auxEne = std::sqrt(std::pow(pxmc[j], 2) +
                                  std::pow(pymc[j], 2) +
                                  std::pow(pzmc[j], 2) +
                                  std::pow(PhysicsConstants::mPiCh, 2));
      std::vector<Double_t> auxiliaryVec = {pxmc[j], pymc[j], pzmc[j], auxEne, 16};

      trkMC.push_back(auxiliaryVec);

      CurvMC.push_back(1000. / std::sqrt(std::pow(pxmc[j], 2) + std::pow(pymc[j], 2)));
      PhivMC.push_back(atan2(pymc[j], pxmc[j]));
      CotvMC.push_back(pzmc[j] / std::sqrt(std::pow(pxmc[j], 2) + std::pow(pymc[j], 2)));

      if (pidmc[j] == 9)
      {
        CurvMC.back() = -CurvMC.back();
      }

      for (Int_t k = 0; k < ntmc; k++)
      {
        if ((mother[vtxmc[k] - 1] == 8) && (pidmc[k] == 5 || pidmc[k] == 6))
        {
          muonAlertPlus = true;
        }

        if ((mother[vtxmc[k] - 1] == 9) && (pidmc[k] == 5 || pidmc[k] == 6))
        {
          muonAlertMinus = true;
        }
      }
    }
  }
}

void GeneratedVariables::KSLGenerated(Int_t nvtxmc, const Int_t *mother, Double_t Kl[9], const Double_t *xvmc, const Double_t *yvmc, const Double_t *zvmc, Double_t Ks[9], Int_t ntmc, const Int_t *pidmc, const Double_t *pxmc, const Double_t *pymc, const Double_t *pzmc)
{
  for (Int_t j = 0; j < nvtxmc; j++)
  {
    if (mother[j] == 10)
    {
      Kl[6] = xvmc[j];
      Kl[7] = yvmc[j];
      Kl[8] = zvmc[j];
    }
  }

  for (Int_t j = 0; j < nvtxmc; j++)
  {
    if (mother[j] == 16)
    {
      Ks[6] = xvmc[j];
      Ks[7] = yvmc[j];
      Ks[8] = zvmc[j];
    }
  }

  for (Int_t j = 0; j < ntmc; j++)
  {
    if (pidmc[j] == 10)
    {
      Kl[0] = pxmc[j];
      Kl[1] = pymc[j];
      Kl[2] = pzmc[j];
      Kl[5] = PhysicsConstants::mK0;
      Kl[4] = std::pow(Kl[0], 2) + std::pow(Kl[1], 2) + std::pow(Kl[2], 2);
      Kl[3] = std::sqrt(Kl[4] + std::pow(Kl[5], 2));
      Kl[4] = std::sqrt(Kl[4]);
    }
  }

  for (Int_t j = 0; j < ntmc; j++)
  {
    if (pidmc[j] == 16)
    {
      Ks[0] = pxmc[j];
      Ks[1] = pymc[j];
      Ks[2] = pzmc[j];
      Ks[5] = PhysicsConstants::mK0;
      Ks[4] = std::pow(Ks[0], 2) + std::pow(Ks[1], 2) + std::pow(Ks[2], 2);
      Ks[3] = std::sqrt(Ks[4] + std::pow(Ks[5], 2));
      Ks[4] = std::sqrt(Ks[4]);
    }
  }
}
void GeneratedVariables::IPGenerated(Int_t nvtxmc, const Int_t *mother, std::vector<Double_t> &ipmc, const Double_t *xvmc, const Double_t *yvmc, const Double_t *zvmc)
{
  for (Int_t j = 0; j < nvtxmc; j++)
  {
    if (mother[j] == 50)
    {
      ipmc[0] = xvmc[j];
      ipmc[1] = yvmc[j];
      ipmc[2] = zvmc[j];
    }
  }
};