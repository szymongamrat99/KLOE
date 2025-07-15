#include <GeneratedVariables.h>

void GeneratedVariables::classifyChannel(
    Int_t ntmc,
    Int_t nvtxmc,
    Int_t* pidmcOld,
    Int_t* vtxmcOld,
    Int_t* motherOld,
    UInt_t mcflag,
    Int_t& mctruth_int
) {
    UInt_t Ks = 0, Kl = 0, Ksregen = 0, piplusks = 0, pipluskl = 0, piminusks = 0, piminuskl = 0,
                 muonplusks = 0, muonpluskl = 0, muonminusks = 0, muonminuskl = 0, electronks = 0, electronkl = 0,
                 positronks = 0, positronkl = 0, pi0ks = 0, pi0kl = 0, pi0phi = 0, piplusphi = 0, piminusphi = 0,
                 otherphi = 0, otherkl = 0, otherks = 0, gammaphi = 0;

    if (mcflag == 1) {
        for (Int_t j = 0; j < ntmc; ++j) {
            if (motherOld[vtxmcOld[j] - 1] == 50) {
                switch (pidmcOld[j]) {
                    case 10: Kl++; break;
                    case 16: Ks++; break;
                    case 7:  pi0phi++; break;
                    case 8:  piplusphi++; break;
                    case 9:  piminusphi++; break;
                    case 1:  gammaphi++; break;
                    default: otherphi++; break;
                }
            } else if (motherOld[vtxmcOld[j] - 1] == 10) {
                switch (pidmcOld[j]) {
                    case 16: Ksregen++; break;
                    case 7:  pi0kl++; break;
                    case 8:  pipluskl++; break;
                    case 9:  piminuskl++; break;
                    case 5:  muonpluskl++; break;
                    case 6:  muonminuskl++; break;
                    case 2:  positronkl++; break;
                    case 3:  electronkl++; break;
                    default: otherkl++; break;
                }
            } else if (motherOld[vtxmcOld[j] - 1] == 16) {
                switch (pidmcOld[j]) {
                    case 7:  pi0ks++; break;
                    case 8:  piplusks++; break;
                    case 9:  piminusks++; break;
                    case 5:  muonplusks++; break;
                    case 6:  muonminusks++; break;
                    case 2:  positronks++; break;
                    case 3:  electronks++; break;
                    default: otherks++; break;
                }
            }
        }

        Bool_t signal_cond = (pi0phi == 0 && piplusphi == 0 && piminusphi == 0 && otherphi == 0 && otherks == 0 && otherkl == 0 &&
                            positronkl + positronks == 0 && electronkl + electronks == 0 && muonminuskl + muonminusks == 0 &&
                            muonpluskl + muonplusks == 0 && Ksregen == 0 &&
                            Ks == 1 && Kl == 1 &&
                            ((pi0ks == 2 && pipluskl == 1 && piminuskl == 1 && pi0kl == 0 && piplusks == 0 && piminusks == 0) ||
                             (pi0kl == 2 && piplusks == 1 && piminusks == 1 && pi0ks == 0 && pipluskl == 0 && piminuskl == 0)));

        Bool_t pipi_cond = (pi0phi == 0 && piplusphi == 0 && piminusphi == 0 && otherphi == 0 && otherks == 0 && otherkl == 0 &&
                          positronkl + positronks == 0 && electronkl + electronks == 0 && muonminuskl + muonminusks == 0 &&
                          muonpluskl + muonplusks == 0 && Ksregen == 0 &&
                          Ks == 1 && Kl == 1 &&
                          ((piminusks == 1 && piplusks == 1 && pipluskl == 1 && piminuskl == 1 && pi0kl == 0 && pi0ks == 0)));

        Bool_t regen_cond = (Ksregen == 1 && Ks == 1 && Kl == 1);

        Bool_t omega_cond = (pi0phi == 2 && piplusphi == 1 && piminusphi == 1 && otherphi == 0 && otherks == 0 && otherkl == 0 &&
                           positronkl + positronks == 0 && electronkl + electronks == 0 && muonminuskl + muonminusks == 0 &&
                           muonpluskl + muonplusks == 0 && Ksregen == 0 &&
                           Ks == 0 && Kl == 0 && pi0ks == 0 && pi0kl == 0 && pipluskl + piplusks == 0 && piminuskl + piminusks == 0);

        Bool_t three_cond = (pi0phi == 0 && piplusphi == 0 && piminusphi == 0 && otherphi == 0 && otherks == 0 && otherkl == 0 &&
                           positronkl + positronks == 0 && electronkl + electronks == 0 && muonminuskl + muonminusks == 0 &&
                           muonpluskl + muonplusks == 0 && Ksregen == 0 &&
                           Ks == 1 && Kl == 1 && (pi0kl == 3 && piplusks == 1 && piminusks == 1 && pi0ks == 0 && pipluskl == 0 && piminuskl == 0));

        Bool_t semi_cond = (pi0phi == 0 && piplusphi == 0 && piminusphi == 0 && otherphi == 0 && otherks == 0 && otherkl == 0 &&
                          Ksregen == 0 && Ks == 1 && Kl == 1 &&
                          ((pi0ks == 2 && positronkl == 1 && piminuskl == 1 && pi0kl == 0) ||
                           (pi0ks == 2 && pipluskl == 1 && electronkl == 1 && pi0kl == 0) ||
                           (pi0ks == 2 && pipluskl == 1 && muonminuskl == 1 && pi0kl == 0) ||
                           (pi0ks == 2 && piminuskl == 1 && muonpluskl == 1 && pi0kl == 0) ||
                           (pi0kl == 2 && positronks == 1 && piminusks == 1 && pi0ks == 0) ||
                           (pi0kl == 2 && piplusks == 1 && electronks == 1 && pi0ks == 0) ||
                           (pi0kl == 2 && piplusks == 1 && muonminusks == 1 && pi0ks == 0) ||
                           (pi0kl == 2 && piminusks == 1 && muonplusks == 1 && pi0ks == 0)));

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
    } else if (mcflag == 0) {
        mctruth_int = 0; // Data event
    }
}

ErrorHandling::ErrorCodes GeneratedVariables::FindNeutralCluster(
    Int_t nclu,
    Int_t ntcl,
    const Int_t* asscl,
    Int_t NCLMIN,
    ErrorHandling::ErrorLogs& logger,
    std::vector<Int_t>& neuclulist
) {
    neuclulist.clear();
    Int_t neucluind = 0;
    for (Int_t i = 1; i <= nclu; ++i) { // Fortran: 1-based
        Bool_t neuclu = true;
        for (Int_t j = 1; j <= ntcl; ++j) {
            if (asscl[j - 1] == i) { // C++: 0-based
                neuclu = false;
                break;
            }
        }
        if (neuclu) {
            ++neucluind;
            neuclulist.push_back(i);
        }
    }
    if (neucluind < NCLMIN) {
        auto err = ErrorHandling::ErrorCodes::NOT_ENOUGH_NEUTRAL_CLUSTERS;
        return err;
    }
    return ErrorHandling::ErrorCodes::NO_ERROR;
}
