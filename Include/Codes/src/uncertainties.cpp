#include "uncertainties.h"
#include "TMath.h"
#include <const.h>

#define sigma_coor 1.3

/**
 * @brief Uncertainty of the Enecl in mega-electronvolts (MeV). Dependent of the energy of photon.
 * @param Enecl Cluster's energy in mega-electronvolts (MeV)
 * @return Returns:
 *          - For Enecl > 0. returns the calculated uncertainty
 *          - For Enecl <= 0 returns -999
 */
Double_t clu_ene_error(Float_t Enecl)
{
    Double_t sigma_E = 0.;

    // Error in MeV

    sigma_E = Enecl * 0.057 / sqrt(0.001 * Enecl);

    if (Enecl > 0.)
        return sigma_E;
    else
        return -999;
};

/**
 * @brief Uncertainty of the Tcl in nanoseconds (ns). Dependent of the energy of photon and the offset.
 * @param Enecl Cluster's energy in mega-electronvolts (MeV)
 * @return Returns:
 *          - For Enecl > 0. returns the calculated uncertainty
 *          - For Enecl <= 0 returns -999
 */
Double_t clu_time_error(Float_t Enecl)
{
    Double_t sigma_time = 0.;

    // Error in ns

    sigma_time = sqrt(pow(0.054 / sqrt(0.001 * Enecl), 2) + pow(0.05, 2));

    if (Enecl > 0.)
        return sigma_time;
    else
        return -999;
}

Double_t clu_x_error(Float_t x_coor, Float_t y_coor, Float_t z_coor, Float_t Enecl)
{
    Double_t sigma_x = 0.;

    if (sqrt(pow(x_coor, 2) + pow(y_coor, 2)) > 200. && abs(z_coor) < 165)
        sigma_x = sigma_coor;
    else
        sigma_x = sigma_coor;

    if (Enecl > 0.)
        return sigma_x;
    else
        return -999;
}

Double_t clu_y_error(Float_t x_coor, Float_t y_coor, Float_t z_coor, Float_t Enecl)
{
    Double_t sigma_y = 0.;

    // Error in ns

    if (sqrt(pow(x_coor, 2) + pow(y_coor, 2)) > 200. && abs(z_coor) < 165)
        sigma_y = sigma_coor;
    else
        sigma_y = sigma_coor / sqrt(0.001 * Enecl);

    if (Enecl > 0.)
        return sigma_y;
    else
        return -999;
}

Double_t clu_z_error(Float_t x_coor, Float_t y_coor, Float_t z_coor, Float_t Enecl)
{
    Double_t sigma_z = 0.;

    // Error in ns

    if (sqrt(pow(x_coor, 2) + pow(y_coor, 2)) > 200. && abs(z_coor) < 165)
        sigma_z = sigma_coor / sqrt(0.001 * Enecl);
    else
        sigma_z = sigma_coor;

    if (Enecl > 0.)
        return sigma_z;
    else
        return -999;
}