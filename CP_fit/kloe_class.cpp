#include "kloe_class.h"
#include <TMath.h>
#include <TLorentzVector.h>

namespace KLOE
{
    pm00::pm00(TLorentzVector *mom_list, TLorentzVector *pos_list)
    {
        phi_mom = mom_list[0];
        phi_pos = pos_list[0];

        kaon_mom[0] = mom_list[1];
        kaon_pos[0] = pos_list[1];

        kaon_mom[1] = mom_list[2];
        kaon_pos[1] = pos_list[2];

        pi_ch_mom[0] = mom_list[3];
        pi_ch_pos[0] = pos_list[3];

        pi_ch_mom[1] = mom_list[4];
        pi_ch_pos[1] = pos_list[4];

        pi_ne_mom[0] = mom_list[5];
        pi_ne_pos[0] = pos_list[5];

        pi_ne_mom[1] = mom_list[6];
        pi_ne_pos[1] = pos_list[6];

        photon_mom[0] = mom_list[7];
        photon_pos[0] = pos_list[7];

        photon_mom[1] = mom_list[8];
        photon_pos[1] = pos_list[8];

        photon_mom[2] = mom_list[9];
        photon_pos[2] = pos_list[9];

        photon_mom[3] = mom_list[10];
        photon_pos[3] = pos_list[10];
    };

    pm00::pm00()
    {
        phi_mom.SetPxPyPzE(0.,0.,0.,0.);
        phi_pos.SetXYZT(0.,0.,0.,0.);

        kaon_mom[0].SetPxPyPzE(0.,0.,0.,0.);
        kaon_pos[0].SetXYZT(0.,0.,0.,0.);

        kaon_mom[1].SetPxPyPzE(0.,0.,0.,0.);
        kaon_pos[1].SetXYZT(0.,0.,0.,0.);

        pi_ch_mom[0].SetPxPyPzE(0.,0.,0.,0.);
        pi_ch_pos[0].SetXYZT(0.,0.,0.,0.);

        pi_ch_mom[1].SetPxPyPzE(0.,0.,0.,0.);
        pi_ch_pos[1].SetXYZT(0.,0.,0.,0.);

        pi_ne_mom[0].SetPxPyPzE(0.,0.,0.,0.);
        pi_ne_pos[0].SetXYZT(0.,0.,0.,0.);

        pi_ne_mom[1].SetPxPyPzE(0.,0.,0.,0.);
        pi_ne_pos[1].SetXYZT(0.,0.,0.,0.);

        photon_mom[0].SetPxPyPzE(0.,0.,0.,0.);
        photon_pos[0].SetXYZT(0.,0.,0.,0.);

        photon_mom[1].SetPxPyPzE(0.,0.,0.,0.);
        photon_pos[1].SetXYZT(0.,0.,0.,0.);

        photon_mom[2].SetPxPyPzE(0.,0.,0.,0.);
        photon_pos[2].SetXYZT(0.,0.,0.,0.);

        photon_mom[3].SetPxPyPzE(0.,0.,0.,0.);
        photon_pos[3].SetXYZT(0.,0.,0.,0.);
    }

    void pm00::inv_mass_calc(TLorentzVector four_mom)
    {
        inv_mass = four_mom.Mag(); 
    };

    void pm00::angle(TLorentzVector vec1, TLorentzVector vec2)
    {
        angle_vec = vec1.Angle(vec2.Vect()); 
    };

    void pm00::cyl_comp(TLorentzVector vec)
    {
        transv = vec.Perp();
        azim_angle = vec.Phi();
    };

    void pm00::boost_vector(TLorentzVector four_mom)
    {
        boost = four_mom.BoostVector();
    };

    void pm00::lorentz_transf(TLorentzVector four_mom)
    {
        four_mom.Boost(-boost);
    };
}