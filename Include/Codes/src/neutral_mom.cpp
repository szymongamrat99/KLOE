#include "neutral_mom.h"
#include "TMath.h"
#include <const.h>

void neutral_mom(Float_t Xcl, Float_t Ycl, Float_t Zcl, Float_t Enecl, Float_t *neu_vtx, Float_t *mom_vec)
{
    Float_t distance = sqrt(pow(Xcl - neu_vtx[0],2) + pow(Ycl - neu_vtx[1],2) + pow(Zcl - neu_vtx[2],2));

    mom_vec[0] = Enecl*(Xcl - neu_vtx[0])/distance;
    mom_vec[1] = Enecl*(Ycl - neu_vtx[1])/distance;
    mom_vec[2] = Enecl*(Zcl - neu_vtx[2])/distance;
    mom_vec[3] = Enecl;

}

void neutral_mom(Double_t Xcl, Double_t Ycl, Double_t Zcl, Double_t Enecl, Double_t *neu_vtx, Double_t *mom_vec)
{
    Double_t distance = sqrt(pow(Xcl - neu_vtx[0],2) + pow(Ycl - neu_vtx[1],2) + pow(Zcl - neu_vtx[2],2));

    mom_vec[0] = Enecl*(Xcl - neu_vtx[0])/distance;
    mom_vec[1] = Enecl*(Ycl - neu_vtx[1])/distance;
    mom_vec[2] = Enecl*(Zcl - neu_vtx[2])/distance;
    mom_vec[3] = Enecl;

}