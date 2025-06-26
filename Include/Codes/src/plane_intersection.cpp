#include "plane_intersection.h"
#include "TMath.h"
#include <const.h>

void plane_intersection(Float_t *pointplane, Float_t *perp_vec/*Float_t *vectorplane1, Float_t *vectorplane2*/,
                        Float_t *pointline, Float_t *vectorline,
                        Float_t *ip_plane)
{
    Float_t n[3], d, l_n, numerator;
    
    n[0] = perp_vec[0];
    n[1] = perp_vec[1];
    n[2] = perp_vec[2];

    l_n = vectorline[0]*n[0] + vectorline[1]*n[1] + vectorline[2]*n[2];

    numerator = (pointplane[0] - pointline[0])*n[0] +
                (pointplane[1] - pointline[1])*n[1] +
                (pointplane[2] - pointline[2])*n[2];

    d = numerator/l_n;

    ip_plane[0] = pointline[0] + vectorline[0]*d;
    ip_plane[1] = pointline[1] + vectorline[1]*d;
    ip_plane[2] = pointline[2] + vectorline[2]*d;
}

void plane_intersection(Double_t *pointplane, Double_t *perp_vec, Double_t *pointline, Double_t *vectorline, Double_t *ip_plane)
{
Double_t n[3], d, l_n, numerator;

n[0] = perp_vec[0];
n[1] = perp_vec[1];
n[2] = perp_vec[2];

l_n = vectorline[0]*n[0] + vectorline[1]*n[1] + vectorline[2]*n[2];

numerator = (pointplane[0] - pointline[0])*n[0] +
(pointplane[1] - pointline[1])*n[1] +
(pointplane[2] - pointline[2])*n[2];

d = numerator/l_n;

ip_plane[0] = pointline[0] + vectorline[0]*d;
ip_plane[1] = pointline[1] + vectorline[1]*d;
ip_plane[2] = pointline[2] + vectorline[2]*d;
}