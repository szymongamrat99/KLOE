#include <closest_approach.h>
#include <TMath.h>
#include <const.h>

void closest_approach(Float_t point1[3], Float_t vector1[3], 
                      Float_t point2[3], Float_t vector2[3],
                      Float_t ip_closest[3])
{
    Float_t n[3], n_len_2, vector1xn[3], vector2xn[3], 
            numerator1, numerator2, t1, t2;

    n[0] = vector1[1]*vector2[2] - vector1[2]*vector2[1];
    n[1] = vector1[2]*vector2[0] - vector1[0]*vector2[2];
    n[2] = vector1[0]*vector2[1] - vector1[1]*vector2[0];

    n_len_2 = pow(n[0],2) + pow(n[1],2) + pow(n[2],2);

    vector1xn[0] = vector1[1]*n[2] - vector1[2]*n[1];
    vector1xn[1] = vector1[2]*n[0] - vector1[0]*n[2];
    vector1xn[2] = vector1[0]*n[1] - vector1[1]*n[0];

    vector2xn[0] = vector2[1]*n[2] - vector2[2]*n[1];
    vector2xn[1] = vector2[2]*n[0] - vector2[0]*n[2];
    vector2xn[2] = vector2[0]*n[1] - vector2[1]*n[0];

    numerator1 = vector2xn[0]*(point2[0] - point1[0]) +
                 vector2xn[1]*(point2[1] - point1[1]) +
                 vector2xn[2]*(point2[2] - point1[2]);

    numerator2 = vector1xn[0]*(point2[0] - point1[0]) +
                 vector1xn[1]*(point2[1] - point1[1]) +
                 vector1xn[2]*(point2[2] - point1[2]);

    t1 = numerator1/n_len_2;
    t2 = numerator2/n_len_2;

    ip_closest[0] = t1*vector1[0] + point1[0];
    ip_closest[1] = t1*vector1[1] + point1[1];
    ip_closest[2] = t1*vector1[2] + point1[2];
}