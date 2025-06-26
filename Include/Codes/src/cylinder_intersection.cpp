#include <cylinder_intersection.h>

namespace KLOE
{
  // Default constructor
  CylinderIntersection::CylinderIntersection(Double_t zmax, Double_t Rmax) : _zmax(zmax), _Rmax(Rmax) { std::cout << "Cylinder intersection finder initialized." << std::endl; };

  Int_t CylinderIntersection::inter_point(Float_t *momentum, Float_t *neu_vtx, Float_t *cluster)
  {
    // Step 1 - checking the intersection of momentum line with barrel

    barrel_inter(momentum, neu_vtx, cluster);

    if (abs(cluster[2]) < _zmax)
    {
      return 1;
    }
    else
    {
      endcap_inter(momentum, neu_vtx, cluster);
      return 2;
    }

    return 0;
  };

  void CylinderIntersection::barrel_inter(Float_t *momentum, Float_t *neu_vtx, Float_t *cluster)
  {
    Float_t transv_vector_pr = 0., mom_val2 = 0., delta = 0., b = 0., a = 0., num[2] = {0., 0.}, den = 0., par[2] = {0., 0.}, tmp[3] = {0., 0., 0.}, clus_mom_scalar_pr = 0., gamma_path = 0., cosangle = -999.;

    transv_vector_pr = momentum[0] * neu_vtx[1] - momentum[1] * neu_vtx[0];
    mom_val2 = pow(momentum[0], 2) + pow(momentum[1], 2);

    b = momentum[0] * neu_vtx[0] + momentum[1] * neu_vtx[1];

    delta = mom_val2 * pow(_Rmax, 2) - pow(transv_vector_pr, 2);

    if (delta < 0)
    {
      cluster[0] = 0.;
      cluster[1] = 0.;
      cluster[2] = 0.;
    }
    else
    {
      delta = sqrt(delta);

      den = mom_val2;

      if (den == 0)
      {
        cluster[0] = 1.;
        cluster[1] = 1.;
        cluster[2] = 1.;
      }
      else
      {
        num[0] = -b - delta;
        num[1] = -b + delta;

        par[0] = num[0] / den;
        par[1] = num[1] / den;

        // Calculation of temporary clusters and cross check with momentum direction

        mom_val2 += pow(momentum[2], 2);
        mom_val2 = sqrt(mom_val2);

        for (Int_t i = 0; i < 2; i++)
        {
          clus_mom_scalar_pr = 0.;
          gamma_path = 0.;

          for (Int_t j = 0; j < 3; j++)
          {
            tmp[j] = neu_vtx[j] + par[i] * momentum[j];

            clus_mom_scalar_pr += (tmp[j] - neu_vtx[j]) * momentum[j];
            gamma_path += pow(tmp[j] - neu_vtx[j], 2);
          }

          gamma_path = sqrt(gamma_path);
          cosangle = clus_mom_scalar_pr / (gamma_path * mom_val2);

          if (cosangle > 0.99)
          {
            cluster[0] = tmp[0];
            cluster[1] = tmp[1];
            cluster[2] = tmp[2];
            break;
          }
        }
      }
    }
  };

  void CylinderIntersection::endcap_inter(Float_t *momentum, Float_t *neu_vtx, Float_t *cluster)
  {
    Float_t z_axis[3] = {0., 0., 1.}, r0[3] = {0., 0., 0.}, num = 0., den = 0., par = 0.;

    if (cluster[2] >= _zmax)
    {
      r0[2] = _zmax;
    }
    else if (cluster[2] <= -_zmax)
    {
      r0[2] = -_zmax;
    }

    num = (r0[0] - neu_vtx[0]) * z_axis[0] +
          (r0[1] - neu_vtx[1]) * z_axis[1] +
          (r0[2] - neu_vtx[2]) * z_axis[2];

    den = momentum[0] * z_axis[0] +
          momentum[1] * z_axis[1] +
          momentum[2] * z_axis[2];

    if (den == 0)
    {
      cluster[0] = 2.;
      cluster[1] = 2.;
      cluster[2] = 2.;
    }
    else
    {
      par = num / den;

      cluster[0] = neu_vtx[0] + momentum[0] * par;
      cluster[1] = neu_vtx[1] + momentum[1] * par;
      cluster[2] = neu_vtx[2] + momentum[2] * par;
    }
  };
}