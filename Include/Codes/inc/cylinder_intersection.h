#ifndef CYLINDER_INTERSECTION_H
#define CYLINDER_INTERSECTION_H

#include <kloe_class.h>

namespace KLOE
{
  class CylinderIntersection : private pm00
  {
    private:
      const Double_t
          _zmax = 165.0,
          _Rmax = 200.0;

      void barrel_inter(Float_t *, Float_t *, Float_t *);
      void endcap_inter(Float_t *, Float_t *, Float_t *);

    public:
      CylinderIntersection(Double_t zmax = 165.0, Double_t rmax = 200.0);

      Double_t getZmax() { return _zmax; };
      Double_t getRmax() { return _Rmax; };

      Int_t inter_point(Float_t *, Float_t *, Float_t *);
  };
}

#endif // !CYLINDER_INTERSECTION_H