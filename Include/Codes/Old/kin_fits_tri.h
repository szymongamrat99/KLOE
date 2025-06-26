#include <TMath.h>
#include <TMatrixD.h>
#include <TVectorT.h>
#include <TF1.h>

#include "kloe_class.h"
#include "reconstructor.h"

namespace KLOE
{

  class kin_fits : public pm00
  {
    public:
      kin_fits(UInt_t N_init, UInt_t M_init, Double_t *par, Double_t *sigmas, Int_t iter_num_init) : _N(N_init), _M(M_init), _iter_num(iter_num_init)
      {
        _P.ResizeTo(_N);
        _P0.ResizeTo(_N);
        _P1.ResizeTo(_N);
        _C.ResizeTo(_M);
        _Lambda.ResizeTo(_M);
        _Corr.ResizeTo(_N);

        _V.ResizeTo(_N, _N);
        _V0.ResizeTo(_N, _N);

        for(Int_t i = 0; i < _N; i++)
        {
          _P(i) = par[i];
          _P0(i) = par[i];

          for(Int_t j = 0; j < _N; j++)
          {
            if(i == j) _V(i,j) = pow(sigmas[i],2);
            else _V(i,j) = 0.;

            if(i == j) _V0(i,j) = pow(sigmas[i],2);
            else _V0(i,j) = 0.;
          }
        }

        _V_updated.ResizeTo(_N, _N);

        _D.ResizeTo(_M, _N);
        _D_T.ResizeTo(_N, _M);

        _Aux.ResizeTo(_M, _M);

        ////////////////////////////////////////////////////////////////////////////////////////////////////////

        _constraints = new TF1*[_M];

        _constraints[0] = new TF1("Total energy conservation", this, &kin_fits::ene_consv_pi, 0, 1, _N);
        _constraints[1] = new TF1("Total x-mom conservation", this, &kin_fits::mom_x_consv_pi, 0, 1, _N);
        _constraints[2] = new TF1("Total y-mom conservation", this, &kin_fits::mom_y_consv_pi, 0, 1, _N);
        _constraints[3] = new TF1("Total z-mom conservation", this, &kin_fits::mom_z_consv_pi, 0, 1, _N);

        _constraints[4] = new TF1("Cluster time 1 conservation", this, &kin_fits::cluster_time_cons_first, 0, 1, _N);
        _constraints[5] = new TF1("Cluster time 2 conservation", this, &kin_fits::cluster_time_cons_second, 0, 1, _N);
        _constraints[6] = new TF1("Cluster time 3 conservation", this, &kin_fits::cluster_time_cons_third, 0, 1, _N);
        _constraints[7] = new TF1("Cluster time 4 conservation", this, &kin_fits::cluster_time_cons_fourth, 0, 1, _N);

        _constraints[8] = new TF1("Neutral inv mass conservation", this, &kin_fits::minv_neu_cons, 0, 1, _N);
        _constraints[9] = new TF1("Charged inv mass conservation", this, &kin_fits::minv_ch_cons, 0, 1, _N);

        /////////////////////////////////////////////////////////////////////////////////////////////////////////

      };

      kin_fits(UInt_t N_init, UInt_t M_init, Int_t iter_num_init) : _N(N_init), _M(M_init), _iter_num(iter_num_init)
      {
        _P.ResizeTo(_N);
        _P0.ResizeTo(_N);
        _P1.ResizeTo(_N);
        _C.ResizeTo(_M);
        _Lambda.ResizeTo(_M);
        _Corr.ResizeTo(_N);

        _V.ResizeTo(_N, _N);
        _V0.ResizeTo(_N, _N);

        _V_updated.ResizeTo(_N, _N);

        _D.ResizeTo(_M, _N);
        _D_T.ResizeTo(_N, _M);

        _Aux.ResizeTo(_M, _M);

        ////////////////////////////////////////////////////////////////////////////////////////////////////////

        _constraints = new TF1*[_M];

        _constraints[0] = new TF1("Total energy conservation", this, &kin_fits::ene_consv_pi, 0, 1, _N);
        _constraints[1] = new TF1("Total x-mom conservation", this, &kin_fits::mom_x_consv_pi, 0, 1, _N);
        _constraints[2] = new TF1("Total y-mom conservation", this, &kin_fits::mom_y_consv_pi, 0, 1, _N);
        _constraints[3] = new TF1("Total z-mom conservation", this, &kin_fits::mom_z_consv_pi, 0, 1, _N);

        _constraints[4] = new TF1("Cluster time 1 conservation", this, &kin_fits::cluster_time_cons_first, 0, 1, _N);
        _constraints[5] = new TF1("Cluster time 2 conservation", this, &kin_fits::cluster_time_cons_second, 0, 1, _N);
        _constraints[6] = new TF1("Cluster time 3 conservation", this, &kin_fits::cluster_time_cons_third, 0, 1, _N);
        _constraints[7] = new TF1("Cluster time 4 conservation", this, &kin_fits::cluster_time_cons_fourth, 0, 1, _N);

        _constraints[8] = new TF1("Neutral inv mass conservation", this, &kin_fits::minv_neu_cons, 0, 1, _N);
        _constraints[9] = new TF1("Charged inv mass conservation", this, &kin_fits::minv_ch_cons, 0, 1, _N);

        /////////////////////////////////////////////////////////////////////////////////////////////////////////

      };

      kin_fits(UInt_t N_init, UInt_t M_init, Int_t iter_num_init, Float_t *neu_vtx[2]) : _N(N_init), _M(M_init), _iter_num(iter_num_init)
      {
        _P.ResizeTo(_N);
        _P0.ResizeTo(_N);
        _P1.ResizeTo(_N);
        _C.ResizeTo(_M);
        _Lambda.ResizeTo(_M);
        _Corr.ResizeTo(_N);

        _V.ResizeTo(_N, _N);
        _V0.ResizeTo(_N, _N);

        _V_updated.ResizeTo(_N, _N);

        _D.ResizeTo(_M, _N);
        _D_T.ResizeTo(_N, _M);

        _Aux.ResizeTo(_M, _M);

        ////////////////////////////////////////////////////////////////////////////////////////////////////////

        _constraints = new TF1*[_M];

        _constraints[0] = new TF1("Total energy conservation", this, &kin_fits::ene_consv_pi, 0, 1, _N);
        _constraints[1] = new TF1("Total x-mom conservation", this, &kin_fits::mom_x_consv_pi, 0, 1, _N);
        _constraints[2] = new TF1("Total y-mom conservation", this, &kin_fits::mom_y_consv_pi, 0, 1, _N);
        _constraints[3] = new TF1("Total z-mom conservation", this, &kin_fits::mom_z_consv_pi, 0, 1, _N);

        _constraints[4] = new TF1("Cluster time 1 conservation", this, &kin_fits::cluster_time_cons_first, 0, 1, _N);
        _constraints[5] = new TF1("Cluster time 2 conservation", this, &kin_fits::cluster_time_cons_second, 0, 1, _N);
        _constraints[6] = new TF1("Cluster time 3 conservation", this, &kin_fits::cluster_time_cons_third, 0, 1, _N);
        _constraints[7] = new TF1("Cluster time 4 conservation", this, &kin_fits::cluster_time_cons_fourth, 0, 1, _N);

        _constraints[8] = new TF1("Neutral inv mass conservation", this, &kin_fits::minv_neu_cons, 0, 1, _N);
        _constraints[9] = new TF1("Charged inv mass conservation", this, &kin_fits::minv_ch_cons, 0, 1, _N);

        /////////////////////////////////////////////////////////////////////////////////////////////////////////

      };

      void FillConstructors(Double_t *par, Double_t *sigmas);

      void solution();

      TVectorD FinalPars() { return _P; };
      TMatrixD FinalCovMtx() { return _V; };

      Double_t Chi2Final() { return _chi2; };

      ////////////////////////////////////////////////////////////////////////////////////////////////////
      
      void cons_vect();
      void cons_diff();

      Double_t ene_consv_pi(const Double_t *xx, const Double_t *pp);
      Double_t mom_x_consv_pi(const Double_t *xx, const Double_t *pp);
      Double_t mom_y_consv_pi(const Double_t *xx, const Double_t *pp);
      Double_t mom_z_consv_pi(const Double_t *xx, const Double_t *pp);

      Double_t cluster_time_cons_first(const Double_t *xx, const Double_t *pp);
      Double_t cluster_time_cons_second(const Double_t *xx, const Double_t *pp);
      Double_t cluster_time_cons_third(const Double_t *xx, const Double_t *pp);
      Double_t cluster_time_cons_fourth(const Double_t *xx, const Double_t *pp);

      Double_t minv_ch_cons(const Double_t *xx, const Double_t *pp);
      Double_t minv_neu_cons(const Double_t *xx, const Double_t *pp);

      Double_t tri_kaon_path(const Double_t *xx, const Double_t *pp);

      ////////////////////////////////////////////////////////////////////////////////////////////////////

    private:
      UInt_t _N, _M; //Num of parameters, constraints
      TVectorD _P, _P0, _P1, _C, _Lambda, _Corr;
      TMatrixD _V, _V0, _D, _D_T, _V_updated, _Aux;

      TF1 **_constraints;

      Double_t _det, _chi2;
      Int_t _iter_num, _fail;
  };

}