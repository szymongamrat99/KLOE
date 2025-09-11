#ifndef KINFITTER_H
#define KINFITTER_H

#include <vector> // for std::vector and other dynamic stuff
#include <algorithm>
#include <utility> // for std::swap

#include <TMath.h>
#include <TMatrixD.h>
#include <TVectorD.h>
#include <TError.h>
#include <TF1.h>
#include <TLorentzVector.h>

#include <const.h>

// Inclusion of base classes

#include <ConstraintsOmega.h>
#include <ConstraintsTrilateration.h>
#include <ConstraintsSignal.h>

namespace KLOE
{
  class KinFitter : public ConstraintsOmega, public ConstraintsTrilateration
  {
  private:
    Double_t
        _CHISQR,
        _FUNVAL,
        _CHISQRTMP,
        _FUNVALTMP,
        _CHISQRSTEP,
        _det;

    TMatrixD
        _V,
        _D,
        _D_T,
        _V_T,
        _V_final,
        _V_aux,
        _V_min,
        _Aux,
        _V_invert,
        _V_init;

    TVectorD
        _X,
        _C,
        _X_final,
        _L,
        _CORR,
        _X_init,
        _X_min,
        _C_min,
        _L_min,
        _C_aux,
        _L_aux,
        _X_init_min,
        _X_init_aux;

    std::vector<ChPart>
        _PiCh;

    std::vector<NeuPart>
        _Photon,
        _PiNeu,
        _Kaon;
    Phi
        _PhiMeson;

    KinFit *_baseObj;

    std::map<std::string, Double_t (KinFit::*)(Double_t *, Double_t *)>
        constraintMap = {
            {"energyconsvlab", &KinFit::EnergyConsvLAB},
            {"pxconsvlab", &KinFit::PxConsvLAB},
            {"pyconsvlab", &KinFit::PyConsvLAB},
            {"pzconsvlab", &KinFit::PzConsvLAB},
            {"photon1pathlab", &KinFit::Photon1PathConsvLAB},
            {"photon2pathlab", &KinFit::Photon2PathConsvLAB},
            {"photon3pathlab", &KinFit::Photon3PathConsvLAB},
            {"photon4pathlab", &KinFit::Photon4PathConsvLAB},
            {"energyconsvcm", &KinFit::EnergyConsvCM},
            {"minvconsv", &KinFit::MinvConsv},
            {"neutralxpathconsvlab", &KinFit::NeutralXPathConsvLAB},
            {"neutralypathconsvlab", &KinFit::NeutralYPathConsvLAB},
            {"neutralzpathconsvlab", &KinFit::NeutralZPathConsvLAB}};

  protected:
    Int_t
        _N_free,
        _N_const,
        _M,
        _M_act;
    Bool_t
        _err_flag,
        _fail;
    Int_t
        *_selected,
        _counter = 0,
        _jmin,
        _jmax,
        _loopcount,
        _N_clus;
    ErrorHandling::ErrorLogs &_logger;
    std::vector<TF1 *> _constraints;
    Double_t _value_min;

    std::vector<Int_t> _chosen;

    std::string
        _mode;

  public:
    KinFitter(std::string mode, Int_t N_free, Int_t N_const, Int_t M, Int_t M_active, Int_t loopcount, Double_t chisqrstep, ErrorHandling::ErrorLogs &logger);

    KinFitter(std::string mode, Int_t N_free, Int_t N_const, Int_t M, Int_t M_active, Int_t loopcount, Int_t jmin, Int_t jmax, Double_t chisqrstep, ErrorHandling::ErrorLogs &logger);

    Int_t ParameterInitialization(Float_t *Params, Float_t *Errors);
    Int_t ParameterInitialization(Double_t *Params, Double_t *Errors);

    Int_t ConstraintSet(std::vector<std::string> ConstrSet);
    Int_t ConstraintSet(std::vector<TF1 *> ConstrSet);

    void PhotonPairing(std::vector<NeuPart> _Photons);

    Double_t FitFunction(Double_t bunchCorr = 0);

    void GetResults(TVectorD &X, TMatrixD &V, TVectorD &X_init, TMatrixD &V_init, TVectorD &C, TVectorD &L);

    Double_t EnergyCalc(Double_t *p, Double_t mass);
    Double_t EnergyCalc(TLorentzVector p, Double_t mass);
  };

}

#endif
