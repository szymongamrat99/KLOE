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
#include <ConstraintsTest.h>

namespace KLOE
{
  class KinFitter
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
        _V_real,
        _D,
        _D_T,
        _D_real,
        _D_T_real,
        _V_T,
        _V_final,
        _V_aux,
        _V_min,
        _Aux,
        _Aux_real,
        _V_invert,
        _V_init;

    TVectorD
        _X,
        _X_real,
        _C,
        _X_final,
        _X_final_real,
        _L,
        _CORR,
        _CORR_real,
        _X_init,
        _X_min,
        _C_min,
        _L_min,
        _C_aux,
        _L_aux,
        _X_init_min,
        _X_init_aux;



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
            {"minvconsvneutralkaon", &KinFit::MinvConsvNeuKaon},
            {"minvconsvchargedkaon", &KinFit::MinvConsvChKaon},
            {"minvconsvomega", &KinFit::MinvConsvOmega},
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

    ErrorHandling::ErrorCodes _err_code = ErrorHandling::ErrorCodes::NO_ERROR;

    std::vector<TF1 *> _constraints;
    Double_t _value_min;

    std::vector<Int_t> _chosen;

    std::string
        _mode;

    Double_t AdjustCyclicalVar(Double_t angleCorrected, Double_t angleOriginal);

  public:
    KinFitter(std::string mode, Int_t N_free, Int_t N_const, Int_t M, Int_t M_active, Int_t loopcount, Double_t chisqrstep, ErrorHandling::ErrorLogs &logger);

    KinFitter(std::string mode, Int_t N_free, Int_t N_const, Int_t M, Int_t M_active, Int_t loopcount, Int_t jmin, Int_t jmax, Double_t chisqrstep, ErrorHandling::ErrorLogs &logger);

    Int_t ParameterInitialization(Float_t *Params, Float_t *Errors);
    Int_t ParameterInitialization(Double_t *Params, Double_t *Errors);

    Int_t ConstraintSet(std::vector<std::string> ConstrSet);
    Int_t ConstraintSet(std::vector<TF1 *> ConstrSet);

    Double_t FitFunction(Double_t bunchCorr = 0);

    void GetResults(TVectorD &X, TMatrixD &V, TVectorD &X_init, TMatrixD &V_init, std::vector<Float_t> &ipFit, std::vector<Float_t> photonFit[4], std::vector<Float_t> &KnerecFit, std::vector<Float_t> &phiFit);
    void GetResults(TVectorD &X, TMatrixD &V, TVectorD &X_init, TMatrixD &V_init, std::vector<Float_t> trkFit[2], std::vector<Float_t> &KchrecFit, std::vector<Float_t> &KchboostFit, std::vector<Float_t> &ipFit, std::vector<Float_t> photonFit[4], std::vector<Float_t> &KnerecFit, std::vector<Float_t> &KnereclorFit);
    void GetResults(TVectorD &X, TMatrixD &V, TVectorD &X_init, TMatrixD &V_init, std::vector<Float_t> trkFit[2], std::vector<Float_t> &OmegaFit, std::vector<Float_t> &ipFit, std::vector<Float_t> photonFit[4], std::vector<Float_t> Pi0OmegaFit[2], std::vector<Float_t> &PhiMomFit);
    void GetResults(TVectorD &X, TMatrixD &V, TVectorD &X_init, TMatrixD &V_init);


    Double_t EnergyCalc(Double_t *p, Double_t mass);

    Double_t DerivativeCalc(Int_t i, Int_t j);

    Double_t EnergyCalc(TLorentzVector p, Double_t mass);

    ErrorHandling::ErrorCodes GetErrorCode() { return _err_code; }
  };

}

#endif
