// File with methods from KinFit class for KLOE
// Date: 09/02/2025
// Author: Szymon Gamrat

#include <KinFitter.h>
#include <TROOT.h>
#include <TDecompSVD.h>

#include <chrono>

using namespace KLOE;

KinFitter::KinFitter(std::string mode, Int_t N_free, Int_t N_const, Int_t M, Int_t M_active, Int_t loopcount, Double_t chisqrstep, ErrorHandling::ErrorLogs &logger) : _N_free(N_free), _N_const(N_const), _M(M), _M_act(M_active), _loopcount(loopcount), _jmin(0), _jmax(0), _CHISQRSTEP(chisqrstep), _logger(logger), _mode(mode), _V(N_free + N_const, N_free + N_const), _V_T(N_free + N_const, N_free + N_const), _V_init(N_free + N_const, N_free + N_const), _V_invert(N_free + N_const, N_free + N_const), _V_final(N_free + N_const, N_free + N_const), _V_aux(N_free + N_const, N_free + N_const), _D(M, N_free + N_const), _D_T(N_free + N_const, M), _Aux(M, M), _C(M), _L(M), _CORR(N_free + N_const), _X(N_free + N_const), _X_init(N_free + N_const), _X_final(N_free + N_const), _X_init_aux(N_free + N_const), _C_aux(M), _L_aux(M)
{
  if (_mode == "Omega")
    _objOmega = std::make_unique<ConstraintsOmega>(logger);
  else if (_mode == "SignalGlobal")
    _objSignal = std::make_unique<ConstraintsSignal>(logger);
  else if (_mode == "Trilateration")
    _objTrilateration = std::make_unique<ConstraintsTrilateration>(logger);
  else if (_mode == "Test")
    _baseObj = std::make_unique<ConstraintsTest>(logger);

  _D_real.ResizeTo(M, N_free);
  _D_T_real.ResizeTo(N_free, M);

  _CORR_real.ResizeTo(N_free);
  _Aux_real.ResizeTo(M, M);

  _Aux.Zero();
  _D.Zero();
  _D_T.Zero();
  _C.Zero();
  _L.Zero();
};

KinFitter::KinFitter(std::string mode, Int_t N_free, Int_t N_const, Int_t M, Int_t M_active, Int_t loopcount, Int_t jmin, Int_t jmax, Double_t chisqrstep, ErrorHandling::ErrorLogs &logger) : _N_free(N_free), _N_const(N_const), _M(M), _M_act(M_active), _CHISQRSTEP(chisqrstep), _loopcount(loopcount), _jmin(jmin), _jmax(jmax), _logger(logger), _mode(mode), _V(N_free + N_const, N_free + N_const), _V_T(N_free + N_const, N_free + N_const), _V_init(N_free + N_const, N_free + N_const), _V_invert(N_free + N_const, N_free + N_const), _V_final(N_free + N_const, N_free + N_const), _V_aux(N_free + N_const, N_free + N_const), _D(M, N_free + N_const), _D_T(N_free + N_const, M), _Aux(M, M), _C(M), _L(M), _CORR(N_free + N_const), _X(N_free + N_const), _X_init(N_free + N_const), _X_final(N_free + N_const), _X_init_aux(N_free + N_const), _C_aux(M), _L_aux(M)
{

  if (_mode == "Omega")
    _objOmega = std::make_unique<ConstraintsOmega>(logger);
  else if (_mode == "SignalGlobal")
    _objSignal = std::make_unique<ConstraintsSignal>(logger);
  else if (_mode == "Trilateration")
    _objTrilateration = std::make_unique<ConstraintsTrilateration>(logger);
  else if (_mode == "Test")
    _baseObj = std::make_unique<ConstraintsTest>(logger);

  _D_real.ResizeTo(M, N_free);
  _D_T_real.ResizeTo(N_free, M);

  _Aux_real.ResizeTo(M, M);
  _CORR_real.ResizeTo(N_free);
};

KinFitter::~KinFitter()
{
  // No need to manually delete unique_ptr objects
}

Int_t KinFitter::ParameterInitialization(Float_t *Params, Float_t *Errors)
{
  TMatrixD SubV(_N_free, _N_free);
  SubV.Zero();

  for (Int_t i = 0; i < _N_free + _N_const; i++)
  {
    _X_init(i) = Params[i];

    // Under assumption that there are no cross-correlations between variables
    _V_init(i, i) = pow(Errors[i], 2);

    if (std::isnan(_V_init(i, i)))
      _err_flag = true;
  }

  _V_init.GetSub(0, _N_free - 1, 0, _N_free - 1, SubV);

  SubV.Invert();

  _X = _X_init;
  _V = _V_init;
  _V_invert.SetSub(0, 0, SubV);

  return 0;
};

Int_t KinFitter::ParameterInitialization(Double_t *Params, Double_t *Errors)
{
  TMatrixD SubV(_N_free, _N_free);
  SubV.Zero();

  for (Int_t i = 0; i < _N_free + _N_const; i++)
  {
    _X_init(i) = Params[i];

    // Under assumption that there are now cross-correlations between variables
    _V_init(i, i) = pow(Errors[i], 2);

    if (std::isnan(_V_init(i, i)))
      _err_flag = true;
  }

  _V_init.GetSub(0, _N_free - 1, 0, _N_free - 1, SubV);

  SubV.Invert();

  _X = _X_init;
  _V = _V_init;
  _V_invert.SetSub(0, 0, SubV);

  return 0;
};

Double_t KinFitter::FitFunction(Double_t bunchCorr)
{
  pm00 objAux;

  _CHISQR = 999999.;
  _CHISQRTMP = 999999.;
  _Aux.Zero();

  Int_t iter = -1;

  // Correction of cluster time - [ns]
  if (bunchCorr != 0 && _mode == "Trilateration")
  {
    for (Int_t i = 0; i < 4; i++)
      _X[i * 5 + 3] += bunchCorr;
  }
  // ------------------------------------------

  _err_flag = 0;

  _err_code = ErrorHandling::ErrorCodes::NO_ERROR;

  Double_t CHISQRTOT = 0.;

  for (Int_t i = 0; i < _loopcount; i++)
  {
    try
    {
      // Enforce positive energies
      if (_X(4) < 0)
        _X(4) = MIN_CLU_ENE;
      if (_X(9) < 0)
        _X(9) = MIN_CLU_ENE;
      if (_X(14) < 0)
        _X(14) = MIN_CLU_ENE;
      if (_X(19) < 0)
        _X(19) = MIN_CLU_ENE;

      for (Int_t l = 0; l < _M; l++)
      {
        _constraints[l]->SetParameters(_X.GetMatrixArray());
        _C(l) = _constraints[l]->EvalPar(0, _X.GetMatrixArray());
        for (Int_t m = 0; m < _N_free + _N_const; m++)
        {
          Double_t auxVal = 0.;

          if (m < _N_free)
          {
            _constraints[l]->SetParError(m, sqrt(_V_init(m, m)));
            auxVal = _constraints[l]->GradientPar(m, 0);

            _D(l, m) = auxVal;
          }
          else
            _D(l, m) = 0;

          std::cout << _D(l, m) << " ";
        }
      }

      std::cout << std::endl;

      _D_T.Transpose(_D);

      _Aux = (_D * _V * _D_T);

      TDecompSVD svd(_Aux);
      svd.SetTol(1e-17); // Ustawienie tolerancji dla wartości singularnych, aby uniknąć problemów z macierzą osobliwą
      Bool_t svdStatus = false;
      Double_t condition = svd.Condition();

      if (condition > 1e17)
      {
        throw ErrorHandling::ErrorCodes::DET_ZERO;
      }

      _Aux = svd.Invert(svdStatus);

      if (!svdStatus)
        throw ErrorHandling::ErrorCodes::DET_ZERO;

      _L = (_Aux * _C);

      _CORR = _V * _D_T * _L;

      _X = _X - _CORR;

      _V_final = _V - _V * _D_T * _Aux * _D * _V;

      TVectorD deltaX = _X - _X_init;

      _CHISQR = _V_invert.Similarity(deltaX);

      iter++;

      Double_t maxDeltaX = 0;
      for (int j = 0; j < _CORR.GetNrows(); j++)
      {
        double relChange = abs(_CORR(j) / (_X(j) + 1e-9)); // 1e-9 zapobiega /0
        if (relChange > maxDeltaX)
          maxDeltaX = relChange;
      }

      if (abs(_CHISQR - _CHISQRTMP) < _CHISQRSTEP || maxDeltaX < 1e-5)
        break;

      _L_aux = _L;
      _C_aux = _C;
      _FUNVALTMP = _FUNVAL;
      _CHISQRTMP = _CHISQR;

      _Aux.Zero();
      _C.Zero();
      _L.Zero();
      _D.Zero();
      _D_T.Zero();
      _CORR.Zero();
    }
    catch (ErrorHandling::ErrorCodes err)
    {
      // _err_code = err;
      // _logger.getErrLog(err, "iteration no. " + std::to_string(i));
      break;
    }
  }

  if (_mode == "SignalGlobal")
  {
    _objSignal->IntermediateReconstruction(_X.GetMatrixArray());
  }
  else if (_mode == "Trilateration")
  {
    _objTrilateration->IntermediateReconstruction(_X.GetMatrixArray());
  }
  else if (_mode == "Omega")
  {
    _objOmega->IntermediateReconstruction(_X.GetMatrixArray());
  }
  else
  {
    _baseObj->SetParameters(_X.GetMatrixArray());
    _baseObj->IntermediateReconstruction();
  }

  if (iter != -1)
  {
    _CHISQR = Dot((_X - _X_init), _V_invert * (_X - _X_init));
    _V = _V_final;
  }

  return _CHISQR;
};

Double_t KinFitter::EnergyCalc(Double_t *p, Double_t mass)
{
  Double_t energy = 0.;

  for (Int_t i = 0; i < 3; i++)
    energy += pow(p[i], 2);

  energy += pow(mass, 2);
  energy = sqrt(energy);

  return energy;
};

Double_t KinFitter::EnergyCalc(TLorentzVector p, Double_t mass)
{
  Double_t energy = 0.;

  for (Int_t i = 0; i < 3; i++)
    energy += pow(p[i], 2);

  energy += pow(mass, 2);
  energy = sqrt(energy);

  return energy;
};

Int_t KinFitter::ConstraintSet(std::vector<std::string> ConstSet)
{
  for (Size_t i = 0; i < ConstSet.size(); i++)
  {
    std::transform(ConstSet[i].begin(),
                   ConstSet[i].end(),
                   ConstSet[i].begin(),
                   ::tolower);

    if (_mode == "SignalGlobal")
      _constraints.push_back(std::make_unique<TF1>(ConstSet[i].c_str(), _objSignal.get(), constraintMapSignal[ConstSet[i]], 0, 1, _N_free + _N_const));
    else if (_mode == "Trilateration")
      _constraints.push_back(std::make_unique<TF1>(ConstSet[i].c_str(), _objTrilateration.get(), constraintMapTrilateration[ConstSet[i]], 0, 1, _N_free + _N_const));
    else if (_mode == "Omega")
      _constraints.push_back(std::make_unique<TF1>(ConstSet[i].c_str(), _objOmega.get(), constraintMapOmega[ConstSet[i]], 0, 1, _N_free + _N_const));
    else
      _constraints.push_back(std::make_unique<TF1>(ConstSet[i].c_str(), _baseObj.get(), constraintMap[ConstSet[i]], 0, 1, _N_free + _N_const));
  }

  return 0;
}

void KinFitter::GetResults(TVectorD &X, TMatrixD &V, TVectorD &X_init, TMatrixD &V_init, std::vector<Double_t> &ipFit, std::vector<Double_t> photonFit[4], std::vector<Double_t> &KnerecFit, std::vector<Double_t> &phiFit)
{
  X = _X;
  V = _V;
  X_init = _X_init;
  V_init = _V_init;

  ipFit = _objTrilateration->fip;

  for (Int_t i = 0; i < 4; i++)
    photonFit[i] = _objTrilateration->fphoton[i].total;

  KnerecFit = _objTrilateration->fKnerec.total;
  phiFit = _objTrilateration->fphi.total;
}

void KinFitter::GetResults(TVectorD &X, TMatrixD &V, TVectorD &X_init, TMatrixD &V_init, std::vector<Double_t> trkFit[2], std::vector<Double_t> &KchrecFit, std::vector<Double_t> &KchboostFit, std::vector<Double_t> &ipFit, std::vector<Double_t> photonFit[4], std::vector<Double_t> &KnerecFit, std::vector<Double_t> &KnereclorFit)
{
  X = _X;
  V = _V;
  X_init = _X_init;
  V_init = _V_init;

  for (Int_t i = 0; i < 2; i++)
  {
    trkFit[i] = _objSignal->fpionCh[i].fourMom;
  }

  KchrecFit = _objSignal->fKchrec.total;
  KchboostFit = _objSignal->fKchboost.total;

  ipFit = _objSignal->fip;

  for (Int_t i = 0; i < 4; i++)
    photonFit[i] = _objSignal->fphoton[i].total;

  KnerecFit = _objSignal->fKnerec.total;
  KnereclorFit = _objSignal->fKnereclor.total;
}

void KinFitter::GetResults(TVectorD &X, TMatrixD &V, TVectorD &X_init, TMatrixD &V_init, std::vector<Double_t> trkFit[2], std::vector<Double_t> &OmegaFit, std::vector<Double_t> &ipFit, std::vector<Double_t> photonFit[4], std::vector<Double_t> Pi0OmegaFit[2], std::vector<Double_t> &PhiMomFit)
{
  X = _X;
  V = _V;
  X_init = _X_init;
  V_init = _V_init;

  for (Int_t i = 0; i < 2; i++)
  {
    if (trkFit[i].size() != 4)
      trkFit[i].resize(4);
    trkFit[i] = _objOmega->fpionCh[i].fourMom;
  }

  OmegaFit = _objOmega->fomega.total;

  ipFit = _objOmega->fip;

  if (photonFit[0].size() != 8)
    photonFit[0].resize(8);
  if (photonFit[1].size() != 8)
    photonFit[1].resize(8);
  for (Int_t i = 0; i < 4; i++)
    photonFit[i] = _objOmega->fphoton[i].total;

  for (Int_t i = 0; i < 2; i++)
  {
    if (Pi0OmegaFit[i].size() != 5)
      Pi0OmegaFit[i].resize(5);
    Pi0OmegaFit[i] = _objOmega->fpionNe[i].total;
  }

  PhiMomFit = _objOmega->fphi.total;
}

void KinFitter::GetResults(TVectorD &X, TMatrixD &V, TVectorD &X_init, TMatrixD &V_init)
{
  for (Int_t i = 0; i < _X.GetNrows(); i++)
    X[i] = _X[i];

  // X = _X;
  V = _V;
  X_init = _X_init;
  V_init = _V_init;
}

Double_t KinFitter::AdjustCyclicalVar(Double_t angleCorrected, Double_t angleOriginal)
{
  Double_t diff = angleCorrected - angleOriginal;

  Long64_t n = round(diff / M_PI);

  Double_t angleAdjusted = angleCorrected + n * M_PI;

  return angleAdjusted;
}

Double_t KinFitter::DerivativeCalc(Int_t i, Int_t j)
{
  Double_t derivative = 0., pplus = 0., pminus = 0.;
  Double_t step = 0.;

  Double_t derivativeAux = 0.;

  TVectorD X_aux_plus = _X,
           X_aux_minus = _X;
  TMatrixD V_aux = _V;

  if (V_aux(j, j) == 0)
  {
    derivative = 999.;
  }
  else
  {
    step = 0.01 * sqrt(V_aux(j, j));

    Int_t L = 0;

    while (abs(derivativeAux - derivative) > 0.01 * abs(derivative) || L < 1)
    {
      derivativeAux = derivative;

      X_aux_plus[j] = _X[j] + step;
      X_aux_minus[j] = _X[j] - step;

      Double_t f_plus = _constraints[i]->EvalPar(0, X_aux_plus.GetMatrixArray());
      Double_t f_minus = _constraints[i]->EvalPar(0, X_aux_minus.GetMatrixArray());

      derivative = (f_plus - f_minus) / (2 * step);

      if (derivative == 0)
        break;

      step = step / 2.;
      L++;
    }
  }

  return derivativeAux;
};

Bool_t KinFitter::CheckMatrixNan(TMatrixD &matrix)
{
  for (Int_t i = 0; i < matrix.GetNrows(); i++)
  {
    for (Int_t j = 0; j < matrix.GetNcols(); j++)
    {
      if (std::isnan(matrix(i, j)) || std::isinf(matrix(i, j)))
        return true;
    }
  }
  return false;
}
