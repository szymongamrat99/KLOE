// File with methods from KinFit class for KLOE
// Date: 09/02/2025
// Author: Szymon Gamrat

#include <KinFitter.h>

using namespace KLOE;

KinFitter::KinFitter(std::string mode, Int_t N_free, Int_t N_const, Int_t M, Int_t M_active, Int_t loopcount, Double_t chisqrstep, ErrorHandling::ErrorLogs &logger) : _N_free(N_free), _N_const(N_const), _M(M), _M_act(M_active), _loopcount(loopcount), _jmin(0), _jmax(0), _CHISQRSTEP(chisqrstep), _logger(logger), _mode(mode), _V(N_free + N_const, N_free + N_const), _V_T(N_free + N_const, N_free + N_const), _V_init(N_free + N_const, N_free + N_const), _V_invert(N_free + N_const, N_free + N_const), _V_final(N_free + N_const, N_free + N_const), _V_aux(N_free + N_const, N_free + N_const), _D(M, N_free + N_const), _D_T(N_free + N_const, M), _Aux(M, M), _C(M), _L(M), _CORR(N_free + N_const), _X(N_free + N_const), _X_init(N_free + N_const), _X_final(N_free + N_const), _X_init_aux(N_free + N_const), _C_aux(M), _L_aux(M)
{
  if (_mode == "Omega")
    _baseObj = new ConstraintsOmega();
  else if (_mode == "SignalGlobal")
    _baseObj = new ConstraintsSignal();
  else if (_mode == "Trilateration")
    _baseObj = new ConstraintsTrilateration();
  else if (_mode == "Test")
    _baseObj = new ConstraintsTest();

  _D_real.ResizeTo(M, N_free);
  _D_T_real.ResizeTo(N_free, M);

  _CORR_real.ResizeTo(N_free);
  _Aux_real.ResizeTo(M, M);
};

KinFitter::KinFitter(std::string mode, Int_t N_free, Int_t N_const, Int_t M, Int_t M_active, Int_t loopcount, Int_t jmin, Int_t jmax, Double_t chisqrstep, ErrorHandling::ErrorLogs &logger) : _N_free(N_free), _N_const(N_const), _M(M), _M_act(M_active), _CHISQRSTEP(chisqrstep), _loopcount(loopcount), _jmin(jmin), _jmax(jmax), _logger(logger), _mode(mode), _V(N_free + N_const, N_free + N_const), _V_T(N_free + N_const, N_free + N_const), _V_init(N_free + N_const, N_free + N_const), _V_invert(N_free + N_const, N_free + N_const), _V_final(N_free + N_const, N_free + N_const), _V_aux(N_free + N_const, N_free + N_const), _D(M, N_free + N_const), _D_T(N_free + N_const, M), _Aux(M, M), _C(M), _L(M), _CORR(N_free + N_const), _X(N_free + N_const), _X_init(N_free + N_const), _X_final(N_free + N_const), _X_init_aux(N_free + N_const), _C_aux(M), _L_aux(M)
{
  if (_mode == "Omega")
    _baseObj = new ConstraintsOmega();
  else if (_mode == "SignalGlobal")
    _baseObj = new ConstraintsSignal();
  else if (_mode == "Trilateration")
    _baseObj = new ConstraintsTrilateration();
  else if (_mode == "Test")
    _baseObj = new ConstraintsTest();

  _D_real.ResizeTo(M, N_free);
  _D_T_real.ResizeTo(N_free, M);

  _Aux_real.ResizeTo(M, M);
  _CORR_real.ResizeTo(N_free);
};

Int_t KinFitter::ParameterInitialization(Float_t *Params, Float_t *Errors)
{
  TMatrixD SubV(_N_free, _N_free);

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
      if (_mode == "SignalGlobal")
      {
        if (_X(4) < 0)
          _X(4) = MIN_CLU_ENE;
        if (_X(9) < 0)
          _X(9) = MIN_CLU_ENE;
        if (_X(14) < 0)
          _X(14) = MIN_CLU_ENE;
        if (_X(19) < 0)
          _X(19) = MIN_CLU_ENE;

        // if(_X(3) < 0)
        //   _X(3) = 0.;
        // if(_X(8) < 0)
        //   _X(8) = 0.;
        // if(_X(13) < 0)
        //   _X(13) = 0.;
        // if(_X(18) < 0)
        //   _X(18) = 0.;

        _baseObj->SetParameters(_X.GetMatrixArray());
        _baseObj->IntermediateReconstruction();
      }

      for (Int_t l = 0; l < _M; l++)
      {
        _C(l) = _constraints[l]->EvalPar(0, _X.GetMatrixArray());
        for (Int_t m = 0; m < _N_free + _N_const; m++)
        {
          _constraints[l]->SetParameters(_X.GetMatrixArray());
          if (m < _N_free)
          {
            _D(l, m) = _constraints[l]->GradientPar(m, 0, 0.01 * sqrt(_V_init(m, m)));
          }
          else
            _D(l, m) = 0;
        }
      }

      _D_T = _D_T.Transpose(_D);

      _Aux = (_D * _V * _D_T);

      _Aux = _Aux.Invert(&_det);

      if (_det == 0)
        throw ErrorHandling::ErrorCodes::DET_ZERO;
      else if (TMath::IsNaN(_det))
        throw ErrorHandling::ErrorCodes::NAN_VAL;

      _L = (_Aux * _C);

      _CORR = _V * _D_T * _L;

      _X = _X - _CORR;

      _V_final = _V - _V * _D_T * _Aux * _D * _V;

      _CHISQR = Dot((_X - _X_init), _V_invert * (_X - _X_init));

      if (_CHISQR - _CHISQRTMP < _CHISQRSTEP && _CHISQR - _CHISQRTMP > 0)
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

  _baseObj->SetParameters(_X.GetMatrixArray());
  _baseObj->IntermediateReconstruction();

  _V = _V_final;

  return _CHISQRTMP;
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

    _constraints.push_back(new TF1(ConstSet[i].c_str(), _baseObj, constraintMap[ConstSet[i]], 0, 1, _N_free + _N_const));
  }

  return 0;
}

void KinFitter::GetResults(TVectorD &X, TMatrixD &V, TVectorD &X_init, TMatrixD &V_init, TVectorD &C, TVectorD &L)
{
  X = _X;
  V = _V;
  X_init = _X_init;
  V_init = _V_init;
  C = _C_aux;
  L = _L_aux;
}

void KinFitter::GetResults(TVectorD &X, TMatrixD &V, TVectorD &X_init, TMatrixD &V_init, std::vector<Float_t> trkFit[2], std::vector<Float_t> &KchrecFit, std::vector<Float_t> &KchboostFit, std::vector<Float_t> &ipFit, std::vector<Float_t> photonFit[4], std::vector<Float_t> &KnerecFit, std::vector<Float_t> &KnereclorFit)
{
  X = _X;
  V = _V;
  X_init = _X_init;
  V_init = _V_init;

  for (Int_t i = 0; i < 2; i++)
  {
    trkFit[i] = _baseObj->pionCh[i].fourMom;
  }

  KchrecFit = _baseObj->Kchrec.total;
  KchboostFit = _baseObj->Kchboost.total;

  ipFit = _baseObj->ip;

  for (Int_t i = 0; i < 4; i++)
    photonFit[i] = _baseObj->photon[i].total;

  KnerecFit = _baseObj->Knerec.total;
  KnereclorFit = _baseObj->Knereclor.total;
}

void KinFitter::GetResults(TVectorD &X, TMatrixD &V, TVectorD &X_init, TMatrixD &V_init, std::vector<Float_t> trkFit[2], std::vector<Float_t> &OmegaFit, std::vector<Float_t> &ipFit, std::vector<Float_t> photonFit[4], std::vector<Float_t> &Pi0OmegaFit, std::vector<Float_t> &PhiMomFit)
{
  X = _X;
  V = _V;
  X_init = _X_init;
  V_init = _V_init;

  for (Int_t i = 0; i < 2; i++)
  {
    trkFit[i] = _baseObj->pionCh[i].fourMom;
  }

  OmegaFit = _baseObj->omega.total;

  ipFit = _baseObj->ip;

  for (Int_t i = 0; i < 4; i++)
    photonFit[i] = _baseObj->photon[i].total;

  Pi0OmegaFit =  _baseObj->pionNe[0].total;
  PhiMomFit = _baseObj->phi.total;
}

void KinFitter::GetResults(TVectorD &X, TMatrixD &V, TVectorD &X_init, TMatrixD &V_init)
{
  X = _X;
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

    while (abs(derivativeAux - derivative) > 0.001 * abs(derivative) || L < 1)
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

// void KinFitter::PhotonPairing(std::vector<NeuPart> _Photons)
// {
//   Int_t PhotonsNum = _Photons.size();
//   if (PhotonsNum % 2 != 0)
//   {
//     std::cerr << "Liczba fotonów musi być parzysta." << std::endl;
//     return;
//   }
//   const Int_t div = Int_t(PhotonsNum / 2.);
//   Double_t AuxMinv[div];
//   Double_t AuxChi2 = 0., AuxChi2Min = 99999.;
//   do
//   {
//     for (Int_t i = 0; i < div; i++)
//     {
//       AuxMinv[i] += pow(_Photons[i * div].FourMom[3] + _Photons[i * div + 1].FourMom[3], 2);
//       for (Int_t j = 0; j < 3; j++)
//       {
//         AuxMinv[i] -= pow(_Photons[i * div].FourMom[j] + _Photons[i * div + 1].FourMom[j], 2);
//       }
//       AuxMinv[i] = sqrt(AuxMinv[i]);
//       AuxChi2 += pow(AuxMinv[i] - mPi0, 2);
//     }
//     AuxChi2 = sqrt(AuxChi2);
//     if (AuxChi2 < AuxChi2Min)
//     {
//       AuxChi2Min = AuxChi2;
//     }
//   } while (std::next_permutation(_Photons.begin(), _Photons.end(), [](const auto & lhs, const auto & rhs)
//                                  { return lhs.index < rhs.index; }));
// }