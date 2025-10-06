#include <charged_mom.h>

namespace KLOE
{
  template <typename F, typename T>
  ChargedVtxRec<F, T>::ChargedVtxRec(Int_t &nv, Int_t &ntv, T *ivOld, F *IP, F *CurV, F *PhiV, F *CotV, F *xvOld, F *yvOld, F *zvOld, Int_t &mode) : _iv(ivOld), _nv(nv), _ntv(ntv), _IP(IP), _CurV(CurV), _PhiV(PhiV), _CotV(CotV), _xv(xvOld), _yv(yvOld), _zv(zvOld), _mode(mode){};

  template <typename F, typename T>
  ChargedVtxRec<F, T>::ChargedVtxRec() : _iv(nullptr), _nv(def), _ntv(def), _IP(nullptr), _CurV(nullptr), _PhiV(nullptr), _CotV(nullptr), _xv(nullptr), _yv(nullptr), _zv(nullptr), _mode(def){};

  template <typename F, typename T>
  void ChargedVtxRec<F, T>::charged_mom(F CurvOld, F PhivOld, F CotvOld, F *mom_vec, Int_t mode)
  {
    mom_vec[0] = cos(PhivOld) * 1000. / abs(CurvOld);
    mom_vec[1] = sin(PhivOld) * 1000. / abs(CurvOld);
    mom_vec[2] = CotvOld * 1000. / abs(CurvOld);

    if (mode == 1)
      mom_vec[3] = sqrt(pow(mom_vec[0], 2) + pow(mom_vec[1], 2) + pow(mom_vec[2], 2) + pow(mPiCh, 2));
    else if (mode == 2)
      mom_vec[3] = sqrt(pow(mom_vec[0], 2) + pow(mom_vec[1], 2) + pow(mom_vec[2], 2) + pow(mElec, 2));
    else if (mode == 3)
      mom_vec[3] = sqrt(pow(mom_vec[0], 2) + pow(mom_vec[1], 2) + pow(mom_vec[2], 2) + pow(mMuon, 2));
  }

  template <typename F, typename T>
  void ChargedVtxRec<F, T>::charged_mom(Int_t &i, F *mom_vec)
  {
    mom_vec[0] = cos(_PhiV[i]) * 1000. / abs(_CurV[i]);
    mom_vec[1] = sin(_PhiV[i]) * 1000. / abs(_CurV[i]);
    mom_vec[2] = _CotV[i] * 1000. / abs(_CurV[i]);

    if (_mode == 1)
      mom_vec[3] = sqrt(pow(mom_vec[0], 2) + pow(mom_vec[1], 2) + pow(mom_vec[2], 2) + pow(mPiCh, 2));
    else if (_mode == 2)
      mom_vec[3] = sqrt(pow(mom_vec[0], 2) + pow(mom_vec[1], 2) + pow(mom_vec[2], 2) + pow(mElec, 2));
    else if (_mode == 3)
      mom_vec[3] = sqrt(pow(mom_vec[0], 2) + pow(mom_vec[1], 2) + pow(mom_vec[2], 2) + pow(mMuon, 2));
  }

  template <typename F, typename T>
  ErrorHandling::ErrorCodes ChargedVtxRec<F, T>::findKchRec(T mcflag, F *KchRec, F *trk1, F *trk2, Int_t *vtaken, ErrorHandling::ErrorLogs &logger)
  {
    F mom_vec1Tmp[4], mom_vec2Tmp[4], KchTmp[9];
    std::vector<Int_t> ivTmp(_iv, _iv + MaxNumTrkV);
    std::map<Int_t, Int_t> mapTmp = pm00::CountRepeatingElements(ivTmp);
    pm00::Clear1DArray(3, vtaken);
    bool found = false;

    // Initialization of momentum smearing
    // -------------------------------------------------------------

    const Int_t numberOfMomenta = 2;

    TVectorT<Double_t>
        momVecMC(numberOfMomenta * 3),
        momVecSmeared(numberOfMomenta * 3);

    std::vector<double> elems = properties["momSmearing"]["covarianceMatrix"]["fElements"].get<std::vector<Double_t>>();

    Int_t nRows = properties["momSmearing"]["covarianceMatrix"]["fNrows"],
          nCols = properties["momSmearing"]["covarianceMatrix"]["fNcols"];

    TMatrixT<Double_t>
        covMatrix(nRows, nCols, elems.data());

    KLOE::MomentumSmearing<Double_t> CovMatrixCalcObj(momVecMC, covMatrix);
    // -------------------------------------------------------------

    if (_nv == 0 || _ntv == 0)
    {
      ErrorHandling::ErrorCodes err = ErrorHandling::ErrorCodes::NOT_ENOUGH_CHARGED_TRACKS;
      return err;
    }

    for (Int_t i = 0; i < _nv; i++)
    {
      if (mapTmp[_iv[i]] == 2)
      {
        for (Int_t j1 = 0; j1 < _ntv - 1; j1++)
          for (Int_t j2 = j1 + 1; j2 < _ntv; j2++)
          {
            if (_iv[j1] - 1 == i && _iv[j2] - 1 == i)
            {
              if (pm00::signum(_CurV[j1]) != pm00::signum(_CurV[j2]))
              {
                ChargedVtxRec::charged_mom(j1, mom_vec1Tmp);
                ChargedVtxRec::charged_mom(j2, mom_vec2Tmp);

                Float_t KchrecSmeared[9], KchboostSmeared[9], energyPion[2];

                if (mcflag == 1)
                {
                  momVecMC[0] = mom_vec1Tmp[0];
                  momVecMC[1] = mom_vec1Tmp[1];
                  momVecMC[2] = mom_vec1Tmp[2];
                  momVecMC[3] = mom_vec2Tmp[0];
                  momVecMC[4] = mom_vec2Tmp[1];
                  momVecMC[5] = mom_vec2Tmp[2];

                  CovMatrixCalcObj.SetMCVector(momVecMC);
                  CovMatrixCalcObj.SmearMomentum();
                  CovMatrixCalcObj.GetSmearedMomentum(momVecSmeared);
                }
                else
                {
                  momVecSmeared[0] = mom_vec1Tmp[0];
                  momVecSmeared[1] = mom_vec1Tmp[1];
                  momVecSmeared[2] = mom_vec1Tmp[2];
                  momVecSmeared[3] = mom_vec2Tmp[0];
                  momVecSmeared[4] = mom_vec2Tmp[1];
                  momVecSmeared[5] = mom_vec2Tmp[2];
                }

                energyPion[0] = sqrt(pow(momVecSmeared[0], 2) + pow(momVecSmeared[1], 2) + pow(momVecSmeared[2], 2) + pow(mPiCh, 2));
                energyPion[1] = sqrt(pow(momVecSmeared[3], 2) + pow(momVecSmeared[4], 2) + pow(momVecSmeared[5], 2) + pow(mPiCh, 2));

                for (Int_t k = 0; k < 3; k++)
                  KchTmp[k] = momVecSmeared[k] + momVecSmeared[k + 3];

                KchTmp[3] = energyPion[0] + energyPion[1];

                KchTmp[4] = pow(KchTmp[0], 2) + pow(KchTmp[1], 2) + pow(KchTmp[2], 2);
                KchTmp[5] = sqrt(pow(KchTmp[3], 2) - KchTmp[4]);
                KchTmp[4] = sqrt(KchTmp[4]);

                if (vtaken[0] <= -1 || abs(KchTmp[5] - mK0) < abs(KchRec[5] - mK0))
                {
                  vtaken[0] = i;
                  vtaken[1] = j1;
                  vtaken[2] = j2;
                  for (Int_t k = 0; k < 6; k++)
                    KchRec[k] = KchTmp[k];
                  KchRec[6] = _xv[vtaken[0]];
                  KchRec[7] = _yv[vtaken[0]];
                  KchRec[8] = _zv[vtaken[0]];
                  for (Int_t k = 0; k < 4; k++)
                  {
                    trk1[k] = mom_vec1Tmp[k];
                    trk2[k] = mom_vec2Tmp[k];
                  }
                  found = true;
                }
              }
            }
          }
      }
    }
    if (!found)
    {
      ErrorHandling::ErrorCodes err = ErrorHandling::ErrorCodes::NOT_ENOUGH_CHARGED_TRACKS;
      return err;
    }
    return ErrorHandling::ErrorCodes::NO_ERROR;
  }

  template <typename F, typename T>
  ErrorHandling::ErrorCodes ChargedVtxRec<F, T>::findKchRec(T mcflag, Bool_t smearingFlag, TMatrixT<Double_t> covMatrix, std::vector<F> &KchRec, std::vector<F> &trk1, std::vector<F> &trk2, std::vector<Int_t> &vtaken, ErrorHandling::ErrorLogs &logger)
  {
    F mom_vec1Tmp[4], mom_vec2Tmp[4], KchTmp[9];
    std::vector<Int_t> ivTmp(_iv, _iv + MaxNumTrkV);
    std::map<Int_t, Int_t> mapTmp = pm00::CountRepeatingElements(ivTmp);
    vtaken.clear();
    vtaken.resize(3);
    bool found = false;

    // Initialization of momentum smearing
    // -------------------------------------------------------------

    const Int_t numberOfMomenta = 2;

    TVectorT<Double_t>
        momVecMC(numberOfMomenta * 3),
        momVecSmeared(numberOfMomenta * 3);

    KLOE::MomentumSmearing<Double_t> CovMatrixCalcObj(momVecMC, covMatrix);
    // -------------------------------------------------------------

    if (_nv == 0 || _ntv == 0)
    {
      ErrorHandling::ErrorCodes err = ErrorHandling::ErrorCodes::NOT_ENOUGH_CHARGED_TRACKS;
      return err;
    }

    for (Int_t i = 0; i < _nv; i++)
    {
      if (mapTmp[_iv[i]] == 2)
      {
        for (Int_t j1 = 0; j1 < _ntv - 1; j1++)
          for (Int_t j2 = j1 + 1; j2 < _ntv; j2++)
          {
            if (_iv[j1] - 1 == i && _iv[j2] - 1 == i)
            {
              if (pm00::signum(_CurV[j1]) != pm00::signum(_CurV[j2]))
              {
                ChargedVtxRec::charged_mom(j1, mom_vec1Tmp);
                ChargedVtxRec::charged_mom(j2, mom_vec2Tmp);

                Float_t KchrecSmeared[9], KchboostSmeared[9], energyPion[2];

                if (mcflag == 1 && smearingFlag == 1)
                {
                  momVecMC[0] = mom_vec1Tmp[0];
                  momVecMC[1] = mom_vec1Tmp[1];
                  momVecMC[2] = mom_vec1Tmp[2];
                  momVecMC[3] = mom_vec2Tmp[0];
                  momVecMC[4] = mom_vec2Tmp[1];
                  momVecMC[5] = mom_vec2Tmp[2];

                  CovMatrixCalcObj.SetMCVector(momVecMC);
                  CovMatrixCalcObj.SmearMomentum();
                  CovMatrixCalcObj.GetSmearedMomentum(momVecSmeared);
                }
                else
                {
                  momVecSmeared[0] = mom_vec1Tmp[0];
                  momVecSmeared[1] = mom_vec1Tmp[1];
                  momVecSmeared[2] = mom_vec1Tmp[2];
                  momVecSmeared[3] = mom_vec2Tmp[0];
                  momVecSmeared[4] = mom_vec2Tmp[1];
                  momVecSmeared[5] = mom_vec2Tmp[2];
                }

                energyPion[0] = sqrt(pow(momVecSmeared[0], 2) + pow(momVecSmeared[1], 2) + pow(momVecSmeared[2], 2) + pow(mPiCh, 2));
                energyPion[1] = sqrt(pow(momVecSmeared[3], 2) + pow(momVecSmeared[4], 2) + pow(momVecSmeared[5], 2) + pow(mPiCh, 2));

                for (Int_t k = 0; k < 3; k++)
                  KchTmp[k] = momVecSmeared[k] + momVecSmeared[k + 3];

                KchTmp[3] = energyPion[0] + energyPion[1];

                KchTmp[4] = pow(KchTmp[0], 2) + pow(KchTmp[1], 2) + pow(KchTmp[2], 2);
                KchTmp[5] = sqrt(pow(KchTmp[3], 2) - KchTmp[4]);
                KchTmp[4] = sqrt(KchTmp[4]);
                if (vtaken[0] <= -1 || abs(KchTmp[5] - mK0) < abs(KchRec[5] - mK0))
                {
                  vtaken[0] = i;
                  vtaken[1] = j1;
                  vtaken[2] = j2;
                  for (Int_t k = 0; k < 6; k++)
                    KchRec[k] = KchTmp[k];
                  KchRec[6] = _xv[vtaken[0]];
                  KchRec[7] = _yv[vtaken[0]];
                  KchRec[8] = _zv[vtaken[0]];

                  for (Int_t k = 0; k < 3; k++)
                  {
                    trk1[k] = momVecSmeared[k];
                    trk2[k] = momVecSmeared[k + 3];
                  }

                  trk1[3] = sqrt(pow(trk1[0], 2) + pow(trk1[1], 2) + pow(trk1[2], 2) + pow(mPiCh, 2));
                  trk2[3] = sqrt(pow(trk2[0], 2) + pow(trk2[1], 2) + pow(trk2[2], 2) + pow(mPiCh, 2));

                  found = true;
                }
              }
            }
          }
      }
    }
    if (!found)
    {
      ErrorHandling::ErrorCodes err = ErrorHandling::ErrorCodes::NOT_ENOUGH_CHARGED_TRACKS;
      return err;
    }
    return ErrorHandling::ErrorCodes::NO_ERROR;
  }

  template <typename F, typename T>
  ErrorHandling::ErrorCodes ChargedVtxRec<F, T>::findKSLRec(Int_t kaonFlag, Int_t KSvtx, F *KchRec, F *trk1, F *trk2, Int_t *vtaken, ErrorHandling::ErrorLogs &logger)
  {
    F mom_vec1Tmp[4], mom_vec2Tmp[4], KchTmp[9], cyl_vol[2];
    std::vector<Int_t> ivTmp(_iv, _iv + MaxNumTrkV);
    std::map<Int_t, Int_t> mapTmp = pm00::CountRepeatingElements(ivTmp);
    pm00::Clear1DArray(3, vtaken);
    vtaken[0] = -1;
    _dist = 999.;
    bool found = false;
    Bool_t distFlag = false;
    for (Int_t i = 0; i < _nv; i++)
    {
      cyl_vol[0] = sqrt(pow(_xv[i] - _IP[0], 2) + pow(_yv[i] - _IP[1], 2));
      cyl_vol[1] = abs(_zv[i] - _IP[2]);
      if (kaonFlag == 16)
      {
        distFlag = (cyl_vol[0] < 10. && cyl_vol[1] < 20.);
      }
      else if (kaonFlag == 10)
      {
        distFlag = (KSvtx != i);
      }
      else
      {
        ErrorHandling::ErrorCodes err = ErrorHandling::ErrorCodes::NOT_RECOGNIZED;

        return err;
      }
      if (mapTmp[_iv[i]] == 2 && distFlag)
      {
        for (Int_t j1 = 0; j1 < _ntv - 1; j1++)
          for (Int_t j2 = j1 + 1; j2 < _ntv; j2++)
          {
            if (_iv[j1] - 1 == i && _iv[j2] - 1 == i)
            {
              if (pm00::signum(_CurV[j1]) != pm00::signum(_CurV[j2]))
              {
                ChargedVtxRec::charged_mom(j1, mom_vec1Tmp);
                ChargedVtxRec::charged_mom(j2, mom_vec2Tmp);
                for (Int_t k = 0; k < 4; k++)
                  KchTmp[k] = mom_vec1Tmp[k] + mom_vec2Tmp[k];
                KchTmp[4] = pow(KchTmp[0], 2) + pow(KchTmp[1], 2) + pow(KchTmp[2], 2);
                KchTmp[5] = sqrt(pow(KchTmp[3], 2) - KchTmp[4]);
                KchTmp[4] = sqrt(KchTmp[4]);
                if (vtaken[0] <= -1 || abs(KchTmp[5] - mK0) < abs(KchRec[5] - mK0))
                {
                  vtaken[0] = i;
                  vtaken[1] = j1;
                  vtaken[2] = j2;
                  for (Int_t k = 0; k < 6; k++)
                    KchRec[k] = KchTmp[k];
                  KchRec[6] = _xv[vtaken[0]];
                  KchRec[7] = _yv[vtaken[0]];
                  KchRec[8] = _zv[vtaken[0]];
                  for (Int_t k = 0; k < 4; k++)
                  {
                    trk1[k] = mom_vec1Tmp[k];
                    trk2[k] = mom_vec2Tmp[k];
                  }
                  found = true;
                }
              }
            }
          }
      }
    }
    if (!found)
    {
      ErrorHandling::ErrorCodes err = ErrorHandling::ErrorCodes::NOT_ENOUGH_CHARGED_TRACKS;
      return err;
    }
    return ErrorHandling::ErrorCodes::NO_ERROR;
  }

  template <typename F, typename T>
  ErrorHandling::ErrorCodes ChargedVtxRec<F, T>::findKSLRec(Int_t kaonFlag, Int_t KSvtx, std::vector<F> &KchRec, std::vector<F> &trk1, std::vector<F> &trk2, std::vector<Int_t> &vtaken, ErrorHandling::ErrorLogs &logger)
  {
    F mom_vec1Tmp[4], mom_vec2Tmp[4], KchTmp[9], cyl_vol[2];
    std::vector<Int_t> ivTmp(_iv, _iv + MaxNumTrkV);
    std::map<Int_t, Int_t> mapTmp = pm00::CountRepeatingElements(ivTmp);
    vtaken.clear();
    vtaken.resize(3);
    vtaken[0] = -1;
    _dist = 999.;
    bool found = false;
    Bool_t distFlag = false;
    for (Int_t i = 0; i < _nv; i++)
    {
      cyl_vol[0] = sqrt(pow(_xv[i] - _IP[0], 2) + pow(_yv[i] - _IP[1], 2));
      cyl_vol[1] = abs(_zv[i] - _IP[2]);
      if (kaonFlag == 16)
      {
        distFlag = (cyl_vol[0] < 10. && cyl_vol[1] < 20.);
      }
      else if (kaonFlag == 10)
      {
        distFlag = (KSvtx != i);
      }
      else
      {
        ErrorHandling::ErrorCodes err = ErrorHandling::ErrorCodes::NOT_RECOGNIZED;

        return err;
      }
      if (mapTmp[_iv[i]] == 2 && distFlag)
      {
        for (Int_t j1 = 0; j1 < _ntv - 1; j1++)
          for (Int_t j2 = j1 + 1; j2 < _ntv; j2++)
          {
            if (_iv[j1] - 1 == i && _iv[j2] - 1 == i)
            {
              if (pm00::signum(_CurV[j1]) != pm00::signum(_CurV[j2]))
              {
                ChargedVtxRec::charged_mom(j1, mom_vec1Tmp);
                ChargedVtxRec::charged_mom(j2, mom_vec2Tmp);
                for (Int_t k = 0; k < 4; k++)
                  KchTmp[k] = mom_vec1Tmp[k] + mom_vec2Tmp[k];
                KchTmp[4] = pow(KchTmp[0], 2) + pow(KchTmp[1], 2) + pow(KchTmp[2], 2);
                KchTmp[5] = sqrt(pow(KchTmp[3], 2) - KchTmp[4]);
                KchTmp[4] = sqrt(KchTmp[4]);
                if (vtaken[0] <= -1 || abs(KchTmp[5] - mK0) < abs(KchRec[5] - mK0))
                {
                  vtaken[0] = i;
                  vtaken[1] = j1;
                  vtaken[2] = j2;
                  for (Int_t k = 0; k < 6; k++)
                    KchRec[k] = KchTmp[k];
                  KchRec[6] = _xv[vtaken[0]];
                  KchRec[7] = _yv[vtaken[0]];
                  KchRec[8] = _zv[vtaken[0]];
                  for (Int_t k = 0; k < 4; k++)
                  {
                    trk1[k] = mom_vec1Tmp[k];
                    trk2[k] = mom_vec2Tmp[k];
                  }
                  found = true;
                }
              }
            }
          }
      }
    }
    if (!found)
    {
      ErrorHandling::ErrorCodes err = ErrorHandling::ErrorCodes::NOT_ENOUGH_CHARGED_TRACKS;
      return err;
    }
    return ErrorHandling::ErrorCodes::NO_ERROR;
  }

  template <typename F, typename T>
  ErrorHandling::ErrorCodes ChargedVtxRec<F, T>::findKClosestRec(F *KchRec, F *trk1, F *trk2, Int_t *vtaken, ErrorHandling::ErrorLogs &logger)
  {
    F mom_vec1Tmp[4], mom_vec2Tmp[4], KchTmp[9], cyl_vol[2], distCyl = 999.;
    std::vector<Int_t> ivTmp(_iv, _iv + MaxNumTrkV);
    std::map<Int_t, Int_t> mapTmp = pm00::CountRepeatingElements(ivTmp);
    pm00::Clear1DArray(3, vtaken);
    bool found = false;
    for (Int_t i = 0; i < _nv; i++)
    {
      cyl_vol[0] = sqrt(pow(_xv[i] - _IP[0], 2) + pow(_yv[i] - _IP[1], 2) + pow(_zv[i] - _IP[2], 2));
      cyl_vol[1] = abs(_zv[i] - _IP[2]);
      if (mapTmp[_iv[i]] == 2)
      {
        for (Int_t j1 = 0; j1 < _ntv - 1; j1++)
          for (Int_t j2 = j1 + 1; j2 < _ntv; j2++)
          {
            if (_iv[j1] - 1 == i && _iv[j2] - 1 == i)
            {
              if (pm00::signum(_CurV[j1]) != pm00::signum(_CurV[j2]))
              {
                ChargedVtxRec::charged_mom(j1, mom_vec1Tmp);
                ChargedVtxRec::charged_mom(j2, mom_vec2Tmp);
                for (Int_t k = 0; k < 4; k++)
                  KchTmp[k] = mom_vec1Tmp[k] + mom_vec2Tmp[k];
                KchTmp[4] = pow(KchTmp[0], 2) + pow(KchTmp[1], 2) + pow(KchTmp[2], 2);
                KchTmp[5] = sqrt(pow(KchTmp[3], 2) - KchTmp[4]);
                KchTmp[4] = sqrt(KchTmp[4]);
                if (vtaken[0] <= -1 || cyl_vol[0] < distCyl)
                {
                  distCyl = cyl_vol[0];
                  vtaken[0] = i;
                  vtaken[1] = j1;
                  vtaken[2] = j2;
                  for (Int_t k = 0; k < 6; k++)
                    KchRec[k] = KchTmp[k];
                  KchRec[6] = _xv[vtaken[0]];
                  KchRec[7] = _yv[vtaken[0]];
                  KchRec[8] = _zv[vtaken[0]];
                  for (Int_t k = 0; k < 4; k++)
                  {
                    trk1[k] = mom_vec1Tmp[k];
                    trk2[k] = mom_vec2Tmp[k];
                  }
                  found = true;
                }
              }
            }
          }
      }
    }
    if (!found)
    {
      ErrorHandling::ErrorCodes err = ErrorHandling::ErrorCodes::NOT_ENOUGH_CHARGED_TRACKS;
      return err;
    }
    return ErrorHandling::ErrorCodes::NO_ERROR;
  }

  template <typename F, typename T>
  ErrorHandling::ErrorCodes ChargedVtxRec<F, T>::findKClosestRec(std::vector<F> &KchRec, std::vector<F> &trk1, std::vector<F> &trk2, std::vector<Int_t> &vtaken, ErrorHandling::ErrorLogs &logger)
  {
    F mom_vec1Tmp[4], mom_vec2Tmp[4], KchTmp[9], cyl_vol[2], distCyl = 999.;
    std::vector<Int_t> ivTmp(_iv, _iv + MaxNumTrkV);
    std::map<Int_t, Int_t> mapTmp = pm00::CountRepeatingElements(ivTmp);
    vtaken.clear();
    vtaken.resize(3);
    bool found = false;
    for (Int_t i = 0; i < _nv; i++)
    {
      cyl_vol[0] = sqrt(pow(_xv[i] - _IP[0], 2) + pow(_yv[i] - _IP[1], 2) + pow(_zv[i] - _IP[2], 2));
      cyl_vol[1] = abs(_zv[i] - _IP[2]);
      if (mapTmp[_iv[i]] == 2)
      {
        for (Int_t j1 = 0; j1 < _ntv - 1; j1++)
          for (Int_t j2 = j1 + 1; j2 < _ntv; j2++)
          {
            if (_iv[j1] - 1 == i && _iv[j2] - 1 == i)
            {
              if (pm00::signum(_CurV[j1]) != pm00::signum(_CurV[j2]))
              {
                ChargedVtxRec::charged_mom(j1, mom_vec1Tmp);
                ChargedVtxRec::charged_mom(j2, mom_vec2Tmp);
                for (Int_t k = 0; k < 4; k++)
                  KchTmp[k] = mom_vec1Tmp[k] + mom_vec2Tmp[k];
                KchTmp[4] = pow(KchTmp[0], 2) + pow(KchTmp[1], 2) + pow(KchTmp[2], 2);
                KchTmp[5] = sqrt(pow(KchTmp[3], 2) - KchTmp[4]);
                KchTmp[4] = sqrt(KchTmp[4]);
                if (vtaken[0] <= -1 || cyl_vol[0] < distCyl)
                {
                  distCyl = cyl_vol[0];
                  vtaken[0] = i;
                  vtaken[1] = j1;
                  vtaken[2] = j2;
                  for (Int_t k = 0; k < 6; k++)
                    KchRec[k] = KchTmp[k];
                  KchRec[6] = _xv[vtaken[0]];
                  KchRec[7] = _yv[vtaken[0]];
                  KchRec[8] = _zv[vtaken[0]];
                  for (Int_t k = 0; k < 4; k++)
                  {
                    trk1[k] = mom_vec1Tmp[k];
                    trk2[k] = mom_vec2Tmp[k];
                  }
                  found = true;
                }
              }
            }
          }
      }
    }
    if (!found)
    {
      ErrorHandling::ErrorCodes err = ErrorHandling::ErrorCodes::NOT_ENOUGH_CHARGED_TRACKS;
      return err;
    }
    return ErrorHandling::ErrorCodes::NO_ERROR;
  }

  template <typename F, typename T>
  Int_t ChargedVtxRec<F, T>::KaonMomFromBoost(F *pKaon, F *pboost, F *pKaonBoost) const
  {
    std::string
        name = "";
    name = base_path + logs_dir + "KchFromBoost_" + pm00::getCurrentDate() + ".log";

    ErrorHandling::ErrorLogs logger(name);

    F
        pK_from_boost = 0.,
        pb_mod = 0.,
        pK_mod = 0.,
        dot = 0.,
        cosb_sq = 0.,
        beta_sq = 0.,
        beta_gamma_sq = 0.,
        gamma_sq = 0.,
        pcm_sq = 0.,
        C = 0.,
        A = 0.,
        B = 0.;

    for (Int_t i = 0; i < 3; i++)
    {
      pb_mod += pow(pboost[i], 2);
      pK_mod += pow(pKaon[i], 2);
      dot += pKaon[i] * pboost[i];
    }

    cosb_sq = pow(dot, 2) / (pb_mod * pK_mod);

    pb_mod = sqrt(pb_mod);
    pK_mod = sqrt(pK_mod);

    beta_sq = pow(pb_mod / pboost[3], 2);
    beta_gamma_sq = pow(pb_mod, 2) / (pow(pboost[3], 2) - pow(pb_mod, 2));
    gamma_sq = 1. + beta_gamma_sq;

    pcm_sq = 0.25 * (pow(pboost[3], 2) - pow(pb_mod, 2)) - pow(mK0, 2);

    C = pow(beta_gamma_sq * pow(mK0, 2) - pcm_sq, 2);
    A = pow(gamma_sq * (1. - beta_sq * cosb_sq), 2);
    B = gamma_sq * ((1. + beta_sq * cosb_sq) * (beta_gamma_sq * pow(mK0, 2) - pcm_sq) - 2. * beta_gamma_sq * pow(mK0, 2) * cosb_sq);

    F
        disc = pow(B, 2) - A * C;

    try
    {
      // Check, if this is a quadratic equation
      if (A == 0)
        throw ErrorHandling::ErrorCodes::DENOM_EQ_ZERO;

      // Check, if discriminant is good
      if (disc < 0. && disc > -100.)
        disc = 0.;
      else if (disc < -100.)
        throw ErrorHandling::ErrorCodes::DELTA_LT_ZERO;
    }
    catch (ErrorHandling::ErrorCodes err)
    {

      return Int_t(err);
    }

    disc = sqrt(disc);

    F
        P1 = (-B + disc) / A,
        P2 = (-B - disc) / A;

    if (P1 > 0. && P2 < 0.)
      pK_from_boost = sqrt(P1);
    else if (P2 > 0. && P1 < 0.)
      pK_from_boost = sqrt(P2);
    else if (P1 > 0. && P2 > 0.)
    {
      if (dot < 0.)
        pK_from_boost = sqrt(std::min(P1, P2));
      else
        pK_from_boost = sqrt(std::max(P1, P2));
    }

    for (Int_t i = 0; i < 3; i++)
      pKaonBoost[i] = pKaon[i] * (pK_from_boost / pK_mod);

    pKaonBoost[3] = sqrt(pow(pK_from_boost, 2) + pow(mK0, 2));
    pKaonBoost[4] = pK_from_boost;
    pKaonBoost[5] = mK0;
    pKaonBoost[6] = pKaon[6];
    pKaonBoost[7] = pKaon[7];
    pKaonBoost[8] = pKaon[8];

    return 0;
  }
  template <typename F, typename T>
  Int_t ChargedVtxRec<F, T>::KaonMomFromBoost(std::vector<F> &pKaon, F *pboost, std::vector<F> &pKaonBoost) const
  {
    std::string
        name = "";
    name = base_path + logs_dir + "KchFromBoost_" + pm00::getCurrentDate() + ".log";

    ErrorHandling::ErrorLogs logger(name);

    F
        pK_from_boost = 0.,
        pb_mod = 0.,
        pK_mod = 0.,
        dot = 0.,
        cosb_sq = 0.,
        beta_sq = 0.,
        beta_gamma_sq = 0.,
        gamma_sq = 0.,
        pcm_sq = 0.,
        C = 0.,
        A = 0.,
        B = 0.;

    for (Int_t i = 0; i < 3; i++)
    {
      pb_mod += pow(pboost[i], 2);
      pK_mod += pow(pKaon[i], 2);
      dot += pKaon[i] * pboost[i];
    }

    cosb_sq = pow(dot, 2) / (pb_mod * pK_mod);

    pb_mod = sqrt(pb_mod);
    pK_mod = sqrt(pK_mod);

    beta_sq = pow(pb_mod / pboost[3], 2);
    beta_gamma_sq = pow(pb_mod, 2) / (pow(pboost[3], 2) - pow(pb_mod, 2));
    gamma_sq = 1. + beta_gamma_sq;

    pcm_sq = 0.25 * (pow(pboost[3], 2) - pow(pb_mod, 2)) - pow(mK0, 2);

    C = pow(beta_gamma_sq * pow(mK0, 2) - pcm_sq, 2);
    A = pow(gamma_sq * (1. - beta_sq * cosb_sq), 2);
    B = gamma_sq * ((1. + beta_sq * cosb_sq) * (beta_gamma_sq * pow(mK0, 2) - pcm_sq) - 2. * beta_gamma_sq * pow(mK0, 2) * cosb_sq);

    F
        disc = pow(B, 2) - A * C;

    try
    {
      // Check, if this is a quadratic equation
      if (A == 0)
        throw ErrorHandling::ErrorCodes::DENOM_EQ_ZERO;

      // Check, if discriminant is good
      if (disc < 0. && disc > -100.)
        disc = 0.;
      else if (disc < -100.)
        throw ErrorHandling::ErrorCodes::DELTA_LT_ZERO;
    }
    catch (ErrorHandling::ErrorCodes err)
    {

      return Int_t(err);
    }

    disc = sqrt(disc);

    F
        P1 = (-B + disc) / A,
        P2 = (-B - disc) / A;

    if (P1 > 0. && P2 < 0.)
      pK_from_boost = sqrt(P1);
    else if (P2 > 0. && P1 < 0.)
      pK_from_boost = sqrt(P2);
    else if (P1 > 0. && P2 > 0.)
    {
      if (dot < 0.)
        pK_from_boost = sqrt(std::min(P1, P2));
      else
        pK_from_boost = sqrt(std::max(P1, P2));
    }

    for (Int_t i = 0; i < 3; i++)
      pKaonBoost[i] = pKaon[i] * (pK_from_boost / pK_mod);

    pKaonBoost[3] = sqrt(pow(pK_from_boost, 2) + pow(mK0, 2));
    pKaonBoost[4] = pK_from_boost;
    pKaonBoost[5] = mK0;
    pKaonBoost[6] = pKaon[6];
    pKaonBoost[7] = pKaon[7];
    pKaonBoost[8] = pKaon[8];

    return 0;
  }

  template <typename F, typename T>
  Int_t ChargedVtxRec<F, T>::IPBoostCorr(F *X_line, F *vec_line, F *X_plane, F *vec_plane, F *int_point) const
  {
    std::string
        name = "";
    name = base_path + logs_dir + "IPBoostCorrection_" + pm00::getCurrentDate() + ".log";

    ErrorHandling::ErrorLogs logger(name);

    F dot_prod_up = 0.,
      dot_prod_down = 0.;

    for (Int_t i = 0; i < 3; i++)
    {
      dot_prod_up += (X_line[i] - X_plane[i]) * vec_plane[i];
      dot_prod_down += vec_line[i] * vec_plane[i];
    }

    try
    {
      if (dot_prod_down == 0)
        throw ErrorHandling::ErrorCodes::DENOM_EQ_ZERO;

      for (Int_t i = 0; i < 3; i++)
      {
        int_point[i] = X_line[i] + (dot_prod_up / dot_prod_down) * vec_line[i];
      }

      return 0;
    }
    catch (ErrorHandling::ErrorCodes err)
    {
      return Int_t(err);
    }
  }

  template <typename F, typename T>
  Int_t ChargedVtxRec<F, T>::IPBoostCorr(F *X_line, F *vec_line, F *X_plane, F *vec_plane, std::vector<F> &int_point) const
  {
    std::string
        name = "";
    name = base_path + logs_dir + "IPBoostCorrection_" + pm00::getCurrentDate() + ".log";

    ErrorHandling::ErrorLogs logger(name);

    F dot_prod_up = 0.,
      dot_prod_down = 0.;

    for (Int_t i = 0; i < 3; i++)
    {
      dot_prod_up += (X_line[i] - X_plane[i]) * vec_plane[i];
      dot_prod_down += vec_line[i] * vec_plane[i];
    }

    try
    {
      if (dot_prod_down == 0)
        throw ErrorHandling::ErrorCodes::DENOM_EQ_ZERO;

      for (Int_t i = 0; i < 3; i++)
      {
        int_point[i] = X_line[i] + (dot_prod_up / dot_prod_down) * vec_line[i];
      }

      return 0;
    }
    catch (ErrorHandling::ErrorCodes err)
    {
      return Int_t(err);
    }
  }

  // Explicit definition of template class
  template class ChargedVtxRec<Float_t, Int_t>;
  template class ChargedVtxRec<Double_t, Int_t>;
  template class ChargedVtxRec<Float_t, UChar_t>;
  template class ChargedVtxRec<Double_t, UChar_t>;
}