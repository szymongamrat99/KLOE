#ifndef FRACTION_FITTER_H
#define FRACTION_FITTER_H

#include <iostream>
#include <vector>

#include <TMath.h>
#include <TH1D.h>
#include <TFractionFitter.h>

#include <const.h>

class FractionFitter
{
private:
  TFractionFitter
      *_fit;
  std::vector<TH1 *>
      _MC;
  TH1
      *_data,
      *_result;
  TObjArray
      *_obj;

  TFitResultPtr
      _resultPtr;

  Double_t
      _fitStats;

  UInt_t
      _MCNum;

  Int_t
      _status;

public:
  FractionFitter(std::vector<TH1 *> MC, TH1 *data);

  void SetRange(Int_t binMin, Int_t binMax);
  void ReleaseRange();
  void Constrain(std::vector<Int_t> lower, std::vector<Int_t> upper);
  TH1 *GetResultHisto();
  void GetResults(std::vector<Double_t> &norm, std::vector<Double_t> &error, Double_t &chi2);

  void FitFractions();
};

#endif // !FRACTION_FITTER_H