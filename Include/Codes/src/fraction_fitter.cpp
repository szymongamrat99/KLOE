#include <fraction_fitter.h>

FractionFitter::FractionFitter(std::vector<TH1 *> MC, TH1 *data) : _MCNum(MC.size()), _MC(MC), _data(data)
{
    _obj = new TObjArray(_MCNum);

    for (Int_t i = 0; i < _MCNum; i++)
        _obj->Add(_MC[i]);

    _fit = new TFractionFitter(_data, _obj);
}

void FractionFitter::SetRange(Int_t binMin, Int_t binMax)
{
    _fit->SetRangeX(binMin, binMax);
}

void FractionFitter::ReleaseRange()
{
    _fit->ReleaseRangeX();
}

void FractionFitter::Constrain(std::vector<Int_t> lower, std::vector<Int_t> upper)
{
    for (Int_t i = 0; i < _MCNum; i++)
        _fit->Constrain(i, lower[i], upper[i]);
}

TH1* FractionFitter::GetResultHisto()
{
    return _result;
}

void FractionFitter::FitFractions()
{
    _status = _fit->Fit();
    _result = _fit->GetPlot();
}

void FractionFitter::GetResults(std::vector<Double_t> &norm, std::vector<Double_t> &error, Double_t &chi2)
{
    for (Int_t i = 0; i < _MCNum; i++)
        _fit->GetResult(i, norm[i], error[i]);

    chi2 = _fit->GetChisquare();
}