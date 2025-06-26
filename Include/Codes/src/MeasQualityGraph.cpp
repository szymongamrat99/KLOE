#include <MeasQualityGraph.h>

namespace KLOE
{
  MeasQualityGraph::MeasQualityGraph(TString mode, Int_t NumOfPts, std::vector<Double_t> &x, std::vector<std::vector<Double_t>> &y, TString xTitle, TString yTitle, TString yRightTitle, Double_t *legendPos) : _mode(mode), _NumOfPts(NumOfPts), _x(x), _y(y), _xTitle(xTitle), _yTitle(yTitle), _yRightTitle(yRightTitle)
  {
    _multiGraph = new TMultiGraph();
    _legend = new TLegend(legendPos[0], legendPos[2], legendPos[1], legendPos[3]);

    if (ToLower(_mode) == "efficiency" && _y.size() == channNum)
    {
      for (Int_t i = 0; i < channNum; i++)
      {
        _effGraphs.push_back(new TGraph(_NumOfPts, _x.data(), _y[i].data()));
        _effGraphs[i]->SetMarkerColor(channColor[i]);
        _multiGraph->Add(_effGraphs[i]);
      }

      _minYLimitL = 0.0;
      _maxYLimitL = 1.0;
      _minYLimitR = 0.0;
      _maxYLimitR = 1.0;
    }
    else if (ToLower(_mode) == "fitresulterr" && _y.size() == 2)
    {
      _ReErrGraph = new TGraph(_NumOfPts, _x.data(), _y[0].data());
      _ImErrGraph = new TGraph(_NumOfPts, _x.data(), _y[1].data());

      _ReErrGraph->SetMarkerStyle(20);
      _ImErrGraph->SetMarkerStyle(22);

      _ReErrGraph->SetMarkerSize(2);
      _ImErrGraph->SetMarkerSize(2);

      _multiGraph->Add(_ReErrGraph);
      _multiGraph->Add(_ImErrGraph);

      _legend->AddEntry(_ReErrGraph, "Re(#varepsilon'/#varepsilon) relative error", "p");
      _legend->AddEntry(_ImErrGraph, "Im(#varepsilon'/#varepsilon) relative error", "p");

      _minYLimitL = 0.0;
      _maxYLimitL = 1.0;
      _minYLimitR = 0.0;
      _maxYLimitR = 1.0;
    }
    else if (ToLower(_mode) == "errvsvalue" && _y.size() == 2)
    {
      _ReErrGraph = new TGraph(_NumOfPts, _x.data(), _y[0].data());
      _ImErrGraph = new TGraph(_NumOfPts, _x.data(), _y[1].data());

      _ReErrGraph->SetMarkerStyle(20);
      _ImErrGraph->SetMarkerStyle(22);

      _ReErrGraph->SetMarkerSize(2);
      _ImErrGraph->SetMarkerSize(2);

      _legend->AddEntry(_ReErrGraph, "Parameter value", "p");
      _legend->AddEntry(_ImErrGraph, "Parameter error", "p");

      _minYLimitL = *std::min_element(_y[0].begin(), _y[0].end());
      _maxYLimitL = *std::max_element(_y[0].begin(), _y[0].end());
      _minYLimitR = *std::min_element(_y[1].begin(), _y[1].end());
      _maxYLimitR = *std::max_element(_y[1].begin(), _y[1].end());
    }
    else
    {
      std::cout << "Wrong Parameters!" << std::endl;
    }

    _minLimit = *std::min_element(_x.begin(), _x.end());
    _maxLimit = *std::max_element(_x.begin(), _x.end());

    _canva = new TCanvas();
    _canva->SetRightMargin(0.15);
    _canva->Range(_minLimit - 0.1 * _minLimit, 0, _maxLimit + 0.1 * _maxLimit, 1);

    SetGraphStyle();
  }

  void MeasQualityGraph::SetGraphStyle()
  {
    TAxis *auxAxis[2];

    if (ToLower(_mode) == "fitresulterr")
    {
      auxAxis[0] = _multiGraph->GetXaxis();
      auxAxis[1] = _multiGraph->GetYaxis();
    }
    if (ToLower(_mode) == "errvsvalue")
    {
      auxAxis[0] = _ReErrGraph->GetXaxis();
      auxAxis[1] = _ImErrGraph->GetYaxis();
    }

    auxAxis[0]->SetTitle(_xTitle);
    auxAxis[1]->SetTitle(_yTitle);

    auxAxis[0]->SetLimits(_minLimit - 0.1 * _minLimit, _maxLimit + 0.1 * _maxLimit);

    auxAxis[1]->SetRangeUser(_minYLimitL - 0.1 * _minYLimitL, _maxYLimitL + 0.1 * _maxLimit);

    for (Int_t i = 0; i < 2; i++)
    {
      auxAxis[i]->CenterTitle(1);
      auxAxis[i]->SetMaxDigits(3);
    }

    _legend->SetFillColor(kWhite);
  }

  void MeasQualityGraph::DrawGraphs(TString imageName)
  {
    _canva->cd();

    if (ToLower(_mode) == "fitresulterr")
    {
     _multiGraph->Draw("AP");
    }
    if (ToLower(_mode) == "errvsvalue")
    {
      _ReErrGraph->Draw("AP");
      _ImErrGraph->Draw("SAME");

      _axis = new TGaxis(gPad->GetUxmax(), gPad->GetUymin(), gPad->GetUxmax(), gPad->GetUymax(), _minYLimitR - 0.1 * _minYLimitR, _maxYLimitR + 0.1 * _maxYLimitR, 110, "+L");

      _axis->SetTitle(_yRightTitle);
      _axis->CenterTitle(1);
      _axis->SetVertical(1);
      _axis->SetTitleOffset(1.5);
      _axis->SetTitleFont(62);
      _axis->SetTitleSize(0.05);

      _axis->Draw();
    }

    _legend->Draw();

    _canva->Print(imageName);
  }
}