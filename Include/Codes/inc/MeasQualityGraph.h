#ifndef MEASQUALITYGRAPH_H
#define MEASQUALITYGRAPH_H

#include <algorithm>

#include <TGraph.h>
#include <TMultiGraph.h>
#include <TGaxis.h>
#include <TPad.h>
#include <TCanvas.h>
#include <TLegend.h>

#include <kloe_class.h>

namespace KLOE
{
    class MeasQualityGraph : protected pm00
    {
    private:
        std::vector<TGraph *>
            _effGraphs;

        std::vector<Double_t>
            _x;

        std::vector<std::vector<Double_t>>
            _y;

        TGraph
            *_ReErrGraph,
            *_ImErrGraph;

        TMultiGraph
            *_multiGraph;

        TLegend 
            *_legend;

        TString
            _mode,
            _xTitle,
            _yTitle,
            _yRightTitle;

        TGaxis
            *_axis;

        TCanvas
            *_canva;

        TPad
            *_pad,
            *_overlay;

        Int_t
            _NumOfPts;

        Double_t
            _minLimit,
            _maxLimit,
            _minYLimitL,
            _maxYLimitL,
            _minYLimitR,
            _maxYLimitR;

        void SetGraphStyle();

        // Modes:
        // Efficiency
        // FitResultErr

    public:
        MeasQualityGraph(TString mode, Int_t NumOfPts, std::vector<Double_t> &x, std::vector<std::vector<Double_t>> &y, TString xTitle, TString yTitle, TString yRightTitle, Double_t *legendPos);

        void DrawGraphs(TString imageName);
    };
}

#endif // !MEASQUALITYGRAPH_H