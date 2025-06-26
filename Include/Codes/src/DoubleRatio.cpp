#include <const.h>
#include <DoubleRatio.h>

namespace KLOE
{
    DoubleRatio::DoubleRatio()
    {
        _NeventsKL.assign(2,0);
        _NeventsKS.assign(2,0);
    }

    DoubleRatio::DoubleRatio(std::vector<Int_t> NeventsKL, std::vector<Int_t> NeventsKS) : _NeventsKL(NeventsKL), _NeventsKS(NeventsKS)
    {
    }

    void DoubleRatio::SetTimeKch(Double_t tKch)
    {
        _tKch = tKch;
    }

    void DoubleRatio::SetTimeKne(Double_t tKne)
    {
        _tKne = tKne;
    }

    void DoubleRatio::GetNevents()
    {
        if(_tKch < _tKne)
        {
            _NeventsKS[0]++;
            _NeventsKL[1]++;
        }
        else
        {
            _NeventsKS[1]++;
            _NeventsKL[0]++;
        }
    }
}