#ifndef DOUBLE_RATIO_H
#define DOUBLE_RATIO_H

#include <vector>
#include <TMath.h>

namespace KLOE
{
  class DoubleRatio
  {
  private:
    std::vector<Int_t> // 0 - pm; 1 - 00
                    _NeventsKL,
                    _NeventsKS;

    std::vector<Double_t> // 0 - pm; 1 - 00
                      _EffKL,
                      _EffKS,
                      _ErrEffKL,
                      _ErrEffKS;
    Double_t
          _tKch,
          _tKne;

    void GetNevents();
    
  public:
    DoubleRatio();
    DoubleRatio(std::vector<Int_t> NeventsKL, std::vector<Int_t> NeventsKS);

    void SetTimeKch(Double_t tKch);
    void SetTimeKne(Double_t tKne);
  };
  

} // namespace KLOE


#endif