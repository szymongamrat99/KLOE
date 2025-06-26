#ifndef EFFICIENCY_ANALYSIS_H
#define EFFICIENCY_ANALYSIS_H

#include <vector>
#include <fstream>

#include <TMath.h>
#include <TEfficiency.h>
#include <TTree.h>
#include <TH1.h>
#include <TGraphAsymmErrors.h>
#include <TF1.h>

#include <const.h>

namespace KLOE
{
  class EfficiencyAnalysis
  {
    private:
      TH1 
        *sigTotal;
      TTree
          *treeCuts;
      std::vector<Float_t>
                        *vars;
      std::ifstream 
                  *cuts_file;

      TEfficiency
                *efficiencyCalc;

      Bool_t
          bPassed;
      
    public:
      EfficiencyAnalysis(TH1 *sigTotal, std::ifstream *cuts_file);


  };
};

#endif // !EFFICIENCY_ANALYSIS_H