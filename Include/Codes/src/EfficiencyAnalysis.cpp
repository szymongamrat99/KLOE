#include <EfficiencyAnalysis.h>

KLOE::EfficiencyAnalysis::EfficiencyAnalysis(TH1 *sigTotal, std::ifstream *cuts_file) : sigTotal(sigTotal), cuts_file(cuts_file)
{
  efficiencyCalc = new TEfficiency("Efficiency", "SignalEfficiency", 181, -90, 90);

  
}
