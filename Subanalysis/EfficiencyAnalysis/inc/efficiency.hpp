#ifndef EFFICIENCY_HPP
#define EFFICIENCY_HPP

#include <iostream>
#include <fstream>
#include <vector>

#include <TTree.h>
#include <TFile.h>
#include <TLegend.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TMultiGraph.h>
#include <TAxis.h>
#include <TGaxis.h>
#include <TObjArray.h>
#include <TObjString.h>

#include <const.h>
#include <ErrorLogs.h>
#include <MainMenu.h>

int efficiencyScan(UInt_t, UInt_t);

#endif //! EFFICIENCY_HPP