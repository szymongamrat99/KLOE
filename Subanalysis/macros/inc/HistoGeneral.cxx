#include "HistoGeneral.h"


namespace HistoGeneral
{
  const std::map<Int_t, TString> channelTypes = {
      {0, "Data"},
      {1, "Signal"},
      {2, "Regeneration"},
      {3, "Omega-pi0"},
      {4, "3pi0"},
      {5, "Semileptonic"},
      {6, "Other"},
      {7, "pi+pi-pi+pi-"},
      {8, "MC sum"}};

  const std::map<TString, Color_t> channelColors = {
      {"Data", kBlack},
      {"Signal", kRed},
      {"Regeneration", kGreen},
      {"Omega-pi0", kViolet},
      {"3pi0", kCyan},
      {"Semileptonic", kBlue},
      {"Other", kGreen - 1},
      {"pi+pi-pi+pi-", kYellow},
      {"MC sum", kOrange + 7}};

  const std::map<TString, TString> histDrawOption = {
      {"Data", "PE"},   // Data
      {"Signal", "HIST"}, // Signal
      {"Regeneration", "HIST"}, // Regeneration
      {"Omega-pi0", "HIST"}, // Omega-pi0
      {"3pi0", "HIST"}, // 3pi0
      {"Semileptonic", "HIST"}, // Semileptonic
      {"Other", "HIST"}, // Other
      {"pi+pi-pi+pi-", "HIST"}, // pi+pi-pi+pi-
      {"MC sum", "HIST"}  // MC sum
  };
}
