#pragma once

#include <map>

#include <TString.h>
#include <const.h>

namespace HistoGeneral
{
  extern const std::map<Int_t, TString> channelTypes;
  extern const std::map<TString, Color_t> channelColors;
  extern const std::map<TString, TString> histDrawOption;
}