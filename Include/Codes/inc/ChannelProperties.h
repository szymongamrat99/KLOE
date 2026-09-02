#pragma once

namespace KLOE
{
  enum ChannelFlags : unsigned int {
    kNone             = 0,
    kHasPM            = 1u << 0,  // pi+ pi-
    kHas00            = 1u << 1,  // pi0 pi0
    kHas000           = 1u << 2,  // 3 pi0
    kHasSemiLeptonic  = 1u << 3,  // semileptonic from Kaon
    kHasISR           = 1u << 4   // ISR photon in the event
  };
}