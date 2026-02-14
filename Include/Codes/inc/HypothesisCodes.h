#pragma once

namespace KLOE
{
  enum class HypothesisCode
  {
    SIGNAL = 1,  // For Pi+Pi-Pi0Pi0
    OMEGAPI = 2, // For Omega-Pi0
    FOUR_PI = 3, // For Pi+Pi-Pi+Pi-
    SIMONA_ANALYSIS = 4, // For Simona's analysis
    THREE_PI0 = 5, // For 3Pi0

    INVALID_VALUE // For invalid values
  };

  enum class TrilaterationCode
  {
    THREE_PI0 = 1, // For 3Pi0
    TWO_PI0 = 2,  // For 2Pi0

    INVALID_VALUE // For invalid values
  };

}