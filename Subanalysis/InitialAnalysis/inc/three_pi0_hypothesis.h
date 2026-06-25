#ifndef THREE_PI0_HYPOTHESIS_H
#define THREE_PI0_HYPOTHESIS_H

#include <base_hypothesis.h>

class ThreePi0Hypothesis : public BaseHypothesis
{
public:
    ThreePi0Hypothesis() = default;
    ~ThreePi0Hypothesis() override = default;

    ErrorHandling::ErrorCodes Process(EventContext& ctx) override;

    private:
      int minClusters = 6;
      double minClusterEnergy = 20.0;

      bool noError = true;

      ErrorHandling::ErrorCodes errorCode = ErrorHandling::ErrorCodes::NO_ERROR;

      bool CheckErrorCode(ErrorHandling::ErrorCodes code)
      {
        if (code != ErrorHandling::ErrorCodes::NO_ERROR)
        {
          // Handle the error, e.g., log it or take appropriate action
          return false;
        }
        return true;
      }

};

#endif // THREE_PI0_HYPOTHESIS_H