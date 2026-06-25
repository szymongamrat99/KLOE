#ifndef BASE_HYPOTHESIS_H
#define BASE_HYPOTHESIS_H

#include <ErrorLogs.h>
#include <event_context.h>

class BaseHypothesis
{
public:
    BaseHypothesis() = default;

    virtual ~BaseHypothesis() = default;

    virtual ErrorHandling::ErrorCodes Process(
        EventContext& ctx
    ) = 0;
};

#endif // BASE_HYPOTHESIS_H