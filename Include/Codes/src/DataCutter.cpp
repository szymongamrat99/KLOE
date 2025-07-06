#include <DataCutter.h>
#include <cmath> // For std::abs

namespace KLOE {

/**
 * @class BaseCut
 * @brief Base class for all types of cuts.
 *
 * This class provides the interface for evaluating cuts and tracking statistics.
 */
// No implementation for evaluate() here, as it is a pure virtual method.

/**
 * @class AsymmetricCut
 * @brief Implements an asymmetric cut using a relational operator and a threshold.
 *
 * The cut passes if the value compared to the threshold with the given operator is true.
 */
AsymmetricCut::AsymmetricCut(
    const std::string& cutName,
    RelationalOperator op,
    double threshold,
    bool active
) : BaseCut(cutName, active), m_op(op), m_threshold(threshold) {}

/**
 * @brief Helper function to evaluate the relational operator.
 * @param value Value to compare
 * @param threshold Threshold value
 * @param op Relational operator
 * @return true if the comparison passes, false otherwise
 */
bool AsymmetricCut::evaluateOperator(double value, double threshold, RelationalOperator op) const {
    switch (op) {
        case RelationalOperator::LESS:       return value < threshold;
        case RelationalOperator::LESS_EQUAL: return value <= threshold;
        case RelationalOperator::GREATER:    return value > threshold;
        case RelationalOperator::GREATER_EQUAL: return value >= threshold;
        case RelationalOperator::EQUAL:      return value == threshold;
        case RelationalOperator::NOT_EQUAL:  return value != threshold;
        default: return false; // Should not happen
    }
}

/**
 * @brief Evaluates the asymmetric cut for a given value.
 * @param value Value to evaluate
 * @return true if the value passes the cut, false otherwise
 */
bool AsymmetricCut::evaluate(double value) const {
    return evaluateOperator(value, m_threshold, m_op);
}

/**
 * @class SymmetricCut
 * @brief Implements a symmetric cut around a center value with an absolute threshold.
 *
 * The cut passes if the absolute difference between the value and the center is less than or equal to the threshold.
 */
SymmetricCut::SymmetricCut(
    const std::string& cutName,
    double center,
    double absThreshold,
    bool active
) : BaseCut(cutName, active), m_center(center), m_absThreshold(absThreshold) {}

/**
 * @brief Evaluates the symmetric cut for a given value.
 * @param value Value to evaluate
 * @return true if the value passes the cut, false otherwise
 */
bool SymmetricCut::evaluate(double value) const {
    return std::abs(value - m_center) <= m_absThreshold;
}

/**
 * @class DataCutter
 * @brief Main class for managing and applying multiple data cuts.
 *
 * Provides methods to add cuts, apply them to data, activate/deactivate cuts, print statistics, and reset statistics.
 */
DataCutter::DataCutter()
    : m_eventsAfterAllCuts(0), m_totalProcessedEvents(0) {}

/**
 * @brief Adds an asymmetric cut to the list.
 * @param cutName Name of the cut
 * @param op Relational operator
 * @param threshold Threshold value
 * @param active Whether the cut is active
 */
void DataCutter::addAsymmetricCut(
    const std::string& cutName,
    RelationalOperator op,
    double threshold,
    bool active
) {
    m_cuts.push_back(std::make_unique<AsymmetricCut>(cutName, op, threshold, active));
}

/**
 * @brief Adds a symmetric cut to the list.
 * @param cutName Name of the cut
 * @param center Center value
 * @param absThreshold Absolute threshold
 * @param active Whether the cut is active
 */
void DataCutter::addSymmetricCut(
    const std::string& cutName,
    double center,
    double absThreshold,
    bool active
) {
    m_cuts.push_back(std::make_unique<SymmetricCut>(cutName, center, absThreshold, active));
}

/**
 * @brief Applies all active cuts to the given value.
 * @param value Value to process through the cuts
 * @return true if the value passes ALL active cuts, false otherwise
 */
bool DataCutter::applyCuts(double value) {
    m_totalProcessedEvents++;
    bool allCutsPassed = true;
    for (const auto& cutPtr : m_cuts) {
        if (cutPtr->isActive) {
            if (!cutPtr->evaluate(value)) {
                allCutsPassed = false;
                break;
            } else {
                cutPtr->passedEvents++;
            }
        }
    }

    if (allCutsPassed) {
        m_eventsAfterAllCuts++;
    }
    return allCutsPassed;
}

/**
 * @brief Activates or deactivates a cut by name.
 * @param cutName Name of the cut
 * @param activate True to activate; false to deactivate
 * @return true if the cut was found and its status changed, false otherwise
 */
bool DataCutter::setCutActive(const std::string& cutName, bool activate) {
    for (auto& cutPtr : m_cuts) {
        if (cutPtr->name == cutName) {
            cutPtr->isActive = activate;
            return true;
        }
    }
    std::cerr << "Warning: Cut '" << cutName << "' not found." << std::endl;
    return false;
}

/**
 * @brief Generates and displays cut statistics.
 * @param totalInitialEvents Total number of events before any cuts
 * @param signalEvents Number of signal events in initial events (optional, for S/B metrics)
 * @param backgroundEvents Number of background events in initial events (optional)
 * @param outputFilePath Path to the file where statistics will be saved
 */
void DataCutter::printStatistics(
    size_t totalInitialEvents,
    size_t signalEvents,
    size_t backgroundEvents,
    const std::string& outputFilePath
) const {
    std::ostream* out = &std::cout;
    std::ofstream outFile;

    if (!outputFilePath.empty()) {
        outFile.open(outputFilePath);
        if (outFile.is_open()) {
            out = &outFile;
        } else {
            std::cerr << "Warning: Could not open output file '" << outputFilePath << "'. Printing to console only." << std::endl;
        }
    }

    *out << "--- Data Cutting Statistics ---" << std::endl;
    *out << std::fixed << std::setprecision(2);

    *out << "Total Events Processed: " << m_totalProcessedEvents << std::endl;
    *out << "Total Events After All Active Cuts: " << m_eventsAfterAllCuts << std::endl;

    if (totalInitialEvents > 0) {
        double overallEfficiency = (double)m_eventsAfterAllCuts / totalInitialEvents * 100.0;
        *out << "Overall Efficiency (relative to initial events): " << overallEfficiency << "%" << std::endl;
    } else {
        *out << "Overall Efficiency (relative to initial events): N/A (Total initial events not provided or is zero)" << std::endl;
    }

    *out << "\n--- Individual Cut Statistics (active cuts only) ---" << std::endl;
    for (const auto& cutPtr : m_cuts) {
        if (cutPtr->isActive) {
            *out << "Cut Name: " << cutPtr->name << std::endl;
            *out << "  Passed Events: " << cutPtr->passedEvents << std::endl;
            if (m_totalProcessedEvents > 0) {
                double cutEfficiency = (double)cutPtr->passedEvents / m_totalProcessedEvents * 100.0;
                *out << "  Efficiency (relative to processed events): " << cutEfficiency << "%" << std::endl;
            } else {
                *out << "  Efficiency (relative to processed events): N/A (No events processed yet)" << std::endl;
            }
        }
    }

    if (signalEvents > 0 && backgroundEvents > 0) {
        *out << "\n--- Signal vs Background Estimation (REQUIRES CAREFUL SETUP) ---" << std::endl;
        *out << "Initial Signal Events: " << signalEvents << std::endl;
        *out << "Initial Background Events: " << backgroundEvents << std::endl;
        *out << "NOTE: To get proper S/B and signal/background efficiencies, 'applyCuts' needs" << std::endl;
        *out << "to be modified to distinguish between signal and background events during processing." << std::endl;
    } else if (signalEvents > 0 || backgroundEvents > 0) {
        *out << "\nNOTE: Both signalEvents and backgroundEvents must be provided for S/B statistics." << std::endl;
    }

    *out << "------------------------------" << std::endl;

    if (outFile.is_open()) {
        outFile.close();
    }
}

/**
 * @brief Resets all cut statistics.
 */
void DataCutter::resetStatistics() {
    m_eventsAfterAllCuts = 0;
    m_totalProcessedEvents = 0;
    for (auto& cutPtr : m_cuts) {
        cutPtr->passedEvents = 0;
    }
    std::cout << "DataCutter statistics have been reset." << std::endl;
}

} // namespace KLOE