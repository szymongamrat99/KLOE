#ifndef DATA_CUTTER_H
#define DATA_CUTTER_H

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip> // For std::setprecision, std::fixed
#include <memory>  // For std::unique_ptr

namespace KLOE {

/**
 * @brief Enumeration for relational operator types.
 */
enum class RelationalOperator {
    LESS,           ///< <
    LESS_EQUAL,     ///< <=
    GREATER,        ///< >
    GREATER_EQUAL,  ///< >=
    EQUAL,          ///< ==
    NOT_EQUAL       ///< !=
};

/**
 * @brief Base class for all cuts.
 */
class BaseCut {
public:
    std::string name;           ///< Name of the cut
    bool isActive;              ///< Whether the cut is active
    size_t passedEvents;        ///< Number of events that passed this cut

    /**
     * @brief Constructor for BaseCut.
     * @param n Name of the cut
     * @param active Whether the cut is active
     */
    BaseCut(const std::string& n, bool active = true)
        : name(n), isActive(active), passedEvents(0) {}

    /**
     * @brief Virtual method to evaluate the cut for a given value.
     * @param value Value to evaluate
     * @return true if the value passes the cut, false otherwise
     */
    virtual bool evaluate(double value) const = 0;

    /**
     * @brief Virtual destructor for safe inheritance.
     */
    virtual ~BaseCut() = default;
};

/**
 * @brief Class for asymmetric cuts.
 */
class AsymmetricCut : public BaseCut {
private:
    RelationalOperator m_op;    ///< Relational operator
    double m_threshold;         ///< Threshold value

    /**
     * @brief Helper function to evaluate the operator.
     * @param value Value to compare
     * @param threshold Threshold value
     * @param op Relational operator
     * @return true if the comparison passes, false otherwise
     */
    bool evaluateOperator(double value, double threshold, RelationalOperator op) const;

public:
    /**
     * @brief Constructor for AsymmetricCut.
     * @param cutName Name of the cut
     * @param op Relational operator
     * @param threshold Threshold value
     * @param active Whether the cut is active
     */
    AsymmetricCut(
        const std::string& cutName,
        RelationalOperator op,
        double threshold,
        bool active = true
    );

    /**
     * @brief Evaluates the cut for a given value.
     * @param value Value to evaluate
     * @return true if the value passes the cut, false otherwise
     */
    bool evaluate(double value) const override;
};

/**
 * @brief Class for symmetric cuts.
 */
class SymmetricCut : public BaseCut {
private:
    double m_center;        ///< Center value
    double m_absThreshold;  ///< Absolute threshold

public:
    /**
     * @brief Constructor for SymmetricCut.
     * @param cutName Name of the cut
     * @param center Center value
     * @param absThreshold Absolute threshold
     * @param active Whether the cut is active
     */
    SymmetricCut(
        const std::string& cutName,
        double center,
        double absThreshold,
        bool active = true
    );

    /**
     * @brief Evaluates the cut for a given value.
     * @param value Value to evaluate
     * @return true if the value passes the cut, false otherwise
     */
    bool evaluate(double value) const override;
};

/**
 * @brief Main class for managing and applying data cuts.
 */
class DataCutter {
public:
    /**
     * @brief Constructor for DataCutter.
     */
    DataCutter();

    /**
     * @brief Adds an asymmetric cut.
     * @param cutName Name of the cut
     * @param op Relational operator
     * @param threshold Threshold value
     * @param active Whether the cut is active
     */
    void addAsymmetricCut(
        const std::string& cutName,
        RelationalOperator op,
        double threshold,
        bool active = true
    );

    /**
     * @brief Adds a symmetric cut.
     * @param cutName Name of the cut
     * @param center Center value
     * @param absThreshold Absolute threshold
     * @param active Whether the cut is active
     */
    void addSymmetricCut(
        const std::string& cutName,
        double center,
        double absThreshold,
        bool active = true
    );

    /**
     * @brief Applies all active cuts to the given value.
     * @param value Value to process through the cuts
     * @return true if the value passes ALL active cuts, false otherwise
     */
    bool applyCuts(double value);

    /**
     * @brief Activates or deactivates a cut by name.
     * @param cutName Name of the cut
     * @param activate True to activate; false to deactivate
     * @return true if the cut was found and its status changed, false otherwise
     */
    bool setCutActive(const std::string& cutName, bool activate);

    /**
     * @brief Generates and displays cut statistics.
     * @param totalInitialEvents Total number of events before any cuts
     * @param signalEvents Number of signal events in initial events (optional, for S/B metrics)
     * @param backgroundEvents Number of background events in initial events (optional)
     * @param outputFilePath Path to the file where statistics will be saved
     */
    void printStatistics(
        size_t totalInitialEvents,
        size_t signalEvents = 0,
        size_t backgroundEvents = 0,
        const std::string& outputFilePath = "cut_statistics.txt"
    ) const;

    /**
     * @brief Resets all cut statistics.
     */
    void resetStatistics();

private:
    std::vector<std::unique_ptr<BaseCut>> m_cuts; ///< Vector of unique_ptr to BaseCut objects
    size_t m_eventsAfterAllCuts;                  ///< Number of events that passed all cuts
    size_t m_totalProcessedEvents;                ///< Total number of events processed by applyCuts
};

} // namespace KLOE

#endif // DATA_CUTTER_H