#ifndef FILE_MANAGER_H
#define FILE_MANAGER_H

#include <TChain.h>
#include <ErrorLogs.h>
#include <vector>
#include <string>

namespace KLOE {
struct RunStats {
    int minRun;
    int maxRun;
    size_t fileCount;
    double totalGB;
    std::vector<int> runList;
};

class FileManager {
public:
    /**
     * @brief Method to initialize a TChain with the list of files. Defined per branch, which files (Prod2ntu, Prod2root, Old analysis) are to be taken into account.
     * @param chain_init address to the externally initialized TChain object
     * @param dataTypeOpt address to the variable, which stores Controls::DataType value
     * @param firstData first data file index
     * @param lastData last data file index
     * @param firstMC first MC file index
     * @param lastMC last MC file index
     * @param logger address to the ErrorHandling::ErrorLogs object to log the errors / infos during analysis operation
     * @param csFlag control sample flag
     */
    void chainInit(TChain &chain_init, Controls::DataType &dataTypeOpt, UInt_t &firstData, UInt_t &lastData, UInt_t &firstMC, UInt_t &lastMC, ErrorHandling::ErrorLogs &logger, Int_t csFlag);

    /**
     * @brief Method to initialize a TChain with the list of files. Defined per branch, which files (Prod2ntu, Prod2root, Old analysis) are to be taken into account.
     * @param chain_init address to the externally initialized TChain object
     * @param logger address to the ErrorHandling::ErrorLogs object to log the errors / infos during analysis operation
     */
    void chainInit(TChain &chain_init, ErrorHandling::ErrorLogs &logger);

    /**
     * @brief Get statistics about runs in a directory matching a regex pattern.
     * @param directory Path to the directory to search.
     * @param regex_pattern Regex pattern with a capturing group for the run number.
     * @return RunStats structure with min/max run, file count, total GB, and run list.
     */
    RunStats getRunStats(const std::string& directory, const std::string& regex_pattern);

    /**
     * @brief Initialize a TChain with files in a given run range (using regex and directory).
     * @param chain_init Reference to the TChain object.
     * @param logger Reference to the ErrorHandling::ErrorLogs object.
     * @param directory Directory to search for files.
     * @param regex_pattern Regex pattern to extract run number.
     * @param minRun Minimal run number to include.
     * @param maxRun Maximal run number to include.
     */
    void chainInit(TChain &chain_init, ErrorHandling::ErrorLogs &logger,
                  const std::string& directory, const std::string& regex_pattern,
                  int minRun, int maxRun);
};
}

#endif // FILE_MANAGER_H
