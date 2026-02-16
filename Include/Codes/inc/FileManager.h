#ifndef FILE_MANAGER_H
#define FILE_MANAGER_H

#include <TChain.h>
#include <ErrorLogs.h>
#include <vector>
#include <string>
#include <MainMenu.h>
#include <TChainElement.h>

namespace KLOE
{
  struct RunStats
  {
    int minRun;
    int maxRun;
    size_t fileCount;
    Long64_t totalEvents; // Całkowita liczba zdarzeń
    std::vector<int> runList;
    double totalLuminosity; // Całkowita luminozność w nb^-1
  };

  class FileManager
  {

  private:
    ErrorHandling::ErrorLogs &_logger;

  public:
    FileManager(ErrorHandling::ErrorLogs &logger) : _logger(logger) {}

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
    void chainInit(TChain &chain_init, Controls::DataType &dataTypeOpt, UInt_t &firstData, UInt_t &lastData, UInt_t &firstMC, UInt_t &lastMC, Int_t csFlag, Int_t oldAnaFlag);

    /**
     * @brief Method to initialize a TChain with the list of files. Defined per branch, which files (Prod2ntu, Prod2root, Old analysis) are to be taken into account.
     * @param chain_init address to the externally initialized TChain object
     * @param logger address to the ErrorHandling::ErrorLogs object to log the errors / infos during analysis operation
     */
    void chainInit(TChain &chain_init);

    /**
     * @brief Get statistics about runs in a directory matching a regex pattern.
     * @param directory Path to the directory to search.
     * @param regex_pattern Regex pattern with a capturing group for the run number.
     * @return RunStats structure with min/max run, file count, total events, and run list.
     * @note This function only scans filenames (fast). Use UpdateRunStatsFromChain() to calculate events and luminosity.
     */
    RunStats getRunStats(const std::string &directory, const std::string &regex_pattern);

    /**
     * @brief Update RunStats with actual event counts and luminosity from a TChain.
     * @param stats RunStats structure to update.
     * @param chain TChain containing the files.
     * @note Call this after chainInit() if you need accurate event counts and luminosity.
     */
    static void UpdateRunStatsFromChain(RunStats &stats, TChain &chain);

    /**
     * @brief Initialize a TChain with files in a given run range (using regex and directory).
     * @param chain_init Reference to the TChain object.
     * @param logger Reference to the ErrorHandling::ErrorLogs object.
     * @param directory Directory to search for files.
     * @param regex_pattern Regex pattern to extract run number.
     * @param minRun Minimal run number to include.
     * @param maxRun Maximal run number to include.
     */
    void chainInit(TChain &chain_init,
                   const std::string &directory, const std::string &regex_pattern,
                   int minRun, int maxRun);

    void chainInit(TChain &chain_init,
                   const std::vector<std::string> &fileList, const std::string &regex_pattern,
                   int minRun, int maxRun);

    /**
     * @brief Przelicz liczbę zdarzeń na luminozność.
     * @param nEvents Liczba zdarzeń.
     * @return Luminozność w nb^-1.
     */
    static double EventsToLuminosity(Long64_t nEvents);

    /**
     * @brief Loguj informacje o plikach w TChain wraz z luminozością.
     * @param chain Reference do TChain.
     * @param logFile Ścieżka do pliku logu.
     */
    static void LogChainLuminosity(TChain &chain, ErrorHandling::ErrorLogs &logger, const std::string &logFile = "input_files.log");

    /**
     * @brief Wczytaj listę ścieżek plików z pliku tekstowego.
     * @param filePath Ścieżka do pliku zawierającego listę plików.
     * @return Wektor ścieżek plików.
     * @throws std::runtime_error jeśli plik nie może być otwarty.
     */
    static std::vector<std::string> LoadFileListFromFile(const std::string &filePath);

    /**
     * @brief Inicjalizuj TChain z listy ścieżek plików.
     * @param chain Reference do TChain.
     * @param logger Reference do ErrorHandling::ErrorLogs.
     * @param fileList Wektor ścieżek plików.
     */
    void chainInit(TChain &chain,
                   const std::vector<std::string> &fileList);

    /**
     * @brief Sprawdź czy nazwa pliku ma odpowiedni format: job_v{wersja}_{typ}_{luminosity}_inv_pb_{numer}.txt
     * @param filename Nazwa pliku do sprawdzenia.
     * @return true jeśli plik ma prawidłowy format, false w przeciwnym razie.
     */
    static bool ValidateJobListFilename(const std::string &filename);
  };
}

#endif // FILE_MANAGER_H
