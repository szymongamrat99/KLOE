#ifndef ERROR_LOGS_H
#define ERROR_LOGS_H

#include <fstream>
#include <iostream>
#include <map>
#include <unordered_map>
#include <chrono>
#include <boost/filesystem.hpp>

#include <TString.h>

/**
 * @namespace ErrorHandling
 * @brief Subspace for Error Handling methods
 */
namespace ErrorHandling
{

  // Funkcja pomocnicza do wycinania samej nazwy pliku
  // Używamy constexpr, aby kompilator mógł to zoptymalizować
  inline const char *trimFilePath(const char *path)
  {
    const char *file = strrchr(path, '/');
#ifdef _WIN32 // Jeśli pracujesz też na Windowsie
    if (!file)
      file = strrchr(path, '\\');
#endif
    return file ? file + 1 : path;
  }

  // Preprocessor macros
  // 1. Standardowy log (np. do pętli zdarzeń, limit 10)
#define LOG_EVENT(logger, code, info, logType) \
  logger.getErrLog(code, info, ErrorHandling::trimFilePath(__FILE__), __LINE__, 10, -1, logType)

// 2. Log modułowy (wyższy limit dla ważniejszych kroków, np. 100)
#define LOG_MODULE(logger, code, info, logType) \
  logger.getErrLog(code, info, ErrorHandling::trimFilePath(__FILE__), __LINE__, 100, -1, logType)

// 3. Log krytyczny/zawsze (limit 0 = wyłączony)
#define LOG_CRITICAL(logger, code, info, logType) \
  logger.getErrLog(code, info, ErrorHandling::trimFilePath(__FILE__), __LINE__, 0, -1, logType)

// 4. Log błędu fizycznego (zliczany osobno per mctruth)
#define LOG_PHYSICS_ERROR(logger, code, mctruth, logType) \
  logger.getErrLog(code, "", ErrorHandling::trimFilePath(__FILE__), __LINE__, 0, mctruth, logType)

  /**
   * @enum ErrorCodes
   * @brief Error codes for the analysis. Include technical problems, as well as, the mathematical errors.
   */
  enum class ErrorCodes
  {
    DATA_TYPE = 100,             /*!< Improper data type*/
    RANGE = 101,                 /*!< Number out of range*/
    MENU_RANGE = 102,            /*!< Choice out of menu range*/
    FILE_NOT_EXIST = 103,        /*!< File does not exist*/
    TREE_NOT_EXIST = 104,        /*!< TTree does not exist*/
    NULL_POINTER = 105,          /*!< Null Pointer Exception*/
    INITIALIZATION_FAILED = 106, /*!< Initialization failed*/

    DELTA_LT_ZERO = 200,         /*!< Negative delta of a quadratic equation*/
    DENOM_EQ_ZERO = 201,         /*!< Denominator of a fraction equal to zero*/
    DET_ZERO = 202,              /*!< Matrix determinant equal to zero*/
    NAN_VAL = 203,               /*!< NaN value in return*/
    CHI_SQR_STEP = 204,          /*!< Step of minimization too small*/
    INVALID_PHOTON_NUMBER = 205, /*!< Invalid number of photons for neutral reconstruction*/
    FOUR_MOM_NOT_FILLED = 206,   /*!< Four momentum of a particle not filled*/

    // Physics-related errors
    NO_VTX_WITH_TWO_TRACKS = 300,                   /*!< No charged vertices or tracks available for analysis*/
    LESS_THAN_FOUR_NEUTRAL_CLUSTERS = 301,          /*!< Less than four neutral clusters available for analysis*/
    LESS_THAN_SIX_NEUTRAL_CLUSTERS = 302,           /*!< Less than six neutral clusters available for analysis*/
    NO_VTX_WITH_OPPOSITE_TRACKS = 303,              /*!< Not enough charged tracks for analysis*/
    LESS_THAN_FOUR_CLUSTERS_WITH_GOOD_ENERGY = 304, /*!< Less than four neutral clusters with good energy available for analysis*/
    NO_TWO_VTX_WITH_TWO_TRACKS = 305,               /*!< No charged vertices or tracks available for analysis*/
    CHARGED_KAON_MASS_PRE = 306,                    /*!< Did not pass charged kaon invariant mass in preselection*/
    TRILATERATION_KIN_FIT = 307,                    /*!< Did not pass the trilateration kin fit*/
    NO_VALID_SIX_GAMMA_SOLUTION = 308,              /*!< Did not find a valid six photon solution*/
    TRIANGLE_REC = 309,                             /*!< Did not pass the triangle reconstruction*/
    SIGNAL_KIN_FIT = 310,                           /*!< Did not pass the signal kin fit*/
    OMEGA_KIN_FIT = 311,                            /*!< Did not pass the omega kin fit*/

    CUT_CHI2_SIGNAL = 400,  /*!< Did not pass chi2 cut for signal hypothesis*/
    CUT_COMB_MPI0 = 401,    /*!< Did not pass combined mpi0 cut*/
    CUT_TRCV = 402,         /*!< Did not pass TrcSum cut*/
    CUT_INV_MASS_KCH = 403, /*!< Did not pass invariant mass cut for charged kaon*/
    CUT_INV_MASS_KNE = 404, /*!< Did not pass invariant mass cut for neutral kaon*/
    CUT_QMISS = 405,        /*!< Did not pass Qmiss cut*/

    NOT_RECOGNIZED = 666, /*!< Unexpected exception*/

    NO_ERROR = 0, /*!< No error occurred*/

    NUM_CODES = 16 /*!< Number of codes for loops*/
  };

  /**
   * @enum InfoCodes
   * @brief Info codes for the analysis. Include function execution, execution time, etc.
   */
  enum class InfoCodes
  {
    FILE_ADDED = 303,            /*!< File added to the chain*/
    FUNC_EXEC_TIME = 304,        /*!< Execution time of a single function*/
    FUNC_EXECUTED = 305,         /*!< Analysis step of a given kind executed*/
    VARIABLES_INITIALIZED = 306, /*!< Variables initialized*/
    CONFIG_LOADED = 307,         /*!< Configuration loaded from file*/
    PRETTY_PRINTED = 308,       /*!< Pretty print of configuration to screen*/

    NUM_CODES = 5 /*!< Number of codes for loops*/
  };

  struct LogFiles
  {
    enum class LogType
    {
      CUT_LIMITS,
      GENERAL,
      ANALYSIS_CONFIG,
      WEB_SERVICE,
      ERROR,
      PHYSICS_CONSTANTS
    };

    enum class LogLevel
    {
      INFO,
      WARNING,
      ERROR
    };

    const std::map<LogType, std::string> logFileNames = {
        {LogType::CUT_LIMITS, "cut.limits.log"},
        {LogType::GENERAL, "general.log"},
        {LogType::ANALYSIS_CONFIG, "analysis.config.log"},
        {LogType::WEB_SERVICE, "rti.web-service.log"},
        {LogType::ERROR, "error.log"},
        {LogType::PHYSICS_CONSTANTS, "physics.constants.log"}};
  };

  /**
   * @class ErrorLogs
   * @brief Class for info and error logs output - both to screen and log files.
   */
  class ErrorLogs
  {
  private:
    std::map<ErrorCodes, TString>
        Logs; /*!< Dictionary of error logs messages*/

    std::map<ErrorCodes, int>
        _errCount; /*!< Counter to get the summary of thrown exceptions*/

    std::string
        _logDirectory,
        _currentLogDir; /*!< Base name for log files, to which the type-specific suffixes will be added*/

    LogFiles::LogLevel _logLevel; /*!< Log level for screen output*/

    std::map<LogFiles::LogType, std::ofstream>
        _logFile; /*!< Log file object*/

    std::unordered_map<std::string, Long64_t>
        _spamCounter; /*!< Number of errors*/

    bool _printToScreen = false; /*!< If true, print logs to screen */

    // Counter for physics-related errors per mctruth value (codes 300-399)
    std::map<int, std::map<ErrorCodes, int>> _physicsErrCountPerMctruth;

    /**
     * @brief private method to get the error message with ErrorCodes parameter.
     * @param code code from ErrorCodes enumeration
     * @returns the error message tied with a given code
     */
    std::string _getErrorMessage(ErrorCodes code) const
    {
      switch (code)
      {
      // General logs
      case ErrorCodes::DATA_TYPE:
        return "Invalid input data type.";
      case ErrorCodes::RANGE:
        return "File number outside available range.";
      case ErrorCodes::MENU_RANGE:
        return "Chosen menu option outside range.";
      case ErrorCodes::FILE_NOT_EXIST:
        return "File does not exist. Check the name and path.";
      case ErrorCodes::TREE_NOT_EXIST:
        return "Tree does not exist. Check the name.";
      case ErrorCodes::NULL_POINTER:
        return "Null pointer exception. Initialize the object.";

      // Math logs
      case ErrorCodes::DELTA_LT_ZERO:
        return "Quadratic delta less than 0!";
      case ErrorCodes::DENOM_EQ_ZERO:
        return "Denominator in fraction equal to 0!";
      case ErrorCodes::DET_ZERO:
        return "Determinant of a matrix equal to 0!";
      case ErrorCodes::NAN_VAL:
        return "Value is NaN!";
      case ErrorCodes::CHI_SQR_STEP:
        return "Chi-squared step less than threshold!";
      case ErrorCodes::INVALID_PHOTON_NUMBER:
        return "Invalid number of photons for neutral reconstruction.";
      case ErrorCodes::FOUR_MOM_NOT_FILLED:
        return "Four momentum of a particle not filled.";

      // Physics-related logs
      case ErrorCodes::NO_VTX_WITH_OPPOSITE_TRACKS:
        return "No vertices with two opposite tracks available for analysis.";
      case ErrorCodes::LESS_THAN_FOUR_NEUTRAL_CLUSTERS:
        return "Less than four neutral clusters available for analysis.";
      case ErrorCodes::LESS_THAN_SIX_NEUTRAL_CLUSTERS:
        return "Less than six neutral clusters available for analysis.";
      case ErrorCodes::CHARGED_KAON_MASS_PRE:
        return "Did not pass charged kaon ivariant mass in preselection.";
      case ErrorCodes::TRILATERATION_KIN_FIT:
        return "Did not pass the trilateration kin fit.";
      case ErrorCodes::NO_VALID_SIX_GAMMA_SOLUTION:
        return "Did not find a valid six photon solution.";
      case ErrorCodes::TRIANGLE_REC:
        return "Did not pass the triangle reconstruction.";
      case ErrorCodes::SIGNAL_KIN_FIT:
        return "Did not pass the signal global kin fit.";
      case ErrorCodes::NO_VTX_WITH_TWO_TRACKS:
        return "No vertices with two tracks available for analysis.";
      case ErrorCodes::NO_TWO_VTX_WITH_TWO_TRACKS:
        return "No two vertices with two tracks available for analysis.";
      case ErrorCodes::LESS_THAN_FOUR_CLUSTERS_WITH_GOOD_ENERGY:
        return "Less than four neutral clusters with good energy available for analysis.";
      case ErrorCodes::OMEGA_KIN_FIT:
        return "Did not pass the omega kin fit.";

      // Cut-related logs
      case ErrorCodes::CUT_CHI2_SIGNAL:
        return "Did not pass chi2 cut for signal hypothesis.";
      case ErrorCodes::CUT_COMB_MPI0:
        return "Did not pass combined mpi0 cut.";
      case ErrorCodes::CUT_TRCV:
        return "Did not pass TrcSum cut.";
      case ErrorCodes::CUT_INV_MASS_KCH:
        return "Did not pass invariant mass cut for charged kaon.";
      case ErrorCodes::CUT_INV_MASS_KNE:
        return "Did not pass invariant mass cut for neutral kaon.";
      case ErrorCodes::CUT_QMISS:
        return "Did not pass Qmiss cut.";

      // Not recognized logs
      case ErrorCodes::NOT_RECOGNIZED:
        return "Error not recognized.";

      default:
        return "Improper error code.";
      }
    }

    /**
     * @brief private method to get the info message with InfoCodes parameter.
     * @param code code from InfoCodes enumeration
     * @returns the info message tied with a given code
     */
    std::string _getInfoMessage(InfoCodes code) const
    {
      switch (code)
      {
      // Informative logs
      case InfoCodes::FILE_ADDED:
        return "File added to the ongoing analysis";
      case InfoCodes::FUNC_EXEC_TIME:
        return "Execution time";
      case InfoCodes::FUNC_EXECUTED:
        return "Executed part of the analysis";
      case InfoCodes::VARIABLES_INITIALIZED:
        return "Variables initialized";
      case InfoCodes::CONFIG_LOADED:
        return "Configuration loaded";
      case InfoCodes::PRETTY_PRINTED:
        return "";

      default:
        return "Improper info code.";
      }
    }

    std::string _createDatedDir()
    {
      auto now = std::chrono::system_clock::now();
      std::time_t currentTime = std::chrono::system_clock::to_time_t(now);
      char buffer[100];
      std::strftime(buffer, sizeof(buffer), "%Y-%m-%d", std::localtime(&currentTime));
      std::string datedDir = _logDirectory + std::string(buffer) + "/";
      return datedDir;
    }

    /**
     * @brief private method to get the timestamp for log generation
     * @returns the timestamp in the format: %Y-%m-%d %H\:%M\:%S
     */
    std::string _getTimestamp()
    {
      auto now = std::chrono::system_clock::now();
      std::time_t currentTime = std::chrono::system_clock::to_time_t(now);
      char buffer[100];
      std::strftime(buffer, sizeof(buffer), "%Y-%m-%d %H:%M:%S", std::localtime(&currentTime));
      return std::string(buffer);
    }

    void _OpenLogFile(LogFiles::LogType logType)
    {
      std::string logFilePath = _currentLogDir + LogFiles().logFileNames.at(logType);
      if (_logFile.find(logType) == _logFile.end())
      {
        _logFile[logType].open(logFilePath, std::ios::out | std::ios::app);
        if (!_logFile[logType].is_open())
        {
          std::cerr << "Failed to open log file: " << logFilePath << std::endl;
          return;
        }
      }
    }

    void _PrintStatistics()
    {
      _OpenLogFile(LogFiles::LogType::ERROR);

      _logFile[LogFiles::LogType::ERROR] << "Physics-related error statistics per mctruth value:" << std::endl;
      for (const auto &mctruthPair : _physicsErrCountPerMctruth)
      {
        int mctruth = mctruthPair.first;
        const auto &errCounts = mctruthPair.second;
        _logFile[LogFiles::LogType::ERROR] << "MC Truth Value: " << mctruth << std::endl;

        std::cout << "MC Truth Value: " << mctruth << std::endl;
        for (const auto &errPair : errCounts)
        {
          ErrorCodes errCode = errPair.first;
          int count = errPair.second;
          _logFile[LogFiles::LogType::ERROR] << "  " << _getErrorMessage(errCode) << ": " << count << " occurrences" << std::endl
                                             << std::endl;
          std::cout << "  " << _getErrorMessage(errCode) << ": " << count << " occurrences" << std::endl;
        }
      }

      _logFile[LogFiles::LogType::ERROR] << "Error counts (all):" << std::endl;
      std::cout << "Error counts (all):" << std::endl;
      for (const auto &pair : _spamCounter)
      {
        _logFile[LogFiles::LogType::ERROR] << "  " << pair.first << ": " << pair.second << " occurrences" << std::endl;
        std::cout << "  " << pair.first << ": " << pair.second << " occurrences" << std::endl;
      }

      _logFile[LogFiles::LogType::ERROR].flush();

      _logFile[LogFiles::LogType::ERROR] << "End of error statistics." << std::endl
                                         << std::endl
                                         << std::endl;
      std::cout << "End of error statistics." << std::endl;
    }

  public:
    // Constructor
    ErrorLogs(const std::string &logDirectory = "") : _logDirectory(logDirectory)
    {
      if (_logDirectory.empty())
      {
        std::cerr << "Log directory not provided, ending." << std::endl;
        return;
      }

      _currentLogDir = _createDatedDir(); /*!< Create a directory with the current date */

      // Create the directory if it does not exist
      if (!boost::filesystem::exists(_currentLogDir))
      {
        if (!boost::filesystem::create_directories(_currentLogDir))
        {
          std::cerr << "Failed to create log directory: " << _currentLogDir << std::endl;
          return;
        }
      }
    }

    ~ErrorLogs()
    {
      _PrintStatistics();

      for (auto &pair : _logFile)
      {
        if (pair.second.is_open())
        {
          pair.second.close();
        }
      }
    }

    /**
     * @brief Increment physics-related error counter for a given mctruth value.
     * @param errCode Error code (should be in 300-399)
     * @param mctruth MC truth value
     */
    void countPhysicsError(ErrorCodes errCode, int mctruth)
    {
      if ((Int_t)errCode >= 300 && (Int_t)errCode <= 399)
        _physicsErrCountPerMctruth[mctruth][errCode]++;
    }

    void getErrLog(ErrorCodes &errCode, const std::string &additionalInfo = "", const std::string file = "Unknown", int line = 0, int limit = 10, int mctruth = -1, LogFiles::LogType logType = LogFiles::LogType::ERROR)
    {
      if (errCode == ErrorCodes::NO_ERROR)
        return; // No error, nothing to log

      _OpenLogFile(logType);

      std::string key = file + ":" + std::to_string(line) + ":" + std::to_string((Int_t)errCode);
      _spamCounter[key]++;

      Bool_t shouldLog = (limit <= 0) || (_spamCounter[key] <= limit);

      if (shouldLog)
      {
        std::string
            timestamp = _getTimestamp(),
            errorMessage = _getErrorMessage(errCode),
            concatMessage;

        if (additionalInfo.empty())
          concatMessage = "[" + timestamp + "] (" + file + ":" + std::to_string(line) + ") Error: " + errorMessage;
        else
          concatMessage = "[" + timestamp + "] (" + file + ":" + std::to_string(line) + ") Error: " + errorMessage + " - " + additionalInfo;

        _errCount[errCode]++;

        Bool_t isPhysicsError = ((Int_t)errCode >= 300 && (Int_t)errCode <= 399);
        // Zliczaj błędy fizyczne dla każdego wywołania getErrLog jeśli mctruth podany
        if (mctruth != -1)
          countPhysicsError(errCode, mctruth);

        if (_printToScreen && !isPhysicsError)
          std::cerr << concatMessage << std::endl;

        if (_logFile.at(logType).is_open() && !isPhysicsError)
        {
          _logFile.at(logType) << concatMessage << std::endl;
        }
      }
      else if (limit > 0 && _spamCounter[key] == limit + 1)
      {
        std::string timestamp = _getTimestamp();
        std::string logMessage = "[" + timestamp + "] (" + file + ":" + std::to_string(line) + ") Error: Further occurrences of this error will be suppressed. Original message: " + _getErrorMessage(errCode);
        if (_printToScreen)
          std::cerr << logMessage << std::endl;
        if (_logFile.at(logType).is_open())
        {
          _logFile.at(logType) << logMessage << std::endl;
        }
      }
    };

    void setPrintToScreen(bool print) { _printToScreen = print; }
    bool getPrintToScreen() const { return _printToScreen; }

    std::string getCurrentLogDir() const { return _currentLogDir; }

    void getLog(InfoCodes &infoCode, const std::string &additionalInfo = "", LogFiles::LogType logType = LogFiles::LogType::GENERAL)
    {
      _OpenLogFile(logType);

      std::string
          timestamp = _getTimestamp(),
          infoMessage = _getInfoMessage(infoCode),
          concatMessage;

      if (additionalInfo.empty())
        concatMessage = "[" + timestamp + "] Info: " + infoMessage;
      else
        concatMessage = "[" + timestamp + "] Info: " + infoMessage + " - " + additionalInfo;

      if (_printToScreen)
        std::cerr << concatMessage << std::endl;

      if (_logFile.at(logType).is_open())
      {
        _logFile.at(logType) << concatMessage << std::endl;
      }
    };

    void prettyPrint(const std::string &additionalInfo = "", const std::string &beginMsg = "", const std::string &endMsg = "", LogFiles::LogType logType = LogFiles::LogType::ANALYSIS_CONFIG)
    {
      _OpenLogFile(logType);

      std::string
          timestamp = _getTimestamp(),
          concatMessageBegin,
          concatMessageEnd;

      concatMessageBegin = "[" + timestamp + "] Info: " + beginMsg + "\n";
      concatMessageEnd = "\n[" + timestamp + "] Info: " + endMsg;

      if (_printToScreen)
        std::cerr << concatMessageBegin << additionalInfo << concatMessageEnd << std::endl;

      if (_logFile.at(logType).is_open())
      {
        _logFile.at(logType) << concatMessageBegin << additionalInfo << concatMessageEnd << std::endl;
      }
    };

    // Funkcja wypisująca statystyki
    void printErrStats() const
    {
      std::cout << "Error statistics:" << std::endl;
      for (const auto &pair : _errCount)
      {
        std::cout << _getErrorMessage(pair.first) << ": " << pair.second << " occurrences" << std::endl;
      }
    }

    /**
     * @brief Print physics-related error statistics per mctruth to the log file and optionally to screen.
     * @param printToScreen If true, also print to screen (default: false)
     */
    void printPhysicsErrorStatsPerMctruth(bool printToScreen = false, LogFiles::LogType logType = LogFiles::LogType::ERROR)
    {
      _OpenLogFile(logType);

      if (!_logFile.at(logType).is_open())
        return;

      _logFile.at(logType) << "Physics-related error statistics per mctruth:\n";
      for (const auto &mctruth_pair : _physicsErrCountPerMctruth)
      {
        int mctruth = mctruth_pair.first;
        _logFile.at(logType) << "mctruth = " << mctruth << ":\n";
        if (printToScreen)
          std::cout << "mctruth = " << mctruth << ":\n";
        for (const auto &err_pair : mctruth_pair.second)
        {
          _logFile.at(logType) << "  " << _getErrorMessage(err_pair.first) << ": " << err_pair.second << " occurrences\n";
          if (printToScreen)
            std::cout << "  " << _getErrorMessage(err_pair.first) << ": " << err_pair.second << " occurrences\n";
        }
      }
      _logFile.at(logType).flush();
    }

    /**
     * @brief Get map of physics error counts for a given mctruth value.
     * @param mctruth MC truth value
     * @return Map of ErrorCodes to error counts for the given mctruth (empty if not present)
     */
    std::map<ErrorCodes, int> getPhysicsErrorCountsForMctruth(int mctruth) const
    {
      auto it = _physicsErrCountPerMctruth.find(mctruth);
      if (it != _physicsErrCountPerMctruth.end())
        return it->second;
      else
        return {};
    }
  };
}

#endif