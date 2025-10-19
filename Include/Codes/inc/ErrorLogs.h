#ifndef ERROR_LOGS_H
#define ERROR_LOGS_H

#include <fstream>
#include <iostream>
#include <map>
#include <chrono>

#include <TString.h>

/**
 * @namespace ErrorHandling
 * @brief Subspace for Error Handling methods
 */
namespace ErrorHandling
{
  /**
   * @enum ErrorCodes
   * @brief Error codes for the analysis. Include technical problems, as well as, the mathematical errors.
   */
  enum class ErrorCodes
  {
    DATA_TYPE = 100,      /*!< Improper data type*/
    RANGE = 101,          /*!< Number out of range*/
    MENU_RANGE = 102,     /*!< Choice out of menu range*/
    FILE_NOT_EXIST = 103, /*!< File does not exist*/
    TREE_NOT_EXIST = 104, /*!< TTree does not exist*/
    NULL_POINTER = 105,   /*!< Null Pointer Exception*/

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
    TRIANGLE_REC = 308,                             /*!< Did not pass the triangle reconstruction*/
    SIGNAL_KIN_FIT = 309,                           /*!< Did not pass the signal kin fit*/
    OMEGA_KIN_FIT = 310,                            /*!< Did not pass the omega kin fit*/

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
    FILE_ADDED = 303,     /*!< File added to the chain*/
    FUNC_EXEC_TIME = 304, /*!< Execution time of a single function*/
    FUNC_EXECUTED = 305,  /*!< Analysis step of a given kind executed*/

    NUM_CODES = 3 /*!< Number of codes for loops*/
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

    std::ofstream
        _logFile; /*!< Log file object*/

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

      default:
        return "Improper info code.";
      }
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

  public:
    // Constructor
    ErrorLogs(const std::string &logFilename)
    {
      _logFile.open(logFilename, std::ios::app);

      if (!_logFile.is_open())
      {
        std::cerr << "Failed to open Log file" << std::endl;
      }
    }

    ~ErrorLogs()
    {
      if (_logFile.is_open())
      {
        _logFile.close();
      }
    }

    /**
     * @brief Increment physics-related error counter for a given mctruth value.
     * @param errCode Error code (should be in 300-399)
     * @param mctruth MC truth value
     */
    void countPhysicsError(ErrorCodes errCode, int mctruth)
    {
      if (static_cast<int>(errCode) >= 300 && static_cast<int>(errCode) <= 399)
        _physicsErrCountPerMctruth[mctruth][errCode]++;
    }

    void getErrLog(ErrorCodes &errCode, const std::string &additionalInfo = "", int mctruth = -1)
    {
      if (errCode == ErrorCodes::NO_ERROR)
        return; // No error, nothing to log
      else
      {
        std::string
            timestamp = _getTimestamp(),
            errorMessage = _getErrorMessage(errCode),
            concatMessage = "[" + timestamp + "] Error: " + errorMessage + " - " + additionalInfo;

        _errCount[errCode]++;

        Bool_t isPhysicsError = (static_cast<int>(errCode) >= 300 && static_cast<int>(errCode) <= 399);
        // Zliczaj błędy fizyczne dla każdego wywołania getErrLog jeśli mctruth podany
        if (mctruth != -1)
          countPhysicsError(errCode, mctruth);
        if (_printToScreen && !isPhysicsError)
          std::cerr << concatMessage << std::endl;

        if (_logFile.is_open() && !isPhysicsError)
        {
          _logFile << concatMessage << std::endl;
        }
      }
    };

    void setPrintToScreen(bool print) { _printToScreen = print; }
    bool getPrintToScreen() const { return _printToScreen; }

    void getLog(InfoCodes &infoCode, const std::string &additionalInfo = "")
    {
      std::string
          timestamp = _getTimestamp(),
          infoMessage = _getInfoMessage(infoCode),
          concatMessage = "[" + timestamp + "] Info: " + infoMessage + " - " + additionalInfo;

      if (_printToScreen)
        std::cerr << concatMessage << std::endl;

      if (_logFile.is_open())
      {
        _logFile << concatMessage << std::endl;
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
    void printPhysicsErrorStatsPerMctruth(bool printToScreen = false)
    {
      if (!_logFile.is_open())
        return;
      _logFile << "Physics-related error statistics per mctruth:\n";
      for (const auto &mctruth_pair : _physicsErrCountPerMctruth)
      {
        int mctruth = mctruth_pair.first;
        _logFile << "mctruth = " << mctruth << ":\n";
        if (printToScreen)
          std::cout << "mctruth = " << mctruth << ":\n";
        for (const auto &err_pair : mctruth_pair.second)
        {
          _logFile << "  " << _getErrorMessage(err_pair.first) << ": " << err_pair.second << " occurrences\n";
          if (printToScreen)
            std::cout << "  " << _getErrorMessage(err_pair.first) << ": " << err_pair.second << " occurrences\n";
        }
      }
      _logFile.flush();
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