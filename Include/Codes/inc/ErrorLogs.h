#ifndef ERROR_LOGS_H
#define ERROR_LOGS_H

#include <fstream>
#include <iostream>
#include <map>
#include <chrono>

#include <TString.h>

#include <const.h>

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

    DELTA_LT_ZERO = 200, /*!< Negative delta of a quadratic equation*/
    DENOM_EQ_ZERO = 201, /*!< Denominator of a fraction equal to zero*/
    DET_ZERO = 202,      /*!< Matrix determinant equal to zero*/
    NAN_VAL = 203,       /*!< NaN value in return*/
    CHI_SQR_STEP = 204,  /*!< Step of minimization too small*/

    NOT_RECOGNIZED = 666, /*!< Unexpected exception*/

    NUM_CODES = 11 /*!< Number of codes for loops*/
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

    void getErrLog(ErrorCodes &errCode, const std::string &additionalInfo = "")
    {
      std::string
          timestamp = _getTimestamp(),
          errorMessage = _getErrorMessage(errCode),
          concatMessage = "[" + timestamp + "] Error: " + errorMessage + " - " + additionalInfo;

      _errCount[errCode]++;
      std::cerr << concatMessage << std::endl;

      if (_logFile.is_open())
      {
        _logFile << concatMessage << std::endl;
      }
    };

    void getLog(InfoCodes &infoCode, const std::string &additionalInfo = "")
    {
      std::string
          timestamp = _getTimestamp(),
          infoMessage = _getInfoMessage(infoCode),
          concatMessage = "[" + timestamp + "] Info: " + infoMessage + " - " + additionalInfo;

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
  };
}

#endif