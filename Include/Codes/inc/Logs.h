#ifndef LOGS_H
#define LOGS_H

#include <fstream>
#include <iostream>
#include <map>
#include <ctime>

#include <TString.h>

namespace LogsHandling
{
  enum class GeneralLogs
  {
    TIME_STAMP,
    USED_FILES,
    USED_DATA_TYPE,
    GENERAL_SETTINGS,
    OPERATION_TIME,

    NEUTREC_KIN_FIT_SETTINGS,

    NUM
  };

  class Logs
  {
  private:
    GeneralLogs generalLogsList[int(GeneralLogs::NUM)];
    std::map<GeneralLogs, TString> generalLogs;
    std::map<GeneralLogs, int> errCount;
    const time_t timestamp = time(0);
    tm *now = localtime(&timestamp); 

  public:
    Logs()
    {
      // Creation of iterable list of general logs
      generalLogsList[0] = GeneralLogs::TIME_STAMP;
      generalLogsList[1] = GeneralLogs::USED_FILES;
      generalLogsList[2] = GeneralLogs::USED_DATA_TYPE;
      generalLogsList[3] = GeneralLogs::GENERAL_SETTINGS;
      generalLogsList[4] = GeneralLogs::OPERATION_TIME;
      generalLogsList[5] = GeneralLogs::NEUTREC_KIN_FIT_SETTINGS;

      // General Logs
      generalLogs[GeneralLogs::TIME_STAMP] = "Time of init: " + std::to_string(now->tm_year + 1900) + "-" + std::to_string(now->tm_mon + 1) + "-" + std::to_string(now->tm_mday) + " " + std::to_string(now->tm_hour) + ":" + std::to_string(now->tm_min) + ":" + std::to_string(now->tm_sec);
      generalLogs[GeneralLogs::USED_FILES] = "Files range used for the generation: ";
      generalLogs[GeneralLogs::USED_DATA_TYPE] = "Data type used for the generation: ";
      generalLogs[GeneralLogs::GENERAL_SETTINGS] = "General settings of the program: ";
      generalLogs[GeneralLogs::OPERATION_TIME] = "Operation times of the program: ";
      generalLogs[GeneralLogs::NEUTREC_KIN_FIT_SETTINGS] = "Neutrec kin fit settings: ";

    }

    void setErrCount(GeneralLogs logsCode)
    {
      errCount[logsCode]++;
    }

    void getErrLog(GeneralLogs logsCode)
    {
      std::cerr << "--------------------------------------------------" << std::endl;
      std::cerr << generalLogs[logsCode] << std::endl;
      std::cerr << "--------------------------------------------------" << std::endl;
    };

    void getErrLog(GeneralLogs logsCode, std::ofstream& LogFile)
    {
      LogFile << "--------------------------------------------------" << std::endl;
      LogFile << generalLogs[logsCode] << std::endl;
      LogFile << "--------------------------------------------------" << std::endl;
    };

    void errLogStats()
    {
      std::cout << "Logs of program operation:" << std::endl;
      for(int i = 0; i < int(GeneralLogs::NUM); i++)
      {
        std::cout << "Error code " << int(generalLogsList[i]) << ": " << errCount[generalLogsList[i]] << " exceptions" << std::endl;
      }
    };

    void errLogStats(std::ofstream& LogFile)
    {
      LogFile << "Logs of program operation:" << std::endl;
      for(int i = 0; i < int(GeneralLogs::NUM); i++)
      {
        LogFile << "--------------------------------------------------" << std::endl;
        LogFile << "Error code " << int(generalLogsList[i]) << ": " << errCount[generalLogsList[i]] << " exceptions" << std::endl;
        LogFile << "--------------------------------------------------" << std::endl;
      }
    };
  };
}

#endif