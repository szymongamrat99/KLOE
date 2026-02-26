#ifndef SPLITFILEWRITER_H
#define SPLITFILEWRITER_H

#include <iostream>
#include <boost/filesystem.hpp>

#include <string>
#include <vector>
#include <map>
#include <fstream>
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <TDirectory.h>
#include <TRandom.h>
#include <kloe_class.h>
#include <MainMenu.h>

class SplitFileWriter
{
private:
  // Structure to track information about generated files
  struct FileInfo
  {
    std::string filename;
    Long64_t events;
    Double_t luminosity;
    Double_t sizeMB;
    std::map<std::string, Double_t> inputFiles; // input file -> luminosity contribution
  };

public:
  SplitFileWriter(const std::string &baseName,
                  Long64_t maxSizeBytes = 1500000000,
                  bool splitByRun = false,
                  const std::string &outputDir = "output",
                  const std::string &logFile = "file_weights.log",
                  Controls::FileType fileType = Controls::FileType::DATA,
                  bool singleFile = false,
                  Int_t currentFileNumber = 1);

  ~SplitFileWriter();

  void Fill(const std::map<std::string, Int_t> &intVars,
            const std::map<std::string, Float_t> &floatVars,
            const std::map<std::string, std::vector<Int_t>> &intArrays,
            const std::map<std::string, std::vector<Float_t>> &floatArrays);

  void Close();

  // Funkcja do aktualizacji luminozności z bieżącego pliku wejściowego
  void AddInputFileLuminosity(const std::string &inputFileName, Double_t luminosityPerEvent);

  // Ustawienie aktualnie przetwarzanego pliku wejściowego
  void SetCurrentInputFile(const std::string &inputFileName, Long64_t fileEvents, Double_t luminosityPerEvent);

  // Funkcja do przeliczania liczby zdarzeń na luminozność
  static Double_t CalculateLuminosity(Long64_t nEvents, Double_t luminosityPerEvent);

private:
  void OpenNewFile();
  void WriteAndClose();
  void EnsureOutputDirExists();
  void SetBranches();
  void LogFileInfo(const std::string &filename, Long64_t fileSize, Long64_t nEvents);
  std::string FileTypeToString(Controls::FileType fileType) const;
  void WriteFinalSummaryTable();
  Int_t GetMaxFileNumber() const;

  std::string _baseName;
  Long64_t _maxSizeBytes;
  bool _splitByRun;
  bool _singleFile;
  std::string _outputDir;
  std::string _logFile;
  Controls::FileType _fileType;
  Int_t _currentRun;

  TFile *_file;
  TTree *_tree;
  Int_t _fileCounter;
  Int_t _start;
  Long64_t _eventsInCurrentFile;
  std::ofstream _logStream;

  // Śledzenie luminozności
  std::string _currentInputFile;
  Long64_t _currentInputFileEvents;
  Double_t _currentLuminosityPerEvent;
  std::map<std::string, Double_t> _inputFilesLuminosity; // Mapa: nazwa pliku wejściowego -> jego wkład do bieżącego pliku wyjściowego [nb^-1]
  Double_t _totalLuminosityInCurrentFile;

  // Tracking all generated files for final summary
  std::vector<FileInfo> _generatedFiles;
  bool _summaryWritten;

  // Zmienne do zapisania
  std::map<std::string, Int_t> _intVars;                    // Zmienne Int_t
  std::map<std::string, Float_t> _floatVars;                // Zmienne Float_t
  std::map<std::string, std::vector<Int_t>> _intArrays;     // Tablice Int_t
  std::map<std::string, std::vector<Float_t>> _floatArrays; // Tablice Float_t
};

#endif // !SPLITFILEWRITER_H
