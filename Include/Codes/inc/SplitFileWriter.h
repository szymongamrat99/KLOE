#ifndef SPLITFILEWRITER_H
#define SPLITFILEWRITER_H

#include <iostream>
#include <boost/filesystem.hpp>

#include <string>
#include <vector>
#include <map>
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <TDirectory.h>
#include <TRandom.h>

class SplitFileWriter
{
public:
  SplitFileWriter(const std::string &baseName,
                  Long64_t maxSizeBytes = 1500000000,
                  bool splitByRun = false,
                  const std::string &outputDir = "output");

  ~SplitFileWriter();

  void Fill(const std::map<std::string, Int_t> &intVars,
            const std::map<std::string, Float_t> &floatVars,
            const std::map<std::string, std::vector<Int_t>> &intArrays,
            const std::map<std::string, std::vector<Float_t>> &floatArrays);
  void Close();

private:
  void OpenNewFile();
  void WriteAndClose();
  void EnsureOutputDirExists();
  void SetBranches();

  std::string _baseName;
  Long64_t _maxSizeBytes;
  bool _splitByRun;
  std::string _outputDir;
  Int_t _currentRun;

  TFile *_file;
  TTree *_tree;
  Int_t _fileCounter;
  Int_t _start;

  // Zmienne do zapisania
  std::map<std::string, Int_t> _intVars;         // Zmienne Int_t
  std::map<std::string, Float_t> _floatVars;     // Zmienne Float_t
  std::map<std::string, std::vector<Int_t>> _intArrays;     // Zmienne Int_t
  std::map<std::string, std::vector<Float_t>> _floatArrays; // Tablice Float_t

};

#endif // !SPLITFILEWRITER_H
