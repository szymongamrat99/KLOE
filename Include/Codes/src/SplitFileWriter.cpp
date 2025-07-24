// SplitFileWriter.cpp
#include <SplitFileWriter.h>

SplitFileWriter::SplitFileWriter(const std::string &baseName,
                                 Long64_t maxSizeBytes,
                                 bool splitByRun,
                                 const std::string &outputDir)
    : _baseName(baseName),
      _maxSizeBytes(maxSizeBytes),
      _splitByRun(splitByRun),
      _outputDir(outputDir),
      _currentRun(-1),
      _file(nullptr),
      _tree(nullptr),
      _fileCounter(1),
      _start(1)
{
  EnsureOutputDirExists();
  OpenNewFile();
}

SplitFileWriter::~SplitFileWriter()
{
  Close();
}

void SplitFileWriter::EnsureOutputDirExists()
{
  boost::filesystem::create_directories(_outputDir);
}

void SplitFileWriter::SetBranches()
{
  // Dodanie zmiennych typu Int_t (z konkretnymi nazwami)
  for (auto &var : _intVars)
  {
    _tree->Branch(var.first.c_str(), &(var.second), (var.first + "/I").c_str());
  }

  // Dodanie zmiennych typu Float_t (z konkretnymi nazwami)
  for (auto &var : _floatVars)
  {
    _tree->Branch(var.first.c_str(), &(var.second), (var.first + "/F").c_str());
  }

  // Dodanie tablic (z konkretnymi nazwami)
  for (auto &arr : _intArrays)
  {
    _tree->Branch(arr.first.c_str(), &(arr.second));
  }

  // Dodanie tablic (z konkretnymi nazwami)
  for (auto &arr : _floatArrays)
  {
    _tree->Branch(arr.first.c_str(), &(arr.second));
  }

}

void SplitFileWriter::OpenNewFile()
{
  if (_file)
    _file->Close();

  std::string filename = _outputDir + "/" + _baseName + "_" + std::to_string(_fileCounter++) + ".root";
  _file = new TFile(filename.c_str(), "RECREATE");
  _tree = new TTree("h1", "Split data");

  // Auto-save co ~10MB
  _tree->SetAutoSave(10000000);
}

void SplitFileWriter::WriteAndClose()
{
  if (_file)
  {
    _tree->Write();
    _file->Close();
    delete _file;
    _file = nullptr;
    _tree = nullptr;
  }
}

void SplitFileWriter::Fill(const std::map<std::string, Int_t> &intVars,
                           const std::map<std::string, Float_t> &floatVars,
                           const std::map<std::string, std::vector<Int_t>> &intArrays,
                           const std::map<std::string, std::vector<Float_t>> &floatArrays)
{
  // For integer variables
  for (const auto &pair : intVars)
  {
    _intVars[pair.first] = pair.second;
  }

  // For float variables
  for (const auto &pair : floatVars)
  {
    _floatVars[pair.first] = pair.second;
  }

  // For integer arrays
  for (const auto &pair : intArrays)
  {
    // If the vector already exists, update its contents; otherwise, emplace a new one
    auto it = _intArrays.find(pair.first);
    if (it != _intArrays.end())
      it->second = pair.second;
    else
      _intArrays.emplace(pair.first, pair.second);
  }

  // For float arrays
  for (const auto &pair : floatArrays)
  {
    auto it = _floatArrays.find(pair.first);
    if (it != _floatArrays.end())
      it->second = pair.second;
    else
      _floatArrays.emplace(pair.first, pair.second);
  }

  if (_splitByRun && _currentRun != -1 && intVars.at("nrun") != _currentRun)
  {
    WriteAndClose();
    OpenNewFile();
    SetBranches();
  }
  else if (!_splitByRun && _file->GetSize() >= _maxSizeBytes)
  {
    WriteAndClose();
    OpenNewFile();
    SetBranches();
  }

  if (_start == 1)
  {
    SetBranches();
    _start = 0;
  }

  _currentRun = intVars.at("nrun");
  _tree->Fill();
}

void SplitFileWriter::Close()
{
  WriteAndClose();
}