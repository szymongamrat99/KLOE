// SplitFileWriter.cpp
#include <SplitFileWriter.h>
#include <iomanip>
#include <ctime>

SplitFileWriter::SplitFileWriter(const std::string &baseName,
                                 Long64_t maxSizeBytes,
                                 bool splitByRun,
                                 const std::string &outputDir,
                                 const std::string &logFile,
                                 Controls::FileType fileType,
                                 bool singleFile)
    : _baseName(baseName),
      _maxSizeBytes(maxSizeBytes),
      _splitByRun(splitByRun),
      _singleFile(singleFile),
      _outputDir(outputDir),
      _logFile(outputDir + "/" + logFile),
      _fileType(fileType),
      _currentRun(-1),
      _file(nullptr),
      _tree(nullptr),
      _fileCounter(1),
      _start(1),
      _eventsInCurrentFile(0),
      _currentInputFile(""),
      _currentInputFileEvents(0),
      _currentLuminosityPerEvent(0.0),
      _totalLuminosityInCurrentFile(0.0),
      _summaryWritten(false)
{
  EnsureOutputDirExists();
  
  // Jeśli włączony jest tryb single file, ustaw _fileCounter na max istniejący numer + 1
  if (_singleFile)
  {
    _fileCounter = GetMaxFileNumber() + 1;
  }
  
  _logStream.open(_logFile, std::ios::out | std::ios::app);
  if (_logStream.is_open())
  {
    auto t = std::time(nullptr);
    auto tm = *std::localtime(&t);
    std::string fileTypeStr = FileTypeToString(_fileType);
    _logStream << "\n=== Split File Writer Log (" << fileTypeStr << ") ===" << std::endl;
    _logStream << "Start time: " << std::put_time(&tm, "%Y-%m-%d %H:%M:%S") << std::endl;
    _logStream << "Base name: " << _baseName << std::endl;
    _logStream << "File type: " << fileTypeStr << std::endl;
    // _logStream << "Max file size: " << _maxSizeBytes / (1024.0 * 1024.0) << " MB" << std::endl;
    _logStream << "Split by run: " << (_splitByRun ? "Yes" : "No") << std::endl;
    _logStream << "=============================\n" << std::endl;
  }
  
  OpenNewFile();
}

SplitFileWriter::~SplitFileWriter()
{
  Close();
  WriteFinalSummaryTable();
  if (_logStream.is_open())
  {
    _logStream.close();
  }
}

void SplitFileWriter::EnsureOutputDirExists()
{
  boost::filesystem::create_directories(_outputDir);
}

std::string SplitFileWriter::FileTypeToString(Controls::FileType fileType) const
{
  switch (fileType)
  {
  case Controls::FileType::DATA:
    return "DATA";
  case Controls::FileType::ALL_PHYS:
    return "ALL_PHYS";
  case Controls::FileType::ALL_PHYS2:
    return "ALL_PHYS2";
  case Controls::FileType::ALL_PHYS3:
    return "ALL_PHYS3";
  default:
    return "UNKNOWN";
  }
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

  std::string fileTypeStr = FileTypeToString(_fileType);
  std::string filename = _outputDir + "/" + _baseName + "_" + std::to_string(_fileCounter++) + ".root";
  _file = new TFile(filename.c_str(), "RECREATE");
  _tree = new TTree("h1", "Split data");
  _eventsInCurrentFile = 0;
  _totalLuminosityInCurrentFile = 0.0;
  _inputFilesLuminosity.clear();

  // Auto-save co ~10MB
  _tree->SetAutoSave(10000000);
}

Double_t SplitFileWriter::CalculateLuminosity(Long64_t nEvents, Double_t luminosityPerEvent)
{
  // Przeliczenie: Luminozność [nb^-1] = liczba zdarzeń × współczynnik [nb^-1/event]
  return nEvents * luminosityPerEvent + 1.38; // stała przesunięcia
}

void SplitFileWriter::SetCurrentInputFile(const std::string &inputFileName, Long64_t fileEvents, Double_t luminosityPerEvent)
{
  _currentInputFile = inputFileName;
  _currentInputFileEvents = fileEvents;
  _currentLuminosityPerEvent = luminosityPerEvent;
}

void SplitFileWriter::AddInputFileLuminosity(const std::string &inputFileName, Double_t luminosityPerEvent)
{
  // Ta funkcja jest wywoływana przy każdym Fill()
  // Dodaje proporcjonalny wkład luminozności z bieżącego pliku wejściowego
  
  Double_t lumiContribution = CalculateLuminosity(_currentInputFileEvents, luminosityPerEvent);
  
  // Jeśli to nowy plik wejściowy, dodaj go do mapy
  if (_inputFilesLuminosity.find(inputFileName) == _inputFilesLuminosity.end())
  {
    _inputFilesLuminosity[inputFileName] = lumiContribution;
    _totalLuminosityInCurrentFile += lumiContribution;
  }
}

void SplitFileWriter::LogFileInfo(const std::string &filename, Long64_t fileSize, Long64_t nEvents)
{
  if (!_logStream.is_open())
    return;

  boost::filesystem::path p(filename);
  std::string baseName = p.filename().string();

  Double_t sizeMB = fileSize / (1024.0 * 1024.0);
  Double_t sizeGB = fileSize / (1024.0 * 1024.0 * 1024.0);

  std::string fileTypeStr = FileTypeToString(_fileType);

  _logStream << std::fixed << std::setprecision(2);
  _logStream << "File closed (" << fileTypeStr << "): " << baseName << std::endl;
  _logStream << "  Size: " << sizeMB << " MB (" << sizeGB << " GB)" << std::endl;
  _logStream << "  Events: " << nEvents << std::endl;
  _logStream << "  Total luminosity: " << std::setprecision(4) 
             << _totalLuminosityInCurrentFile << " nb^-1" << std::endl;
  
  _logStream << "  Input files contribution:" << std::endl;
  for (const auto &entry : _inputFilesLuminosity)
  {
    boost::filesystem::path inputPath(entry.first);
    _logStream << "    - " << inputPath.filename().string() 
               << ": " << std::setprecision(4) << entry.second << " nb^-1" << std::endl;
  }
  _logStream << std::endl;
}

void SplitFileWriter::WriteFinalSummaryTable()
{
  if (!_logStream.is_open() || _summaryWritten)
    return;

  _summaryWritten = true;

  auto t = std::time(nullptr);
  auto tm = *std::localtime(&t);
  std::string fileTypeStr = FileTypeToString(_fileType);

  _logStream << "\n===============================================" << std::endl;
  _logStream << "FINAL SUMMARY TABLE" << std::endl;
  _logStream << "Analysis completed: " << std::put_time(&tm, "%Y-%m-%d %H:%M:%S") << std::endl;
  _logStream << "===============================================" << std::endl;

  // Parse base filename to extract analysis parameters
  std::string baseName = _baseName;
  
  // Extract file type
  _logStream << "File Type: " << fileTypeStr << std::endl;
  
  // Try to extract other parameters from base filename
  // Example: mk0_initial_analysis_all_phys3_SIGNAL_MIXED_Signal
  if (baseName.find("SIGNAL") != std::string::npos) {
    _logStream << "Hypothesis: SIGNAL" << std::endl;
  } else if (baseName.find("FOUR_PI") != std::string::npos) {
    _logStream << "Hypothesis: FOUR_PI" << std::endl;
  } else if (baseName.find("OMEGAPI") != std::string::npos) {
    _logStream << "Hypothesis: OMEGAPI" << std::endl;
  } else {
    _logStream << "Hypothesis: Unknown" << std::endl;
  }
  
  if (baseName.find("MIXED") != std::string::npos) {
    _logStream << "Smearing: MIXED" << std::endl;
  } else if (baseName.find("NoSmearing") != std::string::npos) {
    _logStream << "Smearing: No Smearing" << std::endl;
  } else {
    _logStream << "Smearing: Unknown" << std::endl;
  }
  
  // Check for signal-only generation and extract channel
  if (baseName.find("Signal") != std::string::npos) {
    _logStream << "Generation: Signal Only" << std::endl;
    
    // Extract signal channel - it's the last part after the last underscore
    std::string signalChannel = "Unknown";
    if (baseName.find("_Signal") != std::string::npos) {
      signalChannel = "Signal";
    } else if (baseName.find("_Regeneration") != std::string::npos) {
      signalChannel = "Regeneration";
    } else if (baseName.find("_Omega") != std::string::npos) {
      signalChannel = "Omega";
    } else if (baseName.find("_3pi0") != std::string::npos) {
      signalChannel = "3pi0";
    } else if (baseName.find("_Semileptonic") != std::string::npos) {
      signalChannel = "Semileptonic";
    } else if (baseName.find("_Other") != std::string::npos) {
      signalChannel = "Other";
    } else if (baseName.find("_pi+pi-pi+pi-") != std::string::npos) {
      signalChannel = "pi+pi-pi+pi-";
    }
    _logStream << "Signal Channel: " << signalChannel << std::endl;
  } else {
    _logStream << "Generation: All" << std::endl;
    _logStream << "Signal Channel: All" << std::endl;
  }
  
  _logStream << "===============================================" << std::endl;

  if (_generatedFiles.empty())
  {
    _logStream << "No files were generated." << std::endl;
    _logStream << "===============================================\n" << std::endl;
    return;
  }

  // Calculate totals
  Long64_t totalEvents = 0;
  Double_t totalLuminosity = 0.0;
  Double_t totalSizeMB = 0.0;

  for (const auto &fileInfo : _generatedFiles)
  {
    totalEvents += fileInfo.events;
    totalLuminosity += fileInfo.luminosity;
    totalSizeMB += fileInfo.sizeMB;
  }

  // Write table header
  _logStream << std::left;
  _logStream << std::setw(10) << "File #"
             << std::setw(12) << "Events"
             << std::setw(15) << "Luminosity"
             << std::setw(12) << "Size (MB)" << std::endl;
  _logStream << std::string(49, '-') << std::endl;

  // Write file information
  for (size_t i = 0; i < _generatedFiles.size(); ++i)
  {
    const auto &fileInfo = _generatedFiles[i];
    
    // Extract file number from filename
    std::string fileNumber = std::to_string(i + 1);
    
    _logStream << std::setw(10) << fileNumber
               << std::setw(12) << fileInfo.events
               << std::setw(15) << std::fixed << std::setprecision(4) << fileInfo.luminosity
               << std::setw(12) << std::fixed << std::setprecision(2) << fileInfo.sizeMB << std::endl;
  }

  _logStream << std::string(49, '-') << std::endl;
  _logStream << std::setw(10) << "TOTAL:"
             << std::setw(12) << totalEvents
             << std::setw(15) << std::fixed << std::setprecision(4) << totalLuminosity
             << std::setw(12) << std::fixed << std::setprecision(2) << totalSizeMB << std::endl;

  // Write detailed input file contributions
  _logStream << "\nDetailed Input File Contributions:" << std::endl;
  _logStream << std::string(49, '=') << std::endl;

  for (size_t i = 0; i < _generatedFiles.size(); ++i)
  {
    const auto &fileInfo = _generatedFiles[i];
    _logStream << "File " << (i + 1) << ":" << std::endl;
    
    if (fileInfo.inputFiles.empty())
    {
      _logStream << "  No input file contributions recorded." << std::endl;
    }
    else
    {
      for (const auto &input : fileInfo.inputFiles)
      {
        boost::filesystem::path inputPath(input.first);
        _logStream << "  - " << inputPath.filename().string()
                   << ": " << std::fixed << std::setprecision(4) << input.second << " nb^-1" << std::endl;
      }
    }
    _logStream << std::endl;
  }

  _logStream << "===============================================\n" << std::endl;
}

void SplitFileWriter::WriteAndClose()
{
  if (_file)
  {
    _tree->Write();
    
    std::string filename = _file->GetName();
    Long64_t fileSize = _file->GetSize();
    Long64_t nEvents = _eventsInCurrentFile;
    
    // Store file information for final summary
    FileInfo fileInfo;
    fileInfo.filename = boost::filesystem::path(filename).filename().string();
    fileInfo.events = nEvents;
    fileInfo.luminosity = _totalLuminosityInCurrentFile;
    fileInfo.sizeMB = fileSize / (1024.0 * 1024.0);
    fileInfo.inputFiles = _inputFilesLuminosity;
    _generatedFiles.push_back(fileInfo);
    
    _file->Close();
    
    LogFileInfo(filename, fileSize, nEvents);
    
    delete _file;
    _file = nullptr;
    _tree = nullptr;
  }
}

Int_t SplitFileWriter::GetMaxFileNumber() const
{
  Int_t maxNumber = 0;
  
  try
  {
    boost::filesystem::path outputPath(_outputDir);
    if (!boost::filesystem::exists(outputPath) || !boost::filesystem::is_directory(outputPath))
    {
      return 0; // Brak katalogu lub nie jest katalogiem
    }
    
    // Szukamy plików pasujących do wzorca: baseName_<number>.root
    std::string pattern = _baseName + "_";
    
    for (const auto &entry : boost::filesystem::directory_iterator(outputPath))
    {
      if (boost::filesystem::is_regular_file(entry.path()))
      {
        std::string filename = entry.path().filename().string();
        
        // Sprawdzamy czy plik pasuje do wzorca
        if (filename.find(pattern) == 0 && filename.find(".root") != std::string::npos)
        {
          // Wyodrębniamy numer z nazwy pliku
          // Format: baseName_<number>.root
          size_t startPos = pattern.length();
          size_t endPos = filename.find(".root", startPos);
          
          if (endPos != std::string::npos)
          {
            std::string numberStr = filename.substr(startPos, endPos - startPos);
            try
            {
              Int_t number = std::stoi(numberStr);
              if (number > maxNumber)
              {
                maxNumber = number;
              }
            }
            catch (const std::exception &)
            {
              // Ignore files that don't have valid number format
            }
          }
        }
      }
    }
  }
  catch (const std::exception &e)
  {
    std::cerr << "Error in GetMaxFileNumber(): " << e.what() << std::endl;
  }
  
  return maxNumber;
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

  // Dodaj wkład luminozności z bieżącego pliku wejściowego (jeśli ustawiony)
  if (!_currentInputFile.empty() && _currentLuminosityPerEvent > 0)
  {
    AddInputFileLuminosity(_currentInputFile, _currentLuminosityPerEvent);
  }

  // Logika split - pomijana gdy tryb single file jest aktywny
  if (!_singleFile)
  {
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
  }

  if (_start == 1)
  {
    SetBranches();
    _start = 0;
  }

  _currentRun = intVars.at("nrun");
  _tree->Fill();
  _eventsInCurrentFile++;
}

void SplitFileWriter::Close()
{
  WriteAndClose();
  WriteFinalSummaryTable();
}