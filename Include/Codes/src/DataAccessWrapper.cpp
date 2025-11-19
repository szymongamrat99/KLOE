#include "DataAccessWrapper.h"
#include <TFile.h>
#include <TChainElement.h>
#include <TObjArray.h>
#include <iostream>
#include <algorithm>
#include <type_traits>
#include <regex>

using namespace KLOE;

// Definicja static const zmiennej
const Int_t DataAccessWrapper::kMaxArraySize;

DataAccessWrapper::DataAccessWrapper(TChain &chain, Bool_t useTTreeReader)
    : fChain(chain), fUseTTreeReader(useTTreeReader)
{

  // Regex do rozpoznania plików Data
  std::regex dataFileRegex(R"(prod2root_dk0[2]?_\d+(_v[12])?\.root$)");

  // Pobierz listę gałęzi MC z konfiguracji
  std::vector<TString> mcBranches = fVariableConfig.GetMCBranches();

  // Przejrzyj wszystkie pliki w chain'ie i wyłącz gałęzie MC dla plików Data
  TObjArray *fileList = fChain.GetListOfFiles();
  for (Int_t i = 0; i < fileList->GetEntries(); ++i)
  {
    TChainElement *element = (TChainElement *)fileList->At(i);
    TString filename = element->GetTitle();
    std::string filenameStr = filename.Data();

    // Sprawdź czy to plik Data
    if (std::regex_search(filenameStr, dataFileRegex))
    {
      // Wyłącz gałęzie MC dla tego pliku Data
      for (const auto &branchName : mcBranches)
      {
        fChain.SetBranchStatus(branchName, 0);
      }
      std::cout << "INFO: Disabled MC branches for Data file: " << filename << std::endl;
    }
  }

  // Inicjalizacja map dla v1/v2 arrays (SetBranchAddress)
  for (const auto &key : fVariableConfig.GetAllKeys())
  {
    const auto *info = fVariableConfig.GetVariableInfo(key);
    if (info && info->isArray)
    {
      if (info->type == "Int_t")
      {
        fIntArrays_v1[key].resize(kMaxArraySize, 0);
      }
      else if (info->type == "Float_t")
      {
        fFloatArrays_v1[key].resize(kMaxArraySize, 0.0f);
      }
      else if (info->type == "UInt_t")
      {
        fUIntArrays_v1[key].resize(kMaxArraySize, 0);
      }
    }
  }
}

DataAccessWrapper::~DataAccessWrapper()
{
  // Cleanup jest automatyczny dzięki unique_ptr i std::map
}

Bool_t DataAccessWrapper::Initialize()
{
  if (fChain.GetEntries() == 0)
  {
    std::cerr << "ERROR: DataAccessWrapper - Chain is empty!" << std::endl;
    return false;
  }

  // Mapuj wersje plików
  MapFileVersions();

  // Zawsze inicjalizuj SetBranchAddress (dla v1 i domyślnie v2)
  Bool_t branchOK = InitializeBranchAddress();
  if (!branchOK)
  {
    std::cerr << "ERROR: DataAccessWrapper - Failed to initialize SetBranchAddress" << std::endl;
    return false;
  }

  // Opcjonalnie inicjalizuj TTreeReader (tylko dla v2 jeśli włączone)
  Bool_t readerOK = true;
  if (fUseTTreeReader)
  {
    readerOK = InitializeTreeReader();
    if (!readerOK)
    {
      std::cerr << "ERROR: DataAccessWrapper - Failed to initialize TTreeReader" << std::endl;
      return false;
    }
  }

  fCurrentEntry = -1;

  std::cout << "INFO: DataAccessWrapper initialized successfully" << std::endl;
  std::cout << "INFO: Using " << (fUseTTreeReader ? "TTreeReader" : "SetBranchAddress")
            << " for v2 files" << std::endl;
  PrintFileTypeStats();

  return true;
}

Bool_t DataAccessWrapper::Next()
{
  if (fCurrentEntry + 1 >= fChain.GetEntries())
  {
    return false; // Koniec danych
  }

  fCurrentEntry++;

  // Sprawdź czy zmienił się plik i zaktualizuj wersję
  FileVersion previousVersion = fCurrentFileVersion;
  TString previousFile = GetCurrentFileName();

  // NAJPIERW załaduj entry do chain'a
  Int_t bytesRead = fChain.GetEntry(fCurrentEntry);
  if (bytesRead <= 0)
  {
    std::cerr << "ERROR: Failed to read entry " << fCurrentEntry << " from TChain" << std::endl;
    return false;
  }

  UpdateCurrentFileVersion();
  TString currentFile = GetCurrentFileName();

  // Sprawdź czy zmienił się plik lub wersja
  Bool_t fileChanged = (previousFile != currentFile);
  Bool_t versionChanged = (previousVersion != fCurrentFileVersion);

  if (fileChanged || versionChanged)
  {
    // std::cout << "INFO: File/Version change detected at entry " << fCurrentEntry << std::endl;
    // std::cout << "  Previous: " << previousFile.Data() << " ("
    //           << (previousVersion == FileVersion::V1 ? "v1" : previousVersion == FileVersion::V2 ? "v2" : "unknown") << ")" << std::endl;
    // std::cout << "  Current: " << currentFile.Data() << " ("
    //           << (fCurrentFileVersion == FileVersion::V1 ? "v1" : fCurrentFileVersion == FileVersion::V2 ? "v2" : "unknown") << ")" << std::endl;

    // REINICJALIZUJ TTreeReader tylko jeśli używamy go dla v2
    if (fUseTTreeReader && fCurrentFileVersion == FileVersion::V2)
    {
      // std::cout << "INFO: Reinitializing TTreeReader for v2 file..." << std::endl;
      if (!ReinitializeTreeReader())
      {
        std::cerr << "ERROR: Failed to reinitialize TTreeReader for entry " << fCurrentEntry << std::endl;
        return false;
      }
      // Ustaw TTreeReader na aktualny entry
      fReader->SetEntry(fCurrentEntry);
    }
  }

  // Aktualizuj statystyki w zależności od używanej metody dostępu
  if (fCurrentFileVersion == FileVersion::V2)
  {
    fV2Events++;
    if (fUseTTreeReader)
    {
      // Dla TTreeReader synchronizuj entry
      if (fReader)
      {
        fReader->SetEntry(fCurrentEntry);
        if (fReader->GetEntryStatus() != TTreeReader::kEntryValid)
        {
          std::cerr << "ERROR: TTreeReader entry " << fCurrentEntry << " is not valid. Status: "
                    << fReader->GetEntryStatus() << std::endl;
          return false;
        }
      }
      fTreeReaderEvents++;
    }
    else
    {
      fSetBranchEvents++;
    }
  }
  else if (fCurrentFileVersion == FileVersion::V1)
  {
    fV1Events++;
    fSetBranchEvents++;
  }

  return true;
}

DataAccessWrapper::FileVersion DataAccessWrapper::DetermineFileVersion(const TString &filename)
{
  if (filename.Contains("_v2.root"))
  {
    return FileVersion::V2;
  }
  else if (filename.Contains("_v1.root") || filename.EndsWith(".root"))
  {
    // Zakładamy, że pliki bez oznaczenia wersji to v1
    return FileVersion::V1;
  }
  return FileVersion::UNKNOWN;
}

void DataAccessWrapper::MapFileVersions()
{
  TObjArray *fileList = fChain.GetListOfFiles();

  for (Int_t i = 0; i < fileList->GetEntries(); ++i)
  {
    TChainElement *element = (TChainElement *)fileList->At(i);
    TString filename = element->GetTitle();

    FileVersion version = DetermineFileVersion(filename);
    fFileVersionMap[filename] = version;

    // std::cout << "INFO: File " << filename << " -> Version "
    //           << (version == FileVersion::V1 ? "v1" : version == FileVersion::V2 ? "v2"
    //                                                                              : "unknown")
    //           << std::endl;
  }
}

Bool_t DataAccessWrapper::InitializeTreeReader()
{
  if (!fUseTTreeReader)
  {
    std::cout << "INFO: TTreeReader disabled - using SetBranchAddress for all files" << std::endl;
    return true; // Sukces, ale nic nie robimy
  }

  try
  {
    fReader = std::make_unique<TTreeReader>(&fChain);

    // Inicjalizuj wszystkie zmienne skalarne i tablice dla v2
    for (const auto &key : fVariableConfig.GetAllKeys())
    {
      const auto *info = fVariableConfig.GetVariableInfo(key);
      if (!info)
        continue;

      TString branchName = info->nameV2;

      if (info->isArray)
      {
        // Tablice
        if (info->type == "Int_t")
        {
          fIntArrays_v2[key] = std::make_unique<TTreeReaderArray<Int_t>>(*fReader, branchName.Data());
        }
        else if (info->type == "Float_t")
        {
          fFloatArrays_v2[key] = std::make_unique<TTreeReaderArray<Float_t>>(*fReader, branchName.Data());
        }
        else if (info->type == "UInt_t")
        {
          fUIntArrays_v2[key] = std::make_unique<TTreeReaderArray<UInt_t>>(*fReader, branchName.Data());
        }
      }
      else
      {
        // Zmienne skalarne
        if (info->type == "Int_t")
        {
          fIntValues_v2[key] = std::make_unique<TTreeReaderValue<Int_t>>(*fReader, branchName.Data());
        }
        else if (info->type == "Float_t")
        {
          fFloatValues_v2[key] = std::make_unique<TTreeReaderValue<Float_t>>(*fReader, branchName.Data());
        }
        else if (info->type == "UInt_t")
        {
          fUIntValues_v2[key] = std::make_unique<TTreeReaderValue<UInt_t>>(*fReader, branchName.Data());
        }
      }
    }

    return true;
  }
  catch (const std::exception &e)
  {
    std::cerr << "ERROR: Failed to initialize TTreeReader: " << e.what() << std::endl;
    return false;
  }
}

Bool_t DataAccessWrapper::ReinitializeTreeReader()
{
  if (!fUseTTreeReader)
  {
    return true; // Sukces, ale nic nie robimy
  }

  try
  {
    // Wyczyść stare readery
    fIntValues_v2.clear();
    fFloatValues_v2.clear();
    fUIntValues_v2.clear();
    fIntArrays_v2.clear();
    fFloatArrays_v2.clear();
    fUIntArrays_v2.clear();

    // Stwórz nowy reader dla aktualnego pliku
    fReader.reset();
    fReader = std::make_unique<TTreeReader>(&fChain);

    // Reinicjalizuj wszystkie zmienne skalarne i tablice dla v2
    for (const auto &key : fVariableConfig.GetAllKeys())
    {
      const auto *info = fVariableConfig.GetVariableInfo(key);
      if (!info)
        continue;

      TString branchName = info->nameV2;

      try
      {
        if (info->isArray)
        {
          // Tablice
          if (info->type == "Int_t")
          {
            fIntArrays_v2[key] = std::make_unique<TTreeReaderArray<Int_t>>(*fReader, branchName.Data());
          }
          else if (info->type == "Float_t")
          {
            fFloatArrays_v2[key] = std::make_unique<TTreeReaderArray<Float_t>>(*fReader, branchName.Data());
          }
          else if (info->type == "UInt_t")
          {
            fUIntArrays_v2[key] = std::make_unique<TTreeReaderArray<UInt_t>>(*fReader, branchName.Data());
          }
        }
        else
        {
          // Zmienne skalarne
          if (info->type == "Int_t")
          {
            fIntValues_v2[key] = std::make_unique<TTreeReaderValue<Int_t>>(*fReader, branchName.Data());
          }
          else if (info->type == "Float_t")
          {
            fFloatValues_v2[key] = std::make_unique<TTreeReaderValue<Float_t>>(*fReader, branchName.Data());
          }
          else if (info->type == "UInt_t")
          {
            fUIntValues_v2[key] = std::make_unique<TTreeReaderValue<UInt_t>>(*fReader, branchName.Data());
          }
        }
      }
      catch (const std::exception &e)
      {
        std::cerr << "WARNING: Failed to reinitialize branch " << branchName.Data()
                  << " for key " << key.Data() << ": " << e.what() << std::endl;
        // Kontynuuj z następnym branch'em
      }
    }

    // std::cout << "INFO: TTreeReader successfully reinitialized" << std::endl;
    return true;
  }
  catch (const std::exception &e)
  {
    std::cerr << "ERROR: Failed to reinitialize TTreeReader: " << e.what() << std::endl;
    return false;
  }
}

Bool_t DataAccessWrapper::InitializeBranchAddress()
{
  try
  {
    // Inicjalizuj wszystkie zmienne dla v1 i v2 (gdy TTreeReader wyłączony)
    for (const auto &key : fVariableConfig.GetAllKeys())
    {
      const auto *info = fVariableConfig.GetVariableInfo(key);
      if (!info)
        continue;

      // Użyj nazw v1 dla plików v1, v2 dla plików v2
      TString branchNameV1 = info->nameV1;
      TString branchNameV2 = info->nameV2;

      if (info->isArray)
      {
        // Tablice - ustaw branch address na C arrays dla obu wersji
        if (info->type == "Int_t")
        {
          fChain.SetBranchAddress(branchNameV1.Data(), fIntArrays_v1[key].data());
          fChain.SetBranchAddress(branchNameV2.Data(), fIntArrays_v1[key].data());
        }
        else if (info->type == "Float_t")
        {
          fChain.SetBranchAddress(branchNameV1.Data(), fFloatArrays_v1[key].data());
          fChain.SetBranchAddress(branchNameV2.Data(), fFloatArrays_v1[key].data());
        }
        else if (info->type == "UInt_t")
        {
          fChain.SetBranchAddress(branchNameV1.Data(), fUIntArrays_v1[key].data());
          fChain.SetBranchAddress(branchNameV2.Data(), fUIntArrays_v1[key].data());
        }
      }
      else
      {
        // Zmienne skalarne dla obu wersji
        if (info->type == "Int_t")
        {
          fIntValues_v1[key] = 0;
          fChain.SetBranchAddress(branchNameV1.Data(), &fIntValues_v1[key]);
          fChain.SetBranchAddress(branchNameV2.Data(), &fIntValues_v1[key]);
        }
        else if (info->type == "Float_t")
        {
          fFloatValues_v1[key] = 0.0f;
          fChain.SetBranchAddress(branchNameV1.Data(), &fFloatValues_v1[key]);
          fChain.SetBranchAddress(branchNameV2.Data(), &fFloatValues_v1[key]);
        }
        else if (info->type == "UInt_t")
        {
          fUIntValues_v1[key] = 0;
          fChain.SetBranchAddress(branchNameV1.Data(), &fUIntValues_v1[key]);
          fChain.SetBranchAddress(branchNameV2.Data(), &fUIntValues_v1[key]);
        }
      }
    }

    return true;
  }
  catch (const std::exception &e)
  {
    std::cerr << "ERROR: Failed to initialize SetBranchAddress: " << e.what() << std::endl;
    return false;
  }
}

TString DataAccessWrapper::GetCurrentFileName() const
{
  if (fCurrentEntry < 0)
    return "";

  TFile *currentFile = fChain.GetFile();
  if (!currentFile)
    return "";

  return currentFile->GetName();
}

void DataAccessWrapper::UpdateCurrentFileVersion()
{
  TString currentFile = GetCurrentFileName();
  auto it = fFileVersionMap.find(currentFile);

  if (it != fFileVersionMap.end())
  {
    fCurrentFileVersion = it->second;
  }
  else
  {
    fCurrentFileVersion = DetermineFileVersion(currentFile);
    fFileVersionMap[currentFile] = fCurrentFileVersion;
  }
}

void DataAccessWrapper::PrintFileTypeStats() const
{
  std::cout << "\n=== DataAccessWrapper Statistics ===" << std::endl;
  std::cout << "Total files in chain: " << fFileVersionMap.size() << std::endl;

  Int_t v1Count = 0, v2Count = 0, unknownCount = 0;
  for (const auto &pair : fFileVersionMap)
  {
    switch (pair.second)
    {
    case FileVersion::V1:
      v1Count++;
      break;
    case FileVersion::V2:
      v2Count++;
      break;
    default:
      unknownCount++;
      break;
    }
  }

  std::cout << "v1.root files: " << v1Count << std::endl;
  std::cout << "v2.root files: " << v2Count << std::endl;
  std::cout << "Unknown files: " << unknownCount << std::endl;
  std::cout << "Events processed - v1: " << fV1Events << ", v2: " << fV2Events << std::endl;
  std::cout << "Access method - TTreeReader: " << fTreeReaderEvents
            << ", SetBranchAddress: " << fSetBranchEvents << std::endl;
  std::cout << "TTreeReader " << (fUseTTreeReader ? "ENABLED" : "DISABLED")
            << " for v2 files" << std::endl;
  std::cout << "======================================" << std::endl;
}

// ==================== TEMPLATE IMPLEMENTATIONS ====================

// Template specializations for ArrayConverter - C++14 compatible (musi być przed ArrayValueGetter)
namespace
{
  template <typename T>
  struct ArrayConverter
  {
    static const std::vector<T> &Convert(const DataAccessWrapper *wrapper, const TString &key, Int_t size)
    {
      static std::vector<T> empty;
      return empty; // Default implementation
    }
  };

  template <>
  struct ArrayConverter<Int_t>
  {
    static const std::vector<Int_t> &Convert(const DataAccessWrapper *wrapper, const TString &key, Int_t size)
    {
      size = std::min(size, wrapper->kMaxArraySize);
      auto &cache = wrapper->fIntVectorCache[key];
      auto it = wrapper->fIntArrays_v1.find(key);
      if (it != wrapper->fIntArrays_v1.end())
      {
        cache.clear();
        cache.reserve(size);
        for (Int_t i = 0; i < size; ++i)
        {
          cache.push_back(it->second[i]);
        }
      }
      return cache;
    }
  };

  template <>
  struct ArrayConverter<Float_t>
  {
    static const std::vector<Float_t> &Convert(const DataAccessWrapper *wrapper, const TString &key, Int_t size)
    {
      size = std::min(size, wrapper->kMaxArraySize);
      auto &cache = wrapper->fFloatVectorCache[key];
      auto it = wrapper->fFloatArrays_v1.find(key);
      if (it != wrapper->fFloatArrays_v1.end())
      {
        cache.clear();
        cache.reserve(size);
        for (Int_t i = 0; i < size; ++i)
        {
          cache.push_back(it->second[i]);
        }
      }
      return cache;
    }
  };

  template <>
  struct ArrayConverter<UInt_t>
  {
    static const std::vector<UInt_t> &Convert(const DataAccessWrapper *wrapper, const TString &key, Int_t size)
    {
      size = std::min(size, wrapper->kMaxArraySize);
      auto &cache = wrapper->fUIntVectorCache[key];
      auto it = wrapper->fUIntArrays_v1.find(key);
      if (it != wrapper->fUIntArrays_v1.end())
      {
        cache.clear();
        cache.reserve(size);
        for (Int_t i = 0; i < size; ++i)
        {
          cache.push_back(it->second[i]);
        }
      }
      return cache;
    }
  };
}

// Template specializations for GetScalarValue - C++14 compatible
namespace
{
  template <typename T>
  struct ScalarValueGetter
  {
    static T Get(const DataAccessWrapper *wrapper, const TString &key)
    {
      return T{}; // Default implementation
    }
  };

  template <>
  struct ScalarValueGetter<Int_t>
  {
    static Int_t Get(const DataAccessWrapper *wrapper, const TString &key)
    {
      if (wrapper->GetCurrentFileVersion() == DataAccessWrapper::FileVersion::V2 && wrapper->fUseTTreeReader)
      {
        // Użyj TTreeReader dla v2 gdy włączony
        auto it = wrapper->fIntValues_v2.find(key);
        if (it != wrapper->fIntValues_v2.end())
        {
          return **it->second;
        }
      }

      // Użyj SetBranchAddress dla v1 i domyślnie v2
      auto it = wrapper->fIntValues_v1.find(key);
      if (it != wrapper->fIntValues_v1.end())
      {
        return it->second;
      }

      return Int_t{};
    }
  };

  template <>
  struct ScalarValueGetter<Float_t>
  {
    static Float_t Get(const DataAccessWrapper *wrapper, const TString &key)
    {
      if (wrapper->GetCurrentFileVersion() == DataAccessWrapper::FileVersion::V2 && wrapper->fUseTTreeReader)
      {
        // Użyj TTreeReader dla v2 gdy włączony
        auto it = wrapper->fFloatValues_v2.find(key);
        if (it != wrapper->fFloatValues_v2.end())
        {
          return **it->second;
        }
      }

      // Użyj SetBranchAddress dla v1 i domyślnie v2
      auto it = wrapper->fFloatValues_v1.find(key);
      if (it != wrapper->fFloatValues_v1.end())
      {
        return it->second;
      }

      return Float_t{};
    }
  };

  template <>
  struct ScalarValueGetter<UInt_t>
  {
    static UInt_t Get(const DataAccessWrapper *wrapper, const TString &key)
    {
      if (wrapper->GetCurrentFileVersion() == DataAccessWrapper::FileVersion::V2 && wrapper->fUseTTreeReader)
      {
        // Użyj TTreeReader dla v2 gdy włączony
        auto it = wrapper->fUIntValues_v2.find(key);
        if (it != wrapper->fUIntValues_v2.end())
        {
          return **it->second;
        }
      }

      // Użyj SetBranchAddress dla v1 i domyślnie v2
      auto it = wrapper->fUIntValues_v1.find(key);
      if (it != wrapper->fUIntValues_v1.end())
      {
        return it->second;
      }

      return UInt_t{};
    }
  };
}

template <typename T>
T DataAccessWrapper::GetScalarValue(const TString &key) const
{
  return ScalarValueGetter<T>::Get(this, key);
}

// Template specializations for GetArrayValue - C++14 compatible
namespace
{
  template <typename T>
  struct ArrayValueGetter
  {
    static const std::vector<T> &Get(const DataAccessWrapper *wrapper, const TString &key)
    {
      static std::vector<T> empty;
      return empty; // Default implementation
    }
  };

  template <>
  struct ArrayValueGetter<Int_t>
  {
    static const std::vector<Int_t> &Get(const DataAccessWrapper *wrapper, const TString &key)
    {
      if (wrapper->GetCurrentFileVersion() == DataAccessWrapper::FileVersion::V2 && wrapper->fUseTTreeReader)
      {
        // Użyj TTreeReader dla v2 gdy włączony
        auto it = wrapper->fIntArrays_v2.find(key);
        if (it != wrapper->fIntArrays_v2.end())
        {
          auto &cache = wrapper->fIntVectorCache[key];
          cache.assign(it->second->begin(), it->second->end());
          return cache;
        }
      }

      // Użyj SetBranchAddress dla v1 i domyślnie v2
      const auto *info = wrapper->fVariableConfig.GetVariableInfo(key);
      if (info && info->sizeVariable.Length() > 0)
      {
        Int_t size = wrapper->GetScalarValue<Int_t>(info->sizeVariable);
        return ArrayConverter<Int_t>::Convert(wrapper, key, size);
      }

      static std::vector<Int_t> empty;
      return empty;
    }
  };

  template <>
  struct ArrayValueGetter<Float_t>
  {
    static const std::vector<Float_t> &Get(const DataAccessWrapper *wrapper, const TString &key)
    {
      if (wrapper->GetCurrentFileVersion() == DataAccessWrapper::FileVersion::V2 && wrapper->fUseTTreeReader)
      {
        // Użyj TTreeReader dla v2 gdy włączony
        auto it = wrapper->fFloatArrays_v2.find(key);
        if (it != wrapper->fFloatArrays_v2.end())
        {
          auto &cache = wrapper->fFloatVectorCache[key];
          cache.assign(it->second->begin(), it->second->end());
          return cache;
        }
      }

      // Użyj SetBranchAddress dla v1 i domyślnie v2
      const auto *info = wrapper->fVariableConfig.GetVariableInfo(key);
      if (info && info->sizeVariable.Length() > 0)
      {
        Int_t size = wrapper->GetScalarValue<Int_t>(info->sizeVariable);
        return ArrayConverter<Float_t>::Convert(wrapper, key, size);
      }

      static std::vector<Float_t> empty;
      return empty;
    }
  };

  template <>
  struct ArrayValueGetter<UInt_t>
  {
    static const std::vector<UInt_t> &Get(const DataAccessWrapper *wrapper, const TString &key)
    {
      if (wrapper->GetCurrentFileVersion() == DataAccessWrapper::FileVersion::V2 && wrapper->fUseTTreeReader)
      {
        // Użyj TTreeReader dla v2 gdy włączony
        auto it = wrapper->fUIntArrays_v2.find(key);
        if (it != wrapper->fUIntArrays_v2.end())
        {
          auto &cache = wrapper->fUIntVectorCache[key];
          cache.assign(it->second->begin(), it->second->end());
          return cache;
        }
      }

      // Użyj SetBranchAddress dla v1 i domyślnie v2
      const auto *info = wrapper->fVariableConfig.GetVariableInfo(key);
      if (info && info->sizeVariable.Length() > 0)
      {
        Int_t size = wrapper->GetScalarValue<Int_t>(info->sizeVariable);
        return ArrayConverter<UInt_t>::Convert(wrapper, key, size);
      }

      static std::vector<UInt_t> empty;
      return empty;
    }
  };
}

template <typename T>
const std::vector<T> &DataAccessWrapper::GetArrayValue(const TString &key) const
{
  return ArrayValueGetter<T>::Get(this, key);
}

template <typename T>
const std::vector<T> &DataAccessWrapper::ConvertArrayToVector(const TString &key, Int_t size) const
{
  return ArrayConverter<T>::Convert(this, key, size);
}

// ==================== PUBLIC GETTERS ====================

// Zmienne skalarne Int_t
Int_t DataAccessWrapper::GetNRun() const { return GetScalarValue<Int_t>("nrun"); }
Int_t DataAccessWrapper::GetNEv() const { return GetScalarValue<Int_t>("nev"); }
Int_t DataAccessWrapper::GetNClu() const { return GetScalarValue<Int_t>("nclu"); }
Int_t DataAccessWrapper::GetNCluMC() const { return GetScalarValue<Int_t>("nclumc"); }
Int_t DataAccessWrapper::GetNTCl() const { return GetScalarValue<Int_t>("ntcl"); }
Int_t DataAccessWrapper::GetNV() const { return GetScalarValue<Int_t>("nv"); }
Int_t DataAccessWrapper::GetNTV() const { return GetScalarValue<Int_t>("ntv"); }
Int_t DataAccessWrapper::GetNTMC() const { return GetScalarValue<Int_t>("ntmc"); }
Int_t DataAccessWrapper::GetNVtxMC() const { return GetScalarValue<Int_t>("nvtxmc"); }
Int_t DataAccessWrapper::GetEclFilfo() const { return GetScalarValue<Int_t>("eclfilfo"); }
Int_t DataAccessWrapper::GetEclFilfoWord() const { return GetScalarValue<Int_t>("eclfilfoword"); }
Int_t DataAccessWrapper::GetBunchNum() const { return GetScalarValue<Int_t>("bunchnum"); }
Int_t DataAccessWrapper::GetNECls() const { return GetScalarValue<Int_t>("necls"); }

// Zmienne skalarne Float_t
Float_t DataAccessWrapper::GetT0Step1() const { return GetScalarValue<Float_t>("t0step1"); }
Float_t DataAccessWrapper::GetBx() const { return GetScalarValue<Float_t>("bx"); }
Float_t DataAccessWrapper::GetBy() const { return GetScalarValue<Float_t>("by"); }
Float_t DataAccessWrapper::GetBz() const { return GetScalarValue<Float_t>("bz"); }
Float_t DataAccessWrapper::GetBxErr() const { return GetScalarValue<Float_t>("bxerr"); }
Float_t DataAccessWrapper::GetByErr() const { return GetScalarValue<Float_t>("byerr"); }
Float_t DataAccessWrapper::GetBzErr() const { return GetScalarValue<Float_t>("bzerr"); }
Float_t DataAccessWrapper::GetBlumx() const { return GetScalarValue<Float_t>("blumx"); }
Float_t DataAccessWrapper::GetBlumz() const { return GetScalarValue<Float_t>("blumz"); }
Float_t DataAccessWrapper::GetBpx() const { return GetScalarValue<Float_t>("bpx"); }
Float_t DataAccessWrapper::GetBpy() const { return GetScalarValue<Float_t>("bpy"); }
Float_t DataAccessWrapper::GetBpz() const { return GetScalarValue<Float_t>("bpz"); }
Float_t DataAccessWrapper::GetBpxErr() const { return GetScalarValue<Float_t>("bpxerr"); }
Float_t DataAccessWrapper::GetBpyErr() const { return GetScalarValue<Float_t>("bpyerr"); }
Float_t DataAccessWrapper::GetBpzErr() const { return GetScalarValue<Float_t>("bpzerr"); }
Float_t DataAccessWrapper::GetBRoots() const { return GetScalarValue<Float_t>("broots"); }
Float_t DataAccessWrapper::GetBRootsErr() const { return GetScalarValue<Float_t>("brootserr"); }

// Tablice Int_t
const std::vector<Int_t> &DataAccessWrapper::GetEclStream() const { return GetArrayValue<Int_t>("eclstream"); }
const std::vector<Int_t> &DataAccessWrapper::GetAssCl() const { return GetArrayValue<Int_t>("asscl"); }
const std::vector<Int_t> &DataAccessWrapper::GetIv() const { return GetArrayValue<Int_t>("iv"); }

// Tablice Int_t
const std::vector<Int_t> &DataAccessWrapper::GetVtxMC() const { return GetArrayValue<Int_t>("vtxmc"); }
const std::vector<Int_t> &DataAccessWrapper::GetPidMC() const { return GetArrayValue<Int_t>("pidmc"); }
const std::vector<Int_t> &DataAccessWrapper::GetMother() const { return GetArrayValue<Int_t>("mother"); }
const std::vector<Int_t> &DataAccessWrapper::GetKine() const { return GetArrayValue<Int_t>("kine"); }
const std::vector<Int_t> &DataAccessWrapper::GetKinMom() const { return GetArrayValue<Int_t>("kinmom"); }
const std::vector<Int_t> &DataAccessWrapper::GetPNum1() const { return GetArrayValue<Int_t>("pnum1"); }
const std::vector<Int_t> &DataAccessWrapper::GetPNum2() const { return GetArrayValue<Int_t>("pnum2"); }
const std::vector<Int_t> &DataAccessWrapper::GetPNum3() const { return GetArrayValue<Int_t>("pnum3"); }

// Tablice Float_t
const std::vector<Float_t> &DataAccessWrapper::GetXCl() const { return GetArrayValue<Float_t>("xcl"); }
const std::vector<Float_t> &DataAccessWrapper::GetYCl() const { return GetArrayValue<Float_t>("ycl"); }
const std::vector<Float_t> &DataAccessWrapper::GetZCl() const { return GetArrayValue<Float_t>("zcl"); }
const std::vector<Float_t> &DataAccessWrapper::GetTCl() const { return GetArrayValue<Float_t>("tcl"); }
const std::vector<Float_t> &DataAccessWrapper::GetEneCl() const { return GetArrayValue<Float_t>("enecl"); }
const std::vector<Float_t> &DataAccessWrapper::GetCurv() const { return GetArrayValue<Float_t>("curv"); }
const std::vector<Float_t> &DataAccessWrapper::GetPhiv() const { return GetArrayValue<Float_t>("phiv"); }
const std::vector<Float_t> &DataAccessWrapper::GetCotv() const { return GetArrayValue<Float_t>("cotv"); }
const std::vector<Float_t> &DataAccessWrapper::GetPxtv() const { return GetArrayValue<Float_t>("pxtv"); }
const std::vector<Float_t> &DataAccessWrapper::GetPytv() const { return GetArrayValue<Float_t>("pytv"); }
const std::vector<Float_t> &DataAccessWrapper::GetPztv() const { return GetArrayValue<Float_t>("pztv"); }
const std::vector<Float_t> &DataAccessWrapper::GetVtxCov1() const { return GetArrayValue<Float_t>("vtxcov1"); }
const std::vector<Float_t> &DataAccessWrapper::GetVtxCov2() const { return GetArrayValue<Float_t>("vtxcov2"); }
const std::vector<Float_t> &DataAccessWrapper::GetVtxCov3() const { return GetArrayValue<Float_t>("vtxcov3"); }
const std::vector<Float_t> &DataAccessWrapper::GetVtxCov4() const { return GetArrayValue<Float_t>("vtxcov4"); }
const std::vector<Float_t> &DataAccessWrapper::GetVtxCov5() const { return GetArrayValue<Float_t>("vtxcov5"); }
const std::vector<Float_t> &DataAccessWrapper::GetVtxCov6() const { return GetArrayValue<Float_t>("vtxcov6"); }
const std::vector<Float_t> &DataAccessWrapper::GetXv() const { return GetArrayValue<Float_t>("xv"); }
const std::vector<Float_t> &DataAccessWrapper::GetYv() const { return GetArrayValue<Float_t>("yv"); }
const std::vector<Float_t> &DataAccessWrapper::GetZv() const { return GetArrayValue<Float_t>("zv"); }
const std::vector<Float_t> &DataAccessWrapper::GetXvMC() const { return GetArrayValue<Float_t>("xvmc"); }
const std::vector<Float_t> &DataAccessWrapper::GetYvMC() const { return GetArrayValue<Float_t>("yvmc"); }
const std::vector<Float_t> &DataAccessWrapper::GetZvMC() const { return GetArrayValue<Float_t>("zvmc"); }
const std::vector<Float_t> &DataAccessWrapper::GetPxMC() const { return GetArrayValue<Float_t>("pxmc"); }
const std::vector<Float_t> &DataAccessWrapper::GetPyMC() const { return GetArrayValue<Float_t>("pymc"); }
const std::vector<Float_t> &DataAccessWrapper::GetPzMC() const { return GetArrayValue<Float_t>("pzmc"); }