#include "VariableConfig.h"
#include <iostream>

using namespace KLOE;

VariableConfig::VariableConfig()
{
  InitializeDefaultMappings();
}

void VariableConfig::AddVariable(const TString &key, const VariableInfo &info)
{
  fVariableMap[key] = info;
}

const VariableConfig::VariableInfo *VariableConfig::GetVariableInfo(const TString &key) const
{
  auto it = fVariableMap.find(key);
  if (it != fVariableMap.end())
  {
    return &(it->second);
  }
  return nullptr;
}

TString VariableConfig::GetVariableName(const TString &key, Bool_t isV2) const
{
  const VariableInfo *info = GetVariableInfo(key);
  if (info)
  {
    return isV2 ? info->nameV2 : info->nameV1;
  }
  return "";
}

Bool_t VariableConfig::HasVariable(const TString &key) const
{
  return fVariableMap.find(key) != fVariableMap.end();
}

std::vector<TString> VariableConfig::GetAllKeys() const
{
  std::vector<TString> keys;
  keys.reserve(fVariableMap.size());

  for (const auto &pair : fVariableMap)
  {
    keys.push_back(pair.first);
  }

  return keys;
}

void VariableConfig::PrintConfiguration() const
{
  std::cout << "\n=== Variable Configuration ===" << std::endl;
  std::cout << "Total variables: " << fVariableMap.size() << std::endl;

  for (const auto &pair : fVariableMap)
  {
    const TString &key = pair.first;
    const VariableInfo &info = pair.second;

    std::cout << "Key: " << key.Data() << std::endl;
    std::cout << "  v1: " << info.nameV1.Data() << " -> v2: " << info.nameV2.Data() << std::endl;
    std::cout << "  Type: " << info.type.Data();
    if (info.isArray)
    {
      std::cout << "[] (size: " << info.sizeVariable.Data() << ")";
    }
    std::cout << std::endl
              << std::endl;
  }
  std::cout << "==============================" << std::endl;
}

void VariableConfig::InitializeDefaultMappings()
{
  // Skalarne zmienne Int_t
  AddVariable("nrun", VariableInfo("nRun", "nRun", "Int_t"));
  AddVariable("nev", VariableInfo("nEv", "nEv", "Int_t"));
  AddVariable("nclu", VariableInfo("nClu", "nClu", "Int_t"));
  AddVariable("ntcl", VariableInfo("nTcl", "nTcl", "Int_t"));
  AddVariable("nv", VariableInfo("nV", "nV", "Int_t"));
  AddVariable("ntv", VariableInfo("nTv", "nTv", "Int_t"));
  AddVariable("ntmc", VariableInfo("nTMC", "nTMC", "Int_t"));
  AddVariable("nvtxmc", VariableInfo("nVtxMC", "nVtxMC", "Int_t"));
  AddVariable("eclfilfo", VariableInfo("EclFilfo", "EclFilfo", "Int_t"));
  // AddVariable("eclfilfoword", VariableInfo("EclFilfoword", "EclFilfoWord", "Int_t"));
  AddVariable("necls", VariableInfo("NEcls", "NEcls", "Int_t"));

  // Skalarne zmienne Float_t
  AddVariable("t0step1", VariableInfo("T0Step1", "T0Step1", "Float_t"));
  AddVariable("bx", VariableInfo("Bx", "Bx", "Float_t"));
  AddVariable("by", VariableInfo("By", "By", "Float_t"));
  AddVariable("bz", VariableInfo("Bz", "Bz", "Float_t"));
  AddVariable("bxerr", VariableInfo("BSx", "BSx", "Float_t"));
  AddVariable("byerr", VariableInfo("BSy", "BSy", "Float_t"));
  AddVariable("bzerr", VariableInfo("BSz", "BSz", "Float_t"));
  AddVariable("bpx", VariableInfo("BPx", "BPx", "Float_t"));
  AddVariable("bpy", VariableInfo("BPy", "BPy", "Float_t"));
  AddVariable("bpz", VariableInfo("BPz", "BPz", "Float_t"));
  AddVariable("bpxerr", VariableInfo("BWidPx", "BWidPx", "Float_t"));
  AddVariable("bpyerr", VariableInfo("BWidPy", "BWidPy", "Float_t"));
  AddVariable("bpzerr", VariableInfo("BWidPz", "BWidPz", "Float_t"));
  AddVariable("broots", VariableInfo("Broots", "Broots", "Float_t"));
  AddVariable("brootserr", VariableInfo("BrootsErr", "BrootsErr", "Float_t"));

  // Tablice Int_t
  AddVariable("eclstream", VariableInfo("EclStream", "EclStream", "Int_t", true, "necls"));
  AddVariable("asscl", VariableInfo("AssCl", "AssCl", "Int_t", true, "nclu"));
  AddVariable("iv", VariableInfo("iV", "iV", "Int_t", true, "ntv"));

  // Tablice UInt_t (MC)
  AddVariable("vtxmc", VariableInfo("VtxMC", "VtxMC", "Int_t", true, "ntmc"));
  AddVariable("pidmc", VariableInfo("PidMC", "PidMC", "Int_t", true, "ntmc"));
  AddVariable("mother", VariableInfo("Mother", "Mother", "Int_t", true, "ntmc"));

  // Tablice Float_t - clustery
  AddVariable("xcl", VariableInfo("XCl", "XCl", "Float_t", true, "nclu"));
  AddVariable("ycl", VariableInfo("YCl", "YCl", "Float_t", true, "nclu"));
  AddVariable("zcl", VariableInfo("ZCl", "ZCl", "Float_t", true, "nclu"));
  AddVariable("tcl", VariableInfo("TCl", "TCl", "Float_t", true, "ntcl"));
  AddVariable("enecl", VariableInfo("EneCl", "EneCl", "Float_t", true, "nclu"));

  // Tablice Float_t - tracki
  AddVariable("curv", VariableInfo("CurV", "CurV", "Float_t", true, "ntv"));
  AddVariable("phiv", VariableInfo("PhiV", "PhiV", "Float_t", true, "ntv"));
  AddVariable("cotv", VariableInfo("CoTv", "CoTv", "Float_t", true, "ntv"));

  // Tablice Float_t - vertex
  AddVariable("xv", VariableInfo("xV", "xV", "Float_t", true, "nv"));
  AddVariable("yv", VariableInfo("yV", "yV", "Float_t", true, "nv"));
  AddVariable("zv", VariableInfo("zV", "zV", "Float_t", true, "nv"));

  // Tablice Float_t - MC vertex
  AddVariable("xvmc", VariableInfo("xVMC", "xVMC", "Float_t", true, "nvtxmc"));
  AddVariable("yvmc", VariableInfo("yVMC", "yVMC", "Float_t", true, "nvtxmc"));
  AddVariable("zvmc", VariableInfo("zVMC", "zVMC", "Float_t", true, "nvtxmc"));

  // Tablice Float_t - MC momentum
  AddVariable("pxmc", VariableInfo("PxMC", "PxMC", "Float_t", true, "ntmc"));
  AddVariable("pymc", VariableInfo("PyMC", "PyMC", "Float_t", true, "ntmc"));
  AddVariable("pzmc", VariableInfo("PzMC", "PzMC", "Float_t", true, "ntmc"));
}