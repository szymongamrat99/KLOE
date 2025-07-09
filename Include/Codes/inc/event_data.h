#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TVector3.h"

#include <iostream>
#include <vector> // Dla std::vector<float>

struct BhabhaIP
{
  BhabhaIP(TTreeReader &reader) : px(reader, "Bpx"),
                                  py(reader, "Bpy"),
                                  pz(reader, "Bpz"),
                                  energy(reader, "Broots"),
                                  x(reader, "Bx"),
                                  y(reader, "By"),
                                  z(reader, "Bz") {}

  TTreeReaderValue<Float_t> px;
  TTreeReaderValue<Float_t> py;
  TTreeReaderValue<Float_t> pz;
  TTreeReaderValue<Float_t> energy;
  TTreeReaderValue<Float_t> x;
  TTreeReaderValue<Float_t> y;
  TTreeReaderValue<Float_t> z;

  // TVector3 getMomentum() const
  // {
  //   return TVector3(*px, *py, *pz);
  // };

  bool isValid() const
  {
    return px.IsValid() && py.IsValid() && pz.IsValid();
  }
};

struct GeneralEventPropertiesMC
{
  GeneralEventPropertiesMC(TTreeReader &reader) : ntmc(reader, "nTMC"),
                                           nvtxmc(reader, "nVtxMC"),
                                           pidmc(reader, "PidMC"),
                                           mother(reader, "Mother"),
                                           vtxmc(reader, "VtxMC")
                                           {}

  TTreeReaderValue<Int_t> ntmc;
  TTreeReaderValue<Int_t> nvtxmc;

  TTreeReaderArray<Int_t> pidmc;
  TTreeReaderArray<Int_t> mother;
  TTreeReaderArray<Int_t> vtxmc;

  bool isValid() const
  {
    return ntmc.IsValid() && nvtxmc.IsValid() && pidmc.IsValid() && mother.IsValid() && vtxmc.IsValid();
  }

};