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

  TVector3 getMomentum() const
  {
    return TVector3(*px, *py, *pz);
  };

  bool isValid() const
  {
    return px.IsValid() && py.IsValid() && pz.IsValid();
  }
};

struct TrackProperties
{
  TrackProperties(TTreeReader &reader) : curv(reader, "CurV"),
                                  phiv(reader, "PhiV"),
                                  cotv(reader, "CotV"),
                                  energy(reader, "Broots"),
                                  x(reader, "Bx"),
                                  y(reader, "By"),
                                  z(reader, "Bz") {}

  TTreeReaderArray<Float_t> curv;
  TTreeReaderArray<Float_t> phiv;
  TTreeReaderArray<Float_t> cotv;
  TTreeReaderValue<Float_t> energy;
  TTreeReaderValue<Float_t> x;
  TTreeReaderValue<Float_t> y;
  TTreeReaderValue<Float_t> z;

  TVector3 getMomentum() const
  {
    return TVector3(*px, *py, *pz);
  };

  bool isValid() const
  {
    return px.IsValid() && py.IsValid() && pz.IsValid();
  }
};