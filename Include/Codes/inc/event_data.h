#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TVector3.h"

#include <iostream>
#include <vector> // Dla std::vector<float>

struct BhabhaIP
{
  BhabhaIP(TTreeReader &reader) : px(reader, "BPx"),
                                  py(reader, "BPy"),
                                  pz(reader, "BPz"),
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
                                                  vtxmc(reader, "VtxMC"),
                                                  xvmc(reader, "xVMC"),
                                                  yvmc(reader, "yVMC"),
                                                  zvmc(reader, "zVMC"),
                                                  pxmc(reader, "PxMC"),
                                                  pymc(reader, "PyMC"),
                                                  pzmc(reader, "PzMC")
  {
  }

  TTreeReaderValue<Int_t> ntmc;
  TTreeReaderValue<Int_t> nvtxmc;

  TTreeReaderArray<Int_t> pidmc;
  TTreeReaderArray<Int_t> mother;
  TTreeReaderArray<Int_t> vtxmc;

  TTreeReaderArray<Float_t> xvmc;
  TTreeReaderArray<Float_t> yvmc;
  TTreeReaderArray<Float_t> zvmc;
  TTreeReaderArray<Float_t> pxmc;
  TTreeReaderArray<Float_t> pymc;
  TTreeReaderArray<Float_t> pzmc;

  bool isValid() const
  {
    return ntmc.IsValid() && nvtxmc.IsValid() && pidmc.IsValid() && mother.IsValid() && vtxmc.IsValid();
  }
};

struct ClusterProperties
{
  ClusterProperties(TTreeReader &reader) : nclu(reader, "nClu"),
                                           ntcl(reader, "nTcl"),
                                           asscl(reader, "AssCl"),
                                           t0step1(reader, "T0Step1"),
                                           xcl(reader, "XCl"),
                                           ycl(reader, "YCl"),
                                           zcl(reader, "ZCl"),
                                           tcl(reader, "TCl"),
                                           enecl(reader, "EneCl")
  {
  }

  TTreeReaderValue<Int_t> nclu;
  TTreeReaderValue<Int_t> ntcl;

  TTreeReaderArray<Int_t> asscl;

  TTreeReaderValue<Float_t> t0step1;

  TTreeReaderArray<Float_t> xcl;
  TTreeReaderArray<Float_t> ycl;
  TTreeReaderArray<Float_t> zcl;
  TTreeReaderArray<Float_t> tcl;
  TTreeReaderArray<Float_t> enecl;

  bool isValid() const
  {
    return nclu.IsValid() && ntcl.IsValid() && asscl.IsValid();
  }
};

struct ChargedVertexProperties
{
  ChargedVertexProperties(TTreeReader &reader) : nv(reader, "nV"),
                                                 ntv(reader, "nTv"),
                                                 iv(reader, "iV"),
                                                 Curv(reader, "CurV"),
                                                 Phiv(reader, "PhiV"),
                                                 Cotv(reader, "CoTv"),
                                                 xv(reader, "xV"),
                                                 yv(reader, "yV"),
                                                 zv(reader, "zV")
  {
  }

  TTreeReaderValue<Int_t> nv;
  TTreeReaderValue<Int_t> ntv;

  TTreeReaderArray<Int_t> iv;

  TTreeReaderArray<Float_t> Curv;
  TTreeReaderArray<Float_t> Phiv;
  TTreeReaderArray<Float_t> Cotv;
  TTreeReaderArray<Float_t> xv;
  TTreeReaderArray<Float_t> yv;
  TTreeReaderArray<Float_t> zv;
};

struct GeneralEventProperties
{
  GeneralEventProperties(TTreeReader &reader) : nev(reader, "nEv"),
                                                nrun(reader, "nRun"),
                                                necls(reader, "NEcls"),
                                                eclfilfo(reader, "EclFilfo"),
                                                // eclfilfoword(reader, "EclFilfoWord"),
                                                eclstream(reader, "EclStream")
  {
  }

  TTreeReaderValue<Int_t> nev;
  TTreeReaderValue<Int_t> nrun;
  TTreeReaderValue<Int_t> necls;
  TTreeReaderValue<Int_t> eclfilfo;
  // TTreeReaderValue<Int_t> eclfilfoword;

  TTreeReaderArray<Int_t> eclstream;
};