// Author: Szymon Gamrat
// Date of last update: 26.05.2024

#ifndef FORT_COMMON_H
#define FORT_COMMON_H

#ifdef __cplusplus
extern "C"
{
#endif

  const int NeleCluMax = 2000;
  const int MaxNumClu = 100;
  const int MaxNumVtx = 20;
  const int MaxNumTrkV = 30;
  const int MaxNumTrk = 100;
  const int MaxNumDHSP = 1000;
  const int nMaxDC = 1500;
  const int NQihiMax = 1000;
  const int NQcalMax = 32;
  const int MaxNumFirstHit = 300;
  const int MaxNtrkGen = 50;
  const int MaxNvtxGen = 50;
  const int TriggerElements = 300;

//======================== Include constans.inc =========================
      const float       Mpip = 139.57039;             
      const float       Mpio = 134.9768;            
      const float       Mko = 497.611;            
      const float       Momega = 782.66;          
      const float       Cvel = 29.9792458;         
      const float       TauKs = 0.08954;           
      const float       Mphi = 1019.461;          
//========================= Include klspm00.cin =========================
      const int NCLMIN = 4;
      const float ECLMIN = 20.;
      const float vtxdist = 500.;
      const int testcounter = 100;
      const int MaxNumFitPar = 38;
      const int MaxNumConstr = 10;
      const int MaxNumComega = 8;
      const int cutN = 6;
      const int errN = 17;

      const int MaxNumOverlapStream = 8;
      const int MaxNumCLINF = 100;
      const int NTCLOMax = 40;

  extern struct
  {
    int nev;
    int pileup;
    int gcod;
    int phid;
    int a1typ;
    int a2typ;
    int a3typ;
    int b1typ;
    int b2typ;
    int b3typ;
    int nrundata;
    float tphased_mc;
    float t0dc0;
    float t0hit0;
    float t0clu0;
    float T0step1;
    float DelayCable;
    float Tbunch;
    int TimeSec;
    int TimeMusec;
    int mcflag;
    float Bpx;
    float Bpy;
    float Bpz;
    float Bx;
    float By;
    float Bz;
    float Bwidpx;
    float Bwidpy;
    float Bwidpz;
    float Bsx;
    float Bsy;
    float Bsz;
    float Blumx;
    float Blumz;
    float Broots;
    float BrootsErr;
    int necls2;
    int ECLtrgw2;
    int ECLfilfo2;
    int ECLword2[MaxNumOverlapStream];
    int ECLstream2[MaxNumOverlapStream];
    int ECLtagnum2[MaxNumOverlapStream];
    int ECLevtype2[MaxNumOverlapStream];
    int nclu;
    float EneCl[MaxNumClu];
    float Tcl[MaxNumClu];
    float Xcl[MaxNumClu];
    float Ycl[MaxNumClu];
    float Zcl[MaxNumClu];
    float Xacl[MaxNumClu];
    float Yacl[MaxNumClu];
    float Zacl[MaxNumClu];
    float XRmCl[MaxNumClu];
    float YRmsCl[MaxNumClu];
    float ZrmsCl[MaxNumClu];
    float TrmsCl[MaxNumClu];
    int FlagCl[MaxNumClu];
    int nclumc;
    int Npar[MaxNumClu];
    int Pnum1[MaxNumClu];
    int Pid1[MaxNumClu];
    int Pnum2[MaxNumClu];
    int Pid2[MaxNumClu];
    int Pnum3[MaxNumClu];
    int Pid3[MaxNumClu];
    int ntv;
    int iv[MaxNumTrkV];
    int trknumv[MaxNumTrkV];
    float CurV[MaxNumTrkV];
    float PhiV[MaxNumTrkV];
    float CotV[MaxNumTrkV];
    float PxTV[MaxNumTrkV];
    float PyTV[MaxNumTrkV];
    float PzTV[MaxNumTrkV];
    int nv;
    int vtx[MaxNumVtx];
    float xv[MaxNumVtx];
    float yv[MaxNumVtx];
    float zv[MaxNumVtx];
    float chivtx[MaxNumVtx];
    int qualv[MaxNumVtx];
    int fitidv[MaxNumVtx];
    float VTXcov1[MaxNumVtx];
    float VTXcov2[MaxNumVtx];
    float VTXcov3[MaxNumVtx];
    float VTXcov4[MaxNumVtx];
    float VTXcov5[MaxNumVtx];
    float VTXcov6[MaxNumVtx];
    int nt;
    int trkind[MaxNumTrk];
    float chi2fit[MaxNumTrk];
    float chi2ms[MaxNumTrk];
    int ntfmc;
    int trkine1[MaxNumTrk];
    int trtype1[MaxNumTrk];
    int trhits1[MaxNumTrk];
    int trkine2[MaxNumTrk];
    int trtype2[MaxNumTrk];
    int trhits2[MaxNumTrk];
    int trkine3[MaxNumTrk];
    int trtype3[MaxNumTrk];
    int trhits3[MaxNumTrk];
    int ntmc;
    int kine[MaxNtrkGen];
    int pidmc[MaxNtrkGen];
    int virmom[MaxNtrkGen];
    float pxmc[MaxNtrkGen];
    float pymc[MaxNtrkGen];
    float pzmc[MaxNtrkGen];
    float themc[MaxNtrkGen];
    float phimc[MaxNtrkGen];
    int vtxmc[MaxNtrkGen];
    int nvtxmc;
    int kinmom[MaxNvtxGen];
    int mother[MaxNvtxGen];
    float xvmc[MaxNvtxGen];
    float yvmc[MaxNvtxGen];
    float zvmc[MaxNvtxGen];
    float ntvtx[MaxNvtxGen];
    int ntcl;
    int Asstr[NTCLOMax];
    int Asscl[NTCLOMax];
    int verver[NTCLOMax];
    float xext[NTCLOMax];
    float yext[NTCLOMax];
    float zext[NTCLOMax];
    float Assleng[NTCLOMax];
    float AssChi[NTCLOMax];
    float extPx[NTCLOMax];
    float extPy[NTCLOMax];
    float extPz[NTCLOMax];
    int ncli2;
    int ECLOword2[MaxNumCLINF];
    int idpart2[MaxNumCLINF];
    int dtclpo2[MaxNumCLINF];
    int dvvnpo2[MaxNumCLINF];
    int stre2[MaxNumCLINF];
    int algo2[MaxNumCLINF];
    float KchMC[9];
    float KneMC[9];
    float KchRec[9];
    float KchRecKS[9];
    float KchRecKL[9];
    float KchRecClose[9];
    float KneRec[9];
    float KchBoost[9];
    float ip[3];
    float ip_closest[3];
    float ip_plane[3];
    float ipmc[3];
    int vtaken[3];
    int vtakenks[3];
    int vtakenkl[3];
    int vtakenclose[3];
    int mcISR;
    int mctruth;
    int ncl;
    int ncll[MaxNumClu];
    int nclwrong;
    int ncllwrong[MaxNumClu];
    float DlMC;
    float DtMC;
    float DlBoostLor;
    float DtBoostLor;
    float DlBoostRec;
    float DtBoostRec;
    float DlRec;
    float DtRec;
    float cldist;
    float KneRecLor[9];
    float KneRecLorFit[9];
    float pi0[2];
    float minv4gam;
    float Rc;
    float Rtc;
    float Rn;
    float Rtn;
    float RcMC;
    float RtcMC;
    float RnMC;
    float RtnMC;
    float ominv[2];
    float chdist;
    float trc[MaxNumClu];
    float trcFit[MaxNumClu];
    float trcv[MaxNumClu];
    float trcvFit[MaxNumClu];
    int ErrId;
    int CutId;
    int g4taken[4];
    int g4takenFit[4];
    float g4vtxerr[3];
    float g4vtxerrFit[3];
    float KchFit[9];
    float KneFit[9];
    float chdistFit;
    float cldistFit;
    float ipFit[3];
    int nclwrongFit;
    int ncllwrongFit[MaxNumClu];
    float DlFit;
    float DtFit;
    float RcFit;
    float RtcFit;
    float RnFit;
    float RtnFit;
    float pi0Fit[2];
    float minv4gamFit;
    float ominvFit[2];
    float Chi2;
    float Chi2_w;
    int Niter;
    int Niter_w;
    float FitPar[MaxNumFitPar];
    float ErrPar[MaxNumFitPar];
    float BkgFitPar[MaxNumFitPar];
    float BkgErrPar[MaxNumFitPar];
    float FitParStart[MaxNumFitPar];
    float ErrParStart[MaxNumFitPar];
    float PpioOmega[4];
    float PpioOmegaFit[4];
    float P4PriRest[4];
    float P4PriRestFit[4];
    float trk1[4];
    float trk2[4];
    float trk1KS[4];
    float trk2KS[4];
    float trk1KL[4];
    float trk2KL[4];
    float trk1Close[4];
    float trk2Close[4];
    float cosTrk;
    float cosTrkKS;
    float cosTrkKL;
    float cosTrkClose;
    float cosTrkCM;
    float Qmiss;
    int Nconstr;
    int Nconstr_w;
    int nparfit;
    int gpairtaken[2];
    int gpairtakenFit[2];
    float PgamRec1fit[4];
    float PgamRec2fit[4];
    float PgamRec3fit[4];
    float PgamRec4fit[4];
    float PgamRec1[4];
    float PgamRec2[4];
    float PgamRec3[4];
    float PgamRec4[4];
    int ChVtxId;
    float Trkk1[3];
    float Trkk2[3];
    float ChaVtx[3];
    float NeuVtx[3];
    float PhiVtx[3];
    float TrkIdx[2];
    int CluIdx[4];
    int simok;
    float test[testcounter];
    int vtx_with_two_tracks;
  } interfcommon_;

#ifdef __cplusplus
}
#endif

#endif // !FORT_COMMON_H