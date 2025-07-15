      MODULE ANALYSISMODULE

      IMPLICIT NONE

C= Include /kloe/soft/off/offline/inc/development/tls/maxstructdim.cin =
      INTEGER, PARAMETER :: NeleCluMax = 2000
      INTEGER, PARAMETER ::     MaxNumClu = 100
      INTEGER, PARAMETER ::    MaxNumVtx = 20
      INTEGER, PARAMETER ::    MaxNumTrkV = 30
      INTEGER, PARAMETER ::    MaxNumTrk = 100
      INTEGER, PARAMETER ::    MaxNumDHSP = 1000
      INTEGER, PARAMETER ::    nMaxDC = 1500
      INTEGER, PARAMETER ::    NQihiMax = 1000
      INTEGER, PARAMETER ::   NQcalMax = 32
      INTEGER, PARAMETER ::    MaxNumFirstHit = 300
      INTEGER, PARAMETER ::    MaxNtrkGen = 50
      INTEGER, PARAMETER ::    MaxNvtxGen = 50
      INTEGER, PARAMETER :: TriggerElements = 300
C======================== Include constans.inc =========================
      REAL, PARAMETER ::       Mpip = 139.57039             
      REAL, PARAMETER ::       Mpio = 134.9768              
      REAL, PARAMETER ::       Mko = 497.611             
      REAL, PARAMETER ::       Momega = 782.66          
      REAL, PARAMETER ::      Cvel = 29.9792458         
      REAL, PARAMETER ::       TauKs = 0.08954           
      REAL, PARAMETER ::       Mphi = 1019.461          
C========================= Include klspm00.cin =========================
      INTEGER, PARAMETER :: NCLMIN = 4
      REAL, PARAMETER :: ECLMIN = 20.
      REAL, PARAMETER :: vtxdist = 500.
      INTEGER, PARAMETER :: testcounter = 100
      INTEGER, PARAMETER :: MaxNumFitPar = 38
      INTEGER, PARAMETER :: MaxNumConstr = 10
      INTEGER, PARAMETER :: MaxNumComega = 8
      INTEGER, PARAMETER :: cutN = 6
      INTEGER, PARAMETER :: errN = 17
      CHARACTER*40 messageCut(testcounter),messageErr(testcounter)
      LOGICAL    McFlag,EventOK,ErrFlag,truthmc,truthregg,truthomegaa,
     &  truthsemii,truththreee,truthelsee,                                
     &           truthdoublee
      INTEGER    nmctruth0,ntruth,nmctruth(8),ntruthregg,ntruthsemii,
     &           ntruththreee,ntruthomegaa,ntruthelsee,ntruthdoublee                
      INTEGER    selected,analised
      INTEGER    counter(testcounter),counterBcg(testcounter)
      INTEGER    ErrFlagCount(testcounter),ErrFlagCountBcg(testcounter)
      LOGICAL    kinfitsig,kinfitomega,planerec,closestrec,doublepipi,
     &  onlytruemc,makecuts                                               
      COMMON / flagcommon / McFlag,EventOK,ErrFlag,truthmc
      COMMON / truthcommon / nmctruth0,ntruth,nmctruth
      COMMON / counterscommon / selected,analised,
     &                          counter,counterBcg
      COMMON / errcommon / ErrFlagCount,ErrFlagCountBcg,
     &                     messageCut,messageErr
      COMMON / talk / kinfitsig,kinfitomega,planerec,closestrec,
     &                     doublepipi,onlytruemc,makecuts

      SAVE /flagcommon/, /truthcommon/, /counterscommon/
      SAVE /errcommon/, /talk/    
C== Include /kloe/soft/off/offline/inc/development/tls/evtstruct.cin ===
        TYPE EventInfo
        SEQUENCE
          INTEGER RunNumber
          INTEGER EventNumber
          INTEGER McFlag
          INTEGER EvFlag
          INTEGER Pileup
          INTEGER GenCod
          INTEGER PhiDecay
          INTEGER A1type
          INTEGER A2type
          INTEGER A3type
          INTEGER B1type
          INTEGER B2type
          INTEGER B3type
          INTEGER T3DOWN
          INTEGER T3FLAG
          REAL ECAP(2)
          REAL DCNOISE(4)
        END TYPE
        Integer len_eventinfostru
        Parameter (len_eventinfostru=21)
        TYPE (EventInfo) INFO
C== Include /kloe/soft/off/offline/inc/development/tls/emcstruct.cin ===
        TYPE EmcCluster
          SEQUENCE
         INTEGER n
         INTEGER nmc
         REAL    E(MaxNumCLu)
         REAL    T(MaxNumCLu)
         REAL    X(MaxNumCLu)
         REAL    Y(MaxNumCLu)
         REAL    Z(MaxNumCLu)
         REAL    XA(MaxNumCLu)
         REAL    YA(MaxNumCLu)
         REAL    ZA(MaxNumCLu)
         REAL    Xrms(MaxNumCLu)
         REAL    Yrms(MaxNumCLu)
         REAL    Zrms(MaxNumCLu)
         REAL    XArms(MaxNumCLu)
         REAL    YArms(MaxNumCLu)
         REAL    ZArms(MaxNumCLu)
         REAL    Trms(MaxNumClu)
         INTEGER Flag(MaxNumCLu)
         INTEGER npart(MaxNumCLu)
         INTEGER part1(MaxNumCLu)
         INTEGER pid1 (MaxNumClu)
         INTEGER part2(MaxNumCLu)
         INTEGER pid2 (MaxNumClu)
         INTEGER part3(MaxNumCLu)
         INTEGER pid3 (MaxNumClu)
      END TYPE
        TYPE (EmcCluster) CLUSTER
C=== Include /kloe/soft/off/offline/inc/development/tls/vtxstru.cin ====
        TYPE Vertex
         SEQUENCE
         INTEGER n
         INTEGER Ntrk(MaxNumvtx)
         REAL    X(MaxNumVtx)
         REAL    Y(MaxNumVtx)
         REAL    Z(MaxNumVtx)
         REAL    COV1(MaxNumVtx)
         REAL    COV2(MaxNumVtx)
         REAL    COV3(MaxNumVtx)
         REAL    COV4(MaxNumVtx)
         REAL    COV5(MaxNumVtx)
         REAL    COV6(MaxNumVtx)
         REAL    CHI2(MaxNumVtx)
         INTEGER QUAL(MaxNumVtx)
         INTEGER FITID(MaxNumVtx)
        END TYPE
        TYPE (Vertex) VTX
        TYPE TracksVertex
         SEQUENCE
         INTEGER n
         INTEGER ivOld(MaxNumtrkv)
         INTEGER TrkPoi(MaxNumtrkv)
         REAL    cur(MaxNumtrkv)
         REAL    phi(MaxNumtrkv)
         REAL    cot(MaxNumtrkv)
         REAL    px(MaxNumtrkv)
         REAL    py(MaxNumtrkv)
         REAL    pz(MaxNumtrkv)
         REAL    pmod(MaxNumtrkv)
         REAL    Length(MaxNumtrkv)
         REAL    CHI2(MaxNumTrkV)
         INTEGER ipid(MaxNumtrkv)
         REAL cov11(MaxNumTrkV)
         REAL cov12(MaxNumTrkV)
         REAL cov13(MaxNumTrkV)
         REAL cov22(MaxNumTrkV)
         REAL cov23(MaxNumTrkV)
         REAL cov33(MaxNumTrkV)
        END TYPE
        TYPE (TracksVertex) TRKV
C=== Include /kloe/soft/off/offline/inc/development/tls/trkstru.cin ====
        TYPE AllTracks
        SEQUENCE
         INTEGER n
         INTEGER trkind(MaxNumTrk)
         INTEGER version(MaxNumTrk)
         REAL    cur(MaxNumtrk)
         REAL    phi(MaxNumTrk)
         REAL    cot(MaxNumTrk)
         REAL    px(MaxNumTrk)
         REAL    py(MaxNumTrk)
         REAL    pz(MaxNumTrk)
         REAL    Pmod(MaxNumTrk)
         REAL    x(MaxNumTrk)
         REAL    y(MaxNumTrk)
         REAL    z(MaxNumTrk)
         REAL    Length(MaxNumTrk)
         REAL    curlast(MaxNumTrk)
         REAL    philast(MaxNumTrk)
         REAL    cotlast(MaxNumTrk)
         REAL    pxlast(MaxNumTrk)
         REAL    pylast(MaxNumTrk)
         REAL    pzlast(MaxNumTrk)
         REAL    Pmodlast(MaxNumTrk)
         REAL    xlast (MaxNumTrk)
         REAL    ylast (MaxNumTrk)
         REAL    zlast (MaxNumTrk)
         REAL    xpca  (MaxNumTrk)
         REAL    ypca  (MaxNumTrk)
         REAL    zpca  (MaxNumTrk)
         REAL    Qtrk  (MaxNumTrk)
         REAL    CotPca(MaxNumTrk)
         REAL    PhiPca(MaxNumTrk)
         INTEGER NumPRhit(MaxNumTrk)
         INTEGER NumFitHit(MaxNumTrk)
         REAL    CHI2FIT(MaxNumTrk)
         REAL    CHI2MS(MaxNumTrk)
         REAL    SigPCA(MaxNumTrk)
         REAL    SigZeta(MaxNumTrk)
         REAL    SigCurv(MaxNumTrk)
         REAL    SigCot(MaxNumTrk)
         REAL    SigPhi(MaxNumTrk)
         INTEGER NMSkink(MaxNumTrk)
        END TYPE
        TYPE (AllTracks) TRK
        TYPE AllTracksMC
        SEQUENCE
         INTEGER n
         INTEGER ncontr(MaxNumTrk)
         INTEGER kine1(MaxNumTrk)
         INTEGER type1(MaxNumTrk)
         INTEGER hits1(MaxNumTrk)
         INTEGER kine2(MaxNumTrk)
         INTEGER type2(MaxNumTrk)
         INTEGER hits2(MaxNumTrk)
         INTEGER kine3(MaxNumTrk)
         INTEGER type3(MaxNumTrk)
         INTEGER hits3(MaxNumTrk)
         REAL xfirst(MaxNumTrk)
         REAL yfirst(MaxNumTrk)
         REAL zfirst(MaxNumTrk)
         REAL pxfirst(MaxNumTrk)
         REAL pyfirst(MaxNumTrk)
         REAL pzfirst(MaxNumTrk)
         REAL xlast(MaxNumTrk)
         REAL ylast(MaxNumTrk)
         REAL zlast(MaxNumTrk)
         REAL pxlast(MaxNumTrk)
         REAL pylast(MaxNumTrk)
         REAL pzlast(MaxNumTrk)
         REAL xmcfirst(MaxNumTrk)
         REAL ymcfirst(MaxNumTrk)
         REAL zmcfirst(MaxNumTrk)
         REAL pxmcfirst(MaxNumTrk)
         REAL pymcfirst(MaxNumTrk)
         REAL pzmcfirst(MaxNumTrk)
         REAL xmclast(MaxNumTrk)
         REAL ymclast(MaxNumTrk)
         REAL zmclast(MaxNumTrk)
         REAL pxmclast(MaxNumTrk)
         REAL pymclast(MaxNumTrk)
         REAL pzmclast(MaxNumTrk)
        END TYPE
        TYPE (AllTracksMC) TRKMC
C== Include /kloe/soft/off/offline/inc/development/tls/tclostruct.cin ==
      INTEGER NTCLOMax
      PARAMETER (NTCLOMax = 40)
      TYPE TrackCluster
        SEQUENCE
         INTEGER nt
         INTEGER verver(NTCLOMax)
         INTEGER trknum(NTCLOMax)
         INTEGER clunum(NTCLOMax)
         REAL    xext(NTCLOMax)
         REAL    yext(NTCLOMax)
         REAL    zext(NTCLOMax)
         REAL    leng(NTCLOMax)
         REAL    chi(NTCLOMax)
         REAL    px(NTCLOmax)
         REAL    py(NTCLOmax)
         REAL    pz(NTCLOmax)
      END TYPE
      TYPE(TrackCluster) TCLO
C= Include /kloe/soft/off/offline/inc/development/tls/geanfistruct.cin =
        TYPE GeanfiInformation
         SEQUENCE
         INTEGER Ntrk
         INTEGER kin(maxNtrkGen)
         INTEGER Pid(MaxNtrkGen)
         INTEGER virmom(maxNtrkGen)
         INTEGER Indv(MaxNtrkGen)
         REAL    Px(MaxNtrkGen)
         REAL    Py(MaxNtrkGen)
         REAL    Pz(MaxNtrkGen)
         REAL    Xcv(MaxNtrkGen)
         REAL    ycv(MaxNtrkGen)
         REAL    zcv(MaxNtrkGen)
         REAL    tofcv(MaxNtrkGen)
         REAL    Theta(MaxNtrkGen)
         REAL    Phi(MaxNtrkGen)
         INTEGER ndchmc (MaxNtrkGen)
         INTEGER nlaymc (MaxNtrkGen)
         INTEGER TrkFlag(MaxNtrkGen)
         REAL    Tofmc  (MaxNtrkGen)
         REAL    TrkLen (MaxNtrkGen)
         REAL    xfhmc(MaxNtrkGen)
         REAL    yfhmc(MaxNtrkGen)
         REAL    zfhmc(MaxNtrkGen)
         REAL    pxfhmc(MaxNtrkGen)
         REAL    pyfhmc(MaxNtrkGen)
         REAL    pzfhmc(MaxNtrkGen)
         REAL    xlhmc(MaxNtrkGen)
         REAL    ylhmc(MaxNtrkGen)
         REAL    zlhmc(MaxNtrkGen)
         REAL    pxlhmc(MaxNtrkGen)
         REAL    pylhmc(MaxNtrkGen)
         REAL    pzlhmc(MaxNtrkGen)
         INTEGER NumVtx
         INTEGER Kinmom(MaxNvtxGen)
         INTEGER motherOld(MaxNvtxGen)
         REAL    Tof(MaxNvtxGen)
         REAL    Xv(MaxNvtxGen)
         REAL    Yv(MaxNvtxGen)
         REAL    Zv(MaxNvtxGen)
         REAL    TrkVtx(MaxNvtxGen)
        END TYPE
        TYPE( GeanfiInformation) MC
C=== Include /kloe/soft/off/offline/inc/development/tls/eclostru.cin ===
      INTEGER MaxNumCLINF
      PARAMETER (MaxNumCLINF = 100)     
      TYPE ECLOStru
        SEQUENCE
          INTEGER  n
          INTEGER  TotWord(MaxNumCLINF)
          INTEGER  idpart(MaxNumCLINF)
          INTEGER  dtclpo(MaxNumCLINF)
          INTEGER  dvvnpo(MaxNumCLINF)
          INTEGER  stre(MaxNumCLINF)
          INTEGER  algo(MaxNumCLINF)
          INTEGER  n2
          INTEGER  TotWord2(MaxNumCLINF)
          INTEGER  idpart2(MaxNumCLINF)
          INTEGER  dtclpo2(MaxNumCLINF)
          INTEGER  dvvnpo2(MaxNumCLINF)
          INTEGER  stre2(MaxNumCLINF)
          INTEGER  algo2(MaxNumCLINF)
        END TYPE
      TYPE (ECLOStru) CLINF
C=== Include /kloe/soft/off/offline/inc/development/tls/t0struct.cin ===
      TYPE t0struct
        SEQUENCE
          REAL     dc_step0
          REAL     hit_step0
          REAL     clus_step0
          REAL     step1
          REAL     cable
          REAL     tbunch
          REAL     tphased_mc
        END TYPE
      TYPE (t0struct) T0STRU
C=== Include /kloe/soft/off/offline/inc/development/tls/eclsstru.cin ===
      INTEGER MaxNumOverlapStream
      PARAMETER (MaxNumOverlapStream = 8)       
      TYPE ECLSStru
        SEQUENCE
          INTEGER  n
          INTEGER  Trigger
          INTEGER  Filfo
          INTEGER  totword(MaxNumOverlapStream)
          INTEGER  stream(MaxNumOverlapStream)
          INTEGER  tagnum(MaxNumOverlapStream)
          INTEGER  evntyp(MaxNumOverlapStream)
          INTEGER  n2
          INTEGER  Trigger2
          INTEGER  Filfo2
          INTEGER  totword2(MaxNumOverlapStream)
          INTEGER  stream2(MaxNumOverlapStream)
          INTEGER  tagnum2(MaxNumOverlapStream)
          INTEGER  evntyp2(MaxNumOverlapStream)
        END TYPE
      TYPE (ECLSStru) ECLS
C== Include /kloe/soft/off/offline/inc/development/tls/bposstruct.cin ==
      TYPE  Bposition
        sequence
        REAL px
        REAL py
        REAL pz
        REAL errpx
        REAL errpy
        REAL errpz
        REAL larpx
        REAL elarpx
        REAL larpy
        REAL elarpy
        REAL larpz
        REAL elarpz
        REAL x
        REAL y
        REAL z
        REAL errX
        REAL errY
        REAL errz
        REAL lumx
        REAL elumx
        REAL lumz
        REAL elumz
        REAL Ene
        REAL ErrEne
        REAL Dum1
        REAL ErrDum1
        REAL Dum2
        REAL ErrDum2
      END TYPE
      TYPE (Bposition) BPOS
C====================== Include interfstruct.cin =======================
      TYPE NtupleEvent
        SEQUENCE
          TYPE (EventInfo)         INFO
          TYPE (EmcCluster)        CLU
          TYPE (Vertex)            VTX
          TYPE (TracksVertex)      TRKV
          TYPE (AllTracks)         TRK
          TYPE (AllTracksMC)       TRKMC
          TYPE (TrackCluster)      TCLO
          TYPE (GeanfiInformation) MC
          TYPE (t0struct)          T0STRU
          TYPE (ECLSStru)          ECLS
          TYPE (Bposition)         BPOS
      END TYPE
      TYPE ( NtupleEvent) Evt
      TYPE interfstru
        SEQUENCE
          INTEGER nev
          INTEGER pileup
          INTEGER gcod
          INTEGER phid
          INTEGER a1typ
          INTEGER a2typ
          INTEGER a3typ
          INTEGER b1typ
          INTEGER b2typ
          INTEGER b3typ
          INTEGER nrundata
          REAL    tphased_mc
          REAL    t0dc0
          REAL    t0hit0
          REAL    t0clu0
          REAL    T0step1
          REAL    DelayCable
          REAL    Tbunch
          INTEGER TimeSec
          INTEGER TimeMusec
          INTEGER mcflag
          REAL    Bpx
          REAL    Bpy
          REAL    Bpz
          REAL    Bx
          REAL    By
          REAL    Bz
          REAL    Bwidpx
          REAL    Bwidpy
          REAL    Bwidpz
          REAL    Bsx
          REAL    Bsy
          REAL    Bsz
          REAL    Blumx
          REAL    Blumz
          REAL    Broots
          REAL    BrootsErr
          INTEGER necls2
          INTEGER ECLtrgw2
          INTEGER ECLfilfo2
          INTEGER ECLword2(MaxNumOverlapStream)
          INTEGER ECLstream2(MaxNumOverlapStream)
          INTEGER ECLtagnum2(MaxNumOverlapStream)
          INTEGER ECLevtype2(MaxNumOverlapStream)
          INTEGER nclu
          REAL    EneCl(MaxNumClu)
          REAL    TclOld(MaxNumClu)
          REAL    Xcl(MaxNumClu)
          REAL    Ycl(MaxNumClu)
          REAL    Zcl(MaxNumClu)
          REAL    Xacl(MaxNumClu)
          REAL    Yacl(MaxNumClu)
          REAL    Zacl(MaxNumClu)
          REAL    XRmCl(MaxNumClu)
          REAL    YRmsCl(MaxNumClu)
          REAL    ZrmsCl(MaxNumClu)
          REAL    TrmsCl(MaxNumClu)
          INTEGER FlagCl(MaxNumClu)
          INTEGER nclumc
          INTEGER Npar(MaxNumClu)
          INTEGER Pnum1(MaxNumClu)
          INTEGER Pid1(MaxNumClu)
          INTEGER Pnum2(MaxNumClu)
          INTEGER Pid2(MaxNumClu)
          INTEGER Pnum3(MaxNumClu)
          INTEGER Pid3(MaxNumClu)
          INTEGER ntv
          INTEGER ivOld(MaxNumtrkv)
          INTEGER trknumv(MaxNumtrkv)
          REAL    CurV(MaxNumtrkv)
          REAL    PhiV(MaxNumtrkv)
          REAL    CotV(MaxNumtrkv)
          REAL    PxTV(MaxNumtrkv)
          REAL    PyTV(MaxNumtrkv)
          REAL    PzTV(MaxNumtrkv)
          INTEGER nv
          INTEGER vtx(MaxNumVtx)
          REAL    xvOld(MaxNumVtx)
          REAL    yvOld(MaxNumVtx)
          REAL    zvOld(MaxNumVtx)
          REAL    chivtx(MaxNumVtx)
          INTEGER qualv(MaxNumVtx)
          INTEGER fitidv(MaxNumVtx)
          REAL    VTXcov1(MaxNumVtx)
          REAL    VTXcov2(MaxNumVtx)
          REAL    VTXcov3(MaxNumVtx)
          REAL    VTXcov4(MaxNumVtx)
          REAL    VTXcov5(MaxNumVtx)
          REAL    VTXcov6(MaxNumVtx)
          INTEGER nt
          INTEGER trkind(MaxNumTrk)
          REAL    chi2fit(MaxNumTrk)
          REAL    chi2ms(MaxNumTrk)
          INTEGER ntfmc
          INTEGER trkine1(MaxNumTrk)
          INTEGER trtype1(MaxNumTrk)
          INTEGER trhits1(MaxNumTrk)
          INTEGER trkine2(MaxNumTrk)
          INTEGER trtype2(MaxNumTrk)
          INTEGER trhits2(MaxNumTrk)
          INTEGER trkine3(MaxNumTrk)
          INTEGER trtype3(MaxNumTrk)
          INTEGER trhits3(MaxNumTrk)
          INTEGER ntmc
          INTEGER kine(MaxNtrkGen)
          INTEGER pidmcOld(MaxNtrkGen)
          INTEGER virmom(MaxNtrkGen)
          REAL    pxmc(MaxNtrkGen)
          REAL    pymc(MaxNtrkGen)
          REAL    pzmc(MaxNtrkGen)
          REAL    themc(MaxNtrkGen)
          REAL    phimc(MaxNtrkGen)
          INTEGER vtxmcOld(MaxNtrkGen)
          INTEGER nvtxmc
          INTEGER kinmom(MaxNvtxGen)
          INTEGER motherOld(MaxNvtxGen)
          REAL    xvmc(MaxNvtxGen)
          REAL    yvmc(MaxNvtxGen)
          REAL    zvmc(MaxNvtxGen)
          REAL    ntvtx(MaxNvtxGen)
          INTEGER ntcl
          INTEGER Asstr(NTCLOMax)
          INTEGER Asscl(NTCLOMax)
          INTEGER verver(NTCLOMax)
          REAL    xext(NTCLOMax)
          REAL    yext(NTCLOMax)
          REAL    zext(NTCLOMax)
          REAL    Assleng(NTCLOMax)
          REAL    AssChi(NTCLOMax)
          REAL    extPx(NTCLOMax)
          REAL    extPy(NTCLOMax)
          REAL    extPz(NTCLOMax)
          INTEGER ncli2
          INTEGER ECLOword2(MaxNumCLINF)
          INTEGER idpart2(MaxNumCLINF)
          INTEGER dtclpo2(MaxNumCLINF)
          INTEGER dvvnpo2(MaxNumCLINF)
          INTEGER stre2(MaxNumCLINF)
          INTEGER algo2(MaxNumCLINF)
          REAL    KchMC(9)
          REAL    KneMC(9)
          REAL    KchRec(9)
          REAL    KchRecKS(9)
          REAL    KchRecKL(9)
          REAL    KchRecClose(9)
          REAL    KneRec(9)
          REAL    KchBoost(9)
          REAL    ip(3)
          REAL    ip_closest(3)
          REAL    ip_plane(3)
          REAL    ipmc(3)
          INTEGER vtaken(3)
          INTEGER vtakenks(3)
          INTEGER vtakenkl(3)
          INTEGER vtakenclose(3)
          INTEGER mcISR
          INTEGER mctruth
          INTEGER ncl
          INTEGER ncll(MaxNumClu)
          INTEGER nclwrong
          INTEGER ncllwrong(MaxNumClu)
          REAL    DlMC
          REAL    DtMC
          REAL    DlBoostLor
          REAL    DtBoostLor
          REAL    DlBoostRec
          REAL    DtBoostRec
          REAL    DlRec
          REAL    DtRec
          REAL    cldist
          REAL    KneRecLor(9)
          REAL    KneRecLorFit(9)
          REAL    pi0(2)
          REAL    minv4gam
          REAL    Rc
          REAL    Rtc
          REAL    Rn
          REAL    Rtn
          REAL    RcMC
          REAL    RtcMC
          REAL    RnMC
          REAL    RtnMC
          REAL    ominv(2)
          REAL    chdist
          REAL    trc(MaxNumClu)
          REAL    trcFit(MaxNumClu)
          REAL    trcv(MaxNumClu)
          REAL    trcvFit(MaxNumClu)
          INTEGER ErrId
          INTEGER CutId
          INTEGER g4taken(4)
          INTEGER g4takenFit(4)
          REAL    g4vtxerr(3)
          REAL    g4vtxerrFit(3)
          REAL    KchFit(9)
          REAL    KneFit(9)
          REAL    chdistFit
          REAL    cldistFit
          REAL    ipFit(3)
          INTEGER nclwrongFit
          INTEGER ncllwrongFit(MaxNumClu)
          REAL    DlFit
          REAL    DtFit
          REAL    RcFit
          REAL    RtcFit
          REAL    RnFit
          REAL    RtnFit
          REAL    pi0Fit(2)
          REAL    minv4gamFit
          REAL    ominvFit(2)
          REAL    Chi2
          REAL    Chi2_w
          INTEGER Niter
          INTEGER Niter_w
          REAL    FitPar(MaxNumFitPar)
          REAL    ErrPar(MaxNumFitPar)
          REAL    BkgFitPar(MaxNumFitPar)
          REAL    BkgErrPar(MaxNumFitPar)
          REAL    FitParStart(MaxNumFitPar)
          REAL    ErrParStart(MaxNumFitPar)
          REAL    PpioOmega(4)
          REAL    PpioOmegaFit(4)
          REAL    P4PriRest(4)
          REAL    P4PriRestFit(4)
          REAL    trk1(4)
          REAL    trk2(4)
          REAL    trk1KS(4)
          REAL    trk2KS(4)
          REAL    trk1KL(4)
          REAL    trk2KL(4)
          REAL    trk1Close(4)
          REAL    trk2Close(4)
          REAL    cosTrk
          REAL    cosTrkKS
          REAL    cosTrkKL
          REAL    cosTrkClose
          REAL    cosTrkCM
          REAL    Qmiss
          INTEGER Nconstr
          INTEGER Nconstr_w
          INTEGER nparfit
          INTEGER gpairtaken(2)
          INTEGER gpairtakenFit(2)
          REAL    PgamRec1fit(4)
          REAL    PgamRec2fit(4)
          REAL    PgamRec3fit(4)
          REAL    PgamRec4fit(4)
          REAL    PgamRec1(4)
          REAL    PgamRec2(4)
          REAL    PgamRec3(4)
          REAL    PgamRec4(4)
          INTEGER ChVtxId
          REAL    Trkk1(3)
          REAL    Trkk2(3)
          REAL    ChaVtx(3)
          REAL    NeuVtx(3)
          REAL    PhiVtx(3)
          REAL    TrkIdx(2)
          INTEGER CluIdx(4)
          INTEGER simok
          REAL    test(testcounter)
          INTEGER vtx_with_two_tracks
      END TYPE
      TYPE ( interfstru ) Interf
      COMMON / InterfCommon / Interf,Evt

      END MODULE ANALYSISMODULE