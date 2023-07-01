C=======================================================================
C=======================================================================
      SUBROUTINE clearstruct
C
        IMPLICIT NONE
C========================= Include klspm00.inc =========================
      INTEGER NCLMIN
      PARAMETER( NCLMIN = 4 )
      REAL ECLMIN
      PARAMETER( ECLMIN = 20. )
      REAL vtxdist
      PARAMETER( vtxdist = 50. )
      INTEGER testcounter
      PARAMETER( testcounter = 100 )
      INTEGER MaxNumFitPar
      PARAMETER( MaxNumFitPar = 38 )
      INTEGER MaxNumConstr
      PARAMETER( MaxNumConstr = 10 )
      INTEGER MaxNumComega
      PARAMETER( MaxNumComega = 8 )
      INTEGER cutN
      PARAMETER( cutN = 6 )
      INTEGER errN
      PARAMETER( errN = 17 )
      CHARACTER*40 messageCut(testcounter),messageErr(testcounter)
      LOGICAL    McFlag,EventOK,ErrFlag,truthmc,truthregg,truthomegaa,tr
     1  uthsemii,truththreee,truthelsee                                 
      INTEGER    nmctruth0,ntruth,nmctruth(7),ntruthregg,ntruthsemii,ntr
     1  uththreee,ntruthomegaa,ntruthelsee                              
      INTEGER    selected,analised
      INTEGER    counter(testcounter),counterBcg(testcounter)
      INTEGER    ErrFlagCount(testcounter),ErrFlagCountBcg(testcounter)
      LOGICAL    cutchar,cutneu,onlytruemc,makecuts
      COMMON / flagcommon / McFlag,EventOK,ErrFlag,truthmc
      COMMON / truthcommon / nmctruth0,ntruth,nmctruth
      COMMON / counterscommon / selected,analised,
     &                          counter,counterBcg
      COMMON / errcommon / ErrFlagCount,ErrFlagCountBcg,
     &                     messageCut,messageErr
      COMMON / talk / cutchar,cutneu,onlytruemc,makecuts
      SAVE /flagcommon/, /truthcommon/, /counterscommon/
      SAVE /errcommon/, /talk/
C-----------------------------------------------------------------------
C= Include /kloe/soft/off/offline/inc/development/tls/maxstructdim.cin =
      INTEGER     NeleCluMax
      PARAMETER ( NeleCluMax = 2000 )
      INTEGER     MaxNumClu
      PARAMETER ( MaxNumClu = 100 )
      INTEGER    MaxNumVtx
      PARAMETER (MaxNumVtx = 20)
      INTEGER    MaxNumTrkV
      PARAMETER (MaxNumTrkV = 30)
      INTEGER    MaxNumTrk
      PARAMETER (MaxNumTrk = 100)
      INTEGER    MaxNumDHSP
      PARAMETER (MaxNumDHSP = 1000)
      Integer    nMaxDC
      Parameter (nMaxDC=1500)
      INTEGER    NQihiMax
      PARAMETER (NQihiMax=1000)
      INTEGER    NQcalMax
      PARAMETER (NQcalMax= 32 )
      INTEGER    MaxNumFirstHit
      PARAMETER (MaxNumFirstHit = 300 )
      INTEGER    MaxNtrkGen
      PARAMETER (MaxNtrkGen =50)
      INTEGER    MaxNvtxGen
      PARAMETER (MaxNvTxGen =30)
      integer TriggerElements
      parameter (TriggerElements = 300 )
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
         INTEGER iv(MaxNumtrkv)
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
         INTEGER mother(MaxNvtxGen)
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
C====================== Include interfstruct.inc =======================
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
          REAL    Tcl(MaxNumClu)
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
          INTEGER iv(MaxNumtrkv)
          INTEGER trknumv(MaxNumtrkv)
          REAL    CurV(MaxNumtrkv)
          REAL    PhiV(MaxNumtrkv)
          REAL    CotV(MaxNumtrkv)
          REAL    PxTV(MaxNumtrkv)
          REAL    PyTV(MaxNumtrkv)
          REAL    PzTV(MaxNumtrkv)
          INTEGER nv
          INTEGER vtx(MaxNumVtx)
          REAL    xv(MaxNumVtx)
          REAL    yv(MaxNumVtx)
          REAL    zv(MaxNumVtx)
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
          INTEGER pidmc(MaxNtrkGen)
          INTEGER virmom(MaxNtrkGen)
          REAL    pxmc(MaxNtrkGen)
          REAL    pymc(MaxNtrkGen)
          REAL    pzmc(MaxNtrkGen)
          REAL    themc(MaxNtrkGen)
          REAL    phimc(MaxNtrkGen)
          INTEGER vtxmc(MaxNtrkGen)
          INTEGER nvtxmc
          INTEGER kinmom(MaxNvtxGen)
          INTEGER mother(MaxNvtxGen)
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
          REAL    KchMC(9)
          REAL    KneMC(9)
          REAL    KchRec(9)
          REAL    KneRec(9)
          REAL    KchBoost(9)
          REAL    ip(3)
          REAL    ipmc(3)
          INTEGER vtaken(3)
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
          REAL    cosTrk
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
          INTEGER omegatruth_top
          INTEGER clu_index(4)
          INTEGER chtrk_index(2)
          INTEGER chvtx_index
          REAL P4gam_w(4,4)
          REAL P4ch_w(2,4)
          REAL P4ch_tot_w(4)
          REAL om_dist(4)
      END TYPE
      TYPE ( interfstru) Interf
      COMMON / InterfCommon / Interf,Evt
C-----------------------------------------------------------------------
      INTEGER i
C integers
      Interf%nev=0
      Interf%pileup=0
      Interf%gcod=0
      Interf%phid=0
      Interf%a1typ=0
      Interf%a2typ=0
      Interf%a3typ=0
      Interf%b1typ=0
      Interf%b2typ=0
      Interf%b3typ=0
      Interf%nrundata=0
      Interf%TimeSec=0
      Interf%TimeMusec=0
      Interf%mcflag=0
      Interf%necls2=0
      Interf%ECLtrgw2=0
      Interf%ECLfilfo2=0
      Interf%nclu=0
      Interf%nclumc=0
      Interf%ntv=0
      Interf%nv=0
      Interf%nt=0
      Interf%ntfmc=0
      Interf%ntmc=0
      Interf%nvtxmc=0
      Interf%ntcl=0
C     Interf%ncli2=0
      Interf%vtaken=0
      Interf%mcISR=0
      Interf%mctruth=0
      Interf%ncl=0
      Interf%nclwrong=0
      Interf%Niter=0
      Interf%Niter_w=0
      Interf%Nconstr=0
      Interf%Nconstr_w=0
      Interf%nparfit=0
      Interf%ErrId=0
      Interf%CutId=0
      Interf%ChVtxId=0
      Interf%simok=0
C reals
      Interf%tphased_mc=0.
      Interf%t0dc0=0.
      Interf%t0hit0=0.
      Interf%t0clu0=0.
      Interf%T0step1=0.
      Interf%DelayCable=0.
      Interf%Tbunch=0.
      Interf%Bpx=0.
      Interf%Bpy=0.
      Interf%Bpz=0.
      Interf%Bx=0.
      Interf%By=0.
      Interf%Bz=0.
      Interf%Bwidpx=0.
      Interf%Bwidpy=0.
      Interf%Bwidpz=0.
      Interf%Bsx=0.
      Interf%Bsy=0.
      Interf%Bsz=0.
      Interf%Blumx=0.
      Interf%Blumz=0.
      Interf%Broots=0.
      Interf%BrootsErr=0.
      Interf%DlMC=0.
      Interf%DtMC=0.
      Interf%DlBoostLor=0.
      Interf%DtBoostLor=0.
      Interf%DlBoostRec=0.
      Interf%DtBoostRec=0.
      Interf%DlRec=0.
      Interf%DtRec=0.
      Interf%DlFit=0.
      Interf%DtFit=0.
      Interf%cldist=0.
      Interf%cldistFit=0.
      Interf%minv4gam=0.
      Interf%minv4gamFit=0.
      Interf%RcFit=0.
      Interf%RtcFit=0.
      Interf%RnFit=0.
      Interf%RtnFit=0.
      Interf%Rc=0.
      Interf%Rtc=0.
      Interf%Rn=0.
      Interf%Rtn=0.
      Interf%RcMC=0.
      Interf%RtcMC=0.
      Interf%RnMC=0.
      Interf%RtnMC=0.
      Interf%ominvFit=0.
      Interf%chdist=0.
      Interf%chdistFit=0.
      Interf%Chi2=0.
      Interf%Chi2_w=0.
      Interf%cosTrk=0.
      Interf%Qmiss=0.
      Interf%omegatruth_top=-999
      Interf%chvtx_index=-999
C loops
      DO i = 1,MaxNumOverlapStream
          Interf%ECLword2(i)=0
          Interf%ECLstream2(i)=0
          Interf%ECLtagnum2(i)=0
          Interf%ECLevtype2(i)=0
      ENDDO
      DO i = 1,MaxNumClu
          Interf%EneCl(i)=0.
          Interf%Tcl(i)=0.
          Interf%Xcl(i)=0.
          Interf%Ycl(i)=0.
          Interf%Zcl(i)=0.
          Interf%Xacl(i)=0.
          Interf%Yacl(i)=0.
          Interf%Zacl(i)=0.
          Interf%XRmCl(i)=0.
          Interf%YRmsCl(i)=0.
          Interf%ZrmsCl(i)=0.
          Interf%TrmsCl(i)=0.
          Interf%trc(i)=0.
          Interf%trcFit(i)=0.
          Interf%trcv(i)=0.
          Interf%trcvFit(i)=0.
          Interf%FlagCl(i)=0
          Interf%Npar(i)=0
          Interf%Pnum1(i)=0
          Interf%Pid1(i)=0
          Interf%Pnum2(i)=0
          Interf%Pid2(i)=0
          Interf%Pnum3(i)=0
          Interf%Pid3(i)=0
          Interf%ncll(i)=0
          Interf%ncllwrong(i)=0
          Interf%ncllwrongFit(i)=0
      ENDDO
      DO i = 1,MaxNumtrkv
          Interf%iv(i)=0
          Interf%trknumv(i)=0
          Interf%CurV(i)=0.
          Interf%PhiV(i)=0.
          Interf%CotV(i)=0.
          Interf%PxTV(i)=0.
          Interf%PyTV(i)=0.
          Interf%PzTV(i)=0.
      ENDDO
      DO i = 1,MaxNumVtx
          Interf%vtx(i)=0
          Interf%xv(i)=0.
          Interf%yv(i)=0.
          Interf%zv(i)=0.
          Interf%chivtx(i)=0.
          Interf%qualv(i)=0
          Interf%fitidv(i)=0
          Interf%VTXcov1(i)=0.
          Interf%VTXcov2(i)=0.
          Interf%VTXcov3(i)=0.
          Interf%VTXcov4(i)=0.
          Interf%VTXcov5(i)=0.
          Interf%VTXcov6(i)=0.
      ENDDO
      DO i = 1,MaxNumTrk
          Interf%chi2fit(i)=0.
          Interf%chi2ms(i)=0.
          Interf%trkind(i)=0
          Interf%trkine1(i)=0
          Interf%trtype1(i)=0
          Interf%trhits1(i)=0
          Interf%trkine2(i)=0
          Interf%trtype2(i)=0
          Interf%trhits2(i)=0
          Interf%trkine3(i)=0
          Interf%trtype3(i)=0
          Interf%trhits3(i)=0
      ENDDO
      DO i = 1,MaxNtrkGen
          Interf%kine(i)=0
          Interf%pidmc(i)=0
          Interf%virmom(i)=0
          Interf%vtxmc(i)=0
          Interf%kinmom(i)=0
          Interf%mother(i)=0
          Interf%pxmc(i)=0.
          Interf%pymc(i)=0.
          Interf%pzmc(i)=0.
          Interf%themc(i)=0.
          Interf%phimc(i)=0.
          Interf%xvmc(i)=0.
          Interf%yvmc(i)=0.
          Interf%zvmc(i)=0.
          Interf%ntvtx(i)=0.
      ENDDO
      DO i = 1,NTCLOMax
          Interf%Asstr(i)=0
          Interf%Asscl(i)=0
          Interf%verver(i)=0
          Interf%xext(i)=0.
          Interf%yext(i)=0.
          Interf%zext(i)=0.
          Interf%Assleng(i)=0.
          Interf%AssChi(i)=0.
          Interf%extPx(i)=0.
          Interf%extPy(i)=0.
          Interf%extPz(i)=0.
      ENDDO
C     DO i = 1,MaxNumCLINF
C         Interf%ECLOword2(i)=0
C         Interf%idpart2(i)=0
C         Interf%dtclpo2(i)=0
C         Interf%dvvnpo2(i)=0
C         Interf%stre2(i)=0
C         Interf%algo2(i)=0
C     ENDDO
      DO i = 1,MaxNumFitPar
          Interf%FitPar(i)=0.
          Interf%ErrPar(i)=0.
          Interf%BkgFitPar(i)=0.
          Interf%BkgErrPar(i)=0.
          Interf%FitParStart(i)=0.
          Interf%ErrParStart(i)=0.
      ENDDO
      DO i = 1,9
          Interf%KchMC(i)=0.
          Interf%KchRec(i)=0.
          Interf%KchFit(i)=0.
          Interf%KneRec(i)=0.
          Interf%KneFit(i)=0.
          Interf%KchBoost(i)=0.
          Interf%KneRecLor(i)=0.
          Interf%KneRecLorFit(i)=0.
          Interf%KneMC(i)=0.
      ENDDO
      DO i = 1,3
          Interf%ip(i)=0.
          Interf%ipFit(i)=0.
          Interf%ipmc(i)=0.
          Interf%g4vtxerr(i)=0.
          Interf%g4vtxerrFit(i)=0.
          Interf%Trkk1(i)=0.
          Interf%Trkk2(i)=0.
          Interf%ChaVtx(i)=0.
          Interf%NeuVtx(i)=0.
          Interf%PhiVtx(i)=0.
      ENDDO
      DO i = 1,2
          Interf%pi0(i)=0.
          Interf%pi0Fit(i)=0.
          Interf%ominv(i)=0.
          Interf%gpairtaken(i)=0
          Interf%gpairtakenFit(i)=0
          Interf%TrkIdx(i)=0
          Interf%chtrk_index(i)=-999
      ENDDO
      DO i = 1,4
          Interf%PpioOmega(i)=0.
          Interf%PpioOmegaFit(i)=0.
          Interf%P4PriRest(i)=0.
          Interf%P4PriRestFit(i)=0.
          Interf%trk1(i)=0.
          Interf%trk2(i)=0.
          Interf%g4taken(i)=0
          Interf%g4takenFit(i)=0
          Interf%PgamRec1(i)=0.
          Interf%PgamRec2(i)=0.
          Interf%PgamRec3(i)=0.
          Interf%PgamRec4(i)=0.
          Interf%PgamRec1Fit(i)=0.
          Interf%PgamRec2Fit(i)=0.
          Interf%PgamRec3Fit(i)=0.
          Interf%PgamRec4Fit(i)=0.
          Interf%CluIdx(i)=0
          Interf%clu_index(i)=-999
      ENDDO
      DO i = 1,testcounter
          Interf%test(i)=0.
      ENDDO
      END
C=======================================================================
C=======================================================================
C=======================================================================
      SUBROUTINE fillstruct
C
        IMPLICIT NONE
C====== Include /kloe/soft/off/offline/inc/development/runtyp.cin ======
        INTEGER    EXKPHY               
        PARAMETER (EXKPHY       = 0)
        INTEGER    EXCOSR               
        PARAMETER (EXCOSR       = 1)
        INTEGER    EXOFSI               
        PARAMETER (EXOFSI       = 2)
        INTEGER    EXTEST               
        PARAMETER (EXTEST       = 3)
        INTEGER    RUKPHY               
        PARAMETER (RUKPHY       = 1)
        INTEGER    RUCOSR               
        PARAMETER (RUCOSR       = 2)
        INTEGER    RUCALI               
        PARAMETER (RUCALI       = 3)
C====== Include /kloe/soft/off/offline/inc/development/jobsta.cin ======
        COMMON /JOBSTI/ NRUN, NEV, NRUNA, NEVA, IBEGRN, IENDRN,
     1                  IUDBUG, NEVP, NEVPR, NEVPA, NREC, NRECF
        INTEGER         NRUN, NEV, NRUNA, NEVA, IBEGRN, IENDRN
        INTEGER         IUDBUG, NEVP, NEVPR, NEVPA, NREC, NRECF
        COMMON /JOBLRI/ EXPTYP, PARTNO, RUNTYP, LRECNO, NRCTHS,
     1                  CDFCID, IRTYPA, VERNUM, RELNUM, CRPRID,
     2                  CRYEAR,  CRMON,  CRDAY,
     3                  CRHOUR,  CRMIN,  CRSEC
        INTEGER         EXPTYP, PARTNO, RUNTYP, LRECNO, NRCTHS
        INTEGER         CDFCID, IRTYPA, VERNUM, RELNUM, CRPRID
        INTEGER         CRYEAR,  CRMON,  CRDAY
        INTEGER         CRHOUR,  CRMIN,  CRSEC
        COMMON /JOBEVC/ TRIGNO, TRGSUM(4),
     1                  STEVTI, SDONET, RPTMOD, RESERV, SOPTYP, SBUFNO,
     2                  SREADO(2), SRESIO(2), CREADO(2), TERSUM(2)
        INTEGER         TRIGNO,    TRGSUM
        INTEGER         STEVTI, SDONET, RPTMOD, RESERV, SOPTYP, SBUFNO
        INTEGER         SREADO,    SRESIO,    CREADO,    TERSUM
        COMMON /JOBSTR/ ROOTS, BMAGNT, BMAGFT, BMAGBT
        REAL            ROOTS, BMAGNT, BMAGFT, BMAGBT
C==== Include /kloe/soft/off/a_c/production/aix/library/anerror.cin ====
      INTEGER    ANSUCC
      PARAMETER (ANSUCC = '8198009'X)
      INTEGER    ANBEGN
      PARAMETER (ANBEGN = '8198019'X)
      INTEGER    ANPDUN
      PARAMETER (ANPDUN = '8198021'X)
      INTEGER    ANAST
      PARAMETER (ANAST   = '81980A3'X)
      INTEGER    ANILPS
      PARAMETER (ANILPS = '8198122'X)
      INTEGER    ANILMD
      PARAMETER (ANILMD = '819812A'X)
      INTEGER    ANLSTR
      PARAMETER (ANLSTR = '81980C0'X)
      INTEGER    ANNINI
      PARAMETER (ANNINI = '8198102'X)
      INTEGER    ANOVRF
      PARAMETER (ANOVRF = '819810A'X)
      INTEGER    ANNOME
      PARAMETER (ANNOME = '8198112'X)
      INTEGER    ANNOLU
      PARAMETER (ANNOLU = '819811A'X)
      INTEGER    ANNOEP
      PARAMETER (ANNOEP = '8198011'X)
      INTEGER    ANILGP
      PARAMETER (ANILGP = '8198132'X)
      INTEGER    ANILPA
      PARAMETER (ANILPA = '819813A'X)
      INTEGER    ANILQU
      PARAMETER (ANILQU = '8198142'X)
      INTEGER    ANILSY
      PARAMETER (ANILSY = '819814A'X)
      INTEGER    ANILLU
      PARAMETER (ANILLU = '8198152'X)
      INTEGER    ANILVB
      PARAMETER (ANILVB = '819815A'X)
      INTEGER    ANILLI
      PARAMETER (ANILLI = '8198162'X)
      INTEGER    ANDUMD
      PARAMETER (ANDUMD = '819816A'X)
      INTEGER    ANYBER
      PARAMETER (ANYBER = '8198172'X)
      INTEGER    ANYBFE
      PARAMETER (ANYBFE = '081981B2'X)
      INTEGER    ANEFOP
      PARAMETER (ANEFOP = '819817A'X)
      INTEGER    ANEOF
      PARAMETER (ANEOF = '8198083'X)
      INTEGER    ANIFNO
      PARAMETER (ANIFNO = '819808B'X)
      INTEGER    ANCLOS
      PARAMETER (ANCLOS = '8198093'X)
      INTEGER    ANOFNO
      PARAMETER (ANOFNO = '819809B'X)
      INTEGER    ANNOLR
      PARAMETER (ANNOLR = '8198182'X)
      INTEGER    ANILST
      PARAMETER (ANILST = '819818A'X)
      INTEGER    ANNOEV
      PARAMETER (ANNOEV = '8198192'X)
      INTEGER    ANILBL
      PARAMETER (ANILBL = '819819A'X)
      INTEGER    ANNEIF
      PARAMETER (ANNEIF = '81981A2'X)
      INTEGER    ANNOBK
      PARAMETER (ANNOBK = '81980C8'X)
      INTEGER    ANILBK
      PARAMETER (ANILBK = '81981AA'X)
C==== Include /kloe/soft/off/ybos/production/aix/library/errcod.cin ====
      INTEGER    YESUCC
      PARAMETER (YESUCC = '08508009'X)
      INTEGER    YEWLOK
      PARAMETER (YEWLOK = '08508083'X)
      INTEGER    YEWRNS
      PARAMETER (YEWRNS = '0850808B'X)
      INTEGER    YEINDX
      PARAMETER (YEINDX = '08508093'X)
      INTEGER    YEEOF
      PARAMETER (YEEOF  = '08508100'X)
      INTEGER    YEPRBD
      PARAMETER (YEPRBD = '08508108'X)
      INTEGER    YERSTF
      PARAMETER (YERSTF = '08508110'X)
      INTEGER    YEOVRF
      PARAMETER (YEOVRF = '08508182'X)
      INTEGER    YEBKNG
      PARAMETER (YEBKNG = '0850818A'X)
      INTEGER    YEBKTS
      PARAMETER (YEBKTS = '08508192'X)
      INTEGER    YEBKOV
      PARAMETER (YEBKOV = '0850819A'X)
      INTEGER    YEBKSP
      PARAMETER (YEBKSP = '085081A2'X)
      INTEGER    YENOBK
      PARAMETER (YENOBK = '085081AA'X)
      INTEGER    YENONR
      PARAMETER (YENONR = '085081B2'X)
      INTEGER    YEMANY
      PARAMETER (YEMANY = '085081BA'X)
      INTEGER    YEDUPB
      PARAMETER (YEDUPB = '085081C2'X)
      INTEGER    YEBSIL
      PARAMETER (YEBSIL = '085081CA'X)
      INTEGER    YEILOP
      PARAMETER (YEILOP = '085081D2'X)
      INTEGER    YEAROF
      PARAMETER (YEAROF = '085081DA'X)
      INTEGER    YETRUN
      PARAMETER (YETRUN = '085081E2'X)
      INTEGER    YECORR
      PARAMETER (YECORR = '085081EA'X)
      INTEGER    YEWRGI
      PARAMETER (YEWRGI = '085081F2'X)
      INTEGER    YELWIU
      PARAMETER (YELWIU = '085081FA'X)
      INTEGER    YEILUN
      PARAMETER (YEILUN = '08508202'X)
      INTEGER    YEIOUS
      PARAMETER (YEIOUS = '0850820A'X)
      INTEGER    YELROF
      PARAMETER (YELROF = '08508212'X)
      INTEGER    YELRLN
      PARAMETER (YELRLN = '0850821A'X)
      INTEGER    YEPRXL
      PARAMETER (YEPRXL = '08508222'X)
      INTEGER    YEBKOS
      PARAMETER (YEBKOS = '0850822A'X)
      INTEGER    YEBKXD
      PARAMETER (YEBKXD = '08508232'X)
      INTEGER    YELRZL
      PARAMETER (YELRZL = '0850823A'X)
      INTEGER    YELRWT
      PARAMETER (YELRWT = '08508242'X)
      INTEGER    YEDKOP
      PARAMETER (YEDKOP = '0850824A'X)
      INTEGER    YEDKWT
      PARAMETER (YEDKWT = '08508252'X)
      INTEGER    YEDKRD
      PARAMETER (YEDKRD = '0850825A'X)
      INTEGER    YEILUS
      PARAMETER (YEILUS = '08508282'X)
      INTEGER    YERSUS
      PARAMETER (YERSUS = '0850828A'X)
      INTEGER    YENRUS
      PARAMETER (YENRUS = '08508292'X)
      INTEGER    YEILLA
      PARAMETER (YEILLA = '08508302'X)
      INTEGER    YEDUPA
      PARAMETER (YEDUPA = '0850830A'X)
      INTEGER    YEBTIL
      PARAMETER (YEBTIL = '08508382'X)
      INTEGER    YEGRIL
      PARAMETER (YEGRIL = '0850838A'X)
      INTEGER    YECHKF
      PARAMETER (YECHKF = '08508392'X)
      INTEGER    YESTOF
      PARAMETER (YESTOF = '0850839A'X)
      INTEGER    YENOTA
      PARAMETER (YENOTA = '085083A2'X)
      INTEGER    YEBNKP
      PARAMETER (YEBNKP = '085083AA'X)
      INTEGER    YEBDCH
      PARAMETER (YEBDCH = '085083C2'X)
      INTEGER    YEBDBK
      PARAMETER (YEBDBK = '085083CA'X)
      INTEGER    YEBDLU
      PARAMETER (YEBDLU = '085083D2'X)
      INTEGER    YEBDBS
      PARAMETER (YEBDBS = '085083DA'X)
C====== Include /kloe/soft/off/offline/inc/development/erlevl.cin ======
        INTEGER    ERMINI
        PARAMETER (ERMINI = 0)
        INTEGER    ERMAXI
        PARAMETER (ERMAXI = 7)
        INTEGER    ERSUCC
        PARAMETER (ERSUCC = 0)
        INTEGER    ERINFO
        PARAMETER (ERINFO = 1)
        INTEGER    ERWARN
        PARAMETER (ERWARN = 2)
        INTEGER    EREROR
        PARAMETER (EREROR = 3)
        INTEGER    ERSEVR
        PARAMETER (ERSEVR = 4)
        INTEGER         ELSUCC
        PARAMETER       (ELSUCC = 10)
        INTEGER         ELINFO
        PARAMETER       (ELINFO = 11)
        INTEGER         ELWARN
        PARAMETER       (ELWARN = 12)
        INTEGER         ELWARN2
        PARAMETER       (ELWARN2 = 13)
        INTEGER         ELEROR
        PARAMETER       (ELEROR = 14)
        INTEGER         ELEROR2
        PARAMETER       (ELEROR2 = 15)
        INTEGER         ELNEXT
        PARAMETER       (ELNEXT = 16)
        INTEGER         ELNEXT2
        PARAMETER       (ELNEXT2 = 17)
        INTEGER         ELSEVR
        PARAMETER       (ELSEVR = 18)
        INTEGER         ELABOR
        PARAMETER       (ELABOR = 19)
C====== Include /kloe/soft/off/offline/inc/development/bpoybs.cin ======
      INTEGER RUNGID
      PARAMETER (RUNGID = 0)
      INTEGER RUNVER
      PARAMETER (RUNVER = 1)
      INTEGER RUNRND
      PARAMETER (RUNRND = 3)
      INTEGER RUNEVS
      PARAMETER (RUNEVS = 5)
      INTEGER PARNUM
      PARAMETER (PARNUM = 0)
      INTEGER PARNAM
      PARAMETER (PARNAM = 1)
      INTEGER PARTRT
      PARAMETER (PARTRT = 6)
      INTEGER PARMAS
      PARAMETER (PARMAS = 7)
      INTEGER PARCHR
      PARAMETER (PARCHR = 8)
      INTEGER PARLTI
      PARAMETER (PARLTI = 9)
      INTEGER PARBRT
      PARAMETER (PARBRT = 10)
      INTEGER PARDMO
      PARAMETER (PARDMO = 16)
      INTEGER PAUSER
      PARAMETER (PAUSER = 22)
      INTEGER MATNUM
      PARAMETER (MATNUM = 0)
      INTEGER MATNAM
      PARAMETER (MATNAM = 1)
      INTEGER MATATN
      PARAMETER (MATATN = 6)
      INTEGER MATMAS
      PARAMETER (MATMAS = 7)
      INTEGER MATDEN
      PARAMETER (MATDEN = 8)
      INTEGER MATRLE
      PARAMETER (MATRLE = 9)
      INTEGER MATABL
      PARAMETER (MATABL = 10)
      INTEGER MAUSER
      PARAMETER (MAUSER = 11)
      INTEGER TMENUM
      PARAMETER (TMENUM = 0)
      INTEGER TMENAM
      PARAMETER (TMENAM = 1)
      INTEGER TMEMAT
      PARAMETER (TMEMAT = 6)
      INTEGER TMESTF
      PARAMETER (TMESTF = 7)
      INTEGER TMEMFI
      PARAMETER (TMEMFI = 8)
      INTEGER TMEMFM
      PARAMETER (TMEMFM = 9)
      INTEGER TMEMAA
      PARAMETER (TMEMAA = 10)
      INTEGER TMEMAS
      PARAMETER (TMEMAS = 11)
      INTEGER TMEMAE
      PARAMETER (TMEMAE = 12)
      INTEGER TMEEPS
      PARAMETER (TMEEPS = 13)
      INTEGER TMEMIS
      PARAMETER (TMEMIS = 14)
      INTEGER TMUSER
      PARAMETER (TMUSER = 15)
      INTEGER BRIHVB
      PARAMETER (BRIHVB = 0)
      INTEGER BRIHVS
      PARAMETER (BRIHVS = 1)
      INTEGER BRIVD1
      PARAMETER (BRIVD1 = 2)
      INTEGER BRIPRS
      PARAMETER (BRIPRS = 10)
      INTEGER BRITMP
      PARAMETER (BRITMP = 12)
      INTEGER BRIISO
      PARAMETER (BRIISO = 14)
      INTEGER BRIMGC
      PARAMETER (BRIMGC = 17)
      INTEGER BRIVD2
      PARAMETER (BRIVD2 = 18)
      INTEGER BRIVD3
      PARAMETER (BRIVD3 = 32)
      INTEGER BRINCR
      PARAMETER (BRINCR = 64)
      INTEGER BRINDB
      PARAMETER (BRINDB = 65)
      INTEGER BRINCB
      PARAMETER (BRINCB = 66)
      INTEGER BRINDH
      PARAMETER (BRINDH = 67)
      INTEGER HEARUN
      PARAMETER (HEARUN = 0)
      INTEGER HEAEVT
      PARAMETER (HEAEVT = 1)
      INTEGER HEARND
      PARAMETER (HEARND = 2)
      INTEGER HEAQMI
      PARAMETER (HEAQMI = 4)
      INTEGER HEAKIN
      PARAMETER (HEAKIN = 5)
      INTEGER HEAVER
      PARAMETER (HEAVER = 6)
      INTEGER HEAUSE
      PARAMETER (HEAUSE = 7)
      INTEGER KINPPX
      PARAMETER (KINPPX = 0)
      INTEGER KINPPY
      PARAMETER (KINPPY = 1)
      INTEGER KINPPZ
      PARAMETER (KINPPZ = 2)
      INTEGER KINPEN
      PARAMETER (KINPEN = 3)
      INTEGER KINPTY
      PARAMETER (KINPTY = 4)
      INTEGER KINVXO
      PARAMETER (KINVXO = 5)
      INTEGER KINNVX
      PARAMETER (KINNVX = 6)
      INTEGER KINVLS
      PARAMETER (KINVLS = 7)
      INTEGER VERXCO
      PARAMETER (VERXCO = 0)
      INTEGER VERYCO
      PARAMETER (VERYCO = 1)
      INTEGER VERZCO
      PARAMETER (VERZCO = 2)
      INTEGER VERTOF
      PARAMETER (VERTOF = 3)
      INTEGER VERTNO
      PARAMETER (VERTNO = 4)
      INTEGER VERTAR
      PARAMETER (VERTAR = 5)
      INTEGER VERNTR
      PARAMETER (VERNTR = 6)
      INTEGER VERTLS
      PARAMETER (VERTLS = 7)
      INTEGER DHIHDS
      PARAMETER (DHIHDS = 0)
      INTEGER DHIVRN
      PARAMETER (DHIVRN = 1)
      INTEGER DHINRO
      PARAMETER (DHINRO = 2)
      INTEGER DHINCO
      PARAMETER (DHINCO = 3)
      INTEGER DHIPTY
      PARAMETER (DHIPTY = 0)
      INTEGER DHITRA
      PARAMETER (DHITRA = 1)
      INTEGER DHIADR
      PARAMETER (DHIADR = 2)
      INTEGER DHIXCO
      PARAMETER (DHIXCO = 3)
      INTEGER DHIYCO
      PARAMETER (DHIYCO = 4)
      INTEGER DHIZCO
      PARAMETER (DHIZCO = 5)
      INTEGER DHIPPX
      PARAMETER (DHIPPX = 6)
      INTEGER DHIPPY
      PARAMETER (DHIPPY = 7)
      INTEGER DHIPPZ
      PARAMETER (DHIPPZ = 8)
      INTEGER DHITOF
      PARAMETER (DHITOF = 9)
      INTEGER DHIELO
      PARAMETER (DHIELO =10)
      INTEGER DHILEN
      PARAMETER (DHILEN =11)
      INTEGER DHITIM
      PARAMETER (DHITIM =12)
      INTEGER DHIDIS
      PARAMETER (DHIDIS =13)
      INTEGER DHIFLA
      PARAMETER (DHIFLA =14)
      INTEGER VHIHDS
      PARAMETER (VHIHDS = 0)
      INTEGER VHIVRN
      PARAMETER (VHIVRN = 1)
      INTEGER VHINRO
      PARAMETER (VHINRO = 2)
      INTEGER VHINCO
      PARAMETER (VHINCO = 3)
      INTEGER VHIPTY
      PARAMETER (VHIPTY = 0)
      INTEGER VHITRA
      PARAMETER (VHITRA = 1)
      INTEGER VHIADR
      PARAMETER (VHIADR = 2)
      INTEGER VHIXEN
      PARAMETER (VHIXEN = 3)
      INTEGER VHIYEN
      PARAMETER (VHIYEN = 4)
      INTEGER VHIZEN
      PARAMETER (VHIZEN = 5)
      INTEGER VHIXEX
      PARAMETER (VHIXEX = 6)
      INTEGER VHIYEX
      PARAMETER (VHIYEX = 7)
      INTEGER VHIZEX
      PARAMETER (VHIZEX = 8)
      INTEGER VHIPMO
      PARAMETER (VHIPMO = 9)
      INTEGER VHITOF
      PARAMETER (VHITOF =10)
      INTEGER VHIELO
      PARAMETER (VHIELO =11)
      INTEGER VHITXC
      PARAMETER (VHITXC =12)
      INTEGER VHITYC
      PARAMETER (VHITYC =13)
      INTEGER CHIHDS
      PARAMETER (CHIHDS = 0)
      INTEGER CHIVRN
      PARAMETER (CHIVRN = 1)
      INTEGER CHINRO
      PARAMETER (CHINRO = 2)
      INTEGER CHINCO
      PARAMETER (CHINCO = 3)
      INTEGER CHIPTY
      PARAMETER (CHIPTY = 0)
      INTEGER CHITRA
      PARAMETER (CHITRA = 1)
      INTEGER CHIADR
      PARAMETER (CHIADR = 2)
      INTEGER CHIXCO
      PARAMETER (CHIXCO = 3)
      INTEGER CHIYCO
      PARAMETER (CHIYCO = 4)
      INTEGER CHIZCO
      PARAMETER (CHIZCO = 5)
      INTEGER CHITOF
      PARAMETER (CHITOF = 6)
      INTEGER CHIELO
      PARAMETER (CHIELO = 7)
      INTEGER CHILEN
      PARAMETER (CHILEN = 8)
      INTEGER CFHHDS
      PARAMETER (CFHHDS = 0)
      INTEGER CFHVRN
      PARAMETER (CFHVRN = 1)
      INTEGER CFHNRO
      PARAMETER (CFHNRO = 2)
      INTEGER CFHNCO
      PARAMETER (CFHNCO = 3)
      INTEGER CFHPTY
      PARAMETER (CFHPTY = 0)
      INTEGER CFHTRA
      PARAMETER (CFHTRA = 1)
      INTEGER CFHADR
      PARAMETER (CFHADR = 2)
      INTEGER CFHXCO
      PARAMETER (CFHXCO = 3)
      INTEGER CFHYCO
      PARAMETER (CFHYCO = 4)
      INTEGER CFHZCO
      PARAMETER (CFHZCO = 5)
      INTEGER CFHPPX
      PARAMETER (CFHPPX = 6)
      INTEGER CFHPPY
      PARAMETER (CFHPPY = 7)
      INTEGER CFHPPZ
      PARAMETER (CFHPPZ = 8)
      INTEGER CFHTOF
      PARAMETER (CFHTOF = 9)
      INTEGER CFHLEN
      PARAMETER (CFHLEN = 10)
      INTEGER AHIHDS
      PARAMETER (AHIHDS = 0)
      INTEGER AHIVRN
      PARAMETER (AHIVRN = 1)
      INTEGER AHINRO
      PARAMETER (AHINRO = 2)
      INTEGER AHINCO
      PARAMETER (AHINCO = 3)
      INTEGER AHIPTY
      PARAMETER (AHIPTY = 0)
      INTEGER AHITRA
      PARAMETER (AHITRA = 1)
      INTEGER AHIADR
      PARAMETER (AHIADR = 2)
      INTEGER AHIXCO
      PARAMETER (AHIXCO = 3)
      INTEGER AHIYCO
      PARAMETER (AHIYCO = 4)
      INTEGER AHIZCO
      PARAMETER (AHIZCO = 5)
      INTEGER AHIPPX
      PARAMETER (AHIPPX = 6)
      INTEGER AHIPPY
      PARAMETER (AHIPPY = 7)
      INTEGER AHIPPZ
      PARAMETER (AHIPPZ = 8)
      INTEGER AHITOF
      PARAMETER (AHITOF = 9)
      INTEGER AHIELO
      PARAMETER (AHIELO =10)
      INTEGER AHILEN
      PARAMETER (AHILEN =11)
      INTEGER QIHHDS
      PARAMETER (QIHHDS = 0)
      INTEGER QIHVRN
      PARAMETER (QIHVRN = 1)
      INTEGER QIHNRO
      PARAMETER (QIHNRO = 2)
      INTEGER QIHNCO
      PARAMETER (QIHNCO = 3)
      INTEGER QIHPTY
      PARAMETER (QIHPTY = 0)
      INTEGER QIHTRA
      PARAMETER (QIHTRA = 1)
      INTEGER QIHADR
      PARAMETER (QIHADR = 2)
      INTEGER QIHXCO
      PARAMETER (QIHXCO = 3)
      INTEGER QIHYCO
      PARAMETER (QIHYCO = 4)
      INTEGER QIHZCO
      PARAMETER (QIHZCO = 5)
      INTEGER QIHPPX
      PARAMETER (QIHPPX = 6)
      INTEGER QIHPPY
      PARAMETER (QIHPPY = 7)
      INTEGER QIHPPZ
      PARAMETER (QIHPPZ = 8)
      INTEGER QIHTOF
      PARAMETER (QIHTOF = 9)
      INTEGER QIHELO
      PARAMETER (QIHELO =10)
      INTEGER QIHLEN
      PARAMETER (QIHLEN =11)
      INTEGER QHIHDS
      PARAMETER (QHIHDS = 0)
      INTEGER QHIVRN
      PARAMETER (QHIVRN = 1)
      INTEGER QHINRO
      PARAMETER (QHINRO = 2)
      INTEGER QHINCO
      PARAMETER (QHINCO = 3)
      INTEGER QHIPTY
      PARAMETER (QHIPTY = 0)
      INTEGER QHITRA
      PARAMETER (QHITRA = 1)
      INTEGER QHIADR
      PARAMETER (QHIADR = 2)
      INTEGER QHIXCO
      PARAMETER (QHIXCO = 3)
      INTEGER QHIYCO
      PARAMETER (QHIYCO = 4)
      INTEGER QHIZCO
      PARAMETER (QHIZCO = 5)
      INTEGER QHIPPX
      PARAMETER (QHIPPX = 6)
      INTEGER QHIPPY
      PARAMETER (QHIPPY = 7)
      INTEGER QHIPPZ
      PARAMETER (QHIPPZ = 8)
      INTEGER QHITOF
      PARAMETER (QHITOF = 9)
      INTEGER QHIELO
      PARAMETER (QHIELO =10)
      INTEGER QHILEN
      PARAMETER (QHILEN =11)
      INTEGER ITHHDS
      PARAMETER (ITHHDS = 0)
      INTEGER ITHVRN
      PARAMETER (ITHVRN = 1)
      INTEGER ITHNRO
      PARAMETER (ITHNRO = 2)
      INTEGER ITHNCO
      PARAMETER (ITHNCO = 3)
      INTEGER ITHINX
      PARAMETER (ITHINX = 0)
      INTEGER ITHINY
      PARAMETER (ITHINY = 1)
      INTEGER ITHINZ
      PARAMETER (ITHINZ = 2)
      INTEGER ITHOUX
      PARAMETER (ITHOUX = 3)
      INTEGER ITHOUY
      PARAMETER (ITHOUY = 4)
      INTEGER ITHOUZ
      PARAMETER (ITHOUZ = 5)
      INTEGER ITHPIX
      PARAMETER (ITHPIX = 6)
      INTEGER ITHPIY
      PARAMETER (ITHPIY = 7)
      INTEGER ITHPIZ
      PARAMETER (ITHPIZ = 8)
      INTEGER ITHTRL
      PARAMETER (ITHTRL = 9)
      INTEGER ITHDER
      PARAMETER (ITHDER = 10)
      INTEGER ITHTOF
      PARAMETER (ITHTOF = 11)
      INTEGER ITHPLA
      PARAMETER (ITHPLA = 12)
      INTEGER ITHTYP
      PARAMETER (ITHTYP = 13)
      INTEGER ITHTRA
      PARAMETER (ITHTRA = 14)
      INTEGER QTHHDS
      PARAMETER (QTHHDS = 0)
      INTEGER QTHVRN
      PARAMETER (QTHVRN = 1)
      INTEGER QTHNRO
      PARAMETER (QTHNRO = 2)
      INTEGER QTHNCO
      PARAMETER (QTHNCO = 3)
      INTEGER QTHPTY
      PARAMETER (QTHPTY = 0)
      INTEGER QTHTRA
      PARAMETER (QTHTRA = 1)
      INTEGER QTHADR
      PARAMETER (QTHADR = 2)
      INTEGER QTHXCO
      PARAMETER (QTHXCO = 3)
      INTEGER QTHYCO
      PARAMETER (QTHYCO = 4)
      INTEGER QTHZCO
      PARAMETER (QTHZCO = 5)
      INTEGER QTHTOF
      PARAMETER (QTHTOF = 6)
      INTEGER QTHELO
      PARAMETER (QTHELO = 7)
      INTEGER QTHLEN
      PARAMETER (QTHLEN = 8)
      INTEGER CTHHDS
      PARAMETER (CTHHDS = 0)
      INTEGER CTHVRN
      PARAMETER (CTHVRN = 1)
      INTEGER CTHNRO
      PARAMETER (CTHNRO = 2)
      INTEGER CTHNCO
      PARAMETER (CTHNCO = 3)
      INTEGER CTHPTY
      PARAMETER (CTHPTY = 0)
      INTEGER CTHTRA
      PARAMETER (CTHTRA = 1)
      INTEGER CTHADR
      PARAMETER (CTHADR = 2)
      INTEGER CTHXCO
      PARAMETER (CTHXCO = 3)
      INTEGER CTHYCO
      PARAMETER (CTHYCO = 4)
      INTEGER CTHZCO
      PARAMETER (CTHZCO = 5)
      INTEGER CTHTOF
      PARAMETER (CTHTOF = 6)
      INTEGER CTHELO
      PARAMETER (CTHELO = 7)
      INTEGER CTHLEN
      PARAMETER (CTHLEN = 8)
      INTEGER LEHHDS
      PARAMETER (LEHHDS = 0)
      INTEGER LEHVRN
      PARAMETER (LEHVRN = 1)
      INTEGER LEHNRO
      PARAMETER (LEHNRO = 2)
      INTEGER LEHNCO
      PARAMETER (LEHNCO = 3)
      INTEGER LEHPTY
      PARAMETER (LEHPTY = 0)
      INTEGER LEHTRA
      PARAMETER (LEHTRA = 1)
      INTEGER LEHADR
      PARAMETER (LEHADR = 2)
      INTEGER LEHXCO
      PARAMETER (LEHXCO = 3)
      INTEGER LEHYCO
      PARAMETER (LEHYCO = 4)
      INTEGER LEHZCO
      PARAMETER (LEHZCO = 5)
      INTEGER LEHXOU
      PARAMETER (LEHXOU = 6)
      INTEGER LEHYOU
      PARAMETER (LEHYOU = 7)
      INTEGER LEHZOU
      PARAMETER (LEHZOU = 8)
      INTEGER LEHTOF
      PARAMETER (LEHTOF = 9)
      INTEGER LEHELO
      PARAMETER (LEHELO = 10)
      INTEGER LEHLEN
      PARAMETER (LEHLEN = 11)
      INTEGER DTCHDS
      PARAMETER (DTCHDS = 0)
      INTEGER DTCVRN
      PARAMETER (DTCVRN = 1)
      INTEGER DTCNRO
      PARAMETER (DTCNRO = 2)
      INTEGER DTCNCO
      PARAMETER (DTCNCO = 3)
      INTEGER DTCCAL
      PARAMETER (DTCCAL = 4)
      INTEGER DTCADR
      PARAMETER (DTCADR = 0)
      INTEGER DTCTIM
      PARAMETER (DTCTIM = 1)
      INTEGER DTKHDS
      PARAMETER (DTKHDS = 0)
      INTEGER DTKVRN
      PARAMETER (DTKVRN = 1)
      INTEGER DTKNRO
      PARAMETER (DTKNRO = 2)
      INTEGER DTKNTR
      PARAMETER (DTKNTR = 0)
      INTEGER DTKADR
      PARAMETER (DTKADR = 1)
      INTEGER DTHHDS
      PARAMETER (DTHHDS = 0)
      INTEGER DTHVRN
      PARAMETER (DTHVRN = 1)
      INTEGER DTHNRO
      PARAMETER (DTHNRO = 2)
      INTEGER DTHNHI
      PARAMETER (DTHNHI = 0)
      INTEGER DHNHDS
      PARAMETER (DHNHDS = 0)
      INTEGER DHNVRN
      PARAMETER (DHNVRN = 1)
      INTEGER DHNNRO
      PARAMETER (DHNNRO = 2)
      INTEGER DHNNHI
      PARAMETER (DHNNHI = 0)
      INTEGER CELHDS
      PARAMETER (CELHDS = 0)
      INTEGER CELVRN
      PARAMETER (CELVRN = 1)
      INTEGER CELCAL
      PARAMETER (CELCAL = 2)
      INTEGER CELNRO
      PARAMETER (CELNRO = 3)
      INTEGER CELNCO
      PARAMETER (CELNCO = 4)
      INTEGER CELADR
      PARAMETER (CELADR = 0)
      INTEGER CELEA
      PARAMETER (CELEA = 1)
      INTEGER CELEB
      PARAMETER (CELEB = 2)
      INTEGER CELTA
      PARAMETER (CELTA = 3)
      INTEGER CELTB
      PARAMETER (CELTB = 4)
      INTEGER CEKHDS
      PARAMETER (CEKHDS = 0)
      INTEGER CEKVRN
      PARAMETER (CEKVRN = 1)
      INTEGER CEKNRO
      PARAMETER (CEKNRO = 2)
      INTEGER CEKNTR
      PARAMETER (CEKNTR = 0)
      INTEGER CEKADR
      PARAMETER (CEKADR = 1)
      INTEGER CHNHDS
      PARAMETER (CHNHDS = 0)
      INTEGER CHNVRN
      PARAMETER (CHNVRN = 1)
      INTEGER CHNNRO
      PARAMETER (CHNNRO = 2)
      INTEGER CHNNHI
      PARAMETER (CHNNHI = 0)
      INTEGER QCAHDS
      PARAMETER (QCAHDS = 0)
      INTEGER QCAVRN
      PARAMETER (QCAVRN = 1)
      INTEGER QCACAL
      PARAMETER (QCACAL = 2)
      INTEGER QCANRO
      PARAMETER (QCANRO = 3)
      INTEGER QCANCO
      PARAMETER (QCANCO = 4)
      INTEGER QCADR
      PARAMETER (QCADR = 0)
      INTEGER QCAE
      PARAMETER (QCAE = 1)
      INTEGER QCAT
      PARAMETER (QCAT = 2)
      INTEGER QCKHDS
      PARAMETER (QCKHDS = 0)
      INTEGER QCKVRN
      PARAMETER (QCKVRN = 1)
      INTEGER QCKNRO
      PARAMETER (QCKNRO = 2)
      INTEGER QCKNTR
      PARAMETER (QCKNTR = 0)
      INTEGER QCKADR
      PARAMETER (QCKADR = 1)
      INTEGER ITCHDS
      PARAMETER (ITCHDS = 0)
      INTEGER ITCVRN
      PARAMETER (ITCVRN = 1)
      INTEGER ITCNRO
      PARAMETER (ITCNRO = 2)
      INTEGER ITCNCO
      PARAMETER (ITCNCO = 3)
      INTEGER ITCCAL
      PARAMETER (ITCCAL = 4)
      INTEGER ITCADR
      PARAMETER (ITCADR = 0)
      INTEGER ITCTIM
      PARAMETER (ITCTIM = 1)
      INTEGER ITKHDS
      PARAMETER (ITKHDS = 0)
      INTEGER ITKVRN
      PARAMETER (ITKVRN = 1)
      INTEGER ITKNRO
      PARAMETER (ITKNRO = 2)
      INTEGER ITKNTR
      PARAMETER (ITKNTR = 0)
      INTEGER ITKADR
      PARAMETER (ITKADR = 1)
      INTEGER ITLHDS
      PARAMETER (ITLHDS = 0)
      INTEGER ITLVRN
      PARAMETER (ITLVRN = 1)
      INTEGER ITHLRO
      PARAMETER (ITHLRO = 2)
      INTEGER ITHLHI
      PARAMETER (ITHLHI = 0)
      INTEGER IHNHDS
      PARAMETER (IHNHDS = 0)
      INTEGER IHNVRN
      PARAMETER (IHNVRN = 1)
      INTEGER IHNNRO
      PARAMETER (IHNNRO = 2)
      INTEGER IHNNHI
      PARAMETER (IHNNHI = 0)
      INTEGER QCTHDS
      PARAMETER (QCTHDS = 0)
      INTEGER QCTVRN
      PARAMETER (QCTVRN = 1)
      INTEGER QCTCAL
      PARAMETER (QCTCAL = 2)
      INTEGER QCTNRO
      PARAMETER (QCTNRO = 3)
      INTEGER QCTADR
      PARAMETER (QCTADR = 1)
      INTEGER QCTNTD
      PARAMETER (QCTNTD = 0)
      INTEGER QCTTIM
      PARAMETER (QCTTIM = 0)
      INTEGER QCTWID
      PARAMETER (QCTWID = 1)
      INTEGER QTKHDS
      PARAMETER (QTKHDS = 0)
      INTEGER QTKVRN
      PARAMETER (QTKVRN = 1)
      INTEGER QTKNRO
      PARAMETER (QTKNRO = 2)
      INTEGER QTKNTR
      PARAMETER (QTKNTR = 0)
      INTEGER QTKADR
      PARAMETER (QTKADR = 1)
      INTEGER QTLHDS
      PARAMETER (QTLHDS = 0)
      INTEGER QTLVRN
      PARAMETER (QTLVRN = 1)
      INTEGER QTLNRO
      PARAMETER (QTLNRO = 2)
      INTEGER QTLNHI
      PARAMETER (QTLNHI = 0)
      INTEGER THNHDS
      PARAMETER (THNHDS = 0)
      INTEGER THNVRN
      PARAMETER (THNVRN = 1)
      INTEGER THNNRO
      PARAMETER (THNNRO = 2)
      INTEGER THNNHI
      PARAMETER (THNNHI = 0)
      INTEGER CCTHDS
      PARAMETER (CCTHDS = 0)
      INTEGER CCTVRN
      PARAMETER (CCTVRN = 1)
      INTEGER CCTCAL
      PARAMETER (CCTCAL = 2)
      INTEGER CCTNRO
      PARAMETER (CCTNRO = 3)
      INTEGER CCTNCO
      PARAMETER (CCTNCO = 4)
      INTEGER CCTADR
      PARAMETER (CCTADR = 0)
      INTEGER CCTTIM
      PARAMETER (CCTTIM = 2)
      INTEGER CCTELO
      PARAMETER (CCTELO = 1)
      INTEGER CCKHDS
      PARAMETER (CCKHDS = 0)
      INTEGER CCKVRN
      PARAMETER (CCKVRN = 1)
      INTEGER CCKNRO
      PARAMETER (CCKNRO = 2)
      INTEGER CCKNTR
      PARAMETER (CCKNTR = 0)
      INTEGER CCKADR
      PARAMETER (CCKADR = 1)
      INTEGER CCLHDS
      PARAMETER (CCLHDS = 0)
      INTEGER CCLVRN
      PARAMETER (CCLVRN = 1)
      INTEGER CCLNRO
      PARAMETER (CCLNRO = 2)
      INTEGER CCLNHI
      PARAMETER (CCLNHI = 0)
      INTEGER YHNHDS
      PARAMETER (YHNHDS = 0)
      INTEGER YHNVRN
      PARAMETER (YHNVRN = 1)
      INTEGER YHNNRO
      PARAMETER (YHNNRO = 2)
      INTEGER YHNNHI
      PARAMETER (YHNNHI = 0)
      INTEGER LETHDS
      PARAMETER (LETHDS = 0)
      INTEGER LETVRN
      PARAMETER (LETVRN = 1)
      INTEGER LETCAL
      PARAMETER (LETCAL = 2)
      INTEGER LETNRO
      PARAMETER (LETNRO = 3)
      INTEGER LETNCO
      PARAMETER (LETNCO = 4)
      INTEGER LETADR
      PARAMETER (LETADR = 0)
      INTEGER LETTIM
      PARAMETER (LETTIM = 2)
      INTEGER LETELO
      PARAMETER (LETELO = 1)
      INTEGER LEKHDS
      PARAMETER (LEKHDS = 0)
      INTEGER LEKVRN
      PARAMETER (LEKVRN = 1)
      INTEGER LEKNRO
      PARAMETER (LEKNRO = 2)
      INTEGER LEKNTR
      PARAMETER (LEKNTR = 0)
      INTEGER LEKADR
      PARAMETER (LEKADR = 1)
      INTEGER LELHDS
      PARAMETER (LELHDS = 0)
      INTEGER LELVRN
      PARAMETER (LELVRN = 1)
      INTEGER LELNRO
      PARAMETER (LELNRO = 2)
      INTEGER LELNHI
      PARAMETER (LELNHI = 0)
      INTEGER LENHDS
      PARAMETER (LENHDS = 0)
      INTEGER LENVRN
      PARAMETER (LENVRN = 1)
      INTEGER LENNRO
      PARAMETER (LENNRO = 2)
      INTEGER LENNHI
      PARAMETER (LENNHI = 0)
      INTEGER DTDHDS
      PARAMETER (DTDHDS = 0)
      INTEGER DTDVRN
      PARAMETER (DTDVRN = 1)
      INTEGER DTDNRO
      PARAMETER (DTDNRO = 2)
      INTEGER DTDNCO
      PARAMETER (DTDNCO = 3)
      INTEGER DTDZES
      PARAMETER (DTDZES = 4)
      INTEGER DTDADR
      PARAMETER (DTDADR = 0)
      INTEGER DADHDS
      PARAMETER (DADHDS = 0)
      INTEGER DADVRN
      PARAMETER (DADVRN = 1)
      INTEGER DADNRO
      PARAMETER (DADNRO = 2)
      INTEGER DADNCO
      PARAMETER (DADNCO = 3)
      INTEGER DADPES
      PARAMETER (DADPES = 4)
      INTEGER DADADR
      PARAMETER (DADADR = 0)
      INTEGER CTDHDS
      PARAMETER (CTDHDS = 0)
      INTEGER CTDVRN
      PARAMETER (CTDVRN = 1)
      INTEGER CTDNRO
      PARAMETER (CTDNRO = 2)
      INTEGER CTDNCO
      PARAMETER (CTDNCO = 3)
      INTEGER CTDADR
      PARAMETER (CTDADR = 0)
      INTEGER CADHDS
      PARAMETER (CADHDS = 0)
      INTEGER CADVRN
      PARAMETER (CADVRN = 1)
      INTEGER CADNRO
      PARAMETER (CADNRO = 2)
      INTEGER CADNCO
      PARAMETER (CADNCO = 3)
      INTEGER CADADR
      PARAMETER (CADADR = 0)
C======= Include /kloe/soft/off/offline/inc/development/bcs.cin ========
        INTEGER          IW
        REAL             RW(200)
        DOUBLE PRECISION DW(100)
        CHARACTER*4      AW(200)
        INTEGER*2        IW2(400)
        EQUIVALENCE (IW,IW2,RW,DW)
        EQUIVALENCE (IW,AW)
        COMMON /BCS/ IW(200)
C========================= Include klspm00.inc =========================
      INTEGER NCLMIN
      PARAMETER( NCLMIN = 4 )
      REAL ECLMIN
      PARAMETER( ECLMIN = 20. )
      REAL vtxdist
      PARAMETER( vtxdist = 50. )
      INTEGER testcounter
      PARAMETER( testcounter = 100 )
      INTEGER MaxNumFitPar
      PARAMETER( MaxNumFitPar = 38 )
      INTEGER MaxNumConstr
      PARAMETER( MaxNumConstr = 10 )
      INTEGER MaxNumComega
      PARAMETER( MaxNumComega = 8 )
      INTEGER cutN
      PARAMETER( cutN = 6 )
      INTEGER errN
      PARAMETER( errN = 17 )
      CHARACTER*40 messageCut(testcounter),messageErr(testcounter)
      LOGICAL    McFlag,EventOK,ErrFlag,truthmc,truthregg,truthomegaa,tr
     1  uthsemii,truththreee,truthelsee                                 
      INTEGER    nmctruth0,ntruth,nmctruth(7),ntruthregg,ntruthsemii,ntr
     1  uththreee,ntruthomegaa,ntruthelsee                              
      INTEGER    selected,analised
      INTEGER    counter(testcounter),counterBcg(testcounter)
      INTEGER    ErrFlagCount(testcounter),ErrFlagCountBcg(testcounter)
      LOGICAL    cutchar,cutneu,onlytruemc,makecuts
      COMMON / flagcommon / McFlag,EventOK,ErrFlag,truthmc
      COMMON / truthcommon / nmctruth0,ntruth,nmctruth
      COMMON / counterscommon / selected,analised,
     &                          counter,counterBcg
      COMMON / errcommon / ErrFlagCount,ErrFlagCountBcg,
     &                     messageCut,messageErr
      COMMON / talk / cutchar,cutneu,onlytruemc,makecuts
      SAVE /flagcommon/, /truthcommon/, /counterscommon/
      SAVE /errcommon/, /talk/
C-----------------------------------------------------------------------
C= Include /kloe/soft/off/offline/inc/development/tls/maxstructdim.cin =
      INTEGER     NeleCluMax
      PARAMETER ( NeleCluMax = 2000 )
      INTEGER     MaxNumClu
      PARAMETER ( MaxNumClu = 100 )
      INTEGER    MaxNumVtx
      PARAMETER (MaxNumVtx = 20)
      INTEGER    MaxNumTrkV
      PARAMETER (MaxNumTrkV = 30)
      INTEGER    MaxNumTrk
      PARAMETER (MaxNumTrk = 100)
      INTEGER    MaxNumDHSP
      PARAMETER (MaxNumDHSP = 1000)
      Integer    nMaxDC
      Parameter (nMaxDC=1500)
      INTEGER    NQihiMax
      PARAMETER (NQihiMax=1000)
      INTEGER    NQcalMax
      PARAMETER (NQcalMax= 32 )
      INTEGER    MaxNumFirstHit
      PARAMETER (MaxNumFirstHit = 300 )
      INTEGER    MaxNtrkGen
      PARAMETER (MaxNtrkGen =50)
      INTEGER    MaxNvtxGen
      PARAMETER (MaxNvTxGen =30)
      integer TriggerElements
      parameter (TriggerElements = 300 )
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
         INTEGER iv(MaxNumtrkv)
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
         INTEGER mother(MaxNvtxGen)
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
C====================== Include interfstruct.inc =======================
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
          REAL    Tcl(MaxNumClu)
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
          INTEGER iv(MaxNumtrkv)
          INTEGER trknumv(MaxNumtrkv)
          REAL    CurV(MaxNumtrkv)
          REAL    PhiV(MaxNumtrkv)
          REAL    CotV(MaxNumtrkv)
          REAL    PxTV(MaxNumtrkv)
          REAL    PyTV(MaxNumtrkv)
          REAL    PzTV(MaxNumtrkv)
          INTEGER nv
          INTEGER vtx(MaxNumVtx)
          REAL    xv(MaxNumVtx)
          REAL    yv(MaxNumVtx)
          REAL    zv(MaxNumVtx)
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
          INTEGER pidmc(MaxNtrkGen)
          INTEGER virmom(MaxNtrkGen)
          REAL    pxmc(MaxNtrkGen)
          REAL    pymc(MaxNtrkGen)
          REAL    pzmc(MaxNtrkGen)
          REAL    themc(MaxNtrkGen)
          REAL    phimc(MaxNtrkGen)
          INTEGER vtxmc(MaxNtrkGen)
          INTEGER nvtxmc
          INTEGER kinmom(MaxNvtxGen)
          INTEGER mother(MaxNvtxGen)
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
          REAL    KchMC(9)
          REAL    KneMC(9)
          REAL    KchRec(9)
          REAL    KneRec(9)
          REAL    KchBoost(9)
          REAL    ip(3)
          REAL    ipmc(3)
          INTEGER vtaken(3)
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
          REAL    cosTrk
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
          INTEGER omegatruth_top
          INTEGER clu_index(4)
          INTEGER chtrk_index(2)
          INTEGER chvtx_index
          REAL P4gam_w(4,4)
          REAL P4ch_w(2,4)
          REAL P4ch_tot_w(4)
          REAL om_dist(4)
      END TYPE
      TYPE ( interfstru) Interf
      COMMON / InterfCommon / Interf,Evt
C-----------------------------------------------------------------------
C
C External functions
C
      INTEGER GETECLS,GETEVCL,FillBPOSCommon,GETTIME,BLOCAT
      INTEGER GETCLUSTRU,GETCLUSTRUCORR,Recover_Splitting,getbpos
      INTEGER gettclostru,trkv2stru,getmcstru,T0MCRD,t0glrd,geteclo
C
C Local declarations
C
      INTEGER istat,ind1,Ind,inddat,NrunBPOS,NrunBPOSold
      INTEGER  TimeSec,TimeMusec,mcflag_int,rstat
      LOGICAL enecorflag
      REAL    T0G_PHASED_DC, T0G_HIT_STEP0
      REAL    T0G_CLUSTER, T0G_STEP1, tphased_mc
      REAL    DELTA_CAVI_CALO, BUNCH_PERIOD
      enecorflag=.false.
C-----------------------------------------------------------------------
C Fill BPOS common
C-----------------------------------------------------------------------
      NrunBPOS = Nrun
      mcflag_int = 0
C MC
      IF( mcflag )THEN                 
         mcflag_int = 1
         IStat = BLOCAT(IW,'BRID',1,Ind,IndDat)
         IF( IStat.EQ.YESUCC ) THEN
            NrunBPOS = IW(IndDat)
         ELSE
            CALL ERLOGR('PROD2NTU',ERWARN,0,IStat,
     &           'Unable locating BRID bank: using LRID info')
         ENDIF
      ENDIF
      IF( NrunBPOS.NE.NrunBPOSold) THEN
         NrunBPOSold = NrunBPOS
         IStat = FillBPOSCommon(NrunBPOS)
         IF( IStat.NE.0 )THEN
            CALL ERLOGR('PROD2NTU',EREROR,0,IStat,
     &           'Failed to load boost information for this run')
         ENDIF
      ENDIF
C Average values for X
      evt%bpos%x     = 0.   
C " "      " "   for Y
      evt%bpos%y     = 0.    
C " "      " "   for Z
      evt%bpos%z     = 0.   
C Error on determination of X
      evt%bpos%errx  = 0.  
      evt%bpos%erry  = 0.
      evt%bpos%errz  = 0.
C Luminous region  in X
      evt%bpos%lumx  = 0.   
C Error on Lregion in X
      evt%bpos%elumx = 0.   
C Luminous region  in Z
      evt%bpos%lumz  = 0.   
C Error on Lregion in Z
      evt%bpos%elumz = 0.  
C ----------  fill Beam Momentum -------------------------------
      evt%bpos%px    = 0.
      evt%bpos%py    = 0.
      evt%bpos%pz    = 0.
      evt%bpos%errpx = 0.
      evt%bpos%errpy = 0.
      evt%bpos%errpz = 0.
      evt%bpos%larpx = 0.
      evt%bpos%larpy = 0.
      evt%bpos%larpz = 0.
      evt%bpos%elarpx= 0.
      evt%bpos%elarpy= 0.
      evt%bpos%elarpz= 0.
      evt%bpos%ene    = 0.
      evt%bpos%errene = 0.
      evt%bpos%dum1   = 0.
      evt%bpos%errdum1= 0.
      evt%bpos%dum2   = 0.
      evt%bpos%ErrDum2= 0.
      rstat =  getbpos( evt%Bpos )
C-----------------------------------------------------------------------
C Fill event classification information
C-----------------------------------------------------------------------
      IF( mcflag )  THEN
         istat = GETEVCL(evt%Info)
         istat = BLOCAT(iw,'EVCL',1,ind1,inddat)
      ELSE
         istat = GETTIME(timesec,timemusec)
C Reduce the timing
         timesec = timesec-946080000  
      ENDIF
C-----------------------------------------------------------------------
C Fill Event Infos
C-----------------------------------------------------------------------
      Evt%Info%RunNumber   = Nrun
      Evt%Info%EventNumber = Nev
      istat = GETECLS(evt%ecls)
C-----------------------------------------------------------------------
C Fill Cluster Information
C-----------------------------------------------------------------------
C     istat = GETCLUSTRUCORR(evt%clu,enecorflag)
      istat = GETCLUSTRU(evt%clu)
      istat = Recover_Splitting(evt%clu)
C-----------------------------------------------------------------------
C Fill Montecarlo Information
C-----------------------------------------------------------------------
      IF( mcflag ) istat = GETMCSTRU(evt%mc)
C-----------------------------------------------------------------------
C Fill tracks/verticies Information
C-----------------------------------------------------------------------
      istat = TRKV2STRU(evt%vtx,evt%trkv,evt%trk,evt%trkmc)
C-----------------------------------------------------------------------
C Fill TCLO
C-----------------------------------------------------------------------
      istat =  GETTCLOSTRU(evt%tclo)
C-----------------------------------------------------------------------
C Fill T0GLOBAL
C-----------------------------------------------------------------------
      istat = T0GLRD( T0G_PHASED_DC, T0G_HIT_STEP0,
     &       T0G_CLUSTER, T0G_STEP1,DELTA_CAVI_CALO, BUNCH_PERIOD )
      Evt%T0stru%dc_step0   = T0G_PHASED_DC
      Evt%T0stru%hit_step0  = T0G_HIT_STEP0
      Evt%T0stru%clus_step0 = T0G_CLUSTER
      Evt%T0stru%step1      = T0G_STEP1
      Evt%T0stru%cable      = DELTA_CAVI_CALO
      Evt%T0stru%tbunch     = BUNCH_PERIOD
      IF( mcflag ) THEN
         istat = T0MCRD(tphased_mc)
         Evt%T0stru%tphased_mc = tphased_mc
      ELSE
         Evt%T0stru%tphased_mc = 0.
      ENDIF
C-----------------------------------------------------------------------
C Fill ECLO
C-----------------------------------------------------------------------
C     istat = GETECLO(evt%eclo)
C-----------------------------------------------------------------------
      Interf%nev        = Evt%Info%EventNumber
      Interf%pileup     = Evt%Info%pileup
      Interf%gcod       = Evt%Info%gencod
      Interf%phid       = Evt%Info%phidecay
      Interf%a1typ      = Evt%Info%a1type
      Interf%a2typ      = Evt%Info%a2type
      Interf%a3typ      = Evt%Info%a3type
      Interf%b1typ      = Evt%Info%b1type
      Interf%b2typ      = Evt%Info%b2type
      Interf%b3typ      = Evt%Info%b3type
      Interf%nrundata   = nrunbpos
      Interf%tphased_mc = Evt%T0stru%tphased_mc
      Interf%t0dc0      = Evt%T0stru%dc_step0
      Interf%t0hit0     = Evt%T0stru%hit_step0
      Interf%t0clu0     = Evt%T0stru%clus_step0
      Interf%T0step1    = Evt%T0stru%step1
      Interf%DelayCable = Evt%T0stru%cable
      Interf%Tbunch     = Evt%T0stru%tbunch
      Interf%TimeSec    = TimeSec
      Interf%TimeMusec  = TimeMusec
      Interf%mcflag     = mcflag_int
      Interf%Bpx        = Evt%Bpos%px
      Interf%Bpy        = Evt%Bpos%py
      Interf%Bpz        = Evt%Bpos%pz
      Interf%Bx         = Evt%Bpos%x
      Interf%By         = Evt%Bpos%y
      Interf%Bz         = Evt%Bpos%z
      Interf%Bwidpx     = Evt%Bpos%errpx
      Interf%Bwidpy     = Evt%Bpos%errpy
      Interf%Bwidpz     = Evt%Bpos%errpz
      Interf%Bsx        = Evt%Bpos%errx
      Interf%Bsy        = Evt%Bpos%erry
      Interf%Bsz        = Evt%Bpos%errz
      Interf%Blumx      = Evt%Bpos%lumx
      Interf%Blumz      = Evt%Bpos%lumz
      Interf%Broots     = Evt%Bpos%Ene
      Interf%BrootsErr  = Evt%Bpos%ErrEne
      Interf%necls2     = Evt%ecls%n2
      Interf%ECLtrgw2   = Evt%ecls%trigger2
      Interf%ECLfilfo2  = Evt%ecls%Filfo2
      Interf%ECLword2   = Evt%ecls%totword2
      Interf%ECLstream2 = Evt%ecls%stream2
      Interf%ECLtagnum2 = Evt%ecls%tagnum2
      Interf%ECLevtype2 = Evt%ecls%evntyp2
      Interf%nclu       = Evt%Clu%n
      Interf%EneCl      = Evt%Clu%E
      Interf%Tcl        = Evt%Clu%T
      Interf%Xcl        = Evt%Clu%X
      Interf%Ycl        = Evt%Clu%Y
      Interf%Zcl        = Evt%Clu%Z
      Interf%Xacl       = Evt%Clu%Xa
      Interf%Yacl       = Evt%Clu%Ya
      Interf%Zacl       = Evt%Clu%Za
      Interf%XRmCl      = Evt%Clu%Xrms
      Interf%YRmsCl     = Evt%Clu%Yrms
      Interf%ZrmsCl     = Evt%Clu%Zrms
      Interf%TrmsCl     = Evt%Clu%Trms
      Interf%FlagCl     = Evt%Clu%Flag
      Interf%nclumc     = Evt%Clu%nmc
      Interf%Npar       = Evt%Clu%Npart
      Interf%Pnum1      = Evt%Clu%part1
      Interf%Pid1       = Evt%Clu%pid1
      Interf%Pnum2      = Evt%Clu%part2
      Interf%Pid2       = Evt%Clu%pid2
      Interf%Pnum3      = Evt%Clu%part3
      Interf%Pid3       = Evt%Clu%pid3
      Interf%ntv        = Evt%Trkv%n
      Interf%iv         = Evt%Trkv%iv
      Interf%trknumv    = Evt%Trkv%trkpoi
      Interf%CurV       = Evt%Trkv%cur
      Interf%PhiV       = Evt%Trkv%phi
      Interf%CotV       = Evt%Trkv%cot
      Interf%PxTV       = Evt%Trkv%px
      Interf%PyTV       = Evt%Trkv%py
      Interf%PzTV       = Evt%Trkv%pz
      Interf%nv         = Evt%vtx%n
      Interf%vtx        = Evt%vtx%Ntrk
      Interf%xv         = Evt%vtx%X
      Interf%yv         = Evt%vtx%Y
      Interf%zv         = Evt%vtx%Z
      Interf%chivtx     = Evt%vtx%Chi2
      Interf%qualv      = Evt%vtx%Qual
      Interf%fitidv     = Evt%vtx%Fitid
      Interf%VTXcov1    = Evt%vtx%Cov1
      Interf%VTXcov2    = Evt%vtx%Cov2
      Interf%VTXcov3    = Evt%vtx%Cov3
      Interf%VTXcov4    = Evt%vtx%Cov4
      Interf%VTXcov5    = Evt%vtx%Cov5
      Interf%VTXcov6    = Evt%vtx%Cov6
      Interf%nt         = Evt%Trk%n
      Interf%trkind     = Evt%Trk%trkind
      Interf%chi2fit    = Evt%Trk%chi2fit
      Interf%chi2ms     = Evt%Trk%chi2ms
      Interf%ntfmc      = Evt%TrkMC%n
      Interf%trkine1    = Evt%TrkMC%kine1
      Interf%trtype1    = Evt%TrkMC%type1
      Interf%trhits1    = Evt%TrkMC%hits1
      Interf%trkine2    = Evt%TrkMC%kine2
      Interf%trtype2    = Evt%TrkMC%type2
      Interf%trhits2    = Evt%TrkMC%hits2
      Interf%trkine3    = Evt%TrkMC%kine3
      Interf%trtype3    = Evt%TrkMC%type3
      Interf%trhits3    = Evt%TrkMC%hits3
      Interf%ntmc       = Evt%MC%ntrk
      Interf%kine       = Evt%MC%kin
      Interf%pidmc      = Evt%MC%pid
      Interf%virmom     = Evt%MC%virmom
      Interf%pxmc       = Evt%MC%px
      Interf%pymc       = Evt%MC%py
      Interf%pzmc       = Evt%MC%pz
      Interf%themc      = Evt%MC%theta
      Interf%phimc      = Evt%MC%phi
      Interf%vtxmc      = Evt%MC%Indv
      Interf%nvtxmc     = Evt%MC%numvtx
      Interf%kinmom     = Evt%MC%kinmom
      Interf%mother     = Evt%MC%mother
      Interf%xvmc       = Evt%MC%xv
      Interf%yvmc       = Evt%MC%yv
      Interf%zvmc       = Evt%MC%zv
      Interf%ntvtx      = Evt%MC%trkvtx
      Interf%ntcl       = Evt%Tclo%nt
      Interf%Asstr      = Evt%Tclo%trknum
      Interf%Asscl      = Evt%Tclo%clunum
      Interf%verver     = Evt%Tclo%verver
      Interf%xext       = Evt%Tclo%xext
      Interf%yext       = Evt%Tclo%yext
      Interf%zext       = Evt%Tclo%zext
      Interf%Assleng    = Evt%Tclo%leng
      Interf%AssChi     = Evt%Tclo%chi
      Interf%extPx      = Evt%Tclo%px
      Interf%extPy      = Evt%Tclo%py
      Interf%extPz      = Evt%Tclo%pz
C     Interf%ncli2      = Evt%eclo%n2
C     Interf%ECLOword2  = Evt%eclo%totword2
C     Interf%idpart2    = Evt%eclo%idpart2
C     Interf%dtclpo2    = Evt%eclo%dtclpo2
C     Interf%dvvnpo2    = Evt%eclo%dvvnpo2
C     Interf%stre2      = Evt%eclo%stre2
C     Interf%algo2      = Evt%eclo%algo2
      END
C=======================================================================
C=======================================================================
C=======================================================================
      SUBROUTINE find_range(Min,Max,String)
C-----------------------------------------------------------------------
        IMPLICIT NONE
C-----------------------------------------------------------------------
      INTEGER Min
      INTEGER Max
      CHARACTER*(*) String
C-----------------------------------------------------------------------
      WRITE(String,'(A,I5,A,I5,A)') '[',Min,',',Max,']'
      END
C==============================================================================
C==============================================================================
C==============================================================================
      SUBROUTINE GetKslEvent(ntmc,mother,vtxmc,pidmc,xvmc,yvmc,zvmc,
     &           pxmc,pymc,pzmc,nvtxmc,ipmc,KchMC,KneMC,DtMC,DlMC,
     &           truth,truthreg,truthsemi,truththree,truthomega,truthels
     &e)
C-----------------------------------------------------------------------
C
C  Description:
C  ------------
C in: ntmc,mother(),vtxmc(),pidmc(),xvmc(),yvmc(),zvmc(),
C in: pxmc(),pymc(),pzmc(),nvtxmc
C out: truth,truthreg,truthsemi,truththree,truthomega,truthelse,ipmc(),KchMC(),K
C      neMC(),DtMC,DlMC                                                         
C
C-----------------------------------------------------------------------
C
        IMPLICIT NONE
C
C-----------------------------------------------------------------------
C= Include /kloe/soft/off/offline/inc/development/tls/maxstructdim.cin =
      INTEGER     NeleCluMax
      PARAMETER ( NeleCluMax = 2000 )
      INTEGER     MaxNumClu
      PARAMETER ( MaxNumClu = 100 )
      INTEGER    MaxNumVtx
      PARAMETER (MaxNumVtx = 20)
      INTEGER    MaxNumTrkV
      PARAMETER (MaxNumTrkV = 30)
      INTEGER    MaxNumTrk
      PARAMETER (MaxNumTrk = 100)
      INTEGER    MaxNumDHSP
      PARAMETER (MaxNumDHSP = 1000)
      Integer    nMaxDC
      Parameter (nMaxDC=1500)
      INTEGER    NQihiMax
      PARAMETER (NQihiMax=1000)
      INTEGER    NQcalMax
      PARAMETER (NQcalMax= 32 )
      INTEGER    MaxNumFirstHit
      PARAMETER (MaxNumFirstHit = 300 )
      INTEGER    MaxNtrkGen
      PARAMETER (MaxNtrkGen =50)
      INTEGER    MaxNvtxGen
      PARAMETER (MaxNvTxGen =30)
      integer TriggerElements
      parameter (TriggerElements = 300 )
C======================== Include constans.inc =========================
      REAL       Mpip
      PARAMETER( Mpip = 139.57 )             
      REAL       Mpio
      PARAMETER( Mpio = 134.98 )             
      REAL       Mko
      PARAMETER( Mko  = 497.61 )             
      REAL       Momega
      PARAMETER( Momega  = 782.65 )          
      REAL       Cvel
      PARAMETER( Cvel = 29.9792458 )         
      REAL       EMCvel
      PARAMETER( EMCvel = 28.17 )            
      REAL       TauKs
      PARAMETER( TauKs = 0.08953 )   
      REAL       Mphi
      PARAMETER( Mphi  = 1019.460 )          
C========================= Include klspm00.inc =========================
      INTEGER NCLMIN
      PARAMETER( NCLMIN = 4 )
      REAL ECLMIN
      PARAMETER( ECLMIN = 20. )
      REAL vtxdist
      PARAMETER( vtxdist = 50. )
      INTEGER testcounter
      PARAMETER( testcounter = 100 )
      INTEGER MaxNumFitPar
      PARAMETER( MaxNumFitPar = 38 )
      INTEGER MaxNumConstr
      PARAMETER( MaxNumConstr = 10 )
      INTEGER MaxNumComega
      PARAMETER( MaxNumComega = 8 )
      INTEGER cutN
      PARAMETER( cutN = 6 )
      INTEGER errN
      PARAMETER( errN = 17 )
      CHARACTER*40 messageCut(testcounter),messageErr(testcounter)
      LOGICAL    McFlag,EventOK,ErrFlag,truthmc,truthregg,truthomegaa,tr
     1  uthsemii,truththreee,truthelsee                                 
      INTEGER    nmctruth0,ntruth,nmctruth(7),ntruthregg,ntruthsemii,ntr
     1  uththreee,ntruthomegaa,ntruthelsee                              
      INTEGER    selected,analised
      INTEGER    counter(testcounter),counterBcg(testcounter)
      INTEGER    ErrFlagCount(testcounter),ErrFlagCountBcg(testcounter)
      LOGICAL    cutchar,cutneu,onlytruemc,makecuts
      COMMON / flagcommon / McFlag,EventOK,ErrFlag,truthmc
      COMMON / truthcommon / nmctruth0,ntruth,nmctruth
      COMMON / counterscommon / selected,analised,
     &                          counter,counterBcg
      COMMON / errcommon / ErrFlagCount,ErrFlagCountBcg,
     &                     messageCut,messageErr
      COMMON / talk / cutchar,cutneu,onlytruemc,makecuts
      SAVE /flagcommon/, /truthcommon/, /counterscommon/
      SAVE /errcommon/, /talk/
C-----------------------------------------------------------------------
C
C External functions
C
C
C Local declarations
C
      INTEGER NpioKs,NpipKs,NpimKs,NepKs,NemKs,NmupKs,NmumKs,NothKs,Npio
     1  Kl,NpipKl,NpimKl,NepKl,                                         
     &        NemKl,NmupKl,NmumKl,NothKl,Nreg,Nother,Npiow,Npipw,Npimw,N
     1  othw,NKS,NKL,Nisr,nvtxmc                                        
      INTEGER ntmc,mother(MaxNvtxGen),vtxmc(MaxNtrkGen),pidmc(MaxNtrkGen
     &)
      INTEGER KLpic(2),KLpio(2),KSpic(2),KSpio(2),i
      REAL    xvmc(MaxNvtxGen),yvmc(MaxNvtxGen),zvmc(MaxNvtxGen),Broots
      REAL    ipmc(3),KchMC(9),KneMC(9),Ks(9),Kl(9),DtMC,DlMC,PhiPMC(3)
      REAL    pxmc(MaxNtrkGen),pymc(MaxNtrkGen),pzmc(MaxNtrkGen)
      LOGICAL truth,truthreg,truthsemi,truththree,truthomega,truthelse
      LOGICAL condsig,condreg,condsemi,condthree,condomega,condelse
C
C-----------------------------------------------------------------------
C
      NpioKs = 0
      NpipKs = 0
      NpimKs = 0
      NepKs = 0
      NemKs = 0
      NmupKs = 0
      NmumKs = 0
      NothKs = 0
      NpioKl = 0
      NpipKl = 0
      NpimKl = 0
      NepKl = 0
      NemKl = 0
      NmupKl = 0
      NmumKl = 0
      NothKl = 0
      Nreg = 0
      Nother = 0
      Npiow = 0
      Npipw = 0
      Npimw = 0
      Nothw = 0
      NKS = 0
      NKL = 0
      Nisr = 0
      truth = .FALSE.
      truthreg = .FALSE.
      truthsemi = .FALSE.
      truththree = .FALSE.
      truthomega = .FALSE.
      truthelse = .FALSE.
C
      DO i = 1,ntmc
C Kl vertex
         IF( mother(vtxmc(i)).EQ.10 )THEN       
C pi0
            IF( pidmc(i).EQ.7 )THEN             
               NpioKl = NpioKl + 1
C pi+
            ELSEIF( pidmc(i).EQ.8 ) THEN   
               NpipKl = NpipKl + 1
C pi-
            ELSEIF( pidmc(i).EQ.9 ) THEN   
               NpimKl = NpimKl + 1
C e+
            ELSEIF( pidmc(i).EQ.2 ) THEN   
               NepKl = NepKl + 1
C e-
            ELSEIF( pidmc(i).EQ.3 ) THEN   
               NemKl = NemKl + 1
C mu+
            ELSEIF( pidmc(i).EQ.5 ) THEN   
               NmupKl = NmupKl + 1
C mu-
            ELSEIF( pidmc(i).EQ.6 ) THEN   
               NmumKl = NmumKl + 1
C regeneration
            ELSEIF( pidmc(i).EQ.16 ) THEN   
               Nreg = Nreg + 1
            ELSE
               NothKl = NothKl + 1
            ENDIF
C Ks vertex
         ELSEIF( mother(vtxmc(i)).EQ.16 )THEN   
C pi0
            IF( pidmc(i).EQ.7 )THEN             
               NpioKs = NpioKs + 1
C pi+
            ELSEIF( pidmc(i).EQ.8 ) THEN   
               NpipKs = NpipKs + 1
C pi-
            ELSEIF( pidmc(i).EQ.9 ) THEN   
               NpimKs = NpimKs + 1
C e+
            ELSEIF( pidmc(i).EQ.2 ) THEN   
               NepKs = NepKs + 1
C e-
            ELSEIF( pidmc(i).EQ.3 ) THEN   
               NemKs = NemKs + 1
C mu+
            ELSEIF( pidmc(i).EQ.5 ) THEN   
               NmupKs = NmupKs + 1
C mu-
            ELSEIF( pidmc(i).EQ.6 ) THEN   
               NmumKs = NmumKs + 1
            ELSE
               NothKs = NothKs + 1
            ENDIF
C Omega vertex
         ELSEIF( mother(vtxmc(i)).EQ.50 )THEN   
C pi0
            IF( pidmc(i).EQ.7 )THEN             
               Npiow = Npiow + 1
C pi+
            ELSEIF( pidmc(i).EQ.8 ) THEN   
               Npipw = Npipw + 1
C pi-
            ELSEIF( pidmc(i).EQ.9 ) THEN   
               Npimw = Npimw + 1
C Kl
            ELSEIF( pidmc(i).EQ.10 ) THEN   
               NKL = NKL + 1
C Ks
            ELSEIF( pidmc(i).EQ.16 ) THEN   
               NKS = NKS + 1
            ELSEIF( pidmc(i).EQ.1 ) THEN
               Nisr = Nisr + 1
            ELSE
               Nothw = Nothw + 1
            ENDIF
         ELSE
            Nother = Nother + 1
         ENDIF
      ENDDO
      condsig = ( ( NpioKs.EQ.2 .AND. NpipKl.EQ.1 .AND. NpimKl.EQ.1 .AND
     &.
     &         NpipKs.EQ.0 .AND. NpimKs.EQ.0 .AND. NpioKl.EQ.0 ).OR.
     &       ( NpipKs.EQ.1 .AND. NpimKs.EQ.1 .AND. NpioKl.EQ.2 .AND.
     &         NpioKs.EQ.0 .AND. NpipKl.EQ.0 .AND. NpimKl.EQ.0  ) ) .AND
     1  . NKS.EQ.1 .AND. NKL.EQ.1 .AND. Nisr.GE.0                       
      condthree = ( ( NpioKs.EQ.3 .AND. NpipKl.EQ.1 .AND. NpimKl.EQ.1 .A
     &ND.
     &         NpipKs.EQ.0 .AND. NpimKs.EQ.0 .AND. NpioKl.EQ.0 ).OR.
     &       ( NpipKs.EQ.1 .AND. NpimKs.EQ.1 .AND. NpioKl.EQ.3 .AND.
     &         NpioKs.EQ.0 .AND. NpipKl.EQ.0 .AND. NpimKl.EQ.0 ) ) .AND.
     1   NKS.EQ.1 .AND. NKL.EQ.1 .AND. Nisr.GE.0                        
      condsemi = ( ( NpioKs.EQ.2 .AND. NpipKl.EQ.1 .AND. (NmumKl.EQ.1 .O
     1  R. NemKl.EQ.1) .AND.                                            
     &         NpipKs.EQ.0 .AND. NmumKs.EQ.0 .AND. NemKs.EQ.0 .AND. Npio
     1  Kl.EQ.0 ).OR.                                                   
     &       ( NpipKs.EQ.1 .AND. (NmumKs.EQ.1 .OR. NemKs.EQ.1) .AND. Npi
     1  oKl.EQ.2 .AND.                                                  
     &         NpioKs.EQ.0 .AND. NpipKl.EQ.0 .AND. NmumKl.EQ.0 .AND. Nem
     1  Kl.EQ.0 ).OR.                                                   
     &       ( NpioKs.EQ.2 .AND. NpimKl.EQ.1 .AND. (NmupKl.EQ.1 .OR. Nep
     1  Kl.EQ.1) .AND.                                                  
     &         NpimKs.EQ.0 .AND. NmupKs.EQ.0 .AND. NepKs.EQ.0 .AND. Npio
     1  Kl.EQ.0 ).OR.                                                   
     &       ( NpimKs.EQ.1 .AND. (NmupKs.EQ.1 .OR. NepKs.EQ.1) .AND. Npi
     1  oKl.EQ.2 .AND.                                                  
     &         NpioKs.EQ.0 .AND. NpimKl.EQ.0 .AND. NmupKl.EQ.0 .AND. Nep
     1  Kl.EQ.0 ) ) .AND. NKS.EQ.1 .AND. NKL.EQ.1 .AND. Nisr.GE.0       
      condomega = ( Npiow.EQ.2 .AND. Npipw.EQ.1 .AND. Npimw.EQ.1 .AND. N
     1  isr.GE.0)                                                       
      condreg = ( Nreg.EQ.1 .AND. NKS.EQ.1 .AND. NKL.EQ.1 .AND. Nisr.GE.
     &0 )
      IF( NothKs.EQ.0 .AND. NothKl.EQ.0 .AND. NmumKs.EQ.0 .AND. NmupKs.E
     1  Q.0 .AND. NepKs.EQ.0 .AND. NemKs.EQ.0 .AND.                     
     &    NmumKl.EQ.0 .AND. NmupKl.EQ.0 .AND. NepKl.EQ.0 .AND. NemKl.EQ.
     1  0 .AND. Nreg.EQ.0 .AND.                                         
     &    Npiow.EQ.0 .AND. Npipw.EQ.0 .AND. Npimw.EQ.0 .AND. Nothw.EQ.0)
     & THEN
         IF( condsig  ) THEN
                truth = .TRUE.
         ELSE IF( condthree  ) THEN
                truththree = .TRUE.
         ENDIF
      ELSE IF( NothKs.EQ.0 .AND. NothKl.EQ.0 .AND. Nreg.EQ.0 .AND. Npiow
     1  .EQ.0 .AND.                                                     
     &    Npipw.EQ.0 .AND. Npimw.EQ.0 .AND. Nothw.EQ.0 ) THEN
         IF( condsemi  ) truthsemi = .TRUE.
      ELSE IF( NothKs.EQ.0 .AND. NothKl.EQ.0 .AND. NmumKs.EQ.0 .AND. Nmu
     1  pKs.EQ.0 .AND. NepKs.EQ.0 .AND. NemKs.EQ.0 .AND.                
     &    NmumKl.EQ.0 .AND. NmupKl.EQ.0 .AND. NepKl.EQ.0 .AND. NemKl.EQ.
     1  0 .AND. Nreg.EQ.0 .AND.                                         
     &    Nothw.EQ.0 .AND. NKL.EQ.0 .AND. NKS.EQ.0) THEN
         IF( condomega  ) truthomega = .TRUE.
      ELSE IF( Npiow.EQ.0 .AND. Npipw.EQ.0 .AND. Npimw.EQ.0 .AND. Nothw.
     1  EQ.0 ) THEN                                                     
         IF( condreg  ) truthreg = .TRUE.
      ELSE
         truthelse = .TRUE.
      ENDIF
      IF( truth ) THEN
        DO i = 1,nvtxmc
           IF( mother(i).eq.50 ) THEN
              ipmc(1) = xvmc(i)
              ipmc(2) = yvmc(i)
              ipmc(3) = zvmc(i)
              EXIT
           ENDIF
        ENDDO
        DO i = 1,nvtxmc
C KL
           IF( mother(i).eq.10 ) THEN 
              Kl(7)=xvmc(i)
              Kl(8)=yvmc(i)
              Kl(9)=zvmc(i)
              EXIT
           ENDIF
        ENDDO
        DO i = 1,nvtxmc
C KS
           IF( mother(i).eq.16 ) THEN 
              Ks(7)=xvmc(i)
              Ks(8)=yvmc(i)
              Ks(9)=zvmc(i)
              EXIT
           ENDIF
        ENDDO
       DO i = 1,ntmc
CKL
           IF( pidmc(i).eq.10 ) THEN 
              Kl(1)=pxmc(i)
              Kl(2)=pymc(i)
              Kl(3)=pzmc(i)
              Kl(6)=Mko
              Kl(5)=Kl(1)**2 + Kl(2)**2 + Kl(3)**2
              Kl(4)=SQRT( Kl(5)+Kl(6)**2 )
              Kl(5)=SQRT( Kl(5) )
              EXIT
           ENDIF
        ENDDO
        DO i = 1,ntmc
CKS
           IF( pidmc(i).eq.16 ) THEN 
              Ks(1)=pxmc(i)
              Ks(2)=pymc(i)
              Ks(3)=pzmc(i)
              Ks(6)=Mko
              Ks(5)=Ks(1)**2 + Ks(2)**2 + Ks(3)**2
              Ks(4)=SQRT( Ks(5)+Ks(6)**2 )
              Ks(5)=SQRT( Ks(5) )
              EXIT
           ENDIF
        ENDDO
        DO i = 1,ntmc
C Kl vertex
           IF( mother(vtxmc(i)).EQ.10 ) THEN      
C pi0
              IF( pidmc(i).EQ.7 ) THEN             
                 KneMC=Kl
                 KchMC=Ks
                 EXIT
              ELSEIF( pidmc(i).EQ.8 .OR. pidmc(i).EQ.9 ) THEN
                 KchMC=Kl
                 KneMC=Ks
                 EXIT
              ENDIF
           ENDIF
        ENDDO
        DO i = 1,3
          PhiPMC(i) = KchMC(i)+KneMC(i)
        ENDDO
        Broots = KchMC(4)+KneMC(4)
        CALL GetTdiff(KchMC,KneMC,ipmc,PhiPMC,Broots,DlMC,DtMC)
      ENDIF
      END
C=======================================================================
C=======================================================================
C=======================================================================
      SUBROUTINE find_kchrec(qualv,nv,ntv,IV,PxTV,PyTV,PzTV,
     &                       xv,yv,zv,vtaken,KchRec,
     &                       trk1,trk2,cosTrk)
C-----------------------------------------------------------------------
        IMPLICIT NONE
C-----------------------------------------------------------------------
C========================= Include klspm00.inc =========================
      INTEGER NCLMIN
      PARAMETER( NCLMIN = 4 )
      REAL ECLMIN
      PARAMETER( ECLMIN = 20. )
      REAL vtxdist
      PARAMETER( vtxdist = 50. )
      INTEGER testcounter
      PARAMETER( testcounter = 100 )
      INTEGER MaxNumFitPar
      PARAMETER( MaxNumFitPar = 38 )
      INTEGER MaxNumConstr
      PARAMETER( MaxNumConstr = 10 )
      INTEGER MaxNumComega
      PARAMETER( MaxNumComega = 8 )
      INTEGER cutN
      PARAMETER( cutN = 6 )
      INTEGER errN
      PARAMETER( errN = 17 )
      CHARACTER*40 messageCut(testcounter),messageErr(testcounter)
      LOGICAL    McFlag,EventOK,ErrFlag,truthmc,truthregg,truthomegaa,tr
     1  uthsemii,truththreee,truthelsee                                 
      INTEGER    nmctruth0,ntruth,nmctruth(7),ntruthregg,ntruthsemii,ntr
     1  uththreee,ntruthomegaa,ntruthelsee                              
      INTEGER    selected,analised
      INTEGER    counter(testcounter),counterBcg(testcounter)
      INTEGER    ErrFlagCount(testcounter),ErrFlagCountBcg(testcounter)
      LOGICAL    cutchar,cutneu,onlytruemc,makecuts
      COMMON / flagcommon / McFlag,EventOK,ErrFlag,truthmc
      COMMON / truthcommon / nmctruth0,ntruth,nmctruth
      COMMON / counterscommon / selected,analised,
     &                          counter,counterBcg
      COMMON / errcommon / ErrFlagCount,ErrFlagCountBcg,
     &                     messageCut,messageErr
      COMMON / talk / cutchar,cutneu,onlytruemc,makecuts
      SAVE /flagcommon/, /truthcommon/, /counterscommon/
      SAVE /errcommon/, /talk/
C======================== Include constans.inc =========================
      REAL       Mpip
      PARAMETER( Mpip = 139.57 )             
      REAL       Mpio
      PARAMETER( Mpio = 134.98 )             
      REAL       Mko
      PARAMETER( Mko  = 497.61 )             
      REAL       Momega
      PARAMETER( Momega  = 782.65 )          
      REAL       Cvel
      PARAMETER( Cvel = 29.9792458 )         
      REAL       EMCvel
      PARAMETER( EMCvel = 28.17 )            
      REAL       TauKs
      PARAMETER( TauKs = 0.08953 )   
      REAL       Mphi
      PARAMETER( Mphi  = 1019.460 )          
C= Include /kloe/soft/off/offline/inc/development/tls/maxstructdim.cin =
      INTEGER     NeleCluMax
      PARAMETER ( NeleCluMax = 2000 )
      INTEGER     MaxNumClu
      PARAMETER ( MaxNumClu = 100 )
      INTEGER    MaxNumVtx
      PARAMETER (MaxNumVtx = 20)
      INTEGER    MaxNumTrkV
      PARAMETER (MaxNumTrkV = 30)
      INTEGER    MaxNumTrk
      PARAMETER (MaxNumTrk = 100)
      INTEGER    MaxNumDHSP
      PARAMETER (MaxNumDHSP = 1000)
      Integer    nMaxDC
      Parameter (nMaxDC=1500)
      INTEGER    NQihiMax
      PARAMETER (NQihiMax=1000)
      INTEGER    NQcalMax
      PARAMETER (NQcalMax= 32 )
      INTEGER    MaxNumFirstHit
      PARAMETER (MaxNumFirstHit = 300 )
      INTEGER    MaxNtrkGen
      PARAMETER (MaxNtrkGen =50)
      INTEGER    MaxNvtxGen
      PARAMETER (MaxNvTxGen =30)
      integer TriggerElements
      parameter (TriggerElements = 300 )
C-----------------------------------------------------------------------
C in:qualv(),nv,ntv,IV(),xv(),yv(),zv()
C in:CurV(),PhiV(),CotV() or PxTV(),PyTV(),PzTV()
C out:vtaken(),Kchrec(),trk1(),trk2(),cosTrk
C 1 px
C 2 py
C 3 pz
C 4 E
C 5 |p|
C 6 minv
C 7 xv
C 8 yv
C 9 zv
C Local variables
      INTEGER vtaken(3),i,j,k,l,nv,ntv,iv(MaxNumtrkv),ii
      REAL    KchRec(9),KchTmp(9),P4trkRec1(4),P4trkRec2(4)
      REAL    CurV(MaxNumtrkv),PhiV(MaxNumtrkv),CotV(MaxNumtrkv)
      REAL    PxTV(MaxNumtrkv),PyTV(MaxNumtrkv),PzTV(MaxNumtrkv)
      REAL    xv(MaxNumVtx),yv(MaxNumVtx),zv(MaxNumVtx)
      INTEGER qualv(MaxNumVtx)
      REAL    P4trkRec1mod,P4trkRec2mod
      REAL    trk1(4),trk2(4),cosTrk
C --- rec
      DO i = 1,3
         vtaken(i)=0
      ENDDO
      KchRec(6) = 100000000.
      DO i = 1,nv
C        IF( qualv(i).eq.1 ) THEN
            DO j = 1,ntv-1
C              IF( IV(j).eq.i.and.curv(j).ne.0 ) THEN
               IF( IV(j).eq.i ) THEN
C                 P4trkRec1(1) = COS( PhiV(j) ) / ABS( CurV(j) )*1000.
C                 P4trkRec1(2) = SIN( PhiV(j) ) / ABS( CurV(j) )*1000.
C                 P4trkRec1(3) =      CotV(j)   / ABS( CurV(j) )*1000.
                  P4trkRec1(1) = PxTV(j)
                  P4trkRec1(2) = PyTV(j)
                  P4trkRec1(3) = PzTV(j)
                  P4trkRec1mod = P4trkRec1(1)**2 + P4trkRec1(2)**2 +
     &                           P4trkRec1(3)**2
                  P4trkRec1(4) = SQRT( P4trkRec1mod + Mpip**2 )
                  P4trkRec1mod = SQRT( P4trkRec1mod )
                  DO k = j+1,ntv
C                    IF( IV(k).eq.i.and.curv(k).ne.0 ) THEN
                     IF( IV(k).eq.i ) THEN
C                       P4trkRec2(1) = COS( PhiV(k) ) / ABS( CurV(k) )*
C    &                                                   1000.
C                       P4trkRec2(2) = SIN( PhiV(k) ) / ABS( CurV(k) )*
C    &                                                   1000.
C                       P4trkRec2(3) =      CotV(k)   / ABS( CurV(k) )*
C    &                                                   1000.
                        P4trkRec2(1) = PxTV(k)
                        P4trkRec2(2) = PyTV(k)
                        P4trkRec2(3) = PzTV(k)
                        P4trkRec2mod = P4trkRec2(1)**2 + P4trkRec2(2)**2
     & +
     &                                 P4trkRec2(3)**2
                        P4trkRec2(4) = SQRT( P4trkRec2mod + Mpip**2 )
                        P4trkRec2mod = SQRT( P4trkRec2mod )
                        DO l = 1,4
                           KchTmp(l) = P4trkRec1(l) + P4trkRec2(l)
                        ENDDO
                        KchTmp(5) = KchTmp(1)**2 + KchTmp(2)**2 +
     &                              KchTmp(3)**2
                        KchTmp(6) = SQRT( KchTmp(4)**2 - KchTmp(5) )
                        KchTmp(5) = sqrt( KchTmp(5) )
                        IF( (vtaken(1).eq.0).or.
     &                      ( ABS( KchTmp(6) - Mko ).lt.
     &                        ABS( KchRec(6) - Mko ) ) ) THEN
C vtx id
                           vtaken(1) = i 
                           vtaken(2) = j
                           vtaken(3) = k
                           DO ii = 1,6
                              KchRec(ii) = KchTmp(ii)
                           ENDDO
                           DO ii = 1,4
                              trk1(ii) = P4trkRec1(ii)
                              trk2(ii) = P4trkRec2(ii)
                           ENDDO
                           cosTrk = ( trk1(1)*trk2(1)+
     &                                trk1(2)*trk2(2)+
     &                                trk1(3)*trk2(3) ) /
     &                                (P4trkRec1mod*P4trkRec2mod)
                        ENDIF
                     ENDIF
                  ENDDO
               ENDIF
            ENDDO
C        ENDIF
      ENDDO
      IF( vtaken(1).ne.0 ) THEN
         KchRec(7) = xv(vtaken(1))
         KchRec(8) = yv(vtaken(1))
         KchRec(9) = zv(vtaken(1))
      ELSE
         ErrFlag = .TRUE.
      ENDIF
      END
C=======================================================================
C=======================================================================
C=======================================================================
      SUBROUTINE cor_ip_boost(KchRec,Bpx,Bpy,Bpz,Bx,By,Bz,Broots,
     &                        trk1,trk2,KchBoost,ip,chdist,Qmiss)
C-----------------------------------------------------------------------
        IMPLICIT NONE
C-----------------------------------------------------------------------
C========================= Include klspm00.inc =========================
      INTEGER NCLMIN
      PARAMETER( NCLMIN = 4 )
      REAL ECLMIN
      PARAMETER( ECLMIN = 20. )
      REAL vtxdist
      PARAMETER( vtxdist = 50. )
      INTEGER testcounter
      PARAMETER( testcounter = 100 )
      INTEGER MaxNumFitPar
      PARAMETER( MaxNumFitPar = 38 )
      INTEGER MaxNumConstr
      PARAMETER( MaxNumConstr = 10 )
      INTEGER MaxNumComega
      PARAMETER( MaxNumComega = 8 )
      INTEGER cutN
      PARAMETER( cutN = 6 )
      INTEGER errN
      PARAMETER( errN = 17 )
      CHARACTER*40 messageCut(testcounter),messageErr(testcounter)
      LOGICAL    McFlag,EventOK,ErrFlag,truthmc,truthregg,truthomegaa,tr
     1  uthsemii,truththreee,truthelsee                                 
      INTEGER    nmctruth0,ntruth,nmctruth(7),ntruthregg,ntruthsemii,ntr
     1  uththreee,ntruthomegaa,ntruthelsee                              
      INTEGER    selected,analised
      INTEGER    counter(testcounter),counterBcg(testcounter)
      INTEGER    ErrFlagCount(testcounter),ErrFlagCountBcg(testcounter)
      LOGICAL    cutchar,cutneu,onlytruemc,makecuts
      COMMON / flagcommon / McFlag,EventOK,ErrFlag,truthmc
      COMMON / truthcommon / nmctruth0,ntruth,nmctruth
      COMMON / counterscommon / selected,analised,
     &                          counter,counterBcg
      COMMON / errcommon / ErrFlagCount,ErrFlagCountBcg,
     &                     messageCut,messageErr
      COMMON / talk / cutchar,cutneu,onlytruemc,makecuts
      SAVE /flagcommon/, /truthcommon/, /counterscommon/
      SAVE /errcommon/, /talk/
C= Include /kloe/soft/off/offline/inc/development/tls/maxstructdim.cin =
      INTEGER     NeleCluMax
      PARAMETER ( NeleCluMax = 2000 )
      INTEGER     MaxNumClu
      PARAMETER ( MaxNumClu = 100 )
      INTEGER    MaxNumVtx
      PARAMETER (MaxNumVtx = 20)
      INTEGER    MaxNumTrkV
      PARAMETER (MaxNumTrkV = 30)
      INTEGER    MaxNumTrk
      PARAMETER (MaxNumTrk = 100)
      INTEGER    MaxNumDHSP
      PARAMETER (MaxNumDHSP = 1000)
      Integer    nMaxDC
      Parameter (nMaxDC=1500)
      INTEGER    NQihiMax
      PARAMETER (NQihiMax=1000)
      INTEGER    NQcalMax
      PARAMETER (NQcalMax= 32 )
      INTEGER    MaxNumFirstHit
      PARAMETER (MaxNumFirstHit = 300 )
      INTEGER    MaxNtrkGen
      PARAMETER (MaxNtrkGen =50)
      INTEGER    MaxNvtxGen
      PARAMETER (MaxNvTxGen =30)
      integer TriggerElements
      parameter (TriggerElements = 300 )
C-----------------------------------------------------------------------
C External function
      REAL pk_from_boost
      INTEGER dist_lines
C 1 px
C 2 py
C 3 pz
C 4 E
C 5 |p|
C 6 minv
C 7 xv
C 8 yv
C 9 zv
C in: KchRec(),Bpx,Bpy,Bpz,Bx,By,Bz,Broots,trk1(),trk2()
C out: KchBoost(),ip,chdist,Qmiss
C Local variables
      INTEGER i,status
      REAL KchRec(9),KchBoost(9),ip(3),chdist,P4KchBoost(4)
      REAL P4PhiBha(4),pboost,xB(3),x(3),beamLine(3),p(3),pktP(3)
      REAL pktB(3),Bpx,Bpy,Bpz,Bx,By,Bz,Broots,trk1(4),trk2(4)
      REAL Emiss,Pmiss,Qmiss
C-----------------------------------------------------------------------
C --- Boost
      DO i = 1,4
         P4KchBoost(i) = KchRec(i)
      ENDDO
      P4PhiBha(1) = Bpx
      P4PhiBha(2) = Bpy
      P4PhiBha(3) = Bpz
      P4PhiBha(4) = Broots
      pboost = pk_from_boost(P4KchBoost,P4PhiBha)
      DO i = 1,4
         KchBoost(i) = P4KchBoost(i)
      ENDDO
      KchBoost(5) = pboost
      KchBoost(6) = SQRT( KchBoost(4)**2 - KchBoost(5)**2 )
      DO i = 7,9
         KchBoost(i) = KchRec(i)
      ENDDO
C --- IP
      xB(1) = Bx
      xB(2) = By
      xB(3) = Bz
      x(1) = KchBoost(7)
      x(2) = KchBoost(8)
      x(3) = KchBoost(9)
      BeamLine(1) = 0.
      BeamLine(2) = 0.
      BeamLine(3) = 1.
      DO i = 1,3
         p(i) = KchBoost(i)
      ENDDO
      status = dist_lines(x,p,xB,BeamLine,pktP,pktB,chdist)
C     integer function dist_lines (x1, p1, x2, p2, x1PCA, x2PCA, dist)
C     x1PCA     point on line 1 of closest approach to line 2
C     x2PCA     point on line 2 of closest approach to line 1
C     dist      distance of closest approach
      IF( status.eq.0) THEN
         ErrFlag = .TRUE.
         RETURN
      ENDIF
      ip(1) = xB(1)
      ip(2) = xB(2)
      ip(3) = pktB(3)
      Emiss = Broots/2 - trk1(4) - trk2(4)
      Pmiss = ( KchBoost(1)-trk1(1)-trk2(1) )**2 +
     &        ( KchBoost(2)-trk1(2)-trk2(2) )**2 +
     &        ( KchBoost(3)-trk1(3)-trk2(3) )**2
      Qmiss = sqrt( Emiss**2 + Pmiss )
      END
C=======================================================================
C=======================================================================
C=======================================================================
      REAL FUNCTION pk_from_boost(pK, pb)
C     Copied from Tommaso with inputs reconfigured and some
C     aesthetic changes.
C-----------------------------------------------------------------------
        IMPLICIT NONE
C     Global specifications
C======================== Include constans.inc =========================
      REAL       Mpip
      PARAMETER( Mpip = 139.57 )             
      REAL       Mpio
      PARAMETER( Mpio = 134.98 )             
      REAL       Mko
      PARAMETER( Mko  = 497.61 )             
      REAL       Momega
      PARAMETER( Momega  = 782.65 )          
      REAL       Cvel
      PARAMETER( Cvel = 29.9792458 )         
      REAL       EMCvel
      PARAMETER( EMCvel = 28.17 )            
      REAL       TauKs
      PARAMETER( TauKs = 0.08953 )   
      REAL       Mphi
      PARAMETER( Mphi  = 1019.460 )          
C     Argument specifications
      real pK(4)
      real pb(4)
C     Local specifications
      integer i
      real pb_mod, pK_mod, dot
      real cosb_sq, b_sq, bg_sq, g_sq, pcm_sq
      real C, A, B, disc
      real P1, P2
      real K0mass
C     Executable statements
      K0mass = Mko
      pk_from_boost = 0.
      pb_mod = 0.
      pK_mod = 0.
      dot = 0.
      do i = 1, 3
        pb_mod = pb_mod + pb(i)**2
        pK_mod = pK_mod + pK(i)**2
        dot = dot + pb(i)*pK(i)
      end do
      cosb_sq = dot**2/(pb_mod*pK_mod)
      pb_mod = sqrt(pb_mod)
      pK_mod = sqrt(pK_mod)
      b_sq = (pb_mod/pb(4))**2
      bg_sq = pb_mod**2/(pb(4)**2 - pb_mod**2)
      g_sq = 1. + bg_sq
      pcm_sq = 0.25*(pb(4)**2 - pb_mod**2) - k0mass**2
      C = (bg_sq*k0mass**2 - pcm_sq)**2
      A = (g_sq*(1. - b_sq*cosb_sq))**2
      B = g_sq *
     &    ((1. + b_sq*cosb_sq)*(bg_sq*k0mass**2 - pcm_sq) -
     &    2.*bg_sq*k0mass**2*cosb_sq)
      disc = B**2 - A*C
      if (disc.lt.0) then
        if (
     &      (disc.lt.-100.) ) then
            WRITE(*,*)'Large negative discriminant, suspect solution'
        end if
        disc = 0.
      end if
      disc = sqrt(disc)
      P1 = (-B + disc)/A
      P2 = (-B - disc)/A
      if (P1.gt.0.and.P2.lt.0) then
        pk_from_boost = sqrt(P1)
      else if (P2.gt.0.and.P1.lt.0) then
        pk_from_boost = sqrt(P2)
      else if (P1.gt.0.and.P2.gt.0) then
C       For boost directed along the NEGATIVE x axis
        if (dot.lt.0) then
          pk_from_boost = sqrt(min(P1, P2))
        else
          pk_from_boost = sqrt(max(P1, P2))
        end if
      end if
C     Fix-up the KS momentum on output
      do i = 1, 3
        pK(i) = pK(i)*pk_from_boost/pK_mod
      end do
      pK(4) = sqrt(pk_from_boost**2 + k0mass**2)
      return
      end
C=======================================================================
C=======================================================================
C=======================================================================
      INTEGER FUNCTION dist_lines (x1, p1, x2, p2, x1PCA, x2PCA, dist)
C     x1PCA     point on line 1 of closest approach to line 2
C     x2PCA     point on line 2 of closest approach to line 1
C     dist      distance of closest approach
        IMPLICIT NONE
C     Argument specifications
      real x1(3), p1(3), x2(3), p2(3)
      real x1PCA(3), x2PCA(3), dist
C     Local specifications
      integer k
      real p1mod, p2mod, dp1(3), dp2(3)
      real D(3), u(3), D1(3), D2(3), umod, D1u, D2u
C     Executable statements
      dist_lines = 0
      p1mod = sqrt(p1(1)**2 + p1(2)**2 + p1(3)**2)
      if (p1mod.eq.0.) goto 13
      p2mod = sqrt(p2(1)**2 + p2(2)**2 + p2(3)**2)
      if (p2mod.eq.0.) goto 13
      do k = 1, 3
        dp1(k) = p1(k)/p1mod
        dp2(k) = p2(k)/p2mod
        D(k) = x2(k) - x1(k)
      end do
      call crossProd (dp1, dp2, u)
      umod = sqrt(u(1)**2 + u(2)**2 + u(3)**2)
      if (umod.eq.0.) goto 13
      call crossProd (D, dp1, D1)
      D1u = D1(1)*u(1) + D1(2)*u(2) + D1(3)*u(3)
      call crossProd (D, dp2, D2)
      D2u = D2(1)*u(1) + D2(2)*u(2) + D2(3)*u(3)
      dist = 0.
      do k = 1, 3
        x1PCA(k) = x1(k) + (D2u/umod)*dp1(k)
        x2PCA(k) = x2(k) + (D1u/umod)*dp2(k)
        dist = dist + (x2PCA(k) - x1PCA(k))**2
      end do
      dist = sqrt(dist)
      dist_lines = 1
 13   return
      end
C=======================================================================
C=======================================================================
C=======================================================================
      SUBROUTINE crossProd (A, B, C)
        IMPLICIT NONE
      real A(3), B(3), C(3)
      C(1) = A(2)*B(3) - B(2)*A(3)
      C(2) = B(1)*A(3) - A(1)*B(3)
      C(3) = A(1)*B(2) - B(1)*A(2)
      return
      end
C=======================================================================
C=======================================================================
C=======================================================================
      SUBROUTINE when()
        IMPLICIT NONE
C
C Local declarations
C
      INTEGER today(3), now(3)
C-----------------------------------------------------------------------
C today(1)=day, (2)=month, (3)=year
      call idate(today) 
C now(1)=hour, (2)=minute, (3)=second
      call itime(now)   
         write( *, 998 ) today(1),today(2),today(3),
     &                    now(1),now(2),now(3)
 998     format ( '  ', i2.2, '.', i2.2, '.', i4.4, ' ',
     &                  i2.2, ':', i2.2, ':', i2.2 )
      end
C=======================================================================
C=======================================================================
C=======================================================================
      SUBROUTINE summary
        IMPLICIT NONE
C========================= Include klspm00.inc =========================
      INTEGER NCLMIN
      PARAMETER( NCLMIN = 4 )
      REAL ECLMIN
      PARAMETER( ECLMIN = 20. )
      REAL vtxdist
      PARAMETER( vtxdist = 50. )
      INTEGER testcounter
      PARAMETER( testcounter = 100 )
      INTEGER MaxNumFitPar
      PARAMETER( MaxNumFitPar = 38 )
      INTEGER MaxNumConstr
      PARAMETER( MaxNumConstr = 10 )
      INTEGER MaxNumComega
      PARAMETER( MaxNumComega = 8 )
      INTEGER cutN
      PARAMETER( cutN = 6 )
      INTEGER errN
      PARAMETER( errN = 17 )
      CHARACTER*40 messageCut(testcounter),messageErr(testcounter)
      LOGICAL    McFlag,EventOK,ErrFlag,truthmc,truthregg,truthomegaa,tr
     1  uthsemii,truththreee,truthelsee                                 
      INTEGER    nmctruth0,ntruth,nmctruth(7),ntruthregg,ntruthsemii,ntr
     1  uththreee,ntruthomegaa,ntruthelsee                              
      INTEGER    selected,analised
      INTEGER    counter(testcounter),counterBcg(testcounter)
      INTEGER    ErrFlagCount(testcounter),ErrFlagCountBcg(testcounter)
      LOGICAL    cutchar,cutneu,onlytruemc,makecuts
      COMMON / flagcommon / McFlag,EventOK,ErrFlag,truthmc
      COMMON / truthcommon / nmctruth0,ntruth,nmctruth
      COMMON / counterscommon / selected,analised,
     &                          counter,counterBcg
      COMMON / errcommon / ErrFlagCount,ErrFlagCountBcg,
     &                     messageCut,messageErr
      COMMON / talk / cutchar,cutneu,onlytruemc,makecuts
      SAVE /flagcommon/, /truthcommon/, /counterscommon/
      SAVE /errcommon/, /talk/
C-----------------------------------------------------------------------
C     INTEGER   sistle     <-- function to count letters in character variable
C------------------------------------------------------------------------------
      INTEGER i
      CALL when
      IF( mcflag ) THEN
         WRITE(*,999) analised,'analised including'
         WRITE(*,999) ntruth,'truth'
      ELSE
         WRITE(*,999) analised,'analised with'
      ENDIF
      WRITE(*,999) selected,'selected including:'
      IF( mcflag ) THEN
         WRITE(*,999) nmctruth(1),'truth cuts passed'
         WRITE(*,999) nmctruth(2),'truth rejected by cuts'
         WRITE(*,999) nmctruth0,  'truth with errors'
         IF( .not.onlytruemc ) THEN
           WRITE(*,999) nmctruth(3),'background'
         ENDIF
         WRITE(*,*) '   ErrFlg rejection TRUTH:'
         DO i=1,errN
           WRITE(*,999) ErrFlagCount(i),messageErr(i)
         ENDDO
         WRITE(*,*) '   Cuts rejection TRUTH:'
         DO i=1,cutN
           WRITE(*,999) counter(i),messageCut(i)
         ENDDO
      ENDIF
      IF( .not.(onlytruemc.and.mcflag) ) THEN
         WRITE(*,*) '  ErrFlg rejection BCG:'
         DO i=2,errN
           WRITE(*,999) ErrFlagCountBcg(i),messageErr(i)
         ENDDO
         WRITE(*,*) '   Cuts rejection BCG:'
         DO i=1,cutN
           WRITE(*,999) counterBcg(i),messageCut(i)
         ENDDO
      ENDIF
  999 format ('    ',i10,' ',a)
      WRITE(*,*)'----------------------------------------'
      END
C=======================================================================
C=======================================================================
C=======================================================================
      SUBROUTINE Loren4(Psys,Pold,Pnew)
C-----------------------------------------------------------------------
C
C Input:  Psys (1:4) 4-momentum of the SYSTEM
C                    in your starting frame (e.g. lab)
C                    1:3 = Px,Py,Pz, 4 = energy
C         Pold(1:4)  4-momentum of particle in old system
C
C Output: Pnew(1:4)  4-momentum of particle in new system
C
C-----------------------------------------------------------------------
C
        IMPLICIT NONE
C========================= Include klspm00.inc =========================
      INTEGER NCLMIN
      PARAMETER( NCLMIN = 4 )
      REAL ECLMIN
      PARAMETER( ECLMIN = 20. )
      REAL vtxdist
      PARAMETER( vtxdist = 50. )
      INTEGER testcounter
      PARAMETER( testcounter = 100 )
      INTEGER MaxNumFitPar
      PARAMETER( MaxNumFitPar = 38 )
      INTEGER MaxNumConstr
      PARAMETER( MaxNumConstr = 10 )
      INTEGER MaxNumComega
      PARAMETER( MaxNumComega = 8 )
      INTEGER cutN
      PARAMETER( cutN = 6 )
      INTEGER errN
      PARAMETER( errN = 17 )
      CHARACTER*40 messageCut(testcounter),messageErr(testcounter)
      LOGICAL    McFlag,EventOK,ErrFlag,truthmc,truthregg,truthomegaa,tr
     1  uthsemii,truththreee,truthelsee                                 
      INTEGER    nmctruth0,ntruth,nmctruth(7),ntruthregg,ntruthsemii,ntr
     1  uththreee,ntruthomegaa,ntruthelsee                              
      INTEGER    selected,analised
      INTEGER    counter(testcounter),counterBcg(testcounter)
      INTEGER    ErrFlagCount(testcounter),ErrFlagCountBcg(testcounter)
      LOGICAL    cutchar,cutneu,onlytruemc,makecuts
      COMMON / flagcommon / McFlag,EventOK,ErrFlag,truthmc
      COMMON / truthcommon / nmctruth0,ntruth,nmctruth
      COMMON / counterscommon / selected,analised,
     &                          counter,counterBcg
      COMMON / errcommon / ErrFlagCount,ErrFlagCountBcg,
     &                     messageCut,messageErr
      COMMON / talk / cutchar,cutneu,onlytruemc,makecuts
      SAVE /flagcommon/, /truthcommon/, /counterscommon/
      SAVE /errcommon/, /talk/
C
      REAL      Psys(4), Pold(4), Pnew(4)
      REAL      Esys, Pmod, Msys, ScalProd
C
C-----------------------------------------------------------------------
C
      CALL VZERO(Pnew,4)
C
      Esys = Psys(4)
      Pmod = SQRT( Psys(1)**2 + Psys(2)**2 + Psys(3)**2 )
C
      IF( Esys.LT.Pmod )THEN
         ErrFlag = .TRUE.
         RETURN
      ELSE
         Msys = SQRT( Esys**2 - Pmod**2 )
      ENDIF
C
      ScalProd = Pold(1)*Psys(1) + Pold(2)*Psys(2) +
     &     Pold(3)*Psys(3)
      Pnew(1) = Pold(1) + Psys(1)*(ScalProd/(Esys+Msys)-Pold(4))/Msys
      Pnew(2) = Pold(2) + Psys(2)*(ScalProd/(Esys+Msys)-Pold(4))/Msys
      Pnew(3) = Pold(3) + Psys(3)*(ScalProd/(Esys+Msys)-Pold(4))/Msys
C
C See CERNLIB U101 lorenb
C
      Pnew(4) = Pold(4)*Esys/Msys - ScalProd/Msys
C
      END
C=======================================================================
C=======================================================================
C=======================================================================
      SUBROUTINE check_isr(ntmc,pidmc,virmom,vtxmc,mother,isr)
        IMPLICIT NONE
C-----------------------------------------------------------------------
C= Include /kloe/soft/off/offline/inc/development/tls/maxstructdim.cin =
      INTEGER     NeleCluMax
      PARAMETER ( NeleCluMax = 2000 )
      INTEGER     MaxNumClu
      PARAMETER ( MaxNumClu = 100 )
      INTEGER    MaxNumVtx
      PARAMETER (MaxNumVtx = 20)
      INTEGER    MaxNumTrkV
      PARAMETER (MaxNumTrkV = 30)
      INTEGER    MaxNumTrk
      PARAMETER (MaxNumTrk = 100)
      INTEGER    MaxNumDHSP
      PARAMETER (MaxNumDHSP = 1000)
      Integer    nMaxDC
      Parameter (nMaxDC=1500)
      INTEGER    NQihiMax
      PARAMETER (NQihiMax=1000)
      INTEGER    NQcalMax
      PARAMETER (NQcalMax= 32 )
      INTEGER    MaxNumFirstHit
      PARAMETER (MaxNumFirstHit = 300 )
      INTEGER    MaxNtrkGen
      PARAMETER (MaxNtrkGen =50)
      INTEGER    MaxNvtxGen
      PARAMETER (MaxNvTxGen =30)
      integer TriggerElements
      parameter (TriggerElements = 300 )
C-----------------------------------------------------------------------
C in: ntmc,pidmc(),virmom(),vtxmc(),mother()
C out: isr
      INTEGER isr,i,ntmc,pidmc(MaxNtrkGen),virmom(MaxNtrkGen)
      INTEGER vtxmc(MaxNtrkGen),mother(MaxNvtxGen)
      ISR = 0
      DO i = 1,ntmc
          IF( pidmc(i).eq.1.and.virmom(i).eq.1.and.
     &        vtxmc(i).eq.1.and.
     &        mother(vtxmc(i)).eq.50 ) THEN
            ISR = 1
            EXIT
          ENDIF
      ENDDO
      END
C=======================================================================
C=======================================================================
C=======================================================================
      SUBROUTINE find_neuclu(nclu,ntcl,Asscl,
     &                       neuclulist,neucluind)
        IMPLICIT NONE
C-----------------------------------------------------------------------
C========================= Include klspm00.inc =========================
      INTEGER NCLMIN
      PARAMETER( NCLMIN = 4 )
      REAL ECLMIN
      PARAMETER( ECLMIN = 20. )
      REAL vtxdist
      PARAMETER( vtxdist = 50. )
      INTEGER testcounter
      PARAMETER( testcounter = 100 )
      INTEGER MaxNumFitPar
      PARAMETER( MaxNumFitPar = 38 )
      INTEGER MaxNumConstr
      PARAMETER( MaxNumConstr = 10 )
      INTEGER MaxNumComega
      PARAMETER( MaxNumComega = 8 )
      INTEGER cutN
      PARAMETER( cutN = 6 )
      INTEGER errN
      PARAMETER( errN = 17 )
      CHARACTER*40 messageCut(testcounter),messageErr(testcounter)
      LOGICAL    McFlag,EventOK,ErrFlag,truthmc,truthregg,truthomegaa,tr
     1  uthsemii,truththreee,truthelsee                                 
      INTEGER    nmctruth0,ntruth,nmctruth(7),ntruthregg,ntruthsemii,ntr
     1  uththreee,ntruthomegaa,ntruthelsee                              
      INTEGER    selected,analised
      INTEGER    counter(testcounter),counterBcg(testcounter)
      INTEGER    ErrFlagCount(testcounter),ErrFlagCountBcg(testcounter)
      LOGICAL    cutchar,cutneu,onlytruemc,makecuts
      COMMON / flagcommon / McFlag,EventOK,ErrFlag,truthmc
      COMMON / truthcommon / nmctruth0,ntruth,nmctruth
      COMMON / counterscommon / selected,analised,
     &                          counter,counterBcg
      COMMON / errcommon / ErrFlagCount,ErrFlagCountBcg,
     &                     messageCut,messageErr
      COMMON / talk / cutchar,cutneu,onlytruemc,makecuts
      SAVE /flagcommon/, /truthcommon/, /counterscommon/
      SAVE /errcommon/, /talk/
C= Include /kloe/soft/off/offline/inc/development/tls/maxstructdim.cin =
      INTEGER     NeleCluMax
      PARAMETER ( NeleCluMax = 2000 )
      INTEGER     MaxNumClu
      PARAMETER ( MaxNumClu = 100 )
      INTEGER    MaxNumVtx
      PARAMETER (MaxNumVtx = 20)
      INTEGER    MaxNumTrkV
      PARAMETER (MaxNumTrkV = 30)
      INTEGER    MaxNumTrk
      PARAMETER (MaxNumTrk = 100)
      INTEGER    MaxNumDHSP
      PARAMETER (MaxNumDHSP = 1000)
      Integer    nMaxDC
      Parameter (nMaxDC=1500)
      INTEGER    NQihiMax
      PARAMETER (NQihiMax=1000)
      INTEGER    NQcalMax
      PARAMETER (NQcalMax= 32 )
      INTEGER    MaxNumFirstHit
      PARAMETER (MaxNumFirstHit = 300 )
      INTEGER    MaxNtrkGen
      PARAMETER (MaxNtrkGen =50)
      INTEGER    MaxNvtxGen
      PARAMETER (MaxNvTxGen =30)
      integer TriggerElements
      parameter (TriggerElements = 300 )
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
C-----------------------------------------------------------------------
C in: nclu,ntcl,asscl()
C out: neuclulist,neucluind
      INTEGER neucluind,neuclulist(MaxNumCLu),i,j,nclu,ntcl
      INTEGER Asscl(NTCLOMax)
      LOGICAL neuclu
C-----------------------------------------------------------------------
      neucluind = 0
      DO i = 1,nclu
         neuclu = .TRUE.
         DO j = 1,ntcl
            IF( asscl(j).eq.i ) neuclu = .FALSE.
         ENDDO
         IF( neuclu ) THEN
            neucluind = neucluind + 1
            neuclulist(neucluind) = i
         ENDIF
      ENDDO
      IF( neucluind.lt.NCLMIN) ErrFlag = .TRUE.
      END
C=======================================================================
C=======================================================================
C=======================================================================
      SUBROUTINE GetTdiff(Kch,Kne,PhiVtx,PhiP,Broots,DL,DT)
C-----------------------------------------------------------------------
C
C Input variables:
C ----------------
C 1 px
C 2 py
C 3 pz
C 4 E
C 5 |p|
C 6 minv
C 7 xv
C 8 yv
C 9 zv
C Kch(9)  = 4-momentum of first kaon (pi+pi- final state)
C Kne(9)  = 4-momentum of second kaon (pi0pi0 final state)
C PhiVtx(3) = Phi vertex coordinate
C PhiP(3) = Phi momentum
C Broots = energy of the system
C
C Output variables:
C -----------------
C DL        = Difference of the decay length between pi+pi- and pi0pi0
C             final states
C DT        = Difference of the decay time between pi+pi- and pi0pi0
C             final states
C
C-----------------------------------------------------------------------
        IMPLICIT NONE
C-----------------------------------------------------------------------
C======================== Include constans.inc =========================
      REAL       Mpip
      PARAMETER( Mpip = 139.57 )             
      REAL       Mpio
      PARAMETER( Mpio = 134.98 )             
      REAL       Mko
      PARAMETER( Mko  = 497.61 )             
      REAL       Momega
      PARAMETER( Momega  = 782.65 )          
      REAL       Cvel
      PARAMETER( Cvel = 29.9792458 )         
      REAL       EMCvel
      PARAMETER( EMCvel = 28.17 )            
      REAL       TauKs
      PARAMETER( TauKs = 0.08953 )   
      REAL       Mphi
      PARAMETER( Mphi  = 1019.460 )          
C========================= Include klspm00.inc =========================
      INTEGER NCLMIN
      PARAMETER( NCLMIN = 4 )
      REAL ECLMIN
      PARAMETER( ECLMIN = 20. )
      REAL vtxdist
      PARAMETER( vtxdist = 50. )
      INTEGER testcounter
      PARAMETER( testcounter = 100 )
      INTEGER MaxNumFitPar
      PARAMETER( MaxNumFitPar = 38 )
      INTEGER MaxNumConstr
      PARAMETER( MaxNumConstr = 10 )
      INTEGER MaxNumComega
      PARAMETER( MaxNumComega = 8 )
      INTEGER cutN
      PARAMETER( cutN = 6 )
      INTEGER errN
      PARAMETER( errN = 17 )
      CHARACTER*40 messageCut(testcounter),messageErr(testcounter)
      LOGICAL    McFlag,EventOK,ErrFlag,truthmc,truthregg,truthomegaa,tr
     1  uthsemii,truththreee,truthelsee                                 
      INTEGER    nmctruth0,ntruth,nmctruth(7),ntruthregg,ntruthsemii,ntr
     1  uththreee,ntruthomegaa,ntruthelsee                              
      INTEGER    selected,analised
      INTEGER    counter(testcounter),counterBcg(testcounter)
      INTEGER    ErrFlagCount(testcounter),ErrFlagCountBcg(testcounter)
      LOGICAL    cutchar,cutneu,onlytruemc,makecuts
      COMMON / flagcommon / McFlag,EventOK,ErrFlag,truthmc
      COMMON / truthcommon / nmctruth0,ntruth,nmctruth
      COMMON / counterscommon / selected,analised,
     &                          counter,counterBcg
      COMMON / errcommon / ErrFlagCount,ErrFlagCountBcg,
     &                     messageCut,messageErr
      COMMON / talk / cutchar,cutneu,onlytruemc,makecuts
      SAVE /flagcommon/, /truthcommon/, /counterscommon/
      SAVE /errcommon/, /talk/
C-----------------------------------------------------------------------
C
      INTEGER    ILoop, i
      REAL       DL, DT, PhiVtx(3)
      REAL       P4Kch(4), P4Kne(4), P4KchNew(4), P4KneNew(4)
      REAL       Path4KchNew(4), Path4KneNew(4)
      REAL       P4Vph(4), P4VphNew(4), L1, L2, Time
      REAL       Kch(9),Kne(9),PhiP(3),Broots,Path4Kch(4),Path4Kne(4)
      REAL       Beta1, Beta2
C
C-----------------------------------------------------------------------
C
C
      DL = 999.
      DT = 999.
      DO i = 1,3
         P4Kch(i) = Kch(i)
         P4Kne(i) = Kne(i)
         Path4Kch(i) = Kch(i+6)
         Path4Kne(i) = Kne(i+6)
      ENDDO
      P4Kch(4) = Kch(4)
      P4Kne(4) = Kne(4)
      L1 = SQRT( (Path4Kch(1)-PhiVtx(1))**2 +
     &     (Path4Kch(2)-PhiVtx(2))**2 + (Path4Kch(3)-PhiVtx(3))**2 )
      Beta1 = SQRT( P4Kch(1)**2 + P4Kch(2)**2 + P4Kch(3)**2 ) / P4Kch(4)
      Path4Kch(4) = L1 / Beta1
      L2 = SQRT( (Path4Kne(1)-PhiVtx(1))**2 +
     &     (Path4Kne(2)-PhiVtx(2))**2 + (Path4Kne(3)-PhiVtx(3))**2 )
      Beta2 = SQRT( P4Kne(1)**2 + P4Kne(2)**2 + P4Kne(3)**2 ) / P4Kne(4)
      Path4Kne(4) = L2 / Beta2
      CALL Loren4(P4Kch,Path4Kch,Path4KchNew)
      IF ( ErrFlag ) RETURN
      CALL Loren4(P4Kne,Path4Kne,Path4KneNew)
      IF ( ErrFlag ) RETURN
      DT = ( Path4KchNew(4)-Path4KneNew(4) ) / ( Cvel * TauKs )
      DL = L1 - L2
      END
C=======================================================================
C=======================================================================
C=======================================================================
      SUBROUTINE find_neuvtx(Bpx,Bpy,Bpz,Broots,KchBoost,ncl,enecl,
     &           ncll,xcl,ycl,zcl,ip,tcl,cldist,KneRecLor,trc,
     &           nclwrong,ncllwrong,KneRec,minv4gam,pi0,
     &           g4taken,trcv,PgamRecTaken,gpairtaken,test,
     &           g4vtxerr)
        IMPLICIT NONE
C-----------------------------------------------------------------------
C========================= Include klspm00.inc =========================
      INTEGER NCLMIN
      PARAMETER( NCLMIN = 4 )
      REAL ECLMIN
      PARAMETER( ECLMIN = 20. )
      REAL vtxdist
      PARAMETER( vtxdist = 50. )
      INTEGER testcounter
      PARAMETER( testcounter = 100 )
      INTEGER MaxNumFitPar
      PARAMETER( MaxNumFitPar = 38 )
      INTEGER MaxNumConstr
      PARAMETER( MaxNumConstr = 10 )
      INTEGER MaxNumComega
      PARAMETER( MaxNumComega = 8 )
      INTEGER cutN
      PARAMETER( cutN = 6 )
      INTEGER errN
      PARAMETER( errN = 17 )
      CHARACTER*40 messageCut(testcounter),messageErr(testcounter)
      LOGICAL    McFlag,EventOK,ErrFlag,truthmc,truthregg,truthomegaa,tr
     1  uthsemii,truththreee,truthelsee                                 
      INTEGER    nmctruth0,ntruth,nmctruth(7),ntruthregg,ntruthsemii,ntr
     1  uththreee,ntruthomegaa,ntruthelsee                              
      INTEGER    selected,analised
      INTEGER    counter(testcounter),counterBcg(testcounter)
      INTEGER    ErrFlagCount(testcounter),ErrFlagCountBcg(testcounter)
      LOGICAL    cutchar,cutneu,onlytruemc,makecuts
      COMMON / flagcommon / McFlag,EventOK,ErrFlag,truthmc
      COMMON / truthcommon / nmctruth0,ntruth,nmctruth
      COMMON / counterscommon / selected,analised,
     &                          counter,counterBcg
      COMMON / errcommon / ErrFlagCount,ErrFlagCountBcg,
     &                     messageCut,messageErr
      COMMON / talk / cutchar,cutneu,onlytruemc,makecuts
      SAVE /flagcommon/, /truthcommon/, /counterscommon/
      SAVE /errcommon/, /talk/
C======================== Include constans.inc =========================
      REAL       Mpip
      PARAMETER( Mpip = 139.57 )             
      REAL       Mpio
      PARAMETER( Mpio = 134.98 )             
      REAL       Mko
      PARAMETER( Mko  = 497.61 )             
      REAL       Momega
      PARAMETER( Momega  = 782.65 )          
      REAL       Cvel
      PARAMETER( Cvel = 29.9792458 )         
      REAL       EMCvel
      PARAMETER( EMCvel = 28.17 )            
      REAL       TauKs
      PARAMETER( TauKs = 0.08953 )   
      REAL       Mphi
      PARAMETER( Mphi  = 1019.460 )          
C======================= Include resolutions.inc =======================
      REAL       sigmaEgg
      PARAMETER( sigmaEgg = 19.348 )           
C= Include /kloe/soft/off/offline/inc/development/tls/maxstructdim.cin =
      INTEGER     NeleCluMax
      PARAMETER ( NeleCluMax = 2000 )
      INTEGER     MaxNumClu
      PARAMETER ( MaxNumClu = 100 )
      INTEGER    MaxNumVtx
      PARAMETER (MaxNumVtx = 20)
      INTEGER    MaxNumTrkV
      PARAMETER (MaxNumTrkV = 30)
      INTEGER    MaxNumTrk
      PARAMETER (MaxNumTrk = 100)
      INTEGER    MaxNumDHSP
      PARAMETER (MaxNumDHSP = 1000)
      Integer    nMaxDC
      Parameter (nMaxDC=1500)
      INTEGER    NQihiMax
      PARAMETER (NQihiMax=1000)
      INTEGER    NQcalMax
      PARAMETER (NQcalMax= 32 )
      INTEGER    MaxNumFirstHit
      PARAMETER (MaxNumFirstHit = 300 )
      INTEGER    MaxNtrkGen
      PARAMETER (MaxNtrkGen =50)
      INTEGER    MaxNvtxGen
      PARAMETER (MaxNvTxGen =30)
      integer TriggerElements
      parameter (TriggerElements = 300 )
C-----------------------------------------------------------------------
C
C External functions
C
      REAL    snrm2,sdot
C
C Local declarations
C
      LOGICAL warning
      REAL    KneRecLor(9)
      REAL    P4PhiLab(4),P4LabInCM(4),P4KneuCM(4),P4KneuLab(4)
      REAL    P4KchCM(4),P4KchLab(4),beta_k(3),ia(3)
      REAL    neuvtxx(MaxNumCLu,3),delta,a,b,c,kpath,kpath_tmp
      REAL    costhe,beta_k_mod,ia_mod,NeuVtx(3)
      REAL    Bpx,Bpy,Bpz,Broots,KchBoost(9),EneCl(MaxNumClu)
      REAL    Xcl(MaxNumClu),Ycl(MaxNumClu),Zcl(MaxNumClu),ip(3)
      REAL    Tcl(MaxNumClu),cldist,trc_tmp
      REAL    trc(MaxNumCLu),trcvtmp(MaxNumClu)
      REAL    KneRec(9),minv4gam,minv4gam_tmp,pi0(2),pi0tmp(2)
      REAL    numerator(3),vcoor(3),x2(3),vdist,denominator(3)
      REAL    vcoorsigm(3),vene,p(3),r(3),test(10),chigamtmp
      REAL    g4vtxerr(3),chigam_tmp(7),chigam,trcv(MaxNumClu)
      REAL    PgamRec(4,4),minv1e(3),minv2e(3),Ltot,PgamRecTaken(4,4)
      INTEGER g4taken(4)
      INTEGER g4ind(4),ii,jj,kk,ll,gpair(3,2),gpairok(2),gpairtaken(2)
      INTEGER i,j,k,ILoop,ncl,ncll(MaxNumCLu)
      INTEGER nclwrong,ncllwrong(MaxNumCLu)
C------------------------------------------------------------------------------
C 1 px
C 2 py
C 3 pz
C 4 E
C 5 |p|
C 6 minv
C 7 xv
C 8 yv
C 9 zv
C in:  Bpx,Bpy,Bpz,Broots,KchBoost(),ncl,enecl(),ncll(),xcl(),ycl(),zcl(),ip(),t
C      cl()                                                                     
C out: cldistKneRecLor(),trc(),nclwrong,ncllwrong
C out: KneRec(),minv4gam,pi0(),g4taken(),trcv(),PgamRecTaken(),gpairtaken,test
C out: g4vtxerr
      P4PhiLab(1) = Bpx
      P4PhiLab(2) = Bpy
      P4PhiLab(3) = Bpz
      P4PhiLab(4) = Broots
      DO i = 1,4
         P4KchLab(i)=KchBoost(i)
      ENDDO
      CALL Loren4(P4PhiLab,P4KchLab,P4KchCM)
      IF ( ErrFlag ) THEN
         RETURN
      ENDIF
      DO i = 1,3
         P4KneuCM(i) = -P4KchCM(i)
         P4LabInCM(i) = -P4PhiLab(i)
      ENDDO
      P4KneuCM(4) = P4KchCM(4)
      P4LabInCM(4) = P4PhiLab(4)
      CALL Loren4(P4LabInCM,P4KneuCM,P4KneuLab)
      IF ( ErrFlag ) THEN
         RETURN
      ENDIF
      DO i = 1,3
         KneRecLor(i) = P4KneuLab(i)
      ENDDO
      KneRecLor(5) = KneRecLor(1)**2 + KneRecLor(2)**2 +
     &               KneRecLor(3)**2
      KneRecLor(4) = SQRT( KneRecLor(5) + Mko**2 )
      IF( ( (KneRecLor(4)**2 - KneRecLor(5)).lt.0. ).or.
     &    ( P4KneuLab(4).eq.0. ) ) THEN
         ErrFlag = .TRUE.
         RETURN
      ENDIF
      KneRecLor(6) = SQRT( KneRecLor(4)**2 - KneRecLor(5) )
      KneRecLor(5) = SQRT( KneRecLor(5) )
      DO i = 1,3
         beta_k(i) = P4KneuLab(i)/P4KneuLab(4)
         KneRecLor(i+6) = 0.
      ENDDO
      beta_k_mod = beta_k(1)**2+beta_k(2)**2+beta_k(3)**2
      beta_k_mod = sqrt(beta_k_mod)
      IF ( ( beta_k_mod.eq.0. ).or.
     &     ((beta_k_mod*Cvel).eq.0).or.((beta_k_mod**2).eq.0) ) THEN
         ErrFlag = .TRUE.
         RETURN
      ENDIF
      nclwrong = 0
      DO i = 1,ncl
        IF( (enecl(ncll(i)).lt.ECLMIN) ) THEN
           nclwrong = nclwrong + 1
           ncllwrong(nclwrong) = ncll(i)
        ELSE
           warning = .FALSE.
           ia(1) = xcl(ncll(i))-ip(1)
           ia(2) = ycl(ncll(i))-ip(2)
           ia(3) = zcl(ncll(i))-ip(3)
           ia_mod = sqrt(ia(1)**2+ia(2)**2+ia(3)**2)
           costhe = 0.
           DO j = 1,3
              costhe = costhe + ia(j)*beta_k(j)
           ENDDO
           IF ( (ia_mod*beta_k_mod).eq.0. ) THEN
              ErrFlag = .TRUE.
              RETURN
           ENDIF
           costhe = costhe/(ia_mod*beta_k_mod)
           a = ( 1. - 1 / (beta_k_mod**2) )
           b = 2*( Cvel*tcl(ncll(i))/beta_k_mod -
     &              ia_mod*costhe )
           c = ia_mod**2 - (Cvel*tcl(ncll(i)))**2
           delta = b**2 - 4*a*c
           IF( (delta.lt.0.).or.(a.eq.0.) ) THEN
              warning = .TRUE.
           ELSEIF( delta.eq.0 ) THEN
              kpath = b/(-2.*a)
              kpath_tmp = kpath
           ELSE
              kpath     = (b-sqrt(delta))/(-2.*a)
              kpath_tmp = (b+sqrt(delta))/(-2.*a)
           ENDIF
           IF( warning ) THEN
              nclwrong = nclwrong + 1
              ncllwrong(nclwrong) = ncll(i)
           ELSE
              DO j = 1,3
                 neuvtxx(i,j) = ip(j) + kpath*beta_k(j)/beta_k_mod
              ENDDO
              trc(i) =  tcl(ncll(i))-sqrt(
     &                 (xcl(ncll(i))-neuvtxx(i,1))**2 +
     &                 (ycl(ncll(i))-neuvtxx(i,2))**2 +
     &                 (zcl(ncll(i))-neuvtxx(i,3))**2  )/Cvel -
     &                 kpath/(beta_k_mod*Cvel)
              DO j = 1,3
                 neuvtxx(i,j) = ip(j) + kpath_tmp*beta_k(j)/beta_k_mod
              ENDDO
              IF( abs(trc(i)).gt.abs(tcl(ncll(i))-sqrt(
     &                 (xcl(ncll(i))-neuvtxx(i,1))**2 +
     &                 (ycl(ncll(i))-neuvtxx(i,2))**2 +
     &                 (zcl(ncll(i))-neuvtxx(i,3))**2  )/Cvel -
     &                 kpath_tmp/(beta_k_mod*Cvel)) ) THEN
                 kpath = kpath_tmp
                 trc(i) = tcl(ncll(i))-sqrt(
     &                    (xcl(ncll(i))-neuvtxx(i,1))**2 +
     &                    (ycl(ncll(i))-neuvtxx(i,2))**2 +
     &                    (zcl(ncll(i))-neuvtxx(i,3))**2  )/Cvel -
     &                    kpath/(beta_k_mod*Cvel)
              ELSE
                DO j = 1,3
                   neuvtxx(i,j) = ip(j) + kpath*beta_k(j)/beta_k_mod
                ENDDO
              ENDIF
              trcvtmp(i)= tcl(ncll(i)) - kpath/(beta_k_mod*Cvel)
           ENDIF
        ENDIF
      ENDDO
      IF( (ncl-nclwrong).lt.NCLMIN ) THEN
          ErrFlag = .TRUE.
          RETURN
      ENDIF
      chigam = 9999999999999.
      DO ii =          1,ncl-3
        DO jj =     ii+1,ncl-2
          DO kk =   jj+1,ncl-1
            DO ll = kk+1,ncl
              DO i = 1,7
                 chigam_tmp(i) = 0.
              ENDDO
              g4ind(1) = ii
              g4ind(2) = jj
              g4ind(3) = kk
              g4ind(4) = ll
              DO k = 1,nclwrong
                 IF( ( ncllwrong(k).eq.ncll(g4ind(1)) ).or.
     &               ( ncllwrong(k).eq.ncll(g4ind(2)) ).or.
     &               ( ncllwrong(k).eq.ncll(g4ind(3)) ).or.
     &               ( ncllwrong(k).eq.ncll(g4ind(4)) ) ) THEN
                        goto 16
                 ENDIF
              ENDDO
              vene = enecl(ncll(g4ind(1))) + enecl(ncll(g4ind(2))) +
     &               enecl(ncll(g4ind(3))) + enecl(ncll(g4ind(4)))
              IF( vene.le.0. ) THEN
                      goto 16
              ENDIF
              DO i = 1,3
                vcoor(i) = ( neuvtxx(g4ind(1),i)*enecl(ncll(g4ind(1))) +
     &                       neuvtxx(g4ind(2),i)*enecl(ncll(g4ind(2))) +
     &                       neuvtxx(g4ind(3),i)*enecl(ncll(g4ind(3))) +
     &                       neuvtxx(g4ind(4),i)*enecl(ncll(g4ind(4))) )
     & / vene
                vcoorsigm(i) = 0.
              ENDDO
              DO i = 1,3
C                vcoorsigm(i) = sqrt(
C     &           ((enecl(ncll(g4ind(1))))*
C     &              (neuvtxx(g4ind(1),i)-vcoor(i)))**2 +
C     &           ((enecl(ncll(g4ind(2))))*
C     &              (neuvtxx(g4ind(2),i)-vcoor(i)))**2 +
C     &           ((enecl(ncll(g4ind(3))))*
C     &              (neuvtxx(g4ind(3),i)-vcoor(i)))**2 +
C     &           ((enecl(ncll(g4ind(4))))*
C     &              (neuvtxx(g4ind(4),i)-vcoor(i)))**2 ) /
C     &           ( vene * 4. * sqrt( 4.-1. ) )
                DO j = 1,4
                   vcoorsigm(i) = vcoorsigm(i) + enecl(ncll(g4ind(j))) *
     &                ( neuvtxx(g4ind(j),i)-vcoor(i) )**2
                ENDDO
                vcoorsigm(i) = sqrt( vcoorsigm(i)/vene )
                IF( vcoorsigm(i).eq.0. ) THEN
                      goto 16
                ENDIF
                x2(i) = ip(i) + KneRecLor(i)
                p(i) = KneRecLor(i)
              ENDDO
              CALL crossProd(vcoor-ip,vcoor-x2,numerator)
              denominator = x2 - ip
C distance between point vcoor and line (defined by point ip and x2)
              IF( snrm2(3,denominator,1).eq.0. ) goto 16
              vdist = snrm2(3,numerator,1)/snrm2(3,denominator,1)
C point at line (defined by point ip and vector p) closest to the point vcoor
              IF( snrm2(3,p,1).eq.0. ) goto 16
              r=(sdot(3,p,1,vcoor-ip,1)/(snrm2(3,p,1))**2)*p
              DO i = 1,3
                 chigam_tmp(1) = chigam_tmp(1) +
     &             ((r(i)-vcoor(i))**2)/(vcoorsigm(i)**2)
              ENDDO
C
C 4-momenta for K --> pi0pi0
C
              DO i = 1,4
C Clust-NeuVtx dist.
                Ltot = SQRT( (Xcl(ncll(g4ind(i)))-vcoor(1))**2 + 
     &                       (Ycl(ncll(g4ind(i)))-vcoor(2))**2 +
     &                       (Zcl(ncll(g4ind(i)))-vcoor(3))**2 )
                IF( Ltot.eq.0. ) THEN
                      goto 16
                ENDIF
                PgamRec(i,1) = EneCl(ncll(g4ind(i))) *
     &                    (Xcl(ncll(g4ind(i)))-vcoor(1)) / Ltot
                PgamRec(i,2) = EneCl(ncll(g4ind(i))) *
     &                    (Ycl(ncll(g4ind(i)))-vcoor(2)) / Ltot
                PgamRec(i,3) = EneCl(ncll(g4ind(i))) *
     &                    (Zcl(ncll(g4ind(i)))-vcoor(3)) / Ltot
                PgamRec(i,4) = EneCl(ncll(g4ind(i)))
              ENDDO
              Minv1e(1) = SQRT(
     &                   ( PgamRec(1,4) + PgamRec(2,4) )**2 -
     &                   ( PgamRec(1,1) + PgamRec(2,1) )**2 -
     &                   ( PgamRec(1,2) + PgamRec(2,2) )**2 -
     &                   ( PgamRec(1,3) + PgamRec(2,3) )**2 )
              Minv2e(1) = SQRT(
     &                   ( PgamRec(3,4) + PgamRec(4,4) )**2 -
     &                   ( PgamRec(3,1) + PgamRec(4,1) )**2 -
     &                   ( PgamRec(3,2) + PgamRec(4,2) )**2 -
     &                   ( PgamRec(3,3) + PgamRec(4,3) )**2 )
              gpair(1,1) = g4ind(1)
              gpair(1,2) = g4ind(2)
              Minv1e(2) = SQRT(
     &                   ( PgamRec(1,4) + PgamRec(3,4) )**2 -
     &                   ( PgamRec(1,1) + PgamRec(3,1) )**2 -
     &                   ( PgamRec(1,2) + PgamRec(3,2) )**2 -
     &                   ( PgamRec(1,3) + PgamRec(3,3) )**2 )
              Minv2e(2) = SQRT(
     &                   ( PgamRec(2,4) + PgamRec(4,4) )**2 -
     &                   ( PgamRec(2,1) + PgamRec(4,1) )**2 -
     &                   ( PgamRec(2,2) + PgamRec(4,2) )**2 -
     &                   ( PgamRec(2,3) + PgamRec(4,3) )**2 )
              gpair(2,1) = g4ind(1)
              gpair(2,2) = g4ind(3)
              Minv1e(3) = SQRT(
     &                   ( PgamRec(1,4) + PgamRec(4,4) )**2 -
     &                   ( PgamRec(1,1) + PgamRec(4,1) )**2 -
     &                   ( PgamRec(1,2) + PgamRec(4,2) )**2 -
     &                   ( PgamRec(1,3) + PgamRec(4,3) )**2 )
              Minv2e(3) = SQRT(
     &                   ( PgamRec(2,4) + PgamRec(3,4) )**2 -
     &                   ( PgamRec(2,1) + PgamRec(3,1) )**2 -
     &                   ( PgamRec(2,2) + PgamRec(3,2) )**2 -
     &                   ( PgamRec(2,3) + PgamRec(3,3) )**2 )
              gpair(3,1) = g4ind(1)
              gpair(3,2) = g4ind(4)
              IF( ( (ABS(Minv1e(1)-Mpio)+ABS(Minv2e(1)-Mpio)).le.
     &              (ABS(Minv1e(2)-Mpio)+ABS(Minv2e(2)-Mpio))).and.
     &            ( (ABS(Minv1e(1)-Mpio)+ABS(Minv2e(1)-Mpio)).le.
     &              (ABS(Minv1e(3)-Mpio)+ABS(Minv2e(3)-Mpio))))THEN
                pi0tmp(1) = Minv1e(1)
                pi0tmp(2) = Minv2e(1)
                gpairok(1) = gpair(1,1)
                gpairok(2) = gpair(1,2)
              ELSEIF( ( (ABS(Minv1e(2)-Mpio)+ABS(Minv2e(2)-Mpio)).le.
     &                  (ABS(Minv1e(1)-Mpio)+ABS(Minv2e(1)-Mpio))).and.
     &                ( (ABS(Minv1e(2)-Mpio)+ABS(Minv2e(2)-Mpio)).le.
     &                  (ABS(Minv1e(3)-Mpio)+ABS(Minv2e(3)-Mpio))))THEN
                pi0tmp(1) = Minv1e(2)
                pi0tmp(2) = Minv2e(2)
                gpairok(1) = gpair(2,1)
                gpairok(2) = gpair(2,2)
              ELSEIF( ( (ABS(Minv1e(3)-Mpio)+ABS(Minv2e(3)-Mpio)).le.
     &                  (ABS(Minv1e(1)-Mpio)+ABS(Minv2e(1)-Mpio))).and.
     &                ( (ABS(Minv1e(3)-Mpio)+ABS(Minv2e(3)-Mpio)).le.
     &                  (ABS(Minv1e(2)-Mpio)+ABS(Minv2e(2)-Mpio))))THEN
                pi0tmp(1) = Minv1e(3)
                pi0tmp(2) = Minv2e(3)
                gpairok(1) = gpair(3,1)
                gpairok(2) = gpair(3,2)
              ELSE
                goto 16
              ENDIF
              minv4gam_tmp = sqrt(
     &              ( PgamRec(1,4) + PgamRec(2,4) +
     &                PgamRec(3,4) + PgamRec(4,4) )**2 -
     &              ( PgamRec(1,1) + PgamRec(2,1) +
     &                PgamRec(3,1) + PgamRec(4,1) )**2 -
     &              ( PgamRec(1,2) + PgamRec(2,2) +
     &                PgamRec(3,2) + PgamRec(4,2) )**2 -
     &              ( PgamRec(1,3) + PgamRec(2,3) +
     &                PgamRec(3,3) + PgamRec(4,3) )**2 )
             IF( .not.(minv4gam_tmp.gt.0.and.
     &                    pi0tmp(1).gt.0.and.
     &                    pi0tmp(2).gt.0.) ) THEN
                   goto 16
              ENDIF
              chigam_tmp(2) = chigam_tmp(2) +
     &            ( (minv4gam_tmp-Mko) / (2*19.348) )**2
              DO i = 1,2
                chigam_tmp(3) = chigam_tmp(3) +
     &           ( ( (pi0tmp(i)-Mpio) / (19.348) )**2 )/2.
              ENDDO
              chigam_tmp(4) = sqrt( vcoorsigm(1)**2 +
     &                              vcoorsigm(2)**2 +
     &                              vcoorsigm(3)**2 )
              DO i = 1,4
                 trcvtmp(g4ind(i)) = trcvtmp(g4ind(i)) - sqrt(
     &             (xcl(ncll(g4ind(i)))-vcoor(1))**2 +
     &             (ycl(ncll(g4ind(i)))-vcoor(2))**2 +
     &             (zcl(ncll(g4ind(i)))-vcoor(3))**2 )/Cvel
                 chigam_tmp(5) = chigam_tmp(5) +
     &           enecl(ncll(g4ind(i)))*abs(trcvtmp(g4ind(i)))
                 chigam_tmp(6) = chigam_tmp(6) + sqrt(
     &             (vcoor(1)-neuvtxx(g4ind(i),1))**2 +
     &             (vcoor(2)-neuvtxx(g4ind(i),2))**2 +
     &             (vcoor(3)-neuvtxx(g4ind(i),3))**2 ) *
     &                 enecl(ncll(g4ind(i)))
                 DO j = 1,3
                    chigam_tmp(7) = chigam_tmp(7) +
     &                 (((vcoor(j)-neuvtxx(g4ind(i),j))/
     &                 (vcoorsigm(j)))**2)*enecl(ncll(g4ind(i)))
                 ENDDO
              ENDDO
              DO i = 5,7
                chigam_tmp(i) = chigam_tmp(i)/vene
              ENDDO
              chigamtmp =   0.
              chigamtmp =   chigamtmp
C    &                    + chigam_tmp(1)
C    &                    + chigam_tmp(2)
C    &                    + chigam_tmp(3)
C    &                    + chigam_tmp(4)
C    &                    + chigam_tmp(5)
     &                    + chigam_tmp(6)
C    &                    + chigam_tmp(7)
              IF( chigamtmp.lt.chigam ) THEN
                 chigam = chigamtmp
C test
                 DO i = 1,7
                    test(i) = chigam_tmp(i)
                 ENDDO
                 test(8)  = r(1)
                 test(9)  = r(2)
                 test(10) = r(3)
                 KneRec(1) = PgamRec(1,1) + PgamRec(2,1) +
     &                       PgamRec(3,1) + PgamRec(4,1)
                 KneRec(2) = PgamRec(1,2) + PgamRec(2,2) +
     &                       PgamRec(3,2) + PgamRec(4,2)
                 KneRec(3) = PgamRec(1,3) + PgamRec(2,3) +
     &                       PgamRec(3,3) + PgamRec(4,3)
                 KneRec(4) = PgamRec(1,4) + PgamRec(2,4) +
     &                       PgamRec(3,4) + PgamRec(4,4)
                 KneRec(5) = KneRec(1)**2+KneRec(2)**2+KneRec(3)**2
                 IF( (KneRec(4)**2 - KneRec(5) ).lt.0. ) THEN
                        WRITE(*,*)' ERROR ERYK 1'
                        ErrFlag = .TRUE.
                        RETURN
                 ENDIF
                 KneRec(6) = SQRT( KneRec(4)**2 - KneRec(5) )
                 KneRec(5) = sqrt( KneRec(5) )
                 DO i = 1,3
                   KneRecLor(i+6) = vcoor(i)
                   KneRec(i+6)    = vcoor(i)
                   g4vtxerr(i) = vcoorsigm(i)
                 ENDDO
                 DO i = 1,4
                   g4taken(i) = g4ind(i)
                   trcv(g4taken(i)) = trcvtmp(g4taken(i))
                 ENDDO
                 cldist = vdist
                 pi0(1) = pi0tmp(1)
                 pi0(2) = pi0tmp(2)
                 minv4gam = minv4gam_tmp
                 gpairtaken = gpairok
                 DO i = 1,4
                   DO j = 1,4
                     PgamRecTaken(i,j) = PgamRec(i,j)
                   ENDDO
                 ENDDO
              ENDIF
 16           continue
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      END
C=======================================================================
C=======================================================================
C=======================================================================
      SUBROUTINE find_omegaminv(FitPar,P4gam,Kch,
     &                          MinvOmega,PpioOmega,P4PriRest,coss)
        IMPLICIT NONE
C-----------------------------------------------------------------------
C======================== Include constans.inc =========================
      REAL       Mpip
      PARAMETER( Mpip = 139.57 )             
      REAL       Mpio
      PARAMETER( Mpio = 134.98 )             
      REAL       Mko
      PARAMETER( Mko  = 497.61 )             
      REAL       Momega
      PARAMETER( Momega  = 782.65 )          
      REAL       Cvel
      PARAMETER( Cvel = 29.9792458 )         
      REAL       EMCvel
      PARAMETER( EMCvel = 28.17 )            
      REAL       TauKs
      PARAMETER( TauKs = 0.08953 )   
      REAL       Mphi
      PARAMETER( Mphi  = 1019.460 )          
C========================= Include klspm00.inc =========================
      INTEGER NCLMIN
      PARAMETER( NCLMIN = 4 )
      REAL ECLMIN
      PARAMETER( ECLMIN = 20. )
      REAL vtxdist
      PARAMETER( vtxdist = 50. )
      INTEGER testcounter
      PARAMETER( testcounter = 100 )
      INTEGER MaxNumFitPar
      PARAMETER( MaxNumFitPar = 38 )
      INTEGER MaxNumConstr
      PARAMETER( MaxNumConstr = 10 )
      INTEGER MaxNumComega
      PARAMETER( MaxNumComega = 8 )
      INTEGER cutN
      PARAMETER( cutN = 6 )
      INTEGER errN
      PARAMETER( errN = 17 )
      CHARACTER*40 messageCut(testcounter),messageErr(testcounter)
      LOGICAL    McFlag,EventOK,ErrFlag,truthmc,truthregg,truthomegaa,tr
     1  uthsemii,truththreee,truthelsee                                 
      INTEGER    nmctruth0,ntruth,nmctruth(7),ntruthregg,ntruthsemii,ntr
     1  uththreee,ntruthomegaa,ntruthelsee                              
      INTEGER    selected,analised
      INTEGER    counter(testcounter),counterBcg(testcounter)
      INTEGER    ErrFlagCount(testcounter),ErrFlagCountBcg(testcounter)
      LOGICAL    cutchar,cutneu,onlytruemc,makecuts
      COMMON / flagcommon / McFlag,EventOK,ErrFlag,truthmc
      COMMON / truthcommon / nmctruth0,ntruth,nmctruth
      COMMON / counterscommon / selected,analised,
     &                          counter,counterBcg
      COMMON / errcommon / ErrFlagCount,ErrFlagCountBcg,
     &                     messageCut,messageErr
      COMMON / talk / cutchar,cutneu,onlytruemc,makecuts
      SAVE /flagcommon/, /truthcommon/, /counterscommon/
      SAVE /errcommon/, /talk/
C= Include /kloe/soft/off/offline/inc/development/tls/maxstructdim.cin =
      INTEGER     NeleCluMax
      PARAMETER ( NeleCluMax = 2000 )
      INTEGER     MaxNumClu
      PARAMETER ( MaxNumClu = 100 )
      INTEGER    MaxNumVtx
      PARAMETER (MaxNumVtx = 20)
      INTEGER    MaxNumTrkV
      PARAMETER (MaxNumTrkV = 30)
      INTEGER    MaxNumTrk
      PARAMETER (MaxNumTrk = 100)
      INTEGER    MaxNumDHSP
      PARAMETER (MaxNumDHSP = 1000)
      Integer    nMaxDC
      Parameter (nMaxDC=1500)
      INTEGER    NQihiMax
      PARAMETER (NQihiMax=1000)
      INTEGER    NQcalMax
      PARAMETER (NQcalMax= 32 )
      INTEGER    MaxNumFirstHit
      PARAMETER (MaxNumFirstHit = 300 )
      INTEGER    MaxNtrkGen
      PARAMETER (MaxNtrkGen =50)
      INTEGER    MaxNvtxGen
      PARAMETER (MaxNvTxGen =30)
      integer TriggerElements
      parameter (TriggerElements = 300 )
C-----------------------------------------------------------------------
      REAL    MinvOmega(2),P4gam(4,4),Kch(9),P4Sys(4),Ptr1Omega(4)
      REAL    Ptr2Omega(4),PpioOmega(4),P4PriRest(4),P4Kne(4)
      REAL    Chi2Min,BestDiff,Chi2Pio,P4Pri(4),P4pio(4)
      REAL    mpio1,mpio2,minv1,minv2,PtrkTmp1(4),PtrkTmp2(4)
      REAl    FitPar(6),P4OmegaRest(4),p1(9),p2(9),coss,P4Omega(4)
      INTEGER ILoop,JLoop,LLoop,KLoop,k,MLoop
C in:  FitPar(),P4gam(),Kch()
C out: MinvOmega,PpioOmega(),P4PriRest(),coss
C-----------------------------------------------------------------------
C
C Best pi0 and omega masses
C
      Chi2Min = 999999.
      BestDiff = 999999.
C
      DO ILoop = 1,1
         DO JLoop = ILoop+1,4
            KLoop = 0
            LLoop = 0
            DO MLoop = 2,4
               IF( MLoop.NE.ILoop .AND. MLoop.NE.JLoop )THEN
                  IF( KLoop.EQ.0 )THEN
                     KLoop = MLoop
                  ELSE
                     LLoop = MLoop
                  ENDIF
               ENDIF
            ENDDO
            Minv1 = ( P4gam(ILoop,4) + P4gam(JLoop,4) )**2 -
     &           ( P4gam(ILoop,1) + P4gam(JLoop,1) )**2 -
     &           ( P4gam(ILoop,2) + P4gam(JLoop,2) )**2 -
     &           ( P4gam(ILoop,3) + P4gam(JLoop,3) )**2
            IF( Minv1.gt.0. ) THEN
                    Minv1 = SQRT(Minv1)
            ELSE
                    Minv1 = 0.
            ENDIF
            Minv2 = ( P4gam(KLoop,4) + P4gam(LLoop,4) )**2 -
     &           ( P4gam(KLoop,1) + P4gam(LLoop,1) )**2 -
     &           ( P4gam(KLoop,2) + P4gam(LLoop,2) )**2 -
     &           ( P4gam(KLoop,3) + P4gam(LLoop,3) )**2
            IF( Minv2.gt.0. ) THEN
                    Minv2 = SQRT(Minv2)
            ELSE
                    Minv2 = 0.
            ENDIF
            Chi2Pio = ((Minv1-Mpio)/4.)**2 + ((Minv2-Mpio)/4.)**2
            IF( Chi2Pio.LT.Chi2Min )THEN
               Chi2Min = Chi2Pio
               Mpio1 = Minv1
               Mpio2 = Minv2
CCC Omega mass...
               Minv1 =
     &              (P4gam(ILoop,4)+P4gam(JLoop,4)+Kch(4))**2 -
     &              (P4gam(ILoop,1)+P4gam(JLoop,1)+Kch(1))**2 -
     &              (P4gam(ILoop,2)+P4gam(JLoop,2)+Kch(2))**2 -
     &              (P4gam(ILoop,3)+P4gam(JLoop,3)+Kch(3))**2
               IF( Minv1.gt.0. ) THEN
                       Minv1 = SQRT(Minv1)
               ELSE
                       Minv1 = 0.
               ENDIF
               Minv2 =
     &              (P4gam(KLoop,4)+P4gam(LLoop,4)+Kch(4))**2 -
     &              (P4gam(KLoop,1)+P4gam(LLoop,1)+Kch(1))**2 -
     &              (P4gam(KLoop,2)+P4gam(LLoop,2)+Kch(2))**2 -
     &              (P4gam(KLoop,3)+P4gam(LLoop,3)+Kch(3))**2
               IF( Minv2.gt.0. ) THEN
                       Minv2 = SQRT(Minv2)
               ELSE
                       Minv2 = 0.
               ENDIF
               IF( ABS(Minv1-Momega).LT.BestDiff )THEN
                  BestDiff = ABS(Minv1-Momega)
                  MinvOmega(1) = Minv1
                  MinvOmega(2) = Minv2
                  P4pio(1) = P4gam(ILoop,1) + P4gam(JLoop,1)
                  P4pio(2) = P4gam(ILoop,2) + P4gam(JLoop,2)
                  P4pio(3) = P4gam(ILoop,3) + P4gam(JLoop,3)
                  P4pio(4) = P4gam(ILoop,4) + P4gam(JLoop,4)
C Primary pion
                  P4Pri(1) = P4gam(KLoop,1) + P4gam(LLoop,1)
                  P4Pri(2) = P4gam(KLoop,2) + P4gam(LLoop,2)
                  P4Pri(3) = P4gam(KLoop,3) + P4gam(LLoop,3)
                  P4Pri(4) = P4gam(KLoop,4) + P4gam(LLoop,4)
               ENDIF
               IF( ABS(Minv2-Momega).LT.BestDiff )THEN
                  BestDiff = ABS(Minv2-Momega)
                  MinvOmega(1) = Minv2
                  MinvOmega(2) = Minv1
                  P4pio(1) = P4gam(KLoop,1) + P4gam(LLoop,1)
                  P4pio(2) = P4gam(KLoop,2) + P4gam(LLoop,2)
                  P4pio(3) = P4gam(KLoop,3) + P4gam(LLoop,3)
                  P4pio(4) = P4gam(KLoop,4) + P4gam(LLoop,4)
C
                  P4Pri(1) = P4gam(ILoop,1) + P4gam(JLoop,1)
                  P4Pri(2) = P4gam(ILoop,2) + P4gam(JLoop,2)
                  P4Pri(3) = P4gam(ILoop,3) + P4gam(JLoop,3)
                  P4Pri(4) = P4gam(ILoop,4) + P4gam(JLoop,4)
               ENDIF
CCC
            ENDIF
         ENDDO
      ENDDO
C
C------------------------------------------------------------------------------
C Find pion energies in omega rest frame
C------------------------------------------------------------------------------
C
      DO ILoop = 1,4
         P4Sys(ILoop) = Kch(ILoop) + P4pio(ILoop)
      ENDDO
C
      PtrkTmp1(1) = FitPar(1)
      PtrkTmp1(2) = FitPar(2)
      PtrkTmp1(3) = FitPar(3)
      PtrkTmp1(4) = SQRT( PtrkTmp1(1)**2 + PtrkTmp1(2)**2 +
     &                    PtrkTmp1(3)**2 + Mpip**2 )
      PtrkTmp2(1) = FitPar(4)
      PtrkTmp2(2) = FitPar(5)
      PtrkTmp2(3) = FitPar(6)
      PtrkTmp1(4) = SQRT( PtrkTmp2(1)**2 + PtrkTmp2(2)**2 +
     &                    PtrkTmp2(3)**2 + Mpip**2 )
      CALL Loren4(P4Sys,PtrkTmp1,Ptr1Omega)
      IF ( ErrFlag ) RETURN
      CALL Loren4(P4Sys,PtrkTmp2,Ptr2Omega)
      IF ( ErrFlag ) RETURN
      CALL Loren4(P4Sys,P4pio,PpioOmega)
      IF ( ErrFlag ) RETURN
C
C------------------------------------------------------------------------------
C Find pion energies in phi rest frame
C------------------------------------------------------------------------------
C
      P4Kne(1) = P4gam(1,1) + P4gam(2,1) + P4gam(3,1) + P4gam(4,1)
      P4Kne(2) = P4gam(1,2) + P4gam(2,2) + P4gam(3,2) + P4gam(4,2)
      P4Kne(3) = P4gam(1,3) + P4gam(2,3) + P4gam(3,3) + P4gam(4,3)
      P4Kne(4) = P4gam(1,4) + P4gam(2,4) + P4gam(3,4) + P4gam(4,4)
      DO ILoop = 1,4
         P4Sys(ILoop) = Kch(ILoop) + P4Kne(ILoop)
      ENDDO
      CALL Loren4(P4Sys,P4Pri,P4PriRest)
      IF ( ErrFlag ) RETURN
      P4Omega = P4Sys-P4Pri
      CALL Loren4(P4Sys,P4Omega,P4OmegaRest)
      IF ( ErrFlag ) RETURN
      DO k=1,4
        p1(k) = P4PriRest(k)
        p2(k) = P4OmegaRest(k)
      ENDDO
      DO k=5,9
        p1(k) = 0.
        p2(k) = 0.
      ENDDO
      CALL do_cosin(p1,p2,coss)
      END
C=======================================================================
C=======================================================================
C=======================================================================
      SUBROUTINE do_cuts
        IMPLICIT NONE
C========================= Include klspm00.inc =========================
      INTEGER NCLMIN
      PARAMETER( NCLMIN = 4 )
      REAL ECLMIN
      PARAMETER( ECLMIN = 20. )
      REAL vtxdist
      PARAMETER( vtxdist = 50. )
      INTEGER testcounter
      PARAMETER( testcounter = 100 )
      INTEGER MaxNumFitPar
      PARAMETER( MaxNumFitPar = 38 )
      INTEGER MaxNumConstr
      PARAMETER( MaxNumConstr = 10 )
      INTEGER MaxNumComega
      PARAMETER( MaxNumComega = 8 )
      INTEGER cutN
      PARAMETER( cutN = 6 )
      INTEGER errN
      PARAMETER( errN = 17 )
      CHARACTER*40 messageCut(testcounter),messageErr(testcounter)
      LOGICAL    McFlag,EventOK,ErrFlag,truthmc,truthregg,truthomegaa,tr
     1  uthsemii,truththreee,truthelsee                                 
      INTEGER    nmctruth0,ntruth,nmctruth(7),ntruthregg,ntruthsemii,ntr
     1  uththreee,ntruthomegaa,ntruthelsee                              
      INTEGER    selected,analised
      INTEGER    counter(testcounter),counterBcg(testcounter)
      INTEGER    ErrFlagCount(testcounter),ErrFlagCountBcg(testcounter)
      LOGICAL    cutchar,cutneu,onlytruemc,makecuts
      COMMON / flagcommon / McFlag,EventOK,ErrFlag,truthmc
      COMMON / truthcommon / nmctruth0,ntruth,nmctruth
      COMMON / counterscommon / selected,analised,
     &                          counter,counterBcg
      COMMON / errcommon / ErrFlagCount,ErrFlagCountBcg,
     &                     messageCut,messageErr
      COMMON / talk / cutchar,cutneu,onlytruemc,makecuts
      SAVE /flagcommon/, /truthcommon/, /counterscommon/
      SAVE /errcommon/, /talk/
C======================== Include constans.inc =========================
      REAL       Mpip
      PARAMETER( Mpip = 139.57 )             
      REAL       Mpio
      PARAMETER( Mpio = 134.98 )             
      REAL       Mko
      PARAMETER( Mko  = 497.61 )             
      REAL       Momega
      PARAMETER( Momega  = 782.65 )          
      REAL       Cvel
      PARAMETER( Cvel = 29.9792458 )         
      REAL       EMCvel
      PARAMETER( EMCvel = 28.17 )            
      REAL       TauKs
      PARAMETER( TauKs = 0.08953 )   
      REAL       Mphi
      PARAMETER( Mphi  = 1019.460 )          
C-----------------------------------------------------------------------
C= Include /kloe/soft/off/offline/inc/development/tls/maxstructdim.cin =
      INTEGER     NeleCluMax
      PARAMETER ( NeleCluMax = 2000 )
      INTEGER     MaxNumClu
      PARAMETER ( MaxNumClu = 100 )
      INTEGER    MaxNumVtx
      PARAMETER (MaxNumVtx = 20)
      INTEGER    MaxNumTrkV
      PARAMETER (MaxNumTrkV = 30)
      INTEGER    MaxNumTrk
      PARAMETER (MaxNumTrk = 100)
      INTEGER    MaxNumDHSP
      PARAMETER (MaxNumDHSP = 1000)
      Integer    nMaxDC
      Parameter (nMaxDC=1500)
      INTEGER    NQihiMax
      PARAMETER (NQihiMax=1000)
      INTEGER    NQcalMax
      PARAMETER (NQcalMax= 32 )
      INTEGER    MaxNumFirstHit
      PARAMETER (MaxNumFirstHit = 300 )
      INTEGER    MaxNtrkGen
      PARAMETER (MaxNtrkGen =50)
      INTEGER    MaxNvtxGen
      PARAMETER (MaxNvTxGen =30)
      integer TriggerElements
      parameter (TriggerElements = 300 )
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
         INTEGER iv(MaxNumtrkv)
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
         INTEGER mother(MaxNvtxGen)
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
C====================== Include interfstruct.inc =======================
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
          REAL    Tcl(MaxNumClu)
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
          INTEGER iv(MaxNumtrkv)
          INTEGER trknumv(MaxNumtrkv)
          REAL    CurV(MaxNumtrkv)
          REAL    PhiV(MaxNumtrkv)
          REAL    CotV(MaxNumtrkv)
          REAL    PxTV(MaxNumtrkv)
          REAL    PyTV(MaxNumtrkv)
          REAL    PzTV(MaxNumtrkv)
          INTEGER nv
          INTEGER vtx(MaxNumVtx)
          REAL    xv(MaxNumVtx)
          REAL    yv(MaxNumVtx)
          REAL    zv(MaxNumVtx)
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
          INTEGER pidmc(MaxNtrkGen)
          INTEGER virmom(MaxNtrkGen)
          REAL    pxmc(MaxNtrkGen)
          REAL    pymc(MaxNtrkGen)
          REAL    pzmc(MaxNtrkGen)
          REAL    themc(MaxNtrkGen)
          REAL    phimc(MaxNtrkGen)
          INTEGER vtxmc(MaxNtrkGen)
          INTEGER nvtxmc
          INTEGER kinmom(MaxNvtxGen)
          INTEGER mother(MaxNvtxGen)
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
          REAL    KchMC(9)
          REAL    KneMC(9)
          REAL    KchRec(9)
          REAL    KneRec(9)
          REAL    KchBoost(9)
          REAL    ip(3)
          REAL    ipmc(3)
          INTEGER vtaken(3)
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
          REAL    cosTrk
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
          INTEGER omegatruth_top
          INTEGER clu_index(4)
          INTEGER chtrk_index(2)
          INTEGER chvtx_index
          REAL P4gam_w(4,4)
          REAL P4ch_w(2,4)
          REAL P4ch_tot_w(4)
          REAL om_dist(4)
      END TYPE
      TYPE ( interfstru) Interf
      COMMON / InterfCommon / Interf,Evt
C-----------------------------------------------------------------------
C Local variables
      INTEGER i
      REAL    tmp
C-----------------------------------------------------------------------
C
      EventOK = .FALSE.
C      IF( Interf%nt.lt.2 ) THEN
C        CALL statisticss(1)
C        RETURN
C      ENDIF
C      IF( sqrt((Interf%pi0(1)-Mpio)**2+
C     &         (Interf%pi0(2)-Mpio)**2).gt.35 ) THEN
C        CALL statisticss(2)
C        RETURN
C      ENDIF
C      tmp=0.
C      DO i=1,4
C         tmp=tmp+Interf%trcv(Interf%g4taken(i))
C      ENDDO
C      IF( tmp.lt.-1 ) THEN
C        CALL statisticss(3)
C        RETURN
C      ENDIF
       IF( abs(Interf%minv4gam-Mko).gt.150 ) THEN
         CALL statisticss(4)
         RETURN
       ENDIF
       IF( abs(Interf%kchrec(6)-Mko).gt.2 ) THEN
         CALL statisticss(5)
         RETURN
       ENDIF
C      IF( Interf%Qmiss.gt.6 ) THEN
C        CALL statisticss(6)
C        RETURN
C      ENDIF
      EventOK = .TRUE.
CC probably do_distance function is wrong, (Eryk, 25.06.2013)
      END
C=======================================================================
C=======================================================================
C=======================================================================
      SUBROUTINE fit_interf_pmoo(P,DP,M,CHISQR,C,D,MTX,Z,loopcount)
C
C Description:
C ------------
C
C Kinematic Fit for Phi --> KlKs --> pi+pi-pi0pi0 events
C
C Number of input parameters: 38
C * Curvature, cotg(theta), phi for each track   (3x2)
C * E, X, Y, Z, T for each photon                (5x4)
C * Charge, neutral and phi vertex coordinates   (3x3)
C * Beam Energies                                (1x2)
C * Phi momentum along x axis                    ( 1 )
C
C Total Constraints: 10
C * Time of flight of photons (4)
C * 4-momentum conservation   (4)
C * M(pi+pi-) = M(ko)         (1)
C * M( gggg ) = M(ko)         (1)
C
C Input parameters:
C -----------------
C M      = Number of constraints
C N      = Number of input parameters
C P (N)  = Measured Parameters
C DP(N)  = Parameter Errors
C
C Output parameters:
C ------------------
C CHISQR  = Chi2 value
C P (N)   = Adjusted Parameters
C loopcount = number of iterations
C
C-----------------------------------------------------------------------
        IMPLICIT NONE
C========================= Include klspm00.inc =========================
      INTEGER NCLMIN
      PARAMETER( NCLMIN = 4 )
      REAL ECLMIN
      PARAMETER( ECLMIN = 20. )
      REAL vtxdist
      PARAMETER( vtxdist = 50. )
      INTEGER testcounter
      PARAMETER( testcounter = 100 )
      INTEGER MaxNumFitPar
      PARAMETER( MaxNumFitPar = 38 )
      INTEGER MaxNumConstr
      PARAMETER( MaxNumConstr = 10 )
      INTEGER MaxNumComega
      PARAMETER( MaxNumComega = 8 )
      INTEGER cutN
      PARAMETER( cutN = 6 )
      INTEGER errN
      PARAMETER( errN = 17 )
      CHARACTER*40 messageCut(testcounter),messageErr(testcounter)
      LOGICAL    McFlag,EventOK,ErrFlag,truthmc,truthregg,truthomegaa,tr
     1  uthsemii,truththreee,truthelsee                                 
      INTEGER    nmctruth0,ntruth,nmctruth(7),ntruthregg,ntruthsemii,ntr
     1  uththreee,ntruthomegaa,ntruthelsee                              
      INTEGER    selected,analised
      INTEGER    counter(testcounter),counterBcg(testcounter)
      INTEGER    ErrFlagCount(testcounter),ErrFlagCountBcg(testcounter)
      LOGICAL    cutchar,cutneu,onlytruemc,makecuts
      COMMON / flagcommon / McFlag,EventOK,ErrFlag,truthmc
      COMMON / truthcommon / nmctruth0,ntruth,nmctruth
      COMMON / counterscommon / selected,analised,
     &                          counter,counterBcg
      COMMON / errcommon / ErrFlagCount,ErrFlagCountBcg,
     &                     messageCut,messageErr
      COMMON / talk / cutchar,cutneu,onlytruemc,makecuts
      SAVE /flagcommon/, /truthcommon/, /counterscommon/
      SAVE /errcommon/, /talk/
C-----------------------------------------------------------------------
C
C External functions
C
        REAL VDOTN, VMOD
C
C Local declarations
C
        REAL XI2
        COMMON/PHIFIT/XI2(50)
        INTEGER I, J, K, L, N
C Number of parameters
        PARAMETER (N=36)                  
C Number of constraints in fit
        INTEGER M                         
C Parameters & errors
        REAL    P(N), DP(N), P0(N), P1(N) 
C Parameter Corrections
        REAL*8  A(N)                      
C Derivative of constraint M w.r.t. par. N
        REAL*8  D(N,M)                    
C Balance of constraints
        REAL*8  C(M)                      
C Covariance Matrix (?)
        REAL*8  MTX(N+M,N+M)              
C (0,-C)
        REAL*8  Z(N+M)                    
        REAL*8  C1(M)
        INTEGER ITERATIONMAX, IFAIL, loopcount
        REAL    CHISQR, DX2, DX2MAX, R(100), C2
        REAL    PTEMP(N)
        LOGICAL FailFlag
C------------------------------------------------------------------------------
C in:
C out:
        FailFlag = .FALSE.
C       ITERATIONMAX = 50
        ITERATIONMAX = 100
        CHISQR = 0.
        DX2 = 1.0e6
        DX2MAX = 0.01
        DO I = 1,N
          A(I)  = 0.
          P0(I) = P(I)
          Z(I)  = 0.
          DO J = 1,M
            D(I,J) = 0.d0
          ENDDO
        ENDDO
C
        CALL CONSTRAINT(M,N,P,C)
C
        loopcount = 0
        DO WHILE( loopcount.LT.ITERATIONMAX .AND. DX2.GT.DX2MAX )
          loopcount = loopcount+1
          DO I = 1,N
            P1(I) = P(I)
          ENDDO
         CALL DERIVATIVE(N,M,D,C,P,DP,PTEMP,C1)
          DO J = 1,M
            Z(N+J) = -C(J)
          ENDDO
          DO I = 1,N+M
            DO J = I,N+M
              IF( I.GT.N )THEN
                MTX(I,J) = 0.
              ELSEIF( J.GT.N )THEN
                MTX(I,J) = D(I,J-N)
              ELSEIF( I.EQ.J .AND. DP(I).GT.0 )THEN
                MTX(I,J) = 1/(DP(I)*DP(I))
              ELSE
                MTX(I,J) = 0.
              ENDIF
              IF( J.GT.I )THEN
                MTX(J,I) = MTX(I,J)
              ENDIF
            ENDDO
          ENDDO
          CALL DINV(N+M,MTX,N+M,R,IFAIL)
          IF(IFAIL.NE.0) FailFlag = .TRUE.
          DX2 = CHISQR
          CHISQR = 0.
          DO I = 1,N
            A(I) = 0.
            DO J = 1,N+M
              A(I) = A(I) + MTX(I,J)*Z(J)
            ENDDO
C            P(I)= P1(I) + A(I)
            IF( ABS(A(I)).LT.1.0E3 )THEN
              P(I) = P1(I) + A(I)
            ENDIF
CC--- To avoid negative energy photons... ------------------
C!            IF( MOD(I+4,5).EQ.0 )THEN
C!               P(I) = ABS(P(I))
C!            ENDIF
            IF( I.EQ.7.or.I.EQ.12.or.I.EQ.17.or.I.EQ.22 )THEN
              IF( P(I).LT. ECLMIN ) THEN
                P(I) = ECLMIN
              ENDIF
            ENDIF
CC----------------------------------------------------------
            IF( DP(I).GT.0 )THEN
              CHISQR = CHISQR + ((P(I)-P0(I))/DP(I))**2
            ENDIF
          ENDDO
          DX2 = ABS(CHISQR-DX2)
          CALL CONSTRAINT(M,N,P,C)
          C2 = 0.
          DO J = 1,M
            C2 = C2 + C(J)*C(J)
          ENDDO
          XI2(loopcount) = CHISQR
        ENDDO
C
        IF( FailFlag ) loopcount = 100 + loopcount
C
        RETURN
        END
C
C
C=======================================================================
C=======================================================================
C=======================================================================
        SUBROUTINE DERIVATIVE(N,M,D,C,P,DP,P1,C1)
C
        IMPLICIT NONE
C
C Local declarations
C
        INTEGER I,J,K,L,M,N
        REAL    P(N), DP(N), P1(N)
        REAL*8  D(N,M), C(M), C1(M), D1, STEP
C------------------------------------------------------------------------------
        DO I=1,N
          P1(I)=P(I)
        ENDDO
        DO J=1,M
          DO I=1,N
            IF(DP(I).EQ.0)THEN
              IF(P1(I).EQ.0)THEN
                STEP=0.1
              ELSE
                STEP=0.001*P1(I)
              ENDIF
            ELSE
              STEP=0.01*DP(I)
            ENDIF
            L=0
            DO WHILE(ABS(D1-D(I,J)).GT.0.001*ABS(D(I,J)).OR.L.LE.1)
              D1=D(I,J)
              P1(I)=P(I)+STEP
              CALL CONSTRAINT(M,N,P1,C1)
              D(I,J)=(C1(J)-C(J))/STEP
              STEP=STEP/2
              L=L+1
            ENDDO
            P1(I)=P(I)
          ENDDO
        ENDDO
C
        RETURN
        END
C
C==============================================================================
C==============================================================================
C==============================================================================
C==============================================================================
        SUBROUTINE CONSTRAINT(M,N,P,C)
C
        IMPLICIT NONE
C-----------------------------------------------------------------------
C========================= Include klspm00.inc =========================
      INTEGER NCLMIN
      PARAMETER( NCLMIN = 4 )
      REAL ECLMIN
      PARAMETER( ECLMIN = 20. )
      REAL vtxdist
      PARAMETER( vtxdist = 50. )
      INTEGER testcounter
      PARAMETER( testcounter = 100 )
      INTEGER MaxNumFitPar
      PARAMETER( MaxNumFitPar = 38 )
      INTEGER MaxNumConstr
      PARAMETER( MaxNumConstr = 10 )
      INTEGER MaxNumComega
      PARAMETER( MaxNumComega = 8 )
      INTEGER cutN
      PARAMETER( cutN = 6 )
      INTEGER errN
      PARAMETER( errN = 17 )
      CHARACTER*40 messageCut(testcounter),messageErr(testcounter)
      LOGICAL    McFlag,EventOK,ErrFlag,truthmc,truthregg,truthomegaa,tr
     1  uthsemii,truththreee,truthelsee                                 
      INTEGER    nmctruth0,ntruth,nmctruth(7),ntruthregg,ntruthsemii,ntr
     1  uththreee,ntruthomegaa,ntruthelsee                              
      INTEGER    selected,analised
      INTEGER    counter(testcounter),counterBcg(testcounter)
      INTEGER    ErrFlagCount(testcounter),ErrFlagCountBcg(testcounter)
      LOGICAL    cutchar,cutneu,onlytruemc,makecuts
      COMMON / flagcommon / McFlag,EventOK,ErrFlag,truthmc
      COMMON / truthcommon / nmctruth0,ntruth,nmctruth
      COMMON / counterscommon / selected,analised,
     &                          counter,counterBcg
      COMMON / errcommon / ErrFlagCount,ErrFlagCountBcg,
     &                     messageCut,messageErr
      COMMON / talk / cutchar,cutneu,onlytruemc,makecuts
      SAVE /flagcommon/, /truthcommon/, /counterscommon/
      SAVE /errcommon/, /talk/
C======================== Include constans.inc =========================
      REAL       Mpip
      PARAMETER( Mpip = 139.57 )             
      REAL       Mpio
      PARAMETER( Mpio = 134.98 )             
      REAL       Mko
      PARAMETER( Mko  = 497.61 )             
      REAL       Momega
      PARAMETER( Momega  = 782.65 )          
      REAL       Cvel
      PARAMETER( Cvel = 29.9792458 )         
      REAL       EMCvel
      PARAMETER( EMCvel = 28.17 )            
      REAL       TauKs
      PARAMETER( TauKs = 0.08953 )   
      REAL       Mphi
      PARAMETER( Mphi  = 1019.460 )          
C= Include /kloe/soft/off/offline/inc/development/tls/maxstructdim.cin =
      INTEGER     NeleCluMax
      PARAMETER ( NeleCluMax = 2000 )
      INTEGER     MaxNumClu
      PARAMETER ( MaxNumClu = 100 )
      INTEGER    MaxNumVtx
      PARAMETER (MaxNumVtx = 20)
      INTEGER    MaxNumTrkV
      PARAMETER (MaxNumTrkV = 30)
      INTEGER    MaxNumTrk
      PARAMETER (MaxNumTrk = 100)
      INTEGER    MaxNumDHSP
      PARAMETER (MaxNumDHSP = 1000)
      Integer    nMaxDC
      Parameter (nMaxDC=1500)
      INTEGER    NQihiMax
      PARAMETER (NQihiMax=1000)
      INTEGER    NQcalMax
      PARAMETER (NQcalMax= 32 )
      INTEGER    MaxNumFirstHit
      PARAMETER (MaxNumFirstHit = 300 )
      INTEGER    MaxNtrkGen
      PARAMETER (MaxNtrkGen =50)
      INTEGER    MaxNvtxGen
      PARAMETER (MaxNvTxGen =30)
      integer TriggerElements
      parameter (TriggerElements = 300 )
C-----------------------------------------------------------------------
C
C External functions
C
C
C Local declarations
C
        INTEGER    M, N, ILoop, OffSet
        REAL       P(N), Ptrk(2,4), Pgam(4,4), Tgam(4)
        REAL       Ltot, PxTot, PyTot, PzTot, EnTot
        REAL       Mgggg, Mpipi, Vkne, Tkne
        REAL       Mkne,Vgam
             REAL       SinThetaBoost, CosThetaBoost
        REAL       PxPhi, PyPhi, PzPhi, EnPhi
        REAL*8     C(M)
C------------------------------------------------------------------------------
        Vgam = Cvel
        Mkne = Mko
C
C------------------------------------------------------------------------------
C Track related parameters
C------------------------------------------------------------------------------
C
C Px, Py, Py and E of the two pion from track parameters
C
        DO ILoop = 1,2
           Ptrk(ILoop,1) = (1000./ABS(P(3*(ILoop-1)+1)))*COS(P(3*(ILoop-
     &1)+2));
           Ptrk(ILoop,2) = (1000./ABS(P(3*(ILoop-1)+1)))*SIN(P(3*(ILoop-
     &1)+2));
           Ptrk(ILoop,3) = (1000./ABS(P(3*(ILoop-1)+1)))*P(3*(ILoop-1)+3
     &);
           Ptrk(ILoop,4) = SQRT( Ptrk(ILoop,1)**2 + Ptrk(ILoop,2)**2 +
     &          Ptrk(ILoop,3)**2 + Mpip**2 )
        ENDDO
C
C Invariant mass of pi+pi-
C
        PxTot = Ptrk(1,1) + Ptrk(2,1)
        PyTot = Ptrk(1,2) + Ptrk(2,2)
        PzTot = Ptrk(1,3) + Ptrk(2,3)
        EnTot = Ptrk(1,4) + Ptrk(2,4)
C
        Mpipi = SQRT( EnTot**2 - (PxTot**2+PyTot**2+PzTot**2) )
C
C------------------------------------------------------------------------------
C Photon related parameters
C------------------------------------------------------------------------------
C
C Px, Py, Pz, E and Time-of-Flight from cluster variables
C
        DO ILoop = 1,4
           OffSet = 5*(ILoop-1) + 6
C Cluster-NeuVtx distance
           Ltot = SQRT( (P(2+OffSet)-P(27))**2 +
     &          (P(3+OffSet)-P(28))**2 + (P(4+OffSet)-P(29))**2 )
           Pgam(ILoop,1) = P(1+OffSet) * (P(2+OffSet)-P(27))/Ltot
           Pgam(ILoop,2) = P(1+OffSet) * (P(3+OffSet)-P(28))/Ltot
           Pgam(ILoop,3) = P(1+OffSet) * (P(4+OffSet)-P(29))/Ltot
           Pgam(ILoop,4) = P(1+OffSet)
           Tgam(ILoop)   = Ltot / Vgam
        ENDDO
C
C Invariant mass of 4 photons
C
        PxTot = Pgam(1,1) + Pgam(2,1) + Pgam(3,1) + Pgam(4,1)
        PyTot = Pgam(1,2) + Pgam(2,2) + Pgam(3,2) + Pgam(4,2)
        PzTot = Pgam(1,3) + Pgam(2,3) + Pgam(3,3) + Pgam(4,3)
        EnTot = Pgam(1,4) + Pgam(2,4) + Pgam(3,4) + Pgam(4,4)
C
        Mgggg = SQRT( EnTot**2 - (PxTot**2+PyTot**2+PzTot**2) )
C
C Time-of-flight of the kaon decaying into photons
C
C NeuVtx-PhiVtx distance
        Ltot = SQRT( (P(27)-P(30))**2 +
     &       (P(28)-P(31))**2 + (P(29)-P(32))**2 )
        Vkne = Vgam * SQRT(PxTot**2+PyTot**2+PzTot**2) / EnTot
        Tkne = Ltot / Vkne
C Centroid-Apex distance - PROBABLY NEEDED, TO CHECK
C        Lcorr = SQRT( (P(27)-P(30))**2 +
C     &       (P(28)-P(31))**2 + (P(29)-P(32))**2 )
C        Vcorr = EMCvel
C        Tcorr = Ltot / Vcorr
C
C------------------------------------------------------------------------------
C Constraints
C------------------------------------------------------------------------------
C
C Total ToF for photons
C
        DO ILoop = 1,4
           OffSet = 5*(ILoop-1) + 6
           C(ILoop) = Tkne + Tgam(ILoop) - P(OffSet+5)
        ENDDO
C
C 4-momentum conservation
C
        PxTot = Ptrk(1,1) + Ptrk(2,1) + Pgam(1,1) + Pgam(2,1) +
     &       Pgam(3,1) + Pgam(4,1)
        PyTot = Ptrk(1,2) + Ptrk(2,2) + Pgam(1,2) + Pgam(2,2) +
     &       Pgam(3,2) + Pgam(4,2)
        PzTot = Ptrk(1,3) + Ptrk(2,3) + Pgam(1,3) + Pgam(2,3) +
     &       Pgam(3,3) + Pgam(4,3)
        EnTot = Ptrk(1,4) + Ptrk(2,4) + Pgam(1,4) + Pgam(2,4) +
     &       Pgam(3,4) + Pgam(4,4)
C
C        SinThetaBoost = P(38) / (P(36)+P(37))
C        CosThetaBoost = SQRT(1-SinThetaBoost**2)
C
C == SinThetaBoost * (P(36)+P(37))
        PxPhi = P(34)      
        PyPhi = P(35)
C eryk ERYK: in MC  PzPhi = 0
        PzPhi = P(36)
        EnPhi = P(33)
C
        C(5) = PxTot - PxPhi
        C(6) = PyTot - PyPhi
        C(7) = PzTot - PzPhi
        C(8) = EnTot - EnPhi
C
C Invariant masses
C
        C(9)  = Mgggg-Mkne
        C(10) = Mpipi-Mkne
C
        RETURN
        END
C==============================================================================
      subroutine fit_interf_omega(P,DP,M,CHISQR,C,D,MTX,Z,loopcount)
C==============================================================================
C
C==============================================================================
C
C Description:
C ------------
C
C Kinematic Fit for Phi --> KlKs --> pi+pi-pi0pi0 events
C
C Number of input parameters: 38
C * Curvature, cotg(theta), phi for each track   (3x2)
C * E, X, Y, Z, T for each photon                (5x4)
C * Charge, neutral and phi vertex coordinates   (3x3)
C * Beam Energies                                (1x2)
C * Phi momentum along x axis                    ( 1 )
C
C Total Constraints: 8
C * Time of flight of photons (4)
C * 4-momentum conservation   (4)
C
C Input parameters:
C -----------------
C M      = Number of constraints
C N      = Number of input parameters
C P (N)  = Measured Parameters
C DP(N)  = Parameter Errors
C
C Output parameters:
C ------------------
C CHISQR  = Chi2 value
C P (N)   = Adjusted Parameters
C
C==============================================================================
C
        IMPLICIT NONE
C
C========================= Include klspm00.inc =========================
      INTEGER NCLMIN
      PARAMETER( NCLMIN = 4 )
      REAL ECLMIN
      PARAMETER( ECLMIN = 20. )
      REAL vtxdist
      PARAMETER( vtxdist = 50. )
      INTEGER testcounter
      PARAMETER( testcounter = 100 )
      INTEGER MaxNumFitPar
      PARAMETER( MaxNumFitPar = 38 )
      INTEGER MaxNumConstr
      PARAMETER( MaxNumConstr = 10 )
      INTEGER MaxNumComega
      PARAMETER( MaxNumComega = 8 )
      INTEGER cutN
      PARAMETER( cutN = 6 )
      INTEGER errN
      PARAMETER( errN = 17 )
      CHARACTER*40 messageCut(testcounter),messageErr(testcounter)
      LOGICAL    McFlag,EventOK,ErrFlag,truthmc,truthregg,truthomegaa,tr
     1  uthsemii,truththreee,truthelsee                                 
      INTEGER    nmctruth0,ntruth,nmctruth(7),ntruthregg,ntruthsemii,ntr
     1  uththreee,ntruthomegaa,ntruthelsee                              
      INTEGER    selected,analised
      INTEGER    counter(testcounter),counterBcg(testcounter)
      INTEGER    ErrFlagCount(testcounter),ErrFlagCountBcg(testcounter)
      LOGICAL    cutchar,cutneu,onlytruemc,makecuts
      COMMON / flagcommon / McFlag,EventOK,ErrFlag,truthmc
      COMMON / truthcommon / nmctruth0,ntruth,nmctruth
      COMMON / counterscommon / selected,analised,
     &                          counter,counterBcg
      COMMON / errcommon / ErrFlagCount,ErrFlagCountBcg,
     &                     messageCut,messageErr
      COMMON / talk / cutchar,cutneu,onlytruemc,makecuts
      SAVE /flagcommon/, /truthcommon/, /counterscommon/
      SAVE /errcommon/, /talk/
C======= Include /kloe/soft/off/offline/inc/development/bcs.cin ========
        INTEGER          IW
        REAL             RW(200)
        DOUBLE PRECISION DW(100)
        CHARACTER*4      AW(200)
        INTEGER*2        IW2(400)
        EQUIVALENCE (IW,IW2,RW,DW)
        EQUIVALENCE (IW,AW)
        COMMON /BCS/ IW(200)
C====== Include /kloe/soft/off/offline/inc/development/jobsta.cin ======
        COMMON /JOBSTI/ NRUN, NEV, NRUNA, NEVA, IBEGRN, IENDRN,
     1                  IUDBUG, NEVP, NEVPR, NEVPA, NREC, NRECF
        INTEGER         NRUN, NEV, NRUNA, NEVA, IBEGRN, IENDRN
        INTEGER         IUDBUG, NEVP, NEVPR, NEVPA, NREC, NRECF
        COMMON /JOBLRI/ EXPTYP, PARTNO, RUNTYP, LRECNO, NRCTHS,
     1                  CDFCID, IRTYPA, VERNUM, RELNUM, CRPRID,
     2                  CRYEAR,  CRMON,  CRDAY,
     3                  CRHOUR,  CRMIN,  CRSEC
        INTEGER         EXPTYP, PARTNO, RUNTYP, LRECNO, NRCTHS
        INTEGER         CDFCID, IRTYPA, VERNUM, RELNUM, CRPRID
        INTEGER         CRYEAR,  CRMON,  CRDAY
        INTEGER         CRHOUR,  CRMIN,  CRSEC
        COMMON /JOBEVC/ TRIGNO, TRGSUM(4),
     1                  STEVTI, SDONET, RPTMOD, RESERV, SOPTYP, SBUFNO,
     2                  SREADO(2), SRESIO(2), CREADO(2), TERSUM(2)
        INTEGER         TRIGNO,    TRGSUM
        INTEGER         STEVTI, SDONET, RPTMOD, RESERV, SOPTYP, SBUFNO
        INTEGER         SREADO,    SRESIO,    CREADO,    TERSUM
        COMMON /JOBSTR/ ROOTS, BMAGNT, BMAGFT, BMAGBT
        REAL            ROOTS, BMAGNT, BMAGFT, BMAGBT
C====== Include /kloe/soft/off/offline/inc/development/erlevl.cin ======
        INTEGER    ERMINI
        PARAMETER (ERMINI = 0)
        INTEGER    ERMAXI
        PARAMETER (ERMAXI = 7)
        INTEGER    ERSUCC
        PARAMETER (ERSUCC = 0)
        INTEGER    ERINFO
        PARAMETER (ERINFO = 1)
        INTEGER    ERWARN
        PARAMETER (ERWARN = 2)
        INTEGER    EREROR
        PARAMETER (EREROR = 3)
        INTEGER    ERSEVR
        PARAMETER (ERSEVR = 4)
        INTEGER         ELSUCC
        PARAMETER       (ELSUCC = 10)
        INTEGER         ELINFO
        PARAMETER       (ELINFO = 11)
        INTEGER         ELWARN
        PARAMETER       (ELWARN = 12)
        INTEGER         ELWARN2
        PARAMETER       (ELWARN2 = 13)
        INTEGER         ELEROR
        PARAMETER       (ELEROR = 14)
        INTEGER         ELEROR2
        PARAMETER       (ELEROR2 = 15)
        INTEGER         ELNEXT
        PARAMETER       (ELNEXT = 16)
        INTEGER         ELNEXT2
        PARAMETER       (ELNEXT2 = 17)
        INTEGER         ELSEVR
        PARAMETER       (ELSEVR = 18)
        INTEGER         ELABOR
        PARAMETER       (ELABOR = 19)
C====== Include /kloe/soft/off/offline/inc/development/oferco.cin ======
      INTEGER    OFSUCC                 
      PARAMETER (OFSUCC = 0)
      INTEGER    OFFAIL                 
      PARAMETER (OFFAIL = 1)
      INTEGER    OFFINF                 
      PARAMETER (OFFINF = 2)
      INTEGER    OFERRE                 
      PARAMETER (OFERRE = 3)            
      INTEGER    OFFEOF                 
      PARAMETER (OFFEOF = 4)
C==== Include /kloe/soft/off/ybos/production/aix/library/errcod.cin ====
      INTEGER    YESUCC
      PARAMETER (YESUCC = '08508009'X)
      INTEGER    YEWLOK
      PARAMETER (YEWLOK = '08508083'X)
      INTEGER    YEWRNS
      PARAMETER (YEWRNS = '0850808B'X)
      INTEGER    YEINDX
      PARAMETER (YEINDX = '08508093'X)
      INTEGER    YEEOF
      PARAMETER (YEEOF  = '08508100'X)
      INTEGER    YEPRBD
      PARAMETER (YEPRBD = '08508108'X)
      INTEGER    YERSTF
      PARAMETER (YERSTF = '08508110'X)
      INTEGER    YEOVRF
      PARAMETER (YEOVRF = '08508182'X)
      INTEGER    YEBKNG
      PARAMETER (YEBKNG = '0850818A'X)
      INTEGER    YEBKTS
      PARAMETER (YEBKTS = '08508192'X)
      INTEGER    YEBKOV
      PARAMETER (YEBKOV = '0850819A'X)
      INTEGER    YEBKSP
      PARAMETER (YEBKSP = '085081A2'X)
      INTEGER    YENOBK
      PARAMETER (YENOBK = '085081AA'X)
      INTEGER    YENONR
      PARAMETER (YENONR = '085081B2'X)
      INTEGER    YEMANY
      PARAMETER (YEMANY = '085081BA'X)
      INTEGER    YEDUPB
      PARAMETER (YEDUPB = '085081C2'X)
      INTEGER    YEBSIL
      PARAMETER (YEBSIL = '085081CA'X)
      INTEGER    YEILOP
      PARAMETER (YEILOP = '085081D2'X)
      INTEGER    YEAROF
      PARAMETER (YEAROF = '085081DA'X)
      INTEGER    YETRUN
      PARAMETER (YETRUN = '085081E2'X)
      INTEGER    YECORR
      PARAMETER (YECORR = '085081EA'X)
      INTEGER    YEWRGI
      PARAMETER (YEWRGI = '085081F2'X)
      INTEGER    YELWIU
      PARAMETER (YELWIU = '085081FA'X)
      INTEGER    YEILUN
      PARAMETER (YEILUN = '08508202'X)
      INTEGER    YEIOUS
      PARAMETER (YEIOUS = '0850820A'X)
      INTEGER    YELROF
      PARAMETER (YELROF = '08508212'X)
      INTEGER    YELRLN
      PARAMETER (YELRLN = '0850821A'X)
      INTEGER    YEPRXL
      PARAMETER (YEPRXL = '08508222'X)
      INTEGER    YEBKOS
      PARAMETER (YEBKOS = '0850822A'X)
      INTEGER    YEBKXD
      PARAMETER (YEBKXD = '08508232'X)
      INTEGER    YELRZL
      PARAMETER (YELRZL = '0850823A'X)
      INTEGER    YELRWT
      PARAMETER (YELRWT = '08508242'X)
      INTEGER    YEDKOP
      PARAMETER (YEDKOP = '0850824A'X)
      INTEGER    YEDKWT
      PARAMETER (YEDKWT = '08508252'X)
      INTEGER    YEDKRD
      PARAMETER (YEDKRD = '0850825A'X)
      INTEGER    YEILUS
      PARAMETER (YEILUS = '08508282'X)
      INTEGER    YERSUS
      PARAMETER (YERSUS = '0850828A'X)
      INTEGER    YENRUS
      PARAMETER (YENRUS = '08508292'X)
      INTEGER    YEILLA
      PARAMETER (YEILLA = '08508302'X)
      INTEGER    YEDUPA
      PARAMETER (YEDUPA = '0850830A'X)
      INTEGER    YEBTIL
      PARAMETER (YEBTIL = '08508382'X)
      INTEGER    YEGRIL
      PARAMETER (YEGRIL = '0850838A'X)
      INTEGER    YECHKF
      PARAMETER (YECHKF = '08508392'X)
      INTEGER    YESTOF
      PARAMETER (YESTOF = '0850839A'X)
      INTEGER    YENOTA
      PARAMETER (YENOTA = '085083A2'X)
      INTEGER    YEBNKP
      PARAMETER (YEBNKP = '085083AA'X)
      INTEGER    YEBDCH
      PARAMETER (YEBDCH = '085083C2'X)
      INTEGER    YEBDBK
      PARAMETER (YEBDBK = '085083CA'X)
      INTEGER    YEBDLU
      PARAMETER (YEBDLU = '085083D2'X)
      INTEGER    YEBDBS
      PARAMETER (YEBDBS = '085083DA'X)
C==== Include /kloe/soft/off/a_c/production/aix/library/anerror.cin ====
      INTEGER    ANSUCC
      PARAMETER (ANSUCC = '8198009'X)
      INTEGER    ANBEGN
      PARAMETER (ANBEGN = '8198019'X)
      INTEGER    ANPDUN
      PARAMETER (ANPDUN = '8198021'X)
      INTEGER    ANAST
      PARAMETER (ANAST   = '81980A3'X)
      INTEGER    ANILPS
      PARAMETER (ANILPS = '8198122'X)
      INTEGER    ANILMD
      PARAMETER (ANILMD = '819812A'X)
      INTEGER    ANLSTR
      PARAMETER (ANLSTR = '81980C0'X)
      INTEGER    ANNINI
      PARAMETER (ANNINI = '8198102'X)
      INTEGER    ANOVRF
      PARAMETER (ANOVRF = '819810A'X)
      INTEGER    ANNOME
      PARAMETER (ANNOME = '8198112'X)
      INTEGER    ANNOLU
      PARAMETER (ANNOLU = '819811A'X)
      INTEGER    ANNOEP
      PARAMETER (ANNOEP = '8198011'X)
      INTEGER    ANILGP
      PARAMETER (ANILGP = '8198132'X)
      INTEGER    ANILPA
      PARAMETER (ANILPA = '819813A'X)
      INTEGER    ANILQU
      PARAMETER (ANILQU = '8198142'X)
      INTEGER    ANILSY
      PARAMETER (ANILSY = '819814A'X)
      INTEGER    ANILLU
      PARAMETER (ANILLU = '8198152'X)
      INTEGER    ANILVB
      PARAMETER (ANILVB = '819815A'X)
      INTEGER    ANILLI
      PARAMETER (ANILLI = '8198162'X)
      INTEGER    ANDUMD
      PARAMETER (ANDUMD = '819816A'X)
      INTEGER    ANYBER
      PARAMETER (ANYBER = '8198172'X)
      INTEGER    ANYBFE
      PARAMETER (ANYBFE = '081981B2'X)
      INTEGER    ANEFOP
      PARAMETER (ANEFOP = '819817A'X)
      INTEGER    ANEOF
      PARAMETER (ANEOF = '8198083'X)
      INTEGER    ANIFNO
      PARAMETER (ANIFNO = '819808B'X)
      INTEGER    ANCLOS
      PARAMETER (ANCLOS = '8198093'X)
      INTEGER    ANOFNO
      PARAMETER (ANOFNO = '819809B'X)
      INTEGER    ANNOLR
      PARAMETER (ANNOLR = '8198182'X)
      INTEGER    ANILST
      PARAMETER (ANILST = '819818A'X)
      INTEGER    ANNOEV
      PARAMETER (ANNOEV = '8198192'X)
      INTEGER    ANILBL
      PARAMETER (ANILBL = '819819A'X)
      INTEGER    ANNEIF
      PARAMETER (ANNEIF = '81981A2'X)
      INTEGER    ANNOBK
      PARAMETER (ANNOBK = '81980C8'X)
      INTEGER    ANILBK
      PARAMETER (ANILBK = '81981AA'X)
C
C External functions
C
        REAL VDOTN, VMOD
C
C Local declarations
C
        REAL XI2
        COMMON/PHIFIT/XI2(50)
        INTEGER I, J, K, L, N
C Number of parameters
        PARAMETER (N=33)                  
C Number of constraints in fit
        INTEGER M                         
C Parameters & errors
        REAL    P(N), DP(N), P0(N), P1(N) 
C Parameter Corrections
        REAL*8  A(N)                      
C Derivative of constraint M w.r.t. par. N
        REAL*8  D(N,M)                    
C Balance of constraints
        REAL*8  C(M)                      
C Covariance Matrix (?)
        REAL*8  MTX(N+M,N+M)              
C (0,-C)
        REAL*8  Z(N+M)                    
        REAL*8  C1(M)
        INTEGER ITERATIONMAX, IFAIL,loopcount
        REAL    CHISQR, DX2, DX2MAX, R(100), C2
        REAL    PTEMP(N)
        LOGICAL FailFlag
C
C==============================================================================
C
        FailFlag = .FALSE.
        ITERATIONMAX = 50
        CHISQR = 0.
        DX2 = 1.0e6
        DX2MAX = 0.01
        DO I = 1,N
          A(I)  = 0.
          P0(I) = P(I)
          Z(I)  = 0.
          DO J = 1,M
            D(I,J) = 0.d0
          ENDDO
        ENDDO
C
        CALL CONSTRAINT_OMEGA(M,N,P,C)
C
        loopcount = 0
        DO WHILE( loopcount.LT.ITERATIONMAX .AND. DX2.GT.DX2MAX )
          loopcount = loopcount+1
          DO I = 1,N
            P1(I) = P(I)
          ENDDO
         CALL DERIVATIVE_OMEGA(N,M,D,C,P,DP,PTEMP,C1)
          DO J = 1,M
            Z(N+J) = -C(J)
          ENDDO
          DO I = 1,N+M
            DO J = I,N+M
              IF( I.GT.N )THEN
                MTX(I,J) = 0.
              ELSEIF( J.GT.N )THEN
                MTX(I,J) = D(I,J-N)
              ELSEIF( I.EQ.J .AND. DP(I).GT.0 )THEN
                MTX(I,J) = 1/(DP(I)*DP(I))
              ELSE
                MTX(I,J) = 0.
              ENDIF
              IF( J.GT.I )THEN
                MTX(J,I) = MTX(I,J)
              ENDIF
            ENDDO
          ENDDO
          CALL DINV(N+M,MTX,N+M,R,IFAIL)
          IF(IFAIL.NE.0) FailFlag = .TRUE.
          DX2 = CHISQR
          CHISQR = 0.
          DO I = 1,N
            A(I) = 0.
            DO J = 1,N+M
              A(I) = A(I) + MTX(I,J)*Z(J)
            ENDDO
C            P(I)= P1(I) + A(I)
            IF( ABS(A(I)).LT.1.0E3 )THEN
              P(I) = P1(I) + A(I)
            ENDIF
C--- To avoid negative energy photons... ------------------
Cc           IF( MOD(I+4,5).EQ.0 )THEN
Cc             P(I) = ABS(P(I))
Cc           ENDIF
            IF( I.EQ.7.or.I.EQ.12.or.I.EQ.17.or.I.EQ.22 )THEN
              IF( P(I).LT. ECLMIN ) THEN
                P(I) = ECLMIN
              ENDIF
            ENDIF
C----------------------------------------------------------
            IF( DP(I).GT.0 )THEN
              CHISQR = CHISQR + ((P(I)-P0(I))/DP(I))**2
            ENDIF
          ENDDO
          DX2 = ABS(CHISQR-DX2)
          CALL CONSTRAINT_OMEGA(M,N,P,C)
          C2 = 0.
          DO J = 1,M
            C2 = C2 + C(J)*C(J)
          ENDDO
          XI2(loopcount) = CHISQR
        ENDDO
C
        IF( FailFlag ) loopcount = 100 + loopcount
C
        RETURN
        END
C
C==============================================================================
        SUBROUTINE DERIVATIVE_OMEGA(N,M,D,C,P,DP,P1,C1)
C==============================================================================
C
        IMPLICIT NONE
C
C Local declarations
C
        INTEGER I,J,K,L,M,N
        REAL    P(N), DP(N), P1(N)
        REAL*8  D(N,M), C(M), C1(M), D1, STEP
C
C==============================================================================
C
        DO I=1,N
          P1(I)=P(I)
        ENDDO
        DO J=1,M
          DO I=1,N
            IF(DP(I).EQ.0)THEN
              IF(P1(I).EQ.0)THEN
                STEP=0.1
              ELSE
                STEP=0.001*P1(I)
              ENDIF
            ELSE
              STEP=0.01*DP(I)
            ENDIF
            L=0
            DO WHILE(ABS(D1-D(I,J)).GT.0.001*ABS(D(I,J)).OR.L.LE.1)
              D1=D(I,J)
              P1(I)=P(I)+STEP
              CALL CONSTRAINT_OMEGA(M,N,P1,C1)
              D(I,J)=(C1(J)-C(J))/STEP
              STEP=STEP/2
              L=L+1
            ENDDO
            P1(I)=P(I)
          ENDDO
        ENDDO
C
        RETURN
        END
C
C==============================================================================
        SUBROUTINE CONSTRAINT_OMEGA(M,N,P,C)
C==============================================================================
C-----------------------------------------------------------------------
C========================= Include klspm00.inc =========================
      INTEGER NCLMIN
      PARAMETER( NCLMIN = 4 )
      REAL ECLMIN
      PARAMETER( ECLMIN = 20. )
      REAL vtxdist
      PARAMETER( vtxdist = 50. )
      INTEGER testcounter
      PARAMETER( testcounter = 100 )
      INTEGER MaxNumFitPar
      PARAMETER( MaxNumFitPar = 38 )
      INTEGER MaxNumConstr
      PARAMETER( MaxNumConstr = 10 )
      INTEGER MaxNumComega
      PARAMETER( MaxNumComega = 8 )
      INTEGER cutN
      PARAMETER( cutN = 6 )
      INTEGER errN
      PARAMETER( errN = 17 )
      CHARACTER*40 messageCut(testcounter),messageErr(testcounter)
      LOGICAL    McFlag,EventOK,ErrFlag,truthmc,truthregg,truthomegaa,tr
     1  uthsemii,truththreee,truthelsee                                 
      INTEGER    nmctruth0,ntruth,nmctruth(7),ntruthregg,ntruthsemii,ntr
     1  uththreee,ntruthomegaa,ntruthelsee                              
      INTEGER    selected,analised
      INTEGER    counter(testcounter),counterBcg(testcounter)
      INTEGER    ErrFlagCount(testcounter),ErrFlagCountBcg(testcounter)
      LOGICAL    cutchar,cutneu,onlytruemc,makecuts
      COMMON / flagcommon / McFlag,EventOK,ErrFlag,truthmc
      COMMON / truthcommon / nmctruth0,ntruth,nmctruth
      COMMON / counterscommon / selected,analised,
     &                          counter,counterBcg
      COMMON / errcommon / ErrFlagCount,ErrFlagCountBcg,
     &                     messageCut,messageErr
      COMMON / talk / cutchar,cutneu,onlytruemc,makecuts
      SAVE /flagcommon/, /truthcommon/, /counterscommon/
      SAVE /errcommon/, /talk/
C======================== Include constans.inc =========================
      REAL       Mpip
      PARAMETER( Mpip = 139.57 )             
      REAL       Mpio
      PARAMETER( Mpio = 134.98 )             
      REAL       Mko
      PARAMETER( Mko  = 497.61 )             
      REAL       Momega
      PARAMETER( Momega  = 782.65 )          
      REAL       Cvel
      PARAMETER( Cvel = 29.9792458 )         
      REAL       EMCvel
      PARAMETER( EMCvel = 28.17 )            
      REAL       TauKs
      PARAMETER( TauKs = 0.08953 )   
      REAL       Mphi
      PARAMETER( Mphi  = 1019.460 )          
C= Include /kloe/soft/off/offline/inc/development/tls/maxstructdim.cin =
      INTEGER     NeleCluMax
      PARAMETER ( NeleCluMax = 2000 )
      INTEGER     MaxNumClu
      PARAMETER ( MaxNumClu = 100 )
      INTEGER    MaxNumVtx
      PARAMETER (MaxNumVtx = 20)
      INTEGER    MaxNumTrkV
      PARAMETER (MaxNumTrkV = 30)
      INTEGER    MaxNumTrk
      PARAMETER (MaxNumTrk = 100)
      INTEGER    MaxNumDHSP
      PARAMETER (MaxNumDHSP = 1000)
      Integer    nMaxDC
      Parameter (nMaxDC=1500)
      INTEGER    NQihiMax
      PARAMETER (NQihiMax=1000)
      INTEGER    NQcalMax
      PARAMETER (NQcalMax= 32 )
      INTEGER    MaxNumFirstHit
      PARAMETER (MaxNumFirstHit = 300 )
      INTEGER    MaxNtrkGen
      PARAMETER (MaxNtrkGen =50)
      INTEGER    MaxNvtxGen
      PARAMETER (MaxNvTxGen =30)
      integer TriggerElements
      parameter (TriggerElements = 300 )
C-----------------------------------------------------------------------
C
        IMPLICIT NONE
C
C External functions
C
C
C Local declarations
C
        INTEGER    M, N, ILoop, OffSet
        REAL       P(N), Ptrk(2,4), Pgam(4,4), Tgam(4)
        REAL       Ltot, PxTot, PyTot, PzTot, EnTot
        REAL       Mgggg, Mpipi, Vkne, Tkne, Vgam
        REAL       Mkne
             REAL       SinThetaBoost, CosThetaBoost
        REAL       PxPhi, PyPhi, PzPhi, EnPhi
        REAL*8     C(M)
C
C==============================================================================
        Vgam = Cvel
        Mkne = Mko
C
C------------------------------------------------------------------------------
C Track related parameters
C------------------------------------------------------------------------------
C
C Px, Py, Py and E of the two pion from track parameters
C
        DO ILoop = 1,2
           Ptrk(ILoop,1) = (1000./ABS(P(3*(ILoop-1)+1)))*COS(P(3*(ILoop-
     &1)+2));
           Ptrk(ILoop,2) = (1000./ABS(P(3*(ILoop-1)+1)))*SIN(P(3*(ILoop-
     &1)+2));
           Ptrk(ILoop,3) = (1000./ABS(P(3*(ILoop-1)+1)))*P(3*(ILoop-1)+3
     &);
           Ptrk(ILoop,4) = SQRT( Ptrk(ILoop,1)**2 + Ptrk(ILoop,2)**2 +
     &          Ptrk(ILoop,3)**2 + Mpip**2 )
        ENDDO
C
C------------------------------------------------------------------------------
C Photon related parameters
C------------------------------------------------------------------------------
C
C Px, Py, Pz, E and Time-of-Flight from cluster variables
C
        DO ILoop = 1,4
           OffSet = 5*(ILoop-1) + 6
C Cluster-PhiVtx distance
           Ltot = SQRT( (P(2+OffSet)-P(27))**2 +
     &          (P(3+OffSet)-P(28))**2 + (P(4+OffSet)-P(29))**2 )
           Pgam(ILoop,1) = P(1+OffSet) * (P(2+OffSet)-P(27))/Ltot
           Pgam(ILoop,2) = P(1+OffSet) * (P(3+OffSet)-P(28))/Ltot
           Pgam(ILoop,3) = P(1+OffSet) * (P(4+OffSet)-P(29))/Ltot
           Pgam(ILoop,4) = P(1+OffSet)
           Tgam(ILoop)   = Ltot / Vgam
        ENDDO
C
C------------------------------------------------------------------------------
C Constraints
C------------------------------------------------------------------------------
C
C ToF for photons
C
        DO ILoop = 1,4
           OffSet = 5*(ILoop-1) + 10
           C(ILoop) = Tgam(ILoop) - P(OffSet+1)
        ENDDO
C
C 4-momentum conservation
C
        PxTot = Ptrk(1,1) + Ptrk(2,1) + Pgam(1,1) + Pgam(2,1) +
     &       Pgam(3,1) + Pgam(4,1)
        PyTot = Ptrk(1,2) + Ptrk(2,2) + Pgam(1,2) + Pgam(2,2) +
     &       Pgam(3,2) + Pgam(4,2)
        PzTot = Ptrk(1,3) + Ptrk(2,3) + Pgam(1,3) + Pgam(2,3) +
     &       Pgam(3,3) + Pgam(4,3)
        EnTot = Ptrk(1,4) + Ptrk(2,4) + Pgam(1,4) + Pgam(2,4) +
     &       Pgam(3,4) + Pgam(4,4)
C
C        SinThetaBoost = P(38) / (P(36)+P(37))
C        CosThetaBoost = SQRT(1-SinThetaBoost**2)
C
C == SinThetaBoost * (P(36)+P(37))
        PxPhi = P(31)      
        PyPhi = P(32)
        PzPhi = P(33)
        EnPhi = P(30)
C
        C(5) = PxTot - PxPhi
        C(6) = PyTot - PyPhi
        C(7) = PzTot - PzPhi
        C(8) = EnTot - EnPhi
C
        RETURN
        END
C=======================================================================
C=======================================================================
C=======================================================================
      SUBROUTINE do_cosinus(p,ix,iy,iz,cosinus)
C-----------------------------------------------------------------------
C returns cosinus of angle between two vectors:
C one from the momenta of 9-vector,
C second one from two points (decay and IP)
C-----------------------------------------------------------------------
        IMPLICIT NONE
C-----------------------------------------------------------------------
C
C Local declarations
C
      REAL p(9),ix,iy,iz,cosinus
      REAL d(3),pmod,ddmod
      INTEGER i
C-----------------------------------------------------------------------
C in: p(9),ix,iy,iz
C out: cosinus
C-----------------------------------------------------------------------
      d(1) = p(7) - ix
      d(2) = p(8) - iy
      d(3) = p(9) - iz
      pmod =  sqrt( p(1)**2 + p(2)**2 + p(3)**2 )
      ddmod = sqrt( d(1)**2 + d(2)**2 + d(3)**2 )
      if ((pmod*ddmod).ne.0.) then
         cosinus = ( p(1)*d(1) + p(2)*d(2) + p(3)*d(3)) / (pmod*ddmod)
      else
         cosinus = -2.
      endif
      END
C=======================================================================
C=======================================================================
C=======================================================================
      SUBROUTINE do_distance(pp,ix,iy,iz,distance)
C-----------------------------------------------------------------------
C distance between point i-xyz (IP) and
C a vector (momentum) going through decay point
C-----------------------------------------------------------------------
        IMPLICIT NONE
C-----------------------------------------------------------------------
C
C Local declarations
C
      REAL pp(9),ix,iy,iz,x(3),p(3),distance
      REAL v1,v2,v3,VxP1,VxP2,VxP3,modvxp
C-----------------------------------------------------------------------
C in: pp(9),ix,iy,iz
C out: distance
C-----------------------------------------------------------------------
      x(1)=pp(7)
      x(2)=pp(8)
      x(3)=pp(9)
      p(1)=pp(1)
      p(2)=pp(2)
      p(3)=pp(3)
      v1=x(1)-ix
      v2=x(2)-iy
      v3=x(3)-iz
      VxP1=v2*p(3)-v3*p(2)
      VxP2=v3*p(1)-v1*p(3)
      VxP3=v1*p(2)+v2*p(1)
      modvxp=sqrt(VxP1**2+VxP2**2+VxP3**2)
      IF( sqrt(p(1)**2+p(2)**2+p(3)**2).ne.0 ) THEN
        distance = modvxp/sqrt(p(1)**2+p(2)**2+p(3)**2)
      ELSE
        distance = -2.
      ENDIF
      END
C=======================================================================
C=======================================================================
C=======================================================================
      SUBROUTINE do_cosin(p1,p2,cosin)
C-----------------------------------------------------------------------
C returns cosinus of angle between two momentum vectors
C-----------------------------------------------------------------------
        IMPLICIT NONE
C-----------------------------------------------------------------------
C
C Local declarations
C
      REAL p1(9),p2(9),cosin
      REAL pmod,ddmod
      INTEGER i
C-----------------------------------------------------------------------
C in: p1(9),p2(9)
C out: cosin
C-----------------------------------------------------------------------
      pmod =  sqrt( p1(1)**2 + p1(2)**2 + p1(3)**2 )
      ddmod = sqrt( p2(1)**2 + p2(2)**2 + p2(3)**2 )
      if ((pmod*ddmod).ne.0.) then
         cosin = ( p1(1)*p2(1) + p1(2)*p2(2) + p1(3)*p2(3)) / (pmod*ddmo
     &d)
      else
         cosin = -2.
      endif
      END
C=======================================================================
C=======================================================================
C=======================================================================
      SUBROUTINE OmegaTopologySel(nv,nclu,ntv,Bx,By,Bz,CurV,PhiV,CotV,
     &          EneCl,Xcl,Ycl,Zcl,Tcl,Xacl,Yacl,Zacl,xv,yv,zv,vtx,iv,
     &          chtrk_index,clu_index,chvtx_index,omegatruthtop)
C=======================================================================
        IMPLICIT NONE
C======================== Include constans.inc =========================
      REAL       Mpip
      PARAMETER( Mpip = 139.57 )             
      REAL       Mpio
      PARAMETER( Mpio = 134.98 )             
      REAL       Mko
      PARAMETER( Mko  = 497.61 )             
      REAL       Momega
      PARAMETER( Momega  = 782.65 )          
      REAL       Cvel
      PARAMETER( Cvel = 29.9792458 )         
      REAL       EMCvel
      PARAMETER( EMCvel = 28.17 )            
      REAL       TauKs
      PARAMETER( TauKs = 0.08953 )   
      REAL       Mphi
      PARAMETER( Mphi  = 1019.460 )          
C=======================================================================
      INTEGER nv, nclu, ntv, vtx(nv), iv(nv),I,J
      REAL CurV(nv),PhiV(nv),CotV(nv),EneCl(nclu),Xcl(nclu),Ycl(nclu),
     &          Zcl(nclu),Tcl(nclu),Xacl(nclu),Yacl(nclu),Zacl(nclu),
     &          Bx,By,Bz,xv(nv),yv(nv),zv(nv)
      REAL chtrk(3,nv),clvars(8,nclu),possys(3)
      REAL fid_radius(nv),fid_height(nv),time_res(nclu),
     &         clus_time(nclu),time_corr(nclu)
      LOGICAL fail
      INTEGER omegatruthtop,chosen,index,
     &         proper_clusters,clu_index_tmp(50),clu_index(4),
     &         chtrk_index(2),chtrk_index_tmp(50),chvtx_index
      fail = .FALSE.
      omegatruthtop = 0
      chosen = 0
      index = 0
      proper_clusters = 0
      DO I=1,nv
        chtrk(1,I) = xv(I)
        chtrk(2,I) = yv(I)
        chtrk(3,I) = zv(I)
      ENDDO
      DO I=1,nclu
        clvars(1,I) = EneCl(I)
        clvars(2,I) = Xcl(I)
        clvars(3,I) = Ycl(I)
        clvars(4,I) = Zcl(I)
        clvars(5,I) = Tcl(I)
        clvars(6,I) = Xacl(I)
        clvars(7,I) = Yacl(I)
        clvars(8,I) = Zacl(I)
      ENDDO
      possys(1) = Bx
      possys(2) = By
      possys(3) = Bz
C=======================================================================
C=======================================================================
      DO I=1,nv
         fid_radius(I) = SQRT( (chtrk(1,I) - possys(1))**2
     &             + (chtrk(2,I) - possys(2))**2 )
         fid_height(I) = ABS( chtrk(3,I) - possys(3) )
         IF(fid_radius(I).LT.2..AND.fid_height(I).LT.1..AND.vtx(I).EQ.2)
     &THEN
            DO J=1,ntv
               IF(iv(I).EQ.J+1)THEN
                  chtrk_index_tmp(I) = I
               ENDIF
            ENDDO
            index = I
            chosen = chosen + 1
         ENDIF
      ENDDO
      IF(chosen.NE.1)THEN
         fail = .TRUE.
      ENDIF
      IF(fail.NEQV..TRUE.)THEN
         DO I=1,nclu
            time_res(I) = SQRT( (0.054/SQRT(clvars(1,I)/1000.))**2 + (0.
     1  147)**2 )                                                       
            clus_time(I) = SQRT( (clvars(6,I) - chtrk(1,index))**2 +
     &                  (clvars(7,I) - chtrk(2,index))**2 +
     &                  (clvars(8,I) - chtrk(3,index))**2)
            clus_time(I) = clus_time(I)/Cvel
            time_corr(I) = SQRT( (clvars(2,I) - clvars(6,I))**2 +
     &                  (clvars(3,I) - clvars(7,I))**2 +
     &                  (clvars(4,I) - clvars(8,I))**2)
            time_corr(I) = time_corr(I)/EMCvel
            IF(clus_time(I) + time_corr(I) - clvars(5,I) <
     &          MIN(3.*time_res(I), 1.).AND.clvars(1,I) > 20.)THEN
               proper_clusters = proper_clusters + 1
               clu_index_tmp(I) = I
            ENDIF
         ENDDO
         IF(proper_clusters.NE.4)THEN
            fail = .TRUE.
         ENDIF
      ENDIF
      IF(fail .NEQV. .TRUE.) THEN
         DO I=1,4
            clu_index(I) = clu_index_tmp(I)
         ENDDO
         DO I=1,2
            chtrk_index(I) = chtrk_index_tmp(I)
         ENDDO
         chvtx_index = index
         omegatruthtop = 1
      ELSE
         DO I=1,4
            clu_index(I) = -999
         ENDDO
         DO I=1,2
            chtrk_index(I) = -999
         ENDDO
         chvtx_index = -999
         omegatruthtop = 0
      ENDIF
      END
C=======================================================================
C=======================================================================
C=======================================================================
      SUBROUTINE statisticss(id)
        IMPLICIT NONE
C========================= Include klspm00.inc =========================
      INTEGER NCLMIN
      PARAMETER( NCLMIN = 4 )
      REAL ECLMIN
      PARAMETER( ECLMIN = 20. )
      REAL vtxdist
      PARAMETER( vtxdist = 50. )
      INTEGER testcounter
      PARAMETER( testcounter = 100 )
      INTEGER MaxNumFitPar
      PARAMETER( MaxNumFitPar = 38 )
      INTEGER MaxNumConstr
      PARAMETER( MaxNumConstr = 10 )
      INTEGER MaxNumComega
      PARAMETER( MaxNumComega = 8 )
      INTEGER cutN
      PARAMETER( cutN = 6 )
      INTEGER errN
      PARAMETER( errN = 17 )
      CHARACTER*40 messageCut(testcounter),messageErr(testcounter)
      LOGICAL    McFlag,EventOK,ErrFlag,truthmc,truthregg,truthomegaa,tr
     1  uthsemii,truththreee,truthelsee                                 
      INTEGER    nmctruth0,ntruth,nmctruth(7),ntruthregg,ntruthsemii,ntr
     1  uththreee,ntruthomegaa,ntruthelsee                              
      INTEGER    selected,analised
      INTEGER    counter(testcounter),counterBcg(testcounter)
      INTEGER    ErrFlagCount(testcounter),ErrFlagCountBcg(testcounter)
      LOGICAL    cutchar,cutneu,onlytruemc,makecuts
      COMMON / flagcommon / McFlag,EventOK,ErrFlag,truthmc
      COMMON / truthcommon / nmctruth0,ntruth,nmctruth
      COMMON / counterscommon / selected,analised,
     &                          counter,counterBcg
      COMMON / errcommon / ErrFlagCount,ErrFlagCountBcg,
     &                     messageCut,messageErr
      COMMON / talk / cutchar,cutneu,onlytruemc,makecuts
      SAVE /flagcommon/, /truthcommon/, /counterscommon/
      SAVE /errcommon/, /talk/
C======================== Include constans.inc =========================
      REAL       Mpip
      PARAMETER( Mpip = 139.57 )             
      REAL       Mpio
      PARAMETER( Mpio = 134.98 )             
      REAL       Mko
      PARAMETER( Mko  = 497.61 )             
      REAL       Momega
      PARAMETER( Momega  = 782.65 )          
      REAL       Cvel
      PARAMETER( Cvel = 29.9792458 )         
      REAL       EMCvel
      PARAMETER( EMCvel = 28.17 )            
      REAL       TauKs
      PARAMETER( TauKs = 0.08953 )   
      REAL       Mphi
      PARAMETER( Mphi  = 1019.460 )          
C-----------------------------------------------------------------------
C= Include /kloe/soft/off/offline/inc/development/tls/maxstructdim.cin =
      INTEGER     NeleCluMax
      PARAMETER ( NeleCluMax = 2000 )
      INTEGER     MaxNumClu
      PARAMETER ( MaxNumClu = 100 )
      INTEGER    MaxNumVtx
      PARAMETER (MaxNumVtx = 20)
      INTEGER    MaxNumTrkV
      PARAMETER (MaxNumTrkV = 30)
      INTEGER    MaxNumTrk
      PARAMETER (MaxNumTrk = 100)
      INTEGER    MaxNumDHSP
      PARAMETER (MaxNumDHSP = 1000)
      Integer    nMaxDC
      Parameter (nMaxDC=1500)
      INTEGER    NQihiMax
      PARAMETER (NQihiMax=1000)
      INTEGER    NQcalMax
      PARAMETER (NQcalMax= 32 )
      INTEGER    MaxNumFirstHit
      PARAMETER (MaxNumFirstHit = 300 )
      INTEGER    MaxNtrkGen
      PARAMETER (MaxNtrkGen =50)
      INTEGER    MaxNvtxGen
      PARAMETER (MaxNvTxGen =30)
      integer TriggerElements
      parameter (TriggerElements = 300 )
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
         INTEGER iv(MaxNumtrkv)
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
         INTEGER mother(MaxNvtxGen)
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
C====================== Include interfstruct.inc =======================
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
          REAL    Tcl(MaxNumClu)
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
          INTEGER iv(MaxNumtrkv)
          INTEGER trknumv(MaxNumtrkv)
          REAL    CurV(MaxNumtrkv)
          REAL    PhiV(MaxNumtrkv)
          REAL    CotV(MaxNumtrkv)
          REAL    PxTV(MaxNumtrkv)
          REAL    PyTV(MaxNumtrkv)
          REAL    PzTV(MaxNumtrkv)
          INTEGER nv
          INTEGER vtx(MaxNumVtx)
          REAL    xv(MaxNumVtx)
          REAL    yv(MaxNumVtx)
          REAL    zv(MaxNumVtx)
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
          INTEGER pidmc(MaxNtrkGen)
          INTEGER virmom(MaxNtrkGen)
          REAL    pxmc(MaxNtrkGen)
          REAL    pymc(MaxNtrkGen)
          REAL    pzmc(MaxNtrkGen)
          REAL    themc(MaxNtrkGen)
          REAL    phimc(MaxNtrkGen)
          INTEGER vtxmc(MaxNtrkGen)
          INTEGER nvtxmc
          INTEGER kinmom(MaxNvtxGen)
          INTEGER mother(MaxNvtxGen)
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
          REAL    KchMC(9)
          REAL    KneMC(9)
          REAL    KchRec(9)
          REAL    KneRec(9)
          REAL    KchBoost(9)
          REAL    ip(3)
          REAL    ipmc(3)
          INTEGER vtaken(3)
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
          REAL    cosTrk
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
          INTEGER omegatruth_top
          INTEGER clu_index(4)
          INTEGER chtrk_index(2)
          INTEGER chvtx_index
          REAL P4gam_w(4,4)
          REAL P4ch_w(2,4)
          REAL P4ch_tot_w(4)
          REAL om_dist(4)
      END TYPE
      TYPE ( interfstru) Interf
      COMMON / InterfCommon / Interf,Evt
C-----------------------------------------------------------------------
C
C Local declarations
C
      INTEGER id
C-----------------------------------------------------------------------
C in: id
C-----------------------------------------------------------------------
      Interf%CutId = id
      IF( truthmc) THEN
        counter(id) = counter(id) + 1
      ELSE
        counterBcg(id) = counterBcg(id) + 1
      ENDIF
      END
C=======================================================================
C=======================================================================
C=======================================================================
      SUBROUTINE statisticsErr(id)
        IMPLICIT NONE
C========================= Include klspm00.inc =========================
      INTEGER NCLMIN
      PARAMETER( NCLMIN = 4 )
      REAL ECLMIN
      PARAMETER( ECLMIN = 20. )
      REAL vtxdist
      PARAMETER( vtxdist = 50. )
      INTEGER testcounter
      PARAMETER( testcounter = 100 )
      INTEGER MaxNumFitPar
      PARAMETER( MaxNumFitPar = 38 )
      INTEGER MaxNumConstr
      PARAMETER( MaxNumConstr = 10 )
      INTEGER MaxNumComega
      PARAMETER( MaxNumComega = 8 )
      INTEGER cutN
      PARAMETER( cutN = 6 )
      INTEGER errN
      PARAMETER( errN = 17 )
      CHARACTER*40 messageCut(testcounter),messageErr(testcounter)
      LOGICAL    McFlag,EventOK,ErrFlag,truthmc,truthregg,truthomegaa,tr
     1  uthsemii,truththreee,truthelsee                                 
      INTEGER    nmctruth0,ntruth,nmctruth(7),ntruthregg,ntruthsemii,ntr
     1  uththreee,ntruthomegaa,ntruthelsee                              
      INTEGER    selected,analised
      INTEGER    counter(testcounter),counterBcg(testcounter)
      INTEGER    ErrFlagCount(testcounter),ErrFlagCountBcg(testcounter)
      LOGICAL    cutchar,cutneu,onlytruemc,makecuts
      COMMON / flagcommon / McFlag,EventOK,ErrFlag,truthmc
      COMMON / truthcommon / nmctruth0,ntruth,nmctruth
      COMMON / counterscommon / selected,analised,
     &                          counter,counterBcg
      COMMON / errcommon / ErrFlagCount,ErrFlagCountBcg,
     &                     messageCut,messageErr
      COMMON / talk / cutchar,cutneu,onlytruemc,makecuts
      SAVE /flagcommon/, /truthcommon/, /counterscommon/
      SAVE /errcommon/, /talk/
C======================== Include constans.inc =========================
      REAL       Mpip
      PARAMETER( Mpip = 139.57 )             
      REAL       Mpio
      PARAMETER( Mpio = 134.98 )             
      REAL       Mko
      PARAMETER( Mko  = 497.61 )             
      REAL       Momega
      PARAMETER( Momega  = 782.65 )          
      REAL       Cvel
      PARAMETER( Cvel = 29.9792458 )         
      REAL       EMCvel
      PARAMETER( EMCvel = 28.17 )            
      REAL       TauKs
      PARAMETER( TauKs = 0.08953 )   
      REAL       Mphi
      PARAMETER( Mphi  = 1019.460 )          
C-----------------------------------------------------------------------
C= Include /kloe/soft/off/offline/inc/development/tls/maxstructdim.cin =
      INTEGER     NeleCluMax
      PARAMETER ( NeleCluMax = 2000 )
      INTEGER     MaxNumClu
      PARAMETER ( MaxNumClu = 100 )
      INTEGER    MaxNumVtx
      PARAMETER (MaxNumVtx = 20)
      INTEGER    MaxNumTrkV
      PARAMETER (MaxNumTrkV = 30)
      INTEGER    MaxNumTrk
      PARAMETER (MaxNumTrk = 100)
      INTEGER    MaxNumDHSP
      PARAMETER (MaxNumDHSP = 1000)
      Integer    nMaxDC
      Parameter (nMaxDC=1500)
      INTEGER    NQihiMax
      PARAMETER (NQihiMax=1000)
      INTEGER    NQcalMax
      PARAMETER (NQcalMax= 32 )
      INTEGER    MaxNumFirstHit
      PARAMETER (MaxNumFirstHit = 300 )
      INTEGER    MaxNtrkGen
      PARAMETER (MaxNtrkGen =50)
      INTEGER    MaxNvtxGen
      PARAMETER (MaxNvTxGen =30)
      integer TriggerElements
      parameter (TriggerElements = 300 )
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
         INTEGER iv(MaxNumtrkv)
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
         INTEGER mother(MaxNvtxGen)
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
C====================== Include interfstruct.inc =======================
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
          REAL    Tcl(MaxNumClu)
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
          INTEGER iv(MaxNumtrkv)
          INTEGER trknumv(MaxNumtrkv)
          REAL    CurV(MaxNumtrkv)
          REAL    PhiV(MaxNumtrkv)
          REAL    CotV(MaxNumtrkv)
          REAL    PxTV(MaxNumtrkv)
          REAL    PyTV(MaxNumtrkv)
          REAL    PzTV(MaxNumtrkv)
          INTEGER nv
          INTEGER vtx(MaxNumVtx)
          REAL    xv(MaxNumVtx)
          REAL    yv(MaxNumVtx)
          REAL    zv(MaxNumVtx)
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
          INTEGER pidmc(MaxNtrkGen)
          INTEGER virmom(MaxNtrkGen)
          REAL    pxmc(MaxNtrkGen)
          REAL    pymc(MaxNtrkGen)
          REAL    pzmc(MaxNtrkGen)
          REAL    themc(MaxNtrkGen)
          REAL    phimc(MaxNtrkGen)
          INTEGER vtxmc(MaxNtrkGen)
          INTEGER nvtxmc
          INTEGER kinmom(MaxNvtxGen)
          INTEGER mother(MaxNvtxGen)
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
          REAL    KchMC(9)
          REAL    KneMC(9)
          REAL    KchRec(9)
          REAL    KneRec(9)
          REAL    KchBoost(9)
          REAL    ip(3)
          REAL    ipmc(3)
          INTEGER vtaken(3)
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
          REAL    cosTrk
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
          INTEGER omegatruth_top
          INTEGER clu_index(4)
          INTEGER chtrk_index(2)
          INTEGER chvtx_index
          REAL P4gam_w(4,4)
          REAL P4ch_w(2,4)
          REAL P4ch_tot_w(4)
          REAL om_dist(4)
      END TYPE
      TYPE ( interfstru) Interf
      COMMON / InterfCommon / Interf,Evt
C-----------------------------------------------------------------------
C
C Local declarations
C
      INTEGER id
C-----------------------------------------------------------------------
C in: id
C-----------------------------------------------------------------------
      Interf%ErrId=id
      IF( truthmc) THEN
        ErrFlagCount(id) = ErrFlagCount(id) + 1
      ELSE
        ErrFlagCountBcg(id) = ErrFlagCountBcg(id) + 1
      ENDIF
      END
C=======================================================================
C=======================================================================
C=======================================================================
      SUBROUTINE getsim(ChVtxId,TrkIdx,Trk1,Trk2,ChaVtx,CluIdx,NeuVtx,
     &                  PhiVtx,simok)
        IMPLICIT NONE
C========================= Include klspm00.inc =========================
      INTEGER NCLMIN
      PARAMETER( NCLMIN = 4 )
      REAL ECLMIN
      PARAMETER( ECLMIN = 20. )
      REAL vtxdist
      PARAMETER( vtxdist = 50. )
      INTEGER testcounter
      PARAMETER( testcounter = 100 )
      INTEGER MaxNumFitPar
      PARAMETER( MaxNumFitPar = 38 )
      INTEGER MaxNumConstr
      PARAMETER( MaxNumConstr = 10 )
      INTEGER MaxNumComega
      PARAMETER( MaxNumComega = 8 )
      INTEGER cutN
      PARAMETER( cutN = 6 )
      INTEGER errN
      PARAMETER( errN = 17 )
      CHARACTER*40 messageCut(testcounter),messageErr(testcounter)
      LOGICAL    McFlag,EventOK,ErrFlag,truthmc,truthregg,truthomegaa,tr
     1  uthsemii,truththreee,truthelsee                                 
      INTEGER    nmctruth0,ntruth,nmctruth(7),ntruthregg,ntruthsemii,ntr
     1  uththreee,ntruthomegaa,ntruthelsee                              
      INTEGER    selected,analised
      INTEGER    counter(testcounter),counterBcg(testcounter)
      INTEGER    ErrFlagCount(testcounter),ErrFlagCountBcg(testcounter)
      LOGICAL    cutchar,cutneu,onlytruemc,makecuts
      COMMON / flagcommon / McFlag,EventOK,ErrFlag,truthmc
      COMMON / truthcommon / nmctruth0,ntruth,nmctruth
      COMMON / counterscommon / selected,analised,
     &                          counter,counterBcg
      COMMON / errcommon / ErrFlagCount,ErrFlagCountBcg,
     &                     messageCut,messageErr
      COMMON / talk / cutchar,cutneu,onlytruemc,makecuts
      SAVE /flagcommon/, /truthcommon/, /counterscommon/
      SAVE /errcommon/, /talk/
C======================== Include constans.inc =========================
      REAL       Mpip
      PARAMETER( Mpip = 139.57 )             
      REAL       Mpio
      PARAMETER( Mpio = 134.98 )             
      REAL       Mko
      PARAMETER( Mko  = 497.61 )             
      REAL       Momega
      PARAMETER( Momega  = 782.65 )          
      REAL       Cvel
      PARAMETER( Cvel = 29.9792458 )         
      REAL       EMCvel
      PARAMETER( EMCvel = 28.17 )            
      REAL       TauKs
      PARAMETER( TauKs = 0.08953 )   
      REAL       Mphi
      PARAMETER( Mphi  = 1019.460 )          
C====== Include /kloe/soft/off/offline/inc/development/runtyp.cin ======
        INTEGER    EXKPHY               
        PARAMETER (EXKPHY       = 0)
        INTEGER    EXCOSR               
        PARAMETER (EXCOSR       = 1)
        INTEGER    EXOFSI               
        PARAMETER (EXOFSI       = 2)
        INTEGER    EXTEST               
        PARAMETER (EXTEST       = 3)
        INTEGER    RUKPHY               
        PARAMETER (RUKPHY       = 1)
        INTEGER    RUCOSR               
        PARAMETER (RUCOSR       = 2)
        INTEGER    RUCALI               
        PARAMETER (RUCALI       = 3)
C====== Include /kloe/soft/off/offline/inc/development/jobsta.cin ======
        COMMON /JOBSTI/ NRUN, NEV, NRUNA, NEVA, IBEGRN, IENDRN,
     1                  IUDBUG, NEVP, NEVPR, NEVPA, NREC, NRECF
        INTEGER         NRUN, NEV, NRUNA, NEVA, IBEGRN, IENDRN
        INTEGER         IUDBUG, NEVP, NEVPR, NEVPA, NREC, NRECF
        COMMON /JOBLRI/ EXPTYP, PARTNO, RUNTYP, LRECNO, NRCTHS,
     1                  CDFCID, IRTYPA, VERNUM, RELNUM, CRPRID,
     2                  CRYEAR,  CRMON,  CRDAY,
     3                  CRHOUR,  CRMIN,  CRSEC
        INTEGER         EXPTYP, PARTNO, RUNTYP, LRECNO, NRCTHS
        INTEGER         CDFCID, IRTYPA, VERNUM, RELNUM, CRPRID
        INTEGER         CRYEAR,  CRMON,  CRDAY
        INTEGER         CRHOUR,  CRMIN,  CRSEC
        COMMON /JOBEVC/ TRIGNO, TRGSUM(4),
     1                  STEVTI, SDONET, RPTMOD, RESERV, SOPTYP, SBUFNO,
     2                  SREADO(2), SRESIO(2), CREADO(2), TERSUM(2)
        INTEGER         TRIGNO,    TRGSUM
        INTEGER         STEVTI, SDONET, RPTMOD, RESERV, SOPTYP, SBUFNO
        INTEGER         SREADO,    SRESIO,    CREADO,    TERSUM
        COMMON /JOBSTR/ ROOTS, BMAGNT, BMAGFT, BMAGBT
        REAL            ROOTS, BMAGNT, BMAGFT, BMAGBT
C==== Include /kloe/soft/off/a_c/production/aix/library/anerror.cin ====
      INTEGER    ANSUCC
      PARAMETER (ANSUCC = '8198009'X)
      INTEGER    ANBEGN
      PARAMETER (ANBEGN = '8198019'X)
      INTEGER    ANPDUN
      PARAMETER (ANPDUN = '8198021'X)
      INTEGER    ANAST
      PARAMETER (ANAST   = '81980A3'X)
      INTEGER    ANILPS
      PARAMETER (ANILPS = '8198122'X)
      INTEGER    ANILMD
      PARAMETER (ANILMD = '819812A'X)
      INTEGER    ANLSTR
      PARAMETER (ANLSTR = '81980C0'X)
      INTEGER    ANNINI
      PARAMETER (ANNINI = '8198102'X)
      INTEGER    ANOVRF
      PARAMETER (ANOVRF = '819810A'X)
      INTEGER    ANNOME
      PARAMETER (ANNOME = '8198112'X)
      INTEGER    ANNOLU
      PARAMETER (ANNOLU = '819811A'X)
      INTEGER    ANNOEP
      PARAMETER (ANNOEP = '8198011'X)
      INTEGER    ANILGP
      PARAMETER (ANILGP = '8198132'X)
      INTEGER    ANILPA
      PARAMETER (ANILPA = '819813A'X)
      INTEGER    ANILQU
      PARAMETER (ANILQU = '8198142'X)
      INTEGER    ANILSY
      PARAMETER (ANILSY = '819814A'X)
      INTEGER    ANILLU
      PARAMETER (ANILLU = '8198152'X)
      INTEGER    ANILVB
      PARAMETER (ANILVB = '819815A'X)
      INTEGER    ANILLI
      PARAMETER (ANILLI = '8198162'X)
      INTEGER    ANDUMD
      PARAMETER (ANDUMD = '819816A'X)
      INTEGER    ANYBER
      PARAMETER (ANYBER = '8198172'X)
      INTEGER    ANYBFE
      PARAMETER (ANYBFE = '081981B2'X)
      INTEGER    ANEFOP
      PARAMETER (ANEFOP = '819817A'X)
      INTEGER    ANEOF
      PARAMETER (ANEOF = '8198083'X)
      INTEGER    ANIFNO
      PARAMETER (ANIFNO = '819808B'X)
      INTEGER    ANCLOS
      PARAMETER (ANCLOS = '8198093'X)
      INTEGER    ANOFNO
      PARAMETER (ANOFNO = '819809B'X)
      INTEGER    ANNOLR
      PARAMETER (ANNOLR = '8198182'X)
      INTEGER    ANILST
      PARAMETER (ANILST = '819818A'X)
      INTEGER    ANNOEV
      PARAMETER (ANNOEV = '8198192'X)
      INTEGER    ANILBL
      PARAMETER (ANILBL = '819819A'X)
      INTEGER    ANNEIF
      PARAMETER (ANNEIF = '81981A2'X)
      INTEGER    ANNOBK
      PARAMETER (ANNOBK = '81980C8'X)
      INTEGER    ANILBK
      PARAMETER (ANILBK = '81981AA'X)
C==== Include /kloe/soft/off/ybos/production/aix/library/errcod.cin ====
      INTEGER    YESUCC
      PARAMETER (YESUCC = '08508009'X)
      INTEGER    YEWLOK
      PARAMETER (YEWLOK = '08508083'X)
      INTEGER    YEWRNS
      PARAMETER (YEWRNS = '0850808B'X)
      INTEGER    YEINDX
      PARAMETER (YEINDX = '08508093'X)
      INTEGER    YEEOF
      PARAMETER (YEEOF  = '08508100'X)
      INTEGER    YEPRBD
      PARAMETER (YEPRBD = '08508108'X)
      INTEGER    YERSTF
      PARAMETER (YERSTF = '08508110'X)
      INTEGER    YEOVRF
      PARAMETER (YEOVRF = '08508182'X)
      INTEGER    YEBKNG
      PARAMETER (YEBKNG = '0850818A'X)
      INTEGER    YEBKTS
      PARAMETER (YEBKTS = '08508192'X)
      INTEGER    YEBKOV
      PARAMETER (YEBKOV = '0850819A'X)
      INTEGER    YEBKSP
      PARAMETER (YEBKSP = '085081A2'X)
      INTEGER    YENOBK
      PARAMETER (YENOBK = '085081AA'X)
      INTEGER    YENONR
      PARAMETER (YENONR = '085081B2'X)
      INTEGER    YEMANY
      PARAMETER (YEMANY = '085081BA'X)
      INTEGER    YEDUPB
      PARAMETER (YEDUPB = '085081C2'X)
      INTEGER    YEBSIL
      PARAMETER (YEBSIL = '085081CA'X)
      INTEGER    YEILOP
      PARAMETER (YEILOP = '085081D2'X)
      INTEGER    YEAROF
      PARAMETER (YEAROF = '085081DA'X)
      INTEGER    YETRUN
      PARAMETER (YETRUN = '085081E2'X)
      INTEGER    YECORR
      PARAMETER (YECORR = '085081EA'X)
      INTEGER    YEWRGI
      PARAMETER (YEWRGI = '085081F2'X)
      INTEGER    YELWIU
      PARAMETER (YELWIU = '085081FA'X)
      INTEGER    YEILUN
      PARAMETER (YEILUN = '08508202'X)
      INTEGER    YEIOUS
      PARAMETER (YEIOUS = '0850820A'X)
      INTEGER    YELROF
      PARAMETER (YELROF = '08508212'X)
      INTEGER    YELRLN
      PARAMETER (YELRLN = '0850821A'X)
      INTEGER    YEPRXL
      PARAMETER (YEPRXL = '08508222'X)
      INTEGER    YEBKOS
      PARAMETER (YEBKOS = '0850822A'X)
      INTEGER    YEBKXD
      PARAMETER (YEBKXD = '08508232'X)
      INTEGER    YELRZL
      PARAMETER (YELRZL = '0850823A'X)
      INTEGER    YELRWT
      PARAMETER (YELRWT = '08508242'X)
      INTEGER    YEDKOP
      PARAMETER (YEDKOP = '0850824A'X)
      INTEGER    YEDKWT
      PARAMETER (YEDKWT = '08508252'X)
      INTEGER    YEDKRD
      PARAMETER (YEDKRD = '0850825A'X)
      INTEGER    YEILUS
      PARAMETER (YEILUS = '08508282'X)
      INTEGER    YERSUS
      PARAMETER (YERSUS = '0850828A'X)
      INTEGER    YENRUS
      PARAMETER (YENRUS = '08508292'X)
      INTEGER    YEILLA
      PARAMETER (YEILLA = '08508302'X)
      INTEGER    YEDUPA
      PARAMETER (YEDUPA = '0850830A'X)
      INTEGER    YEBTIL
      PARAMETER (YEBTIL = '08508382'X)
      INTEGER    YEGRIL
      PARAMETER (YEGRIL = '0850838A'X)
      INTEGER    YECHKF
      PARAMETER (YECHKF = '08508392'X)
      INTEGER    YESTOF
      PARAMETER (YESTOF = '0850839A'X)
      INTEGER    YENOTA
      PARAMETER (YENOTA = '085083A2'X)
      INTEGER    YEBNKP
      PARAMETER (YEBNKP = '085083AA'X)
      INTEGER    YEBDCH
      PARAMETER (YEBDCH = '085083C2'X)
      INTEGER    YEBDBK
      PARAMETER (YEBDBK = '085083CA'X)
      INTEGER    YEBDLU
      PARAMETER (YEBDLU = '085083D2'X)
      INTEGER    YEBDBS
      PARAMETER (YEBDBS = '085083DA'X)
C====== Include /kloe/soft/off/offline/inc/development/erlevl.cin ======
        INTEGER    ERMINI
        PARAMETER (ERMINI = 0)
        INTEGER    ERMAXI
        PARAMETER (ERMAXI = 7)
        INTEGER    ERSUCC
        PARAMETER (ERSUCC = 0)
        INTEGER    ERINFO
        PARAMETER (ERINFO = 1)
        INTEGER    ERWARN
        PARAMETER (ERWARN = 2)
        INTEGER    EREROR
        PARAMETER (EREROR = 3)
        INTEGER    ERSEVR
        PARAMETER (ERSEVR = 4)
        INTEGER         ELSUCC
        PARAMETER       (ELSUCC = 10)
        INTEGER         ELINFO
        PARAMETER       (ELINFO = 11)
        INTEGER         ELWARN
        PARAMETER       (ELWARN = 12)
        INTEGER         ELWARN2
        PARAMETER       (ELWARN2 = 13)
        INTEGER         ELEROR
        PARAMETER       (ELEROR = 14)
        INTEGER         ELEROR2
        PARAMETER       (ELEROR2 = 15)
        INTEGER         ELNEXT
        PARAMETER       (ELNEXT = 16)
        INTEGER         ELNEXT2
        PARAMETER       (ELNEXT2 = 17)
        INTEGER         ELSEVR
        PARAMETER       (ELSEVR = 18)
        INTEGER         ELABOR
        PARAMETER       (ELABOR = 19)
C====== Include /kloe/soft/off/offline/inc/development/bpoybs.cin ======
      INTEGER RUNGID
      PARAMETER (RUNGID = 0)
      INTEGER RUNVER
      PARAMETER (RUNVER = 1)
      INTEGER RUNRND
      PARAMETER (RUNRND = 3)
      INTEGER RUNEVS
      PARAMETER (RUNEVS = 5)
      INTEGER PARNUM
      PARAMETER (PARNUM = 0)
      INTEGER PARNAM
      PARAMETER (PARNAM = 1)
      INTEGER PARTRT
      PARAMETER (PARTRT = 6)
      INTEGER PARMAS
      PARAMETER (PARMAS = 7)
      INTEGER PARCHR
      PARAMETER (PARCHR = 8)
      INTEGER PARLTI
      PARAMETER (PARLTI = 9)
      INTEGER PARBRT
      PARAMETER (PARBRT = 10)
      INTEGER PARDMO
      PARAMETER (PARDMO = 16)
      INTEGER PAUSER
      PARAMETER (PAUSER = 22)
      INTEGER MATNUM
      PARAMETER (MATNUM = 0)
      INTEGER MATNAM
      PARAMETER (MATNAM = 1)
      INTEGER MATATN
      PARAMETER (MATATN = 6)
      INTEGER MATMAS
      PARAMETER (MATMAS = 7)
      INTEGER MATDEN
      PARAMETER (MATDEN = 8)
      INTEGER MATRLE
      PARAMETER (MATRLE = 9)
      INTEGER MATABL
      PARAMETER (MATABL = 10)
      INTEGER MAUSER
      PARAMETER (MAUSER = 11)
      INTEGER TMENUM
      PARAMETER (TMENUM = 0)
      INTEGER TMENAM
      PARAMETER (TMENAM = 1)
      INTEGER TMEMAT
      PARAMETER (TMEMAT = 6)
      INTEGER TMESTF
      PARAMETER (TMESTF = 7)
      INTEGER TMEMFI
      PARAMETER (TMEMFI = 8)
      INTEGER TMEMFM
      PARAMETER (TMEMFM = 9)
      INTEGER TMEMAA
      PARAMETER (TMEMAA = 10)
      INTEGER TMEMAS
      PARAMETER (TMEMAS = 11)
      INTEGER TMEMAE
      PARAMETER (TMEMAE = 12)
      INTEGER TMEEPS
      PARAMETER (TMEEPS = 13)
      INTEGER TMEMIS
      PARAMETER (TMEMIS = 14)
      INTEGER TMUSER
      PARAMETER (TMUSER = 15)
      INTEGER BRIHVB
      PARAMETER (BRIHVB = 0)
      INTEGER BRIHVS
      PARAMETER (BRIHVS = 1)
      INTEGER BRIVD1
      PARAMETER (BRIVD1 = 2)
      INTEGER BRIPRS
      PARAMETER (BRIPRS = 10)
      INTEGER BRITMP
      PARAMETER (BRITMP = 12)
      INTEGER BRIISO
      PARAMETER (BRIISO = 14)
      INTEGER BRIMGC
      PARAMETER (BRIMGC = 17)
      INTEGER BRIVD2
      PARAMETER (BRIVD2 = 18)
      INTEGER BRIVD3
      PARAMETER (BRIVD3 = 32)
      INTEGER BRINCR
      PARAMETER (BRINCR = 64)
      INTEGER BRINDB
      PARAMETER (BRINDB = 65)
      INTEGER BRINCB
      PARAMETER (BRINCB = 66)
      INTEGER BRINDH
      PARAMETER (BRINDH = 67)
      INTEGER HEARUN
      PARAMETER (HEARUN = 0)
      INTEGER HEAEVT
      PARAMETER (HEAEVT = 1)
      INTEGER HEARND
      PARAMETER (HEARND = 2)
      INTEGER HEAQMI
      PARAMETER (HEAQMI = 4)
      INTEGER HEAKIN
      PARAMETER (HEAKIN = 5)
      INTEGER HEAVER
      PARAMETER (HEAVER = 6)
      INTEGER HEAUSE
      PARAMETER (HEAUSE = 7)
      INTEGER KINPPX
      PARAMETER (KINPPX = 0)
      INTEGER KINPPY
      PARAMETER (KINPPY = 1)
      INTEGER KINPPZ
      PARAMETER (KINPPZ = 2)
      INTEGER KINPEN
      PARAMETER (KINPEN = 3)
      INTEGER KINPTY
      PARAMETER (KINPTY = 4)
      INTEGER KINVXO
      PARAMETER (KINVXO = 5)
      INTEGER KINNVX
      PARAMETER (KINNVX = 6)
      INTEGER KINVLS
      PARAMETER (KINVLS = 7)
      INTEGER VERXCO
      PARAMETER (VERXCO = 0)
      INTEGER VERYCO
      PARAMETER (VERYCO = 1)
      INTEGER VERZCO
      PARAMETER (VERZCO = 2)
      INTEGER VERTOF
      PARAMETER (VERTOF = 3)
      INTEGER VERTNO
      PARAMETER (VERTNO = 4)
      INTEGER VERTAR
      PARAMETER (VERTAR = 5)
      INTEGER VERNTR
      PARAMETER (VERNTR = 6)
      INTEGER VERTLS
      PARAMETER (VERTLS = 7)
      INTEGER DHIHDS
      PARAMETER (DHIHDS = 0)
      INTEGER DHIVRN
      PARAMETER (DHIVRN = 1)
      INTEGER DHINRO
      PARAMETER (DHINRO = 2)
      INTEGER DHINCO
      PARAMETER (DHINCO = 3)
      INTEGER DHIPTY
      PARAMETER (DHIPTY = 0)
      INTEGER DHITRA
      PARAMETER (DHITRA = 1)
      INTEGER DHIADR
      PARAMETER (DHIADR = 2)
      INTEGER DHIXCO
      PARAMETER (DHIXCO = 3)
      INTEGER DHIYCO
      PARAMETER (DHIYCO = 4)
      INTEGER DHIZCO
      PARAMETER (DHIZCO = 5)
      INTEGER DHIPPX
      PARAMETER (DHIPPX = 6)
      INTEGER DHIPPY
      PARAMETER (DHIPPY = 7)
      INTEGER DHIPPZ
      PARAMETER (DHIPPZ = 8)
      INTEGER DHITOF
      PARAMETER (DHITOF = 9)
      INTEGER DHIELO
      PARAMETER (DHIELO =10)
      INTEGER DHILEN
      PARAMETER (DHILEN =11)
      INTEGER DHITIM
      PARAMETER (DHITIM =12)
      INTEGER DHIDIS
      PARAMETER (DHIDIS =13)
      INTEGER DHIFLA
      PARAMETER (DHIFLA =14)
      INTEGER VHIHDS
      PARAMETER (VHIHDS = 0)
      INTEGER VHIVRN
      PARAMETER (VHIVRN = 1)
      INTEGER VHINRO
      PARAMETER (VHINRO = 2)
      INTEGER VHINCO
      PARAMETER (VHINCO = 3)
      INTEGER VHIPTY
      PARAMETER (VHIPTY = 0)
      INTEGER VHITRA
      PARAMETER (VHITRA = 1)
      INTEGER VHIADR
      PARAMETER (VHIADR = 2)
      INTEGER VHIXEN
      PARAMETER (VHIXEN = 3)
      INTEGER VHIYEN
      PARAMETER (VHIYEN = 4)
      INTEGER VHIZEN
      PARAMETER (VHIZEN = 5)
      INTEGER VHIXEX
      PARAMETER (VHIXEX = 6)
      INTEGER VHIYEX
      PARAMETER (VHIYEX = 7)
      INTEGER VHIZEX
      PARAMETER (VHIZEX = 8)
      INTEGER VHIPMO
      PARAMETER (VHIPMO = 9)
      INTEGER VHITOF
      PARAMETER (VHITOF =10)
      INTEGER VHIELO
      PARAMETER (VHIELO =11)
      INTEGER VHITXC
      PARAMETER (VHITXC =12)
      INTEGER VHITYC
      PARAMETER (VHITYC =13)
      INTEGER CHIHDS
      PARAMETER (CHIHDS = 0)
      INTEGER CHIVRN
      PARAMETER (CHIVRN = 1)
      INTEGER CHINRO
      PARAMETER (CHINRO = 2)
      INTEGER CHINCO
      PARAMETER (CHINCO = 3)
      INTEGER CHIPTY
      PARAMETER (CHIPTY = 0)
      INTEGER CHITRA
      PARAMETER (CHITRA = 1)
      INTEGER CHIADR
      PARAMETER (CHIADR = 2)
      INTEGER CHIXCO
      PARAMETER (CHIXCO = 3)
      INTEGER CHIYCO
      PARAMETER (CHIYCO = 4)
      INTEGER CHIZCO
      PARAMETER (CHIZCO = 5)
      INTEGER CHITOF
      PARAMETER (CHITOF = 6)
      INTEGER CHIELO
      PARAMETER (CHIELO = 7)
      INTEGER CHILEN
      PARAMETER (CHILEN = 8)
      INTEGER CFHHDS
      PARAMETER (CFHHDS = 0)
      INTEGER CFHVRN
      PARAMETER (CFHVRN = 1)
      INTEGER CFHNRO
      PARAMETER (CFHNRO = 2)
      INTEGER CFHNCO
      PARAMETER (CFHNCO = 3)
      INTEGER CFHPTY
      PARAMETER (CFHPTY = 0)
      INTEGER CFHTRA
      PARAMETER (CFHTRA = 1)
      INTEGER CFHADR
      PARAMETER (CFHADR = 2)
      INTEGER CFHXCO
      PARAMETER (CFHXCO = 3)
      INTEGER CFHYCO
      PARAMETER (CFHYCO = 4)
      INTEGER CFHZCO
      PARAMETER (CFHZCO = 5)
      INTEGER CFHPPX
      PARAMETER (CFHPPX = 6)
      INTEGER CFHPPY
      PARAMETER (CFHPPY = 7)
      INTEGER CFHPPZ
      PARAMETER (CFHPPZ = 8)
      INTEGER CFHTOF
      PARAMETER (CFHTOF = 9)
      INTEGER CFHLEN
      PARAMETER (CFHLEN = 10)
      INTEGER AHIHDS
      PARAMETER (AHIHDS = 0)
      INTEGER AHIVRN
      PARAMETER (AHIVRN = 1)
      INTEGER AHINRO
      PARAMETER (AHINRO = 2)
      INTEGER AHINCO
      PARAMETER (AHINCO = 3)
      INTEGER AHIPTY
      PARAMETER (AHIPTY = 0)
      INTEGER AHITRA
      PARAMETER (AHITRA = 1)
      INTEGER AHIADR
      PARAMETER (AHIADR = 2)
      INTEGER AHIXCO
      PARAMETER (AHIXCO = 3)
      INTEGER AHIYCO
      PARAMETER (AHIYCO = 4)
      INTEGER AHIZCO
      PARAMETER (AHIZCO = 5)
      INTEGER AHIPPX
      PARAMETER (AHIPPX = 6)
      INTEGER AHIPPY
      PARAMETER (AHIPPY = 7)
      INTEGER AHIPPZ
      PARAMETER (AHIPPZ = 8)
      INTEGER AHITOF
      PARAMETER (AHITOF = 9)
      INTEGER AHIELO
      PARAMETER (AHIELO =10)
      INTEGER AHILEN
      PARAMETER (AHILEN =11)
      INTEGER QIHHDS
      PARAMETER (QIHHDS = 0)
      INTEGER QIHVRN
      PARAMETER (QIHVRN = 1)
      INTEGER QIHNRO
      PARAMETER (QIHNRO = 2)
      INTEGER QIHNCO
      PARAMETER (QIHNCO = 3)
      INTEGER QIHPTY
      PARAMETER (QIHPTY = 0)
      INTEGER QIHTRA
      PARAMETER (QIHTRA = 1)
      INTEGER QIHADR
      PARAMETER (QIHADR = 2)
      INTEGER QIHXCO
      PARAMETER (QIHXCO = 3)
      INTEGER QIHYCO
      PARAMETER (QIHYCO = 4)
      INTEGER QIHZCO
      PARAMETER (QIHZCO = 5)
      INTEGER QIHPPX
      PARAMETER (QIHPPX = 6)
      INTEGER QIHPPY
      PARAMETER (QIHPPY = 7)
      INTEGER QIHPPZ
      PARAMETER (QIHPPZ = 8)
      INTEGER QIHTOF
      PARAMETER (QIHTOF = 9)
      INTEGER QIHELO
      PARAMETER (QIHELO =10)
      INTEGER QIHLEN
      PARAMETER (QIHLEN =11)
      INTEGER QHIHDS
      PARAMETER (QHIHDS = 0)
      INTEGER QHIVRN
      PARAMETER (QHIVRN = 1)
      INTEGER QHINRO
      PARAMETER (QHINRO = 2)
      INTEGER QHINCO
      PARAMETER (QHINCO = 3)
      INTEGER QHIPTY
      PARAMETER (QHIPTY = 0)
      INTEGER QHITRA
      PARAMETER (QHITRA = 1)
      INTEGER QHIADR
      PARAMETER (QHIADR = 2)
      INTEGER QHIXCO
      PARAMETER (QHIXCO = 3)
      INTEGER QHIYCO
      PARAMETER (QHIYCO = 4)
      INTEGER QHIZCO
      PARAMETER (QHIZCO = 5)
      INTEGER QHIPPX
      PARAMETER (QHIPPX = 6)
      INTEGER QHIPPY
      PARAMETER (QHIPPY = 7)
      INTEGER QHIPPZ
      PARAMETER (QHIPPZ = 8)
      INTEGER QHITOF
      PARAMETER (QHITOF = 9)
      INTEGER QHIELO
      PARAMETER (QHIELO =10)
      INTEGER QHILEN
      PARAMETER (QHILEN =11)
      INTEGER ITHHDS
      PARAMETER (ITHHDS = 0)
      INTEGER ITHVRN
      PARAMETER (ITHVRN = 1)
      INTEGER ITHNRO
      PARAMETER (ITHNRO = 2)
      INTEGER ITHNCO
      PARAMETER (ITHNCO = 3)
      INTEGER ITHINX
      PARAMETER (ITHINX = 0)
      INTEGER ITHINY
      PARAMETER (ITHINY = 1)
      INTEGER ITHINZ
      PARAMETER (ITHINZ = 2)
      INTEGER ITHOUX
      PARAMETER (ITHOUX = 3)
      INTEGER ITHOUY
      PARAMETER (ITHOUY = 4)
      INTEGER ITHOUZ
      PARAMETER (ITHOUZ = 5)
      INTEGER ITHPIX
      PARAMETER (ITHPIX = 6)
      INTEGER ITHPIY
      PARAMETER (ITHPIY = 7)
      INTEGER ITHPIZ
      PARAMETER (ITHPIZ = 8)
      INTEGER ITHTRL
      PARAMETER (ITHTRL = 9)
      INTEGER ITHDER
      PARAMETER (ITHDER = 10)
      INTEGER ITHTOF
      PARAMETER (ITHTOF = 11)
      INTEGER ITHPLA
      PARAMETER (ITHPLA = 12)
      INTEGER ITHTYP
      PARAMETER (ITHTYP = 13)
      INTEGER ITHTRA
      PARAMETER (ITHTRA = 14)
      INTEGER QTHHDS
      PARAMETER (QTHHDS = 0)
      INTEGER QTHVRN
      PARAMETER (QTHVRN = 1)
      INTEGER QTHNRO
      PARAMETER (QTHNRO = 2)
      INTEGER QTHNCO
      PARAMETER (QTHNCO = 3)
      INTEGER QTHPTY
      PARAMETER (QTHPTY = 0)
      INTEGER QTHTRA
      PARAMETER (QTHTRA = 1)
      INTEGER QTHADR
      PARAMETER (QTHADR = 2)
      INTEGER QTHXCO
      PARAMETER (QTHXCO = 3)
      INTEGER QTHYCO
      PARAMETER (QTHYCO = 4)
      INTEGER QTHZCO
      PARAMETER (QTHZCO = 5)
      INTEGER QTHTOF
      PARAMETER (QTHTOF = 6)
      INTEGER QTHELO
      PARAMETER (QTHELO = 7)
      INTEGER QTHLEN
      PARAMETER (QTHLEN = 8)
      INTEGER CTHHDS
      PARAMETER (CTHHDS = 0)
      INTEGER CTHVRN
      PARAMETER (CTHVRN = 1)
      INTEGER CTHNRO
      PARAMETER (CTHNRO = 2)
      INTEGER CTHNCO
      PARAMETER (CTHNCO = 3)
      INTEGER CTHPTY
      PARAMETER (CTHPTY = 0)
      INTEGER CTHTRA
      PARAMETER (CTHTRA = 1)
      INTEGER CTHADR
      PARAMETER (CTHADR = 2)
      INTEGER CTHXCO
      PARAMETER (CTHXCO = 3)
      INTEGER CTHYCO
      PARAMETER (CTHYCO = 4)
      INTEGER CTHZCO
      PARAMETER (CTHZCO = 5)
      INTEGER CTHTOF
      PARAMETER (CTHTOF = 6)
      INTEGER CTHELO
      PARAMETER (CTHELO = 7)
      INTEGER CTHLEN
      PARAMETER (CTHLEN = 8)
      INTEGER LEHHDS
      PARAMETER (LEHHDS = 0)
      INTEGER LEHVRN
      PARAMETER (LEHVRN = 1)
      INTEGER LEHNRO
      PARAMETER (LEHNRO = 2)
      INTEGER LEHNCO
      PARAMETER (LEHNCO = 3)
      INTEGER LEHPTY
      PARAMETER (LEHPTY = 0)
      INTEGER LEHTRA
      PARAMETER (LEHTRA = 1)
      INTEGER LEHADR
      PARAMETER (LEHADR = 2)
      INTEGER LEHXCO
      PARAMETER (LEHXCO = 3)
      INTEGER LEHYCO
      PARAMETER (LEHYCO = 4)
      INTEGER LEHZCO
      PARAMETER (LEHZCO = 5)
      INTEGER LEHXOU
      PARAMETER (LEHXOU = 6)
      INTEGER LEHYOU
      PARAMETER (LEHYOU = 7)
      INTEGER LEHZOU
      PARAMETER (LEHZOU = 8)
      INTEGER LEHTOF
      PARAMETER (LEHTOF = 9)
      INTEGER LEHELO
      PARAMETER (LEHELO = 10)
      INTEGER LEHLEN
      PARAMETER (LEHLEN = 11)
      INTEGER DTCHDS
      PARAMETER (DTCHDS = 0)
      INTEGER DTCVRN
      PARAMETER (DTCVRN = 1)
      INTEGER DTCNRO
      PARAMETER (DTCNRO = 2)
      INTEGER DTCNCO
      PARAMETER (DTCNCO = 3)
      INTEGER DTCCAL
      PARAMETER (DTCCAL = 4)
      INTEGER DTCADR
      PARAMETER (DTCADR = 0)
      INTEGER DTCTIM
      PARAMETER (DTCTIM = 1)
      INTEGER DTKHDS
      PARAMETER (DTKHDS = 0)
      INTEGER DTKVRN
      PARAMETER (DTKVRN = 1)
      INTEGER DTKNRO
      PARAMETER (DTKNRO = 2)
      INTEGER DTKNTR
      PARAMETER (DTKNTR = 0)
      INTEGER DTKADR
      PARAMETER (DTKADR = 1)
      INTEGER DTHHDS
      PARAMETER (DTHHDS = 0)
      INTEGER DTHVRN
      PARAMETER (DTHVRN = 1)
      INTEGER DTHNRO
      PARAMETER (DTHNRO = 2)
      INTEGER DTHNHI
      PARAMETER (DTHNHI = 0)
      INTEGER DHNHDS
      PARAMETER (DHNHDS = 0)
      INTEGER DHNVRN
      PARAMETER (DHNVRN = 1)
      INTEGER DHNNRO
      PARAMETER (DHNNRO = 2)
      INTEGER DHNNHI
      PARAMETER (DHNNHI = 0)
      INTEGER CELHDS
      PARAMETER (CELHDS = 0)
      INTEGER CELVRN
      PARAMETER (CELVRN = 1)
      INTEGER CELCAL
      PARAMETER (CELCAL = 2)
      INTEGER CELNRO
      PARAMETER (CELNRO = 3)
      INTEGER CELNCO
      PARAMETER (CELNCO = 4)
      INTEGER CELADR
      PARAMETER (CELADR = 0)
      INTEGER CELEA
      PARAMETER (CELEA = 1)
      INTEGER CELEB
      PARAMETER (CELEB = 2)
      INTEGER CELTA
      PARAMETER (CELTA = 3)
      INTEGER CELTB
      PARAMETER (CELTB = 4)
      INTEGER CEKHDS
      PARAMETER (CEKHDS = 0)
      INTEGER CEKVRN
      PARAMETER (CEKVRN = 1)
      INTEGER CEKNRO
      PARAMETER (CEKNRO = 2)
      INTEGER CEKNTR
      PARAMETER (CEKNTR = 0)
      INTEGER CEKADR
      PARAMETER (CEKADR = 1)
      INTEGER CHNHDS
      PARAMETER (CHNHDS = 0)
      INTEGER CHNVRN
      PARAMETER (CHNVRN = 1)
      INTEGER CHNNRO
      PARAMETER (CHNNRO = 2)
      INTEGER CHNNHI
      PARAMETER (CHNNHI = 0)
      INTEGER QCAHDS
      PARAMETER (QCAHDS = 0)
      INTEGER QCAVRN
      PARAMETER (QCAVRN = 1)
      INTEGER QCACAL
      PARAMETER (QCACAL = 2)
      INTEGER QCANRO
      PARAMETER (QCANRO = 3)
      INTEGER QCANCO
      PARAMETER (QCANCO = 4)
      INTEGER QCADR
      PARAMETER (QCADR = 0)
      INTEGER QCAE
      PARAMETER (QCAE = 1)
      INTEGER QCAT
      PARAMETER (QCAT = 2)
      INTEGER QCKHDS
      PARAMETER (QCKHDS = 0)
      INTEGER QCKVRN
      PARAMETER (QCKVRN = 1)
      INTEGER QCKNRO
      PARAMETER (QCKNRO = 2)
      INTEGER QCKNTR
      PARAMETER (QCKNTR = 0)
      INTEGER QCKADR
      PARAMETER (QCKADR = 1)
      INTEGER ITCHDS
      PARAMETER (ITCHDS = 0)
      INTEGER ITCVRN
      PARAMETER (ITCVRN = 1)
      INTEGER ITCNRO
      PARAMETER (ITCNRO = 2)
      INTEGER ITCNCO
      PARAMETER (ITCNCO = 3)
      INTEGER ITCCAL
      PARAMETER (ITCCAL = 4)
      INTEGER ITCADR
      PARAMETER (ITCADR = 0)
      INTEGER ITCTIM
      PARAMETER (ITCTIM = 1)
      INTEGER ITKHDS
      PARAMETER (ITKHDS = 0)
      INTEGER ITKVRN
      PARAMETER (ITKVRN = 1)
      INTEGER ITKNRO
      PARAMETER (ITKNRO = 2)
      INTEGER ITKNTR
      PARAMETER (ITKNTR = 0)
      INTEGER ITKADR
      PARAMETER (ITKADR = 1)
      INTEGER ITLHDS
      PARAMETER (ITLHDS = 0)
      INTEGER ITLVRN
      PARAMETER (ITLVRN = 1)
      INTEGER ITHLRO
      PARAMETER (ITHLRO = 2)
      INTEGER ITHLHI
      PARAMETER (ITHLHI = 0)
      INTEGER IHNHDS
      PARAMETER (IHNHDS = 0)
      INTEGER IHNVRN
      PARAMETER (IHNVRN = 1)
      INTEGER IHNNRO
      PARAMETER (IHNNRO = 2)
      INTEGER IHNNHI
      PARAMETER (IHNNHI = 0)
      INTEGER QCTHDS
      PARAMETER (QCTHDS = 0)
      INTEGER QCTVRN
      PARAMETER (QCTVRN = 1)
      INTEGER QCTCAL
      PARAMETER (QCTCAL = 2)
      INTEGER QCTNRO
      PARAMETER (QCTNRO = 3)
      INTEGER QCTADR
      PARAMETER (QCTADR = 1)
      INTEGER QCTNTD
      PARAMETER (QCTNTD = 0)
      INTEGER QCTTIM
      PARAMETER (QCTTIM = 0)
      INTEGER QCTWID
      PARAMETER (QCTWID = 1)
      INTEGER QTKHDS
      PARAMETER (QTKHDS = 0)
      INTEGER QTKVRN
      PARAMETER (QTKVRN = 1)
      INTEGER QTKNRO
      PARAMETER (QTKNRO = 2)
      INTEGER QTKNTR
      PARAMETER (QTKNTR = 0)
      INTEGER QTKADR
      PARAMETER (QTKADR = 1)
      INTEGER QTLHDS
      PARAMETER (QTLHDS = 0)
      INTEGER QTLVRN
      PARAMETER (QTLVRN = 1)
      INTEGER QTLNRO
      PARAMETER (QTLNRO = 2)
      INTEGER QTLNHI
      PARAMETER (QTLNHI = 0)
      INTEGER THNHDS
      PARAMETER (THNHDS = 0)
      INTEGER THNVRN
      PARAMETER (THNVRN = 1)
      INTEGER THNNRO
      PARAMETER (THNNRO = 2)
      INTEGER THNNHI
      PARAMETER (THNNHI = 0)
      INTEGER CCTHDS
      PARAMETER (CCTHDS = 0)
      INTEGER CCTVRN
      PARAMETER (CCTVRN = 1)
      INTEGER CCTCAL
      PARAMETER (CCTCAL = 2)
      INTEGER CCTNRO
      PARAMETER (CCTNRO = 3)
      INTEGER CCTNCO
      PARAMETER (CCTNCO = 4)
      INTEGER CCTADR
      PARAMETER (CCTADR = 0)
      INTEGER CCTTIM
      PARAMETER (CCTTIM = 2)
      INTEGER CCTELO
      PARAMETER (CCTELO = 1)
      INTEGER CCKHDS
      PARAMETER (CCKHDS = 0)
      INTEGER CCKVRN
      PARAMETER (CCKVRN = 1)
      INTEGER CCKNRO
      PARAMETER (CCKNRO = 2)
      INTEGER CCKNTR
      PARAMETER (CCKNTR = 0)
      INTEGER CCKADR
      PARAMETER (CCKADR = 1)
      INTEGER CCLHDS
      PARAMETER (CCLHDS = 0)
      INTEGER CCLVRN
      PARAMETER (CCLVRN = 1)
      INTEGER CCLNRO
      PARAMETER (CCLNRO = 2)
      INTEGER CCLNHI
      PARAMETER (CCLNHI = 0)
      INTEGER YHNHDS
      PARAMETER (YHNHDS = 0)
      INTEGER YHNVRN
      PARAMETER (YHNVRN = 1)
      INTEGER YHNNRO
      PARAMETER (YHNNRO = 2)
      INTEGER YHNNHI
      PARAMETER (YHNNHI = 0)
      INTEGER LETHDS
      PARAMETER (LETHDS = 0)
      INTEGER LETVRN
      PARAMETER (LETVRN = 1)
      INTEGER LETCAL
      PARAMETER (LETCAL = 2)
      INTEGER LETNRO
      PARAMETER (LETNRO = 3)
      INTEGER LETNCO
      PARAMETER (LETNCO = 4)
      INTEGER LETADR
      PARAMETER (LETADR = 0)
      INTEGER LETTIM
      PARAMETER (LETTIM = 2)
      INTEGER LETELO
      PARAMETER (LETELO = 1)
      INTEGER LEKHDS
      PARAMETER (LEKHDS = 0)
      INTEGER LEKVRN
      PARAMETER (LEKVRN = 1)
      INTEGER LEKNRO
      PARAMETER (LEKNRO = 2)
      INTEGER LEKNTR
      PARAMETER (LEKNTR = 0)
      INTEGER LEKADR
      PARAMETER (LEKADR = 1)
      INTEGER LELHDS
      PARAMETER (LELHDS = 0)
      INTEGER LELVRN
      PARAMETER (LELVRN = 1)
      INTEGER LELNRO
      PARAMETER (LELNRO = 2)
      INTEGER LELNHI
      PARAMETER (LELNHI = 0)
      INTEGER LENHDS
      PARAMETER (LENHDS = 0)
      INTEGER LENVRN
      PARAMETER (LENVRN = 1)
      INTEGER LENNRO
      PARAMETER (LENNRO = 2)
      INTEGER LENNHI
      PARAMETER (LENNHI = 0)
      INTEGER DTDHDS
      PARAMETER (DTDHDS = 0)
      INTEGER DTDVRN
      PARAMETER (DTDVRN = 1)
      INTEGER DTDNRO
      PARAMETER (DTDNRO = 2)
      INTEGER DTDNCO
      PARAMETER (DTDNCO = 3)
      INTEGER DTDZES
      PARAMETER (DTDZES = 4)
      INTEGER DTDADR
      PARAMETER (DTDADR = 0)
      INTEGER DADHDS
      PARAMETER (DADHDS = 0)
      INTEGER DADVRN
      PARAMETER (DADVRN = 1)
      INTEGER DADNRO
      PARAMETER (DADNRO = 2)
      INTEGER DADNCO
      PARAMETER (DADNCO = 3)
      INTEGER DADPES
      PARAMETER (DADPES = 4)
      INTEGER DADADR
      PARAMETER (DADADR = 0)
      INTEGER CTDHDS
      PARAMETER (CTDHDS = 0)
      INTEGER CTDVRN
      PARAMETER (CTDVRN = 1)
      INTEGER CTDNRO
      PARAMETER (CTDNRO = 2)
      INTEGER CTDNCO
      PARAMETER (CTDNCO = 3)
      INTEGER CTDADR
      PARAMETER (CTDADR = 0)
      INTEGER CADHDS
      PARAMETER (CADHDS = 0)
      INTEGER CADVRN
      PARAMETER (CADVRN = 1)
      INTEGER CADNRO
      PARAMETER (CADNRO = 2)
      INTEGER CADNCO
      PARAMETER (CADNCO = 3)
      INTEGER CADADR
      PARAMETER (CADADR = 0)
C======= Include /kloe/soft/off/offline/inc/development/bcs.cin ========
        INTEGER          IW
        REAL             RW(200)
        DOUBLE PRECISION DW(100)
        CHARACTER*4      AW(200)
        INTEGER*2        IW2(400)
        EQUIVALENCE (IW,IW2,RW,DW)
        EQUIVALENCE (IW,AW)
        COMMON /BCS/ IW(200)
C-----------------------------------------------------------------------
C= Include /kloe/soft/off/offline/inc/development/tls/maxstructdim.cin =
      INTEGER     NeleCluMax
      PARAMETER ( NeleCluMax = 2000 )
      INTEGER     MaxNumClu
      PARAMETER ( MaxNumClu = 100 )
      INTEGER    MaxNumVtx
      PARAMETER (MaxNumVtx = 20)
      INTEGER    MaxNumTrkV
      PARAMETER (MaxNumTrkV = 30)
      INTEGER    MaxNumTrk
      PARAMETER (MaxNumTrk = 100)
      INTEGER    MaxNumDHSP
      PARAMETER (MaxNumDHSP = 1000)
      Integer    nMaxDC
      Parameter (nMaxDC=1500)
      INTEGER    NQihiMax
      PARAMETER (NQihiMax=1000)
      INTEGER    NQcalMax
      PARAMETER (NQcalMax= 32 )
      INTEGER    MaxNumFirstHit
      PARAMETER (MaxNumFirstHit = 300 )
      INTEGER    MaxNtrkGen
      PARAMETER (MaxNtrkGen =50)
      INTEGER    MaxNvtxGen
      PARAMETER (MaxNvTxGen =30)
      integer TriggerElements
      parameter (TriggerElements = 300 )
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
         INTEGER iv(MaxNumtrkv)
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
         INTEGER mother(MaxNvtxGen)
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
C=== Include /kloe/soft/off/offline/inc/development/tls/nvostru.cin ====
      INTEGER MaxNumKNVO
      PARAMETER (MaxNumKNVO = 40)       
      TYPE KNVOStru
      SEQUENCE
      integer n                
      integer iknvo(MaxNumKNVO) 
      real    px(MaxNumKNVO)   
      real    py(MaxNumKNVO)
      real    pz(MaxNumKNVO)
      integer pid(MaxNumKNVO)  
      integer bank(MaxNumKNVO) 
      integer vlinked(MaxNumKNVO) 
      END TYPE
      TYPE (KNVOStru) KNVO
      INTEGER MaxNumVNVO
      PARAMETER (MaxNumVNVO = 40)       
      TYPE VNVOStru
      SEQUENCE
      integer n                 
      integer ivnvo(MaxNumVNVO) 
      real vx(MaxNumVNVO)  
      real vy(MaxNumVNVO)
      real vz(MaxNumVNVO)
      integer kori(MaxNumVNVO) 
      integer idvfs(MaxNumVNVO) 
      integer nknv(MaxNumVNVO) 
      integer fknv(MaxNumVNVO) 
      END TYPE
      TYPE (VNVOStru) VNVO
      INTEGER MaxNumVNVOb
      PARAMETER (MaxNumVNVOb = 40)      
      TYPE VNVBStru
      SEQUENCE
      integer n                 
      integer ibank(MaxNumVNVOb) 
      END TYPE
      TYPE (VNVBStru) VNVB
      INTEGER MaxNumINVO
      PARAMETER (MaxNumINVO = 40)       
      TYPE INVOStru
      SEQUENCE
      integer n                 
      integer iinvo(MaxNumINVO) 
      integer iclps(MaxNumINVO) 
      real xi(MaxNumINVO)     
      real yi(MaxNumINVO)
      real zi(MaxNumINVO)
      real ti(MaxNumINVO)     
      real lk(MaxNumINVO)     
      real sigmalk(MaxNumINVO) 
      END TYPE
      TYPE (INVOStru) INVO
C==== Include /kloe/soft/off/offline/inc/development/ecl/nvbnk.cin =====
        integer    KNVO_HDSZ    
        parameter (KNVO_HDSZ=3)
        integer    KNVHSZ       
        parameter (KNVHSZ=0)
        integer    KNVVER       
        parameter (KNVVER=1)
        integer    KNVDSZ       
        parameter (KNVDSZ=2)
        integer     KNVPCO
        parameter  (KNVPCO = 0)
        integer     KNVERR
        parameter  (KNVERR = 3)
        integer     KNVTYP
        parameter  (KNVTYP = 12)
        integer     KNVBNK
        parameter  (KNVBNK = 13)
        integer     KNVLNK
        parameter  (KNVLNK = 14)
        integer    VNVO_VERS    
        parameter (VNVO_VERS=1)
        integer    VNVO_HDSZ    
        parameter (VNVO_HDSZ=3)
        integer    VNVHSZ       
        parameter (VNVHSZ=0)
        integer    VNVVER       
        parameter (VNVVER=1)
        integer    VNVDSZ       
        parameter (VNVDSZ=2)
        integer     VNVXYZ
        parameter  (VNVXYZ = 0)
        integer     VNVERR
        parameter  (VNVERR = 3)
        integer     VNVKOR
        parameter  (VNVKOR = 12)
        integer     VNVDVF
        parameter  (VNVDVF = 13)
        integer     VNVPGE
        parameter  (VNVPGE = 14)
        integer    INVO_VERS    
        parameter (INVO_VERS=1)
        integer    INVO_HDSZ    
        parameter (INVO_HDSZ=3)
        integer    INVHSZ
        parameter (INVHSZ=0)
        integer    INVVER       
        parameter (INVVER=1)
        integer    INVDSZ       
        parameter (INVDSZ=2)
        integer     INVGAM
        parameter  (INVGAM = 0)
        integer     INVERR
        parameter  (INVERR = 4)
        integer     INVCLP
        parameter  (INVCLP = 8)
        integer     INVLKI
        parameter  (INVLKI = 9)
        integer     INVSLK
        parameter  (INVSLK = 10)
        integer NVO_VERS
        COMMON /NVOVRS/  NVO_VERS
C====== Include /kloe/soft/off/offline/inc/development/oferco.cin ======
      INTEGER    OFSUCC                 
      PARAMETER (OFSUCC = 0)
      INTEGER    OFFAIL                 
      PARAMETER (OFFAIL = 1)
      INTEGER    OFFINF                 
      PARAMETER (OFFINF = 2)
      INTEGER    OFERRE                 
      PARAMETER (OFERRE = 3)            
      INTEGER    OFFEOF                 
      PARAMETER (OFFEOF = 4)
C====================== Include interfstruct.inc =======================
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
          REAL    Tcl(MaxNumClu)
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
          INTEGER iv(MaxNumtrkv)
          INTEGER trknumv(MaxNumtrkv)
          REAL    CurV(MaxNumtrkv)
          REAL    PhiV(MaxNumtrkv)
          REAL    CotV(MaxNumtrkv)
          REAL    PxTV(MaxNumtrkv)
          REAL    PyTV(MaxNumtrkv)
          REAL    PzTV(MaxNumtrkv)
          INTEGER nv
          INTEGER vtx(MaxNumVtx)
          REAL    xv(MaxNumVtx)
          REAL    yv(MaxNumVtx)
          REAL    zv(MaxNumVtx)
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
          INTEGER pidmc(MaxNtrkGen)
          INTEGER virmom(MaxNtrkGen)
          REAL    pxmc(MaxNtrkGen)
          REAL    pymc(MaxNtrkGen)
          REAL    pzmc(MaxNtrkGen)
          REAL    themc(MaxNtrkGen)
          REAL    phimc(MaxNtrkGen)
          INTEGER vtxmc(MaxNtrkGen)
          INTEGER nvtxmc
          INTEGER kinmom(MaxNvtxGen)
          INTEGER mother(MaxNvtxGen)
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
          REAL    KchMC(9)
          REAL    KneMC(9)
          REAL    KchRec(9)
          REAL    KneRec(9)
          REAL    KchBoost(9)
          REAL    ip(3)
          REAL    ipmc(3)
          INTEGER vtaken(3)
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
          REAL    cosTrk
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
          INTEGER omegatruth_top
          INTEGER clu_index(4)
          INTEGER chtrk_index(2)
          INTEGER chvtx_index
          REAL P4gam_w(4,4)
          REAL P4ch_w(2,4)
          REAL P4ch_tot_w(4)
          REAL om_dist(4)
      END TYPE
      TYPE ( interfstru) Interf
      COMMON / InterfCommon / Interf,Evt
C-----------------------------------------------------------------------
C External functions
      INTEGER    ANGIST, ANGOHS, ANGPAR
      INTEGER    BLOCAT, BNUMB, BDATA, BNEXT
      INTEGER    GetCluStru, TrkV2Stru, GetNvo
      INTEGER    Recover_Splitting
      INTEGER    GetEmcPassed, fit_interf_pmoo
      REAL       VMOD, VDOTN, VDIST, ANPTRG
C-----------------------------------------------------------------------
C Exported variables
      INTEGER ChVtxId,TrkIdx(2),CluIdx(4),simok
      REAL    Trk1(3),Trk2(3),ChaVtx(3),NeuVtx(3),PhiVtx(3)
C-----------------------------------------------------------------------
C Local variables
C
      INTEGER    Status, InStat, ParSet, MinHisId, MaxHisId,VtxId
      INTEGER    VNVOIdx, VNVONxt, INDDAT, BnkNum, Nknvo
      INTEGER    ILoop, JLoop, Ntracks, Ngammas, FstBank, Nrunold
      INTEGER    TrkId(MaxNumVNVOb), CluId(MaxNumVNVOb),i
      REAL       CovPhiVtx(9), CovNeuVtx(9)
      CHARACTER  ErrTxt*100
C
C-----------------------------------------------------------------------
      Status = GETCLUSTRU(cluster)
      Status = Recover_Splitting(cluster)
      Status = TRKV2STRU(vtx,trkv,trk,trkmc)
      Status = GetNvo(KNVO,VNVO,VNVB,INVO)
C
C Unpack VNVO bank to get covariance matrix for vertex
C
      Nknvo = 0
      Status = BLOCAT(IW,'VNVO',-1,VNVOIdx,INDDAT)
      IF( Status.EQ.YESUCC )THEN
         DO WHILE( Status.EQ.YESUCC .AND. VNVOIdx.GT.0 .AND.
     &      Nknvo.LT.MaxNumVnvo )
C Pointer to data
            Status = BDATA(IW,VNVOIdx,INDDAT)  
C Number of current VNVO bank
            Status = BNUMB(IW,VNVOIdx,BnkNum)  
            IF( Status.EQ.YESUCC )THEN
               Nknvo = Nknvo + 1
C Phi vertex bank number
               IF( BnkNum.EQ.1 )THEN           
                  DO ILoop = 4,12
                     CovPhiVtx(ILoop-3) =
     &                    RW(INDDAT+VNVO_HDSZ+VNVXYZ-1+ILoop)
                  ENDDO
               ENDIF
C Neutral vertex bank number
               IF( BnkNum.EQ.3 )THEN           
                  DO ILoop = 4,12
                     CovNeuVtx(ILoop-3) =
     &                    RW(INDDAT+VNVO_HDSZ+VNVXYZ-1+ILoop)
                  ENDDO
               ENDIF
            ENDIF
C Point to next bank
            Status = BNEXT(IW,VNVOIdx,VNVONxt) 
            VNVOIdx = VNVONxt
         ENDDO
      ENDIF
C
C------------------------------------------------------------------------------
C Loop on VNVO to get verticies information
C------------------------------------------------------------------------------
C
      VtxId   = 0
      Ntracks = 0
      Ngammas = 0
      CALL VZERO(PhiVtx,3)
      CALL VZERO(ChaVtx,3)
      CALL VZERO(NeuVtx,3)
      CALL VZERO(CluId,MaxNumVNVOb)
      CALL VZERO(TrkId,MaxNumVNVOb)
C
      DO ILoop = 1,VNVO%N
C
C Phi vertex
C
         IF( VNVO%Ivnvo(ILoop).EQ.1 )THEN
            PhiVtx(1) = VNVO%Vx(ILoop)
            PhiVtx(2) = VNVO%Vy(ILoop)
            PhiVtx(3) = VNVO%Vz(ILoop)
C
C Charged vertex
C
         ELSEIF( VNVO%Ivnvo(ILoop).EQ.2 )THEN
C Vertex index
            VtxId = VNVO%Idvfs(ILoop)           
C Vertex coordinates
            ChaVtx(1) = VNVO%Vx(ILoop)          
            ChaVtx(2) = VNVO%Vy(ILoop)
            ChaVtx(3) = VNVO%Vz(ILoop)
C Get DVFS index of tracks connected to the vertex
            DO JLoop = 1,TrkV%N
               IF( TrkV%Iv(JLoop).EQ.VtxId )THEN
                  Ntracks = Ntracks + 1
                  TrkId(Ntracks) = JLoop
               ENDIF
            ENDDO
C Check number of connected
            IF( Ntracks.NE.VNVO%Nknv(ILoop) )THEN 
C tracks from VNVO and TrkV
               WRITE(ErrTxt,'(a,i3,a,i3,a)')      
     &              'Different number of connected tracks:',
     &              VNVO%Nknv(ILoop),' from VNVO',Ntracks,' from TrkV'
               CALL ERLOGR('PMOOEV',ERWARN,0,Status,ErrTxt)
            ENDIF
C
C Neutral vertex
C
         ELSEIF( VNVO%Ivnvo(ILoop).EQ.3 )THEN
C Vertex coordinates
            NeuVtx(1) = VNVO%Vx(ILoop)          
            NeuVtx(2) = VNVO%Vy(ILoop)
            NeuVtx(3) = VNVO%Vz(ILoop)
C Number of connected clusters
            Ngammas = VNVO%Nknv(ILoop)          
            FstBank = VNVO%Fknv(ILoop)
C Get CLPS index of clusters connected to the vertex
            DO JLoop = FstBank,FstBank+Ngammas-1
              CluId(JLoop-FstBank+1) = VNVB%Ibank(JLoop)
            ENDDO
         ENDIF
      ENDDO
C
      simok = 0
      IF( Ntracks.EQ.2 )THEN
         IF( TrkV%Cur(TrkId(1)).NE.TrkV%Cur(TrkId(2)) )THEN
            IF( Ngammas.EQ.4 )THEN
               simok = 1
               DO i= 1,4
                 CluIdx(i) = CluId(i)
               ENDDO
               ChVtxId = VtxId
               DO i= 1,2
                 TrkIdx(i) = TrkId(i)
               ENDDO
               Trk1(1) = COS( TrkV%Phi(TrkIdx(1) )) / ABS( TrkV%Cur(TrkI
     1  dx(1) ))*1000.                                                  
               Trk1(2) = SIN( TrkV%Phi(TrkIdx(1) )) / ABS( TrkV%Cur(TrkI
     1  dx(1) ))*1000.                                                  
               Trk1(3) =      TrkV%Cot(TrkIdx(1) )  / ABS( TrkV%Cur(TrkI
     1  dx(1) ))*1000.                                                  
               Trk2(1) = COS( TrkV%Phi(TrkIdx(2) )) / ABS( TrkV%Cur(TrkI
     1  dx(2) ))*1000.                                                  
               Trk2(2) = SIN( TrkV%Phi(TrkIdx(2) )) / ABS( TrkV%Cur(TrkI
     1  dx(2) ))*1000.                                                  
               Trk2(3) =      TrkV%Cot(TrkIdx(2) )  / ABS( TrkV%Cur(TrkI
     1  dx(2) ))*1000.                                                  
            ENDIF
         ENDIF
      ENDIF
C
      RETURN
      END
