C==============================================================================
C  KLSPM00
C==============================================================================
C
C  Description:
C  ------------
C  A_C module to study kaon interferometry for the pi+pi-pi0pi0 final state.
C  It must be executed after PROD2NTU module, since the output of the
C  analysis is stored in a new block added to the PROD2NTU ntuple.
C
C  Language:
C  ---------
C  KLOE Fortran
C
C  Author:
C  -------
C
C  Creation date:
C  --------------
C
C  Modification history:
C  ---------------------
C
C==============================================================================
C==============================================================================
C==============================================================================
      SUBROUTINE PM00IN
C------------------------------------------------------------------------------
C
C  Description:
C  ------------
C
C------------------------------------------------------------------------------
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
C
C Local declarations
C
      INTEGER i
C------------------------------------------------------------------------------
      nmctruth0   = 0
      nmctruth(1) = 0
      nmctruth(2) = 0
      nmctruth(3) = 0
      nmctruth(4) = 0
      nmctruth(5) = 0
      nmctruth(6) = 0
      nmctruth(7) = 0
      selected = 0
      ntruth = 0
      ntruthregg = 0
      ntruththreee = 0
      ntruthsemii = 0
      ntruthomegaa = 0
      ntruthelsee = 0
      analised = 0
      DO i = 1,testcounter
         ErrFlagCount(i) = 0
         ErrFlagCountBcg(i) = 0
         counter(i) = 0
         counterBcg(i) = 0
      ENDDO
      messageErr(1)  = ' GetKslEvent'
      messageErr(2)  = ' find_kchrec'
      messageErr(3)  = ' cor_ip_boost'
      messageErr(4)  = ' ch vtxdist'
      messageErr(5)  = ' find_neuclu'
      messageErr(6)  = ' find_neuvtx'
      messageErr(7)  = ' ne vtxdist'
      messageErr(8)  = ' GetTdiff KchBoost KneRecLor'
      messageErr(9)  = ' GetTdiff KchBoost KneRec'
      messageErr(10) = ' GetTdiff KchRec   KneRec'
      messageErr(11) = ' Knerec(0,0,0,...)'
      messageErr(12) = ' find_omegaminv'
      messageErr(13) = ' fit_interf_pmoo'
      messageErr(14) = ' find_neuvtx fit'
      messageErr(15) = ' GetTdiff fit'
      messageErr(16) = ' find_omegaminv fit'
      messageErr(17) = ' do_cuts'
      messageErr(testcounter) = ' unknown'
      messageCut(1) = '  chi2'
      messageCut(2) = '  pi0fit'
      messageCut(3) = '  trcv'
      messageCut(4) = '  Mkne'
      messageCut(5) = '  Mkch'
      messageCut(6) = '  Qmiss'
      truthmc    = .FALSE.
      cutneu     = .FALSE.
      cutchar    = .FALSE.
      onlytruemc = .FALSE.
      makecuts   = .TRUE.
C
      CALL when
      WRITE(*,*) 'Program started'
C
      END
C
C==============================================================================
C==============================================================================
C==============================================================================
      SUBROUTINE PM00RI
C------------------------------------------------------------------------------
C
C  Description:
C  ------------
C
C------------------------------------------------------------------------------
C
        IMPLICIT NONE
C
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
      INTEGER    ANPIST
C
C Local declarations
C
      INTEGER    Status
C
C------------------------------------------------------------------------------
C
C Put A_C Error Code to SUCCESS at Run_Init
C
      Status = ANPIST(ANSUCC)
C
C Set initialization parameters
C
      END
C
C==============================================================================
C==============================================================================
C==============================================================================
       SUBROUTINE PM00EV
C------------------------------------------------------------------------------
C
        IMPLICIT NONE
C
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
C------------------------------------------------------------------------------
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
C------------------------------------------------------------------------------
C
C External functions
C
      INTEGER    ANGIST, ANGOHS, ANGPAR, ANPTRG,sistle
C
C Local declarations
C
      REAL       FitPar(MaxNumFitPar),ErrPar(MaxNumFitPar)
      REAL       BkgFitPar(MaxNumFitPar),BkgErrPar(MaxNumFitPar)
      REAL       Rttmp,Chi2fit
      REAL*8     Constr(MaxNumConstr)
      REAL*8     D(MaxNumFitPar,MaxNumConstr)
      REAL*8     MTX(MaxNumFitPar+MaxNumConstr,
     &               MaxNumFitPar+MaxNumConstr)
      REAL*8     Z(MaxNumFitPar+MaxNumConstr)
      REAL*8     ConstrOmega(MaxNumComega)
      REAL*8     DOmega(MaxNumFitPar,MaxNumComega)
      REAL*8     MTXOmega(MaxNumFitPar+MaxNumComega,
     &                    MaxNumFitPar+MaxNumComega)
      REAL*8     ZOmega(MaxNumFitPar+MaxNumComega)
      REAL       PxVtx,PyVtx,PzVtx,ErrXv,ErrYv,ErrZv,Eres,Tres,Tcst,Zres
     1  ,Xv,Yv,Zv                                                       
      REAL       ErrPxVtx,ErrPyVtx,ErrPzVtx,SqrtS,ErrSqrtS
      REAL       clustersE(4),clustersX(4),clustersY(4),clustersZ(4)
      REAL       clustersT(4),FitPar6(6),ipcoor(3),PgamRecTaken4x4(4,4)
      REAL       PgamRecTaken4x4fit(4,4),test(10),PhiP(3)
      INTEGER    clustersN(4),clustersNwrong(4)
      INTEGER    nclusters,nclusterswrong
      INTEGER    InStat, ParSet, MinHisId, MaxHisId
      INTEGER    Status,i,KLoop,ILoop,loopcount
      LOGICAL    HisFlg
      INTEGER    kslNumb,icycle
      PARAMETER( icycle = 1000000 )
      PARAMETER( kslNumb =  2 )
C
C------------------------------------------------------------------------------
C
      Interf%ErrId=0
      Interf%CutId=0
C IF ANPACK is set to ANYBER (i.e.
      Status = ANGIST(InStat)         
C YBos ERror), skip the event routine
      IF( InStat.EQ.ANYBER ) RETURN   
C-----------------------------------------------------------------------
C MC or DATA
C-----------------------------------------------------------------------
C  MC
      IF( ExpTyp.EQ.EXOFSI )THEN    
          mcflag=.TRUE.
      ELSE
          mcflag=.FALSE.
      ENDIF
C------------------------------------------------------------------------------
      analised = analised + 1
      IF( MOD(analised,icycle).eq.0 ) THEN
         CALL summary
      ENDIF
C------------------------------------------------------------------------------
      EventOK  = .FALSE.
      ErrFlag  = .FALSE.
      CALL clearstruct
      CALL fillstruct
C------------------------------------------------------------------------------
C Selection of good events:
C------------------------------------------------------------------------------
C
C Look for MC Phi --> KsKl --> pi+pi-pi0pi0 events
C
      IF( mcflag )THEN
        CALL GetKslEvent(Interf%ntmc,Interf%mother,
     &  Interf%vtxmc,Interf%pidmc,Interf%xvmc,Interf%yvmc,Interf%zvmc,
     &  Interf%pxmc,Interf%pymc,Interf%pzmc,Interf%nvtxmc,
     &  Interf%ipmc,Interf%KchMC,Interf%KneMC,Interf%DtMC,Interf%DlMC,
     &  truthmc,truthregg,truthsemii,truththreee,truthomegaa,truthelsee)
          IF( ErrFlag ) THEN
            CALL statisticsErr(1)
            goto 66
          ENDIF
      ELSE
        truthmc = .FALSE.
      ENDIF
      IF( truthmc ) THEN
        ntruth = ntruth + 1
      ELSE IF( truthregg ) THEN
        ntruthregg = ntruthregg + 1
      ELSE IF( truthsemii ) THEN
        ntruthsemii = ntruthsemii + 1
      ELSE IF( truththreee ) THEN
        ntruththreee = ntruththreee + 1
      ELSE IF( truthomegaa ) THEN
        ntruthomegaa = ntruthomegaa + 1
      ELSE IF( truthelsee ) THEN
      ntruthelsee = ntruthelsee + 1
      ENDIF
      IF(onlytruemc.and.mcflag.and..not.truthmc) RETURN
C------------------------------------------------------------------------------
      IF( truthmc ) THEN
       Interf%RcMC  = SQRT( (Interf%KchMC(7)-Interf%ipmc(1))**2 +
     &                      (Interf%KchMC(8)-Interf%ipmc(2))**2 +
     &                      (Interf%KchMC(9)-Interf%ipmc(3))**2 )
       Interf%RtcMC = SQRT( (Interf%KchMC(7)-Interf%ipmc(1))**2 +
     &                      (Interf%KchMC(8)-Interf%ipmc(2))**2 )
       Interf%RnMC  = SQRT( (Interf%KneMC(7)-Interf%ipmc(1))**2 +
     &                      (Interf%KneMC(8)-Interf%ipmc(2))**2 +
     &                      (Interf%KneMC(9)-Interf%ipmc(3))**2 )
       Interf%RtnMC = SQRT( (Interf%KneMC(7)-Interf%ipmc(1))**2 +
     &                      (Interf%KneMC(8)-Interf%ipmc(2))**2 )
      ENDIF
C------------------------------------------------------------------------------
      CALL getsim(Interf%ChVtxId,Interf%TrkIdx,Interf%Trkk1,
     &            Interf%Trkk2,Interf%ChaVtx,Interf%CluIdx,
     &            Interf%NeuVtx,Interf%PhiVtx,Interf%simok)
C------------------------------------------------------------------------------
      IF( mcflag ) CALL check_isr(Interf%ntmc,Interf%pidmc,
     &        Interf%virmom,Interf%vtxmc,Interf%mother,Interf%mcisr)
C------------------------------------------------------------------------------
C add better selection tracs based on chi2_vtx, chi2_trk, chi2_Kmass
      CALL find_kchrec(Interf%qualv,Interf%nv,Interf%ntv,Interf%IV,
     &                 Interf%PxTV,
     &                 Interf%PyTV,Interf%PzTV,Interf%xv,Interf%yv,
     &                 Interf%zv,Interf%vtaken,Interf%KchRec,
     &                 Interf%trk1,Interf%trk2,Interf%cosTrk)
      IF( ErrFlag ) THEN
         CALL statisticsErr(2)
         IF( truthmc ) THEN
            goto 66
         ELSE
            RETURN
         ENDIF
      ENDIF
C------------------------------------------------------------------------------
C should be gauss with width 0.3
      IF ( mcflag ) Interf%Bpz = 0. 
      CALL cor_ip_boost(Interf%KchRec,Interf%Bpx,Interf%Bpy,Interf%Bpz,
     &            Interf%Bx,Interf%By,Interf%Bz,Interf%Broots,
     &            Interf%trk1,Interf%trk2,
     &            Interf%KchBoost,Interf%ip,Interf%chdist,Interf%Qmiss)
      IF( ErrFlag ) THEN
         CALL statisticsErr(3)
         IF( truthmc ) THEN
            goto 66
         ELSE
            RETURN
         ENDIF
      ENDIF
      IF( abs(Interf%Bz-Interf%ip(3)).gt.2 ) Interf%ip(3)=Interf%Bz
C------------------------------------------------------------------------------
      Interf%Rc  = SQRT( (Interf%KchBoost(7)-Interf%ip(1))**2 +
     &                   (Interf%KchBoost(8)-Interf%ip(2))**2 +
     &                   (Interf%KchBoost(9)-Interf%ip(3))**2 )
      Interf%Rtc = SQRT( (Interf%KchBoost(7)-Interf%ip(1))**2 +
     &                   (Interf%KchBoost(8)-Interf%ip(2))**2 )
      IF( Interf%Rc.gt.vtxdist ) THEN
         ErrFlag = .TRUE.
         CALL statisticsErr(4)
         IF( truthmc ) THEN
            goto 66
         ELSE
            RETURN
         ENDIF
      ENDIF
C------------------------------------------------------------------------------
      CALL find_neuclu(Interf%nclu,Interf%ntcl,Interf%asscl,
     &                 Interf%ncll,Interf%ncl)
      IF( ErrFlag ) THEN
         CALL statisticsErr(5)
         IF( truthmc ) THEN
            goto 66
         ELSE
            RETURN
         ENDIF
      ENDIF
C ERYK: chi2 kinfit
      CALL find_neuvtx(Interf%Bpx,Interf%Bpy,Interf%Bpz,Interf%Broots,
     &                 Interf%KchBoost,Interf%ncl,Interf%enecl,
     &                 Interf%ncll,Interf%xacl,Interf%yacl,Interf%zacl,
     &                 Interf%ip,Interf%tcl,Interf%cldist,
     &                 Interf%KneRecLor,Interf%trc,Interf%nclwrong,
     &                 Interf%ncllwrong,Interf%KneRec,
     &                 Interf%minv4gam,Interf%pi0,Interf%g4taken,
     &                 Interf%trcv,PgamRecTaken4x4,Interf%gpairtaken,
     &                 test,Interf%g4vtxerr)
      IF( ErrFlag ) THEN
         CALL statisticsErr(6)
         IF( truthmc ) THEN
            goto 66
         ELSE
            RETURN
         ENDIF
      ENDIF
      DO i=1,10
         Interf%test(i)=test(i)
      ENDDO
C------------------------------------------------------------------------------
      Interf%Rn  = SQRT( (Interf%KneRecLor(7)-Interf%ip(1))**2 +
     &                   (Interf%KneRecLor(8)-Interf%ip(2))**2 +
     &                   (Interf%KneRecLor(9)-Interf%ip(3))**2 )
      Interf%Rtn = SQRT( (Interf%KneRecLor(7)-Interf%ip(1))**2 +
     &                   (Interf%KneRecLor(8)-Interf%ip(2))**2 )
      IF( Interf%Rn.gt.vtxdist ) THEN
         ErrFlag = .TRUE.
         CALL statisticsErr(7)
         IF( truthmc ) THEN
            goto 66
         ELSE
            RETURN
         ENDIF
      ENDIF
      PhiP(1) = Interf%Bpx
      PhiP(2) = Interf%Bpy
      PhiP(3) = Interf%Bpz
      CALL GetTdiff(Interf%KchBoost,Interf%KneRecLor,
     &        Interf%ip,PhiP,Interf%Broots,
     &        Interf%DlBoostLor,Interf%DtBoostLor)
      IF( ErrFlag ) THEN
         CALL statisticsErr(8)
         IF( truthmc ) THEN
            goto 66
         ELSE
            RETURN
         ENDIF
      ENDIF
      CALL GetTdiff(Interf%KchBoost,Interf%KneRec,
     &        Interf%ip,PhiP,Interf%Broots,
     &        Interf%DlBoostRec,Interf%DtBoostRec)
      IF( ErrFlag ) THEN
         CALL statisticsErr(9)
         IF( truthmc ) THEN
            goto 66
         ELSE
            RETURN
         ENDIF
      ENDIF
      CALL GetTdiff(Interf%KchRec,Interf%KneRec,
     &        Interf%ip,PhiP,Interf%Broots,
     &        Interf%DlRec,Interf%DtRec)
      IF( ErrFlag ) THEN
         CALL statisticsErr(10)
         IF( truthmc ) THEN
            goto 66
         ELSE
            RETURN
         ENDIF
      ENDIF
C------------------------------------------------------------------------------
      DO i=1,4
         Interf%PgamRec1(i)=PgamRecTaken4x4(1,i)
         Interf%PgamRec2(i)=PgamRecTaken4x4(2,i)
         Interf%PgamRec3(i)=PgamRecTaken4x4(3,i)
         Interf%PgamRec4(i)=PgamRecTaken4x4(4,i)
      ENDDO
      IF( Interf%KneRec(1).eq.0..and.Interf%KneRec(2).eq.0..and.
     &    Interf%KneRec(3).eq.0. ) THEN
         ErrFlag = .TRUE.
         CALL statisticsErr(11)
         IF( truthmc ) THEN
            goto 66
         ELSE
            RETURN
         ENDIF
      ENDIF
C------------------------------------------------------------------------------
C Omega Topology Selection
C------------------------------------------------------------------------------
      CALL OmegaTopologySel(Interf%nv,Interf%nclu,Interf%ntv,Interf%Bx,I
     1  nterf%By,Interf%Bz,                                             
     &          Interf%CurV,Interf%PhiV,Interf%CotV,Interf%EneCl,Interf%
     1  Xcl,Interf%Ycl,Interf%Zcl,                                      
     &      Interf%Tcl,Interf%Xacl,Interf%Yacl,Interf%Zacl,Interf%xv,Int
     1  erf%yv,Interf%zv,                                               
     &      Interf%vtx,Interf%iv,
     &          Interf%chtrk_index,Interf%clu_index,Interf%chvtx_index,I
     1  nterf%omegatruth_top)                                           
      DO I=1,4
         Interf%om_dist(I) = SQRT( (Interf%Xcl(Interf%clu_index(I)) - In
     1  terf%xv(Interf%chvtx_index))**2 +                               
     &               (Interf%Ycl(Interf%clu_index(I)) - Interf%yv(Interf
     1  %chvtx_index))**2 +                                             
     &               (Interf%Zcl(Interf%clu_index(I)) - Interf%zv(Interf
     1  %chvtx_index))**2)                                              
         Interf%P4gam_w(I,1) = Interf%EneCl(Interf%clu_index(I))*
     &                  (Interf%Xcl(Interf%clu_index(I)) - Interf%xv(Int
     1  erf%chvtx_index))/Interf%om_dist(I)                             
         Interf%P4gam_w(I,2) = Interf%EneCl(Interf%clu_index(I))*
     &                  (Interf%Ycl(Interf%clu_index(I)) - Interf%yv(Int
     1  erf%chvtx_index))/Interf%om_dist(I)                             
         Interf%P4gam_w(I,3) = Interf%EneCl(Interf%clu_index(I))*
     &                  (Interf%Zcl(Interf%clu_index(I)) - Interf%zv(Int
     1  erf%chvtx_index))/Interf%om_dist(I)                             
         Interf%P4gam_w(I,4) = Interf%EneCl(Interf%clu_index(I))
      ENDDO
      DO I=1,2
         Interf%P4ch_w(I,1) = 1000.*COS(Interf%PhiV(Interf%chtrk_index(I
     1  )))/ABS(Interf%CurV(Interf%chtrk_index(I)))                     
         Interf%P4ch_w(I,2) = 1000.*SIN(Interf%PhiV(Interf%chtrk_index(I
     1  )))/ABS(Interf%CurV(Interf%chtrk_index(I)))                     
         Interf%P4ch_w(I,3) = 1000.*Interf%CotV(Interf%chtrk_index(I))/A
     1  BS(Interf%CurV(Interf%chtrk_index(I)))                          
         Interf%P4ch_w(I,4) = SQRT( Interf%P4ch_w(I,1)**2 + Interf%P4ch_
     1  w(I,2)**2 + Interf%P4ch_w(I,3)**2 + Mpip**2)                    
      ENDDO
      Interf%P4ch_tot_w(1) = Interf%P4ch_w(1,1) + Interf%P4ch_w(2,1)
      Interf%P4ch_tot_w(2) = Interf%P4ch_w(1,2) + Interf%P4ch_w(2,2)
      Interf%P4ch_tot_w(3) = Interf%P4ch_w(1,3) + Interf%P4ch_w(2,3)
      Interf%P4ch_tot_w(4) = Interf%P4ch_w(1,4) + Interf%P4ch_w(2,4)
C------------------------------------------------------------------------------
C Track parameters
      DO ILoop = 1,2                     
CEvt%Trkv%px(Interf%vtaken(ILoop+1))
         FitPar6(3*(ILoop-1)+1) = Interf%P4ch_w(ILoop,1) 
CEvt%Trkv%py(Interf%vtaken(ILoop+1))
         FitPar6(3*(ILoop-1)+2) = Interf%P4ch_w(ILoop,2) 
CEvt%Trkv%pz(Interf%vtaken(ILoop+1))
         FitPar6(3*(ILoop-1)+3) = Interf%P4ch_w(ILoop,3) 
      ENDDO
      CALL find_omegaminv(FitPar6,Interf%P4gam_w,Interf%P4ch_tot_w,
     &          Interf%ominv,Interf%PpioOmega,Interf%P4PriRest,
     &          Interf%test(21))
      IF( ErrFlag ) THEN
         CALL statisticsErr(12)
         IF( truthmc ) THEN
            goto 66
         ELSE
            RETURN
         ENDIF
      ENDIF
C------------------------------------------------------------------------------
C Kinematic Fit
C------------------------------------------------------------------------------
      DO ILoop = 1,MaxNumFitPar
         FitPar(ILoop) = 0.
         ErrPar(ILoop) = 0.
      ENDDO
CC Set input values
CC
CC Set MC or DATA parameters
CC
CCCCCCCCCCCCC eryk
CC      IF( mcflag )THEN    ! MC
CC         Eres = 0.05
CC         Tres = 0.05
CC         Tcst = 0.00
CC         Zres = 0.85
CC         Xv = 0
CC         Yv = 0.
CC         Zv = 0.
CC         ErrXv = 0.2
CC         ErrYv = 0.002
CC         ErrZv = 2.0
CC         PxVtx = -12.75
CC         ErrPxVtx = 0.006
CC         SqrtS    = 1020.0
CC         ErrSqrtS =   0.45
CC      ELSE                          ! Data
CC         Eres = 0.060
CC         Tres = 0.054
CC         Tcst = 0.150
CC         Zres = 1.400
CC         Xv = Interf%Bx
CC         ErrXv = SQRT( Interf%Bwidpx**2 + Interf%Blumx**2 )
CC         Yv = Interf%By
CC         ErrYv = Interf%Bsy
CC         Zv = Interf%Bz
CC         ErrZv = SQRT( Interf%Bsz**2 + Interf%Blumz**2 )
CC         PxVtx = Interf%Bpx
CCc           ErrPxVtx = BPos%ElarPx
CC         ErrPxVtx = Interf%Bwidpx
CC         SqrtS    = Interf%Broots
CC         ErrSqrtS = Interf%BrootsErr
CC      ENDIF
CCCCCCCCCCCCC eryk
         IF( mcflag ) THEN
            Tcst = 0.
         ELSE
            Tcst = 0.150
         ENDIF
         Eres = 0.060
         Tres = 0.054
         Zres = 1.400
         Xv = Interf%ip(1)
         ErrXv = SQRT( Interf%Bwidpx**2 + Interf%Blumx**2 )
         Yv = Interf%ip(2)
         ErrYv = Interf%Bsy
         Zv = Interf%ip(3)
C        ErrZv = SQRT( Interf%Bsz**2 + Interf%Blumz**2 )
C        ErrZv = SQRT( Interf%Bwidpz**2 + Interf%Blumz**2 )
         ErrZv = 0.6
         PxVtx = Interf%Bpx
         ErrPxVtx = Interf%Bwidpx
         PyVtx = Interf%Bpy
         ErrPyVtx = Interf%Bwidpy
         PzVtx = Interf%Bpz
         ErrPzVtx = Interf%Bwidpz
C! eryk ERYK dodac PzVtx, ErrPzVtx do kinfit
         SqrtS    = Interf%Broots
         ErrSqrtS = Interf%BrootsErr
CCCCCCCCCCCC eryk
C Track parameters
      DO ILoop = 1,2                     
         FitPar(3*(ILoop-1)+1) = Interf%CurV(Interf%vtaken(ILoop+1))
         FitPar(3*(ILoop-1)+2) = Interf%PhiV(Interf%vtaken(ILoop+1))
         FitPar(3*(ILoop-1)+3) = Interf%CotV(Interf%vtaken(ILoop+1))
C(1.5**2)/2.
         ErrPar(3*(ILoop-1)+1) = 0.025   
C(1.5**2)/2.
         ErrPar(3*(ILoop-1)+2) = 0.013   
C(1.8**2)/2.
         ErrPar(3*(ILoop-1)+3) = 0.012   
      ENDDO
CC
C Cluster parameters
      DO ILoop = 1,4                
         FitPar(6+5*(ILoop-1)+1) =
C Energy
     &       Interf%EneCl(Interf%ncll(Interf%g4taken(ILoop)))  
         FitPar(6+5*(ILoop-1)+2) =
C X coord.
     &        Interf%Xcl(Interf%ncll(Interf%g4taken(ILoop)))  
         FitPar(6+5*(ILoop-1)+3) =
C Y coord.
     &        Interf%Ycl(Interf%ncll(Interf%g4taken(ILoop)))  
         FitPar(6+5*(ILoop-1)+4) =
C Z coord.
     &        Interf%Zcl(Interf%ncll(Interf%g4taken(ILoop)))  
         FitPar(6+5*(ILoop-1)+5) =
C Timing
     &         Interf%Tcl(Interf%ncll(Interf%g4taken(ILoop)))  
         Rttmp = SQRT(
     &         Interf%Xacl(Interf%ncll(Interf%g4taken(ILoop)))**2 +
     &         Interf%Yacl(Interf%ncll(Interf%g4taken(ILoop)))**2 )
         ErrPar(6+5*(ILoop-1)+1) = Eres * SQRT( 1.0E3 *
     &         Interf%EneCl(Interf%ncll(Interf%g4taken(ILoop))) )
         ErrPar(6+5*(ILoop-1)+2) = 1.2
         IF( Rttmp.GT.200. )THEN
            ErrPar(6+5*(ILoop-1)+3) = 1.2
            ErrPar(6+5*(ILoop-1)+4) = Zres * SQRT( 1.0E3 /
     &         Interf%EneCl(Interf%ncll(Interf%g4taken(ILoop))) )
         ELSE
            ErrPar(6+5*(ILoop-1)+3) = Zres * SQRT( 1.0E3 /
     &         Interf%EneCl(Interf%ncll(Interf%g4taken(ILoop))) )
            ErrPar(6+5*(ILoop-1)+4) = 1.2
         ENDIF
         ErrPar(6+5*(ILoop-1)+5) = SQRT( ( Tres * SQRT( 1.0E3 /
     &         Interf%EneCl(Interf%ncll(Interf%g4taken(ILoop)))
     &                                 ) )**2 + Tcst**2 )
      ENDDO
C
C      FitPar(27) = Interf%xv(Interf%vtaken(1))          ! Charged vertex coordi
C      nates                                                                    
C      FitPar(28) = Interf%yv(Interf%vtaken(1))          ! Not used actually
C      FitPar(29) = Interf%zv(Interf%vtaken(1))
C      ErrPar(27) = SQRT( ABS( Interf%VTXcov1(Interf%vtaken(1) ) ) )
C      ErrPar(28) = SQRT( ABS( Interf%VTXcov4(Interf%vtaken(1) ) ) )
C      ErrPar(29) = SQRT( ABS( Interf%VTXcov6(Interf%vtaken(1) ) ) )
CC
C Neutral vertex coordinates
      FitPar(27) = Interf%KneRecLor(7)             
C Improved
      FitPar(28) = Interf%KneRecLor(8)             
      FitPar(29) = Interf%KneRecLor(9)
      ErrPar(27) = 0.8
      ErrPar(28) = 0.8
      ErrPar(29) = 1.1
CC
C Phi vertex coordinates
      FitPar(30) = Xv                    
C Improved
      FitPar(31) = Yv                    
      FitPar(32) = Zv
      ErrPar(30) = ErrXv
      ErrPar(31) = ErrYv
      ErrPar(32) = ErrZv
CC
C      FitPar(36) = SqrtS*0.5                  ! Energy of E+ E- beams
C      FitPar(37) = SqrtS*0.5
C      ErrPar(36) = ErrSqrtS / SQRT(2.)
C      ErrPar(37) = ErrSqrtS / SQRT(2.)
C
C 4-mom Phi meson
      FitPar(33) = SqrtS       
      FitPar(34) = PxVtx
      FitPar(35) = PyVtx
      FitPar(36) = PzVtx
      ErrPar(33) = ErrSqrtS
      ErrPar(34) = ErrPxVtx
      ErrPar(35) = ErrPyVtx
      ErrPar(36) = ErrPzVtx
C      FitPar(38) = PxVtx                 ! Phi momentum along x axis
C      ErrPar(38) = ErrPxVtx
CC
C----------------------------------------------------------------------------
C Set input values for background fit
C SZYMON Using variables found by omega selection from MEMO
C----------------------------------------------------------------------------
      DO i = 1,MaxNumFitPar
         Interf%FitParStart(i) = FitPar(i)
         Interf%ErrParStart(i) = ErrPar(i)
      ENDDO
C
      IF(Interf%omegatruth_top.EQ.1)THEN
C Track parameters
         DO ILoop = 1,2                     
            BkgFitPar(3*(ILoop-1)+1) = Interf%CurV(Interf%chtrk_index(IL
     &oop))
            BkgFitPar(3*(ILoop-1)+2) = Interf%PhiV(Interf%chtrk_index(IL
     &oop))
            BkgFitPar(3*(ILoop-1)+3) = Interf%CotV(Interf%chtrk_index(IL
     &oop))
C(1.5**2)/2.
            BkgErrPar(3*(ILoop-1)+1) = 0.025   
C(1.5**2)/2.
            BkgErrPar(3*(ILoop-1)+2) = 0.013   
C(1.8**2)/2.
            BkgErrPar(3*(ILoop-1)+3) = 0.007   
         ENDDO
C Cluster parameters
         DO ILoop = 1,4                
            BkgFitPar(6+5*(ILoop-1)+1) =
C Energy
     &       Interf%EneCl(Interf%clu_index(ILoop))  
            BkgFitPar(6+5*(ILoop-1)+2) =
C X coord.
     &        Interf%Xcl(Interf%clu_index(ILoop))  
            BkgFitPar(6+5*(ILoop-1)+3) =
C Y coord.
     &        Interf%Ycl(Interf%clu_index(ILoop))  
            BkgFitPar(6+5*(ILoop-1)+4) =
C Z coord.
     &        Interf%Zcl(Interf%clu_index(ILoop))  
            BkgFitPar(6+5*(ILoop-1)+5) =
C Timing
     &         Interf%Tcl(Interf%clu_index(ILoop))  
            Rttmp = SQRT(
     &         Interf%Xacl(Interf%clu_index(ILoop))**2 +
     &         Interf%Yacl(Interf%clu_index(ILoop))**2 )
            BkgErrPar(6+5*(ILoop-1)+1) = Eres * SQRT( 1.0E3 *
     &         Interf%EneCl(Interf%clu_index(ILoop)) )
            BkgErrPar(6+5*(ILoop-1)+2) = 1.2
            IF( Rttmp.GT.200. )THEN
               BkgErrPar(6+5*(ILoop-1)+3) = 1.2
               BkgErrPar(6+5*(ILoop-1)+4) = Zres * SQRT( 1.0E3 /
     &         Interf%EneCl(Interf%clu_index(ILoop)) )
            ELSE
               BkgErrPar(6+5*(ILoop-1)+3) = Zres * SQRT( 1.0E3 /
     &         Interf%EneCl(Interf%clu_index(ILoop)) )
               BkgErrPar(6+5*(ILoop-1)+4) = 1.2
            ENDIF
            BkgErrPar(6+5*(ILoop-1)+5) = SQRT( ( Tres * SQRT( 1.0E3 /
     &         Interf%EneCl(Interf%clu_index(ILoop))
     &                                 ) )**2 + Tcst**2 )
         ENDDO
C Since the background is a prompt event,
         BkgFitPar(27) = Interf%Bx        
C charged and neutral vertex are equal to
         BkgFitPar(28) = Interf%By        
C phi vertex
         BkgFitPar(29) = Interf%Bz        
         BkgErrPar(27) = ErrXv
         BkgErrPar(28) = ErrYv
         BkgErrPar(29) = ErrZv
         BkgFitPar(30) = SqrtS
         BkgFitPar(31) = PxVtx
         BkgFitPar(32) = PyVtx
         BkgFitPar(33) = PzVtx
         BkgErrPar(30) = ErrSqrtS
         BkgErrPar(31) = ErrPxVtx
         BkgErrPar(32) = ErrPyVtx
         BkgErrPar(33) = ErrPzVtx
      ENDIF
CC
CC Perform fit for the signal
CC
      CALL fit_interf_pmoo(FitPar,ErrPar,MaxNumConstr,Chi2Fit,Constr,
     &                       D,MTX,Z,loopcount)
      IF( ErrFlag ) THEN
         CALL statisticsErr(13)
         IF( truthmc ) THEN
            goto 66
         ELSE
            RETURN
         ENDIF
      ENDIF
CC
CC Insert fit results into structure
CC
      Interf%nparfit = MaxNumFitPar
      DO i=1,MaxNumFitPar
         Interf%FitPar(i) = FitPar(i)
         Interf%ErrPar(i) = ErrPar(i)
      ENDDO
      Interf%Niter  = loopcount
      Interf%Chi2   = Chi2Fit
      Interf%Nconstr   = MaxNumConstr
      Interf%KchFit(1) = FitPar(1)+FitPar(4)
      Interf%KchFit(2) = FitPar(2)+FitPar(5)
      Interf%KchFit(3) = FitPar(3)+FitPar(6)
      Interf%KchFit(5) = Interf%KchFit(1)**2 +
     &                   Interf%KchFit(2)**2 +
     &                   Interf%KchFit(3)**2
      Interf%KchFit(4) = sqrt( Interf%KchFit(5) + Mko**2 )
      Interf%KchFit(5) = sqrt( Interf%KchFit(5) )
      Interf%KchFit(6) = Mko
      Interf%KchFit(7) = FitPar(27)
      Interf%KchFit(8) = FitPar(28)
      Interf%KchFit(9) = FitPar(29)
      DO i = 1,3
         ipcoor(i) = FitPar(32+i)
      ENDDO
C
      nclusters=4
      DO i=1,nclusters
         clustersE(i) = FitPar(7 +(i-1)*5)
         clustersX(i) = FitPar(8 +(i-1)*5)
         clustersY(i) = FitPar(9 +(i-1)*5)
         clustersZ(i) = FitPar(10+(i-1)*5)
         clustersT(i) = FitPar(11+(i-1)*5)
         clustersN(i) = i
         clustersNwrong(i) = 0
      ENDDO
      CALL find_neuvtx(FitPar(38),Interf%Bpy,Interf%Bpz,FitPar(36)*2,
     &                 Interf%KchFit,nclusters,clustersE,clustersN,
     &                 clustersX,clustersY,clustersZ,ipcoor,clustersT,
     &                 Interf%cldistFit,Interf%KneRecLorFit,
     &                 Interf%trcFit,Interf%nclwrongFit,
     &                 Interf%ncllwrongFit,Interf%KneFit,
     &                 Interf%minv4gamFit,Interf%pi0Fit,
     &                 Interf%g4takenFit,Interf%trcvFit,
     &                 PgamRecTaken4x4fit,Interf%gpairtakenFit,
     &                 test,Interf%g4vtxerrFit)
      DO i=1,10
         Interf%test(i+10)=test(i)
      ENDDO
      IF( ErrFlag ) THEN
         CALL statisticsErr(14)
         IF( truthmc ) THEN
            goto 66
         ELSE
            RETURN
         ENDIF
      ENDIF
      DO i=1,4
         Interf%PgamRec1fit(i)=PgamRecTaken4x4fit(1,i)
         Interf%PgamRec2fit(i)=PgamRecTaken4x4fit(2,i)
         Interf%PgamRec3fit(i)=PgamRecTaken4x4fit(3,i)
         Interf%PgamRec4fit(i)=PgamRecTaken4x4fit(4,i)
      ENDDO
      Interf%KneFit(7) = FitPar(30)
      Interf%KneFit(8) = FitPar(31)
      Interf%KneFit(9) = FitPar(32)
C chdistFit
      Interf%ipFit(1) = FitPar(33)
      Interf%ipFit(2) = FitPar(34)
      Interf%ipFit(3) = FitPar(35)
C eryk
      CALL GetTdiff(Interf%KchBoost,Interf%KneFit,
     &        Interf%ip,PhiP,Interf%Broots,
     &        Interf%DlFit,Interf%DtFit)
      IF( ErrFlag ) THEN
         CALL statisticsErr(15)
         IF( truthmc ) THEN
            goto 66
         ELSE
            RETURN
         ENDIF
      ENDIF
      Interf%RcFit = SQRT( (Interf%KchBoost(7)-Interf%ipFit(1))**2 +
     &                     (Interf%KchBoost(8)-Interf%ipFit(2))**2 +
     &                     (Interf%KchBoost(9)-Interf%ipFit(3))**2 )
      Interf%RtcFit = SQRT( (Interf%KchBoost(7)-Interf%ipFit(1))**2 +
     &                      (Interf%KchBoost(8)-Interf%ipFit(2))**2 )
      Interf%RnFit = SQRT( (Interf%KneRecLorFit(7)-Interf%ipFit(1))**2 +
     &                     (Interf%KneRecLorFit(8)-Interf%ipFit(2))**2 +
     &                     (Interf%KneRecLorFit(9)-Interf%ipFit(3))**2 )
      Interf%RtnFit = SQRT( (Interf%KneRecLorFit(7)-Interf%ipFit(1))**2+
     &                      (Interf%KneRecLorFit(8)-Interf%ipFit(2))**2)
C
C------------------------------------------------------------------------------
C
C Perform fit for e+e- --> omegapi0 --> pi+pi-pi0pi0 background
C
      IF(Interf%omegatruth_top.EQ.1)THEN
         CALL fit_interf_omega(BkgFitPar,BkgErrPar,MaxNumComega,
     &        Chi2Fit,ConstrOmega,DOmega,MTXOmega,ZOmega,loopcount)
CC
CC Insert fit results into structure
CC
         DO i=1,MaxNumFitPar
            Interf%BkgFitPar(i) = BkgFitPar(i)
            Interf%BkgErrPar(i) = BkgErrPar(i)
         ENDDO
         Interf%Niter_w  = loopcount
         Interf%Chi2_w   = Chi2Fit
         Interf%Nconstr_w   = MaxNumComega
C
C------------------------------------------------------------------------------
C
C Track parameters
         DO ILoop = 1,2                     
            FitPar6(3*(ILoop-1)+1) = Evt%Trkv%px(Interf%chtrk_index(ILoo
     &p))
            FitPar6(3*(ILoop-1)+2) = Evt%Trkv%py(Interf%chtrk_index(ILoo
     &p))
            FitPar6(3*(ILoop-1)+3) = Evt%Trkv%pz(Interf%chtrk_index(ILoo
     &p))
         ENDDO
C
         CALL find_omegaminv(FitPar6,PgamRecTaken4x4fit,Interf%KchFit,
     &          Interf%ominvFit,Interf%PpioOmegaFit,Interf%P4PriRestFit,
     &          Interf%test(22))
         IF( ErrFlag ) THEN
            CALL statisticsErr(16)
            IF( truthmc ) THEN
               goto 66
            ELSE
               RETURN
            ENDIF
         ENDIF
      ELSE
         DO i=1,MaxNumFitPar
            Interf%BkgFitPar(i) = -999.
            Interf%BkgErrPar(i) = -999.
         ENDDO
         Interf%Niter_w  = -999
         Interf%Chi2_w   = -999.
         Interf%Nconstr_w   = -999
      ENDIF
C
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
C do cuts
C------------------------------------------------------------------------------
C
      CALL do_cuts
      IF( ErrFlag ) THEN
        CALL statisticsErr(17)
        IF( truthmc ) THEN
           goto 66
        ELSE
           RETURN
        ENDIF
      ENDIF
      IF( .not.makecuts ) EventOK = .TRUE.
C------------------------------------------------------------------------------
C Fill Ntuple
C------------------------------------------------------------------------------
C
 66   continue
      IF( EventOK.or.truthmc ) THEN
         selected=selected+1
C value of mctruth            signal       cuts passed    error
C                0              +            n/a          +
C                1              +            +            -
C                2              +            -            -
C                3            regen          +            -
C                4            omega          +            -
C                5            three          +            -
C                6            semi           +            -
C                7            else           +            -
         IF( truthmc ) THEN
            IF( ErrFlag  ) THEN
              Interf%mctruth = 0
              nmctruth0 = nmctruth0 + 1
            ELSEIF( EventOK ) THEN
              Interf%mctruth = 1
              nmctruth(1) = nmctruth(1) + 1
            ELSE
               Interf%mctruth = 2
               nmctruth(2) = nmctruth(2) + 1
            ENDIF
         ELSE IF( truthregg ) THEN
            Interf%mctruth = 3
            nmctruth(3) = nmctruth(3) + 1
         ELSE IF( truthomegaa ) THEN
            Interf%mctruth = 4
            nmctruth(4) = nmctruth(4) + 1
         ELSE IF( truththreee ) THEN
            Interf%mctruth = 5
            nmctruth(5) = nmctruth(5) + 1
         ELSE IF( truthsemii ) THEN
            Interf%mctruth = 6
            nmctruth(6) = nmctruth(6) + 1
         ELSE IF( truthelsee ) THEN
            Interf%mctruth = 7
            nmctruth(7) = nmctruth(7) + 1
         ENDIF
         Status = ANGOHS('INTERF',HisFlg,MinHisId,MaxHisId)
         IF( Status.NE.ANSUCC )THEN
            CALL ERLOGR('PM00EV',ERWARN,0,Status,
     &           'Error getting HBook info for INTERF module')
         ELSE
           IF( HisFlg )THEN
               ParSet = ANGPAR()
               IF( ParSet.EQ.1 )THEN
                  CALL HCDir('//PAWC/INTERF',' ')
                  CALL HFNT(1)
               ELSE
                  CALL ERLOGR('PM00EV',ERWARN,0,0,
     &                 'INTERF ntuple not active')
               ENDIF
            ENDIF
         ENDIF
C
         Status = ANPTRG(.TRUE.)
      ELSE
         Status = ANPTRG(.FALSE.)
      ENDIF
C
      END
C
C==============================================================================
C==============================================================================
C==============================================================================
      SUBROUTINE PM00RE
C------------------------------------------------------------------------------
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
C------------------------------------------------------------------------------
      END
C
C
C==============================================================================
C==============================================================================
C==============================================================================
      SUBROUTINE PM00FI
C------------------------------------------------------------------------------
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
C------------------------------------------------------------------------------
      WRITE(6,*)' '
      WRITE(6,*)'======================================================'
      WRITE(6,*)'Summary:'
      WRITE(6,*)'         cutchar =',cutchar,' cutneu =',cutneu,
     &                  ' onlytruemc =',onlytruemc,
     &                  ' makecuts =',makecuts
      CALL summary
      WRITE(6,*)'======================================================'
C
      END
C
C==============================================================================
C==============================================================================
C==============================================================================
      SUBROUTINE PM00HB
C------------------------------------------------------------------------------
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
C------------------------------------------------------------------------------
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
C------------------------------------------------------------------------------
C$$include 'c$inc:erlevl.inc'
C$$include 'a_c$library:anerror.inc'
C
C External functions
C
C     INTEGER    ANGOHS,UIGTLU
C
C Local declarations
C
      INTEGER    LunId
      CHARACTER*16 Range
C
C------------------------------------------------------------------------------
C
C
C Define ntuple id and ntuple name
      LunId=1
      CALL HBNt(LunId,'INTERF','D')
      CALL HCDir('//PAWC/INTERF',' ')
C-----------------------------------------------------------------------
C Event info
C-----------------------------------------------------------------------
C  first 30330 and last 41902 run nr for 2004-2005 at phi mass
C     CALL find_range(30300,41902,Range)
      CALL HBName(LunId,'INFO',Interf%nev,'nev:i')
      CALL HBName(LunId,'INFO',Interf%pileup,'pileup[0,1]:i')
      CALL find_range(0,100,Range)
      CALL HBName(LunId,'INFO',Interf%gcod,'gcod:i::'//Range)
      CALL HBName(LunId,'INFO',Interf%phid,'phid:i::'//Range)
      CALL HBName(LunId,'INFO',Interf%a1typ,'a1typ:i::'//Range)
      CALL HBName(LunId,'INFO',Interf%a2typ,'a2typ:i::'//Range)
      CALL HBName(LunId,'INFO',Interf%a3typ,'a3typ:i::'//Range)
      CALL HBName(LunId,'INFO',Interf%b1typ,'b1typ:i::'//Range)
      CALL HBName(LunId,'INFO',Interf%b2typ,'b2typ:i::'//Range)
      CALL HBName(LunId,'INFO',Interf%b3typ,'b3typ:i::'//Range)
C-----------------------------------------------------------
C *** T0 Finder
C-----------------------------------------------------------
C     CALL HBName(LunId,'INFO',Interf%tphased_mc,'tphased_mc:r')
C     CALL HBName(LunId,'INFO',Interf%t0dc0,'t0dc0:r')
C     CALL HBName(LunId,'INFO',Interf%t0hit0,'t0hit0:r')
C     CALL HBName(LunId,'INFO',Interf%t0clu0,'t0clu0:r')
      CALL HBName(LunId,'INFO',Interf%T0step1,'T0step1:r')
C     CALL HBName(LunId,'INFO',Interf%DelayCable,'DelayCable:r')
C     CALL HBName(LunId,'INFO',Interf%Tbunch,'Tbunch:r')
C------------------------------------------------------------------------------
C     CALL HBName(LunId,'DATA',Interf%TimeSec,'TimeSec:i')
C     CALL HBName(LunId,'DATA',Interf%TimeMusec,'TimeMusec:i')
      CALL HBName(LunId,'DATA',mcflag,'mcflag:i::[0,1]')
C------------------------------------------------------------------------------
      CALL HBName(LunId,'BPOS',Interf%Bx,'Bx:r')
      CALL HBName(LunId,'BPOS',Interf%By,'By:r')
      CALL HBName(LunId,'BPOS',Interf%Bz,'Bz:r')
      CALL HBName(LunId,'BPOS',Interf%Bsx,'Bsx:r')
      CALL HBName(LunId,'BPOS',Interf%Bsy,'Bsy:r')
      CALL HBName(LunId,'BPOS',Interf%Bsz,'Bsz:r')
C------------------------------------------------------------------------------
      CALL HBName(LunId,'BMOM',Interf%Bpx,'Bpx:r')
      CALL HBName(LunId,'BMOM',Interf%Bpy,'Bpy:r')
      CALL HBName(LunId,'BMOM',Interf%Bpz,'Bpz:r')
      CALL HBName(LunId,'BMOM',Interf%Bwidpx,'Bwidpx:r')
      CALL HBName(LunId,'BMOM',Interf%Bwidpy,'Bwidpy:r')
      CALL HBName(LunId,'BMOM',Interf%Bwidpz,'Bwidpz:r')
      CALL HBName(LunId,'BMOM',Interf%Blumx,'Blumx:r')
      CALL HBName(LunId,'BMOM',Interf%Blumz,'Blumz:r')
      CALL HBName(LunId,'BMOM',Interf%Broots,'Broots:r')
      CALL HBName(LunId,'BMOM',Interf%BrootsErr,'BrootsErr:r')
C------------------------------------------------------------------------------
      CALL HBName(LunId,'ECLS2',Interf%necls2,'necls2:i::[0,8]')
      CALL HBName(LunId,'ECLS2',Interf%ECLtrgw2,'ECLtrgw2:i')
      CALL HBName(LunId,'ECLS2',Interf%ECLfilfo2,'ECLfilfo2:i')
      CALL HBName(LunId,'ECLS2',Interf%ECLword2,'ECLword2(necls2):i')
      CALL HBName(LunId,'ECLS2',Interf%ECLstream2,
     &                              'ECLstream2(necls2):i')
      CALL HBName(LunId,'ECLS2',Interf%ECLtagnum2,
     &                              'ECLtagnum2(necls2):i')
      CALL HBName(LunId,'ECLS2',Interf%ECLevtype2,
     &                              'ECLevtype2(necls2):i')
C-----------------------------------------------------------------------
C EMC Clusters
C-----------------------------------------------------------------------
      CALL find_range(0,MaxNumClu,Range)
      CALL HBName(LunId,'CLU',Interf%nclu,'nclu:i::'//Range)
      CALL HBName(LunId,'CLU',Interf%EneCl,'EneCl(nclu):r')
      CALL HBName(LunId,'CLU',Interf%Tcl,'Tcl(nclu):r')
      CALL HBName(LunId,'CLU',Interf%Xcl,'Xcl(nclu):r')
      CALL HBName(LunId,'CLU',Interf%Ycl,'Ycl(nclu):r')
      CALL HBName(LunId,'CLU',Interf%Zcl,'Zcl(nclu):r')
      CALL HBName(LunId,'CLU',Interf%Xacl,'Xacl(nclu):r')
      CALL HBName(LunId,'CLU',Interf%Yacl,'Yacl(nclu):r')
      CALL HBName(LunId,'CLU',Interf%Zacl,'Zacl(nclu):r')
      CALL HBName(LunId,'CLU',Interf%XRmCl,'XRmCl(nclu):r')
      CALL HBName(LunId,'CLU',Interf%YRmsCl,'YRmsCl(nclu):r')
      CALL HBName(LunId,'CLU',Interf%ZrmsCl,'ZrmsCl(nclu):r')
      CALL HBName(LunId,'CLU',Interf%TrmsCl,'TrmsCl(nclu):r')
      CALL HBName(LunId,'CLU',Interf%FlagCl,'FlagCl(nclu):i')
C------------------------------------------------------------------------------
C     CALL HBName(LunId,'CLUMC',Interf%nclumc,'nclumc:i::'//Range)
C     CALL HBName(LunId,'CLUMC',Interf%Npar,'Npar(nclumc)[0,10]:i')
C     CALL find_range(0,100,Range)
C     CALL HBName(LunId,'CLUMC',Interf%Pnum1,'Pnum1(nclumc):i::'//Range)
C     CALL HBName(LunId,'CLUMC',Interf%Pid1,'Pid1(nclumc):i::'//Range)
C     CALL HBName(LunId,'CLUMC',Interf%Pnum2,'Pnum2(nclumc):i::'//Range)
C     CALL HBName(LunId,'CLUMC',Interf%Pid2,'Pid2(nclumc):i::'//Range)
C     CALL HBName(LunId,'CLUMC',Interf%Pnum3,'Pnum3(nclumc):i::'//Range)
C     CALL HBName(LunId,'CLUMC',Interf%Pid3,'Pid3(nclumc):i::'//Range)
C-----------------------------------------------------------------------
C Tracks after Vertex reconstruction
C-----------------------------------------------------------------------
      CALL find_range(0,MaxNumTrkV,Range)
      CALL HBName(LunId,'TRKV',Interf%ntv,'ntv:i::'//Range)
      CALL find_range(0,MaxNumVtx,Range)
      CALL HBName(LunId,'TRKV',Interf%iv,'iv(ntv):i::'//Range)
      CALL HBName(LunId,'TRKV',Interf%trknumv,'trknumv(ntv):i::[0,999]')
      CALL HBName(LunId,'TRKV',Interf%CurV,'CurV(ntv):r')
      CALL HBName(LunId,'TRKV',Interf%PhiV,'PhiV(ntv):r')
      CALL HBName(LunId,'TRKV',Interf%CotV,'CotV(ntv):r')
C-----------------------------------------------------------------------
C Verticies
C-----------------------------------------------------------------------
      CALL find_range(0,MaxNumVtx,Range)
      CALL HBName(LunId,'VTX',Interf%nv,'nv:i::'//Range)
      CALL HBName(LunId,'VTX',Interf%vtx,'vtx(nv):i::[0,10]')
      CALL HBName(LunId,'VTX',Interf%xv,'xv(nv):r')
      CALL HBName(LunId,'VTX',Interf%yv,'yv(nv):r')
      CALL HBName(LunId,'VTX',Interf%zv,'zv(nv):r')
      CALL HBName(LunId,'VTX',Interf%chivtx,'chivtx(nv):r')
      CALL HBName(LunId,'VTX',Interf%qualv,'qualv(nv):i')
      CALL HBName(LunId,'VTX',Interf%fitidv,'fitidv(nv):i')
C     CALL HBName(LunId,'VTX',Interf%VTXcov1,'VTXcov1(nv):r')
C     CALL HBName(LunId,'VTX',Interf%VTXcov2,'VTXcov2(nv):r')
C     CALL HBName(LunId,'VTX',Interf%VTXcov3,'VTXcov3(nv):r')
C     CALL HBName(LunId,'VTX',Interf%VTXcov4,'VTXcov4(nv):r')
C     CALL HBName(LunId,'VTX',Interf%VTXcov5,'VTXcov5(nv):r')
C     CALL HBName(LunId,'VTX',Interf%VTXcov6,'VTXcov6(nv):r')
C-----------------------------------------------------------------------
C Trks before Vertex reconstruction
C-----------------------------------------------------------------------
      CALL find_range(0,MaxNumTrk,Range)
      CALL HBName(LunId,'TRK',Interf%nt,'nt:i::'//Range)
C     CALL HBName(LunId,'TRK',Interf%trkind,'trkind(nt):i')
      CALL HBName(LunId,'TRK',Interf%chi2fit,'chi2fit(nt):r')
      CALL HBName(LunId,'TRK',Interf%chi2ms,'chi2ms(nt):r')
C------------------------------------------------------------------------------
C     CALL HBName(LunId,'TRKMC',Interf%ntfmc,'ntfmc:i::'//Range)
C     CALL HBName(LunId,'TRKMC',Interf%trkine1,'trkine1(ntfmc):i')
C     CALL HBName(LunId,'TRKMC',Interf%trtype1,'trtype1(ntfmc):i')
C     CALL HBName(LunId,'TRKMC',Interf%trhits1,'trhits1(ntfmc):i')
C     CALL HBName(LunId,'TRKMC',Interf%trkine2,'trkine2(ntfmc):i')
C     CALL HBName(LunId,'TRKMC',Interf%trtype2,'trtype2(ntfmc):i')
C     CALL HBName(LunId,'TRKMC',Interf%trhits2,'trhits2(ntfmc):i')
C     CALL HBName(LunId,'TRKMC',Interf%trkine3,'trkine3(ntfmc):i')
C     CALL HBName(LunId,'TRKMC',Interf%trtype3,'trtype3(ntfmc):i')
C     CALL HBName(LunId,'TRKMC',Interf%trhits3,'trhits3(ntfmc):i')
C-----------------------------------------------------------------------
C GEANFI MC Information
C-----------------------------------------------------------------------
      CALL find_range(0,MaxNtrkGen,Range)
      CALL HBName(LunId,'MC',Interf%ntmc,'ntmc:i::'//Range)
      CALL find_range(0,100,Range)
      CALL HBName(LunId,'MC',Interf%kine,'kine(ntmc):i::'//Range)
      CALL HBName(LunId,'MC',Interf%pidmc,'pidmc(ntmc):i::'//Range)
      CALL HBName(LunId,'MC',Interf%virmom,'virmom(ntmc):i::'//Range)
      CALL HBName(LunId,'MC',Interf%pxmc,'pxmc(ntmc):r')
      CALL HBName(LunId,'MC',Interf%pymc,'pymc(ntmc):r')
      CALL HBName(LunId,'MC',Interf%pzmc,'pzmc(ntmc):r')
      CALL HBName(LunId,'MC',Interf%themc,'themc(ntmc):r')
      CALL HBName(LunId,'MC',Interf%phimc,'phimc(ntmc):r')
      CALL HBName(LunId,'MC',Interf%vtxmc,'vtxmc(ntmc):i::'//Range)
      CALL find_range(0,MaxNvtxGen,Range)
      CALL HBName(LunId,'MC',Interf%nvtxmc,'nvtxmc:i::'//Range)
      CALL find_range(0,100,Range)
      CALL HBName(LunId,'MC',Interf%kinmom,
     &                             'kinmom(nvtxmc):i::'//Range)
      CALL HBName(LunId,'MC',Interf%mother,
     &                             'mother(nvtxmc):i::'//Range)
      CALL HBName(LunId,'MC',Interf%xvmc,'xvmc(nvtxmc):r')
      CALL HBName(LunId,'MC',Interf%yvmc,'yvmc(nvtxmc):r')
      CALL HBName(LunId,'MC',Interf%zvmc,'zvmc(nvtxmc):r')
      CALL HBName(LunId,'MC',Interf%ntvtx,'ntvtx(nvtxmc):r')
C-----------------------------------------------------------------------
C Tracks vs cluster merging ... (BANK TCLO)
C-----------------------------------------------------------------------
C     CALL find_range(0,NTCLOmax,Range)
C     CALL HBName(LunId,'TCLO',Interf%ntcl,'ntcl:i::'//Range)
C     CALL HBName(LunId,'TCLO',Interf%Asstr,'Asstr(ntcl):i::[0,999]')
C     CALL HBName(LunId,'TCLO',Interf%Asscl,'Asscl(ntcl):i::[0,100]')
C     CALL HBName(LunId,'TCLO',Interf%verver,
C    &               'verver(ntcl):i::[-100,100]')
C     CALL HBName(LunId,'TCLO',Interf%xext,'xext(ntcl):r')
C     CALL HBName(LunId,'TCLO',Interf%yext,'yext(ntcl):r')
C     CALL HBName(LunId,'TCLO',Interf%zext,'zext(ntcl):r')
C     CALL HBName(LunId,'TCLO',Interf%Assleng,'Assleng(ntcl):r')
C     CALL HBName(LunId,'TCLO',Interf%AssChi,'AssChi(ntcl):r')
C     CALL HBName(LunId,'TCLO',Interf%extPx,'extPx(ntcl):r')
C     CALL HBName(LunId,'TCLO',Interf%extPy,'extPy(ntcl):r')
C     CALL HBName(LunId,'TCLO',Interf%extPz,'extPz(ntcl):r')
C-----------------------------------------------------------------------
C     ECLO  - informatins from ECLO bank
C-----------------------------------------------------------------------
C     CALL find_range(0,MaxNumCLINF,Range)
C     CALL HBName(LunId,'ECLO2',Interf%ncli2,'ncli2:i::'//Range)
C     CALL HBName(LunId,'ECLO2',Interf%ECLOword2,'ECLOword2(ncli2):i')
C     CALL HBName(LunId,'ECLO2',Interf%idpart2,'idpart2(ncli2):i')
C     CALL HBName(LunId,'ECLO2',Interf%dtclpo2,'dtclpo2(ncli2):i')
C     CALL HBName(LunId,'ECLO2',Interf%dvvnpo2,'dvvnpo2(ncli2):i')
C     CALL HBName(LunId,'ECLO2',Interf%stre2,'stre2(ncli2):i')
C     CALL HBName(LunId,'ECLO2',Interf%algo2,'algo2(ncli2):i')
C-----------------------------------------------------------------------
C     INTERF
C-----------------------------------------------------------------------
      CALL HBName(LunId,'INTERF',Interf%mctruth,'mctruth:i::[0,7]')
      CALL find_range(0,testcounter,Range)
      CALL HBName(LunId,'INTERF',Interf%ErrId,'ErrId:i::'//Range)
      CALL HBName(LunId,'INTERF',Interf%CutId,'CutId:i::'//Range)
      CALL HBName(LunId,'INTERF',Interf%KchMC,'KchMC(9):r')
      CALL HBName(LunId,'INTERF',Interf%KneMC,'KneMC(9):r')
      CALL HBName(LunId,'INTERF',Interf%ipmc,'ipmc(3):r')
      CALL HBName(LunId,'INTERF',Interf%KchRec,'KchRec(9):r')
      CALL HBName(LunId,'INTERF',Interf%KchBoost,'KchBoost(9):r')
      CALL HBName(LunId,'INTERF',Interf%KneRec,'KneRec(9):r')
      CALL HBName(LunId,'INTERF',Interf%KneRecLor,'KneRecLor(9):r')
      CALL HBName(LunId,'INTERF',Interf%ip,'ip(3):r')
      CALL HBName(LunId,'INTERF',Interf%chdist,'chdist:r:[0.,]')
      CALL HBName(LunId,'INTERF',Interf%cldist,'cldist:r:[0.,]')
      CALL HBName(LunId,'INTERF',Interf%mcISR,'mcISR:i::[0,1]')
      CALL find_range(0,MaxNumVtx,Range)
      CALL HBName(LunId,'INTERF',Interf%vtaken,'vtaken(3):i::'//Range)
      CALL find_range(0,MaxNumCLu,Range)
      CALL HBName(LunId,'INTERF',Interf%ncl,'ncl:i::'//Range)
      CALL HBName(LunId,'INTERF',Interf%ncll,'ncll(ncl):i::'//Range)
      CALL HBName(LunId,'INTERF',Interf%nclwrong,'nclwrong:i::'//Range)
      CALL HBName(LunId,'INTERF',Interf%ncllwrong,
     &                           'ncllwrong(nclwrong):i::'//Range)
      CALL HBName(LunId,'INTERF',Interf%DlBoostRec,'DlBoostRec:r')
      CALL HBName(LunId,'INTERF',Interf%DtBoostRec,'DtBoostRec:r')
      CALL HBName(LunId,'INTERF',Interf%DlBoostLor,'DlBoostLor:r')
      CALL HBName(LunId,'INTERF',Interf%DtBoostLor,'DtBoostLor:r')
      CALL HBName(LunId,'INTERF',Interf%DlRec,'DlRec:r')
      CALL HBName(LunId,'INTERF',Interf%DtRec,'DtRec:r')
      CALL HBName(LunId,'INTERF',Interf%DlMC,'DlMC:r')
      CALL HBName(LunId,'INTERF',Interf%DtMC,'DtMC:r')
      CALL HBName(LunId,'INTERF',Interf%Rc,'Rc:r:[0.,]')
      CALL HBName(LunId,'INTERF',Interf%Rtc,'Rtc:r:[0.,]')
      CALL HBName(LunId,'INTERF',Interf%Rn,'Rn:r:[0.,]')
      CALL HBName(LunId,'INTERF',Interf%Rtn,'Rtn:r:[0.,]')
      CALL HBName(LunId,'INTERF',Interf%RcMC,'RcMC:r:[0.,]')
      CALL HBName(LunId,'INTERF',Interf%RtcMC,'RtcMC:r:[0.,]')
      CALL HBName(LunId,'INTERF',Interf%RnMC,'RnMC:r:[0.,]')
      CALL HBName(LunId,'INTERF',Interf%RtnMC,'RtnMC:r:[0.,]')
      CALL HBName(LunId,'INTERF',Interf%trc,'trc(ncl):r:[0.,]')
      CALL HBName(LunId,'INTERF',Interf%trcv,'trcv(ncl):r:[0.,]')
      CALL HBName(LunId,'INTERF',Interf%pi0,'pi0(2):r:[0.,]')
      CALL HBName(LunId,'INTERF',Interf%minv4gam,'minv4gam:r:[0.,]')
      CALL HBName(LunId,'INTERF',Interf%PgamRec1,'PgamRec1(4):r')
      CALL HBName(LunId,'INTERF',Interf%PgamRec2,'PgamRec2(4):r')
      CALL HBName(LunId,'INTERF',Interf%PgamRec3,'PgamRec3(4):r')
      CALL HBName(LunId,'INTERF',Interf%PgamRec4,'PgamRec4(4):r')
      CALL HBName(LunId,'INTERF',Interf%gpairtaken,
     &                                  'gpairtaken(2):i:[0,]')
      CALL HBName(LunId,'INTERF',Interf%g4vtxerr,'g4vtxerr(3):r')
      CALL HBName(LunId,'INTERF',Interf%ominv,'ominv(2):r:[0.,]')
      CALL HBName(LunId,'INTERF',Interf%PpioOmega,'PpioOmega(4):r')
      CALL HBName(LunId,'INTERF',Interf%P4PriRest,'P4PriRest(4):r')
      CALL HBName(LunId,'INTERF',Interf%trk1,'trk1(4):r')
      CALL HBName(LunId,'INTERF',Interf%trk2,'trk2(4):r')
      CALL HBName(LunId,'INTERF',Interf%cosTrk,'cosTrk:r:[-1.,1.]')
      CALL HBName(LunId,'INTERF',Interf%Qmiss,'Qmiss:r')
      CALL find_range(0,MaxNumCLu,Range)
      CALL HBName(LunId,'INTERF',Interf%g4taken,'g4taken(4):i::'//Range)
      CALL HBName(LunId,'KINFIT',Interf%nparfit,'nparfit:i')
      CALL HBName(LunId,'KINFIT',Interf%KchFit,'KchFit(9):r')
      CALL HBName(LunId,'KINFIT',Interf%KneFit,'KneFit(9):r')
      CALL HBName(LunId,'KINFIT',Interf%KneRecLorFit,
     &                                 'KneRecLorFit(9):r')
      CALL HBName(LunId,'KINFIT',Interf%chdistFit,'chdistFit:r')
      CALL HBName(LunId,'KINFIT',Interf%cldistFit,
     &                                 'cldistFit:r:[0.,]')
      CALL HBName(LunId,'KINFIT',Interf%ipFit,'ipFit(3):r')
      CALL HBName(LunId,'KINFIT',Interf%nclwrongFit,'nclwrongFit:i::'//R
     &ange)
      CALL HBName(LunId,'KINFIT',Interf%ncllwrongFit,
     &                           'ncllwrongFit(nclwrongFit):i::'//Range)
      CALL HBName(LunId,'KINFIT',Interf%DlFit,'DlFit:r')
      CALL HBName(LunId,'KINFIT',Interf%DtFit,'DtFit:r')
      CALL HBName(LunId,'KINFIT',Interf%RcFit,'RcFit:r')
      CALL HBName(LunId,'KINFIT',Interf%RtcFit,'RtcFit:r')
      CALL HBName(LunId,'KINFIT',Interf%RnFit,'RnFit:r')
      CALL HBName(LunId,'KINFIT',Interf%RtnFit,'RtnFit:r')
      CALL HBName(LunId,'KINFIT',Interf%ncl,'nclfit:i::'//Range)
      CALL HBName(LunId,'KINFIT',Interf%trcFit,'trcFit(nclfit):r')
      CALL HBName(LunId,'KINFIT',Interf%trcvFit,'trcvFit(nclfit):r')
      CALL HBName(LunId,'KINFIT',Interf%pi0Fit,'pi0Fit(2):r:[0.,]')
      CALL HBName(LunId,'KINFIT',Interf%PgamRec1fit,'PgamRec1fit(4):r')
      CALL HBName(LunId,'KINFIT',Interf%PgamRec2fit,'PgamRec2fit(4):r')
      CALL HBName(LunId,'KINFIT',Interf%PgamRec3fit,'PgamRec3fit(4):r')
      CALL HBName(LunId,'KINFIT',Interf%PgamRec4fit,'PgamRec4fit(4):r')
      CALL HBName(LunId,'KINFIT',Interf%gpairtakenFit,
     &                                 'gpairtakenFit(2):i:[0,]')
      CALL HBName(LunId,'KINFIT',Interf%g4vtxerrFit,'g4vtxerrFit(3):r')
      CALL HBName(LunId,'KINFIT',Interf%minv4gamFit,
     &                                 'minv4gamFit:r:[0.,]')
      CALL HBName(LunId,'KINFIT',Interf%ominvFit,'ominvFit(2):r:[0.,]')
      CALL HBName(LunId,'KINFIT',Interf%PpioOmegaFit,
     &                                 'PpioOmegaFit(4):r')
      CALL HBName(LunId,'KINFIT',Interf%P4PriRestFit,
     &                                  'P4PriRestFit(4):r')
      CALL HBName(LunId,'KINFIT',Interf%g4takenFit,
     &                                  'g4takenFit(4):i::'//Range)
      CALL HBName(LunId,'KINFIT',Interf%Niter,'Niter:i')
      CALL HBName(LunId,'KINFIT',Interf%Niter_w,'Niterw:i')
      CALL HBName(LunId,'KINFIT',Interf%Chi2,'Chi2:r')
      CALL HBName(LunId,'KINFIT',Interf%Chi2_w,'Chi2w:r:')
      CALL HBName(LunId,'KINFIT',Interf%FitPar,'FitPar(38):r')
      CALL HBName(LunId,'KINFIT',Interf%ErrPar,'ErrPar(38):r')
      CALL HBName(LunId,'KINFIT',Interf%BkgFitPar,'BkgFitPar(38):r')
      CALL HBName(LunId,'KINFIT',Interf%BkgErrPar,'BkgErrPar(38):r')
      CALL HBName(LunId,'KINFIT',Interf%FitParStart,'FitParStart(38):r')
      CALL HBName(LunId,'KINFIT',Interf%ErrParStart,'ErrParStart(38):r')
      CALL HBName(LunId,'KINFIT',Interf%Nconstr,'Nconstr:i')
      CALL HBName(LunId,'KINFIT',Interf%Nconstr_w,'Nconstr_w:i')
      CALL HBName(LunId,'SIMG',Interf%simok,'simok:i:[0,1]')
      CALL find_range(0,MaxNumVtx,Range)
      CALL HBName(LunId,'SIMG',Interf%ChVtxId,'ChVtxId:i::'//Range)
      CALL find_range(0,MaxNumTrk,Range)
      CALL HBName(LunId,'SIMG',Interf%TrkIdx,'TrkIdx(2):i::'//Range)
      CALL find_range(0,MaxNumCLu,Range)
      CALL HBName(LunId,'SIMG',Interf%CluIdx,'CluIdx(4):i::'//Range)
      CALL HBName(LunId,'SIMG',Interf%Trkk1,'Trkk1(3):r')
      CALL HBName(LunId,'SIMG',Interf%Trkk2,'Trkk2(3):r')
      CALL HBName(LunId,'SIMG',Interf%ChaVtx,'ChaVtx(3):r')
      CALL HBName(LunId,'SIMG',Interf%NeuVtx,'NeuVtx(3):r')
      CALL HBName(LunId,'SIMG',Interf%PhiVtx,'PhiVtx(3):r')
      CALL HBName(LunId,'TEST',Interf%test,'test(100):r')
      CALL HBName(LunId,'INTERF',Interf%omegatruth_top,'omegatruth_top:i
     &')
      CALL HBName(LunId,'INTERF',Interf%chtrk_index,'chtrk_index(2):i')
      CALL HBName(LunId,'INTERF',Interf%clu_index,'clu_index(4):i')
      CALL HBName(LunId,'INTERF',Interf%chvtx_index,'chvtx_index:i')
      END
C
C==============================================================================
C==============================================================================
C==============================================================================
      SUBROUTINE PM00TK
C------------------------------------------------------------------------------
C
C  Description:
C  ------------
C
C------------------------------------------------------------------------------
C
        IMPLICIT NONE
C== Include /kloe/soft/off/uipack/production/aix/library/uierror.cin ===
      INTEGER    UISUCC
      PARAMETER (UISUCC = '8128009'X)
      INTEGER    UIDFLT
      PARAMETER (UIDFLT = '8128011'X)
      INTEGER    UITRUN
      PARAMETER (UITRUN = '8128019'X)
      INTEGER    UIQUPR
      PARAMETER (UIQUPR = '8128021'X)
      INTEGER    UIQUNE
      PARAMETER (UIQUNE = '8128029'X)
      INTEGER    UIQUVA
      PARAMETER (UIQUVA = '8128031'X)
      INTEGER    UIQUDF
      PARAMETER (UIQUDF = '8128039'X)
      INTEGER    UINOIN
      PARAMETER (UINOIN = '8128083'X)
      INTEGER    UIABRT
      PARAMETER (UIABRT = '812808B'X)
      INTEGER    UICANC
      PARAMETER (UICANC = '8128093'X)
      INTEGER    UIMPTY
      PARAMETER (UIMPTY = '81280C0'X)
      INTEGER    UIQUAB
      PARAMETER (UIQUAB = '81280C8'X)
      INTEGER    UIEOF
      PARAMETER (UIEOF  = '81280D0'X)
      INTEGER    UININI
      PARAMETER (UININI = '8128102'X)
      INTEGER    UIOVRF
      PARAMETER (UIOVRF = '812810A'X)
      INTEGER    UINOGP
      PARAMETER (UINOGP = '8128112'X)
      INTEGER    UINOME
      PARAMETER (UINOME = '812811A'X)
      INTEGER    UINOLU
      PARAMETER (UINOLU = '8128122'X)
      INTEGER    UINOPC
      PARAMETER (UINOPC = '812824A'X)
      INTEGER    UINOVP
      PARAMETER (UINOVP = '81281DA'X)
      INTEGER    UIILBU
      PARAMETER (UIILBU = '812812A'X)
      INTEGER    UIILCO
      PARAMETER (UIILCO = '81281E2'X)
      INTEGER    UIILDV
      PARAMETER (UIILDV = '8128132'X)
      INTEGER    UIILGP
      PARAMETER (UIILGP = '812813A'X)
      INTEGER    UIILLI
      PARAMETER (UIILLI = '812818A'X)
      INTEGER    UIILME
      PARAMETER (UIILME = '8128142'X)
      INTEGER    UIILMD
      PARAMETER (UIILMD = '812814A'X)
      INTEGER    UIILPA
      PARAMETER (UIILPA = '81281EA'X)
      INTEGER    UIILPC
      PARAMETER (UIILPC = '812823A'X)
      INTEGER    UIILQU
      PARAMETER (UIILQU = '8128152'X)
      INTEGER    UIILSZ
      PARAMETER (UIILSZ = '812815A'X)
      INTEGER    UIILSB
      PARAMETER (UIILSB = '812821A'X)
      INTEGER    UIILSY
      PARAMETER (UIILSY = '8128162'X)
      INTEGER    UIILTP
      PARAMETER (UIILTP = '812816A'X)
      INTEGER    UIILUI
      PARAMETER (UIILUI = '8128262'X)
      INTEGER    UIILLU
      PARAMETER (UIILLU = '8128172'X)
      INTEGER    UIILVB
      PARAMETER (UIILVB = '812817A'X)
      INTEGER    UIILVP
      PARAMETER (UIILVP = '8128182'X)
      INTEGER    UIDUGP
      PARAMETER (UIDUGP = '8128192'X)
      INTEGER    UIDUME
      PARAMETER (UIDUME = '812819A'X)
      INTEGER    UIDUPA
      PARAMETER (UIDUPA = '81281F2'X)
      INTEGER    UIDUPC
      PARAMETER (UIDUPC = '8128242'X)
      INTEGER    UIDUQU
      PARAMETER (UIDUQU = '81281A2'X)
      INTEGER    UIDUSB
      PARAMETER (UIDUSB = '8128222'X)
      INTEGER    UIDUVB
      PARAMETER (UIDUVB = '81281AA'X)
      INTEGER    UIDUVP
      PARAMETER (UIDUVP = '81281FA'X)
      INTEGER    UIRCSB
      PARAMETER (UIRCSB = '0812822A'X)
      INTEGER    UIRELU
      PARAMETER (UIRELU = '81281B2'X)
      INTEGER    UISTAK
      PARAMETER (UISTAK = '81281BA'X)
      INTEGER    UICIFI
      PARAMETER (UICIFI = '81281C2'X)
      INTEGER    UIRIFI
      PARAMETER (UIRIFI = '0812820A'X)
      INTEGER    UICOFI
      PARAMETER (UICOFI = '81281CA'X)
      INTEGER    UICOWR
      PARAMETER (UICOWR = '8128212'X)
      INTEGER    UIGLOK
      PARAMETER (UIGLOK = '81281D2'X)
      INTEGER    UISLOK
      PARAMETER (UISLOK = '8128232'X)
      INTEGER    UIVLOK
      PARAMETER (UIVLOK = '8128202'X)
      INTEGER    UIEMST
      PARAMETER (UIEMST = '8128252'X)
      INTEGER    UINTME
      PARAMETER (UINTME = '812825A'X)
C=== Include /kloe/soft/off/s_i/production/aix/library/noarginc.cin ====
      EXTERNAL   N$A
      INTEGER    SIMSAR
      PARAMETER (SIMSAR = 0)            
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
C- Local variables
      integer status
      integer GROUP, MENUID1, MENUID2, MENUID3
      character*40 VERB
C
C- External functions
      integer UIUSGP, UIACME, UIGTIN, UIGTRE, UIGTYE, UIGTFI, UIDFFI
      integer SITRFN
C
      status = UIDFFI('/gpfs/user/g/gamrat/analysis/prod/src/klspm00.uid
     &',
     $                 GROUP,MENUID1,N$A,N$A,MENUID3,N$A)
      if (MENUID3.GT.MENUID1) MENUID2=MENUID1+1
 11   continue
      Status = UIUSGP(GROUP,N$A)
      Status = UIACME(MENUID1,VERB,N$A)
      if (Status.EQ.UIABRT) goto 12
C-"Capitalize words" is mandatory
      if (Verb.eq.'CUTCHAR') then
        STATUS = UIGTye("Make cut on charged? Yes/No",cutchar)
      elseif (Verb.eq.'CUTNEU') then
        STATUS = UIGTye("Make cut on neutral? Yes/No",cutneu)
      elseif (Verb.eq.'ONLYTRUEMC') then
        STATUS = UIGTye("Only true MC signal events? Yes/No",onlytruemc)
      elseif (Verb.eq.'MAKECUTS') then
        STATUS = UIGTye("Do cuts? Yes/No",makecuts)
      elseif (Verb.eq.'RETURN') then
         goto 12
      endif
      goto 11
C
 12   continue
C
C
C
C External functions
C
C
C Local declarations
C
C
C------------------------------------------------------------------------------
C
C
      END
