C=======================================================================
C=======================================================================
      SUBROUTINE clearstruct

      USE ANALYSISMODULE
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
      Interf%vtakenks=0
      Interf%vtakenkl=0
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
      Interf%cldist=0.
      Interf%minv4gam=0.
      Interf%Rc=0.
      Interf%Rtc=0.
      Interf%Rn=0.
      Interf%Rtn=0.
      Interf%RcMC=0.
      Interf%RtcMC=0.
      Interf%RnMC=0.
      Interf%RtnMC=0.
      Interf%chdist=0.
      Interf%Chi2=0.
      Interf%Chi2_w=0.
      Interf%cosTrk=0.
      Interf%cosTrkKS=0.
      Interf%cosTrkKL=0.
      Interf%cosTrkCM=0.
      Interf%Qmiss=0.
C loops
      DO i = 1,MaxNumOverlapStream
          Interf%ECLword2(i)=0
          Interf%ECLstream2(i)=0
          Interf%ECLtagnum2(i)=0
          Interf%ECLevtype2(i)=0
      ENDDO
      DO i = 1,MaxNumClu
          Interf%EneCl(i)=0.
          Interf%TclOld(i)=0.
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
          Interf%trcv(i)=0.
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
      ENDDO
      DO i = 1,MaxNumtrkv
          Interf%ivOld(i)=0
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
          Interf%xvOld(i)=0.
          Interf%yvOld(i)=0.
          Interf%zvOld(i)=0.
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
          Interf%pidmcOld(i)=0
          Interf%virmom(i)=0
          Interf%vtxmcOld(i)=0
          Interf%kinmom(i)=0
          Interf%motherOld(i)=0
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
C DO i = 1,MaxNumCLINF
C     Interf%ECLOword2(i)=0
C     Interf%idpart2(i)=0
C     Interf%dtclpo2(i)=0
C     Interf%dvvnpo2(i)=0
C     Interf%stre2(i)=0
C     Interf%algo2(i)=0
C ENDDO
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
          Interf%KchRecKS(i)=0.
          Interf%KchRecKL(i)=0.
          Interf%KneRec(i)=0.
          Interf%KchBoost(i)=0.
          Interf%KneRecLor(i)=0.
          Interf%KneMC(i)=0.
      ENDDO
      DO i = 1,3
          Interf%ip(i)=0.
          Interf%ip_closest(i)=0.
          Interf%ip_plane(i)=0.
          Interf%ipmc(i)=0.
          Interf%g4vtxerr(i)=0.
          Interf%Trkk1(i)=0.
          Interf%Trkk2(i)=0.
          Interf%ChaVtx(i)=0.
          Interf%NeuVtx(i)=0.
          Interf%PhiVtx(i)=0.
      ENDDO
      DO i = 1,2
          Interf%pi0(i)=0.
          Interf%ominv(i)=0.
          Interf%gpairtaken(i)=0
          Interf%TrkIdx(i)=0
      ENDDO
      DO i = 1,4
          Interf%PpioOmega(i)=0.
          Interf%P4PriRest(i)=0.
          Interf%trk1(i)=0.
          Interf%trk2(i)=0.
          Interf%trk1KS(i)=0.
          Interf%trk2KS(i)=0.
          Interf%trk1KL(i)=0.
          Interf%trk2KL(i)=0.
          Interf%g4taken(i)=0
          Interf%PgamRec1(i)=0.
          Interf%PgamRec2(i)=0.
          Interf%PgamRec3(i)=0.
          Interf%PgamRec4(i)=0.
          Interf%CluIdx(i)=0
      ENDDO
      DO i = 1,testcounter
          Interf%test(i)=0.
      ENDDO
      END

C======================================================================
C======================================================================
C======================================================================

      SUBROUTINE GetKslEvent(ntmc,motherOld,vtxmcOld,pidmcOld,xvmc,yvmc,zvmc,
     &           pxmc,pymc,pzmc,nvtxmc,ipmc,KchMC,KneMC,DtMC,DlMC,
     &           truth,truthreg,truthsemi,truththree,truthomega,
     &           truthelse,
     &           truthdouble)
C-----------------------------------------------------------------------
C
C  Description:
C  ------------
C in: ntmc,motherOld(),vtxmcOld(),pidmcOld(),xvmc(),yvmc(),zvmc(),
C in: pxmc(),pymc(),pzmc(),nvtxmc
C out: truth,truthreg,truthsemi,truththree,truthomega,truthelse,ipmc(),KchMC(),K
C      neMC(),DtMC,DlMC                                                         
C
C-----------------------------------------------------------------------
      USE ANALYSISMODULE
      use iso_c_binding
C-----------------------------------------------------------------------
C
C External functions
C
C
C Local declarations
C
      INTEGER NpioKs,NpipKs,NpimKs,NepKs,NemKs,NmupKs,NmumKs,
     &        NothKs,NpioKl,NpipKl,NpimKl,NepKl,
     &        NemKl,NmupKl,NmumKl,NothKl,Nreg,Nother,
     &        Npiow,Npipw,Npimw,Nothw,NKS,NKL,Nisr,nvtxmc
      INTEGER ntmc,motherOld(MaxNvtxGen),vtxmcOld(MaxNtrkGen),
     &        pidmcOld(MaxNtrkGen)
      INTEGER KLpic(2),KLpio(2),KSpic(2),KSpio(2),i
      REAL    xvmc(MaxNvtxGen),yvmc(MaxNvtxGen),zvmc(MaxNvtxGen),Broots
      REAL    ipmc(3),KchMC(9),KneMC(9),Ks(9),Kl(9),DtMC,DlMC,PhiPMC(3)
      REAL    pxmc(MaxNtrkGen),pymc(MaxNtrkGen),pzmc(MaxNtrkGen)
      LOGICAL truth,truthreg,truthsemi,truththree,truthomega,truthelse,
     &        truthdouble
      LOGICAL condsig,condreg,condsemi,condthree,condomega,condelse,
     &        conddouble
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
      truthdouble = .FALSE.
C
      DO i = 1,ntmc
C Kl vertex
         IF( motherOld(vtxmcOld(i)).EQ.10 )THEN  
C pi0
            IF( pidmcOld(i).EQ.7 )THEN        
               NpioKl = NpioKl + 1
C pi+
            ELSEIF( pidmcOld(i).EQ.8 ) THEN   
               NpipKl = NpipKl + 1
C pi-
            ELSEIF( pidmcOld(i).EQ.9 ) THEN   
               NpimKl = NpimKl + 1
C e+
            ELSEIF( pidmcOld(i).EQ.2 ) THEN   
               NepKl = NepKl + 1
C e-
            ELSEIF( pidmcOld(i).EQ.3 ) THEN   
               NemKl = NemKl + 1
C mu+
            ELSEIF( pidmcOld(i).EQ.5 ) THEN   
               NmupKl = NmupKl + 1
C mu-
            ELSEIF( pidmcOld(i).EQ.6 ) THEN   
               NmumKl = NmumKl + 1
C regeneration
            ELSEIF( pidmcOld(i).EQ.16 ) THEN   
               Nreg = Nreg + 1
            ELSE
               NothKl = NothKl + 1
            ENDIF
C Ks vertex
         ELSEIF( motherOld(vtxmcOld(i)).EQ.16 )THEN   
C pi0
            IF( pidmcOld(i).EQ.7 )THEN             
               NpioKs = NpioKs + 1
C pi+
            ELSEIF( pidmcOld(i).EQ.8 ) THEN   
               NpipKs = NpipKs + 1
C pi-
            ELSEIF( pidmcOld(i).EQ.9 ) THEN   
               NpimKs = NpimKs + 1
C e+
            ELSEIF( pidmcOld(i).EQ.2 ) THEN   
               NepKs = NepKs + 1
C e-
            ELSEIF( pidmcOld(i).EQ.3 ) THEN   
               NemKs = NemKs + 1
C mu+
            ELSEIF( pidmcOld(i).EQ.5 ) THEN   
               NmupKs = NmupKs + 1
C mu-
            ELSEIF( pidmcOld(i).EQ.6 ) THEN   
               NmumKs = NmumKs + 1
            ELSE
               NothKs = NothKs + 1
            ENDIF
C Omega vertex
         ELSEIF( motherOld(vtxmcOld(i)).EQ.50 )THEN   
C pi0
            IF( pidmcOld(i).EQ.7 )THEN             
               Npiow = Npiow + 1
C pi+
            ELSEIF( pidmcOld(i).EQ.8 ) THEN   
               Npipw = Npipw + 1
C pi-
            ELSEIF( pidmcOld(i).EQ.9 ) THEN   
               Npimw = Npimw + 1
C Kl
            ELSEIF( pidmcOld(i).EQ.10 ) THEN   
               NKL = NKL + 1
C Ks
            ELSEIF( pidmcOld(i).EQ.16 ) THEN   
               NKS = NKS + 1
            ELSEIF( pidmcOld(i).EQ.1 ) THEN
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
      conddouble = ( NKL.EQ.1 .AND. NKS.EQ.1 .AND. NpipKl.EQ.1 .AND.
     &               NpimKl.EQ.1 .AND. NpipKs.EQ.1 .AND. NpimKs.EQ.1
     &               .AND.Nisr.GE.0 .AND. NpioKl.EQ.0 .AND. NpioKs.EQ.0 
     &)
      IF( NothKs.EQ.0 .AND. NothKl.EQ.0 .AND. NmumKs.EQ.0 .AND. NmupKs.E
     1  Q.0 .AND. NepKs.EQ.0 .AND. NemKs.EQ.0 .AND.                     
     &    NmumKl.EQ.0 .AND. NmupKl.EQ.0 .AND. NepKl.EQ.0 .AND. NemKl.EQ.
     1  0 .AND. Nreg.EQ.0 .AND.                                         
     &    Npiow.EQ.0 .AND. Npipw.EQ.0 .AND. Npimw.EQ.0 .AND. Nothw.EQ.0)
     & THEN
         IF ( conddouble ) THEN
            truthdouble = .TRUE.
         ELSE IF( condsig  ) THEN
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
           IF( motherOld(i).eq.50 ) THEN
              ipmc(1) = xvmc(i)
              ipmc(2) = yvmc(i)
              ipmc(3) = zvmc(i)
              EXIT
           ENDIF
        ENDDO
        DO i = 1,nvtxmc
C KL
           IF( motherOld(i).eq.10 ) THEN 
              Kl(7)=xvmc(i)
              Kl(8)=yvmc(i)
              Kl(9)=zvmc(i)
              EXIT
           ENDIF
        ENDDO
        DO i = 1,nvtxmc
C KS
           IF( motherOld(i).eq.16 ) THEN 
              Ks(7)=xvmc(i)
              Ks(8)=yvmc(i)
              Ks(9)=zvmc(i)
              EXIT
           ENDIF
        ENDDO
       DO i = 1,ntmc
CKL
           IF( pidmcOld(i).eq.10 ) THEN 
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
           IF( pidmcOld(i).eq.16 ) THEN 
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
           IF( motherOld(vtxmcOld(i)).EQ.10 ) THEN      
C pi0
              IF( pidmcOld(i).EQ.7 ) THEN             
                 KneMC=Kl
                 KchMC=Ks
                 EXIT
              ELSEIF( pidmcOld(i).EQ.8 .OR. pidmcOld(i).EQ.9 ) THEN
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
C 7 xvOld
C 8 yvOld
C 9 zvOld
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
      USE ANALYSISMODULE
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
      USE ANALYSISMODULE
C
      REAL      Psys(4), Pold(4), Pnew(4)
      REAL      Esys, Pmod, Msys, ScalProd
C
C-----------------------------------------------------------------------
C
      DO i = 1,4
        Pnew(i) = 0.
      ENDDO
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
      SUBROUTINE find_kchrec(findKS,findKL,last_vtx,findClose,
     &                       Bx,By,Bz,qualv,nv,ntv,
     &                       IV,CurV,PhiV,CotV,
     &                       xvOld,yvOld,zvOld,vtaken,KchRec,
     &                       trk1,trk2,cosTrk)
C-----------------------------------------------------------------------
      USE ANALYSISMODULE
C-----------------------------------------------------------------------
C in:qualv(),nv,ntv,IV(),xvOld(),yvOld(),zvOld()
C in:CurV(),PhiV(),CotV() or PxTV(),PyTV(),PzTV()
C out:vtaken(),Kchrec(),trk1(),trk2(),cosTrk
C 1 px
C 2 py
C 3 pz
C 4 E
C 5 |p|
C 6 minv
C 7 xvOld
C 8 yvOld
C 9 zvOld
C Local variables
      INTEGER findKS,findKL,findClose
      INTEGER vtaken(3),i,j,k,l,nv,ntv,ivOld(MaxNumtrkv),ii,last_vtx
      REAL    KchRec(9),KchTmp(9),P4trkRec1(4),P4trkRec2(4)
      REAL    CurV(MaxNumtrkv),PhiV(MaxNumtrkv),CotV(MaxNumtrkv)
      REAL    PxTV(MaxNumtrkv),PyTV(MaxNumtrkv),PzTV(MaxNumtrkv)
      REAL    xvOld(MaxNumVtx),yvOld(MaxNumVtx),zvOld(MaxNumVtx)
      INTEGER qualv(MaxNumVtx)
      REAL    P4trkRec1mod,P4trkRec2mod
      REAL    trk1(4),trk2(4),cosTrk
      REAL    cyl_vol(2),Bx,By,Bz,dist,dist_tmp
      LOGICAL pipi_bool
C --- standard values for klspm00 selection
C -----------------------------------------
C --- rec

      DO i = 1,3
         vtaken(i)=0
      ENDDO
      dist = 1000000.
      KchRec(6) = 1000000000.
      DO i = 1,nv
         cyl_vol(1) = SQRT( (xvOld(i) - Bx)**2 + (yvOld(i) - By)**2 )
         cyl_vol(2) = ABS( zvOld(i) - Bz )
C --- Find the conditions for KS or KL candidate
         IF(findKS .EQ. 1) THEN
            pipi_bool = (cyl_vol(1).LE.10) .AND. (cyl_vol(2).LE.20)
         ELSE IF(findKL .EQ. 1) THEN
            pipi_bool = i.NE.(last_vtx)
         ELSE
            pipi_bool = .TRUE.
         ENDIF

         IF(pipi_bool) THEN
            DO j = 1,ntv-1
               IF( IV(j).eq.i ) THEN
                  P4trkRec1(1) = COS( PhiV(j) ) / ABS( CurV(j) )*1000.
                  P4trkRec1(2) = SIN( PhiV(j) ) / ABS( CurV(j) )*1000.
                  P4trkRec1(3) =      CotV(j)   / ABS( CurV(j) )*1000.
C
                  P4trkRec1mod = P4trkRec1(1)**2 + P4trkRec1(2)**2 +
     &                           P4trkRec1(3)**2
                  P4trkRec1(4) = SQRT( P4trkRec1mod + Mpip**2 )
                  P4trkRec1mod = SQRT( P4trkRec1mod )
                  DO k = j+1,ntv
                     IF( IV(k).eq.i ) THEN
                        P4trkRec2(1) = COS( PhiV(k) ) / ABS( CurV(k) )*
     &                                                   1000.
                        P4trkRec2(2) = SIN( PhiV(k) ) / ABS( CurV(k) )*
     &                                                   1000.
                        P4trkRec2(3) =      CotV(k)   / ABS( CurV(k) )*
     &                                                   1000.
C
                        P4trkRec2mod = P4trkRec2(1)**2 + P4trkRec2(2)**2
     & +
     &                              P4trkRec2(3)**2
                        P4trkRec2(4) = SQRT( P4trkRec2mod + Mpip**2 )
                        P4trkRec2mod = SQRT( P4trkRec2mod )
                        DO l = 1,4
                           KchTmp(l) = P4trkRec1(l) + P4trkRec2(l)
                        ENDDO
                        KchTmp(5) = KchTmp(1)**2 + KchTmp(2)**2 +
     &                              KchTmp(3)**2
                        KchTmp(6) = SQRT( KchTmp(4)**2 - KchTmp(5) )
                        KchTmp(5) = sqrt( KchTmp(5) )

                        IF( findClose .EQ. 0) THEN
                           IF( (vtaken(1).eq.0).or.
     &                        ( ABS( KchTmp(6) - Mko ).lt.
     &                           ABS( KchRec(6) - Mko ) ) ) THEN
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

                              
                              cosTrk=( trk1(1)*trk2(1)+
     &                                 trk1(2)*trk2(2)+
     &                                 trk1(3)*trk2(3) ) /
     &                                 (P4trkRec1mod*P4trkRec2mod)
                           ENDIF
                        ELSE
                           dist_tmp = SQRT((xvOld(i) - Bx)**2 
     &                                      + (yvOld(i) - By)**2) 
                           IF( (vtaken(1).eq.0).or.
     &                        ( dist_tmp.lt.
     &                           dist ) ) THEN
C vtx id
                              dist = dist_tmp
                        
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
                              cosTrk=( trk1(1)*trk2(1)+
     &                                 trk1(2)*trk2(2)+
     &                                 trk1(3)*trk2(3) ) /
     &                                 (P4trkRec1mod*P4trkRec2mod)
                           ENDIF
                        ENDIF
                     ENDIF
                  ENDDO
               ENDIF
            ENDDO
         ENDIF
      ENDDO
      IF( vtaken(1).ne.0 ) THEN
         KchRec(7) = xvOld(vtaken(1))
         KchRec(8) = yvOld(vtaken(1))
         KchRec(9) = zvOld(vtaken(1))
      ELSE
         ErrFlag = .TRUE.
      ENDIF
      END

C=======================================================================
C=======================================================================
C=======================================================================
      SUBROUTINE cor_ip_boost(KchRec,Bpx,Bpy,Bpz,Bx,By,Bz,Broots,
     &                        trk1,trk2,KchBoost,ip_closest,
     &                        ip_plane,chdist,Qmiss)
C-----------------------------------------------------------------------
      USE ANALYSISMODULE
C-----------------------------------------------------------------------
C External function
      REAL pk_from_boost
      INTEGER dist_lines
      INTEGER plane_intersection
C 1 px
C 2 py
C 3 pz
C 4 E
C 5 |p|
C 6 minv
C 7 xvOld
C 8 yvOld
C 9 zvOld
C in: KchRec(),Bpx,Bpy,Bpz,Bx,By,Bz,Broots,trk1(),trk2()
C out: KchBoost(),ip,chdist,Qmiss
C Local variables
      INTEGER i,status_cl,status_pl
      REAL KchRec(9),KchBoost(9),ip_closest(3),chdist,P4KchBoost(4)
      REAL P4PhiBha(4),pboost,xB(3),x(3),beamLine(3),p(3),pktP(3)
      REAL plane_perp(3),ip_plane(3)
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
C --- IP - both plane and closest approach SG
      xB(1) = Bx
      xB(2) = By
      xB(3) = Bz
      x(1) = KchBoost(7)
      x(2) = KchBoost(8)
      x(3) = KchBoost(9)
      BeamLine(1) = 0.
      BeamLine(2) = 0.
      BeamLine(3) = 1.
      plane_perp(1) = 0.
      plane_perp(2) = P4PhiBha(2)
      plane_perp(3) = 0.
      DO i = 1,3
         p(i) = KchBoost(i)
      ENDDO
      status_cl = dist_lines(x,p,xB,BeamLine,pktP,pktB,chdist)
      status_pl = plane_intersection(x,p,xB,plane_perp,ip_plane)
C     integer function dist_lines (x1, p1, x2, p2, x1PCA, x2PCA, dist)
C     x1PCA     point on line 1 of closest approach to line 2
C     x2PCA     point on line 2 of closest approach to line 1
C     dist      distance of closest approach
      IF( status_cl.eq.0 .and. status_pl.eq.0 ) THEN
         ErrFlag = .TRUE.
         RETURN
      ENDIF
      ip_closest(1) = xB(1)
      ip_closest(2) = xB(2)
      ip_closest(3) = pktB(3)
      Emiss = Kchboost(4) - trk1(4) - trk2(4)
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
      USE ANALYSISMODULE         
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
      INTEGER FUNCTION plane_intersection (X_line, vec_line, X_plane,
     &                                     vec_plane, int_point)
      IMPLICIT NONE

      REAL X_line(3), vec_line(3), X_plane(3), vec_plane(3),
     &     int_point(3), dot_prod_up, dot_prod_down
      plane_intersection = 0
      dot_prod_up = (X_line(1) - X_plane(1)) * vec_plane(1) +
     &              (X_line(2) - X_plane(2)) * vec_plane(2) +
     &              (X_line(3) - X_plane(3)) * vec_plane(3)
      dot_prod_down = vec_line(1) * vec_plane(1) +
     &                vec_line(2) * vec_plane(2) +
     &                vec_line(3) * vec_plane(3)
      if(dot_prod_down.eq.0) goto 13
      int_point(1) = X_line(1) + (dot_prod_up/dot_prod_down) *
     &                                                    vec_line(1)
      int_point(2) = X_line(2) + (dot_prod_up/dot_prod_down) *
     &                                                    vec_line(2)
      int_point(3) = X_line(3) + (dot_prod_up/dot_prod_down) *
     &                                                    vec_line(3)
      plane_intersection = 1
 13   RETURN
      END
C=======================================================================
C=======================================================================
C=======================================================================
      SUBROUTINE cos_pipi_kaonCM (PxPhi,PyPhi,PzPhi,SqrtS,
     &                                  ip_vtx,KchLAB,Ppich1,
     &                                  Ppich2,costrk)

      USE ANALYSISMODULE

      IMPLICIT NONE

C For me event of phi meson decay was on t = 0
C We need the interval, any offset will delete
C Lorentz transf is linear
      REAL PxPhi, PyPhi, PzPhi, SqrtS, ip_vtx(3), KchLAB(9)
      REAL Ppich1(4), Ppich2(4), costrk
      REAL betaKch, dist_Kch, mom_Kch, numerator, denominator
      REAL mom1_len, mom2_len
      REAL PBhabha(4), ChVtxLAB(4), ChVtxCM(4)
      REAL PKchCM(4), Ppich1LAB(4), Ppich2LAB(4)
      REAL Ppich1CM(4), Ppich2CM(4)
      REAL Ppich1CMK(4), Ppich2CMK(4)
      PBhabha(1) = PxPhi
      PBhabha(2) = PyPhi
      PBhabha(3) = PzPhi
      PBhabha(4) = SqrtS
      betaKch = Cvel * KchLAB(5) / KchLAB(4)
      ChVtxLAB(1) = KchLAB(7) - ip_vtx(1)
      ChVtxLAB(2) = KchLAB(8) - ip_vtx(2)
      ChVtxLAB(3) = KchLAB(9) - ip_vtx(3)
      ChVtxLAB(4) = SQRT(ChVtxLAB(1)**2 +
     &                   ChVtxLAB(2)**2 +
     &                   ChVtxLAB(3)**2 ) / betaKch
      CALL Loren4(PBhabha,ChVtxLAB,ChVtxCM)
C Going with pich to CM
      Ppich1LAB(1) = Ppich1(1)
      Ppich1LAB(2) = Ppich1(2)
      Ppich1LAB(3) = Ppich1(3)
      Ppich1LAB(4) = Ppich1(4)
      Ppich2LAB(1) = Ppich2(1)
      Ppich2LAB(2) = Ppich2(2)
      Ppich2LAB(3) = Ppich2(3)
      Ppich2LAB(4) = Ppich2(4)
      CALL Loren4(PBhabha,Ppich1LAB,Ppich1CM)
      CALL Loren4(PBhabha,Ppich2LAB,Ppich2CM)
C Determination of Kaon Momentum using vtx
      dist_Kch = SQRT(ChVtxCM(1)**2 +
     &                ChVtxCM(2)**2 +
     &                ChVtxCM(3)**2 )
      mom_Kch = SQRT( (SqrtS**2 / 4.) - Mko**2 )
      PKchCM(1) = mom_Kch * ChVtxCM(1) / dist_Kch
      PKchCM(2) = mom_Kch * ChVtxCM(2) / dist_Kch
      PKchCM(3) = mom_Kch * ChVtxCM(3) / dist_Kch
      PKchCM(4) = SqrtS / 2.
C Having Kaon Mom in Phi CM we can go to Kaon CM
      CALL Loren4(PKchCM,Ppich1CM,Ppich1CMK)
      CALL Loren4(PKchCM,Ppich2CM,Ppich2CMK)
      mom1_len = SQRT(Ppich1CMK(1)**2 +
     &                Ppich1CMK(2)**2 +
     &                Ppich1CMK(3)**2 )
      mom2_len = SQRT(Ppich2CMK(1)**2 +
     &                Ppich2CMK(2)**2 +
     &                Ppich2CMK(3)**2 )
      numerator = Ppich1CMK(1)*Ppich2CMK(1) +
     &            Ppich1CMK(2)*Ppich2CMK(2) +
     &            Ppich1CMK(3)*Ppich2CMK(3)
      denominator = mom1_len * mom2_len
      IF(denominator.EQ.0) THEN
         ErrFlag = .TRUE.
         costrk = -2.
         GOTO 13
      ENDIF
      costrk = numerator / denominator
      
 13   RETURN
      END
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

      USE ANALYSISMODULE
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
           WRITE(*,999) nmctruth(3),'regeneration'
           WRITE(*,999) nmctruth(4),'omega-pi0'
           WRITE(*,999) nmctruth(5),'3pi0-pi+pi-'
           WRITE(*,999) nmctruth(6),'semileptonic'
           WRITE(*,999) nmctruth(8),'pi+pi-pi+pi-'
           WRITE(*,999) nmctruth(7),'other background'
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

      SUBROUTINE check_isr(ntmc,pidmcOld,virmom,vtxmcOld,motherOld,isr)
         USE ANALYSISMODULE
C-----------------------------------------------------------------------
C in: ntmc,pidmcOld(),virmom(),vtxmcOld(),motherOld()
C out: isr
      INTEGER isr,i,ntmc,pidmcOld(MaxNtrkGen),virmom(MaxNtrkGen)
      INTEGER vtxmcOld(MaxNtrkGen),motherOld(MaxNvtxGen)
      ISR = 0
      DO i = 1,ntmc
          IF( pidmcOld(i).eq.1.and.virmom(i).eq.1.and.
     &        vtxmcOld(i).eq.1.and.
     &        motherOld(vtxmcOld(i)).eq.50 ) THEN
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
      USE ANALYSISMODULE
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
      SUBROUTINE find_neuvtx(Bpx,Bpy,Bpz,Broots,KchBoost,ncl,enecl,
     &           ncll,xcl,ycl,zcl,ip,tcl,cldist,KneRecLor,trc,
     &           nclwrong,ncllwrong,KneRec,minv4gam,pi0,
     &           g4taken,trcv,PgamRecTaken,gpairtaken,test,
     &           g4vtxerr)
   
C-----------------------------------------------------------------------
      USE ANALYSISMODULE
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
      REAL    TclOld(MaxNumClu),cldist,trc_tmp
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
C 7 xvOld
C 8 yvOld
C 9 zvOld
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

      USE ANALYSISMODULE
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

      SUBROUTINE do_cuts

      USE ANALYSISMODULE
C-----------------------------------------------------------------------
C Local variables
      INTEGER i
      REAL    tmp
C-----------------------------------------------------------------------
C
      EventOK = .FALSE.
      IF(kinfitsig) THEN
         IF( Interf%Chi2.gt.40 ) THEN
         CALL statisticss(1)
         RETURN
         ENDIF
         IF( sqrt((Interf%pi0Fit(1)-Mpio)**2+
     &         (Interf%pi0Fit(2)-Mpio)**2).gt.35 ) THEN
         CALL statisticss(2)
         RETURN
         ENDIF
      ENDIF
C      tmp=0.
C      DO i=1,4
C         tmp=tmp+Interf%trcv(Interf%g4taken(i))
C      ENDDO
C      IF( tmp.lt.-1 ) THEN
C        CALL statisticss(3)
C        RETURN
C      ENDIF
       IF( Interf%minv4gam.lt.300 .OR. Interf%minv4gam.gt.1000) THEN
         CALL statisticss(4)
         RETURN
       ENDIF
C       IF( abs(Interf%kchrec(6)-Mko).gt.1.2 ) THEN
C         CALL statisticss(5)
C         RETURN
C       ENDIF
C      IF( Interf%Qmiss.gt.3.75 ) THEN
C        CALL statisticss(6)
C        RETURN
C      ENDIF
C      IF( Interf%cosTrkCM.gt.-0.8 ) THEN
C        CALL statisticss(7)
C        RETURN
C      ENDIF
      EventOK = .TRUE.
CC probably do_distance function is wrong, (Eryk, 25.06.2013)
      END
C================================================================================
C================================================================================
C================================================================================

      SUBROUTINE statisticss(id)

      USE ANALYSISMODULE
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
      
      USE ANALYSISMODULE
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
      USE ANALYSISMODULE
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
        PARAMETER (N=38)                  
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
         USE ANALYSISMODULE
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
           Ptrk(ILoop,1) = P(3*(ILoop-1)+1)
           Ptrk(ILoop,2) = P(3*(ILoop-1)+2)
           Ptrk(ILoop,3) = P(3*(ILoop-1)+3)
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
           Ltot = SQRT( (P(2+OffSet)-P(30))**2 +
     &          (P(3+OffSet)-P(31))**2 + (P(4+OffSet)-P(32))**2 )
           Pgam(ILoop,1) = P(1+OffSet) * (P(2+OffSet)-P(30))/Ltot
           Pgam(ILoop,2) = P(1+OffSet) * (P(3+OffSet)-P(31))/Ltot
           Pgam(ILoop,3) = P(1+OffSet) * (P(4+OffSet)-P(32))/Ltot
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
        Ltot = SQRT( (P(30)-P(33))**2 +
     &       (P(31)-P(34))**2 + (P(32)-P(35))**2 )
        Vkne = Vgam * SQRT(PxTot**2+PyTot**2+PzTot**2) / EnTot
        Tkne = Ltot / Vkne
C
C------------------------------------------------------------------------------
C Constraints
C------------------------------------------------------------------------------
C
C Total ToF for photons
C
        DO ILoop = 1,4
           OffSet = 5*(ILoop-1) + 10
           C(ILoop) = Tkne + Tgam(ILoop) - P(OffSet+1)
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
        SinThetaBoost = P(38) / (P(36)+P(37))
        CosThetaBoost = SQRT(1-SinThetaBoost**2)
C
C == SinThetaBoost * (P(36)+P(37))
        PxPhi = P(38)      
        PyPhi = 0.
C eryk ERYK: in MC  PzPhi = 0
        PzPhi = CosThetaBoost * (P(36)-P(37))
        EnPhi = P(36) + P(37)
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
        USE ANALYSISMODULE
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
        PARAMETER (N=38)                  
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
         USE ANALYSISMODULE
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
           Ptrk(ILoop,1) = P(3*(ILoop-1)+1)
           Ptrk(ILoop,2) = P(3*(ILoop-1)+2)
           Ptrk(ILoop,3) = P(3*(ILoop-1)+3)
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
           Ltot = SQRT( (P(2+OffSet)-P(33))**2 +
     &          (P(3+OffSet)-P(34))**2 + (P(4+OffSet)-P(35))**2 )
           Pgam(ILoop,1) = P(1+OffSet) * (P(2+OffSet)-P(33))/Ltot
           Pgam(ILoop,2) = P(1+OffSet) * (P(3+OffSet)-P(34))/Ltot
           Pgam(ILoop,3) = P(1+OffSet) * (P(4+OffSet)-P(35))/Ltot
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
        SinThetaBoost = P(38) / (P(36)+P(37))
        CosThetaBoost = SQRT(1-SinThetaBoost**2)
C
C == SinThetaBoost * (P(36)+P(37))
        PxPhi = P(38)      
        PyPhi = 0.
        PzPhi = CosThetaBoost * (P(36)-P(37))
        EnPhi = P(36) + P(37)
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

      SUBROUTINE test_array(array, numbers, array1)

      USE ANALYSISMODULE

      REAL*8 array(4,2)
      REAL*8 array1(4)
      REAL numbers

      DO i = 1,4
         DO j = 1,2
            array(i,j) = 111.
         ENDDO

         array1(i) = 999.
      ENDDO

      numbers = 49.

      END




