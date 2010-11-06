C                       *****************
                        SUBROUTINE QTRIA1
C                       *****************
C
     *( F     , XK    , FREQ  , DEPTH , RAISF , GRAVIT, ALFLTA, RFMLTA,
     *  NF    , NPLAN , NPOIN2, TSTOT , TSDER , FTOT  , FMOY  )
C
C**********************************************************************
C  TOMAWAC - V1.1                           (EDF/DER/LNH)  -  26/12/96
C**********************************************************************
C
C  FONCTION : TERME SOURCE LIE AUX INTERACTIONS NON-LINEAIRES ENTRE
C  ********** TRIPLETS DE FREQUENCES.
C
C  ARGUMENTS :
C  ***********
C  +-------------+----+--------------------------------------------+
C  ! NOM         !MODE! SIGNIFICATION - OBSERVATIONS               !
C  +-------------+----+--------------------------------------------+
C  ! F(-,-,-)    !<-->! SPECTRE DIRECTIONNEL DE VARIANCE           !
C  ! XK(-,-)     ! -->! TABLEAU DES NOMBRES D'ONDE                 !
C  ! FREQ(-)     ! -->! TABLEAU DES FREQUENCES DE DISCRETISATION   !
C  ! DEPTH(-)    ! -->! TABLEAU DES PROFONDEURS (METRES)           !
C  ! RAISF       ! -->! RAISON FREQUENTIELLE POUR DISCRETISATION   !
C  ! GRAVIT      ! -->! ACCELERATION DE LA PESANTEUR               !
C  ! NF          ! -->! NOMBRE DE FREQUENCES DE DISCRETISATION     !
C  ! NPLAN       ! -->! NOMBRE DE DIRECTIONS DE DISCRETISATION     !
C  ! NPOIN2      ! -->! NOMBRE DE POINTS DU MAILLAGE SPATIAL       !
C  ! TSTOT(-,-,-)!<-- ! CONTRIBUTION TERME SOURCE - PARTIE TOTALE  !
C  ! TSDER(-,-,-)!<-- ! CONTRIBUTION TERME SOURCE - PARTIE DERIVEE !
C  ! TAUX1(-)    !<-->! TABLEAU DE TRAVAIL (DIMENSION NPOIN2)      !
C  ! TAUX2(-)    !<-->! TABLEAU DE TRAVAIL (DIMENSION NPOIN2)      !
C  +-------------+----+--------------------------------------------+
C  ! MODE   (-> : NON-MODIFIE)  (<-> : MODIFIE)  (<- : INITIALISE) !
C  +---------------------------------------------------------------+
C
C  APPELS :    - PROGRAMME(S) APPELANT  :  SEMIMP
C  ********    - PROGRAMME(S) APPELE(S) :  ------
C
C  REMARQUES :
C  ***********
C
C**********************************************************************
C
      IMPLICIT NONE
C
C.....VARIABLES TRANSMISES
C     """"""""""""""""""""
      INTEGER  NF, NPOIN2, NPLAN
      DOUBLE PRECISION  F(NPOIN2,NPLAN,NF), XK(NPOIN2,NF)
      DOUBLE PRECISION  RAISF , FREQ(NF), DEPTH(NPOIN2)
      DOUBLE PRECISION  TSTOT(NPOIN2,NPLAN,NF), TSDER(NPOIN2,NPLAN,NF)
      DOUBLE PRECISION  GRAVIT, ALFLTA, RFMLTA
      DOUBLE PRECISION  FTOT(NPOIN2) , FMOY(NPOIN2)
C
C.....VARIABLES LOCALES
C     """""""""""""""""
      INTEGER           IIND  , IFMA  , IFF   , IPL   , IPO
      DOUBLE PRECISION  CPH   , CGR   , BIF   , RPS2  , RP    , DEUPI ,
     *                  FPS2  , CPHS2 , XKPS2 , URS   , RIND  , FMAX  ,
     *                  F2P   , CPH2P , XK2P  , CGR2P , E2P   , EPS2  ,
     *                  OMP   , OM2P  , COEF  , D     , DEUKD ,
     *                  SPLUS , SMOIN , DMOIN , DPLUS , FP    , XKP
      PARAMETER(DEUPI=2.D0*3.141592654D0)
C
C.....FONCTION FORMULE
C     """"""""""""""""
      DOUBLE PRECISION  RPP   , DD    , O2P   , CP    , KP    , K2P
      RPP(KP,CP,K2P,O2P,DD) = KP**2*(GRAVIT*DD+2.D0*CP**2)/(K2P*DD)/
     *   (GRAVIT*DD+(2.D0/15.D0)*GRAVIT*DD**3*K2P**2-0.4D0*(O2P*DD)**2)
C
C
      COEF  = GRAVIT*DSQRT(2.D0)
C
      DO IPO=1,NPOIN2
C
        D=DEPTH(IPO)
C
C.......CALCUL DU NOMBRE D'URSELL AU POINT D'ESPACE CONSIDERE
C       """""""""""""""""""""""""""""""""""""""""""""""""""""
        URS = COEF*DSQRT(FTOT(IPO))/(DEUPI*D*FMOY(IPO))**2
C
C.......ON N'ACTIVE LE CALCUL QUE SI URSELL > 0.1
C       """""""""""""""""""""""""""""""""""""""""
        IF (URS.GT.0.1D0) THEN
C
C.........CALCUL DU SINUS DE LA BIPHASE
C         """""""""""""""""""""""""""""
          IF (URS.GT.10.D0) THEN
            BIF=1.D0
          ELSE
            BIF=DSIN(DABS(DEUPI/4.D0*(-1.D0+DTANH(0.2D0/URS))))
          ENDIF
C
C.........CALCUL DE L'INDICE FREQUENTIEL MAXIMUM
C         """"""""""""""""""""""""""""""""""""""
          FMAX = RFMLTA*MAX(FMOY(IPO),FREQ(1))
          RIND = 1.D0 + DLOG(FMAX/FREQ(1))/DLOG(RAISF)
          IFMA = MIN(IDINT(RIND),NF)
C
          DO IFF=1,IFMA
            FP  = FREQ(IFF)
            OMP = DEUPI*FP
            XKP = XK(IPO,IFF)
            CPH = OMP/XKP
C
C
C...........CALCUL DE LA CONTIBUTION S+
C           """""""""""""""""""""""""""
            FPS2 = FP/2.D0
            RIND = 1.D0 + DLOG(FPS2/FREQ(1))/DLOG(RAISF)
            IIND = IDINT(RIND)
            RIND = RIND-DBLE(IIND)
C
            IF (IIND.GT.0) THEN
              DEUKD=2.D0*XKP*D
              IF(DEUKD.LE.7.D2) THEN
                CGR = CPH*(0.5D0+XKP*D/DSINH(2.D0*XKP*D))
              ELSE
                CGR = 0.5D0*CPH
              ENDIF
              CALL WNSCOU(XKPS2,FPS2,D)
              CPHS2 = DEUPI*FPS2/XKPS2
              RPS2 = CPH*CGR*RPP(XKPS2,CPHS2,XKP,OMP,D)**2
C
              DO IPL=1,NPLAN
                EPS2=(1.D0-RIND)*F(IPO,IPL,IIND)+RIND*F(IPO,IPL,IIND+1)
                SPLUS = ALFLTA*RPS2*BIF*(EPS2-2.D0*F(IPO,IPL,IFF))*EPS2
                DPLUS = ALFLTA*RPS2*BIF*(-2.D0*EPS2)
                IF (SPLUS.LT.0.D0) THEN
                  SPLUS = 0.D0
                  DPLUS = 0.D0
                ENDIF
                TSTOT(IPO,IPL,IFF) = TSTOT(IPO,IPL,IFF) + SPLUS
C                TSDER(IPO,IPL,IFF) = TSDER(IPO,IPL,IFF) + DPLUS
              ENDDO
            ENDIF
C
C
C...........CALCUL DE LA CONTIBUTION S-
C           """""""""""""""""""""""""""
            F2P = 2.D0*FP
            RIND = 1.D0 + DLOG(F2P/FREQ(1))/DLOG(RAISF)
            IIND = IDINT(RIND)
            RIND = RIND-DBLE(IIND)
            IF (IIND.LT.IFMA) THEN
              OM2P  = DEUPI*F2P
              CALL WNSCOU(XK2P,F2P,D)
              CPH2P = OM2P/XK2P
              DEUKD=2.D0*XK2P*D
              IF(DEUKD.LE.700D2) THEN
                CGR2P = CPH2P*(0.5D0+XK2P*D/DSINH(2.D0*XK2P*D))
              ELSE
                CGR2P = CPH2P*0.5D0
              ENDIF
              RP    = CPH2P*CGR2P*RPP(XKP,CPH,XK2P,OM2P,D)**2
C
              DO IPL=1,NPLAN
                E2P = (1.D0-RIND)*F(IPO,IPL,IIND)+RIND*F(IPO,IPL,IIND+1)
                SMOIN = 2.D0*ALFLTA*RP*BIF*F(IPO,IPL,IFF)
     *                *(F(IPO,IPL,IFF)-2.D0*E2P)
                DMOIN = 4.D0*ALFLTA*RP*BIF*(F(IPO,IPL,IFF)-E2P)
                IF (SMOIN.LT.0.D0) THEN
                  SMOIN = 0.D0
                  DMOIN = 0.D0
                ENDIF
                TSTOT(IPO,IPL,IFF) = TSTOT(IPO,IPL,IFF) - SMOIN
C                TSDER(IPO,IPL,IFF) = TSDER(IPO,IPL,IFF) - DMOIN
              ENDDO
            ENDIF
C
          ENDDO
        ENDIF
      ENDDO
C
      RETURN
      END
