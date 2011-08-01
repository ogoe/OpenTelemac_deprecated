C                 
C                       *****************
                        SUBROUTINE QWIND1
C                       *****************
C
     *( TSTOT , TSDER , F     , XK    , FREQ  , USOLD , USNEW , TWOLD ,
     *  TWNEW , Z0OLD , Z0NEW , TETA  , ROAIR , ROEAU , BETAM , XKAPPA,
     *  DECAL , GRAVIT, NF    , NPLAN , NPOIN2, CIMPLI, TOLD  , TNEW  ,
     *  CPHAS , USN   , USO   , OMNEW , OMOLD , BETAN , BETAO )
C
C**********************************************************************
C  TOMAWAC - V1P0    P. THELLIER & M. BENOIT (EDF/DER/LNH)  -  11/04/95
C**********************************************************************
C
C  FONCTION : CALCUL DE LA CONTRIBUTION DU TERME SOURCE DE GENERATION
C  ********** PAR LE VENT. UTILISATION DE LA THEORIE DE JANSSEN (1989,
C             1991).
C
C  ARGUMENTS :
C  ***********
C +-------------+----+--------------------------------------------+
C ! NOM         !MODE! SIGNIFICATION - OBSERVATIONS               !
C +-------------+----+--------------------------------------------+
C ! TSTOT(-,-,-)!<-->! CONTRIBUTION TERME SOURCE - PARTIE TOTALE  !
C ! TSDER(-,-,-)!<-->! CONTRIBUTION TERME SOURCE - PARTIE DERIVEE !
C ! F(-,-,-)    ! -->! SPECTRE DIRECTIONNEL                       !
C ! XK(-,-)     ! -->! NOMBRE D'ONDE                              !
C ! FREQ(-)     ! -->! TABLEAU DES FREQUENCES DE DISCRETISATION   !
C ! USOLD(-)    ! -->! VITESSE DE FROTTEMENT A L'INSTANT N        !
C ! USNEW(-)    ! -->! VITESSE DE FROTTEMENT A L'INSTANT N+1      !
C ! TWOLD(-)    ! -->! DIRECTION DU VENT A L'INSTANT N            !
C ! TWNEW(-)    ! -->! DIRECTION DU VENT A L'INSTANT N+1          !
C ! Z0OLD(-)    ! -->! LARGEUR DE RUGOSITE L'INSTANT N            !
C ! Z0NEW(-)    ! -->! LARGEUR DE RUGOSITE L'INSTANT N+1          !
C ! TETA(-)     ! -->! TABLEAU DES DIRECTIONS DE DISCRETISATION   !
C ! ROAIR       ! -->! DENSITE DE L'AIR                           !
C ! ROEAU       ! -->! DENSITE DE L'EAU                           !
C ! BETAM       ! -->! PARAMETRE D'APPORT DU AU VENT              !
C ! XKAPPA      ! -->! CONSTANTE DE VON KARMAN                    !
C ! DECAL       ! -->! CONSTANTE DE DECALAGE DE COURBE CROISSANCE !
C ! GRAVIT      ! -->! ACCELERATION DE LA PESANTEUR               !
C ! NF          ! -->! NOMBRE DE FREQUENCES DE DISCRETISATION     !
C ! NPLAN       ! -->! NOMBRE DE DIRECTIONS DE DISCRETISATION     !
C ! NPOIN2      ! -->! NOMBRE DE POINTS DU MAILLAGE SPATIAL       !
C ! TOLD(-,-)   !<-->! TABLEAU DE TRAVAIL (DIMENSION NPOIN2*NPLAN)!
C ! TNEW(-,-)   !<-->! TABLEAU DE TRAVAIL (DIMENSION NPOIN2*NPLAN)!
C ! CPHAS(-)    !<-->! TABLEAU DE TRAVAIL (DIMENSION NPOIN2)      !
C ! USN(-)      !<-->! TABLEAU DE TRAVAIL (DIMENSION NPOIN2)      !
C ! USO(-)      !<-->! TABLEAU DE TRAVAIL (DIMENSION NPOIN2)      !
C ! OMNEW(-)    !<-->! TABLEAU DE TRAVAIL (DIMENSION NPOIN2)      !
C ! OMOLD(-)    !<-->! TABLEAU DE TRAVAIL (DIMENSION NPOIN2)      !
C ! BETAN(-)    !<-->! TABLEAU DE TRAVAIL (DIMENSION NPOIN2)      !
C ! BETAO(-)    !<-->! TABLEAU DE TRAVAIL (DIMENSION NPOIN2)      !
C +-------------+----+--------------------------------------------+
C ! MODE   (-> : NON-MODIFIE)  (<-> : MODIFIE)  (<- : INITIALISE) !
C +---------------------------------------------------------------+
C
C  APPELS :    - PROGRAMME(S) APPELANT  : SEMIMP
C  ********    - PROGRAMME(S) APPELE(S) :    -
C
C  REMARQUES :
C  ***********
C
C  REFRENCES : - JANSSEN P.A.E.M (1989) : WIND-INDUCED STRESS AND THE
C              DRAG OF AIR FLOW OVER SEA WAVES. JPO, VOl 19, PP 745-754
C              - JANSSEN P.A.E.M (1991) : QUASI-LINEAR THEORY OF WIND-
C              WAVE GENERATION APPLIED TO WAVE FORECASTING. JPO, VOL 21
C              PP 1631-1642.
C
C**********************************************************************
C
      IMPLICIT NONE
C
C.....VARIABLES TRANSMISES
C     """"""""""""""""""""
      INTEGER  NF    , NPLAN , NPOIN2
      DOUBLE PRECISION GRAVIT, ROAIR , ROEAU , BETAM , XKAPPA, DECAL ,
     *                 CIMPLI
      DOUBLE PRECISION  FREQ(NF)    , TETA(NPLAN)  , USOLD(NPOIN2)
      DOUBLE PRECISION TWOLD(NPOIN2), TWNEW(NPOIN2), USNEW(NPOIN2)
      DOUBLE PRECISION Z0OLD(NPOIN2), Z0NEW(NPOIN2), BETAN(NPOIN2)
      DOUBLE PRECISION CPHAS(NPOIN2),   USN(NPOIN2), OMOLD(NPOIN2)
      DOUBLE PRECISION OMNEW(NPOIN2),   USO(NPOIN2), BETAO(NPOIN2)
      DOUBLE PRECISION TSTOT(NPOIN2,NPLAN,NF), TSDER(NPOIN2,NPLAN,NF)
      DOUBLE PRECISION     F(NPOIN2,NPLAN,NF),    XK(NPOIN2,NF)
      DOUBLE PRECISION TOLD(NPOIN2,NPLAN)    ,  TNEW(NPOIN2,NPLAN)
C
C.....VARIABLES LOCALES
C     """""""""""""""""
      INTEGER  JP    , JF    , IP
      DOUBLE PRECISION DEUPI , C1 , DIREC , CONST , DIMPLI
      DOUBLE PRECISION XX    , ZLOGMU
C
C
      DEUPI = 2.D0* 3.14159265358978D0
      C1 = DEUPI * (ROAIR/ROEAU) * (BETAM/XKAPPA**2.D0)
      DIMPLI=1.0D0-CIMPLI
C
C.....PRECALCULS SUR LES DEPENDANCES DIRECTIONNELLES
C     """"""""""""""""""""""""""""""""""""""""""""""
      DO JP=1,NPLAN
        DIREC=TETA(JP)
        DO IP=1,NPOIN2
          TOLD(IP,JP)=COS(DIREC-TWOLD(IP))
          TNEW(IP,JP)=COS(DIREC-TWNEW(IP))
        ENDDO
      ENDDO
C
C.....BOUCLE SUR LES FREQUENCES DE DISCRETISATION.
C     """"""""""""""""""""""""""""""""""""""""""""
      DO JF=1,NF
        CONST=C1*FREQ(JF)
C
C.......CALCUL DE LA VITESSE DE PHASE
C       """""""""""""""""""""""""""""
        DO IP=1,NPOIN2
          CPHAS(IP) = DEUPI * FREQ(JF) / XK(IP,JF)
        ENDDO
C
C.......PRE-CALCUL FREQUENTIELS (OMEGA ET UETOILE/CPHASE)
C       """""""""""""""""""""""""""""""""""""""""""""""""
        DO IP=1,NPOIN2
          OMOLD(IP) = GRAVIT * Z0OLD(IP) / CPHAS(IP)**2.D0
          OMNEW(IP) = GRAVIT * Z0NEW(IP) / CPHAS(IP)**2.D0
          USO(IP) = (USOLD(IP) / CPHAS(IP)) + DECAL
          USN(IP) = (USNEW(IP) / CPHAS(IP)) + DECAL
        ENDDO
C
C.......BOUCLE SUR LES DIRECTIONS DE DISCRETISATION.
C       """"""""""""""""""""""""""""""""""""""""""""
        DO JP=1,NPLAN
C
          DO IP=1,NPOIN2
            BETAO(IP)=0.D0
            BETAN(IP)=0.D0
          ENDDO
C
C.........CALCUL DU TERME SOURCE
C         """"""""""""""""""""""
          DO IP=1,NPOIN2
            IF (TOLD(IP,JP).GT.0.01D0) THEN
              XX = USO(IP) * TOLD(IP,JP)
              ZLOGMU = DLOG(OMOLD(IP)) + XKAPPA/XX
              IF (ZLOGMU.LT.0.D0) THEN
                BETAO(IP) = CONST*OMOLD(IP)*EXP(XKAPPA/XX)*
     *                          ZLOGMU**4.D0*XX**2.D0
              ENDIF
            ENDIF
          ENDDO
          DO IP=1,NPOIN2
            IF (TNEW(IP,JP).GT.0.01D0) THEN
              XX = USN(IP) * TNEW(IP,JP)
              ZLOGMU = DLOG(OMNEW(IP)) + XKAPPA/XX
              IF (ZLOGMU.LT.0.D0) THEN
                BETAN(IP) = CONST*OMNEW(IP)*EXP(XKAPPA/XX)*
     *                          ZLOGMU**4.D0*XX**2.D0
              ENDIF
            ENDIF
          ENDDO
C
C.........PRISE EN COMPTE DU TERME SOURCE.
C         """"""""""""""""""""""""""""""""
          DO IP=1,NPOIN2
            TSTOT(IP,JP,JF) = TSTOT(IP,JP,JF)
     *              + (DIMPLI*BETAO(IP)+CIMPLI*BETAN(IP)) * F(IP,JP,JF)
            TSDER(IP,JP,JF) = TSDER(IP,JP,JF) + BETAN(IP)
          ENDDO
C
        ENDDO
      ENDDO
C
      RETURN
      END
