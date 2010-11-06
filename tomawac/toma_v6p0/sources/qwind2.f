C                       *****************
                        SUBROUTINE QWIND2
C                       *****************
C
     *( TSTOT , TSDER , F     , XK    , FREQ  , USOLD , USNEW , TWOLD ,
     *  TWNEW , TETA  , ROAIR , ROEAU , GRAVIT, NF    , NPLAN , NPOIN2,
     *  CIMPLI, CPHAS , USN   , USO   , BETAN , BETAO )
C
C**********************************************************************
C  TOMAWAC - V1P1    M. BENOIT               (EDF/DER/LNH)  -  26/03/96
C**********************************************************************
C
C  FONCTION : CALCUL DE LA CONTRIBUTION DU TERME SOURCE DE GENERATION
C  ********** PAR LE VENT. UTILISATION DE LA THEORIE DE SNYDER ET AL.
C             (1981), MODIFIEE PAR KOMEN ET AL. (1984) POUR UTILISER
C             LA VITESSE DE FROTTEMENT U* AU LIEU DE LA VITESSE U5
C             (MESUREE A 5 METRES) POUR LE VENT.
C             CETTE THEORIE DE GENERATION EST IDENTIQUE A WAM-CYCLE 3.
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
C ! TETA(-)     ! -->! TABLEAU DES DIRECTIONS DE DISCRETISATION   !
C ! ROAIR       ! -->! DENSITE DE L'AIR                           !
C ! ROEAU       ! -->! DENSITE DE L'EAU                           !
C ! GRAVIT      ! -->! ACCELERATION DE LA PESANTEUR               !
C ! NF          ! -->! NOMBRE DE FREQUENCES DE DISCRETISATION     !
C ! NPLAN       ! -->! NOMBRE DE DIRECTIONS DE DISCRETISATION     !
C ! NPOIN2      ! -->! NOMBRE DE POINTS DU MAILLAGE SPATIAL       !
C ! CPHAS(-)    !<-->! TABLEAU DE TRAVAIL (DIMENSION NPOIN2)      !
C ! USN(-)      !<-->! TABLEAU DE TRAVAIL (DIMENSION NPOIN2)      !
C ! USO(-)      !<-->! TABLEAU DE TRAVAIL (DIMENSION NPOIN2)      !
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
C  REFRENCES : - SNYDER ET AL. (1981) : ARRAY MEASUREMENTS OF ATMOSPHE-
C              RIC PRESSURE FLUCTUATIONS ABOVE SURFACE GRAVITY WAVES.
C              JOURNAL OF FLUID MECH., VOL 102., PP 1-59.
C              - KOMEN ET AL.  (1984) : ON THE EXISTENCE OF A FULLY
C              DEVELOPED WINDSEA SPECTRUM. JPO, VOL 14, PP 1271-1285.
C
C**********************************************************************
C
      IMPLICIT NONE
C
C.....VARIABLES TRANSMISES
C     """"""""""""""""""""
      INTEGER  NF    , NPLAN , NPOIN2
      DOUBLE PRECISION GRAVIT, ROAIR , ROEAU , CIMPLI
      DOUBLE PRECISION  FREQ(NF)    , TETA(NPLAN)
      DOUBLE PRECISION TWOLD(NPOIN2), TWNEW(NPOIN2), USNEW(NPOIN2)
      DOUBLE PRECISION BETAO(NPOIN2), BETAN(NPOIN2), USOLD(NPOIN2)
      DOUBLE PRECISION CPHAS(NPOIN2),   USO(NPOIN2),   USN(NPOIN2)
      DOUBLE PRECISION TSTOT(NPOIN2,NPLAN,NF), TSDER(NPOIN2,NPLAN,NF)
      DOUBLE PRECISION     F(NPOIN2,NPLAN,NF),    XK(NPOIN2,NF)
C
C.....VARIABLES LOCALES
C     """""""""""""""""
      INTEGER  JP    , JF    , IP
      DOUBLE PRECISION DEUPI , C1    , DIREC , CONST , DIMPLI
C
C
      DEUPI = 2.D0* 3.14159265358978D0
      C1 = 0.25D0 * (ROAIR/ROEAU) * DEUPI
      DIMPLI=1.0D0-CIMPLI
C
C.....BOUCLE SUR LES DIRECTIONS DE DISCRETISATION.
C     """"""""""""""""""""""""""""""""""""""""""""
      DO JP=1,NPLAN
C
C.......PRECALCULS SUR LES DEPENDANCES DIRECTIONNELLES
C       """"""""""""""""""""""""""""""""""""""""""""""
        DIREC=TETA(JP)
        DO IP=1,NPOIN2
          USO(IP)=USOLD(IP)*COS(DIREC-TWOLD(IP))
          USN(IP)=USNEW(IP)*COS(DIREC-TWNEW(IP))
        ENDDO
C
C.......BOUCLE SUR LES FREQUENCES DE DISCRETISATION.
C       """"""""""""""""""""""""""""""""""""""""""""
        DO JF=1,NF
          CONST=C1*FREQ(JF)
C
C.........CALCUL DES COEFFICIENTS DE PROPORTIONALITE BETA.
C         """"""""""""""""""""""""""""""""""""""""""""""""
          DO IP=1,NPOIN2
            CPHAS(IP) = DEUPI * FREQ(JF) / XK(IP,JF)
            BETAO(IP)=MAX(28.D0*USO(IP)/CPHAS(IP)-1.D0,0.D0)*CONST
            BETAN(IP)=MAX(28.D0*USN(IP)/CPHAS(IP)-1.D0,0.D0)*CONST
          ENDDO
C
C.........PRISE EN COMPTE DU TERME SOURCE.
C         """"""""""""""""""""""""""""""""
          DO IP=1,NPOIN2
            TSTOT(IP,JP,JF) = TSTOT(IP,JP,JF)
     *             + (DIMPLI*BETAO(IP)+CIMPLI*BETAN(IP)) * F(IP,JP,JF)
            TSDER(IP,JP,JF) = TSDER(IP,JP,JF) + BETAN(IP)
          ENDDO
C
        ENDDO
C
      ENDDO
C
      RETURN
      END
