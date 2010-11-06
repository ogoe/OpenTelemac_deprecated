C                       *****************
                        SUBROUTINE QBREK4
C                       *****************
C
     *( TSTOT , TSDER , F     , FCAR  , VARIAN, DEPTH , BETAIH, EM2SIH,
     *  GRAVIT, NF    , NPLAN , NPOIN2, BETA  )
C
C**********************************************************************
C  TOMAWAC - V1.1    F. BECQ & M. BENOIT     (EDF/DER/LNH)  -  26/03/96
C**********************************************************************
C
C  FONCTION : CALCUL DE LA CONTRIBUTION DU TERME SOURCE DE DEFERLEMENT
C  ********** DU A LA BATHYMETRIE BASE SUR LA METHODE PROPOSEE PAR
C             IZUMIYA ET HORIKAWA (1984).
C
C  ARGUMENTS :
C +-------------+----+--------------------------------------------+
C ! NOM         !MODE! SIGNIFICATION - OBSERVATIONS               !
C +-------------+----+--------------------------------------------+
C ! TSTOT(-,-,-)!<-->! CONTRIBUTION TERME SOURCE - PARTIE TOTALE  !
C ! TSDER(-,-,-)!<-->! CONTRIBUTION TERME SOURCE - PARTIE DERIVEE !
C ! F(-,-,-)    ! -->! SPECTRE DIRECTIONNEL                       !
C ! FCAR(-)     ! -->! FREQUENCE CARACTERISTIQUE                  !
C ! VARIAN(-)   ! -->! VARIANCE TOTALE DU SPECTRE                 !
C ! DEPTH(-)    ! -->! TABLEAU DES PROFONDEURS                    !
C ! BETAIH      ! -->! CONSTANTE BETA0 DU MODELE IH               !
C ! EM2SIH      ! -->! CONSTANTE M2*   DU MODELE IH               !
C ! NF          ! -->! NOMBRE DE FREQUENCES DE DISCRETISATION     !
C ! NPLAN       ! -->! NOMBRE DE DIRECTIONS DE DISCRETISATION     !
C ! NPOIN2      ! -->! NOMBRE DE POINTS DU MAILLAGE SPATIAL       !
C ! BETA(-)     !<-->! TABLEAU DE TRAVAIL (DIMENSION NPOIN2)      !
C +-------------+----+--------------------------------------------+
C ! MODE   (-> : NON-MODIFIE)  (<-> : MODIFIE)  (<- : INITIALISE) !
C +---------------------------------------------------------------+
C
C  APPELS :    - PROGRAMME(S) APPELANT  :  SEMIMP
C  ********    - PROGRAMME(S) APPELE(S) :    -
C
C  REMARQUES :
C  ***********
C  - CE TERME SOURCE EST LINEAIRE EN F(FREQ,TETA) ET LE COEFFICIENT
C    LINEAIRE EST CONSTANT DANS LE TEMPS.
C
C  REFERENCES : - IZUMIYA T., HORIKAWA K. (1984) : WAVE ENERGY EQUATION
C  ************   APPLICABLE IN AND OUTSIDE THE SURF ZONE. COASTAL
C                 ENGINEERING IN JAPAN, VOL 17, PP 119-137.
C
C**********************************************************************
C
      IMPLICIT NONE
C
C.....VARIABLES TRANSMISES
C     """"""""""""""""""""
      INTEGER          NF    , NPLAN , NPOIN2
      DOUBLE PRECISION BETAIH, EM2SIH, GRAVIT
      DOUBLE PRECISION DEPTH(NPOIN2), BETA(NPOIN2), FCAR(NPOIN2)
      DOUBLE PRECISION TSTOT(NPOIN2,NPLAN,NF), TSDER(NPOIN2,NPLAN,NF)
      DOUBLE PRECISION     F(NPOIN2,NPLAN,NF), VARIAN(NPOIN2)
C
C.....VARIABLES LOCALES
C     """""""""""""""""
      INTEGER  JP    , IFF   , IP
      DOUBLE PRECISION COEF  , XKCAR , DEUKD , GG1   , GG2
C
C
      COEF   = -DSQRT(GRAVIT)*BETAIH
C
C.....CALCUL DU COEFFICIENT LINEAIRE BETA : QBREK4 = BETA * F
C     """""""""""""""""""""""""""""""""""""""""""""""""""""""
      DO 40 IP = 1,NPOIN2
        CALL WNSCOU( XKCAR, FCAR(IP), DEPTH(IP) )
        DEUKD=2.D0*XKCAR*DEPTH(IP)
        IF (DEUKD.GT.7.D2) THEN
          GG1 = 0.D0
          GG2 = 0.5D0
        ELSE
          GG1 = DEUKD/SINH(DEUKD)
          GG2 = 0.5D0*(1.D0+GG1)
        ENDIF
        BETA(IP) = COEF/DEPTH(IP)**1.5*DSQRT(VARIAN(IP)*GG1)
     *               *DSQRT(DMAX1(0.D0,GG2*VARIAN(IP)
     *               /(DEPTH(IP)*DEPTH(IP))-EM2SIH))
   40 CONTINUE
C
C
C.....PRISE EN COMPTE DU TERME SOURCE.
C     """"""""""""""""""""""""""""""""
      DO 10 IFF = 1,NF
        DO 20 JP = 1,NPLAN
          DO 30 IP = 1,NPOIN2
            TSTOT(IP,JP,IFF) = TSTOT(IP,JP,IFF)+BETA(IP)*F(IP,JP,IFF)
C            TSDER(IP,JP,IFF) = TSDER(IP,JP,IFF)+BETA(IP)
   30     CONTINUE
   20   CONTINUE
   10 CONTINUE
C
      RETURN
      END
