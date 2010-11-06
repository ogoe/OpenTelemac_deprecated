C                       *****************
                        SUBROUTINE QBREK3
C                       *****************
C
     *( TSTOT , TSDER , F     , FCAR  , VARIAN, DEPTH , ALFARO, GAMARO,
     *  GAM2RO, IEXPRO, IDISRO, NF    , NPLAN , NPOIN2, BETA  )
C
C**********************************************************************
C  TOMAWAC - V1.1    F. BECQ & M. BENOIT     (EDF/DER/LNH)  -  26/03/96
C**********************************************************************
C
C  FONCTION : CALCUL DE LA CONTRIBUTION DU TERME SOURCE DE DEFERLEMENT
C  ********** DU A LA BATHYMETRIE BASE SUR LA METHODE PROPOSEE PAR
C             ROELVINK (1993).
C
C  ARGUMENTS :
C  ***********
C +-------------+----+--------------------------------------------+
C ! NOM         !MODE! SIGNIFICATION - OBSERVATIONS               !
C +-------------+----+--------------------------------------------+
C ! TSTOT(-,-,-)!<-->! CONTRIBUTION TERME SOURCE - PARTIE TOTALE  !
C ! TSDER(-,-,-)!<-->! CONTRIBUTION TERME SOURCE - PARTIE DERIVEE !
C ! F(-,-,-)    ! -->! SPECTRE DIRECTIONNEL                       !
C ! FCAR(-)     ! -->! FREQUENCE CARACTERISTIQUE                  !
C ! VARIAN(-)   ! -->! VARIANCE TOTALE DU SPECTRE                 !
C ! DEPTH(-)    ! -->! TABLEAU DES PROFONDEURS                    !
C ! ALFARO      ! -->! CONSTANTE ALPHA  DU MODELE RO              !
C ! GAMARO      ! -->! CONSTANTE GAMMA  DU MODELE RO              !
C ! GAM2RO      ! -->! CONSTANTE GAMMA2 DU MODELE RO              !
C ! IEXPRO      ! -->! EXPOSANT N DU MODELE RO                    !
C ! IDISRO      ! -->! CHOIX DE LA DISTRIBUTION DES HAUTEURS      !
C ! NF          ! -->! NOMBRE DE FREQUENCES DE DISCRETISATION     !
C ! NPLAN       ! -->! NOMBRE DE DIRECTIONS DE DISCRETISATION     !
C ! NPOIN2      ! -->! NOMBRE DE POINTS DU MAILLAGE SPATIAL       !
C ! BETA(-)     !<-->! TABLEAU DE TRAVAIL (DIMENSION NPOIN2)      !
C +-------------+----+--------------------------------------------+
C ! MODE   (-> : NON-MODIFIE)  (<-> : MODIFIE)  (<- : INITIALISE) !
C +---------------------------------------------------------------+
C
C  APPELS :    - PROGRAMME(S) APPELANT  :  SEMIMP
C  ********    - PROGRAMME(S) APPELE(S) :  GAMMLN
C                                          QGAUSS
C
C  REMARQUES :
C  ***********
C  - CE TERME SOURCE EST LINEAIRE EN F(FREQ,TETA) ET LE COEFFICIENT
C    LINEAIRE EST CONSTANT DANS LE TEMPS.
C
C  REFERENCES : - ROELVINK (1993) : DISSIPATION IN RANDOM WAVE GROUPS
C  ************   INCIDENT ON A BEACH. COASTAL ENG. VOL 19, PP 127-150.
C
C**********************************************************************
C
      IMPLICIT NONE
C
C.....VARIABLES TRANSMISES
C     """"""""""""""""""""
      INTEGER          NF    , NPLAN , NPOIN2, IEXPRO, IDISRO
      DOUBLE PRECISION ALFARO, GAMARO, GAM2RO
      DOUBLE PRECISION DEPTH(NPOIN2), BETA(NPOIN2), FCAR(NPOIN2)
      DOUBLE PRECISION TSTOT(NPOIN2,NPLAN,NF), TSDER(NPOIN2,NPLAN,NF)
      DOUBLE PRECISION     F(NPOIN2,NPLAN,NF), VARIAN(NPOIN2)
C
C.....VARIABLES LOCALES
C     """""""""""""""""
      INTEGER  JP    , IFF   , IP
      DOUBLE PRECISION COEF1 , COEF2 , SEUIL , PIS2  , DEUPI
      DOUBLE PRECISION A     , XM    , SIGMA , BX    , FN
C
C.....FONCTIONS EXTERNES
C     """"""""""""""""""
      DOUBLE PRECISION   GAMMLN, QGAUSS
      EXTERNAL           GAMMLN, QGAUSS
C
      PARAMETER (PIS2 = 1.570796327D0 , DEUPI = 6.283185307D0)
C
C
      SEUIL  = 1.D-6
      COEF1  = -2.D0*ALFARO
      COEF2  = 8.D0/(GAMARO*GAMARO)
C
      IF (IDISRO.EQ.1) THEN
C
C.......CALCUL DU COEFFICIENT LINEAIRE BETA (DISTR. WEIBULL).
C       """""""""""""""""""""""""""""""""""""""""""""""""""""
        DO 40 IP = 1,NPOIN2
          IF (VARIAN(IP).GT.SEUIL) THEN
            BX    = COEF2*VARIAN(IP)/(DEPTH(IP)*DEPTH(IP))
            SIGMA = DSQRT(8.D0*VARIAN(IP))/DEPTH(IP)
            XM    = 1.D0 + 0.7D0*(DTAN(PIS2*SIGMA/GAM2RO))**2
            A     = DEXP(XM*(GAMMLN(1.D0+1.D0/XM,DEUPI)))
            IF (XM.GT.98.D0) THEN
               FN = 1.D0
            ELSE
               FN = QGAUSS(BX,IEXPRO,A,XM)
            ENDIF
            BETA(IP) = COEF1*FCAR(IP)*FN
          ELSE
            BETA(IP) = 0.D0
          ENDIF
   40   CONTINUE
C
      ELSE
C
C.......CALCUL DU COEFFICIENT LINEAIRE BETA (DISTR. RAYLEIGH).
C       """"""""""""""""""""""""""""""""""""""""""""""""""""""
        DO 50 IP = 1,NPOIN2
          BX = COEF2*VARIAN(IP)/(DEPTH(IP)*DEPTH(IP))
          XM = 1.D0
          A  = 1.D0
          FN = QGAUSS(BX,IEXPRO,A,XM)
          BETA(IP) = COEF1*FCAR(IP)*FN
   50   CONTINUE
      ENDIF
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
