C                       *****************
                        SUBROUTINE QBREK2
C                       *****************
C
     *( TSTOT , TSDER , F     , FCAR  , VARIAN, DEPTH , BORETG, GAMATG,
     *  IWHTG , NF    , NPLAN , NPOIN2, BETA  )
C
C**********************************************************************
C  TOMAWAC - V1.1    F. BECQ & M. BENOIT     (EDF/DER/LNH)  -  14/02/96
C**********************************************************************
C
C  FONCTION : CALCUL DE LA CONTRIBUTION DU TERME SOURCE DE DEFERLEMENT
C  ********** DU A LA BATHYMETRIE BASE SUR LA METHODE PROPOSEE PAR
C             THORNTON ET GUZA (1983).
C
C  ARGUMENTS :
C  ***********
C +-------------+----+--------------------------------------------+
C ! NOM         !MODE! SIGNIFICATION - OBSERVATIONS               !
C +-------------+----+--------------------------------------------+
C ! TSTOT(-,-,-)!<-->! CONTRIBUTION TERME SOURCE - PARTIE TOTALE  !
C ! TSDER(-,-,-)!<-->! CONTRIBUTION TERME SOURCE - PARTIE DERIVEE !
C ! F(-,-,-)    ! -->! SPECTRE DIRECTIONNEL                       !
C ! FCAR(-)     ! -->! FREQUENCE CARACTERISTIQUE DU SPECTRE       !
C ! VARIAN(-)   ! -->! VARIANCE TOTALE DU SPECTRE                 !
C ! DEPTH(-)    ! -->! TABLEAU DES PROFONDEURS                    !
C ! BORETG      ! -->! MODELE DEFERLEMENT TG : CONSTANTE B        !
C ! GAMATG      ! -->! MODELE DEFERLEMENT TG : CONSTANTE GAMMA    !
C ! IWHTG       ! -->! MODELE DEFERLEMENT TG : MODE CALCUL DE W(H)!
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
C  REFERENCES : - THORNTON ET GUZA (1983) : TRANSFORMATION OF WAVE
C  ************   HEIGHT DISTRIBUTION.
C
C**********************************************************************
C
      IMPLICIT NONE
C
C.....VARIABLES TRANSMISES
C     """"""""""""""""""""
      INTEGER  NF    , NPLAN , NPOIN2, IWHTG
      DOUBLE PRECISION BORETG, GAMATG
      DOUBLE PRECISION DEPTH(NPOIN2), BETA(NPOIN2), FCAR(NPOIN2)
      DOUBLE PRECISION TSTOT(NPOIN2,NPLAN,NF), TSDER(NPOIN2,NPLAN,NF)
      DOUBLE PRECISION     F(NPOIN2,NPLAN,NF), VARIAN(NPOIN2)
C
C.....VARIABLES LOCALES
C     """""""""""""""""
      INTEGER  JP    , IFF   , IP
      DOUBLE PRECISION COEF  , GAMMA2, SEUIL , DEUPI
C
C
      DEUPI  = 6.283185307D0
      SEUIL  = 1.D-6
      GAMMA2 = GAMATG*GAMATG
      COEF   = -24.D0*DSQRT(DEUPI)*BORETG**3/GAMMA2
C
      IF (IWHTG.EQ.1) THEN
C
C.......CALCUL DU COEFFICIENT LINEAIRE BETA : QBREK2 = BETA * F
C       AVEC FONCTION DE POIDS W(H) = CONSTANTE
C       """""""""""""""""""""""""""""""""""""""""""""""""""""""
        DO 25 IP = 1,NPOIN2
          BETA(IP) = COEF*8.D0*VARIAN(IP)**2.5D0*FCAR(IP)
     &             /(GAMMA2*DEPTH(IP)**5)
   25   CONTINUE
C
      ELSEIF (IWHTG.EQ.2) THEN
C
C.......CALCUL DU COEFFICIENT LINEAIRE BETA : QBREK2 = BETA * F
C       AVEC FONCTION DE POIDS W(H) != CONSTANTE
C       """""""""""""""""""""""""""""""""""""""""""""""""""""""
        DO 35 IP = 1,NPOIN2
          BETA(IP) = (COEF*VARIAN(IP)**1.5D0*FCAR(IP)/
     &                DEPTH(IP)**3)*(1.D0-1.D0/(1.D0+VARIAN(IP)*8.D0
     &                /(GAMMA2*DEPTH(IP)*DEPTH(IP)))**2.5D0)
  35    CONTINUE
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
