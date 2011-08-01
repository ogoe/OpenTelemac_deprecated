C                       *****************
                        SUBROUTINE QFROT1
C                       *****************
C
     *( TSTOT , TSDER , F     , XK    , DEPTH , CFROT1, GRAVIT, NF    ,
     *  NPLAN , NPOIN2, BETA  )
C
C**********************************************************************
C  TOMAWAC - V1.0    P. THELLIER & M. BENOIT (EDF/DER/LNH)  -  03/04/95
C**********************************************************************
C
C  FONCTION : CALCUL DE LA CONTRIBUTION DU TERME SOURCE DE FROTTEMENT
C  ********** SUR LE FOND BASE SUR L'EXPRESSION PROPOSEE PAR HASSELMANN
C             ET AL. (1973) ET MODIFIEE PAR BOUWS ET KOMEN (1983).
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
C ! DEPTH(-)    ! -->! TABLEAU DES PROFONDEURS                    !
C ! CFROT1      ! -->! CONSTANTE DE L'EXPRESSION DE FROTTEMENT    !
C ! GRAVIT      ! -->! ACCELERATION DE LA PESANTEUR               !
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
C  - LA CONSTANTE CFROT1 UTILISEE DANS WAM CYCLE 4 VAUT 0.038 M2.S-3
C
C  REFERENCES : - HASSELMANN ET AL. (1973) : MEASUREMENTS OF WIND-WAVE
C               GROWTH AND SWELL DECAY DURING THE JOINT NORTH SEA
C               WAVE PROJECT (JONSWAP). DEUTSCHEN HYDROGRAPHISVHEN
C               ZEITSCHRIFT, REIHE A(8), NUM 12.
C               - BOUWS E., KOMEN G.J. (1983) : ON THE BALANCE BETWEEN
C               GROWTH AND DISSIPATION IN AN EXTREME DEPTH-LIMITED
C               WIND-SEA IN THE SOUTHERN NORTH-SEA. JPO, VOL 13.
C
C**********************************************************************
C
      IMPLICIT NONE
C
C.....VARIABLES TRANSMISES
C     """"""""""""""""""""
      INTEGER  NF    , NPLAN , NPOIN2 
      DOUBLE PRECISION CFROT1, GRAVIT
      DOUBLE PRECISION DEPTH(NPOIN2), BETA(NPOIN2)
      DOUBLE PRECISION TSTOT(NPOIN2,NPLAN,NF), TSDER(NPOIN2,NPLAN,NF)
      DOUBLE PRECISION     F(NPOIN2,NPLAN,NF),    XK(NPOIN2,NF)
C
C.....VARIABLES LOCALES
C     """""""""""""""""
      INTEGER  JP    , JF    , IP
      DOUBLE PRECISION COEF , DEUKD
C
C
      COEF=-2.D0*CFROT1/GRAVIT
C
C.....BOUCLE SUR LES FREQUENCES DE DISCRETISATION.
C     """"""""""""""""""""""""""""""""""""""""""""
      DO JF=1,NF
C
C.......CALCUL DU COEFFICIENT LINEAIRE BETA : QFROT1 = BETA * F
C       """""""""""""""""""""""""""""""""""""""""""""""""""""""
        DO IP=1,NPOIN2
          DEUKD = MIN(2.D0*DEPTH(IP)*XK(IP,JF),7.D2)
          BETA(IP) = COEF*XK(IP,JF)/SINH(DEUKD)
        ENDDO
C
C.......PRISE EN COMPTE DU TERME SOURCE.
C       """"""""""""""""""""""""""""""""
        DO JP=1,NPLAN
          DO IP=1,NPOIN2
            TSTOT(IP,JP,JF) = TSTOT(IP,JP,JF)+BETA(IP)*F(IP,JP,JF)
            TSDER(IP,JP,JF) = TSDER(IP,JP,JF)+BETA(IP)
          ENDDO
        ENDDO
      ENDDO
C
      RETURN
      END
