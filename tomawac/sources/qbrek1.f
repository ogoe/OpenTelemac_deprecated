C                       *****************
                        SUBROUTINE QBREK1
C                       *****************
C
     *( TSTOT , TSDER , F     , FCAR  , VARIAN, DEPTH , ALFABJ, GAMBJ1,
     *  GAMBJ2, IQBBJ , IHMBJ , NF    , NPLAN , NPOIN2, BETA  )
C
C**********************************************************************
C  TOMAWAC - V1.1        F. BECQ & M. BENOIT (EDF/DER/LNH)  -  14/02/96
C  TOMAWAC - V5.2        OPTIMER                            -  14/06/01
C**********************************************************************
C
C  FONCTION : CALCUL DE LA CONTRIBUTION DU TERME SOURCE DE DEFERLEMENT
C  ********** DU A LA BATHYMETRIE BASE SUR LA METHODE PROPOSEE PAR
C             BATTJES ET JANSSEN (1978).
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
C ! ALFABJ      ! -->! CONSTANTE ALPHA DU MODELE BJ               !
C ! GAMBJ2      ! -->! CONSTANTE GAMMA DU MODELE BJ               !
C ! IQBBJ       ! -->! MODE DE CALCUL DE QB POUR MODELE BJ        !
C ! IHMBJ       ! -->! TYPE DE HAUTEUR DE HOULE MAX POUR MODELE BJ!
C ! NF          ! -->! NOMBRE DE FREQUENCES DE DISCRETISATION     !
C ! NPLAN       ! -->! NOMBRE DE DIRECTIONS DE DISCRETISATION     !
C ! NPOIN2      ! -->! NOMBRE DE POINTS DU MAILLAGE SPATIAL       !
C ! BETA(-)     !<-->! TABLEAU DE TRAVAIL (DIMENSION NPOIN2)      !
C +-------------+----+--------------------------------------------+
C ! MODE   (-> : NON-MODIFIE)  (<-> : MODIFIE)  (<- : INITIALISE) !
C +---------------------------------------------------------------+
C
C  APPELS :    - PROGRAMME(S) APPELANT  :  SEMIMP
C  ********    - PROGRAMME(S) APPELE(S) :  QBBJ78
C                                          WNSCOU
C
C  REMARQUES :
C  ***********
C  - CE TERME SOURCE EST LINEAIRE EN F(FREQ,TETA) ET LE COEFFICIENT
C    LINEAIRE EST CONSTANT DANS LE TEMPS.
C
C  REFERENCES : - BATTJES ET JANSSEN (1978) : ENERGY LOSS AND SET-UP
C  ************   DUE TO BREAKING OF RANDOM WAVES. (ICCE'78)
C
C**********************************************************************
C
      IMPLICIT NONE
C
C.....VARIABLES TRANSMISES
C     """"""""""""""""""""
      INTEGER          NF    , NPLAN  , NPOIN2, IQBBJ , IHMBJ
      DOUBLE PRECISION ALFABJ, GAMBJ1 , GAMBJ2
      DOUBLE PRECISION DEPTH(NPOIN2), BETA(NPOIN2)
      DOUBLE PRECISION TSTOT(NPOIN2,NPLAN,NF), TSDER(NPOIN2,NPLAN,NF)
      DOUBLE PRECISION     F(NPOIN2,NPLAN,NF), VARIAN(NPOIN2)
      DOUBLE PRECISION     FCAR(NPOIN2)
C
C.....VARIABLES LOCALES
C     """""""""""""""""
      INTEGER          JP   , IFF , IP
      DOUBLE PRECISION COEF , HM  , XK8 , XKCAR , B , QB , SEUIL
C
C.....FONCTIONS EXTERNES
C     """"""""""""""""""
      DOUBLE PRECISION   QBBJ78
      EXTERNAL           QBBJ78
C
C
      SEUIL=1.D-6
      COEF =-.25D0*ALFABJ
C
C.....CALCUL DU COEFFICIENT LINEAIRE BETA : QBREK1 = BETA * F
C     """""""""""""""""""""""""""""""""""""""""""""""""""""""
      DO 25 IP = 1,NPOIN2
         IF (VARIAN(IP).GT.SEUIL) THEN
C
C..........CALCUL DE LA HAUTEUR DE HOULE MAXIMALE.
C          """""""""""""""""""""""""""""""""""""""
           IF (IHMBJ.EQ.1) THEN
             HM  = GAMBJ2*DEPTH(IP)
           ELSEIF (IHMBJ.EQ.2) THEN
             CALL WNSCOU(XKCAR,FCAR(IP),DEPTH(IP))
             XK8 = GAMBJ1/XKCAR
             HM  = XK8*DTANH(GAMBJ2*DEPTH(IP)/XK8)
           ENDIF
C
C..........CALCUL DE LA FRACTION DE VAGUES DEFERLANTES.
C          """"""""""""""""""""""""""""""""""""""""""""
           B   = DSQRT(8.D0*VARIAN(IP))/HM
           QB  = QBBJ78(B,IQBBJ)
C
           BETA(IP) = COEF*QB*FCAR(IP)*HM*HM/VARIAN(IP)
         ELSE
           BETA(IP) = 0.D0
         ENDIF
   25 CONTINUE
C
C.....PRISE EN COMPTE DU TERME SOURCE.
C     """"""""""""""""""""""""""""""""
      DO 10 IFF = 1,NF
        DO 20 JP = 1,NPLAN
          DO 30 IP = 1,NPOIN2
            TSTOT(IP,JP,IFF) = TSTOT(IP,JP,IFF)+BETA(IP)*F(IP,JP,IFF)
C           TSDER(IP,JP,IFF) = TSDER(IP,JP,IFF)+BETA(IP)
   30     CONTINUE
   20   CONTINUE
   10 CONTINUE
C
      RETURN
      END
