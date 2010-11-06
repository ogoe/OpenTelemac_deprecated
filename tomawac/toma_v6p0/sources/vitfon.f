C                       *****************
                        SUBROUTINE VITFON
C                       *****************
C
     *( UWBM  , F     , XK    , DEPTH , DFREQ , NF    , NPOIN2, NPLAN ,
     *      GRAVIT, BETA  )
C
C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C
C                                                                     C
C  EDF/DER - LABORATOIRE NATIONAL D'HYDRAULIQUE  -  CHATOU (FRANCE)   C
C                                                                     C
C  CODE TOMAWAC DE MODELISATION DES ETATS DE MER EN ELEMENTS FINIS    C
C                       ---- VERSION 1.2 ----                         C
C                                                                     C
C  VITFON :  CALCUL DE LA VITESSE ORBITALE MAXIMALE SUR LE FOND       C
C  ********  (VITESSE MOYENNE SUR LE SPECTRE)                         C
C                                                                     C
C   - CREE POUR VERSION 1.2  LE 05/07/96 PAR M. BENOIT                C
C                                                                     C
C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C
C                                                                     C
C  ARGUMENTS D'APPEL :                                                C
C  *******************                                                C
C  +-------------+----+--------------------------------------------+  C
C  ! NOM         !MODE! SIGNIFICATION - OBSERVATIONS               !  C
C  +-------------+----+--------------------------------------------+  C
C  ! UWBM(-)     !<-- ! TABLEAU DES VITESSES ORBITALES MAX AU FOND !  C
C  ! F(-,-,-)    ! -->! SPECTRE DIRECTIONNEL DE VARIANCE           !  C
C  ! XK(-,-)     ! -->! TABLEAU DES NOMBRES D'ONDE                 !  C
C  ! DEPTH(-)    ! -->! TABLEAU DES PROFONDEURS (METRES)           !  C
C  ! DFREQ(-)    ! -->! TABLEAU DES PAS DE FREQUENCES              !  C
C  ! NF          ! -->! NOMBRE DE FREQUENCES DE DISCRETISATION     !  C
C  ! NPOIN2      ! -->! NOMBRE DE POINTS DU MAILLAGE SPATIAL 2D    !  C
C  ! NPLAN       ! -->! NOMBRE DE DIRECTIONS DE DISCRETISATION     !  C
C  ! GRAVIT      ! -->! ACCELERATION DE LA PESANTEUR               !  C
C  ! BETA(-)     !<-->! TABLEAU DE TRAVAIL (DIMENSION NPOIN2)      !  C
C  +-------------+----+--------------------------------------------+  C
C  ! MODE   (-> : NON-MODIFIE)  (<-> : MODIFIE)  (<- : INITIALISE) !  C
C  +---------------------------------------------------------------+  C
C                                                                     C
C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C
C                                                                     C
C  APPELS :    - PROGRAMME(S) APPELANT  :  DUMP2D                     C
C  ********    - PROGRAMME(S) APPELE(S) :                             C
C                                                                     C
C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C
C
      IMPLICIT NONE
C
C.....VARIABLES TRANSMISES
C     """"""""""""""""""""
      INTEGER  NF    , NPLAN , NPOIN2
      DOUBLE PRECISION GRAVIT
      DOUBLE PRECISION UWBM(NPOIN2), DEPTH(NPOIN2), BETA(NPOIN2)
      DOUBLE PRECISION F(NPOIN2,NPLAN,NF), XK(NPOIN2,NF), DFREQ(NF)
C
C.....VARIABLES LOCALES
C     """""""""""""""""
      INTEGER  IP    , JP    , JF
      DOUBLE PRECISION DEUPI , DTETAR, DEUKD , COEF
C
C
      DEUPI=2.D0*3.14159265D0
      DTETAR=DEUPI/FLOAT(NPLAN)
C
      DO 30 IP = 1,NPOIN2
        UWBM(IP) = 0.D0
   30 CONTINUE
C
C.....SOMMATIONS SUR LA PARTIE DISCRETISEE DU SPECTRE.
C     """"""""""""""""""""""""""""""""""""""""""""""""
      DO 20 JF = 1,NF
        COEF=2.D0*GRAVIT*DFREQ(JF)*DTETAR
        DO 25 IP = 1,NPOIN2
          DEUKD = MIN(2.D0*DEPTH(IP)*XK(IP,JF),7.D2)
          BETA(IP) = COEF*XK(IP,JF)/SINH(DEUKD)
   25   CONTINUE
        DO 10 JP = 1,NPLAN
          DO 5 IP=1,NPOIN2
            UWBM(IP) = UWBM(IP) + F(IP,JP,JF)*BETA(IP)
    5     CONTINUE
   10   CONTINUE
   20 CONTINUE
C
      DO 35 IP=1,NPOIN2
        UWBM(IP) = DSQRT(UWBM(IP))
   35 CONTINUE
C
      RETURN
      END
