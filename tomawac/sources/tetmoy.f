C                       *****************
                        SUBROUTINE TETMOY
C                       *****************
C
     *( TETAM , F     , COSTET, SINTET, NPLAN , FREQ  , DFREQ , NF    ,
     *  NPOIN2, TAILF , COSMOY, SINMOY, TAUXC , TAUXS )
C
C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C
C                                                                     C
C  EDF/DER - LABORATOIRE NATIONAL D'HYDRAULIQUE  -  CHATOU (FRANCE)   C
C                                                                     C
C  CODE TOMAWAC DE MODELISATION DES ETATS DE MER EN ELEMENTS FINIS    C
C                       ---- VERSION 1.2 ----                         C
C                                                                     C
C  TETMOY : CALCUL DE LA DIRECTION MOYENNE D'UN SPECTRE DIRECTIONNEL  C
C  ******** (EN CALCULANT L'ARCTANGENTE DES SINUS ET COSINUS MOYENS)  C
C           LE RESULTAT EST DONNE EN RADIANS.                         C
C                                                                     C
C   - CREE POUR VERSION 1.0  LE 01/02/95 PAR P. THELIIER ET M. BENOIT C
C   - MOD. POUR VERSION 1.2  LE 05/07/96 PAR M. BENOIT                C
C                                                                     C
C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C
C                                                                     C
C  ARGUMENTS D'APPEL :                                                C
C  *******************                                                C
C  +-------------+----+--------------------------------------------+  C
C  ! NOM         !MODE! SIGNIFICATION - OBSERVATIONS               !  C
C  +-------------+----+--------------------------------------------+  C
C  ! TETAM(-)    !<-- ! VECTEUR DES DIRECTIONS MOYENNES (RADIANS)  !  C
C  ! F(-,-,-)    ! -->! SPECTRE DIRECTIONNEL DE VARIANCE           !  C
C  ! COSTET(-)   ! -->! VECTEUR DES COSINUS DES DIRECTIONS         !  C
C  ! SINTET(-)   ! -->! VECTEUR DES SINUS   DES DIRECTIONS         !  C
C  ! NPLAN       ! -->! NOMBRE DE DIRECTIONS DE DISCRETISATION     !  C
C  ! FREQ(-)     ! -->! VECTEUR DES FREQUENCES DE DISCRETISATION   !  C
C  ! DFREQ(-)    ! -->! VECTEUR DES PAS DE FREQUENCE               !  C
C  ! NF          ! -->! NOMBRE DE FREQUENCES DE DISCRETISATION     !  C
C  ! NPOIN2      ! -->! NOMBRE DE POINTS DU MAILLAGE SPATIAL 2D    !  C
C  ! TAILF       ! -->! FACTEUR DE QUEUE DU SPECTRE                !  C
C  ! COSMOY(-)   !<-->! TABLEAU DE TRAVAIL (DIMENSION NPOIN2)      !  C
C  ! SINMOY(-)   !<-->! TABLEAU DE TRAVAIL (DIMENSION NPOIN2)      !  C
C  ! TAUXC(-)    !<-->! TABLEAU DE TRAVAIL (DIMENSION NPOIN2)      !  C
C  ! TAUXS(-)    !<-->! TABLEAU DE TRAVAIL (DIMENSION NPOIN2)      !  C
C  +-------------+----+--------------------------------------------+  C
C  ! MODE   (-> : NON-MODIFIE)  (<-> : MODIFIE)  (<- : INITIALISE) !  C
C  +---------------------------------------------------------------+  C
C                                                                     C
C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C
C                                                                     C
C  APPELS :    - PROGRAMME(S) APPELANT  :  DUMP2D, TERSOU             C
C  ********    - PROGRAMME(S) APPELE(S) :                             C
C                                                                     C
C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C
C
      IMPLICIT NONE
C
C.....VARIABLES TRANSMISES
C     """"""""""""""""""""
      INTEGER  NF    , NPLAN , NPOIN2
      DOUBLE PRECISION TAILF
      DOUBLE PRECISION TETAM(NPOIN2) , DFREQ(NF)    , COSMOY(NPOIN2)
      DOUBLE PRECISION COSTET(NPLAN) , SINTET(NPLAN), SINMOY(NPOIN2)
      DOUBLE PRECISION TAUXC(NPOIN2) , TAUXS(NPOIN2)
      DOUBLE PRECISION F(NPOIN2,NPLAN,NF), FREQ(NF)
C
C.....VARIABLES LOCALES
C     """""""""""""""""
      INTEGER  IP    , JP    , JF
      DOUBLE PRECISION AUXC  , AUXS  , DEUPI , SEUIL , COEFT , DFDTET
      DOUBLE PRECISION DTETAR
C
C
      DEUPI =2.D0*3.14159265D0
      DTETAR=DEUPI/DBLE(NPLAN)
      SEUIL =1.D-10
      DO 10 IP=1,NPOIN2
        COSMOY(IP)=0.D0
        SINMOY(IP)=0.D0
   10 CONTINUE
C
C-----C-------------------------------------------------------C
C-----C  SOMMATION SUR LA PARTIE DISCRETISEE DU SPECTRE       C
C-----C-------------------------------------------------------C
      DO 30 JF=1,NF
C
        DFDTET=DFREQ(JF)*DTETAR
C
        DO 35 IP=1,NPOIN2
          TAUXC(IP)=0.D0
          TAUXS(IP)=0.D0
   35   CONTINUE
C
        DO 20 JP=1,NPLAN
          AUXC=COSTET(JP)*DFDTET
          AUXS=SINTET(JP)*DFDTET
          DO 40 IP=1,NPOIN2
            TAUXC(IP)=TAUXC(IP)+F(IP,JP,JF)*AUXC
            TAUXS(IP)=TAUXS(IP)+F(IP,JP,JF)*AUXS
   40     CONTINUE
   20   CONTINUE
C
        DO 45 IP=1,NPOIN2
          COSMOY(IP)=COSMOY(IP)+TAUXC(IP)
          SINMOY(IP)=SINMOY(IP)+TAUXS(IP)
   45   CONTINUE
C
   30 CONTINUE
C
C-----C-------------------------------------------------------------C
C-----C  PRISE EN COMPTE EVENTUELLE DE LA PARTIE HAUTES-FREQUENCES  C
C-----C-------------------------------------------------------------C
      IF (TAILF.GT.1.D0) THEN
        COEFT=FREQ(NF)/((TAILF-1.D0)*DFREQ(NF))
        DO 55 IP=1,NPOIN2
          COSMOY(IP)=COSMOY(IP)+TAUXC(IP)*COEFT
          SINMOY(IP)=SINMOY(IP)+TAUXS(IP)*COEFT
   55   CONTINUE
      ENDIF
C
C-----C-------------------------------------------------------------C
C-----C  CALCUL DE LA DIRECTION MOYENNE PROPREMENT DIT              C
C-----C  (EN RADIANS ENTRE 0 ET 2.PI)                               C
C-----C-------------------------------------------------------------C
      DO 60 IP=1,NPOIN2
        IF ((ABS(SINMOY(IP)).LT.SEUIL).AND.
     &      (ABS(COSMOY(IP)).LT.SEUIL)) THEN
          TETAM(IP) = 0.D0
        ELSE
          TETAM(IP)=ATAN2(SINMOY(IP),COSMOY(IP))
          IF (TETAM(IP).LT.0.D0) TETAM(IP)=TETAM(IP)+DEUPI
        ENDIF
   60 CONTINUE
C
      RETURN
      END
