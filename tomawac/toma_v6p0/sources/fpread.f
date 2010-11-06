C                       *****************
                        SUBROUTINE FPREAD
C                       *****************
C
     *( FREAD , F     , FREQ  , DFREQ , NF    , NPLAN , NPOIN2, EXPO  ,
     *  TAILF , DENOM , E     )
C
C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C
C                                                                     C
C  EDF/DER - LABORATOIRE NATIONAL D'HYDRAULIQUE  -  CHATOU (FRANCE)   C
C                                                                     C
C  CODE TOMAWAC DE MODELISATION DES ETATS DE MER EN ELEMENTS FINIS    C
C                       ---- VERSION 1.2 ----                         C
C                                                                     C
C  FPREAD :  CALCUL DE LA FREQUENCE DE PIC DU SPECTRE DE VARIANCE     C
C  ********  PAR LA METHODE DITE DE READ.                             C
C                                                                     C
C   - CREE POUR VERSION 1.1  LE 30/01/96 PAR M. BENOIT                C
C   - MOD. POUR VERSION 1.2  LE 05/07/96 PAR M. BENOIT                C
C                                                                     C
C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C
C                                                                     C
C  ARGUMENTS D'APPEL :                                                C
C  *******************                                                C
C  +-------------+----+--------------------------------------------+  C
C  ! NOM         !MODE! SIGNIFICATION - OBSERVATIONS               !  C
C  +-------------+----+--------------------------------------------+  C
C  ! FREAD(-)    !<-- ! TABLEAU DES FREQUENCES DE PIC (METH. READ) !  C
C  ! F(-,-,-)    ! -->! SPECTRE DIRECTIONNEL DE VARIANCE           !  C
C  ! FREQ(-)     ! -->! TABLEAU DES FREQUENCES DE DISCRETISATION   !  C
C  ! DFREQ(-)    ! -->! TABLEAU DES PAS DE FREQUENCES              !  C
C  ! NF          ! -->! NOMBRE DE FREQUENCES DE DISCRETISATION     !  C
C  ! NPLAN       ! -->! NOMBRE DE DIRECTIONS DE DISCRETISATION     !  C
C  ! NPOIN2      ! -->! NOMBRE DE POINTS DU MAILLAGE SPATIAL 2D    !  C
C  ! EXPO        ! -->! EXPOSANT DE LA METHODE DE READ             !  C
C  ! TAILF       ! -->! FACTEUR DE QUEUE DU SPECTRE                !  C
C  ! DENOM(-)    !<-->! TABLEAU DE TRAVAIL (DIMENSION NPOIN2)      !  C
C  ! E(-)        !<-->! TABLEAU DE TRAVAIL (DIMENSION NPOIN2)      !  C
C  +-------------+----+--------------------------------------------+  C
C  ! MODE   (-> : NON-MODIFIE)  (<-> : MODIFIE)  (<- : INITIALISE) !  C
C  +---------------------------------------------------------------+  C
C                                                                     C
C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C
C                                                                     C
C  APPELS :    - PROGRAMME(S) APPELANT  :  TERSOU, DUMP2D             C
C  ********    - PROGRAMME(S) APPELE(S) :                             C
C                                                                     C
C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C
C
      IMPLICIT NONE
C
C.....VARIABLES TRANSMISES
C     """"""""""""""""""""
      INTEGER  NF    , NPLAN , NPOIN2
      DOUBLE PRECISION EXPO  , TAILF
      DOUBLE PRECISION F(NPOIN2,NPLAN,NF), FREQ(NF), DFREQ(NF)
      DOUBLE PRECISION DENOM(NPOIN2), E(NPOIN2), FREAD(NPOIN2)
C
C.....VARIABLES LOCALES
C     """""""""""""""""
      INTEGER  JP    , JF    , IP
      DOUBLE PRECISION SEUIL , AUXI  , COEFN  , COEFD , DTETAR
C
C
      SEUIL =1.D-20
      DTETAR=2.D0*3.141592654D0/DBLE(NPLAN)
      DO 10 IP = 1,NPOIN2
        FREAD(IP)=0.D0
        DENOM(IP)=0.D0
   10 CONTINUE
C
C-----C-------------------------------------------------------C
C-----C  SOMMATION SUR LA PARTIE DISCRETISEE DU SPECTRE       C
C-----C-------------------------------------------------------C
      DO 20 JF=1,NF
C
C.......INTEGRATION SUR LES DIRECTIONS POUR TROUVER E(F).
C       """""""""""""""""""""""""""""""""""""""""""""""""
        DO 60 IP=1,NPOIN2
          E(IP) = 0.D0
   60   CONTINUE
        DO 30 JP=1,NPLAN
          DO 40 IP=1,NPOIN2
                 E(IP) = E(IP) + F(IP,JP,JF)*DTETAR
   40     CONTINUE
   30   CONTINUE
C
C.......ON SOMME LA CONTRIBUTION DE LA FREQUENCE F.
C       """""""""""""""""""""""""""""""""""""""""""
        DO 50 IP=1,NPOIN2
          IF (E(IP).GT.SEUIL) THEN
            AUXI = E(IP)**EXPO*DFREQ(JF)
            FREAD(IP) = FREAD(IP)+AUXI*FREQ(JF)
            DENOM(IP) = DENOM(IP)+AUXI
          ENDIF
   50   CONTINUE
C
   20 CONTINUE
C
C-----C-------------------------------------------------------------C
C-----C  PRISE EN COMPTE EVENTUELLE DE LA PARTIE HAUTES-FREQUENCES  C
C-----C-------------------------------------------------------------C
      IF (TAILF.GT.1.D0) THEN
        COEFN=FREQ(NF)**2/(TAILF*EXPO-2.D0)
        COEFD=FREQ(NF)   /(TAILF*EXPO-1.D0)
        DO 55 IP=1,NPOIN2
          AUXI=E(IP)**EXPO
          FREAD(IP) = FREAD(IP)+AUXI*COEFN
          DENOM(IP) = DENOM(IP)+AUXI*COEFD
   55   CONTINUE
      ENDIF
C
C-----C-------------------------------------------------------------C
C-----C  CALCUL DE LA FREQUENCE DE PIC PROPREMENT DIT               C
C-----C-------------------------------------------------------------C
      DO 70 IP=1,NPOIN2
        IF (DENOM(IP).LT.1.D-90) THEN
          FREAD(IP) = SEUIL
        ELSE
          FREAD(IP) = FREAD(IP)/DENOM(IP)
        ENDIF 
   70 CONTINUE
C
      RETURN
      END
