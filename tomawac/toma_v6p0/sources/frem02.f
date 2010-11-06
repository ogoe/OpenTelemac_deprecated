C                       *****************
                        SUBROUTINE FREM02
C                       *****************
C
     *( FM02  , F     , FREQ  , DFREQ , TAILF , NF    , NPLAN , NPOIN2,
     *  AUX1  , AUX2  )
C
C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C
C                                                                     C
C  EDF/DER - LABORATOIRE NATIONAL D'HYDRAULIQUE  -  CHATOU (FRANCE)   C
C                                                                     C
C  CODE TOMAWAC DE MODELISATION DES ETATS DE MER EN ELEMENTS FINIS    C
C                       ---- VERSION 1.2 ----                         C
C                                                                     C
C  FREM02 : CALCUL DE LA FREQUENCE MOYENNE FM02 DU SPECTRE EN TOUS    C
C  ******** LES POINTS DU MAILLAGE SPATIAL 2D. CETTE FREQUENCE        C
C           DEFINIE A PARTIR DES MOMENTS M0 ET M2 DU SPECTRE :        C
C                                                                     C
C                     SOMME(  F(FREQ,TETA)*FREQ**2 DFREQ DTETA  )     C
C       FM02 = SQRT(  -------------------------------------------  )  C
C                     SOMME(  F(FREQ,TETA) DFREQ DTETA  )             C
C                                                                     C
C   - CREE POUR VERSION 1.0  LE 09/02/95 PAR P. THELLIER ET M. BENOIT C
C   - MOD. POUR VERSION 1.2  LE 05/07/96 PAR M. BENOIT                C
C                                                                     C
C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C
C                                                                     C
C  ARGUMENTS D'APPEL :                                                C
C  *******************                                                C
C  +-------------+----+--------------------------------------------+  C
C  ! NOM         !MODE! SIGNIFICATION - OBSERVATIONS               !  C
C  +-------------+----+--------------------------------------------+  C
C  ! FM02(-)     !<-- ! TABLEAU DES FREQUENCES MOYENNES F02        !  C
C  ! F(-,-,-)    ! -->! SPECTRE DIRECTIONNEL DE VARIANCE           !  C
C  ! FREQ(-)     ! -->! TABLEAU DES FREQUENCES DE DISCRETISATION   !  C
C  ! DFREQ(-)    ! -->! TABLEAU DES PAS DE FREQUENCE               !  C
C  ! TAILF       ! -->! FACTEUR DE QUEUE DU SPECTRE                !  C
C  ! NF          ! -->! NOMBRE DE FREQUENCES DE DISCRETISATION     !  C
C  ! NPLAN       ! -->! NOMBRE DE DIRECTIONS DE DISCRETISATION     !  C
C  ! NPOIN2      ! -->! NOMBRE DE POINTS DU MAILLAGE SPATIAL 2D    !  C
C  ! AUX1(-)     !<-->! TABLEAU DE TRAVAIL (DIMENSION NPOIN2)      !  C
C  ! AUX2(-)     !<-->! TABLEAU DE TRAVAIL (DIMENSION NPOIN2)      !  C
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
C                                                                     C
C  REMARQUES :                                                        C
C  ***********                                                        C
C  - LA PARTIE HAUTES-FREQUENCES DU SPECTRE N'EST PRISE EN COMPTE QUE C
C    SI LE FACTEUR DE QUEUE (TAILF) EST STRICTEMENT SUPERIEUR A 3.    C
C                                                                     C
C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C
C
      IMPLICIT NONE
C
C.....VARIABLES TRANSMISES
C     """"""""""""""""""""
      INTEGER  NF    , NPLAN , NPOIN2
      DOUBLE PRECISION TAILF
      DOUBLE PRECISION F(NPOIN2,NPLAN,NF)
      DOUBLE PRECISION FREQ(NF), DFREQ(NF), FM02(NPOIN2)
      DOUBLE PRECISION AUX1(NPOIN2), AUX2(NPOIN2)
C
C.....VARIABLES LOCALES
C     """""""""""""""""
      INTEGER  JP    , JF    , IP
      DOUBLE PRECISION SEUIL , DTETAR, AUX3  , AUX4
C
C
      SEUIL = 1.D-20
      DTETAR= 2.D0*3.141592654D0/DBLE(NPLAN)
      DO 30 IP = 1,NPOIN2
        AUX1(IP) = 0.D0
        AUX2(IP) = 0.D0
   30 CONTINUE
C
C-----C-------------------------------------------------------C
C-----C  SOMMATION SUR LA PARTIE DISCRETISEE DU SPECTRE       C
C-----C-------------------------------------------------------C
      DO 20 JF = 1,NF-1
        AUX3=DTETAR*DFREQ(JF)
        AUX4=AUX3*FREQ(JF)**2
        DO 10 JP = 1,NPLAN
          DO 5 IP=1,NPOIN2
            AUX1(IP) = AUX1(IP) + F(IP,JP,JF)*AUX4
            AUX2(IP) = AUX2(IP) + F(IP,JP,JF)*AUX3
    5     CONTINUE
   10   CONTINUE
   20 CONTINUE
C
C-----C-------------------------------------------------------------C
C-----C  PRISE EN COMPTE EVENTUELLE DE LA PARTIE HAUTES-FREQUENCES  C
C-----C-------------------------------------------------------------C
      IF (TAILF.GT.3.D0) THEN
        AUX3=DTETAR*(DFREQ(NF)+FREQ(NF)/(TAILF-1.D0))
        AUX4=DTETAR*(DFREQ(NF)+FREQ(NF)/(TAILF-3.D0))*FREQ(NF)**2
      ELSE
        AUX3=DTETAR*DFREQ(NF)
        AUX4=AUX3*FREQ(NF)**2
      ENDIF
      DO 40 JP = 1,NPLAN
        DO 45 IP=1,NPOIN2
          AUX1(IP) = AUX1(IP) + F(IP,JP,NF)*AUX4
          AUX2(IP) = AUX2(IP) + F(IP,JP,NF)*AUX3
   45   CONTINUE
   40 CONTINUE
C
C-----C-------------------------------------------------------------C
C-----C  CALCUL DE LA FREQUENCE MOYENNE PROPREMENT DIT              C
C-----C-------------------------------------------------------------C
      DO 50 IP=1,NPOIN2
        IF (AUX2(IP).LT.SEUIL) THEN
          FM02(IP) = SEUIL
        ELSE
          FM02(IP) = SQRT(AUX1(IP)/AUX2(IP))
        ENDIF
   50 CONTINUE
C
      RETURN
      END
