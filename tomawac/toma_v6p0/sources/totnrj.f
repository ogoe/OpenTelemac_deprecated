C                       *****************
                        SUBROUTINE TOTNRJ
C                       *****************
C
     *( VARIAN, F     , FREQ  , DFREQ , TAILF , NF    , NPLAN , NPOIN2)
C
C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C
C                                                                     C
C  EDF/DER - LABORATOIRE NATIONAL D'HYDRAULIQUE  -  CHATOU (FRANCE)   C
C                                                                     C
C  CODE TOMAWAC DE MODELISATION DES ETATS DE MER EN ELEMENTS FINIS    C
C                       ---- VERSION 1.2 ----                         C
C                                                                     C
C  TOTNRJ : CALCUL DE LA VARIANCE DU SPECTRE DIRECTIONNEL EN TOUS     C
C  ******** LES POINTS DU MAILLAGE SPATIAL 2D. CETTE VARIANCE EST     C
C           OBTENUE PAR INTEGRATION SUR LES FREQUENCES ET DIRECTIONS  C
C           ET PREND EVENTUELLEMENT EN COMPTE LA PARTIE HAUTES        C
C           FREQUENCES DU SPECTRE.                                    C
C                                                                     C
C   - CREE POUR VERSION 1.0  LE 09/02/95 PAR P. THELLIER ET M. BENOIT C
C                                                                     C
C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C
C                                                                     C
C  ARGUMENTS D'APPEL :                                                C
C  *******************                                                C
C  +-------------+----+--------------------------------------------+  C
C  ! NOM         !MODE! SIGNIFICATION - OBSERVATIONS               !  C
C  +-------------+----+--------------------------------------------+  C
C  ! VARIAN(-)   !<-- ! TABLEAU DES VARIANCES                      !  C
C  ! F(-,-,-)    ! -->! SPECTRE DIRECTIONNEL DE VARIANCE           !  C
C  ! FREQ(-)     ! -->! TABLEAU DES FREQUENCES DE DISCRETISATION   !  C
C  ! DFREQ(-)    ! -->! TABLEAU DES PAS DE FREQUENCE               !  C
C  ! TAILF       ! -->! FACTEUR DE QUEUE DU SPECTRE                !  C
C  ! NF          ! -->! NOMBRE DE FREQUENCES DE DISCRETISATION     !  C
C  ! NPLAN       ! -->! NOMBRE DE DIRECTIONS DE DISCRETISATION     !  C
C  ! NPOIN2      ! -->! NOMBRE DE POINTS DU MAILLAGE SPATIAL 2D    !  C
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
C    SI LE FACTEUR DE QUEUE (TAILF) EST STRICTEMENT SUPERIEUR A 1.    C
C                                                                     C
C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C
C
      IMPLICIT NONE
C
C.....VARIABLES TRANSMISES
C     """"""""""""""""""""
      INTEGER          NF    , NPLAN , NPOIN2
      DOUBLE PRECISION TAILF , VARIAN(NPOIN2), FREQ(NF), DFREQ(NF)
      DOUBLE PRECISION F(NPOIN2,NPLAN,NF)
C
C.....VARIABLES LOCALES
C     """""""""""""""""
      INTEGER          IP    , JP    , JF
      DOUBLE PRECISION AUX1  , DTETAR
C
C
      DTETAR=2.D0*3.14159265D0/FLOAT(NPLAN)
      DO 30 IP = 1,NPOIN2
        VARIAN(IP) = 0.D0
30    CONTINUE
C
C-----C-------------------------------------------------------C
C-----C  SOMMATION SUR LA PARTIE DISCRETISEE DU SPECTRE       C
C-----C-------------------------------------------------------C
      DO 20 JF = 1,NF-1
        AUX1=DFREQ(JF)*DTETAR
        DO 10 JP = 1,NPLAN
          DO 5 IP=1,NPOIN2
            VARIAN(IP) = VARIAN(IP) + F(IP,JP,JF)*AUX1
    5     CONTINUE
   10   CONTINUE
   20 CONTINUE
C
C-----C-------------------------------------------------------------C
C-----C  PRISE EN COMPTE EVENTUELLE DE LA PARTIE HAUTES-FREQUENCES  C
C-----C-------------------------------------------------------------C
      IF (TAILF.GT.1.D0) THEN
        AUX1=DTETAR*(DFREQ(NF)+FREQ(NF)/(TAILF-1.D0))
      ELSE
        AUX1=DTETAR*DFREQ(NF)
      ENDIF
      DO 40 JP = 1,NPLAN
        DO 45 IP=1,NPOIN2
          VARIAN(IP) = VARIAN(IP) + F(IP,JP,NF)*AUX1
   45   CONTINUE
   40 CONTINUE
C
      RETURN
      END
