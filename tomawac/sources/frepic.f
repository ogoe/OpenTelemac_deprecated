C                       *****************
                        SUBROUTINE FREPIC
C                       *****************
C
     *( FPIC  , F     , FREQ  , NF    , NPLAN , NPOIN2, EMAX  , E     )
C
C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C
C                                                                     C
C  EDF/DER - LABORATOIRE NATIONAL D'HYDRAULIQUE  -  CHATOU (FRANCE)   C
C                                                                     C
C  CODE TOMAWAC DE MODELISATION DES ETATS DE MER EN ELEMENTS FINIS    C
C                       ---- VERSION 1.2 ----                         C
C                                                                     C
C  FREM02 : CALCUL DE LA FREQUENCE DE PIC DU SPECTRE E(F) EN TOUS LES C
C  ******** POINTS DU MAILLAGE SPATIAL 2D. CETTE FREQUENCE DE PIC     C
C           EST DEFINIE COMME LA FREQUENCE DE DISCRETISATION PORTANT  C
C           LE MAXIMUM D'ENERGIE.                                     C
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
C  ! FPIC(-)     !<-- ! TABLEAU DES FREQUENCES DE PIC              !  C
C  ! F(-,-,-)    ! -->! SPECTRE DIRECTIONNEL DE VARIANCE           !  C
C  ! FREQ(-)     ! -->! TABLEAU DES FREQUENCES DE DISCRETISATION   !  C
C  ! NF          ! -->! NOMBRE DE FREQUENCES DE DISCRETISATION     !  C
C  ! NPLAN       ! -->! NOMBRE DE DIRECTIONS DE DISCRETISATION     !  C
C  ! NPOIN2      ! -->! NOMBRE DE POINTS DU MAILLAGE SPATIAL 2D    !  C
C  ! EMAX(-)     !<-->! TABLEAU DE TRAVAIL (DIMENSION NPOIN2)      !  C
C  ! E(-)        !<-->! TABLEAU DE TRAVAIL (DIMENSION NPOIN2)      !  C
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
      DOUBLE PRECISION F(NPOIN2,NPLAN,NF), FREQ(NF)  , FPIC(NPOIN2)
      DOUBLE PRECISION EMAX(NPOIN2),E(NPOIN2)
C
C.....VARIABLES LOCALES
C     """""""""""""""""
      INTEGER  JP    , JF    , IP
C
C
      DO 10 IP = 1,NPOIN2
        FPIC(IP) = 1.D-20
        EMAX(IP) = 0.D0
   10 CONTINUE
C
C.....ON PARCOURT LES FREQUENCES DE DISCRETISATION.
C     """""""""""""""""""""""""""""""""""""""""""""
      DO 20 JF = 1,NF
C
C.......INTEGRATION SUR LES DIRECTIONS POUR TROUVER E(F).
C       """""""""""""""""""""""""""""""""""""""""""""""""
        DO 60 IP=1,NPOIN2
          E(IP) = 0.D0
   60   CONTINUE
        DO 30 JP = 1,NPLAN
          DO 40 IP=1,NPOIN2
                 E(IP) = E(IP) + F(IP,JP,JF)
   40     CONTINUE
   30   CONTINUE
C
C.......ON GARDE LE MAXIMUM DE E(F) ET LA FREQUENCE ASSOCIEE.
C       """""""""""""""""""""""""""""""""""""""""""""""""""""
        DO 50 IP=1,NPOIN2
          IF (E(IP).GT.EMAX(IP)) THEN
            EMAX(IP) = E(IP)
            FPIC(IP) = FREQ(JF)
          ENDIF
   50   CONTINUE
C
   20 CONTINUE
C
      RETURN
      END
