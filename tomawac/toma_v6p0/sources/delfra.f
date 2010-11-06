C                       ***************
                        FUNCTION DELFRA
C                       ***************
C
     *( SS    , DEUPI )
C
C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C
C                                                                     C
C  EDF/DER - LABORATOIRE NATIONAL D'HYDRAULIQUE  -  CHATOU (FRANCE)   C
C                                                                     C
C  CODE TOMAWAC DE MODELISATION DES ETATS DE MER EN ELEMENTS FINIS    C
C                       ---- VERSION 1.2 ----                         C
C                                                                     C
C  DELFRA : CALCUL DU COEFFICIENT DE NORMALISATION DE LA FONCTION     C
C  ******** DE REPARTITION ANGULAIRE EN COS **2.S (TETA-TETA0).       C
C                                                                     C
C                               GAMMA( SS + 0.5)                      C
C        DELFRA(SS) = SQRT(PI)  ----------------                      C
C                               GAMMA( SS + 1. )                      C
C                                                                     C
C   - CREE POUR VERSION 1.0  LE 15/11/95 PAR M. BENOIT                C
C   - MOD. POUR VERSION 1.2  LE 07/11/96 PAR M. BENOIT                C
C                                                                     C
C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C
C                                                                     C
C  ARGUMENTS D'APPEL :                                                C
C  *******************                                                C
C  +-------------+----+--------------------------------------------+  C
C  ! NOM         !MODE! SIGNIFICATION - OBSERVATIONS               !  C
C  +-------------+----+--------------------------------------------+  C
C  ! DELFRA      !<-- ! COEFFICIENT DE NORMALISATION DE LA FRA     !  C
C  ! SS          ! -->! EXPOSANT DE LA FONCTION DE REPARTITION ANG.!  C
C  ! DEUPI       ! -->! 2.PI                                       !  C
C  +-------------+----+--------------------------------------------+  C
C  ! MODE   (-> : NON-MODIFIE)  (<-> : MODIFIE)  (<- : INITIALISE) !  C
C  +---------------------------------------------------------------+  C
C                                                                     C
C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C
C                                                                     C
C  APPELS :    - PROGRAMME(S) APPELANT  :  FSPRD1, FSPRD3             C
C  ********    - PROGRAMME(S) APPELE(S) :  GAMMLN                     C
C                                                                     C
C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C
C                                                                     C
C  REMARQUES :                                                        C
C  ***********                                                        C
C                                                                     C
C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C
C
      IMPLICIT NONE
C
C.....VARIABLES TRANSMISES
C     """"""""""""""""""""
      DOUBLE PRECISION DELFRA, SS    , DEUPI
C
C.....FONCTIONS EXTERNES
C     """"""""""""""""""
      DOUBLE PRECISION GAMMLN
      EXTERNAL         GAMMLN
C
C
      DELFRA=SQRT(DEUPI/2.D0)
     *      *EXP(GAMMLN(SS+0.5D0,DEUPI)-GAMMLN(SS+1.D0,DEUPI))
C
      RETURN
      END
