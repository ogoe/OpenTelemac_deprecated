C                       *****************
                        SUBROUTINE FSPRD3
C                       *****************
C
     *( FRA   , DIREC , NPLAN , SPRED1, TETA1 , SPRED2, TETA2 , XLAMDA,
     *  DEUPI )
C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C
C                                                                     C
C  EDF/DER - LABORATOIRE NATIONAL D'HYDRAULIQUE  -  CHATOU (FRANCE)   C
C                                                                     C
C  CODE TOMAWAC DE MODELISATION DES ETATS DE MER EN ELEMENTS FINIS    C
C                       ---- VERSION 1.2 ----                         C
C                                                                     C
C  FSPRD3 : CALCUL DE LA FONCTION DE REPARTITION ANGULAIRE BIMODALE   C
C  ******** POUR UNE SERIE DE DIRECTIONS.                             C
C                                                                     C
C                        2S                                           C
C           MODELE EN COS  ((T-T0)/2.) (DE TYPE MITSUYASU)            C
C                                                                     C
C   - CREE POUR VERSION 1.2  LE 07/11/96 PAR M. BENOIT                C
C                                                                     C
C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C
C                                                                     C
C  ARGUMENTS D'APPEL :                                                C
C  *******************                                                C
C  +-------------+----+--------------------------------------------+  C
C  ! NOM         !MODE! SIGNIFICATION - OBSERVATIONS               !  C
C  +-------------+----+--------------------------------------------+  C
C  ! FRA(-)      !<-- ! VALEURS DE LA FONCTION DE REPARTITION ANG. !  C
C  ! DIREC(-)    ! -->! DIRECTIONS DE DISCRETISATION               !  C
C  ! NPLAN       ! -->! NOMBRE DE DIRECTIONS DE DISCRETISATION     !  C
C  ! SPRED1      ! -->! ETALEMENT DIRECTIONNEL 1 POUR FRA          !  C
C  ! TETA1       ! -->! DIRECTION PRINCIPALE 1 POUR FRA            !  C
C  ! SPRED2      ! -->! ETALEMENT DIRECTIONNEL 2 POUR FRA          !  C
C  ! TETA2       ! -->! DIRECTION PRINCIPALE 2 POUR FRA            !  C
C  ! XLAMDA      ! -->! FACTEUR DE PONDERATION POUR LA FRA         !  C
C  ! DEUPI       ! -->! 2.PI                                       !  C
C  +-------------+----+--------------------------------------------+  C
C  ! MODE   (-> : NON-MODIFIE)  (<-> : MODIFIE)  (<- : INITIALISE) !  C
C  +---------------------------------------------------------------+  C
C                                                                     C
C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C
C                                                                     C
C  APPELS :    - PROGRAMME(S) APPELANT  :  SPEINI, LIMWAC, ...        C
C  ********    - PROGRAMME(S) APPELE(S) :  DELFRA                     C
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
      INTEGER  NPLAN
      DOUBLE PRECISION SPRED1, TETA1 , SPRED2, TETA2 , XLAMDA, DEUPI
      DOUBLE PRECISION FRA(NPLAN)    , DIREC(NPLAN)
C
C.....VARIABLES LOCALES 
C     """""""""""""""""
      INTEGER  JP
      DOUBLE PRECISION DELT1 , DELT2 , FTH   , FRA1  , FRA2  , ARGUM
      DOUBLE PRECISION ARGMI1, ARGMI2
C
C.....FONCTION EXTERNES 
C     """""""""""""""""
      DOUBLE PRECISION DELFRA
      EXTERNAL         DELFRA
C
C
      DELT1 = 0.5D0/DELFRA(SPRED1,DEUPI)
      DELT2 = 0.5D0/DELFRA(SPRED2,DEUPI)
      IF (SPRED1.GT.1.D-1) THEN
        ARGMI1=10.D0**(-4.D0/SPRED1)
      ELSE
        ARGMI1=0.D0
      ENDIF
      IF (SPRED2.GT.1.D-1) THEN
        ARGMI2=10.D0**(-4.D0/SPRED2)
      ELSE
        ARGMI2=0.D0
      ENDIF
C
      DO 110 JP=1,NPLAN
        FTH = DIREC(JP)
C
        ARGUM = ABS(COS(0.5D0*(FTH-TETA1)))
        IF (ARGUM.GT.ARGMI1) THEN
          FRA1=DELT1*ARGUM**(2.D0*SPRED1)
        ELSE
          FRA1=0.D0
        ENDIF
C
        ARGUM = ABS(COS(0.5D0*(FTH-TETA2)))
        IF (ARGUM.GT.ARGMI2) THEN
          FRA2=DELT2*ARGUM**(2.D0*SPRED2)
        ELSE
          FRA2=0.D0
        ENDIF
C
        FRA(JP)=XLAMDA*FRA1+(1.D0-XLAMDA)*FRA2
        IF (FRA(JP).LT.1.D-10) FRA(JP)=0.D0
  110 CONTINUE
C
      RETURN
      END
