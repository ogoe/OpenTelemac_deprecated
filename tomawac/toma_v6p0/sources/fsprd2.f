C                       *****************
                        SUBROUTINE FSPRD2
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
C  FSPRD2 : CALCUL DE LA FONCTION DE REPARTITION ANGULAIRE BIMODALE   C
C  ******** POUR UNE SERIE DE DIRECTIONS.                             C
C                                                                     C
C           MODELE EN EXP -0.5((T-T0)/S)**2  T DANS (T0-PI/2;T0+PI/2) C
C                                                                     C
C   - CREE POUR VERSION 1.0  LE 10/01/96 PAR M. BENOIT                C
C   - MOD. POUR VERSION 1.2  LE 07/11/96 PAR M. BENOIT                C
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
C  ********    - PROGRAMME(S) APPELE(S) :                             C
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
      DOUBLE PRECISION PI    , C1    , C2
C
C
      PI=DEUPI/2.D0
      IF (SPRED1.GT.1.D-4) THEN
        DELT1 = 1.D0/(SPRED1*SQRT(DEUPI))
        C1    = -0.5/(SPRED1*SPRED1)
      ELSE
        DELT1 = 0.D0
        C1    = 0.D0
      ENDIF
      IF (SPRED2.GT.1.D-4) THEN
        DELT2 = 1.D0/(SPRED2*SQRT(DEUPI))
        C2    = -0.5/(SPRED2*SPRED2)
      ELSE
        DELT2 = 0.D0
        C2    = 0.D0
      ENDIF
C
      DO 110 JP=1,NPLAN
        FTH = DIREC(JP)
C
        ARGUM = FTH-TETA1
   50   CONTINUE
        IF (ARGUM.LT.-PI) THEN
          ARGUM=ARGUM+DEUPI
          GOTO 50
        ENDIF
   60   CONTINUE
        IF (ARGUM.GT.PI) THEN
          ARGUM=ARGUM-DEUPI
          GOTO 60
        ENDIF
        FRA1=DELT1*EXP(MAX(-10.D0,C1*ARGUM*ARGUM))
C
        ARGUM = FTH-TETA2
   70   CONTINUE
        IF (ARGUM.LT.-PI) THEN
          ARGUM=ARGUM+DEUPI
          GOTO 70
        ENDIF
   80   CONTINUE
        IF (ARGUM.GT.PI) THEN
          ARGUM=ARGUM-DEUPI
          GOTO 80
        ENDIF
        FRA2=DELT2*EXP(MAX(-10.D0,C2*ARGUM*ARGUM))
C
        FRA(JP)=XLAMDA*FRA1+(1.D0-XLAMDA)*FRA2
        IF (FRA(JP).LT.1.D-10) FRA(JP)=0.D0
  110 CONTINUE
C
      RETURN
      END
