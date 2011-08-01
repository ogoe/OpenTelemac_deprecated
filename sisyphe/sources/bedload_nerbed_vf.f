      ! **************************** !
        SUBROUTINE BEDLOAD_NERBED_VF !
      ! **************************** !

     &(MESH,LIEBOR,KSORT,ELAY,V2DPAR,QSX,QSY,AVA,NPOIN,NSEG,NPTFR,
     & DT,QS,T1,T2,T3,BREACH)

C**********************************************************************
C SISYPHE VERSION 6.0  14/09/04  F. HUVELIN                            
C SISYPHE VERSION 5.3  07/05/02  M. GONZALES DE LINARES                
C                                                                      
C JMH : BUG D'INITIALISATION DE T1 ET T2 CORRIGE LE 31/01/2008         
C JMH : KSORT ADDED (HARDCODED BEFORE !!!!)                                                                  C
C**********************************************************************


                  ! ====================================== !
                  ! Non erodable method for Finite Volumes !
                  ! ====================================== !


C COPYRIGHT EDF-DTMPL-SOGREAH-LHF-GRADIENT
C**********************************************************************C
C                                                                      C
C                 SSSS I   SSSS Y   Y PPPP  H   H EEEEE                C
C                S     I  S      Y Y  P   P H   H E                    C
C                 SSS  I   SSS    Y   PPPP  HHHHH EEEE                 C
C                    S I      S   Y   P     H   H E                    C
C                SSSS  I  SSSS    Y   P     H   H EEEEE                C
C                                                                      C
C----------------------------------------------------------------------C
C                             ARGUMENTS                                C
C .________________.____.______________________________________________C
C |      NOM       |MODE|                   ROLE                       C
C |________________|____|______________________________________________C
C |________________|____|______________________________________________C
C                    <=  Can't be change by the user                   C
C                    =>  Can be changed by the user                    C 
C ---------------------------------------------------------------------C
!                                                                      !
! CALLED BY BEDLOAD_QNEROD                                             !
!                                                                      !
! CALL      ------                                                     !
!                                                                      !
!======================================================================!
!======================================================================!
!                    DECLARATION DES TYPES ET DIMENSIONS               !
!======================================================================!
!======================================================================!

      ! 1/ MODULES
      ! ----------
      USE INTERFACE_SISYPHE, EX_BEDLOAD_NERBED_VF => BEDLOAD_NERBED_VF
      USE BIEF
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      TYPE(BIEF_MESH),  INTENT(INOUT) :: MESH
      TYPE(BIEF_OBJ),   INTENT(IN)    :: LIEBOR
      TYPE(BIEF_OBJ),   INTENT(IN)    :: QSX, QSY
      INTEGER,          INTENT(IN)    :: NPOIN, NSEG, NPTFR,KSORT
      DOUBLE PRECISION, INTENT(IN)    :: DT
      TYPE(BIEF_OBJ),   INTENT(INOUT) :: QS, T1, T2, T3
      TYPE(BIEF_OBJ),   INTENT(INOUT) :: BREACH
      DOUBLE PRECISION, INTENT(IN)    :: ELAY(NPOIN),V2DPAR(NPOIN)
      DOUBLE PRECISION, INTENT(IN)    :: AVA(NPOIN)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER          :: I, K
      INTEGER          :: IEL, IEL1, IEL2, ISEGIN
      DOUBLE PRECISION :: QSP1, QSP2, QSPC
      DOUBLE PRECISION :: XN, YN, TEMP
      DOUBLE PRECISION :: VNOIN1, VNOIN2, RNORM
!
!======================================================================!
!======================================================================!
!                               PROGRAMME                              !
!======================================================================!
!======================================================================!
!
      ! ****************** !
      ! I - INITIALIZATION ! 
      ! ****************** !
!
      ! BREACH indicates if non erodable bed will be reached
      ! during time step for this point
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!     DONNER A T1 ET T2 LA MEME STRUCTURE QUE QS
!      
      CALL CPSTVC(QS,T1)
      CALL CPSTVC(QS,T2)
!      
      DO IEL = 1, NPOIN
        BREACH%I(IEL) = 0
        T1%R(IEL)=0.D0
        T2%R(IEL)=0.D0
      ENDDO
!
      ! ************************************************* !
      ! II - DETERMINER LE FLUX SORTANT DE CHAQUE CELLULE ! (_IMP_)
      ! ************************************************* !
      ! Le principe est que pour chaque segment on calcule QS comme 
      ! demi-somme des QS aux points qui sont les centres des elements 
      ! dont le segment est la frontiere, on le projette sur la normale
      ! du segment, on le multiplie par la longueur du segment, et on 
      ! ajoute (ou retranche) ce terme de flux aux deux elements.
!
      DO ISEGIN = 1, NSEG
!
         IEL1 = MESH%NUBO%I(2*ISEGIN - 1)
         IEL2 = MESH%NUBO%I(2*ISEGIN    )
!
         ! II.1 - Longueur du segment (RNORM) 
         ! ----------------------------------
         VNOIN1 = MESH%VNOIN%R(3*ISEGIN - 2)
         VNOIN2 = MESH%VNOIN%R(3*ISEGIN - 1)
         RNORM  = MESH%VNOIN%R(3*ISEGIN    )
!
         ! II.2 - Projection de QS au segment sur la normale du segment 
         ! ------------------------------------------------------------
         QSP1 = VNOIN1*QSX%R(IEL1) + VNOIN2*QSY%R(IEL1)
         QSP2 = VNOIN1*QSX%R(IEL2) + VNOIN2*QSY%R(IEL2)
         QSPC = (QSP1+QSP2)*0.5D0
!
         ! II.3 - QS tel que le flux sortant soit maximal 
         ! ----------------------------------------------
         T1%R(IEL1) = T1%R(IEL1) + RNORM*MAX(QSPC,QSP1,0.D0)
         T1%R(IEL2) = T1%R(IEL2) - RNORM*MIN(QSPC,QSP2,0.D0)
!
         IF(QSPC > 0.D0) THEN
           T2%R(IEL1) = T2%R(IEL1) + RNORM*QSP1
         ELSEIF(QSPC < 0.D0) THEN
           T2%R(IEL2) = T2%R(IEL2) - RNORM*QSP2
         ENDIF
!
      ENDDO
!
      ! ************************************** !
      ! III - BOUCLE SUR LES POINTS FRONTIERES ! 
      ! ************************************** !
!      
      DO K = 1, NPTFR
         IEL = MESH%NBOR%I(K)
!
         ! III.1 - Evolution libre : on laisse sortir les sediments
         ! ---------------------------------------------------------
         IF (LIEBOR%I(K) == KSORT) THEN
!
            ! XNEBOR(*+NPTFR) et YNEBOR(*+NPTFR)
            ! contiennent le vecteur normal a un point frontiere
            ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            XN   = MESH%XNEBOR%R(K+NPTFR)
            YN   = MESH%YNEBOR%R(K+NPTFR)
            TEMP = QSX%R(IEL)*XN + QSY%R(IEL)*YN

            IF (TEMP > 0.D0) THEN
               T1%R(IEL) = T1%R(IEL) + TEMP
               T2%R(IEL) = T2%R(IEL) + TEMP
            ENDIF
!
         ENDIF
!
         ! III.2 - Pour la paroi solide il n'y a rien a programmer 
         !         car le flux de sediments y est nul               
         !         Pour l'evolution imposee : cf bedload_solvs_vf.f 
         ! --------------------------------------------------------
      ENDDO
!
      IF(NCSIZE > 1) THEN
        CALL PARCOM(T1, 2, MESH)
        CALL PARCOM(T2, 2, MESH)
      ENDIF
!
      ! ************************************************ !
      ! IV - CALCUL DU FLUX MAXIMUM AUTORISE PAR CELLULE ! 
      ! ************************************************ !
!      
      DO I = 1, NPOIN
!
         T3%R(I)=ELAY(I)*V2DPAR(I)*AVA(I)*(1.D0-1.D-6)/DT
         IF (T3%R(I) < 0.D0) T3%R(I) = 0.D0
!
         ! Si flux sortant trop fort, on limite le QS au noeud
         ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         IF(T1%R(I) > T3%R(I)) THEN
            BREACH%I(I) = 1
            IF(T2%R(I) > T3%R(I)) THEN
              QS%R(I) = QS%R(I)*T3%R(I)/T2%R(I)
            ENDIF
         ENDIF
!
      ENDDO
!
!======================================================================!
!======================================================================!
!
      RETURN
      END
