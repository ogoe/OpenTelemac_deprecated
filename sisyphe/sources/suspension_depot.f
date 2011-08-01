cT1: TOB
c À FAIRE  REMPLACER ustar PAR tob
c 
      ! ***************************  !
        SUBROUTINE SUSPENSION_DEPOT  !
      ! ***************************  !

     &(TOB,HN, ACLADM,NPOIN, HMIN,XWC,VITCD,
     & ZERO,KARMAN,XMVE, T1,T2,ZREF,FLUDPT,DEBUG,SEDCO)

C**********************************************************************
C SISYPHE VERSION 5.9  31/07/08               J-M HERVOUET + C VILLARET            
C
C**********************************************************************


          ! ================================================= !
          ! Computation of the flux of deposition and erosion !
          ! ================================================= !


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
C |    HN          | => |  HAUTEUR D'EAU
C |________________|____|______________________________________________C
C                                   C
C ---------------------------------------------------------------------C
! called by suspension_computation                                          !
! call suspension_rouse                                                !   
!                                                                      !
!======================================================================!
!======================================================================!
!                    DECLARATION DES TYPES ET DIMENSIONS               !
!======================================================================!
!======================================================================!

      ! 1/ MODULES
      ! ----------
      USE INTERFACE_SISYPHE,EX_SUSPENSION_DEPOT => SUSPENSION_DEPOT
      USE BIEF
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU         
C
      ! 2/ GLOBAL VARIABLES
      ! -------------------
      TYPE (BIEF_OBJ),  INTENT(IN)    ::  HN, ACLADM,TOB
      INTEGER,          INTENT(IN)    ::  NPOIN,DEBUG
      LOGICAL,          INTENT(IN)    :: SEDCO
      DOUBLE PRECISION, INTENT(IN)    ::  HMIN
      DOUBLE PRECISION, INTENT(IN)    :: XWC
      DOUBLE PRECISION, INTENT(IN)    :: VITCD
      DOUBLE PRECISION, INTENT(IN)    :: ZERO, KARMAN,XMVE
      TYPE (BIEF_OBJ),  INTENT(INOUT) :: T1,T2
      TYPE (BIEF_OBJ),  INTENT(IN)    :: ZREF
      TYPE (BIEF_OBJ),  INTENT(INOUT) :: FLUDPT
C
      ! 3/ LOCAL VARIABLES
      ! ------------------
      INTEGER :: I
      DOUBLE PRECISION:: AUX
!======================================================================!
!                               PROGRAMME                              !
!======================================================================!
!======================================================================!   
!     ! ****************************************            !
      ! THE  TOTAL FRICTION VELOCIT    --> USTAR (T1)              !
      ! remplace par  (V6p0) USTARP (skin friction velocity)
      ! for erosion flux                           
      ! ****************************************            !


      CALL OS('X=CY    ', X=T1, Y=TOB, C=1.D0/XMVE)
      CALL OS('X=+(Y,C)', X=T1, Y=T1, C=ZERO)
      CALL OS('X=SQR(Y)', X=T1, Y=T1) 


      IF(SEDCO) THEN
!      
      ! ************************************************ !
      ! Ia - FORMULATION FOR COHESIVE SEDIMENTS COHESIFS ! 
      !      (WITHOUT BEDLOAD)                           ! 
      ! ************************************************ !      
!
!  COMPUTATION OF THE PROBABILITY FOR DEPOSITION   
!
         DO I = 1, NPOIN
           IF(VITCD.GT.1.D-08) THEN
             AUX = MAX(1.D0-(T1%R(I)/VITCD)**2,ZERO)
           ELSE
             AUX=0.D0
           ENDIF   
!          COMPUTATION OF THE IMPLICIT PART OF DEPOSITION FLUX
           FLUDPT%R(I)= XWC*AUX
         ENDDO
! sediment uniforme sur la verticale
         CALL OS('X=C     ', X=T2, C=1.D0)
!         
      ! ******************************************* !
      ! Ib - FORMULATION FOR NON-COHESIVE SEDIMENTS ! 
      !      (WITH BEDLOAD)                         ! 
      ! ******************************************* !
!      
      ELSE
!       
            ! ***************************************************** !
            !  COMPUTATION OF THE RATIO BETWEEN NEAR BED CONC AND MEAN CONC  ! 
            !                                  -->  T2    (à conserver )    !
            ! ***************************************************** !
        IF (DEBUG > 0) WRITE(LU,*) 'SUSPENSION_ROUSE'
        CALL SUSPENSION_ROUSE(T1,HN,NPOIN, 
     &                        KARMAN,HMIN,ZERO,XWC,ZREF,T2)
        IF (DEBUG > 0) WRITE(LU,*) 'END SUSPENSION_ROUSE'
!
            ! *****************************************************  !
            !  COMPUTATION OF DEPOSITION FLUX --> FLUDPT = XWC * T2  !
            ! *****************************************************  !
!        
         CALL OS('X=CY    ', X=FLUDPT, Y=T2, C=XWC)
!
      ENDIF
!
!======================================================================!
!======================================================================!
!
      RETURN      
      END
