      ! **************************** !
        SUBROUTINE BEDLOAD_DIRECTION ! (_IMP_)
      ! **************************** !

     &  (QU, QV, NPOIN, PI, THETAC)


C**********************************************************************C
C SISYPHE VERSION 5.4  01/10/2003  C. VILLARET                         C
C**********************************************************************C


           ! ============================================= !
           ! Calcul de l'angle THETAC (direction du debit) !
           ! ============================================= !


C COPYRIGHT EDF-BAW-IFH
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
! CALLED BY BEDLOAD_BAILARD                                            !
!           BEDLOAD_DIBWAT                                             !
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
      USE INTERFACE_SISYPHE,
     &    EX_BEDLOAD_DIRECTION => BEDLOAD_DIRECTION
      USE BIEF
      IMPLICIT NONE               
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU                                


      ! 2/ GLOBAL VARIABLES
      ! -------------------
      TYPE(BIEF_OBJ),  INTENT(IN)  :: QU, QV
      INTEGER,          INTENT(IN)  :: NPOIN
      DOUBLE PRECISION, INTENT(IN)  :: PI
      TYPE(BIEF_OBJ),  INTENT(INOUT) :: THETAC


      ! 3/ LOCAL VARIABLES
      ! ------------------
      INTEGER                     :: I
      DOUBLE PRECISION, PARAMETER :: LOCAL_ZERO = 1.D-6


!======================================================================!
!======================================================================!
!                    DECLARATION DES TYPES ET DIMENSIONS               !
!======================================================================!
!======================================================================!

      DO I = 1, NPOIN

         IF (ABS(QU%R(I)) <= LOCAL_ZERO) THEN

            IF (QV%R(I) < = LOCAL_ZERO) THEN
               THETAC%R(I) = -PI*0.5D0
            ELSE
               THETAC%R(I) =  PI*0.5D0
            ENDIF

         ELSE

            THETAC%R(I) = ATAN(QV%R(I) / QU%R(I))

            IF (QU%R(I) < 0.D0) THEN
               THETAC%R(I) = PI + THETAC%R(I)
            ENDIF

         ENDIF

      END DO

!======================================================================!
!======================================================================!

      RETURN                      
      END SUBROUTINE BEDLOAD_DIRECTION
