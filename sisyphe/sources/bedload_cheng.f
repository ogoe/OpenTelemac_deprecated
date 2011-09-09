      ! ************************ !
       SUBROUTINE BEDLOAD_CHENG  
      ! ************************ ! 
 
     &  (TETAP, NPOIN, DENS, GRAV, DM, DSTAR, QSC) 
 
 
C**********************************************************************C 
C                                                                      C
C SISYPHE VERSION 6.1  11/01/2011  O.Goethel                           C 
C**********************************************************************C 
 
 
           ! =========================================== ! 
           !   Bed-load transport formula of Cheng       ! 
           ! =========================================== ! 
 
 
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
! CALLED BY BEDLOAD_FORMULA                                            ! 
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
     %    EX_BEDLOAD_CHENG => BEDLOAD_CHENG 
      USE BIEF 
      IMPLICIT NONE 
      INTEGER LNG,LU 
      COMMON/INFO/LNG,LU 
 
 
      ! 2/ GLOBAL VARIABLES 
      ! ------------------- 
      TYPE(BIEF_OBJ),   INTENT(IN)    :: TETAP
      INTEGER,          INTENT(IN)    :: NPOIN 
      DOUBLE PRECISION, INTENT(IN)    :: DENS, GRAV, DM, DSTAR
      TYPE(BIEF_OBJ),   INTENT(INOUT)   :: QSC
 
 
      ! 3/ LOCAL VARIABLES 
      ! ------------------ 
      !DOUBLE PRECISION :: C2
      INTEGER :: I
 
 
!======================================================================! 
!======================================================================! 
!                               PROGRAMME                              ! 
!======================================================================! 
!======================================================================! 

      ! **************************** !
      ! II - TRANSPORT PAR CHARRIAGE ! (_IMP_)
      ! **************************** !

      DO I = 1, NPOIN

         IF (TETAP%R(I) > 0.D0) THEN
            QSC%R(I) = 13.D0*TETAP%R(I)**1.5D0*exp(-0.05D0/
     &                 (TETAP%R(I)**1.5D0)) * SQRT(DENS*GRAV*(DM**3))
         ELSE
            QSC%R(I) = 0.D0
         ENDIF

      ENDDO

!======================================================================! 
!======================================================================! 
 
      RETURN 
      END SUBROUTINE BEDLOAD_CHENG 
