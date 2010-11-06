      ! ************************ !
        SUBROUTINE BEDLOAD_EINST ! (_IMP)
      ! ************************ !

     &  (TETAP, NPOIN, DENS, GRAV, DM, DSTAR, QSC)


C**********************************************************************C
C SISYPHE VERSION 5.4  --/10/2003   C.VILLARET                         C
C SISYPHE VERSION 5.1  11/09/1995  E. PELTIER                          C
C SISYPHE VERSION 5.1  11/09/1995  C. LENORMANT                        C
C SISYPHE VERSION 5.1  11/09/1995  J.-M. HERVOUET                      C
C**********************************************************************C


           ! ============================================== !
           !   Bed-load transport formula of Einstein-Brown !
           ! ============================================== !


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
     %    EX_BEDLOAD_EINST => BEDLOAD_EINST
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
      INTEGER          :: I
      DOUBLE PRECISION :: CEINST


!======================================================================!
!======================================================================!
!                               PROGRAMME                              !
!======================================================================!
!======================================================================!


 
      ! **************************** !
      ! II - TRANSPORT PAR CHARRIAGE ! (_IMP_)
      ! **************************** !
      CEINST = 36.D0/(DSTAR**3)
      CEINST = SQRT(2.D0/3.D0+CEINST) -  SQRT(CEINST)
      CEINST = CEINST * SQRT(DENS*GRAV*(DM**3))
      DO I = 1, NPOIN

         IF (TETAP%R(I) < 2.5D-3) THEN
            QSC%R(I) = 0.D0
         ELSE IF (TETAP%R(I) < 0.2D0) THEN
            QSC%R(I) = 2.15D0* CEINST * EXP(-0.391D0/TETAP%R(I))
         ELSE
            QSC%R(I) = 40.D0 * CEINST * (TETAP%R(I)**3.D0)
         ENDIF

      ENDDO

!======================================================================!
!======================================================================!

      RETURN
      END SUBROUTINE BEDLOAD_EINST
