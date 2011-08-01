C                     ************************  
                      SUBROUTINE INIT_CONSTANT
C                     ************************  
C
     &(KARIM_HOLLY_YANG,KARMAN,PI)
C
C**********************************************************************
C SISYPHE VERSION 5.7  11/01/04                           F. HUVELIN                            
C**********************************************************************
C
C
C            ! ===================================== !
C            ! Setting the constants used by sisyphe !
C            ! ===================================== !
C
C
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
! CALLED BY SISYPHE                                                    !
!                                                                      !
! CALL                                                                 !
!                                                                      !
!======================================================================!
!======================================================================!
!                    DECLARATION DES TYPES ET DIMENSIONS               !
!======================================================================!
!======================================================================!
!
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      DOUBLE PRECISION, INTENT(INOUT) :: KARIM_HOLLY_YANG
      DOUBLE PRECISION, INTENT(INOUT) :: KARMAN
      DOUBLE PRECISION, INTENT(INOUT) :: PI
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
!
! Karim, Holly & Yang constant
! ----------------------------
!
      KARIM_HOLLY_YANG = 0.85D0
!
! Von Karman constant
! -------------------
!
      KARMAN = 0.4D0
!
! Partheniades constant : doit être exprimée en m/s
! ---------------------
! cette valeur est maintenant donnée dans user_krone_part
!      PARTHENIADES = 2.D-5/XMVS
!
! PI
! --
!
      PI =ACOS(-1.D0)
!
!======================================================================!
!
      RETURN
      END SUBROUTINE
