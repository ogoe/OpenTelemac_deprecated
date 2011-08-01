       SUBROUTINE BEDLOAD_MEYER ! (_IMP_) 
      ! ************************ ! 
 
     &  (TETAP, HIDING, HIDFAC, DENS, GRAV, DM, AC, 
     &   ACP, QSC, SLOPEFF, COEFPN) 
 
 
C**********************************************************************C 
C SISYPHE VERSION 5.4  --/10/2003   C.VILLARET                         C 
C SISYPHE VERSION 5.1  11/09/1995  E. PELTIER                          C 
C SISYPHE VERSION 5.1  11/09/1995  C. LENORMANT                        C 
C SISYPHE VERSION 5.1  11/09/1995  J.-M. HERVOUET                      C 
C**********************************************************************C 
 
 
           ! =========================================== ! 
           !   Bed-load transport formula of Meyer-Peter ! 
           ! =========================================== ! 
 
 
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
     %    EX_BEDLOAD_MEYER => BEDLOAD_MEYER 
      USE BIEF 
      IMPLICIT NONE 
      INTEGER LNG,LU 
      COMMON/INFO/LNG,LU 
 
 
      ! 2/ GLOBAL VARIABLES 
      ! ------------------- 
      TYPE(BIEF_OBJ),   INTENT(IN)    :: TETAP, HIDING 
      INTEGER,          INTENT(IN)    :: HIDFAC, SLOPEFF 
      DOUBLE PRECISION, INTENT(IN)    :: DENS, GRAV, DM, AC 
      TYPE(BIEF_OBJ),   INTENT(INOUT) :: ACP ! work array T1
      TYPE(BIEF_OBJ),   INTENT(INOUT)   :: QSC, COEFPN 
 
 
      ! 3/ LOCAL VARIABLES 
      ! ------------------ 
      DOUBLE PRECISION :: C2 
 
 
!======================================================================! 
!======================================================================! 
!                               PROGRAMME                              ! 
!======================================================================! 
!======================================================================! 
 
      CALL OS('X=C     ', X=ACP, C=AC) 
 
      ! **************************************** ! 
      ! 0 - EFFET DE PENTE : FORMULE DE SOULBY   ! (_IMP_) 
      ! **************************************** ! 
      IF(SLOPEFF == 2) THEN 
        CALL OS('X=XY    ', X=ACP, Y=COEFPN ) 
      ENDIF 

 
      ! **************************************** ! 
      ! III - TRANSPORT PAR CHARRIAGE AVEC       ! (_IMP_) 
      !       CORRECTION POUR LA GRANULO ETENDUE ! (_IMP_) 
      ! **************************************** ! 
      C2 = 8.D0 * SQRT(GRAV*DENS*DM**3) 
      IF ((HIDFAC == 1) .OR. (HIDFAC == 2) ) THEN 
         CALL OS('X=XY    ', X=ACP, Y=HIDING) 
         CALL OS('X=Y-Z   ', X=QSC, Y=TETAP, Z=ACP) 
         CALL OS('X=+(Y,C)', X=QSC, Y=QSC , C=0.D0) 
         CALL OS('X=Y**C  ', X=QSC, Y=QSC , C=1.5D0) 
         CALL OS('X=CX    ', X=QSC, C=C2) 
      ELSE 
          CALL OS('X=Y-Z   ', X=QSC, Y=TETAP, Z=ACP) 
          CALL OS('X=+(Y,C)', X=QSC, Y=QSC, C=0.D0) 
         CALL OS('X=Y**C  ', X=QSC, Y=QSC, C=1.5D0) 
         CALL OS('X=CX    ', X=QSC, C=C2) 
         CALL OS('X=XY    ', X=QSC, Y=HIDING) 
      ENDIF 
 
!======================================================================! 
!======================================================================! 
 
      RETURN 
      END SUBROUTINE BEDLOAD_MEYER 
