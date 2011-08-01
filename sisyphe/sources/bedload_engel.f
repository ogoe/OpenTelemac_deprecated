      ! ************************ !
      ! ************************ !
        SUBROUTINE BEDLOAD_ENGEL ! (_IMP_)
      ! ************************ !

     &  (TOB, CF, DENS, GRAV, DM, XMVE, TETA, QSC)


C**********************************************************************C
C SISYPHE VERSION 5.4  --/10/2003   C.VILLARET                         C
C SISYPHE VERSION 5.1  11/09/1995  E. PELTIER                          C
C SISYPHE VERSION 5.1  11/09/1995  C. LENORMANT                        C
C SISYPHE VERSION 5.1  11/09/1995  J.-M. HERVOUET                      C
C**********************************************************************C


            ! ============================================= !
            ! Bed-load transport formula of Engelund-Hansen !
            ! ============================================= !


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
     %    EX_BEDLOAD_ENGEL => BEDLOAD_ENGEL
      USE BIEF
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU


      ! 2/ GLOBAL VARIABLES
      ! -------------------
      TYPE(BIEF_OBJ),   INTENT(IN)    :: TOB, CF
      DOUBLE PRECISION, INTENT(IN)    :: DENS, GRAV, DM, XMVE
      TYPE(BIEF_OBJ),   INTENT(INOUT) :: TETA ! work array T1
      TYPE(BIEF_OBJ),   INTENT(INOUT) :: QSC


      ! 3/ LOCAL VARIABLES
      ! ------------------
      DOUBLE PRECISION :: CENGEL, C1


!======================================================================!
!======================================================================!
!                               PROGRAMME                              !
!======================================================================!
!======================================================================!

      ! ************************** !
      ! I - CONTRAINTE TOTALE ADIM ! (_IMP_)
      ! ************************** !
      C1 = 1.D0/(DENS*XMVE*GRAV*DM)
      CALL OS('X=CY    ', X=TETA, Y=TOB , C=C1)
      CALL OS('X=Y**C  ', X=TETA, Y=TETA, C=5.D0/2.D0)


      ! *************************** !
      ! II - TRANSPORT PAR CHARRIAGE ! (_IMP_)
      ! *************************** !
      CENGEL = 0.1D0*SQRT(DENS*GRAV*DM**3)
      CALL OS('X=+(Y,C)', X=QSC , Y=CF  , C=1.D-06)
      CALL OS('X=1/Y   ', X=QSC , Y=QSC)
      CALL OS('X=CXY   ', X=QSC , Y=TETA, C=CENGEL)

!======================================================================!
!======================================================================!

      RETURN
      END SUBROUTINE BEDLOAD_ENGEL
