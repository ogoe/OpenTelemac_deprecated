      ! ************************* !
        SUBROUTINE BEDLOAD_CALCDW ! (_IMP_)
      ! ************************* !

     &  (UCW, UW, TW, NPOIN, PI, UW1, UW2, TW1, TW2)

C**********************************************************************C
C SISYPHE VERSION 5.4  Oct   2003  C. VILLARET                         C
C**********************************************************************C


               ! ================================ !
               ! Calcul des vitesse qudratique et !
               !   des periodes en cas de houle   !
               ! ================================ !


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
! CALLED BY BEDLOAD_DIBWAT                                             !
!                                                                      !
! CALL      ------                                                     !
!                                                                      !
!======================================================================!
!======================================================================!
!                    DECLARATION DES TYPES ET DIMENSIONS               !
!======================================================================!
!======================================================================!
!
!     1/ MODULES
! 
      USE INTERFACE_SISYPHE,EX_BEDLOAD_CALCDW => BEDLOAD_CALCDW
      USE BIEF
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
!
!     2/ GLOBAL VARIABLES
!
      TYPE(BIEF_OBJ),   INTENT(IN)    :: UCW, UW, TW
      INTEGER,          INTENT(IN)    :: NPOIN
      DOUBLE PRECISION, INTENT(IN)    :: PI
      TYPE(BIEF_OBJ),   INTENT(INOUT) :: UW1, UW2, TW1, TW2
!
!     3/ LOCAL VARIABLES
! 
      INTEGER                     :: I
      DOUBLE PRECISION            :: UCMOY, RAP
      DOUBLE PRECISION            :: ACOSMRAP, ACOSPRAP, SQRTRAP
      DOUBLE PRECISION, PARAMETER :: ZERO = 1.D-06
!
!======================================================================!
!======================================================================!
!                               PROGRAMME                              !
!======================================================================!
!======================================================================!
!
      DO I = 1,NPOIN
         UCMOY = ABS(UCW%R(I))
         ! ****************** !
         !    I - HOULE SEULE ! (_IMP_)
         ! ****************** !
         IF (UCMOY <= ZERO) THEN
            UW1%R(I) = UW%R(I)
            UW2%R(I) = UW%R(I) 
            TW1%R(I) = TW%R(I) / 2.D0
            TW2%R(I) = TW%R(I) / 2.D0
         ELSE
            RAP = UW%R(I) / UCMOY
            ! ******************** !
            ! II - HOULE DOMINANTE ! (_IMP_)
            ! ******************** !
            IF (RAP > 1.D0) THEN
               ACOSMRAP = ACOS(-1.D0/RAP)
               ACOSPRAP = ACOS( 1.D0/RAP)
               SQRTRAP  = SQRT(1.D0-1.D0/RAP**2)
               TW1%R(I) = TW%R(I)*ACOSMRAP / PI
               TW2%R(I) = TW%R(I)*ACOSPRAP / PI
               UW1%R(I) = 2.D0*UCMOY**2 + UW%R(I)**2
     &                  + 3.D0*UCMOY*UW%R(I)*SQRTRAP/ACOSMRAP
               UW1%R(I) = SQRT(UW1%R(I))
               UW2%R(I) = 2.D0*UCMOY**2 + UW%R(I)**2
     &                  - 3.D0*UCMOY*UW%R(I)*SQRTRAP/ACOSPRAP
               UW2%R(I) = SQRT(UW2%R(I))

            ! ********************** !
            ! III - COURANT DOMINANT ! (_IMP_)
            ! ********************** !
            ELSE
               UW1%R(I) = UCW%R(I)*SQRT(2.D0 + RAP**2)
               UW2%R(I) = ZERO
               TW1%R(I) = TW%R(I)
               TW2%R(I) = ZERO
            ENDIF
         ENDIF
      ENDDO
!
!======================================================================!
!======================================================================!
!
      RETURN
      END SUBROUTINE BEDLOAD_CALCDW
