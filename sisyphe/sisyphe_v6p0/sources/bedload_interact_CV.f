C                      ***************************
                       SUBROUTINE BEDLOAD_INTERACT
C                      ***************************
C
     *(UCMOY,TOBW,TOB,ALPHAW,FW,CF,UW,NPOIN,XMVE,FCW)
C
C**********************************************************************
C SISYPHE VERSION 5.7   01/10/2003             C. VILLARET (LNHE)                 
C**********************************************************************
C
C
C           ! ============================================ !
C           !   Calcul du coefficient de frottement sous   !
C           ! l'action conjugu�e de la houle et du courant !
C           ! ============================================ !
C
C
C COPYRIGHT EDF-BAW-IFH
C----------------------------------------------------------------------!
C                             ARGUMENTS                                !
C .________________.____.______________________________________________!
C |      NOM       |MODE|                   ROLE                       !
C |________________|____|______________________________________________!
C |________________|____|______________________________________________!
C                                                                      !
C                                                                      ! 
C----------------------------------------------------------------------!
C                                                                      !
C CALLED BY BEDLOAD_BAILARD                                            !
C           BEDLOAD_DIBWAT                                             !
C                                                                      !                                                   !
C                                                                      !
C======================================================================!
C======================================================================!
C                    DECLARATION DES TYPES ET DIMENSIONS               !
C======================================================================!
C======================================================================!
C
C 1/ MODULES
C ----------
      USE INTERFACE_SISYPHE,EX_BEDLOAD_INTERACT => BEDLOAD_INTERACT
      USE BIEF
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C 2/ GLOBAL VARIABLES
C -------------------
      TYPE(BIEF_OBJ),   INTENT(IN)  :: UCMOY, TOBW, TOB, ALPHAW
      TYPE(BIEF_OBJ),   INTENT(IN)  :: FW, CF, UW
      INTEGER,          INTENT(IN)  :: NPOIN
      DOUBLE PRECISION, INTENT(IN)  :: XMVE
      TYPE(BIEF_OBJ),   INTENT(INOUT) :: FCW
      ! 3/ LOCAL VARIABLES
      ! ------------------
      INTEGER                     :: I
      DOUBLE PRECISION            :: TX, LOGF
      DOUBLE PRECISION            :: CSAL,CSAL1, CSAL2, CSAL3
      DOUBLE PRECISION            :: AX, MX, NX, BX, PX, QX
      DOUBLE PRECISION            :: UCW2, TAUCW
C LOCAL VARIABLES
      DOUBLE PRECISION             ::ZERO
C
      INTRINSIC MAX
!
!======================================================================!
!======================================================================!
!                               PROGRAMME                              !
!======================================================================!
!======================================================================!
      ZERO = 1.D-6
! 
      DO I = 1, NPOIN
!
            TX = TOB%R(I) / MAX((TOB%R(I) + TOBW%R(I)),ZERO)
   !
         LOGF  = LOG10(2.D0*MAX(FW%R(I),ZERO)/MAX(CF%R(I),ZERO))
         CSAL  = ABS(COS(ALPHAW%R(I)))
         CSAL1 = CSAL**0.82D0
!
! correction le 28/12/2006
!         CSAL2 = CSAL**0.67D0
!
         CSAL3 = CSAL**2.70D0
!
         AX = -0.07D0 + 1.87D0*CSAL1 + (-0.34D0 - 0.12D0*CSAL1)*LOGF
         MX =  0.72D0 - 0.33D0*CSAL1 + ( 0.08D0 + 0.34D0*CSAL1)*LOGF
         NX =  0.78D0 - 0.23D0*CSAL1 + ( 0.12D0 - 0.12D0*CSAL1)*LOGF
!
         BX =  0.27D0 + 0.51D0*CSAL3 + (-0.10D0 - 0.24D0*CSAL3)*LOGF
         PX = -0.75D0 + 0.13D0*CSAL3 + ( 0.12D0 + 0.02D0*CSAL3)*LOGF
         QX =  0.89D0 + 0.40D0*CSAL3 + ( 0.50D0 - 0.28D0*CSAL3)*LOGF
!
         IF(TX <ZERO) THEN
            TAUCW = TOBW%R(I)

         ELSEIF(TX.LT.1.D0.AND.TX.GT.ZERO) THEN
            TAUCW = (1.D0 + BX * TX**PX * (1.D0 - TX)**QX)*TOB%R(I)*TX
     &            + (1.D0 + AX * TX**MX * (1.D0 - TX)**NX)*TOBW%R(I)
         ELSEIF(TX.GT.1.D0) 
            TAUCW = TOB%R(I)
         ENDIF  
!
         UCW2 = (UCMOY%R(I)**2 + 0.5D0 * UW%R(I)**2) * XMVE
         FCW%R(I) = TAUCW / MAX(UCW2,1.D-10)
!    
      ENDDO
!
!======================================================================!
!======================================================================!
!
      RETURN                      
      END
