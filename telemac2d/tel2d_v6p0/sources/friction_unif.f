C                       ************************
                        SUBROUTINE FRICTION_UNIF
C                       ************************
C
     & (MESH, H, U, V, CHESTR, S, KFROT, ITURB, LISRUG, LINDNER,
     &  SB, NDEF, DP, SP, VK, KARMAN, GRAV, T1, T2, CHBORD, CF, CFBOR)
C
C***********************************************************************
C  TELEMAC-2D VERSION 5.5                 J-M HERVOUET (LNH) 30 87 80 18
C***********************************************************************
C
C 20/04/04 : subroutine written by F. Huvelin (from coefro.f)
C
C
          ! ------------------------------------------------ !
          !         Friction calculation for each node       !
          ! when there is only one frictionlaw in the domain !
          ! ------------------------------------------------ !
C
C
C               TTTTT EEEEE L     EEEEE M   M   AA  CCCCC
C                 T   E     L     E     MM MM  A  A C
C                 T   EEE   L     EEE   M M M  AAAA C
C                 T   E     L     E     M   M  A  A C
C                 T   EEEEE LLLLL EEEEE M   M  A  A CCCCC
C
C
C----------------------------------------------------------------------C
C                             ARGUMENTS                                C
C .________________.____.______________________________________________C
C |      NOM       |MODE|                   ROLE                       C
C |________________|____|______________________________________________C
C |                | => |                                              C
C |________________|____|______________________________________________C
C                    <=  input value                                   C
C                    =>  output value                                  C 
C ---------------------------------------------------------------------C
!
!=======================================================================!
!=======================================================================!
!                    DECLARATION DES TYPES ET DIMENSIONS                !
!=======================================================================!
!=======================================================================!
!
      USE INTERFACE_TELEMAC2D, EX_FRICTION_UNIF => FRICTION_UNIF
!
      USE BIEF
!
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      TYPE(BIEF_MESH),  INTENT(IN)      :: MESH
      TYPE(BIEF_OBJ),   INTENT(IN)      :: H, U, V, CHESTR, CHBORD, S
      INTEGER,          INTENT(IN)      :: KFROT, ITURB, LISRUG
      LOGICAL,          INTENT(IN)      :: LINDNER
      DOUBLE PRECISION, INTENT(IN)      :: NDEF, DP, SP
      DOUBLE PRECISION, INTENT(IN)      :: VK, KARMAN, GRAV
!
      DOUBLE PRECISION, INTENT(INOUT)   :: SB
      TYPE(BIEF_OBJ),   INTENT(INOUT)   :: T1, T2
!
      TYPE(BIEF_OBJ),   INTENT(INOUT)   :: CF, CFBOR
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER          :: IELMC, IELMH, I
      DOUBLE PRECISION :: C, CP
!
!=======================================================================!
!=======================================================================!
!                               PROGRAMME                               !
!=======================================================================!
!=======================================================================!
!
! ======================================= !
! INITIALIZATION AND DISCRETIZATION CHECK !
! ======================================= !
!
      ! Element type
      ! ------------
      IELMC = CF%ELM
      IELMH = H%ELM
!
      ! Same discretization for water depth and friction coefficient if needed
      ! ----------------------------------------------------------------------
      IF (KFROT.NE.0.AND.KFROT.NE.2) THEN
!
         ! Maximum between Water depth and 1.D-4
         ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         CALL CPSTVC(H,T1)
         CALL OS('X=Y     ', T1, H, S, C)
         IF(IELMC.NE.IELMH) CALL CHGDIS( T1 , IELMH , IELMC , MESH )
         CALL OS('X=+(Y,C)', T1, T1, S, 1.D-4)
      ENDIF
!
      ! Resultant velocity in T2
      ! ------------------------
      IF ((KFROT==1).OR.(KFROT==6).OR.(KFROT==7)) THEN
         CALL CPSTVC(CF,T2)
         CALL OS('X=N(Y,Z)', T2,  U, V, C)
         CALL OS('X=+(Y,C)', T2, T2, S, 1.D-6)
      ENDIF
!
      ! =============== !
      ! BOTTOM FRICTION !
      ! =============== !
!
      ! Friction coefficient for the bottom
      ! -----------------------------------
      CALL FRICTION_CALC
     &     (1, CF%DIM1, KFROT, NDEF, VK, GRAV,
     &      KARMAN, CHESTR, T1, T1, T2, CF)
!
      ! Friction coefficient for non-submerged vegetation
      ! -------------------------------------------------
      IF(LINDNER) THEN
!
         DO I = 1, CF%DIM1
            CALL FRICTION_LINDNER
     &           (T2%R(I), T1%R(I), CF%R(I),
     &            VK, GRAV, DP, SP, CP)
!
            IF (CP< -0.9) THEN
               CP = 0.75*T1%R(I)*DP/(SP**2)
            ENDIF
!         
            CF%R(I) = (CF%R(I)+2.D0*CP)
!
         ENDDO
      ENDIF
!
! ============= !
! WALL FRICTION !
! ============= !
!
      ! Wall friction computation
      ! -------------------------
      IF ((ITURB == 3).AND.(LISRUG == 2)) THEN
!
         CALL FRICTION_CALC
     &        (1, MESH%NPTFR, KFROT, NDEF, VK, GRAV,
     &         KARMAN, CHBORD, MESH%DISBOR, T1, T2, CFBOR)
      ENDIF
!
!=======================================================================!
!=======================================================================!
!
      RETURN
      END SUBROUTINE FRICTION_UNIF
