C                       *************************
                        SUBROUTINE FRICTION_ZONES
C                       *************************
C
     & (MESH, H, U, V, S, CHESTR, CHBORD, NKFROT, NDEFMA, LINDDP,
     &  LINDSP, KFRO_B, NDEF_B, ITURB, LISRUG, LINDNER, VK,
     &  KARMAN, GRAV, T1, T2, CF, CFBOR)
C
C***********************************************************************
C  TELEMAC-2D VERSION 5.5                 J-M HERVOUET (LNH) 30 87 80 18
C***********************************************************************
C
C 20/04/04 : subroutine written by F. Huvelin
C
C
          ! ----------------------------------------------- !
          !   Friction calculation for each node and zone   !
          ! ----------------------------------------------- !
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
      USE BIEF
      USE INTERFACE_TELEMAC2D, EX_FRICTION_ZONES => FRICTION_ZONES
C
      IMPLICIT NONE      
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      TYPE(BIEF_MESH),    INTENT(IN)    :: MESH
      TYPE(BIEF_OBJ),     INTENT(IN)    :: H, U, V, S
      TYPE(BIEF_OBJ),     INTENT(IN)    :: CHESTR
      TYPE(BIEF_OBJ),     INTENT(IN)    :: CHBORD
      TYPE(BIEF_OBJ),     INTENT(IN)    :: NKFROT
      TYPE(BIEF_OBJ),     INTENT(IN)    :: NDEFMA, LINDDP, LINDSP
      TYPE(BIEF_OBJ),     INTENT(IN)    :: KFRO_B, NDEF_B 
      INTEGER,            INTENT(IN)    :: ITURB, LISRUG
      LOGICAL,            INTENT(IN)    :: LINDNER 
      DOUBLE PRECISION,   INTENT(IN)    :: VK, KARMAN, GRAV
      TYPE(BIEF_OBJ),     INTENT(INOUT) :: CF, CFBOR
      TYPE(BIEF_OBJ),     INTENT(INOUT) :: T1, T2
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER          :: I, J
      INTEGER          :: IELMC, IELMH
      DOUBLE PRECISION :: CP
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
      ! Maximum between Water depth and 1.D-4
      ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      CALL CPSTVC(H,T1)
      CALL OS('X=Y     ', X=T1, Y=H)
      IF(IELMC.NE.IELMH) CALL CHGDIS(T1, IELMH, IELMC, MESH)
      CALL OS('X=+(Y,C)', X=T1, Y=T1, C=1.D-4)
!
      ! Resultant velocity in T2
      ! ------------------------
      CALL CPSTVC(CF,T2)
      CALL OS('X=N(Y,Z)', X=T2, Y=U, Z=V)
      CALL OS('X=+(Y,C)', X=T2, Y=T2, C=1.D-6)
!
!
      ! =============== !
      ! BOTTOM FRICTION !
      ! =============== !
!
      ! Bottom friction calculation
      ! ---------------------------
      DO I = 1, CF%DIM1
!
         ! Friction coefficient for the bottom
         ! -----------------------------------
         CALL FRICTION_CALC
     &        (I, I, NKFROT%I(I), NDEFMA%R(I), VK, GRAV,
     &         KARMAN, CHESTR, T1, T1, T2, CF)
!
         ! Friction coefficient for non-submerged vegetation
         ! -------------------------------------------------
         IF (LINDNER) THEN
            CALL FRICTION_LINDNER
     &           (T2%R(I), T1%R(I), CF%R(I), VK, GRAV, 
     &            LINDDP%R(I),LINDSP%R(I), CP)
!
            IF (CP < -0.9D0) THEN
               CP = 0.75D0*T1%R(I)*LINDDP%R(I) / (LINDSP%R(I))**2
            ENDIF
!         
            CF%R(I) = (CF%R(I)+2.D0*CP)
         ENDIF
!
      END DO
!
!
      ! ============= !
      ! WALL FRICTION !
      ! ============= !
      IF (ITURB.EQ.3.AND.LISRUG.EQ.2) THEN
!
         DO J = 1, MESH%NPTFR
            I = MESH%NBOR%I(J)
!
            ! Bottom friction calculation
            ! ---------------------------
            CALL FRICTION_CALC
     &           (J, J, KFRO_B%I(J), NDEF_B%R(J), VK, GRAV,KARMAN,
     &            CHBORD, MESH%DISBOR , T1, T2, CFBOR)
         END DO
!
      ENDIF
!
!=======================================================================!
!=======================================================================!
      RETURN
      END SUBROUTINE FRICTION_ZONES
