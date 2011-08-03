!                    ************************
                     SUBROUTINE FRICTION_UNIF
!                    ************************
!
     &(MESH,H,U,V,CHESTR,S,KFROT,KFROTL,ITURB,LISRUG,LINDNER,
     & SB,NDEF,DP,SP,VK,KARMAN,GRAV,T1,T2,CHBORD,CF,CFBOR)
!
!***********************************************************************
! TELEMAC2D   V6P0                                   21/08/2010
!***********************************************************************
!
!brief    COMPUTES FRICTION FOR EACH NODE WHEN THERE IS ONLY
!+                ONE FRICTION LAW IN THE DOMAIN.
!
!history  F. HUVELIN
!+        20/04/2004
!+
!+   WRITTEN FROM COEFRO.F
!
!history  J-M HERVOUET (LNHE)
!+
!+        V5P5
!+
!
!history  N.DURAND (HRW), S.E.BOURBAN (HRW)
!+        13/07/2010
!+        V6P0
!+   Translation of French comments within the FORTRAN sources into
!+   English comments
!
!history  N.DURAND (HRW), S.E.BOURBAN (HRW)
!+        21/08/2010
!+        V6P0
!+   Creation of DOXYGEN tags for automated documentation and
!+   cross-referencing of the FORTRAN sources
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!| CF             |<--| ADIMENSIONAL FRICTION COEFFICIENT
!| CFBORD         |<--| ADIMENSIONAL FRICTION COEFFICIENT ON BOUNDARIES
!| CHBORD         |-->| DEFAULT'S MANNING ON BOUNDARY
!| CHESTR         |-->| FRICTION COEFFICIENTS
!| DP             |-->| DIAMETER OF ROUGHNESS ELEMENT  
!| GRAV           |-->| GRAVITY
!| H              |-->| WATER DEPTH
!| ITURB          |---| NOT USED !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!| KARMAN         |-->| VON KARMAN CONSTANT
!| KFROT          |-->| LAW OF BOTTOM FRICTION
!| LINDNER        |-->| IF YES, THERE IS NON-SUBMERGED VEGETATION FRICTION
!| LISRUG         |-->| TURBULENCE REGIME (1: SMOOTH 2: ROUGH)
!| MESH           |-->| MESH STRUCTURE
!| NDEF           |-->| DEFAULT'S MANNING
!| S              |-->| VOID BIEF_OBJ STRUCTURE
!| SB             |---| NOT USED !!!!!!!!!!!!!!!!!!!!!!
!| SP             |-->| SPACING OF ROUGHNESS ELEMENT 
!| T1             |<->| WORK ARRAY IN A BIEF_OBJ STRUCTURE
!| T2             |<->| WORK ARRAY IN A BIEF_OBJ STRUCTURE
!| U              |-->| X-COMPONENT OF VELOCITY
!| V              |-->| Y-COMPONENT OF VELOCITY
!| VK             |-->| KINEMATIC VISCOSITY
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      USE INTERFACE_TELEMAC2D, EX_FRICTION_UNIF => FRICTION_UNIF
!
      USE BIEF
!
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
      TYPE(BIEF_MESH),  INTENT(IN)      :: MESH
      TYPE(BIEF_OBJ),   INTENT(IN)      :: H,U,V,CHESTR,CHBORD,S
      INTEGER,          INTENT(IN)      :: KFROT,KFROTL,ITURB,LISRUG
      LOGICAL,          INTENT(IN)      :: LINDNER
      DOUBLE PRECISION, INTENT(IN)      :: NDEF,DP,SP
      DOUBLE PRECISION, INTENT(IN)      :: VK,KARMAN,GRAV
!
      DOUBLE PRECISION, INTENT(INOUT)   :: SB
      TYPE(BIEF_OBJ),   INTENT(INOUT)   :: T1, T2
!
      TYPE(BIEF_OBJ),   INTENT(INOUT)   :: CF, CFBOR
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
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
      ! ELEMENT TYPE
      ! ------------
      IELMC = CF%ELM
      IELMH = H%ELM
!
      ! SAME DISCRETIZATION FOR WATER DEPTH AND FRICTION COEFFICIENT IF NEEDED
      ! ----------------------------------------------------------------------
      IF (KFROT.NE.0.AND.KFROT.NE.2) THEN
         !
         ! MAXIMUM BETWEEN WATER DEPTH AND 1.D-4
         ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         CALL CPSTVC(H,T1)
         CALL OS('X=Y     ', T1, H, S, C)
         IF(IELMC.NE.IELMH) CALL CHGDIS( T1 , IELMH , IELMC , MESH )
!        NIKURADSE LAW WILL DO ITS OWN CLIPPING
         IF(KFROT.NE.5) THEN
           CALL OS('X=+(Y,C)', T1, T1, S, 1.D-4)
         ENDIF
      ENDIF
!
      ! RESULTANT VELOCITY IN T2
      ! ------------------------
      IF (KFROT.EQ.1.OR.KFROT.EQ.6.OR.KFROT.EQ.7) THEN
         CALL CPSTVC(CF,T2)
         CALL OS('X=N(Y,Z)', T2,  U, V, C)
         CALL OS('X=+(Y,C)', T2, T2, S, 1.D-6)
      ENDIF
!
      ! =============== !
      ! BOTTOM FRICTION !
      ! =============== !
!
      ! FRICTION COEFFICIENT FOR THE BOTTOM
      ! -----------------------------------
      CALL FRICTION_CALC(1, CF%DIM1, KFROT, NDEF, VK, GRAV,
     &                   KARMAN, CHESTR, T1, T1, T2, CF)
!
      ! FRICTION COEFFICIENT FOR NON-SUBMERGED VEGETATION
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
      ! WALL FRICTION COMPUTATION
      ! -------------------------
!
      IF(LISRUG.EQ.2) THEN
        CALL FRICTION_CALC(1,MESH%NPTFR,KFROTL,NDEF,VK,GRAV,
     &                     KARMAN,CHBORD,MESH%DISBOR,T1,T2,CFBOR)
      ENDIF
!
!=======================================================================!
!=======================================================================!
!
      RETURN
      END
