!                    **************************
                     SUBROUTINE FRICTION_CHOICE
!                    **************************
!
     &(FRICTION_PASS,KARMAN)
!
!***********************************************************************
! TELEMAC2D   V6P1                                   21/08/2010
!***********************************************************************
!
!brief    MAIN SUBROUTINE FOR FRICTION COMPUTATION.
!
!history  F. HUVELIN
!+        20/04/2004
!+
!+
!
!history  J-M HERVOUET (LNHE)
!+
!+        V6P0
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
!| FRICTION_PASS  |-->| IF 0, INITIALISATION
!| KARMAN         |-->| VON KARMAN CONSTANT
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      USE BIEF
      USE FRICTION_DEF
      USE DECLARATIONS_TELEMAC
      USE DECLARATIONS_TELEMAC2D
      USE INTERFACE_TELEMAC2D, EX_FRICTION_CHOICE => FRICTION_CHOICE
!
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
      INTEGER,          INTENT(IN) :: FRICTION_PASS
      DOUBLE PRECISION, INTENT(IN) :: KARMAN
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
      INTEGER                     :: I
      DOUBLE PRECISION, PARAMETER :: VK = 1.D-6
!
!======================================================================!
!======================================================================!
!                               PROGRAMME                              !
!======================================================================!
!======================================================================!
!
! INITIALIZATION
! --------------
!
      IF(FRICTION_PASS == 0) THEN
!
! ZONES INITIALIZATION
! --------------------
!
         IF (FRICTB) THEN
            CALL FRICTION_INIT
            CALL STRCHE
! FH : FOR QUASI-BUBBLE
! FH : 2004/03/01
!JAJ FOR QUADRATIC ELEMENTS
! =>
            IF (CF%ELM /= H%ELM) THEN
              IF(CF%ELM==12 .AND. H%ELM==11) THEN
                CALL FRICTION_BUBBLE
     &              (IKLE, NPOIN, NELEM, NELMAX, LINDNER, NKFROT,
     &               CHESTR, NDEFMA, LINDDP, LINDSP)
              ELSEIF (CF%ELM==13 .AND. H%ELM==11) THEN
                CALL FRICTION_QUAD
     &              (IKLE%I, NPOIN, NELEM, NELMAX, LINDNER, NKFROT,
     &               CHESTR, NDEFMA, LINDDP, LINDSP)
!                WRITE (LU,*)
!     &           'FRICTION_CHOICE::QUADRATIC ELEMENTS NOT IMPLEMENTED.'
!                CALL PLANTE(1)
              ELSE
                WRITE (LU,*)
     &           'FRICTION_CHOICE::DISCRETISATION NOT IMPLEMENTED.'
                WRITE (LU,*)
     &           'CF%ELM, H%ELM: ',CF%ELM, H%ELM
                CALL PLANTE(1)
              ENDIF
            ENDIF
! <=
!JAJ FOR QUADRATIC ELEMENTS
! FH : 2004/03/01
! FH : FOR QUASI-BUBBLE
!
!        UNIFORM CASE
!        ------------
!
         ELSE
            ! CHESTR FOR BOUNDARY CONDITIONS INITIALIZATION
            ! -----------------------------------------------
            IF(LISRUG.EQ.2) THEN
              IF(KFROTL.EQ.1) THEN
                DO I = 1, NPTFR
                  CHBORD%R(I) = CHESTR%R(MESH%NBOR%I(I))
                ENDDO
              ELSE
!               JMH 21/12/2010
!               BOUNDARY CONDITIONS FILE DATA IF ANY SUPERSEDE
!               THE KEY-WORD ROUGHNESS COEFFICIENT OF BOUNDARIES
                IF(P_DOTS(CHBORD,CHBORD,MESH).EQ.0.D0) THEN
                  CALL OS('X=C     ', X=CHBORD, C=SB)
                ENDIF
              ENDIF
            ENDIF
            ! TYPE OF FRICTION LAW FOR EACH NODE
            ! ----------------------------------
            DO I=1, CF%DIM1
              NKFROT%I(I) = KFROT
            ENDDO
!
         ENDIF
!
!     COMPUTATION
!     -----------
!
      ELSE
!
         ! FRICTION BY ZONES
         ! -----------------
         IF (FRICTB) THEN
! FH : FOR QUASI-BUBBLE
! FH : 2004/03/01
!JAJ FOR QUADRATIC ELEMENTS
! =>
            IF (CF%ELM /= H%ELM) THEN
              IF (CF%ELM==12 .AND. H%ELM==11) THEN
                CALL FRICTION_BUBBLE
     &              (IKLE, NPOIN, NELEM, NELMAX, LINDNER, NKFROT,
     &               CHESTR, NDEFMA, LINDDP, LINDSP)
              ELSE IF (CF%ELM==13 .AND. H%ELM==11) THEN
                CALL FRICTION_QUAD
     &              (IKLE%I, NPOIN, NELEM, NELMAX, LINDNER, NKFROT,
     &               CHESTR, NDEFMA, LINDDP, LINDSP)
              ELSE
                WRITE (LU,*)
     &           'FRICTION_CHOICE::DISCRETISATION NOT IMPLEMENTED.'
                WRITE (LU,*)
     &           'CF%ELM, H%ELM: ',CF%ELM, H%ELM
                CALL PLANTE(1)
                STOP
              ENDIF
            ENDIF
!
            CALL FRICTION_ZONES
     &           (MESH, H, U, V, S, CHESTR, CHBORD, NKFROT, NDEFMA,
     &            LINDDP, LINDSP, KFRO_B, NDEF_B, ITURB, LISRUG,
     &            LINDNER, VK, KARMAN, GRAV, T1, T2, CF, CFBOR)
! <=
!JAJ FOR QUADRATIC ELEMENTS
! FH : 2004/03/01
! FH : FOR QUASI-BUBBLE
!
         ! UNIFORM FRICTION
         ! ----------------
         ELSE
            CALL FRICTION_UNIF
     &           (MESH,H,U,V,CHESTR,S,KFROT,KFROTL,ITURB, LISRUG,
     &            LINDNER, SB, NDEF, DP, SP, VK, KARMAN, GRAV, T1,
     &            T2, CHBORD, CF, CFBOR)
!
         ENDIF
      ENDIF
!
!======================================================================!
!======================================================================!
!
      RETURN
      END
