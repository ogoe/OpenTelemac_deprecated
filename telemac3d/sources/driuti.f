!                    *****************
                     SUBROUTINE DRIUTI
!                    *****************
!
     & (FRI, RI, ITYP, ITRAC, NPOIN3)
!
!***********************************************************************
! TELEMAC3D   V6P1                                   21/08/2010
!***********************************************************************
!
!brief    USER DEFINED DAMPING FUNCTIONS.
!
!warning  USER SUBROUTINE; DEFAULT VALUE FOR THE DUMPING
!+            FUNCTION (FRI) IS 1
!warning  THE KEYWORD 'DAMPING FUNCTION' MUST BE SET TO 1 IN
!+            THE STEERING FILE (USER-DEFINED)
!
!history  A MALCHEREK (HANOVRE); E PELTIER (LNH)    ; F LEPEINTRE (LNH)
!+        25/11/97
!+        V5P1
!+
!
!history  JACEK A. JANKOWSKI PINXIT
!+        **/03/99
!+
!+   FORTRAN95 VERSION
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
!| FRI            |<->| DAMPING FUNCTION
!| ITRAC          |-->| TRACER NUMBER
!| ITYP           |-->| =1 FOR VELOCITIES
!|                |   | =2 FOR TRACERS
!| NPOIN3         |-->| NUMBER OF POINTS IN THE 3D MESH
!| RI             |<->| RICHARDSON NUMBER
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
      INTEGER, INTENT(IN) :: ITYP, ITRAC, NPOIN3
      DOUBLE PRECISION, INTENT(INOUT) :: FRI(NPOIN3)
      DOUBLE PRECISION, INTENT(INOUT) :: RI(NPOIN3)
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
      IF(ITYP.EQ.1) THEN
!
!        DAMPING FUNCTION FOR VELOCITIES
!
         CALL OV('X=C     ',FRI,FRI,FRI,1.D0,NPOIN3)
!
      ELSEIF(ITYP.EQ.2) THEN
!
!        DAMPING FUNCTION FOR TRACERS
!
         CALL OV('X=C     ',FRI,FRI,FRI,1.D0,NPOIN3)
!
      ELSE
!
         IF (LNG.EQ.1) WRITE(LU,11) ITYP
         IF (LNG.EQ.2) WRITE(LU,12) ITYP
         CALL PLANTE(1)
         STOP
!
      ENDIF
!
!-----------------------------------------------------------------------
!
11    FORMAT('DRIUTI: VARIABLE NON PREVUE ITYP: ',I2)
12    FORMAT('DRIUTI: UNEXPECTED PARAMETER ITYP: ',I2)
!
!-----------------------------------------------------------------------
!
      RETURN
      END SUBROUTINE DRIUTI
