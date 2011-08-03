!                    *****************
                     SUBROUTINE SUMVER
!                    *****************
!
     &(FINT,F,NPOIN2,NPLAN,NPTFR)
!
!***********************************************************************
! TELEMAC3D   V6P1                                   21/08/2010
!***********************************************************************
!
!brief    COMPUTES THE SUM ON THE VERTICAL OF A 3D VARIABLE
!+                DEFINED ON LATERAL BOUNDARIES.
!+
!+            THE RESULT IS A 2D BOUNDARY VECTOR.
!
!history  J-M HERVOUET (LNH)
!+        11/02/08
!+        V5P9
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
!| F              |-->| VARIABLES TO AVERAGE
!| FINT           |<->| MEAN VARIABLE
!| NPLAN          |-->| NUMBER OF HORIZONTAL PLANES
!| NPOIN2         |-->| NUMBER OF POINTS PER HORIZONTAL PLANE
!| NPTFR          |-->| NUMBER OF BOUNDARY POINTS
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
      INTEGER, INTENT(IN) :: NPOIN2,NPLAN,NPTFR
      DOUBLE PRECISION, INTENT(INOUT) :: FINT(NPTFR)
      DOUBLE PRECISION, INTENT(IN)    :: F(NPTFR,NPLAN)
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
      INTEGER IP,IPTFR
!
!-----------------------------------------------------------------------
!
      DO IPTFR=1,NPTFR
        FINT(IPTFR) = F(IPTFR,1)
      ENDDO
!
      DO IP = 2,NPLAN
        DO IPTFR=1,NPTFR
          FINT(IPTFR)=FINT(IPTFR)+F(IPTFR,IP)
        ENDDO
      ENDDO
!
!-----------------------------------------------------------------------
!
      RETURN
      END
