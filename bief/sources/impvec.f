!                    *****************
                     SUBROUTINE IMPVEC
!                    *****************
!
     &(VEC,NOM,NPOIN)
!
!***********************************************************************
! BIEF   V6P1                                   21/08/2010
!***********************************************************************
!
!brief    PRINTS OUT A VECTOR ON THE LISTING.
!
!warning  PRINTS OUT FOR REGULAR GRIDS ARE NOT IMPLEMENTED
!
!history  J-M HERVOUET
!+        17/08/94
!+        V5P1
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
!| NOM            |-->| NAME OF VECTOR (OR SIMPLY A COMMENT).
!| NPOIN          |-->| NUMBER OF POINTS.
!| VEC            |-->| VECTOR TO BE PRINTED.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
      INTEGER, INTENT(IN)           :: NPOIN
      DOUBLE PRECISION, INTENT(IN)  :: VEC(NPOIN)
      CHARACTER(LEN=32), INTENT(IN) :: NOM
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
      INTEGER IPOIN
!
!-----------------------------------------------------------------------
!
      IF(LNG.EQ.1.OR.LNG.EQ.2) WRITE(LU,'(//,1X,A32,//)') NOM
!
      IF(LNG.EQ.1.OR.LNG.EQ.2) THEN
        WRITE(LU,40) (IPOIN,VEC(IPOIN),IPOIN=1,NPOIN)
40      FORMAT(7(1X,1I5,':',G13.5))
      ENDIF
!
!-----------------------------------------------------------------------
!
      RETURN
      END
