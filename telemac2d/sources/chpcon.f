!                    *****************
                     SUBROUTINE CHPCON
!                    *****************
!
     &(UCONV,VCONV,U,V,UN,VN,TETAU)
!
!***********************************************************************
! TELEMAC2D   V6P1                                   21/08/2010
!***********************************************************************
!
!brief    COMPUTES THE ADVECTION VECTOR FIELD UCONV,VCONV.
!
!history  J-M HERVOUET (LNH)
!+        17/08/1994
!+        V5P2
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
!| TETAU          |-->| IMPLICITATION COEFFICIENT ON VELOCITY
!| U              |-->| X-COMPONENT OF VELOCITY AT TIME N+1
!| V              |-->| Y-COMPONENT OF VELOCITY AT TIME N+1
!| UCONV          |-->| X-COMPONENT OF ADVECTION FIELD
!| VCONV          |-->| Y-COMPONENT OF ADVECTION FIELD
!| UN             |-->| X-COMPONENT OF VELOCITY AT TIME N
!| VN             |-->| Y-COMPONENT OF VELOCITY AT TIME N
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      USE BIEF
!
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
      DOUBLE PRECISION, INTENT(IN)  :: TETAU
      TYPE(BIEF_OBJ), INTENT(IN)    :: U,UN,V,VN
      TYPE(BIEF_OBJ), INTENT(INOUT) :: UCONV,VCONV
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
!     CREATES CONVECTION ARRAYS UCONV AND VCONV
!
      CALL OS( 'X=CY    ' , UCONV , UN , U , 1.D0-TETAU )
      CALL OS( 'X=X+CY  ' , UCONV , U  , U ,      TETAU )
      CALL OS( 'X=CY    ' , VCONV , VN , U , 1.D0-TETAU )
      CALL OS( 'X=X+CY  ' , VCONV , V  , U ,      TETAU )
!
!-----------------------------------------------------------------------
!
      RETURN
      END
