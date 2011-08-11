!                    ****************
                     SUBROUTINE OPASS
!                    ****************
!
     &(OP,X,W,IW,Y,IY,LIMVOI,MXPTVS,NPMAX)
!
!***********************************************************************
! BIEF   V6P1                                   21/08/2010
!***********************************************************************
!
!brief    FRONTAL MATRIX-VECTOR PRODUCT FOR ELEMENT 11-11.
!+
!+            OMITS THE DIAGONAL TERMS HERE.
!
!history  J-M HERVOUET (LNH)    ; F  LEPEINTRE (LNH)
!+        05/02/91
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
!| IW             |-->| ADDRESSES IN W
!| IY             |-->| ADDRESSES IN Y
!| LIMVOI         |-->| LIMVOI(N,1): LOWER RANK OF POINTS WITH N NEIGHBOURS
!|                |   | LIMVOI(N,2): UPPER RANK OF POINTS WITH N NEIGHBOURS
!| MXPTVS         |-->| MAXIMUM NUMBER OF NEIGHBOURS OF A POINT
!| NPMAX          |-->| MAXIMUM NUMBER OF POINTS
!| OP             |-->| OPERATION TO BE DONE (SEE ABOVE)
!| W              |-->| OFF-DIAGONAL TERMS OF MATRIX
!| X              |<->| FINAL ASSEMBLED VECTOR
!| Y              |-->| X=AY WHERE OFF-DIAGONAL TERMS OF A ARE HERE W
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
      INTEGER, INTENT(IN) :: NPMAX,MXPTVS
      INTEGER, INTENT(IN) :: IW(NPMAX,*),IY(NPMAX,*),LIMVOI(MXPTVS,2)
      DOUBLE PRECISION, INTENT(INOUT) :: X(*)
      DOUBLE PRECISION, INTENT(IN)    :: W(*),Y(*)
      CHARACTER(LEN=8), INTENT(IN)    :: OP
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
      INTEGER I
!
!-----------------------------------------------------------------------
!
      IF(MXPTVS.GT.11) THEN
       IF(LNG.EQ.1) WRITE(LU,777)
       IF(LNG.EQ.2) WRITE(LU,778)
777    FORMAT(1X,'OPASS (BIEF) : PROGRAMME JUSQU''A 11 VOISINS',/,
     &        1X,'CHOISIR STOCKAGE DES MATRICES : 1')
778    FORMAT(1X,'OPASS (BIEF): IMPLEMENTED UP TO 11 NEIGHBOURS ONLY',/,
     &        1X,'CHOOSE STORAGE OF MATRICES : 1')
       CALL PLANTE(1)
       STOP
      ENDIF
!
!-----------------------------------------------------------------------
!
      IF(OP(1:8).EQ.'X=WY    ') THEN
!
      IF(LIMVOI(1,1).GT.0) THEN
!       THIS CASE IS NOT POSSIBLE IN THEORY
        DO I = LIMVOI(1,1) , LIMVOI(1,2)
           X(I) = W(IW(I,1))*Y(IY(I,1))
        ENDDO
      ENDIF
!
      IF(LIMVOI(2,1).GT.0) THEN
!       THIS CASE ONLY EXISTS IF THERE ARE OVERSTRESSED TRIANGLES
        DO I = LIMVOI(2,1) , LIMVOI(2,2)
           X(I) = W(IW(I,1))*Y(IY(I,1)) + W(IW(I,2))*Y(IY(I,2))
        ENDDO
      ENDIF
!
      IF(MXPTVS.GE.3.AND.LIMVOI(3,1).GT.0) THEN
      DO I = LIMVOI(3,1) , LIMVOI(3,2)
         X(I) = W(IW(I,1))*Y(IY(I,1)) + W(IW(I,2))*Y(IY(I,2))
     &        + W(IW(I,3))*Y(IY(I,3))
      ENDDO
      ENDIF
!
      IF(MXPTVS.GE.4.AND.LIMVOI(4,1).GT.0) THEN
      DO I = LIMVOI(4,1) , LIMVOI(4,2)
         X(I) = W(IW(I,1))*Y(IY(I,1)) + W(IW(I,2))*Y(IY(I,2))
     &        + W(IW(I,3))*Y(IY(I,3)) + W(IW(I,4))*Y(IY(I,4))
      ENDDO
      ENDIF
!
      IF(MXPTVS.GE.5.AND.LIMVOI(5,1).GT.0) THEN
      DO I = LIMVOI(5,1) , LIMVOI(5,2)
         X(I) = W(IW(I,1))*Y(IY(I,1)) + W(IW(I,2))*Y(IY(I,2))
     &        + W(IW(I,3))*Y(IY(I,3)) + W(IW(I,4))*Y(IY(I,4))
     &        + W(IW(I,5))*Y(IY(I,5))
      ENDDO
      ENDIF
!
      IF(MXPTVS.GE.6.AND.LIMVOI(6,1).GT.0) THEN
      DO I = LIMVOI(6,1) , LIMVOI(6,2)
         X(I) = W(IW(I,1))*Y(IY(I,1)) + W(IW(I,2))*Y(IY(I,2))
     &        + W(IW(I,3))*Y(IY(I,3)) + W(IW(I,4))*Y(IY(I,4))
     &        + W(IW(I,5))*Y(IY(I,5)) + W(IW(I,6))*Y(IY(I,6))
      ENDDO
      ENDIF
!
      IF(MXPTVS.GE.7.AND.LIMVOI(7,1).GT.0) THEN
      DO I = LIMVOI(7,1) , LIMVOI(7,2)
         X(I) = W(IW(I,1))*Y(IY(I,1)) + W(IW(I,2))*Y(IY(I,2))
     &        + W(IW(I,3))*Y(IY(I,3)) + W(IW(I,4))*Y(IY(I,4))
     &        + W(IW(I,5))*Y(IY(I,5)) + W(IW(I,6))*Y(IY(I,6))
     &        + W(IW(I,7))*Y(IY(I,7))
      ENDDO
      ENDIF
!
      IF(MXPTVS.GE.8.AND.LIMVOI(8,1).GT.0) THEN
      DO I = LIMVOI(8,1) , LIMVOI(8,2)
         X(I) = W(IW(I,1))*Y(IY(I,1)) + W(IW(I,2))*Y(IY(I,2))
     &        + W(IW(I,3))*Y(IY(I,3)) + W(IW(I,4))*Y(IY(I,4))
     &        + W(IW(I,5))*Y(IY(I,5)) + W(IW(I,6))*Y(IY(I,6))
     &        + W(IW(I,7))*Y(IY(I,7)) + W(IW(I,8))*Y(IY(I,8))
      ENDDO
      ENDIF
!
      IF(MXPTVS.GE.9.AND.LIMVOI(9,1).GT.0) THEN
      DO I = LIMVOI(9,1) , LIMVOI(9,2)
         X(I) = W(IW(I,1))*Y(IY(I,1)) + W(IW(I,2))*Y(IY(I,2))
     &        + W(IW(I,3))*Y(IY(I,3)) + W(IW(I,4))*Y(IY(I,4))
     &        + W(IW(I,5))*Y(IY(I,5)) + W(IW(I,6))*Y(IY(I,6))
     &        + W(IW(I,7))*Y(IY(I,7)) + W(IW(I,8))*Y(IY(I,8))
     &        + W(IW(I,9))*Y(IY(I,9))
      ENDDO
      ENDIF
!
      IF(MXPTVS.GE.10.AND.LIMVOI(10,1).GT.0) THEN
      DO I = LIMVOI(10,1) , LIMVOI(10,2)
         X(I) = W(IW(I,1))*Y(IY(I,1)) + W(IW(I,2))*Y(IY(I,2))
     &        + W(IW(I,3))*Y(IY(I,3)) + W(IW(I,4))*Y(IY(I,4))
     &        + W(IW(I,5))*Y(IY(I,5)) + W(IW(I,6))*Y(IY(I,6))
     &        + W(IW(I,7))*Y(IY(I,7)) + W(IW(I,8))*Y(IY(I,8))
     &        + W(IW(I,9))*Y(IY(I,9)) + W(IW(I,10))*Y(IY(I,10))
      ENDDO
      ENDIF
!
      IF(MXPTVS.GE.11.AND.LIMVOI(11,1).GT.0) THEN
      DO I = LIMVOI(11,1) , LIMVOI(11,2)
         X(I) = W(IW(I,1))*Y(IY(I,1)) + W(IW(I,2))*Y(IY(I,2))
     &        + W(IW(I,3))*Y(IY(I,3)) + W(IW(I,4))*Y(IY(I,4))
     &        + W(IW(I,5))*Y(IY(I,5)) + W(IW(I,6))*Y(IY(I,6))
     &        + W(IW(I,7))*Y(IY(I,7)) + W(IW(I,8))*Y(IY(I,8))
     &        + W(IW(I,9))*Y(IY(I,9)) + W(IW(I,10))*Y(IY(I,10))
     &        + W(IW(I,11))*Y(IY(I,11))
      ENDDO
      ENDIF
!
      ELSE
!
        IF (LNG.EQ.1) WRITE(LU,3000) OP
        IF (LNG.EQ.2) WRITE(LU,3001) OP
3000    FORMAT(1X,'OPASS (BIEF) : OPERATION INCONNUE : ',A8)
3001    FORMAT(1X,'OPASS (BIEF) : UNKNOWN OPERATION : ',A8)
        CALL PLANTE(1)
        STOP
!
      ENDIF
!
!-----------------------------------------------------------------------
!
      RETURN
      END
