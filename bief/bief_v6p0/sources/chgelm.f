!                      *****************
                       SUBROUTINE CHGELM
!                      *****************
     &(V, IELM)
!
!=======================================================================
! BIEF VERSION 5.6    MARCH 1999        JACEK A. JANKOWSKI PINXIT
!=======================================================================
! FUNCTION:
! CHANGES THE TYPE OF ELEMENT AND THE DIMENSION OF A VECTOR
! NOTE: ONLY THE FIRST VECTOR DIMENSION IS TREATED!
!
!-----------------------------------------------------------------------
!
      USE BIEF, EX_CHGELM => CHGELM
      IMPLICIT NONE 
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN)           :: IELM
      TYPE(BIEF_OBJ), INTENT(INOUT) :: V
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      IF(V%TYPE.NE.2) THEN
       IF(LNG.EQ.1) THEN 
        WRITE(LU,*) 'CHGELM : L''OBJET ',V%NAME,' N''EST PAS UN VECTEUR'
        WRITE(LU,*) 'V%TYPE = ',V%TYPE
       ENDIF
       IF(LNG.EQ.2) THEN 
        WRITE(LU,*) 'CHGELM: OBJECT ',V%NAME,' IS NOT A VECTOR.'
        WRITE(LU,*) 'V%TYPE = ',V%TYPE
       ENDIF
       CALL PLANTE(1)
       STOP
      ELSE 
       V%ELM = IELM
       IF(V%MAXDIM1.LT.NBPTS(IELM)) V%MAXDIM1 = NBPTS(IELM)
       V%DIM1 = NBPTS(IELM) 
      ENDIF
!
!-----------------------------------------------------------------------
!
      RETURN
      END
