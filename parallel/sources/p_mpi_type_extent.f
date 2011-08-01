C                       *****************************
                        SUBROUTINE  P_MPI_TYPE_EXTENT
C                       *****************************
C
     *(I1,I2,I3)
C
C***********************************************************************
C  PARA VERSION 6.0         27/10/09        C. DENIS (SINETICS)
C
C***********************************************************************
C
C      FONCTIONS:  CALLING MPI_TYPE_EXTENT
C      ==========
C     
C-----------------------------------------------------------------------
C
C
C APPELE PAR :
C
C SOUS-PROGRAMMES APPELES : NEANT
C
C**********************************************************************
C
      IMPLICIT NONE
      INTEGER LNG,LU,IERR
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(INOUT) :: I1,I2,I3
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      CALL MPI_TYPE_EXTENT(I1,I2,IERR)
C
      IF(IERR.NE.0) THEN
        WRITE(LU,*) 'P_MPI_TYPE_EXTENT:'
        WRITE(LU,*) 'MPI ERROR ',IERR
        STOP
      ENDIF
C     
C-----------------------------------------------------------------------
C
      RETURN
      END
