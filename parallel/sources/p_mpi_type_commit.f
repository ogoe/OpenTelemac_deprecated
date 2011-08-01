C                       *****************************
                        SUBROUTINE  P_MPI_TYPE_COMMIT
C                       *****************************
C
     *(I1,I2)
C
C***********************************************************************
C  PARA VERSION 6.0         27/10/09        C. DENIS (SINETICS)
C
C***********************************************************************
C
C      FONCTIONS: CALLING MPI_TYPE_COMMIT
C      ==========
C
C-----------------------------------------------------------------------
C
C APPELE PAR :
C
C SOUS-PROGRAMMES APPELES : NEANT
C
C**********************************************************************
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN) :: I1,I2
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      CALL MPI_TYPE_COMMIT(I1,I2)
C
      IF(I2.NE.0) THEN
        WRITE(LU,*) 'P_MPI_TYPE_COMMIT:'
        WRITE(LU,*) 'MPI ERROR ',I2
        STOP
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END
