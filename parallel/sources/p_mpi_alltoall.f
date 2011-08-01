C                       *************************
                        SUBROUTINE P_MPI_ALLTOALL
C                       *************************
C
     *(I1,I2,I3,I4,I5,I6,I7,I8)
C
C***********************************************************************
C  PARA VERSION 6.0         27/10/09                 C. DENIS (SINETICS)
C
C***********************************************************************
C
C      FONCTIONS:   CALLING MPI_ALLTOALL
C      ==========
C    
C-----------------------------------------------------------------------
C
C APPELE PAR :
C
C SOUS-PROGRAMMES APPELES : NEANT
C
C***********************************************************************
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN) :: I1(*),I2,I3,I4(*),I5,I6,I7,I8
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      CALL MPI_ALLTOALL(I1,I2,I3,I4,I5,I6,I7,I8)
C
      IF(I8.NE.0) THEN
        WRITE(LU,*) 'P_MPI_ALLTOALL:'
        WRITE(LU,*) 'MPI ERROR ',I8
        STOP
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END
