C                       *************************
                        SUBROUTINE  P_MPI_TYPE_UB
C                       *************************
C
     *(I1,I2,I3)
C
C***********************************************************************
C  PARA VERSION 6.0         27/10/09        C. DENIS (SINETICS)
C
C***********************************************************************
C
C      FONCTIONS:  CALLING MPI_TYPE_UB
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
      INTEGER, INTENT(IN) :: I1,I2,I3
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      CALL MPI_TYPE_UB(I1,I2,I3)
C
      IF(I3.NE.0) THEN
        WRITE(LU,*) 'P_MPI_TYPE_UB:'
        WRITE(LU,*) 'MPI ERROR ',I3
        STOP
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END