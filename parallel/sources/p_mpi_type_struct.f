C                       *****************************
                        SUBROUTINE  P_MPI_TYPE_STRUCT
C                       *****************************
C
     *(I1,I2,I3,I4,I5,I6)
C
C***********************************************************************
C  PARA VERSION 6.0         27/10/09        C. DENIS (SINETICS)
C
C***********************************************************************
C
C      FONCTIONS:  CALLING MPI_TYPE_STRUCT
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
      INTEGER, INTENT(IN) :: I1,I5,I6
      INTEGER, INTENT(IN) :: I2(I1),I3(I1),I4(I1)
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      CALL MPI_TYPE_STRUCT(I1,I2,I3,I4,I5,I6)
C
      IF(I6.NE.0) THEN
        WRITE(LU,*) 'P_MPI_TYPE_STRUCT:'
        WRITE(LU,*) 'MPI ERROR ',I6
        STOP
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END

