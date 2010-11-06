C                       ***************************
                        SUBROUTINE  P_MPI_ALLTOALLV
C                       ***************************
C
     *(I1,I2,I3,I4,I5,I6,I7,I8,I9,I10)
C
C***********************************************************************
C  PARA VERSION 6.0         27/10/09                 C. DENIS (SINETICS)
C
C***********************************************************************
C
C      FONCTIONS:   CALLING MPI_ALLTOALLV
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
      INTEGER, INTENT(IN) :: I1(*),I2(*),I3(*),I4,I5(*),I6(*),I7(*)
      INTEGER, INTENT(IN) :: I8,I9,I10
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      CALL MPI_ALLTOALLV(I1,I2,I3,I4,I5,I6,I7,I8,I9,I10)
C
      IF(I10.NE.0) THEN
        WRITE(LU,*) 'P_MPI_ALLTOALLV:'
        WRITE(LU,*) 'MPI ERROR ',I10
        STOP
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END
      


 
     
 
