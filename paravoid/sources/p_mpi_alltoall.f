C                       *************************
                        SUBROUTINE P_MPI_ALLTOALL
C                       *************************
C
     *(I1,I2,I3,I4,I5,I6,I7,I8)
C
C***********************************************************************
C  PARAVOID VERSION 6.0         27/10/09        C. DENIS (SINETICS)
C
C***********************************************************************
C
C      FONCTIONS:
C      ==========
C
C      FAUX APPEL A LA FONCTION MPI_ALLTOALL
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
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
      INTEGER, INTENT(IN) :: I1(*),I2,I3,I4(*),I5,I6,I7,I8
C
C-----------------------------------------------------------------------
C
      IF(LNG.EQ.1) WRITE(LU,*) 'APPEL DE P_MPI_ALLTOALL VERSION VIDE'
      IF(LNG.EQ.2) WRITE(LU,*) 'CALL OF P_MPI_ALLTOALL VOID VERSION'
C
C-----------------------------------------------------------------------
C
      STOP
      END
 
