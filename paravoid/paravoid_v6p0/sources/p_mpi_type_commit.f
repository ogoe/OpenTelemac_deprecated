C                       *****************************
                        SUBROUTINE  P_MPI_TYPE_COMMIT
C                       *****************************
C
     *(I1,I2)
C
C***********************************************************************
C  PARAVOID VERSION 6.0         27/10/09        C. DENIS (SINETICS)
C
C***********************************************************************
C
C      FONCTIONS:
C      ==========
C      FAUX APPEL A LA FONCTION MPI_TYPE_COMMIT
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
      INTEGER, INTENT(IN) :: I1,I2 
C
C-----------------------------------------------------------------------
C
      IF(LNG.EQ.1) WRITE(LU,*) 'APPEL DE P_MPI_TYPE_COMMIT VERSION VIDE'
      IF(LNG.EQ.2) WRITE(LU,*) 'CALL OF P_MPI_TYPE_COMMIT VOID VERSION' 
C
C-----------------------------------------------------------------------
C
      STOP
      END
 
