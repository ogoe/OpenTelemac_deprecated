C                       *****************************
                        SUBROUTINE GET_MPI_PARAMETERS
C                       *****************************
C
     *(P_INTEGER,P_REAL8,P_UB,P_COMM_WORLD,P_SUCCESS)
C
C***********************************************************************
C  PARALLEL VERSION 6.0    02/02/2009      J-M HERVOUET (01 30 87 80 18)
C             
C***********************************************************************
C
C      FONCTIONS:
C      ==========
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |                | -->|
C |________________|____|______________________________________________|
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C APPELE PAR : PRINCI
C
C SOUS-PROGRAMMES APPELES : NEANT
C
C**********************************************************************
C
      IMPLICIT NONE
C
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(OUT) :: P_INTEGER,P_REAL8,P_UB
      INTEGER, INTENT(OUT) :: P_COMM_WORLD,P_SUCCESS
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INCLUDE 'mpif.h'
C        
      P_INTEGER   =MPI_INTEGER
C     P_REAL8     =MPI_REAL8
C     CHRISTOPHE DENIS + JMH ON 04/12/2009 (FOR BLUE GENE) 
      P_REAL8     =MPI_DOUBLE_PRECISION
      P_UB        =MPI_UB
      P_COMM_WORLD=MPI_COMM_WORLD
      P_SUCCESS   =MPI_SUCCESS           
C
C-----------------------------------------------------------------------
C
      RETURN
      END
