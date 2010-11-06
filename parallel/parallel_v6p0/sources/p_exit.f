C                       *****************
                        SUBROUTINE P_EXIT
C                       *****************
C
C***********************************************************************
C  PARALLEL    VERSION 6.0         16/06/2009       J-M HERVOUET (LNHE)
C  ADAPTED FOR MPI              /10/99       RAINER JOHANNI (SGI MUNICH)
C  VERSION 5.0 MODIFIED       28/12/99    J.A. JANKOWSKI (BAW KARLSRUHE)
C***********************************************************************
C
C      FONCTIONS: FIN DE MPI
C      ==========
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
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
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
      INCLUDE 'mpif.h'
C
      INTEGER IER
C
C-----------------------------------------------------------------------
C
      WRITE(LU,*) ' '
      IF(LNG.EQ.1) WRITE(LU,*) 'SORTIE DE MPI'
      IF(LNG.EQ.2) WRITE(LU,*) 'EXITING MPI'
      WRITE(LU,*) ' '
C
C     TO AVOID EXITING BEFORE EVERYTHING IS DONE IN OTHER PROCESSORS
C
      CALL MPI_BARRIER(MPI_COMM_WORLD,IER)
C
C     EXITING
C
      CALL MPI_FINALIZE(IER)
      IF(IER.NE.0) THEN
        IF(LNG.EQ.1) WRITE(LU,*) 'P_EXIT: ERREUR DANS MPI_FINALIZE'
        IF(LNG.EQ.2) WRITE(LU,*) 'P_EXIT: ERROR IN MPI_FINALIZE'
        WRITE(LU,*) 'MPI ERROR ',IER
        STOP 
      ENDIF 
C
C-----------------------------------------------------------------------
C
      RETURN
      END
