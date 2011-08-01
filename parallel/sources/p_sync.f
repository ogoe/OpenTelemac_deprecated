C                       *****************
                        SUBROUTINE P_SYNC
C                       *****************
C
C
C***********************************************************************
C  PARA       VERSION 5.9         23/06/2008     HANS HERRMANN (HANOVRE)
C             MODIFIED        08/06/96     REINHARD HINKELMANN (HANOVRE)
C             MODIFIED        17/12/96            J-M HERVOUET (LNH)
C  ADAPTED FOR MPI              /10/99       RAINER JOHANNI (SGI MUNICH)
C  VERSION 5.0 MODIFIED       28/12/99    J.A. JANKOWSKI (BAW KARLSRUHE)
C***********************************************************************
C
C      FONCTIONS: SYNCHRONISATION DE TOUS LES PROCESSEURS
C      ==========
C
C      VARIABLEN
C      ---------
C      IPID     : NUMMER DES KNOTENS / NUMERO DU PROCESSEUR
C      NCSIZE   : PROZESSORENANZAHL  / NOMBRE DE PROCESSEURS
C      ITID     : PROZESSNUMMERNFELD (WIRD AUF NCUBE NICHT BENOETIGT)
C
C
C      WICHTIG:
C      DIESE ROUTINE MUSS VON ALLEN KNOTEN/PROZESSEN AUFGERUFEN
C      WERDEN, SONST "STEHT" DAS PROGRAMM.
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |  NPROC         | -->| NOMBRE DE PROCESSEURS OU STATIONS
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
      CALL MPI_BARRIER(MPI_COMM_WORLD,IER)
C
      IF(IER.NE.0) THEN
        IF(LNG.EQ.1) WRITE(LU,*) 'P_SYNC : ERREUR DANS MPI_BARRIER'
        IF(LNG.EQ.2) WRITE(LU,*) 'P_SYNC: ERROR IN MPI_BARRIER'
        WRITE(LU,*) 'MPI ERROR ',IER
        STOP
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END
