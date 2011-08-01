C                       ********************************
                        DOUBLE PRECISION FUNCTION P_DSUM
C                       ********************************
C
     *(MYPART)
C
C***********************************************************************
C  PARALLEL    VERSION 5.9         23/06/2008      J-M HERVOUET (LNHE)
C  ADAPTED FOR MPI              /10/99       RAINER JOHANNI (SGI MUNICH)
C  VERSION 5.0 MODIFIED       28/12/99    J.A. JANKOWSKI (BAW KARLSRUHE)
C***********************************************************************
C
C      FONCTIONS: SOMME D'UNE VALEUR SUR TOUS LES PROCESSEURS.
C      ==========
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE
C |________________|____|______________________________________________|
C |  MYPART        | -->| CONTRIBUTION DU PROCESSEUR APPELANT.
C |________________|____|______________________________________________|
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C
C-----------------------------------------------------------------------
C
C APPELE PAR :
C
C SOUS-PROGRAMMES APPELES : MPI_ALLREDUCE
C
C**********************************************************************
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
      INCLUDE 'mpif.h'
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      DOUBLE PRECISION, INTENT(IN) :: MYPART
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER IER
C
C-----------------------------------------------------------------------
C
      CALL MPI_ALLREDUCE(MYPART,P_DSUM,1,MPI_DOUBLE_PRECISION,            
     &                   MPI_SUM,MPI_COMM_WORLD,IER)
C
      IF(IER.NE.0) THEN
        IF(LNG.EQ.1) WRITE(LU,*) 'P_DSUM: ERREUR DANS MPI_ALLREDUCE'
        IF(LNG.EQ.2) WRITE(LU,*) 'P_DSUM: ERROR IN MPI_ALLREDUCE'
        WRITE(LU,*) 'MPI ERROR: ',IER
        STOP 
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END
