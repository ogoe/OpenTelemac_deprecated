C                       ********************************
                        DOUBLE PRECISION FUNCTION P_DMIN
C                       ********************************
C
     *(MYPART)
C
C***********************************************************************
C  PARA       VERSION 5.9         08/01/97        J-M HERVOUET (LNH)
C  ADAPTED FOR MPI              /10/99       RAINER JOHANNI (SGI MUNICH)
C  VERSION 5.0 MODIFIED       28/12/99    J.A. JANKOWSKI (BAW KARLSRUHE)
C***********************************************************************
C
C      FONCTIONS: MINIMUM D'UNE VALEUR ENTRE TOUS LES PROCESSEURS.
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
      CALL MPI_ALLREDUCE(MYPART,P_DMIN,1,MPI_DOUBLE_PRECISION,MPI_MIN,
     &                   MPI_COMM_WORLD,IER)
C
      IF(IER.NE.0) THEN
        IF(LNG.EQ.1) WRITE(LU,*) 'P_DMIN: ERREUR DANS MPI_ALLREDUCE'
        IF(LNG.EQ.2) WRITE(LU,*) 'P_DMIN: ERROR IN MPI_ALLREDUCE'
        WRITE(LU,*) 'MPI ERROR ',IER
        STOP
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END
