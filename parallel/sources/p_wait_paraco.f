C                       ************************
                        SUBROUTINE P_WAIT_PARACO
C                       ************************
C
     *(IBUF,NB)
C
C***********************************************************************
C  PARA       VERSION 5.9         23/06/2008        PASCAL VEZOLLE (IBM)
C  
C***********************************************************************
C
C      FONCTIONS: ATTENTE A LA FIN DE PARACO
C      ==========
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE
C |________________|____|______________________________________________|
C |                | -->|
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
      INTEGER IBUF(*), NB, IER
      INTEGER WAIT_REQ(MPI_STATUS_SIZE,100)
      SAVE
C
C-----------------------------------------------------------------------
C
      IF(NB.GT.100) THEN
        WRITE(LU,*) 'WAIT_PARACO:'
        WRITE(LU,*) 'DIMENSION OF WAIT_REQ TOO SMALL'
        STOP
      ENDIF
C
      CALL MPI_WAITALL(NB,IBUF,WAIT_REQ,IER)
C
      IF(IER.NE.0) THEN
        WRITE(LU,*) 'WAIT_PARACO:'
        WRITE(LU,*) 'MPI ERROR ',IER
        STOP
      ENDIF
C
C----------------------------------------------------------------------
C
      RETURN
      END
