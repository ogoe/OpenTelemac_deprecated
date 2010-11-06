C                       ********************
                        SUBROUTINE P_IREAD_C
C                       ********************
C
     *(BUFFER,NBYTES,SOURCE,ITAG,IREQ)
C
C***********************************************************************
C  PARA       VERSION 5.9       23/06/2008       PASCAL VEZOLLES (IBM)
C***********************************************************************
C
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |  BUFFER        | -->| ZONE TAMPON POUR LES DONNEES
C |                |    | BUFFER / PUFFERFELD
C |  NBYTES        | -->| NOMBRE DE BYTES A TRANSMETTRE
C |                |    | LENGTH IN BYTES / LAENGE IN BYTES
C |  SOURCE        | -->| ORIGINE DES DONNEES
C |                |    | TID OF THE SENDER / KNOTEN-ID DES SENDER
C |  ITAG          | -->| MESSAGE TAG
C |                |    |
C |  IREQ          | -->| NUMERO DE REQUEST POUR MPI_IRECV
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
      INTEGER NBYTES,SOURCE,ITAG,IREQ,IER
      CHARACTER(LEN=*) BUFFER
C
C-----------------------------------------------------------------------
C     RECEPTION DES DONNEES / RECEIVING DATA / DATEN EMPFANGEN
C-----------------------------------------------------------------------
C
      CALL MPI_IRECV(BUFFER,NBYTES,MPI_BYTE,SOURCE,ITAG,
     &               MPI_COMM_WORLD,IREQ,IER)
C
      IF(IER.NE.0) THEN
        IF(LNG.EQ.1) WRITE(LU,*) 'P_IREAD: ERREUR IN MPI_IRECV'
        IF(LNG.EQ.2) WRITE(LU,*) 'P_IREAD: ERROR IN MPI_IRECV'
        WRITE(LU,*) 'MPI ERROR ',IER
        STOP
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END
