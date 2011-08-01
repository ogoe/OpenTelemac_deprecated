C                       ********************
                        SUBROUTINE P_IWRIT_C
C                       ********************
C
     *(BUFFER,NBYTES,DEST,ITAG,IREQ)
C
C***********************************************************************
C  PARA       VERSION 5.9       23/06/2008       PASCAL VEZOLLES (IBM)
C***********************************************************************
C
C      FONCTIONS: ENVOI DE VALEURS ENTRE PROCESSEURS.
C      ==========
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |  BUFFER        | -->| ZONE TAMPON POUR LES DONNEES
C |                |    | BUFFER / PUFFERFELD
C |  NBYTES        | -->| NOMBRE DE BYTES A TRANSMETTRE
C |                |    | LENGTH IN BYTES / LAENGE IN BYTES
C |  DEST          | -->| DESTINATION DES DONNEES
C |                |    | TID OF THE DEST.  / KNOTEN-ID DES EMPFAENGERS
C |  ITAG          | -->| MESSAGE TAG
C |  IREQ          | -->| NUMERO DE REQUEST POUR MPI_ISEND
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
      INTEGER NBYTES,DEST,ITAG,IREQ,IER
      CHARACTER(LEN=*) BUFFER
C
C-----------------------------------------------------------------------
C
      CALL MPI_ISEND(BUFFER,NBYTES,MPI_BYTE,DEST,ITAG,
     &               MPI_COMM_WORLD,IREQ,IER)
C
      IF(IER.NE.0) THEN
        IF(LNG.EQ.1) WRITE(LU,*) 'P_IWRIT: ERREUR IN MPI_ISEND'
        IF(LNG.EQ.2) WRITE(LU,*) 'P_IWRIT: ERROR IN MPI_ISEND'
        WRITE(LU,*) 'MPI ERROR ',IER
        STOP
      ENDIF
C
C----------------------------------------------------------------------
C
      RETURN
      END
