C                       ********************
                        SUBROUTINE P_IREAD_C
C                       ********************
C
     *(BUFFER,NBYTES,SOURCE,ITAG,IREQ)
C
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
      INTEGER NBYTES,SOURCE,ITAG,IREQ,IER
      CHARACTER(LEN=*) BUFFER
C
C-----------------------------------------------------------------------
C
      IF(LNG.EQ.1) WRITE(LU,*) 'APPEL DE P_IREAD VERSION VIDE'
      IF(LNG.EQ.2) WRITE(LU,*) 'CALL OF P_IREAD IN ITS VOID VERSION'
C
C-----------------------------------------------------------------------
C
      RETURN
      END
