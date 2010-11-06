C                       ******************
                        SUBROUTINE P_IWRIT
C                       ******************
C
     *(BUFFER,NBYTES,DEST,ITAG,IREQ)
C
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
      INTEGER NBYTES,DEST,ITAG,IREQ,IER
      DOUBLE PRECISION BUFFER(*)
C
C-----------------------------------------------------------------------
C
      IF(LNG.EQ.1) WRITE(LU,*) 'APPEL DE P_IWRIT VERSION VIDE'
      IF(LNG.EQ.2) WRITE(LU,*) 'CALL OF P_IWRIT IN ITS VOID VERSION'
C
C----------------------------------------------------------------------
C
      RETURN
      END
