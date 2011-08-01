C                       *****************
                        SUBROUTINE P_WRIT
C                       *****************
C
     *(BUFFER,NBYTES,DEST,TYPE)
C
C***********************************************************************
C  PARA       VERSION 5.9       23/06/2008       HANS HERRMANN (HANOVRE)
C             MODIFIED        08/06/96     REINHARD HINKELMANN (HANOVRE)
C             MODIFIED        17/12/96            J-M HERVOUET (LNH)
C  ADAPTED FOR MPI              /10/99       RAINER JOHANNI (SGI MUNICH)
C  VERSION 5.0 MODIFIED       28/12/99    J.A. JANKOWSKI (BAW KARLSRUHE)
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
C |  TYPE          | -->| TYPE DES DONNEES (MSGTAG DE PVM)
C |                |    |   0 - STRING
C |                |    |   1 - BYTE1
C |                |    |   2 - INTEGER2
C |                |    |   3 - INTEGER4
C |                |    |   4 - REAL4
C |                |    |   5 - COMPLEX8
C |                |    |   6 - REAL8
C |                |    |   7 - COMPLEX16
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
      INTEGER NBYTES,DEST,TYPE,IER
      DOUBLE PRECISION BUFFER(*)
C
C-----------------------------------------------------------------------
C
      CALL MPI_SEND(BUFFER,NBYTES,MPI_BYTE,DEST,4711,
     &              MPI_COMM_WORLD,IER)
C
      IF (IER.NE.0) THEN
        IF(LNG.EQ.1) WRITE(LU,*) 'P_WRITE: ERREUR IN MPI_SEND'
        IF(LNG.EQ.2) WRITE(LU,*) 'P_WRITE: ERROR IN MPI_SEND'
        WRITE(LU,*) 'MPI ERROR ',IER
        STOP
      ENDIF
C
C----------------------------------------------------------------------
C
      RETURN
      END