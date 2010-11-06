C                       *****************
                        SUBROUTINE IMPVEC
C                       *****************
C
     *(VEC,NOM,NPOIN)
C
C-----------------------------------------------------------------------
C BIEF VERSION 5.1          17/08/94    J-M HERVOUET 30 71 80 18
C-----------------------------------------------------------------------
C
C      FONCTION : IMPRESSION D'UN VECTEUR SUR LISTING
C
C                 L'IMPRESSION POUR UN MAILLAGE REGULIER
C                 N'EST PAS PROGRAMMEE.
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|______________________________________________|
C |   VEC          | -->|  VECTEUR A IMPRIMER.
C |   NOM          | -->|  NOM DU VECTEUR OU COMMENTAIRE
C |   NPOIN        | -->|  NOMBRE DE POINTS DU MAILLAGE
C |________________|____|_______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMME APPELANT : DESIMP
C
C**********************************************************************
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER, INTENT(IN)           :: NPOIN
      DOUBLE PRECISION, INTENT(IN)  :: VEC(NPOIN)
      CHARACTER(LEN=32), INTENT(IN) :: NOM
C
C+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C
      INTEGER IPOIN
C
C-----------------------------------------------------------------------
C
      IF(LNG.EQ.1.OR.LNG.EQ.2) WRITE(LU,'(//,1X,A32,//)') NOM
C
      IF(LNG.EQ.1.OR.LNG.EQ.2) THEN
        WRITE(LU,40) (IPOIN,VEC(IPOIN),IPOIN=1,NPOIN)
40      FORMAT(7(1X,1I5,':',G13.5))
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END
 
