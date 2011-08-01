C                       **********************************
                        DOUBLE PRECISION FUNCTION BIEF_SUM
C                       **********************************
C
     *( X )
C
C***********************************************************************
C BIEF VERSION 6.0      11/05/2009    J-M HERVOUET (LNHE) 01 30 87 80 18
C***********************************************************************
C
C   FONCTION : SOMME DES COMPOSANTES D'UN VECTEUR
C
C              X PEUT ETRE UN VECTEUR OU
C
C              UNE STRUCTURE DE BLOC DE VECTEURS EN NOMBRE ET
C
C              CARACTERISTIQUES IDENTIQUES
C
C              ATTENTION ||||
C
C              SI LES VECTEURS ONT UNE DEUXIEME DIMENSION
C              ELLE EST POUR L'INSTANT IGNOREE
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
C .________________.____.______________________________________________.
C |      NOM       |MODE|                   ROLE                       |
C |________________|____|_______________________________________________
C |      X         | -->| LA STRUCTURE DONT ON VEUT LA SOMME
C |________________|____|_______________________________________________
C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
C-----------------------------------------------------------------------
C
C PROGRAMMES APPELES   :
C
C***********************************************************************
C
      USE BIEF, EX_BIEF_SUM => BIEF_SUM
C
      IMPLICIT NONE
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU
C
C-----------------------------------------------------------------------
C
C     STRUCTURES : VECTEURS OU BLOCS
C
      TYPE(BIEF_OBJ), INTENT(IN) :: X
C
C-----------------------------------------------------------------------
C
C  CAS D'UN VECTEUR
C
      IF(X%TYPE.EQ.2) THEN
C
        BIEF_SUM = SOMME(X%R,X%DIM1*X%DIM2)
C
C-----------------------------------------------------------------------
C
C  CAS OU LES STRUCTURES SONT DES BLOCS (A PROGRAMMER)
C
C     ELSEIF(X%TYPE.EQ.4) THEN
C
C-----------------------------------------------------------------------
C
C  ERREUR
C
      ELSE
C
         IF (LNG.EQ.1) WRITE(LU,50) X%NAME,X%TYPE
         IF (LNG.EQ.1) WRITE(LU,53)
         IF (LNG.EQ.2) WRITE(LU,60) X%NAME,X%TYPE
         IF (LNG.EQ.2) WRITE(LU,63)
50       FORMAT(1X,'BIEF_SUM (BIEF) : NOM DE X : ',A6,'  TYPE : ',1I6)
53       FORMAT(1X,'                  CAS NON PREVU')
60       FORMAT(1X,'BIEF_SUM (BIEF): NAME OF X : ',A6,'  TYPE : ',1I6)
63       FORMAT(1X,'                 CASE NOT IMPLEMENTED')
         CALL PLANTE(1)
         STOP
C
      ENDIF
C
C-----------------------------------------------------------------------
C
      RETURN
      END
